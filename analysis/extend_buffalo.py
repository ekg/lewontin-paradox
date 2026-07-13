#!/usr/bin/env python3
"""Extend Buffalo (2021): test gBGC-saturation against his cross-species data.

Tier 1 (n=173): model comparison on log10(pi) vs log10(Nc):
    M0  power law      a + b*x
    M1  quadratic      a + b*x + c*x^2   (concavity test)
    M2  Michaelis-Menten (pure saturation, X=1):  pi = Vmax*Nc/(K+Nc)
    M3  MM + linear tail (X<1):  pi = a*Nc + Vmax*Nc/(K+Nc)
  The gBGC-saturation model predicts M3 > M2 > M0, with a small but
  nonzero linear tail `a` (saturation is big, not the whole story).

Tier 2 (n=18-39): discriminating recombination sign test.
  residual = obs log10(pi) - power-law prediction.
  BGS  predicts residual ~ +recombination  (rec protects diversity).
  gBGC predicts residual ~ -recombination   (biased repair erodes diversity).
  The sign discriminates the two forces. We control for Nc.
"""
import sys, csv, math, os
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

DATA = Path(os.environ.get(
    "BUFFALO_DATA",
    Path(__file__).resolve().parent / "data" / "buffalo" / "combined_data.tsv",
))

def load():
    rows = list(csv.DictReader(open(DATA), delimiter="\t"))
    out = []
    for r in rows:
        try:
            pi = float(r["diversity"]); Nc = float(r["pred_log10_N"])
        except ValueError:
            continue
        if pi <= 0 or Nc <= 0 or not np.isfinite(pi) or not np.isfinite(Nc):
            continue
        rec = None
        try: rec = float(r["ave_rec"])
        except ValueError: pass
        ml = None
        try: ml = float(r["map_length"])
        except ValueError: pass
        gs = None
        try: gs = float(r["genome_size"])
        except ValueError: pass
        out.append(dict(sp=r["species"], pi=pi, logN=Nc, Nc=10**Nc,
                        logpi=math.log10(pi), ave_rec=rec, map_length=ml,
                        genome_size=gs))
    return out

def nloglik(y, yhat, k):
    resid = y - yhat
    s2 = np.var(resid)
    if s2 <= 0: s2 = 1e-12
    n = len(y)
    return 0.5*n*(math.log(2*math.pi*s2)+1) + (0.5*n*s2)/s2 * 0  # noqa
    # use proper gaussian nll with est variance
def aic(y, yhat, k):
    resid = y - yhat
    s2 = np.mean(resid**2)/len(resid)
    if s2 <= 0: s2 = 1e-12
    sigma = math.sqrt(s2)
    n = len(y)
    nll = -np.sum(stats.norm.logpdf(resid, scale=sigma))
    return 2*k + 2*nll, n - k - 1

def fit_power(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    r2 = 1 - np.sum((y-yhat)**2)/np.sum((y-y.mean())**2)
    return coef, yhat, r2

def fit_quad(x, y):
    A = np.vstack([np.ones_like(x), x, x**2]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    r2 = 1 - np.sum((y-yhat)**2)/np.sum((y-y.mean())**2)
    return coef, yhat, r2

def mm(Nc, Vmax, K):
    return Vmax*Nc/(K+Nc)
def mm_lin(Nc, a, Vmax, K):
    return a*Nc + Vmax*Nc/(K+Nc)

def main():
    rows = load()
    logN = np.array([r["logN"] for r in rows])
    logpi = np.array([r["logpi"] for r in rows])
    Nc = np.array([r["Nc"] for r in rows])
    pi = np.array([r["pi"] for r in rows])
    n = len(rows)
    print(f"== Tier 1: model comparison on pi vs Nc (n={n}) ==\n")

    # M0 power law (in log space)
    (a,b), y0, r0 = fit_power(logN, logpi)
    a0,nll0 = aic(logpi, y0, 2)
    print(f"M0 power law : log10 pi = {a:.3f} + {b:.3f}*log10 Nc   R2={r0:.4f}  AIC={a0:.1f}")

    # M1 quadratic
    (qa,qb,qc), y1, r1 = fit_quad(logN, logpi)
    a1,_ = aic(logpi, y1, 3)
    print(f"M1 quadratic : c={qc:.5f} ({'concave' if qc<0 else 'convex'})  R2={r1:.4f}  AIC={a1:.1f}  dAIC(vs M0)={a1-a0:.1f}")

    # M2 Michaelis-Menten (pure saturation) in linear pi
    try:
        popt,_ = curve_fit(mm, Nc, pi, p0=[0.01, 1e6], maxfev=20000)
        Vmax,K = popt
        y2 = np.log10(mm(Nc,*popt))
        a2,_ = aic(logpi, y2, 2)
        r2mm = 1 - np.sum((logpi-y2)**2)/np.sum((logpi-logpi.mean())**2)
        print(f"M2 MM (X=1)  : Vmax={Vmax:.4g}  K={K:.3g} (log10 K={math.log10(K):.2f})  R2={r2mm:.4f}  AIC={a2:.1f}  dAIC(vs M0)={a2-a0:.1f}")
    except Exception as e:
        print("M2 failed:",e); Vmax=K=0; a2=1e9

    # M4 SATURATING POWER LAW (the 2-mechanism model):
    #   pi = A * Nc^gamma / (1 + (Nc/K)^gamma)
    #   low Nc  -> pi ~ A*Nc^gamma      (slope gamma; sublinear => Ne/Nc reduction, Buffalo)
    #   high Nc -> pi ~ A*K^gamma         (plateau => gBGC saturation)
    # Fit in log space (robust over the 12-order Nc span).
    def sat_pow_log(logN, logA, gamma, logK):
        # log10 pi = logA + gamma*logN - log10(1 + (Nc/K)^gamma)
        Nc=10.0**logN; K=10.0**logK
        return logA + gamma*logN - np.log10(1.0 + (Nc/K)**gamma)
    try:
        popt4,_ = curve_fit(sat_pow_log, logN, logpi,
                           p0=[-4.0, 0.12, 9.0], maxfev=60000,
                           bounds=([-12,0,2],[2,2,16]))
        logA,gamma,logK = popt4
        y4 = sat_pow_log(logN,*popt4)
        a4,_ = aic(logpi, y4, 3)
        r4 = 1 - np.sum((logpi-y4)**2)/np.sum((logpi-logpi.mean())**2)
        # low-N slope = gamma; high-N slope -> 0. Report the plateau pi.
        plateau = 10**logA * (10**logK)**gamma
        print(f"M4 sat-pow  : low-N slope gamma={gamma:.3f}  log10 K={logK:.2f}  plateau pi={plateau:.4g}  R2={r4:.4f}  AIC={a4:.1f}  dAIC(vs M0)={a4-a0:.1f}")
        print(f"             : low-Nc slope (Buffalo/Ne-Nc)={gamma:.3f}; high-Nc plateau (gBGC sat) at pi~{plateau:.3g}")
    except Exception as e:
        print("M4 failed:",e); a4=1e9

    # M3 MM + linear tail
    try:
        popt3,_ = curve_fit(mm_lin, Nc, pi, p0=[1e-9, 0.01, 1e6], maxfev=40000)
        aa,VV,KK = popt3
        y3 = np.log10(mm_lin(Nc,*popt3))
        a3,_ = aic(logpi, y3, 3)
        r3 = 1 - np.sum((logpi-y3)**2)/np.sum((logpi-logpi.mean())**2)
        # fraction of pi at high N that is the saturated (plateau) part vs linear tail
        sat_frac = VV*KK/(aa*Nc + VV*KK)  # asymptotic ratio plateau/(plateau+linear) at high N
        print(f"M3 MM+linear : a={aa:.3g}  Vmax={VV:.4g}  K={KK:.3g} (log10 K={math.log10(KK):.2f})  R2={r3:.4f}  AIC={a3:.1f}  dAIC(vs M0)={a3-a0:.1f}")
        print(f"             : linear tail a={aa:.3g} (slope at high Nc) vs Vmax/K={VV/KK:.3g}; saturation dominates at high Nc: {VV*KK>aa*Nc.max()}")
    except Exception as e:
        print("M3 failed:",e)

    print(f"\n  -> best model by AIC: ", end="")
    models={"M0":a0,"M1":a1,"M2":a2,"M3":a3,"M4":a4}
    print(min(models,key=models.get))

    # ---- Tier 2: recombination sign test ----
    print(f"\n== Tier 2: recombination sign test (gBGC predicts -, BGS predicts +) ==\n")
    resid = logpi - y0   # residual from power law (M0)
    # ave_rec
    sub=[(r["ave_rec"],res) for r,res in zip(rows,resid) if r["ave_rec"] is not None]
    if len(sub)>=5:
        x=np.array([s[0] for s in sub]); yv=np.array([s[1] for s in sub])
        rr,pp=stats.pearsonr(x,yv); sr,sp=stats.spearmanr(x,yv)
        print(f"ave_rec  (n={len(sub)}): pearson r={rr:+.3f} p={pp:.3g} | spearman rho={sr:+.3f} p={sp:.3g}")
    # map_length
    sub=[(r["map_length"],res,r["logN"]) for r,res in zip(rows,resid) if r["map_length"] is not None]
    if len(sub)>=5:
        x=np.array([s[0] for s in sub]); yv=np.array([s[1] for s in sub]); z=np.array([s[2] for s in sub])
        rr,pp=stats.pearsonr(x,yv); sr,sp=stats.spearmanr(x,yv)
        print(f"map_len  (n={len(sub)}): pearson r={rr:+.3f} p={pp:.3g} | spearman rho={sr:+.3f} p={sp:.3g}")
        # partial out Nc: regress resid on logN (already done partially) then on maplen
        # partial correlation of resid & maplen controlling for logN
        def partial(x,y,z):
            bx,*_=np.linalg.lstsq(np.vstack([np.ones_like(z),z]).T,x,rcond=None); xr=x-np.vstack([np.ones_like(z),z]).T@bx
            by,*_=np.linalg.lstsq(np.vstack([np.ones_like(z),z]).T,y,rcond=None); yr=y-np.vstack([np.ones_like(z),z]).T@by
            r,p=stats.pearsonr(xr,yr); return r,p
        pr,ppr=partial(x,yv,z)
        print(f"           partial r (control logN) = {pr:+.3f} p={ppr:.3g}")
    # genome_size as proxy (larger genome -> more recombination-tract territory? also BGS proxy)
    sub=[(r["genome_size"],res) for r,res in zip(rows,resid) if r["genome_size"] is not None]
    if len(sub)>=5:
        x=np.array([s[0] for s in sub]); yv=np.array([s[1] for s in sub])
        rr,pp=stats.pearsonr(np.log10(x),yv)
        print(f"genome_sz(n={len(sub)}): pearson r(log10)={rr:+.3f} p={pp:.3g}")

if __name__=="__main__":
    main()
