#!/usr/bin/env python3
"""Figure: extend Buffalo — pi vs Nc with the rising-then-plateau decomposition.

Left: scatter of log10 pi vs log10 Nc (n=172) with M0 power law (dashed) and
       M4 saturating power law (solid); two-half slopes annotated.
Right: residuals from M0 vs recombination (map_length), showing the negative
       (gBGC) direction controlling for Nc.
"""
import csv, math, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

ANALYSIS_DIR = Path(__file__).resolve().parent
DATA = Path(os.environ.get("BUFFALO_DATA", ANALYSIS_DIR / "data" / "buffalo" / "combined_data.tsv"))
OUT = Path(os.environ.get("BUFFALO_FIGURE", ANALYSIS_DIR / "fig_extend_buffalo.pdf"))

rows=[]
for r in csv.DictReader(open(DATA),delimiter="\t"):
    try: pi=float(r["diversity"]); N=float(r["pred_log10_N"])
    except: continue
    if pi>0 and N>0:
        rows.append(dict(logN=N, logpi=math.log10(pi),
                         ml=(float(r["map_length"]) if r["map_length"] not in ("","NA","NaN") else None)))
logN=np.array([r["logN"] for r in rows]); logpi=np.array([r["logpi"] for r in rows])

# M0 power law
A=np.vstack([np.ones_like(logN),logN]).T
b0,*_=np.linalg.lstsq(A,logpi,rcond=None); y0=A@b0
# M4 saturating power law
def sp(lN,lA,g,lK): return lA+g*lN-np.log10(1+(10**lN/10**lK)**g)
p4,_=curve_fit(sp,logN,logpi,p0=[-4,0.12,9],maxfev=60000,bounds=([-12,0,2],[2,2,16]))
# two-half slopes
o=np.argsort(logN); x=logN[o]; y=logpi[o]; h=len(x)//2
sl=stats.linregress(x[:h],y[:h]); sh=stats.linregress(x[h:],y[h:])

fig,(axL,axR)=plt.subplots(1,2,figsize=(11,4.6))
axL.scatter(logN,logpi,s=16,alpha=0.5,edgecolor="none",color="#444")
xs=np.linspace(logN.min(),logN.max(),200)
axL.plot(xs,b0[0]+b0[1]*xs,"--",color="#c0392b",lw=1.6,label=f"M0 power law (slope {b0[1]:.3f}, $R^2$=0.264)")
axL.plot(xs,sp(xs,*p4),"-",color="#2c3e50",lw=1.8,label=f"M4 saturating (low-N $\\gamma$={p4[1]:.2f}, plateau)")
axL.axvline(x[h],color="#888",ls=":",lw=1)
axL.text(x[h]+0.1,axL.get_ylim()[0]+0.15,
         f"low-Nc half slope {sl.slope:+.3f} (p={sl.pvalue:.3f})\nhigh-Nc half slope {sh.slope:+.3f} (p={sh.pvalue:.2f})",
         fontsize=8,color="#555")
axL.set_xlabel("log$_{10}$ census population size $N_c$")
axL.set_ylabel("log$_{10}$ nucleotide diversity $\\pi$")
axL.set_title("Buffalo's 172 metazoans: rising then plateau")
axL.legend(fontsize=8,loc="upper left")

# residual vs map_length, partial out Nc
resid=logpi-y0
sub=[(r["ml"],r["logN"]) for r in rows if r["ml"] is not None]
ml=np.array([s[0] for s in sub]); zN=np.array([s[1] for s in sub]); rsub=resid[[i for i,r in enumerate(rows) if r["ml"] is not None]]
# partial: residualize ml and rsub against logN
Z=np.vstack([np.ones_like(zN),zN]).T
bml,*_=np.linalg.lstsq(Z,ml,rcond=None); mlr=ml-Z@bml
brs,*_=np.linalg.lstsq(Z,rsub,rcond=None); rsr=rsub-Z@brs
pr,pp=stats.pearsonr(mlr,rsr)
axR.scatter(mlr,rsr,s=20,alpha=0.6,edgecolor="none",color="#2c3e50")
m,c=np.polyfit(mlr,rsr,1); axR.plot(mlr,m*mlr+c,"-",color="#c0392b",lw=1.4)
axR.axhline(0,color="#888",lw=0.8)
axR.set_xlabel("map length (residualized for $N_c$)")
axR.set_ylabel("$\\pi$ residual (from power law)")
axR.set_title(f"Recombination sign test (n={len(ml)}): gBGC predicts $-$, BGS predicts $+$\npartial r = {pr:+.2f} (p={pp:.2f})")
plt.tight_layout()
plt.savefig(OUT)
plt.savefig(OUT.with_suffix(".png"),dpi=150)
print("wrote",OUT,"and png")
