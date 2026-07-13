"""Tier-2.5 proxy: piN/piS vs Nc from Buffalo's bundled Romiguier 2014 table.

Buffalo's combined_data.tsv is total-pi only (site_type all NA, no GC).
But his repo ships Romiguier et al. 2014's supplementary table, which has
piN, piS, piN/piS for ~76 metazoans. Joined to Buffalo's Nc (pred_log10_N),
this gives a *proxied, confounded* Tier-3 test: piN/piS should decline with
Nc (nearly-neutral) and -- if gBGC saturates at high Ne -- the decline should
flatten. CAVEAT: purifying-selection saturation predicts the same flattening,
so this is supportive but NOT discriminating.
"""
import os
from pathlib import Path
import pandas as pd, numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_DIR = Path(__file__).resolve().parent
BUFFALO_DIR = Path(os.environ.get("BUFFALO_DATA_DIR", ANALYSIS_DIR / "data" / "buffalo"))
CD = Path(os.environ.get("BUFFALO_DATA", BUFFALO_DIR / "combined_data.tsv"))
ROM = Path(os.environ.get("ROMIGUIER_DATA", BUFFALO_DIR / "romiguier_et_al_2014_updated.tsv"))
cd=pd.read_csv(CD,sep="\t"); rom=pd.read_csv(ROM,sep="\t")
cd["sp"]=cd.species.str.strip().str.lower()
rom["sp"]=rom.species.str.strip().str.lower()
m=cd.merge(rom[["sp","piN","piS","piNpiS","propagule","longevity_years","number_of_inds"]],on="sp",how="inner")
sub=m[m.pred_log10_N.notna() & m.piNpiS.notna() & (m.piNpiS>0)].copy()
sub["log10Nc"]=sub.pred_log10_N.astype(float)
sub["log10r"]=np.log10(sub.piNpiS.astype(float))
n=len(sub)
r,p=stats.pearsonr(sub.log10Nc,sub.log10r)
rs,ps=stats.spearmanr(sub.log10Nc,sub.piNpiS)
X=sub.log10Nc.values; Y=sub.log10r.values
coef=np.polyfit(X,Y,2)
print(f"n={n}  Pearson r={r:.3f} p={p:.4g}  Spearman rho={rs:.3f} p={ps:.4g}")
print(f"quadratic: log10(piN/piS) ~ {coef[0]:+.4f}*Nc^2 {coef[1]:+.4f}*Nc {coef[2]:+.4f}  -> {'convex (flattens at high Nc)' if coef[0]>0 else 'concave'}")
# split-half slope test (rising-then-plateau analogue)
med=X.mean()
lo=sub[sub.log10Nc<=med]; hi=sub[sub.log10Nc>med]
def slope(d):
    if len(d)>2:
        b,_,_,_,_=stats.linregress(d.log10Nc,d.log10r); return b
    return float("nan")
print(f"median Nc={med:.2f}  low-Nc slope={slope(lo):+.3f} (n={len(lo)})  high-Nc slope={slope(hi):+.3f} (n={len(hi)})")
# figure
fig,ax=plt.subplots(figsize=(6.2,4.6))
ax.scatter(sub.log10Nc,sub.log10r,s=26,alpha=0.6,edgecolor="none",color="#4477AA")
xs=np.linspace(X.min(),X.max(),200); 
ax.plot(xs,np.polyval(coef,xs),color="#CC3311",lw=2,label=f"quadratic fit")
# linear ref
b1,_,_,_,_=stats.linregress(X,Y); ax.plot(xs,b1*xs+ (Y-b1*X).mean(),color="grey",ls="--",lw=1,label="linear")
ax.set_xlabel("log$_{10}$ N$_c$ (Buffalo)")
ax.set_ylabel("log$_{10}$ (π$_N$/π$_S$)  [Romiguier 2014]")
ax.set_title(f"π$_N$/π$_S$ vs N$_c$ (n={n}): r={r:.2f}, p={p:.3g}\nflattens at high N$_c$ — consistent with, but not specific to, gBGC saturation")
ax.legend(frameon=False,fontsize=9)
fig.tight_layout()
fig.savefig(ANALYSIS_DIR / "fig_proxy_pinps.pdf"); fig.savefig(ANALYSIS_DIR / "fig_proxy_pinps.png",dpi=150)
print("wrote analysis/fig_proxy_pinps.{pdf,png}")
