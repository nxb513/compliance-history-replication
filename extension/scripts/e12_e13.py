"""E12: air releases (TRI) vs air violations, same plants. E13: stable trait or shared bad years (cross-media co-occurrence decomposition, lead-lag)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
from scipy.stats import spearmanr

# ---------------- E12 ----------------
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet")
frames = []
for y in [2010] + list(range(2017, 2025)):
    t = pd.read_csv(f"{ROOT}/data/raw/tri/TRI_{y}_US.csv", dtype=str, encoding="latin-1", usecols=["1. YEAR", "3. FRS ID", "50. UNIT OF MEASURE", "51. 5.1 - FUGITIVE AIR", "52. 5.2 - STACK AIR"])
    t.columns = ["yr", "REG", "unit", "fug", "stack"]
    for c in ["fug", "stack"]:
        t[c] = pd.to_numeric(t[c], errors="coerce").fillna(0.0); t.loc[t.unit.str.strip().str.lower() == "grams", c] *= 0.00220462
    t["airrel"] = t["fug"] + t["stack"]; frames.append(t[["yr", "REG", "airrel"]])
tri = pd.concat(frames); tri["yr"] = tri.yr.astype(int); tri = tri.groupby(["REG", "yr"]).airrel.sum().reset_index()
m = air.merge(tri, on=["REG", "yr"], how="inner"); m["lair"] = np.log1p(m.airrel)
m = m[m.REG.map(m.groupby("REG").yr.nunique()) >= 2]; m = m[m.pnorm.map(m.groupby("pnorm").REG.nunique()) >= 2].reset_index(drop=True)
print(f"E12 sample: plant-years {len(m)}, plants {m.REG.nunique()}, firms {m.pnorm.nunique()}, share with air releases > 0 {(m.airrel > 0).mean():.3f}, mean any_air {m.any_air.mean():.3f}")
rows = []
for col, lab in [("lair", "log(1 + TRI air releases)"), ("any_air", "any air violation (same plant-years)"), ("any_e90", "any water violation (same plant-years)")]:
    s = partition(m, col); r = reml_best(m, col)
    rows.append(dict(outcome=lab, nested_firm=s["pnorm"], nested_plant=s["REG"], nested_resid=s["resid"], reml_firm=r["firm"], reml_plant=r["plant"], reml_resid=r["resid"],
                     reml_method=r["method"], reml_llf_spread=r["llf_spread"]))
    print(rows[-1], flush=True)
bs = firm_bootstrap(m, lambda d: partition(d, "lair")["pnorm"] - partition(d, "any_air")["pnorm"], B=200, seed=121)
pl = m.groupby("REG").agg(rel=("airrel", "mean"), viol=("any_air", "mean"))
rho = spearmanr(pl.rel, pl.viol).correlation
topv = pl.viol >= pl.viol.quantile(0.9); topr = pl.rel >= pl.rel.quantile(0.9)
viol_any = pl[pl.viol > 0]
e12 = pd.DataFrame(rows)
e12["note"] = ""
e12.loc[0, "note"] = (f"nested firm share, releases minus violations: CI {np.round(ci(bs), 3)}; Spearman(long-run releases, long-run violation rate) {rho:.3f}; "
                      f"top-decile violators that are top-decile emitters {(topv & topr).sum() / topv.sum():.3f}; plants with any air violation reporting zero air releases {(viol_any.rel == 0).mean():.3f}; "
                      f"sample plants {m.REG.nunique()}, firms {m.pnorm.nunique()}, plant-years {len(m)}")
print(e12.round(3).to_string(index=False)); e12.to_csv(f"{OUT_T}/E12_air_releases_vs_violations.csv", index=False)

# ---------------- E13 ----------------
rows, ll = [], []
for prog, V in [("air", "any_air"), ("rcra", "any_rcra")]:
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet").copy()
    P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2].reset_index(drop=True)
    for c in ["any_e90", V]: P[c + "_c"] = P[c] - P.groupby("yr")[c].transform("mean")
    mean_w = P.groupby("REG")["any_e90_c"].transform("mean"); mean_p = P.groupby("REG")[V + "_c"].transform("mean")
    cov_tot = np.mean(P.any_e90_c * P[V + "_c"]); cov_b = np.mean((mean_w - mean_w.mean()) * (mean_p - mean_p.mean())); cov_w = np.mean((P.any_e90_c - mean_w) * (P[V + "_c"] - mean_p))
    R = twoway_demean(P, ["any_e90", V], "REG", "yr")
    rows.append(dict(program=prog, plant_years=len(P), plants=P.REG.nunique(), corr_plant_years_total=np.corrcoef(P.any_e90_c, P[V + "_c"])[0, 1],
                     corr_long_run_plant_means=np.corrcoef(P.groupby("REG").any_e90.mean(), P.groupby("REG")[V].mean())[0, 1],
                     corr_within_plant_deviations=np.corrcoef(R.any_e90, R[V])[0, 1], share_cov_between=cov_b / cov_tot, share_cov_within=cov_w / cov_tot))
    # lead-lag with plant and year FE
    P = P.sort_values(["REG", "yr"]); g = P.groupby("REG")
    for c in ["any_e90", V]: P[c + "_lead"] = np.where(g.yr.shift(-1) == P.yr + 1, g[c].shift(-1), np.nan)
    for dep, x, own in [("any_e90_lead", V, "any_e90"), (V + "_lead", "any_e90", V)]:
        for spec, xs in [("plant + year FE", [x]), ("plant + year FE + own current outcome", [x, own])]:
            Q = P.dropna(subset=[dep]).copy(); D = twoway_demean(Q, [dep] + xs, "REG", "yr")
            f = sm.OLS(D[dep], D[xs]).fit(cov_type="cluster", cov_kwds={"groups": Q.pnorm.astype("category").cat.codes})
            ll.append(dict(program=prog, outcome_t_plus_1=dep.replace("_lead", ""), regressor_t=x, spec=spec, coef=f.params[x], se=f.bse[x], n=int(f.nobs),
                           outcome_mean=Q[dep].mean()))
e13 = pd.DataFrame(rows); e13ll = pd.DataFrame(ll)
print(e13.round(3).to_string(index=False)); print(e13ll.round(4).to_string(index=False))
e13.to_csv(f"{OUT_T}/E13_crossmedia_decomposition.csv", index=False); e13ll.to_csv(f"{OUT_T}/E13_lead_lag.csv", index=False)
