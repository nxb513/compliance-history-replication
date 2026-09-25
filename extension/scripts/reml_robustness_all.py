"""Robustness of every REML variance share in the extension to the optimizer.
Each model is re-fitted with lbfgs, powell and Nelder-Mead (maxiter 2000); the fit with the highest REML log-likelihood is the reference.
Rebuilds each REML input exactly as in the original scripts (e1_e2, e2b_followup, e3_e6, e4_analyze, e5b_followup, cross_analysis)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
only = sys.argv[1] if len(sys.argv) > 1 else None

def fit_all(d, col, three=False):
    d = d.copy(); d["pid"] = d.REG.astype(str)
    vcf = {"fs": "0 + C(fs)", "plant": "0 + C(pid)"} if three else {"plant": "0 + C(pid)"}
    res = []
    for meth in ["lbfgs", "powell", "nm"]:
        try:
            md = smf.mixedlm(f"{col} ~ C(yr)", d, groups=d["pnorm"], re_formula="1", vc_formula=vcf).fit(method=meth, maxiter=2000, reml=True)
            vc = dict(zip(md.model.exog_vc.names, md.vcomp)); vf = float(md.cov_re.iloc[0, 0]); ve = float(md.scale); t = vf + sum(vc.values()) + ve
            res.append(dict(method=meth, llf=md.llf, converged=bool(md.converged), firm=vf / t, firm_within_state=vc.get("fs", np.nan) / t if three else np.nan,
                            plant=vc["plant"] / t, resid=ve / t))
        except Exception as e:
            res.append(dict(method=meth, note=str(e)[:80]))
    return res

jobs = []
pr = load_primary()
jobs.append(("E1/fixed core", "log(1+effluent violations), primary sample", pr, "le90", False))
jobs.append(("E1", "log(1+reporting violations), primary sample", pr, "lrep", False))
pr3 = pr.copy(); pr3["fs"] = pr3.pnorm.astype(str) + "|" + pr3.STATE_CODE.astype(str)
jobs.append(("E2", "nested firm/firm-state/plant, log(1+effluent)", pr3, "le90", True))
jobs.append(("E2", "nested firm/firm-state/plant, any effluent", pr3, "any_e90", True))
# E3 input
frames = []
for y in [2010] + list(range(2017, 2025)):
    t = pd.read_csv(f"{ROOT}/data/raw/tri/TRI_{y}_US.csv", dtype=str, encoding="latin-1", usecols=["1. YEAR", "3. FRS ID", "50. UNIT OF MEASURE", "53. 5.3 - WATER"])
    t.columns = ["yr", "REG", "unit", "water"]; t["water"] = pd.to_numeric(t.water, errors="coerce").fillna(0.0)
    t.loc[t.unit.str.strip().str.lower() == "grams", "water"] *= 0.00220462; frames.append(t)
tri = pd.concat(frames); tri["yr"] = tri.yr.astype(int); tri = tri.groupby(["REG", "yr"])[["water"]].sum().reset_index()
m = pr.merge(tri, on=["REG", "yr"], how="inner"); m["lwater"] = np.log1p(m.water)
m = m[m.REG.map(m.groupby("REG").yr.nunique()) >= 2]; m = m[m.pnorm.map(m.groupby("pnorm").REG.nunique()) >= 2]
jobs.append(("E3", "log(1+TRI water releases), TRI plants", m, "lwater", False))
jobs.append(("E3", "log(1+effluent violations), TRI plants", m, "le90", False))
# E4d and E5b inputs
v = pd.read_parquet(f"{ROOT}/extension/e90_param_primary.parquet"); v["yr"] = v.MONITORING_PERIOD_END_DATE.str[-4:].astype(int)
v["exc"] = pd.to_numeric(v.EXCEEDENCE_PCT, errors="coerce")
v4 = v.merge(pr.drop_duplicates("REG")[["REG", "pnorm"]], left_on="REGISTRY_ID", right_on="REG", how="inner")
sv = v4.dropna(subset=["exc"]).groupby(["REG", "yr"]).agg(med_exc=("exc", "median")).reset_index().merge(pr[["REG", "yr", "pnorm"]], on=["REG", "yr"])
sv["lmed"] = np.log1p(sv.med_exc); sv = sv[sv.pnorm.map(sv.groupby("pnorm").REG.nunique()) >= 2]
jobs.append(("E4d", "log(1+median exceedance %), violating plant-years", sv, "lmed", False))
pl = pd.read_parquet(f"{ROOT}/extension/plant_attributes_primary.parquet")
npv = v.groupby(["REGISTRY_ID", "yr"]).PARAMETER_CODE.nunique().rename("n_params_violated").reset_index().rename(columns={"REGISTRY_ID": "REG"})
d5 = pr.merge(npv, on=["REG", "yr"], how="left"); d5["n_params_violated"] = d5.n_params_violated.fillna(0); d5["n_lim"] = d5.REG.map(pl.n_limited_params)
d5 = d5[d5.n_lim >= 1].copy(); d5["rate"] = (d5.n_params_violated / d5.n_lim).clip(upper=1); d5 = d5[d5.pnorm.map(d5.groupby("pnorm").REG.nunique()) >= 3]
jobs.append(("E5b", "share of limited parameters violated", d5, "rate", False))
# E7 inputs
for prog, V, CLS, keep in [("air", "any_air", "air_class", ["major", "synthetic minor"]), ("rcra", "any_rcra", "rcra_class", ["LQG or TSDF"])]:
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet")
    for sname, S in [("all plants with a program ID", P), ("covered subset", P[P[CLS].isin(keep)])]:
        S = S[S.pnorm.map(S.groupby("pnorm").REG.nunique()) >= 2].reset_index(drop=True)
        jobs.append(("E7", f"{prog}: {V}, {sname}", S, V, False)); jobs.append(("E7", f"{prog}: any_e90 (water), {sname}", S, "any_e90", False))

out = []
for tag, lab, d, col, three in jobs:
    if only and only not in tag: continue
    cols = [col, "pnorm", "REG", "yr"] + (["fs"] if three else [])
    res = fit_all(d[cols].dropna(), col, three)
    ok = [r for r in res if "llf" in r]; best = max(ok, key=lambda r: r["llf"])
    for r in res:
        out.append(dict(analysis=tag, model=lab, n=len(d), **r, best=(r is best), llf_gap_to_best=(best["llf"] - r["llf"]) if "llf" in r else np.nan))
    print(tag, lab, "| best:", best["method"], {k: round(best[k], 3) for k in ["firm", "firm_within_state", "plant", "resid"] if pd.notna(best[k])},
          "| gaps:", {r["method"]: round(best["llf"] - r["llf"], 2) for r in ok}, flush=True)
    pd.DataFrame(out).to_csv(f"{OUT_T}/REML_optimizer_robustness{'_' + only.replace('/', '') if only else ''}.csv", index=False)
