"""E19: water violations at plants covered by a multi-plant federal case (FEC case naming 2+ primary plants of the same firm) and at the firm's
uncovered plants, around the case year; stacked event study with clean controls (as E14)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
pr = load_primary(); grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
F = f"{ROOT}/data/interim/fec"
cf = pd.read_csv(f"{F}/CASE_FACILITIES.csv", dtype=str, encoding="latin-1", usecols=["ACTIVITY_ID", "CASE_NUMBER", "REGISTRY_ID"])
ce = pd.read_csv(f"{F}/CASE_ENFORCEMENTS.csv", dtype=str, encoding="latin-1")
yc = [c for c in ["FISCAL_YEAR", "CASE_FILED_DATE"] if c in ce.columns]
print("CASE_ENFORCEMENTS year columns available:", yc)
ce = ce[["CASE_NUMBER"] + yc].drop_duplicates("CASE_NUMBER"); ce["fy"] = pd.to_numeric(ce["FISCAL_YEAR"], errors="coerce") if "FISCAL_YEAR" in ce else np.nan
cf = cf[cf.REGISTRY_ID.isin(set(grp.index))].merge(ce[["CASE_NUMBER", "fy"]], on="CASE_NUMBER", how="left"); cf["pnorm"] = cf.REGISTRY_ID.map(grp.pnorm)
any_case_years = cf.dropna(subset=["fy"]).groupby("pnorm").fy.apply(lambda s: set(s.astype(int)))
pc = cf.dropna(subset=["fy"]).groupby(["CASE_NUMBER", "pnorm"]).agg(plants=("REGISTRY_ID", "nunique"), fy=("fy", "first"), regs=("REGISTRY_ID", lambda s: set(s))).reset_index()
multi = pc[(pc.plants >= 2) & pc.fy.between(2013, 2022)].sort_values("fy")
ev = multi.groupby("pnorm").first()  # first multi-plant case per firm in 2013-2022
print("event firms:", len(ev), "| by year:", ev.fy.astype(int).value_counts().sort_index().to_dict())
stacks = []
for e in sorted(ev.fy.astype(int).unique()):
    tf = ev[ev.fy == e]
    clean = [f for f in pr.pnorm.unique() if not (set(range(e - 3, e + 4)) & any_case_years.get(f, set()))]
    win = pr[pr.yr.between(e - 3, e + 3)]
    T = win[win.pnorm.isin(tf.index)].copy(); covered = set().union(*tf.regs.values)
    T["role"] = np.where(T.REG.isin(covered), "covered plant", "uncovered sister plant"); C = win[win.pnorm.isin(clean)].copy(); C["role"] = "control"
    S = pd.concat([T, C]); S["stk"] = e; S["k"] = S.yr - e; stacks.append(S)
S = pd.concat(stacks, ignore_index=True); S["ps"] = S.REG.astype(str) + "_" + S.stk.astype(str); S["ys"] = S.yr.astype(str) + "_" + S.stk.astype(str)
print("roles (plant-stacks):", S.drop_duplicates("ps").role.value_counts().to_dict())
out, tests = [], []
for g in ["covered plant", "uncovered sister plant"]:
    D = S[S.role.isin([g, "control"])].copy(); D["tr"] = (D.role == g).astype(float); xs = []
    for k in [-3, -2, 0, 1, 2, 3]: D[f"k{k}"] = ((D.k == k) & (D.tr == 1)).astype(float); xs.append(f"k{k}")
    R = twoway_demean(D, ["any_e90"] + xs, "ps", "ys")
    f = sm.OLS(R.any_e90, R[xs]).fit(cov_type="cluster", cov_kwds={"groups": D.pnorm.astype("category").cat.codes})
    pre = f.wald_test("k-3 = 0, k-2 = 0", scalar=True); post = f.wald_test("k0 = 0, k1 = 0, k2 = 0, k3 = 0", scalar=True)
    for k in xs: out.append(dict(group=g, event_time=int(k[1:]), coef=f.params[k], se=f.bse[k]))
    out.append(dict(group=g, event_time=-1, coef=0.0, se=0.0))
    tests.append(dict(group=g, treated_plant_stacks=D[D.tr == 1].ps.nunique(), firms=D[D.tr == 1].pnorm.nunique(), mean_at_minus1=D[(D.tr == 1) & (D.k == -1)].any_e90.mean(),
                      pretrend_p=float(pre.pvalue), post_joint_p=float(post.pvalue), mean_post=np.mean([f.params[f"k{k}"] for k in [0, 1, 2, 3]])))
out = pd.DataFrame(out).sort_values(["group", "event_time"]); tests = pd.DataFrame(tests)
print(out.round(4).to_string(index=False)); print(tests.round(4).to_string(index=False))
out.to_csv(f"{OUT_T}/E19_event_coefficients.csv", index=False); tests.to_csv(f"{OUT_T}/E19_event_tests.csv", index=False)
