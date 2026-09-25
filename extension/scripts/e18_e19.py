"""E18: air stack-test failures (compliance measured independently of inspector discretion). E19: feasibility of firm-wide settlements (FEC cases naming 2+ plants of a firm)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet")
fac = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FACILITIES.csv", dtype=str, usecols=["PGM_SYS_ID", "REGISTRY_ID"])
stt = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_STACK_TESTS.csv", dtype=str, usecols=["PGM_SYS_ID", "ACTUAL_END_DATE", "AIR_STACK_TEST_STATUS_CODE"]).merge(fac, on="PGM_SYS_ID")
stt = stt[stt.AIR_STACK_TEST_STATUS_CODE.isin(["PSS", "FAI"])]; stt["yr"] = pd.to_numeric(stt.ACTUAL_END_DATE.str[-4:], errors="coerce")
st = stt.groupby(["REGISTRY_ID", "yr"]).agg(n_tests=("AIR_STACK_TEST_STATUS_CODE", "size"), n_fail=("AIR_STACK_TEST_STATUS_CODE", lambda s: (s == "FAI").sum())).reset_index().rename(columns={"REGISTRY_ID": "REG"})
P = air.merge(st, on=["REG", "yr"], how="inner"); P["any_fail"] = (P.n_fail > 0).astype(float)
P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2].reset_index(drop=True)
print(f"E18 plant-years with a completed stack test: {len(P)}, plants {P.REG.nunique()}, firms {P.pnorm.nunique()}, fail rate {P.any_fail.mean():.3f}")
rows = []
for col in ["any_fail", "any_e90"]:
    s = partition(P, col); r = reml_best(P, col)
    rows.append(dict(outcome=col, nested_firm=s["pnorm"], nested_plant=s["REG"], nested_resid=s["resid"], reml_firm=r["firm"], reml_plant=r["plant"], reml_resid=r["resid"], reml_method=r["method"]))
    print(rows[-1], flush=True)
pl = P.groupby("REG").agg(pnorm=("pnorm", "first"), ever_fail=("any_fail", "max"), n_tests=("n_tests", "sum"), cls=("air_class", "first")).join(grp[["chronic", "sib_chronic", "major"]])
pl["c"] = pl.chronic.astype(int); pl["s"] = pl.sib_chronic.astype(int); pl["ltests"] = np.log1p(pl.n_tests); pl["cls"] = pl.cls.fillna("unknown"); pl = pl.reset_index()
cl = dict(cov_type="cluster", cov_kwds={"groups": pl.pnorm.astype("category").cat.codes})
f1 = smf.ols("ever_fail ~ c + s + major + C(cls) + ltests", pl).fit(**cl); f2 = smf.ols("ever_fail ~ c + major + C(cls) + ltests + C(pnorm)", pl).fit(**cl)
g = pl.assign(group=np.select([pl.chronic, ~pl.chronic & pl.sib_chronic], ["chronic water violator", "non-chronic sibling"], "firm with no chronic plant")).groupby("group").agg(plants=("REG", "size"), ever_fail=("ever_fail", "mean"), tests=("n_tests", "mean"))
print(g.round(3).to_string())
lp = dict(own=f1.params.c, own_se=f1.bse.c, sibling=f1.params.s, sib_se=f1.bse.s, own_firmFE=f2.params.c, own_firmFE_se=f2.bse.c, n=int(f1.nobs), mean=pl.ever_fail.mean())
print(lp)
pd.DataFrame(rows).to_csv(f"{OUT_T}/E18_stacktest_partition.csv", index=False); g.to_csv(f"{OUT_T}/E18_stacktest_by_group.csv"); pd.DataFrame([lp]).to_csv(f"{OUT_T}/E18_stacktest_lpm.csv", index=False)

# ---------------- E19 feasibility ----------------
F = f"{ROOT}/data/interim/fec"
cf = pd.read_csv(f"{F}/CASE_FACILITIES.csv", dtype=str, encoding="latin-1")
print("CASE_FACILITIES columns:", list(cf.columns)[:20])
regcol = [c for c in cf.columns if "REGISTRY" in c.upper()][0]; casecol = [c for c in cf.columns if c.upper() in ("CASE_NUMBER", "ACTIVITY_ID")][0]
cf = cf[cf[regcol].isin(set(grp.index))].copy(); cf["pnorm"] = cf[regcol].map(grp.pnorm)
ce = pd.read_csv(f"{F}/CASE_ENFORCEMENTS.csv", dtype=str, encoding="latin-1")
ycol = "FISCAL_YEAR" if "FISCAL_YEAR" in ce.columns else None
if ycol and casecol in ce.columns: cf = cf.merge(ce[[casecol, ycol]].drop_duplicates(casecol), on=casecol, how="left"); cf["fy"] = pd.to_numeric(cf[ycol], errors="coerce")
else: cf["fy"] = np.nan
per_case = cf.groupby([casecol, "pnorm"]).agg(plants=(regcol, "nunique"), fy=("fy", "first")).reset_index()
multi = per_case[(per_case.plants >= 2) & (per_case.fy.between(2010, 2025) | per_case.fy.isna())]
print(f"E19: FEC cases naming >=2 primary plants of the same firm (FY2010-2025 or year missing): {len(multi)}; firms {multi.pnorm.nunique()}; plants covered {int(multi.plants.sum())}")
pd.DataFrame([dict(cases_2plus_same_firm=len(multi), firms=multi.pnorm.nunique(), plants_covered=int(multi.plants.sum()), with_year=int(multi.fy.notna().sum()),
                   cases_2to5=int(multi.plants.between(2, 5).sum()), cases_6plus=int((multi.plants >= 6).sum()))]).to_csv(f"{OUT_T}/E19_firmwide_cases_feasibility.csv", index=False)
multi.to_csv(f"{OUT_T}/E19_firmwide_cases_list.csv", index=False)
