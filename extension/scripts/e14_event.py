"""E14 (pre-specified, descriptive): water violations at sister plants around the first formal NPDES action at a firm (stacked event study, clean controls)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
pr = load_primary()
ids = pd.read_parquet(f"{ROOT}/extension/data_npdes_ids_primary.parquet")
fe = pd.read_parquet(f"{ROOT}/data/interim/formal_enf.parquet", columns=["NPDES_ID", "SETTLEMENT_ENTERED_DATE", "AGENCY"])
fe["yr"] = pd.to_numeric(fe.SETTLEMENT_ENTERED_DATE.str[-4:], errors="coerce")
fe = fe.merge(ids, on="NPDES_ID").rename(columns={"REGISTRY_ID": "REG"})
fe = fe[fe.REG.isin(set(pr.REG))]
firm_of = pr.drop_duplicates("REG").set_index("REG").pnorm; st_of = pr.drop_duplicates("REG").set_index("REG").STATE_CODE
fe["pnorm"] = fe.REG.map(firm_of)
pf = fe.groupby(["pnorm", "yr"]).REG.apply(set)          # plants with a formal action, by firm-year
firm_years = fe.groupby("pnorm").yr.apply(set)
print("formal actions at primary plants:", len(fe), "| plants", fe.REG.nunique(), "| firms", fe.pnorm.nunique(), "| agency", fe.AGENCY.value_counts().to_dict())
events = {}
for f, ys in firm_years.items():
    for e in range(2013, 2023):
        if e in ys and not ({e - 3, e - 2, e - 1} & ys):
            events[f] = e; break
ev = pd.Series(events, name="e"); print("event firms:", len(ev), "| by year:", ev.value_counts().sort_index().to_dict())
stacks = []
for e in range(2013, 2023):
    tf = ev[ev == e].index
    if len(tf) == 0: continue
    clean = [f for f in pr.pnorm.unique() if not (set(range(e - 3, e + 4)) & firm_years.get(f, set()))]
    win = pr[pr.yr.between(e - 3, e + 3)]
    T = win[win.pnorm.isin(tf)].copy(); C = win[win.pnorm.isin(clean)].copy()
    enforced = set().union(*[pf.loc[(f, e)] for f in tf])
    enf_state = {f: {st_of[r] for r in pf.loc[(f, e)]} for f in tf}
    T["role"] = np.where(T.REG.isin(enforced), "enforced plant", np.where([s in enf_state[f] for s, f in zip(T.STATE_CODE, T.pnorm)], "sibling, same state", "sibling, other state"))
    C["role"] = "control"
    S = pd.concat([T, C]); S["stk"] = e; S["k"] = S.yr - e
    stacks.append(S)
S = pd.concat(stacks, ignore_index=True)
S["ps"] = S.REG.astype(str) + "_" + S.stk.astype(str); S["ys"] = S.yr.astype(str) + "_" + S.stk.astype(str)
print("stacked rows:", len(S), "| roles (plant-stacks):", S.drop_duplicates("ps").role.value_counts().to_dict())
out, tests = [], []
for grpname, roles in [("all siblings", ["sibling, same state", "sibling, other state"]), ("siblings, same state", ["sibling, same state"]),
                       ("siblings, other state", ["sibling, other state"]), ("enforced plant", ["enforced plant"])]:
    D = S[S.role.isin(roles + ["control"])].copy(); D["tr"] = D.role.isin(roles).astype(float)
    xs = []
    for k in [-3, -2, 0, 1, 2, 3]:
        D[f"k{k}"] = ((D.k == k) & (D.tr == 1)).astype(float); xs.append(f"k{k}")
    R = twoway_demean(D, ["any_e90"] + xs, "ps", "ys")
    f = sm.OLS(R.any_e90, R[xs]).fit(cov_type="cluster", cov_kwds={"groups": D.pnorm.astype("category").cat.codes})
    pre = f.wald_test("k-3 = 0, k-2 = 0", scalar=True)
    post = f.wald_test("k0 = 0, k1 = 0, k2 = 0, k3 = 0", scalar=True)
    base = D[(D.tr == 1) & (D.k == -1)].any_e90.mean()
    for k in xs: out.append(dict(group=grpname, event_time=int(k[1:]), coef=f.params[k], se=f.bse[k]))
    out.append(dict(group=grpname, event_time=-1, coef=0.0, se=0.0))
    tests.append(dict(group=grpname, treated_plant_stacks=D[D.tr == 1].ps.nunique(), control_plant_stacks=D[D.tr == 0].ps.nunique(), firms_treated=D[D.tr == 1].pnorm.nunique(),
                      mean_at_minus1=base, pretrend_p=float(pre.pvalue), post_joint_p=float(post.pvalue), mean_post=np.mean([f.params[f"k{k}"] for k in [0, 1, 2, 3]])))
out = pd.DataFrame(out).sort_values(["group", "event_time"]); tests = pd.DataFrame(tests)
print(out.round(4).to_string(index=False)); print(tests.round(4).to_string(index=False))
out.to_csv(f"{OUT_T}/E14_event_coefficients.csv", index=False); tests.to_csv(f"{OUT_T}/E14_event_tests.csv", index=False)
