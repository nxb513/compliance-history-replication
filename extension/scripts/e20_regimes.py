"""E20: results within regulatory regimes (2014 HPV Policy for air; 2015 NPDES eReporting rule for water)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
cl = lambda df: dict(cov_type="cluster", cov_kwds={"groups": df.pnorm.astype("category").cat.codes})

def plant_level(Q, V, X, CLS):
    Q = Q[Q.pnorm.map(Q.groupby("pnorm").REG.nunique()) >= 2]
    agg = dict(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), ny=("yr", "nunique"), wr=("any_e90", "mean"), ever=(V, "max"))
    if X: agg["ever_x"] = (X, "max")
    if CLS: agg["cls"] = (CLS, "first")
    pl = Q.groupby("REG").agg(**agg).join(grp[["major"]]); pl = pl[pl.ny >= 2]
    pl["c"] = (pl.wr >= 0.4).astype(int); pl["s"] = ((pl.groupby("pnorm").c.transform("sum") - pl.c) > 0).astype(int)
    pl["cls"] = pl.cls.fillna("unknown") if CLS else "none"
    return pl.reset_index()

rows = []
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet"); rcra = pd.read_parquet(f"{ROOT}/extension/rcra_panel.parquet")
for regime, (a, b), P, V, X, CLS, prog in [("1998 HPV Policy (2010-2014)", (2010, 2014), air, "any_air", "any_hpv", "air_class", "air"),
                                           ("2014 HPV Policy (2015-2025)", (2015, 2025), air, "any_air", "any_hpv", "air_class", "air"),
                                           ("before eDMR (2010-2016)", (2010, 2016), rcra, "any_rcra", "any_rcra_snc", "rcra_class", "rcra"),
                                           ("after eDMR Phase 1 (2017-2025)", (2017, 2025), rcra, "any_rcra", "any_rcra_snc", "rcra_class", "rcra"),
                                           ("before eDMR (2010-2016)", (2010, 2016), air, "any_air", "any_hpv", "air_class", "air"),
                                           ("after eDMR Phase 1 (2017-2025)", (2017, 2025), air, "any_air", "any_hpv", "air_class", "air")]:
    pl = plant_level(P[P.yr.between(a, b)], V, X, CLS)
    for y, lab in [("ever", "ever violated"), ("ever_x", "ever most serious (HPV/SNC)")]:
        f1 = smf.ols(f"{y} ~ c + s + major + C(cls)", pl).fit(**cl(pl)); f2 = smf.ols(f"{y} ~ c + major + C(cls) + C(pnorm)", pl).fit(**cl(pl))
        rows.append(dict(regime=regime, program=prog, outcome=lab, plants=len(pl), chronic_plants=int(pl.c.sum()), mean=pl[y].mean(), own=f1.params.c, own_se=f1.bse.c,
                         sibling=f1.params.s, sib_se=f1.bse.s, own_firmFE=f2.params.c, own_firmFE_se=f2.bse.c))
        print(rows[-1], flush=True)
lpm = pd.DataFrame(rows)

# water decomposition and sibling correlations in each water regime
pr = load_primary(); wrows = []
for regime, (a, b) in [("before eDMR (2010-2016)", (2010, 2016)), ("after eDMR Phase 1 (2017-2025)", (2017, 2025))]:
    Q = pr[pr.yr.between(a, b)]; Q = Q[Q.pnorm.map(Q.groupby("pnorm").REG.nunique()) >= 2].reset_index(drop=True)
    r = reml_best(Q, "any_e90"); s = partition(Q, "any_e90")
    pl = Q.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), w=("any_e90", "mean")); pl["wc"] = pl.w - pl.w.mean()
    A, B, same = [], [], []
    for _, g in pl.groupby("pnorm"):
        if len(g) < 2: continue
        x = g.wc.values; st = g.st.values; i, j = np.triu_indices(len(g), 1); A.append(x[i]); B.append(x[j]); same.append(st[i] == st[j])
    A, B, same = map(np.concatenate, (A, B, same))
    pc = lambda a, b: np.corrcoef(np.concatenate([a, b]), np.concatenate([b, a]))[0, 1]
    wrows.append(dict(regime=regime, plant_years=len(Q), plants=Q.REG.nunique(), firms=Q.pnorm.nunique(), reml_firm=r["firm"], reml_plant=r["plant"], reml_resid=r["resid"],
                      reml_method=r["method"], nested_firm=s["pnorm"], nested_plant=s["REG"], sib_corr_same_state=pc(A[same], B[same]), sib_corr_diff_state=pc(A[~same], B[~same])))
    print(wrows[-1], flush=True)
w = pd.DataFrame(wrows)
pd.set_option("display.width", 250); print(lpm.round(4).to_string(index=False)); print(w.round(3).to_string(index=False))
lpm.to_csv(f"{OUT_T}/E20_regimes_lpm.csv", index=False); w.to_csv(f"{OUT_T}/E20_regimes_water.csv", index=False)
