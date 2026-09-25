"""E16a: budget-constrained targeting simulation. E16b: do regulators already use cross-program information when allocating evaluations?"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm, statsmodels.formula.api as smf
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")

def build(prog):
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet")
    V, X, EV, CLS = ("any_air", "any_hpv", "n_air_eval", "air_class") if prog == "air" else ("any_rcra", "any_rcra_snc", "n_rcra_eval", "rcra_class")
    h = P[P.yr <= 2017]; o = P[P.yr >= 2018]
    H = h.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), nh=("yr", "nunique"), w_h=("any_e90", "mean"), p_h=(V, "mean"), ev_h=(EV, "mean"), cls=(CLS, "first"))
    O = o.groupby("REG").agg(no=("yr", "nunique"), y=(V, "max"), yx=(X, "max"), ev_o=(EV, "mean"))
    d = H.join(O, how="inner"); d = d[(d.nh >= 3) & (d.no >= 3)]; d = d[d.pnorm.map(d.pnorm.value_counts()) >= 2].copy()
    for src in ["w_h", "p_h"]:
        g = d.groupby("pnorm")[src]; n = d.pnorm.map(d.pnorm.value_counts()); d[f"firm_{src}"] = (g.transform("sum") - d[src]) / (n - 1)
    d = d.join(grp[["major"]]); d["cls"] = d.cls.fillna("unknown")
    return d.reset_index()

sim, alloc = [], []
for prog in ["air", "rcra"]:
    d = build(prog)
    feats = {"NPDES major only": ["major"], "own program history": ["p_h"], "own program + own water": ["p_h", "w_h"],
             "own program + own water + sister plants": ["p_h", "w_h", "firm_p_h", "firm_w_h"], "sister plants only": ["firm_p_h", "firm_w_h"]}
    firms = d.pnorm.unique(); rng = np.random.default_rng(171)
    for yc, ylab in [("y", "any violation"), ("yx", "most serious category")]:
        acc = {(k, b): [] for k in feats for b in [0.05, 0.10, 0.20]}
        for _ in range(100):
            tr = set(rng.choice(firms, len(firms) // 2, replace=False)); m = d.pnorm.isin(tr).values; te = d[~m]; ytest = te[yc].values.astype(bool)
            for k, f in feats.items():
                try:
                    fit = sm.Logit(d.loc[m, yc].astype(float), sm.add_constant(d.loc[m, f].astype(float))).fit(disp=0)
                    s = fit.predict(sm.add_constant(te[f].astype(float), has_constant="add")).values + rng.random(len(te)) * 1e-9  # random tie-break
                except Exception: continue
                order = np.argsort(-s)
                for b in [0.05, 0.10, 0.20]:
                    sel = order[: int(round(b * len(te)))]
                    acc[(k, b)].append((ytest[sel].sum() / max(ytest.sum(), 1), ytest[sel].mean(), ytest.mean()))
        for (k, b), v in acc.items():
            v = np.array(v)
            sim.append(dict(program=prog, outcome=ylab, information=k, budget_share=b, recall=v[:, 0].mean(), recall_sd=v[:, 0].std(), precision=v[:, 1].mean(), base_rate=v[:, 2].mean(), splits=len(v)))
    # E16b: allocation of evaluations 2018-2025 vs outcomes 2018-2025, on 2010-2017 information
    cl = dict(cov_type="cluster", cov_kwds={"groups": d.pnorm.astype("category").cat.codes})
    d["lev_o"] = np.log1p(d.ev_o)
    for dep, lab in [("y", "future violation (any)"), ("yx", "future violation (most serious)"), ("lev_o", "future evaluations, log(1+per year)")]:
        for spec, rhs in [("own program + own water + major + class", "p_h + w_h + major + C(cls)"), ("+ sister plants", "p_h + w_h + firm_p_h + firm_w_h + major + C(cls)"),
                          ("+ past evaluations", "p_h + w_h + firm_p_h + firm_w_h + np.log1p(ev_h) + major + C(cls)")]:
            f = smf.ols(f"{dep} ~ {rhs}", d).fit(**cl)
            r = dict(program=prog, outcome=lab, spec=spec, n=int(f.nobs), outcome_mean=d[dep].mean(), sd_w_h=d.w_h.std(), sd_p_h=d.p_h.std())
            for v in ["p_h", "w_h", "firm_p_h", "firm_w_h"]:
                if v in f.params: r[f"b_{v}"] = f.params[v]; r[f"se_{v}"] = f.bse[v]
            alloc.append(r)
sim, alloc = pd.DataFrame(sim), pd.DataFrame(alloc)
pd.set_option("display.width", 250)
print(sim.round(3).to_string(index=False)); print(alloc.round(4).to_string(index=False))
sim.to_csv(f"{OUT_T}/E16a_targeting_simulation.csv", index=False); alloc.to_csv(f"{OUT_T}/E16b_allocation_vs_risk.csv", index=False)
