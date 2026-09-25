"""B1: information value for targeting. Does a plant's own history, its siblings' history, or its firm's history predict violations in 2018-2025?"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
from scipy.stats import rankdata
import statsmodels.api as sm
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")

def auc(y, s):
    y = np.asarray(y, bool); s = np.asarray(s, float); ok = ~np.isnan(s); y, s = y[ok], s[ok]
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0: return np.nan
    r = rankdata(s); return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def top_decile(y, s):
    ok = ~np.isnan(s); y, s = np.asarray(y)[ok], np.asarray(s)[ok]
    cut = np.quantile(s, 0.9); top = s >= cut
    if top.mean() > 0.2:  # many ties at the cut: take the top 10% by random tie-break
        order = np.lexsort((np.random.default_rng(0).random(len(s)), -s)); top = np.zeros(len(s), bool); top[order[: int(round(0.1 * len(s)))]] = True
    return y[top].mean()

def build(prog):
    if prog == "water":
        P = load_primary(); V, X = "any_e90", None
    else:
        P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet"); V, X = ("any_air", "any_hpv") if prog == "air" else ("any_rcra", "any_rcra_snc")
    h = P[P.yr <= 2017]; o = P[P.yr >= 2018]
    H = h.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), nh=("yr", "nunique"), w_h=("any_e90", "mean"), p_h=(V, "mean"))
    O = o.groupby("REG").agg(no=("yr", "nunique"), y=(V, "max"), **({"yx": (X, "max")} if X else {}))
    d = H.join(O, how="inner"); d = d[(d.nh >= 3) & (d.no >= 3)]
    d = d[d.pnorm.map(d.pnorm.value_counts()) >= 2].copy()
    for src in ["w_h", "p_h"]:
        g = d.groupby("pnorm")[src]; n = d.pnorm.map(d.pnorm.value_counts())
        d[f"firm_{src}"] = (g.transform("sum") - d[src]) / (n - 1)
        gs = d.groupby(["pnorm", "st"])[src]; ns = d.groupby(["pnorm", "st"])[src].transform("size")
        d[f"same_{src}"] = np.where(ns > 1, (gs.transform("sum") - d[src]) / (ns - 1), np.nan)
        tot_other = g.transform("sum") - gs.transform("sum"); n_other = n - ns
        d[f"other_{src}"] = np.where(n_other > 0, tot_other / n_other.replace(0, np.nan), np.nan)
    d = d.join(grp[["major"]])
    return d.reset_index()

res, combo, diffs = [], [], []
for prog in ["air", "rcra", "water"]:
    d = build(prog)
    outs = [("y", {"air": "any air violation", "rcra": "any RCRA violation", "water": "any effluent violation"}[prog])]
    if prog != "water": outs.append(("yx", "any high priority air violation" if prog == "air" else "any RCRA significant noncomplier"))
    sets = [("own water history", "w_h"), ("same-state siblings' water history", "same_w_h"), ("other-state siblings' water history", "other_w_h"), ("firm (all siblings') water history", "firm_w_h")]
    if prog != "water":
        sets += [("own program history", "p_h"), ("same-state siblings' program history", "same_p_h"), ("other-state siblings' program history", "other_p_h"), ("firm (all siblings') program history", "firm_p_h")]
    for yc, ylab in outs:
        y = d[yc].astype(bool).values
        for lab, c in sets:
            s = d[c].values; ok = ~np.isnan(s)
            res.append(dict(program=prog, outcome=ylab, information=lab, plants=int(ok.sum()), base_rate=y[ok].mean(), auc=auc(y, s), top_decile_rate=top_decile(y[ok], s[ok]),
                            auc_own_water_same_plants=auc(y[ok], d.w_h.values[ok]), auc_own_program_same_plants=auc(y[ok], d.p_h.values[ok])))
        # combinations: logistic, fitted on half of firms, evaluated on the other half
        feats = {"baseline: NPDES major": ["major"], "own program history": ["p_h"], "own program + own water": ["p_h", "w_h"],
                 "own program + own water + firm histories": ["p_h", "w_h", "firm_p_h", "firm_w_h"], "firm histories only": ["firm_p_h", "firm_w_h"]}
        if prog == "water": feats = {"baseline: NPDES major": ["major"], "own water history": ["w_h"], "own + firm water history": ["w_h", "firm_w_h"], "firm water history only": ["firm_w_h"]}
        firms = d.pnorm.unique(); rng = np.random.default_rng(141); scores = {k: [] for k in feats}
        for _ in range(100):
            tr = set(rng.choice(firms, len(firms) // 2, replace=False)); m = d.pnorm.isin(tr).values
            for k, f in feats.items():
                Xtr = sm.add_constant(d.loc[m, f].astype(float)); Xte = sm.add_constant(d.loc[~m, f].astype(float), has_constant="add")
                try:
                    fit = sm.Logit(d.loc[m, yc].astype(float), Xtr).fit(disp=0); scores[k].append(auc(d.loc[~m, yc].astype(bool), fit.predict(Xte)))
                except Exception: pass
        for k, v in scores.items(): combo.append(dict(program=prog, outcome=ylab, model=k, mean_test_auc=np.nanmean(v), sd=np.nanstd(v), splits=len(v)))
        # firm bootstrap for AUC differences
        F = d.pnorm.astype("category").cat.codes.values; nF = F.max() + 1; idx_by_f = [np.where(F == f)[0] for f in range(nF)]
        base = "p_h" if prog != "water" else "w_h"
        pairs = [("own water minus firm water", "w_h", "firm_w_h"), ("own water minus same-state siblings' water", "w_h", "same_w_h")]
        if prog != "water": pairs += [("own program minus firm program", "p_h", "firm_p_h"), ("own program minus own water", "p_h", "w_h")]
        B = []
        for b in range(200):
            ix = np.concatenate([idx_by_f[f] for f in rng.integers(0, nF, nF)]); yy = y[ix]
            row = []
            for _, a, c in pairs:
                ok = ~np.isnan(d[c].values[ix]); row.append(auc(yy[ok], d[a].values[ix][ok]) - auc(yy[ok], d[c].values[ix][ok]))
            B.append(row)
        B = np.array(B)
        for j, (lab, a, c) in enumerate(pairs):
            ok = ~np.isnan(d[c].values)
            diffs.append(dict(program=prog, outcome=ylab, contrast=lab, estimate=auc(y[ok], d[a].values[ok]) - auc(y[ok], d[c].values[ok]), ci_low=np.nanpercentile(B[:, j], 2.5), ci_high=np.nanpercentile(B[:, j], 97.5)))
    print(prog, "plants", len(d), "firms", d.pnorm.nunique(), flush=True)
res, combo, diffs = pd.DataFrame(res), pd.DataFrame(combo), pd.DataFrame(diffs)
pd.set_option("display.width", 250)
print(res.round(3).to_string(index=False)); print(combo.round(3).to_string(index=False)); print(diffs.round(3).to_string(index=False))
res.to_csv(f"{OUT_T}/B1_auc_single.csv", index=False); combo.to_csv(f"{OUT_T}/B1_auc_combinations.csv", index=False); diffs.to_csv(f"{OUT_T}/B1_auc_differences.csv", index=False)
