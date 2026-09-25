"""E7-E10 for one program panel (air or rcra). Usage: python cross_analysis.py air|rcra"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
prog = sys.argv[1]
P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet")
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
V = "any_air" if prog == "air" else "any_rcra"
EXTRA = ["any_hpv"] if prog == "air" else ["any_rcra_snc"] if "any_rcra_snc" in P else []
EVAL = "n_air_eval" if prog == "air" else "n_rcra_eval"
samples = {"all plants with a program ID": P}
if prog == "air":
    samples["major and synthetic-minor air sources"] = P[P.air_class.isin(["major", "synthetic minor"])]
else:
    samples["large quantity generators or operating TSDFs"] = P[P.rcra_class == "LQG or TSDF"]
CLS = "air_class" if prog == "air" else "rcra_class"
out7, out8, out9 = [], [], []
for sname, S in samples.items():
    S = S.copy(); nf = S.groupby("pnorm").REG.nunique(); S = S[S.pnorm.map(nf) >= 2].reset_index(drop=True)
    print(f"\n##### {prog.upper()} | {sname}: plants {S.REG.nunique()}, firms {S.pnorm.nunique()}, plant-years {len(S)}, mean {V} {S[V].mean():.3f}")
    # ---- E7: partition, program vs water on the same plant-years
    for col in [V] + EXTRA + ["any_e90"]:
        s = partition(S, col)
        out7.append(dict(program=prog, sample=sname, outcome=col, firm=s["pnorm"], plant=s["REG"], resid=s["resid"], D=s["REG"] - s["pnorm"],
                         plants=S.REG.nunique(), firms=S.pnorm.nunique(), plant_years=len(S), mean=S[col].mean()))
    def st(d):
        a = partition(d, V); w = partition(d, "any_e90")
        return (a["REG"] - a["pnorm"]) - (w["REG"] - w["pnorm"]), a["REG"] - a["pnorm"], a["pnorm"] - w["pnorm"]
    bs = firm_bootstrap(S, st, B=200, seed=71)
    out7.append(dict(program=prog, sample=sname, outcome=f"D({V}) minus D(any_e90)", D=np.nan,
                     note=f"diff CI {np.round(ci([b[0] for b in bs]),3)}; D_{V} CI {np.round(ci([b[1] for b in bs]),3)}; firm share diff CI {np.round(ci([b[2] for b in bs]),3)}"))
    try:
        r = reml_shares(S, V); out7.append(dict(program=prog, sample=sname, outcome=f"{V}, REML", firm=r["firm"], plant=r["plant"], resid=r["resid"], D=r["plant"] - r["firm"]))
    except Exception as e:
        print("REML failed:", e)
    # ---- E8: cross-program covariance, plant vs sibling vs unrelated same-state
    pl = S.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), w=("any_e90", "mean"), a=(V, "mean"))
    pl["wc"] = pl.w - pl.w.mean(); pl["ac"] = pl.a - pl.a.mean()
    within = np.corrcoef(pl.wc, pl.ac)[0, 1]
    recs = []
    for f, g in pl.groupby("pnorm"):
        if len(g) < 2: continue
        wv, av, sv = g.wc.values, g.ac.values, g.st.values
        i, j = np.where(~np.eye(len(g), dtype=bool))
        recs.append(pd.DataFrame(dict(pnorm=f, w=wv[i], a=av[j], aa_i=av[i], same=sv[i] == sv[j])))
    sp = pd.concat(recs, ignore_index=True)
    rng = np.random.default_rng(72); stv = pl.st.values; fm = pl.pnorm.values; wv = pl.wc.values; av = pl.ac.values
    groups = pd.Series(np.arange(len(pl))).groupby(stv).apply(np.array); groups = groups[groups.apply(len) >= 2]
    wts = groups.apply(lambda x: len(x) * (len(x) - 1)).values.astype(float); wts /= wts.sum()
    A, B, AA = [], [], []
    for k in rng.choice(len(groups), 300000, p=wts):
        i, j = rng.choice(groups.iloc[k], 2, replace=False)
        if fm[i] != fm[j]: A.append(wv[i]); B.append(av[j]); AA.append(av[i])
    A, B, AA = map(np.array, (A, B, AA))
    c = lambda x, y: np.corrcoef(x, y)[0, 1]
    out8.append(dict(program=prog, sample=sname, water_x_program_same_plant=within,
                     water_x_program_sibling_same_state=c(sp[sp.same].w, sp[sp.same].a), water_x_program_sibling_diff_state=c(sp[~sp.same].w, sp[~sp.same].a),
                     water_x_program_unrelated_same_state=c(A, B),
                     program_x_program_sibling_same_state=c(sp[sp.same].aa_i, sp[sp.same].a), program_x_program_sibling_diff_state=c(sp[~sp.same].aa_i, sp[~sp.same].a),
                     program_x_program_unrelated_same_state=c(AA, B), n_plants=len(pl)))
    # ---- E9 / E10: groups defined on the WATER panel
    g2 = pl.join(grp[["chronic", "sib_chronic", "major"]])
    ever = S.groupby("REG").agg(ever_prog=(V, "max"), rate_prog=(V, "mean"), evals=(EVAL, "mean"), **({"ever_extra": (EXTRA[0], "max")} if EXTRA else {}))
    g2 = g2.join(ever)
    g2["group"] = np.select([g2.chronic, ~g2.chronic & g2.sib_chronic], ["chronic water violator", "non-chronic sibling"], "firm with no chronic plant")
    agg = dict(plants=("ever_prog", "size"), ever_violation=("ever_prog", "mean"), annual_rate=("rate_prog", "mean"), evals_per_year=("evals", "mean"))
    if EXTRA: agg["ever_" + EXTRA[0]] = ("ever_extra", "mean")
    tab = g2.groupby("group").agg(**agg).reset_index(); tab.insert(0, "sample", sname); tab.insert(0, "program", prog)
    print(tab.round(3).to_string(index=False)); out9.append(tab)
    g2 = g2.reset_index(); g2["chronic_i"] = g2.chronic.astype(int); g2["sib_i"] = g2.sib_chronic.astype(int); g2["y"] = g2.ever_prog.astype(float)
    cl = dict(cov_type="cluster", cov_kwds={"groups": g2.pnorm.astype("category").cat.codes})
    g2[CLS] = g2.REG.map(S.drop_duplicates("REG").set_index("REG")[CLS]).fillna("unknown")
    ctrl = " + major" + (f" + C({CLS})" if g2[CLS].nunique() > 1 else "")
    f1 = smf.ols("y ~ chronic_i + sib_i" + ctrl, g2).fit(**cl); f2 = smf.ols("y ~ chronic_i" + ctrl + " + C(pnorm)", g2).fit(**cl)
    lp = pd.DataFrame([dict(program=prog, sample=sname, own_chronic=f1.params["chronic_i"], own_se=f1.bse["chronic_i"], sibling_chronic=f1.params["sib_i"],
                            sib_se=f1.bse["sib_i"], own_chronic_firmFE=f2.params["chronic_i"], own_firmFE_se=f2.bse["chronic_i"], n=int(f1.nobs))])
    print(lp.round(4).to_string(index=False)); out9.append(lp)
pd.DataFrame(out7).to_csv(f"{OUT_T}/E7_{prog}_partition.csv", index=False)
pd.DataFrame(out8).to_csv(f"{OUT_T}/E8_{prog}_crossprogram_corr.csv", index=False)
pd.concat([x for x in out9 if "group" in x], ignore_index=True).to_csv(f"{OUT_T}/E9_E10_{prog}_by_group.csv", index=False)
pd.concat([x for x in out9 if "own_chronic" in x], ignore_index=True).to_csv(f"{OUT_T}/E9_{prog}_lpm.csv", index=False)
print("\nE7:"); print(pd.DataFrame(out7).drop(columns=[c for c in ["plants", "firms", "plant_years"] if c in pd.DataFrame(out7)]).round(3).to_string(index=False))
print("\nE8:"); print(pd.DataFrame(out8).round(3).T.to_string())
