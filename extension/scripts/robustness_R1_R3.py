"""Robustness checks added during the audit: R1 (bootstrap CIs for E8), R2 (chronic thresholds for E9), R3 (linkage coverage)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")

def wcorr(x, y, w):
    w = w / w.sum(); mx, my = (w * x).sum(), (w * y).sum()
    return (w * (x - mx) * (y - my)).sum() / np.sqrt((w * (x - mx) ** 2).sum() * (w * (y - my) ** 2).sum())

# ---------------- R1 ----------------
r1 = []
for prog, V in [("air", "any_air"), ("rcra", "any_rcra")]:
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet"); P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2]
    pl = P.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), w=("any_e90", "mean"), a=(V, "mean")).reset_index()
    firms = pl.pnorm.astype("category"); pl["fid"] = firms.cat.codes; F = pl.fid.max() + 1
    wv, av, sv, fv = pl.w.values, pl.a.values, pl.st.values, pl.fid.values
    I, J = [], []
    for _, g in pl.groupby("fid"):
        idx = g.index.values; ii, jj = np.meshgrid(idx, idx); m = ii != jj; I.append(ii[m]); J.append(jj[m])
    I = np.concatenate(I); J = np.concatenate(J); same = sv[I] == sv[J]
    rng = np.random.default_rng(131); byst = pd.Series(np.arange(len(pl))).groupby(sv).apply(np.array); byst = byst[byst.apply(len) >= 2]
    wts = byst.apply(lambda x: len(x) * (len(x) - 1)).values.astype(float); wts /= wts.sum()
    UI, UJ = [], []
    for k in rng.choice(len(byst), 200000, p=wts):
        i, j = rng.choice(byst.iloc[k], 2, replace=False)
        if fv[i] != fv[j]: UI.append(i); UJ.append(j)
    UI, UJ = np.array(UI), np.array(UJ)
    def stats(fw):
        pw = fw[fv]
        s_plant = wcorr(wv, av, pw)
        s_same = wcorr(wv[I][same], av[J][same], (pw[I] * pw[J])[same] if False else fw[fv[I]][same])
        s_diff = wcorr(wv[I][~same], av[J][~same], fw[fv[I]][~same])
        s_unrel = wcorr(wv[UI], av[UJ], fw[fv[UI]] * fw[fv[UJ]])
        a_same = wcorr(av[I][same], av[J][same], fw[fv[I]][same]); a_diff = wcorr(av[I][~same], av[J][~same], fw[fv[I]][~same]); a_unrel = wcorr(av[UI], av[UJ], fw[fv[UI]] * fw[fv[UJ]])
        return np.array([s_plant, s_same, s_diff, s_unrel, s_plant - s_same, s_diff - s_unrel, a_same, a_diff, a_unrel, a_diff - a_unrel])
    est = stats(np.ones(F))
    B = np.array([stats(np.bincount(rng.integers(0, F, F), minlength=F).astype(float)) for _ in range(200)])
    names = ["water x program, same plant", "water x program, sibling same state", "water x program, sibling other state", "water x program, unrelated same state",
             "contrast: same plant minus sibling same state", "contrast: sibling other state minus unrelated same state",
             "program x program, sibling same state", "program x program, sibling other state", "program x program, unrelated same state",
             "contrast: program sibling other state minus unrelated"]
    for n, e, lo, hi in zip(names, est, np.nanpercentile(B, 2.5, axis=0), np.nanpercentile(B, 97.5, axis=0)):
        r1.append(dict(program=prog, statistic=n, estimate=e, ci_low=lo, ci_high=hi))
r1 = pd.DataFrame(r1); print(r1.round(3).to_string(index=False)); r1.to_csv(f"{OUT_T}/R1_E8_bootstrap_ci.csv", index=False)

# ---------------- R2 ----------------
r2 = []
q90 = grp.water_rate.quantile(0.9)
defs = {"at least 4 violation years": grp.yrs_bad >= 4, "at least 6 years (main)": grp.yrs_bad >= 6, "at least 8 years": grp.yrs_bad >= 8, "top decile of water violation rate": grp.water_rate >= q90}
for prog, V, EX, CLS in [("air", "any_air", "any_hpv", "air_class"), ("rcra", "any_rcra", "any_rcra_snc", "rcra_class")]:
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet"); P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2]
    pl = P.groupby("REG").agg(pnorm=("pnorm", "first"), ever=(V, "max"), ever_x=(EX, "max"), cls=(CLS, "first")).join(grp[["major"]]).reset_index(); pl["cls"] = pl.cls.fillna("unknown")
    for dn, flag in defs.items():
        c = pl.REG.map(flag).astype(int); pl["c"] = c
        n_ch = pl.groupby("pnorm").c.transform("sum"); pl["s"] = ((n_ch - pl.c) > 0).astype(int)
        cl = dict(cov_type="cluster", cov_kwds={"groups": pl.pnorm.astype("category").cat.codes})
        for y in ["ever", "ever_x"]:
            f1 = smf.ols(f"{y} ~ c + s + major + C(cls)", pl).fit(**cl); f2 = smf.ols(f"{y} ~ c + major + C(cls) + C(pnorm)", pl).fit(**cl)
            r2.append(dict(program=prog, outcome={"ever": "ever violated", "ever_x": "ever HPV" if prog == "air" else "ever SNC"}[y], definition=dn, chronic_plants=int(pl.c.sum()),
                           own=f1.params.c, own_se=f1.bse.c, sibling=f1.params.s, sib_se=f1.bse.s, own_firmFE=f2.params.c, own_firmFE_se=f2.bse.c))
r2 = pd.DataFrame(r2); print(r2.round(3).to_string(index=False)); r2.to_csv(f"{OUT_T}/R2_E9_chronic_thresholds.csv", index=False)

# ---------------- R3 ----------------
links = pd.read_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
has = {k: set(links[links.PGM_SYS_ACRNM == a].REGISTRY_ID) for k, a in [("air", "AIR"), ("rcra", "RCRAINFO")]}
g = grp.copy(); g.index.name = "REG"; g = g.reset_index()
r3 = []
for k in ["air", "rcra"]:
    g["linked"] = g.REG.isin(has[k])
    for lk, d in g.groupby("linked"):
        r3.append(dict(program=k, linked=lk, plants=len(d), water_violation_rate=d.water_rate.mean(), chronic_share=d.chronic.mean(), major_share=d.major.mean()))
    f = smf.ols("water_rate ~ linked + C(pnorm)", g.assign(linked=g.linked.astype(int))).fit(cov_type="cluster", cov_kwds={"groups": g.pnorm.astype("category").cat.codes})
    r3.append(dict(program=k, linked="difference within firm", water_violation_rate=f.params.linked, chronic_share=f.bse.linked))
r3 = pd.DataFrame(r3); print(r3.round(3).to_string(index=False)); r3.to_csv(f"{OUT_T}/R3_linkage_coverage.csv", index=False)
