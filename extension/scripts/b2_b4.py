"""B2 (detection source), B3 (heterogeneity of the cross-media gap), B4 (stability across 2010-2017 and 2018-2025)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
attr = pd.read_parquet(f"{ROOT}/extension/plant_attributes_primary.parquet")[["n3"]]
links = pd.read_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
cl = lambda df: dict(cov_type="cluster", cov_kwds={"groups": df.pnorm.astype("category").cat.codes})

def plant_frame(prog):
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet"); P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2]
    V, CLS = ("any_air", "air_class") if prog == "air" else ("any_rcra", "rcra_class")
    pl = P.groupby("REG").agg(pnorm=("pnorm", "first"), ever=(V, "max"), cls=(CLS, "first"), nplant=("nplant", "first")).join(grp[["chronic", "sib_chronic", "major"]])
    pl["c"] = pl.chronic.astype(int); pl["s"] = pl.sib_chronic.astype(int); pl["cls"] = pl.cls.fillna("unknown")
    return P, pl

def lpm(pl, y):
    f1 = smf.ols(f"{y} ~ c + s + major + C(cls)", pl).fit(**cl(pl)); f2 = smf.ols(f"{y} ~ c + major + C(cls) + C(pnorm)", pl).fit(**cl(pl))
    return dict(mean=pl[y].mean(), own=f1.params.c, own_se=f1.bse.c, sibling=f1.params.s, sib_se=f1.bse.s, own_firmFE=f2.params.c, own_firmFE_se=f2.bse.c, n=int(f1.nobs))

# ---------------- B2 ----------------
b2 = []
P, pl = plant_frame("air")
fac = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FACILITIES.csv", dtype=str, usecols=["PGM_SYS_ID", "REGISTRY_ID"]); fac = fac[fac.REGISTRY_ID.isin(set(pl.index))]
vh = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_VIOLATION_HISTORY.csv", dtype=str, usecols=["PGM_SYS_ID", "ENF_RESPONSE_POLICY_CODE", "EARLIEST_FRV_DETERM_DATE", "HPV_DAYZERO_DATE", "AGENCY_TYPE_DESC"]).merge(fac, on="PGM_SYS_ID")
vh["date"] = np.where(vh.ENF_RESPONSE_POLICY_CODE == "HPV", vh.HPV_DAYZERO_DATE.fillna(vh.EARLIEST_FRV_DETERM_DATE), vh.EARLIEST_FRV_DETERM_DATE.fillna(vh.HPV_DAYZERO_DATE))
vh["yr"] = pd.to_numeric(vh.date.str[-4:], errors="coerce"); vh = vh.rename(columns={"REGISTRY_ID": "REG"}).merge(P[["REG", "yr"]], on=["REG", "yr"])  # only panel plant-years
vh["src"] = np.where(vh.AGENCY_TYPE_DESC == "U.S. EPA", "EPA", "state or local")
for src in ["state or local", "EPA"]:
    pl[f"ev_{src[:3]}"] = pl.index.isin(set(vh[vh.src == src].REG)).astype(float)
    b2.append(dict(program="air", detected_by=src, violations=int((vh.src == src).sum()), **lpm(pl, f"ev_{src[:3]}")))
P, pl = plant_frame("rcra")
idmap = links[links.PGM_SYS_ACRNM == "RCRAINFO"][["PGM_SYS_ID", "REGISTRY_ID"]].drop_duplicates().rename(columns={"PGM_SYS_ID": "ID_NUMBER", "REGISTRY_ID": "REG"})
v = pd.read_csv(f"{ROOT}/data/interim/rcra/RCRA_VIOLATIONS.csv", dtype=str, encoding="latin-1", usecols=["ID_NUMBER", "DATE_VIOLATION_DETERMINED", "VIOL_DETERMINED_BY_AGENCY"]).merge(idmap, on="ID_NUMBER")
v["yr"] = pd.to_numeric(v.DATE_VIOLATION_DETERMINED.str[-4:], errors="coerce"); v = v.merge(P[["REG", "yr"]], on=["REG", "yr"])
v["src"] = np.where(v.VIOL_DETERMINED_BY_AGENCY.str.strip() == "E", "EPA", "state or local")
for src in ["state or local", "EPA"]:
    pl[f"ev_{src[:3]}"] = pl.index.isin(set(v[v.src == src].REG)).astype(float)
    b2.append(dict(program="rcra", detected_by=src, violations=int((v.src == src).sum()), **lpm(pl, f"ev_{src[:3]}")))
b2 = pd.DataFrame(b2); print(b2.round(4).to_string(index=False)); b2.to_csv(f"{OUT_T}/B2_detection_source.csv", index=False)

# ---------------- B3 ----------------
b3 = []
for prog in ["air", "rcra"]:
    P, pl = plant_frame(prog); pl = pl.join(attr)
    n3 = pl.n3.fillna("unk").astype(str)
    pl["sector"] = np.select([n3.str[:2].isin(["31", "32", "33"]), n3 == "221", n3.str[:2] == "21"], ["manufacturing", "utilities", "mining and extraction"], "other or unknown")
    pl["size"] = pd.cut(pl.nplant, [2, 5, 15, 10 ** 6], labels=["3-5 plants", "6-15 plants", "16+ plants"]).astype(str)
    pl["majgrp"] = np.where(pl.major == 1, "NPDES major", "NPDES non-major")
    for dim in ["majgrp", "sector", "size"]:
        ctrl = "" if dim == "majgrp" else " + major"  # major is collinear with the major-status groups
        f = smf.ols(f"ever ~ C({dim}) + c:C({dim}) + s:C({dim}){ctrl} + C(cls)", pl).fit(**cl(pl))
        levels = sorted(pl[dim].unique())
        wald = f.wald_test(" = ".join([f"c:C({dim})[{l}]" for l in levels]), scalar=True) if len(levels) > 1 else None
        for l in levels:
            sub = pl[pl[dim] == l]
            b3.append(dict(program=prog, dimension={"majgrp": "major status", "sector": "sector", "size": "firm size"}[dim], group=l, plants=len(sub), chronic_plants=int(sub.c.sum()),
                           mean=sub.ever.mean(), own=f.params[f"c:C({dim})[{l}]"], own_se=f.bse[f"c:C({dim})[{l}]"], sibling=f.params[f"s:C({dim})[{l}]"], sib_se=f.bse[f"s:C({dim})[{l}]"],
                           equality_of_own_p=float(wald.pvalue) if wald is not None else np.nan))
b3 = pd.DataFrame(b3); print(b3.round(4).to_string(index=False)); b3.to_csv(f"{OUT_T}/B3_heterogeneity.csv", index=False)

# ---------------- B4 ----------------
b4, b4c = [], []
rng = np.random.default_rng(151)
for prog, V in [("air", "any_air"), ("rcra", "any_rcra")]:
    P = pd.read_parquet(f"{ROOT}/extension/{prog}_panel.parquet"); CLS = "air_class" if prog == "air" else "rcra_class"
    for per, (a, b) in {"2010-2017": (2010, 2017), "2018-2025": (2018, 2025)}.items():
        Q = P[P.yr.between(a, b)]; Q = Q[Q.pnorm.map(Q.groupby("pnorm").REG.nunique()) >= 2]
        pl = Q.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), yb=("any_e90", "sum"), w=("any_e90", "mean"), ever=(V, "max"), a=(V, "mean"), cls=(CLS, "first")).join(grp[["major"]])
        pl["c"] = (pl.yb >= 3).astype(int); pl["s"] = ((pl.groupby("pnorm").c.transform("sum") - pl.c) > 0).astype(int); pl["cls"] = pl.cls.fillna("unknown")
        b4.append(dict(program=prog, period=per, plants=len(pl), chronic_plants=int(pl.c.sum()), **lpm(pl, "ever")))
        # E8-type correlations within the period
        pl = pl.reset_index(); wc = (pl.w - pl.w.mean()).values; ac = (pl.a - pl.a.mean()).values; st = pl.st.values; fm = pl.pnorm.values
        I, J = [], []
        for _, g in pl.groupby("pnorm"):
            ix = g.index.values; ii, jj = np.meshgrid(ix, ix); m = ii != jj; I.append(ii[m]); J.append(jj[m])
        I, J = np.concatenate(I), np.concatenate(J); same = st[I] == st[J]
        byst = pd.Series(np.arange(len(pl))).groupby(st).apply(np.array); byst = byst[byst.apply(len) >= 2]
        wts = byst.apply(lambda x: len(x) * (len(x) - 1)).values.astype(float); wts /= wts.sum(); UI, UJ = [], []
        for k in rng.choice(len(byst), 150000, p=wts):
            i, j = rng.choice(byst.iloc[k], 2, replace=False)
            if fm[i] != fm[j]: UI.append(i); UJ.append(j)
        c = lambda x, y: np.corrcoef(x, y)[0, 1]
        b4c.append(dict(program=prog, period=per, same_plant=c(wc, ac), sibling_same_state=c(wc[I][same], ac[J][same]), sibling_other_state=c(wc[I][~same], ac[J][~same]),
                        unrelated_same_state=c(wc[UI], ac[UJ])))
b4, b4c = pd.DataFrame(b4), pd.DataFrame(b4c); print(b4.round(4).to_string(index=False)); print(b4c.round(3).to_string(index=False))
b4.to_csv(f"{OUT_T}/B4_periods_lpm.csv", index=False); b4c.to_csv(f"{OUT_T}/B4_periods_correlations.csv", index=False)
