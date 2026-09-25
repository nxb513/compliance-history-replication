"""E21: do assessed penalties move with the violation history of the firm's other plants, holding the plant's own history fixed?
Pre-specified in Addendum 8 of analysis_plan.md. Unit = plant x program x year with at least one formal action, 2013-2025."""
import sys, os, zipfile, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
YR = (2013, 2025)
num = lambda s: pd.to_numeric(s, errors="coerce").fillna(0.0)

# ---- plant-year violation histories (water for all primary plants; air and RCRA for plants in those programs)
P = load_primary()[["REG", "yr", "pnorm", "STATE_CODE", "is_major", "nplant", "any_e90"]].rename(columns={"any_e90": "v_water"})
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet")[["REG", "yr", "any_air", "air_class"]].rename(columns={"any_air": "v_air"})
rc = pd.read_parquet(f"{ROOT}/extension/rcra_panel.parquet")[["REG", "yr", "any_rcra", "rcra_class"]].rename(columns={"any_rcra": "v_rcra"})
H = P.merge(air, on=["REG", "yr"], how="left").merge(rc, on=["REG", "yr"], how="left")
in_air, in_rcra = set(air.REG), set(rc.REG)
H["v_any"] = H[["v_water", "v_air", "v_rcra"]].max(axis=1)
V = {p: H.set_index(["REG", "yr"])[f"v_{p}"] for p in ["water", "air", "rcra", "any"]}
def window(series, r, t):
    vals = [series.get((r, y)) for y in (t - 3, t - 2, t - 1)]
    vals = [v for v in vals if v is not None and not pd.isna(v)]
    return np.mean(vals) if vals else np.nan
plants = H.drop_duplicates("REG").set_index("REG")[["pnorm", "STATE_CODE", "is_major", "nplant", "air_class", "rcra_class"]]
plants["air_class"] = H.dropna(subset=["air_class"]).drop_duplicates("REG").set_index("REG").air_class
plants["rcra_class"] = H.dropna(subset=["rcra_class"]).drop_duplicates("REG").set_index("REG").rcra_class
byfirm = plants.groupby("pnorm").apply(lambda g: list(zip(g.index, g.STATE_CODE))).to_dict()

# ---- formal actions by program -> units (plant, year)
ids = pd.read_parquet(f"{ROOT}/extension/data_npdes_ids_primary.parquet")
fw = pd.read_parquet(f"{ROOT}/data/interim/formal_enf.parquet").merge(ids, on="NPDES_ID")
fw = fw.assign(yr=pd.to_numeric(fw.SETTLEMENT_ENTERED_DATE.str[-4:], errors="coerce"), pen=num(fw.FED_PENALTY_ASSESSED_AMT) + num(fw.STATE_LOCAL_PENALTY_AMT),
               epa=(fw.AGENCY == "EPA").astype(float), jud=(fw.ACTIVITY_TYPE_CODE == "JDC").astype(float))
fac = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FACILITIES.csv", dtype=str, usecols=["PGM_SYS_ID", "REGISTRY_ID"])
fa = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FORMAL_ACTIONS.csv", dtype=str).merge(fac, on="PGM_SYS_ID")
fa = fa.assign(yr=pd.to_numeric(fa.SETTLEMENT_ENTERED_DATE.str[-4:], errors="coerce"), pen=num(fa.PENALTY_AMOUNT),
               epa=(fa.STATE_EPA_FLAG.str.strip() == "E").astype(float), jud=(fa.ACTIVITY_TYPE_CODE == "JDC").astype(float))
links = pd.read_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
rmap = links[links.PGM_SYS_ACRNM == "RCRAINFO"][["PGM_SYS_ID", "REGISTRY_ID"]].drop_duplicates().rename(columns={"PGM_SYS_ID": "ID_NUMBER"})
re_ = pd.read_csv(zipfile.ZipFile(f"{ROOT}/data/raw/echo_crossmedia/rcra_downloads.zip").open("RCRA_ENFORCEMENTS.csv"), dtype=str, encoding="latin-1").merge(rmap, on="ID_NUMBER")
re_["code"] = pd.to_numeric(re_.ENFORCEMENT_TYPE.str.strip().str[-3:], errors="coerce")
re_ = re_[re_.code.between(200, 399) | re_.code.between(500, 699)].copy()
re_ = re_.assign(yr=pd.to_numeric(re_.ENFORCEMENT_ACTION_DATE.str[-4:], errors="coerce"), fmp=pd.to_numeric(re_.FMP_AMOUNT, errors="coerce"), pmp=pd.to_numeric(re_.PMP_AMOUNT, errors="coerce"),
                 epa=(re_.ENFORCEMENT_AGENCY.str.strip() == "E").astype(float), jud=re_.code.between(500, 699).astype(float))
def units(df, prog):
    df = df[df.yr.between(*YR) & df.REGISTRY_ID.isin(plants.index)]
    if prog == "rcra":
        g = df.groupby(["REGISTRY_ID", "yr"]).agg(fmp=("fmp", lambda s: s.sum(min_count=1)), pmp=("pmp", lambda s: s.sum(min_count=1)), epa=("epa", "max"), jud=("jud", "max"), n_act=("code", "size"))
        g["pen"] = g.fmp.where(g.fmp.notna(), g.pmp).fillna(0.0); g = g.drop(columns=["fmp", "pmp"])
    else:
        g = df.groupby(["REGISTRY_ID", "yr"]).agg(pen=("pen", "sum"), epa=("epa", "max"), jud=("jud", "max"), n_act=("pen", "size"))
    g = g.reset_index().rename(columns={"REGISTRY_ID": "REG"}); g["program"] = prog
    return g
U = pd.concat([units(fw, "water"), units(fa, "air"), units(re_, "rcra")], ignore_index=True)
U = U[(U.program == "water") | ((U.program == "air") & U.REG.isin(in_air)) | ((U.program == "rcra") & U.REG.isin(in_rcra))]

# ---- histories
rows = []
for r, t, prog in zip(U.REG, U.yr, U.program):
    f, s = plants.at[r, "pnorm"], plants.at[r, "STATE_CODE"]
    own_same = window(V[prog], r, t)
    others = [p for p in ["water", "air", "rcra"] if p != prog and (p == "water" or (p == "air" and r in in_air) or (p == "rcra" and r in in_rcra))]
    oo = [window(V[p], r, t) for p in others]; oo = [x for x in oo if not pd.isna(x)]
    ss = [window(V["any"], q, t) for q, qs in byfirm[f] if q != r and qs == s]; so = [window(V["any"], q, t) for q, qs in byfirm[f] if q != r and qs != s]
    ss = [x for x in ss if not pd.isna(x)]; so = [x for x in so if not pd.isna(x)]
    rows.append(dict(own_same=own_same, in_other=float(len(others) > 0), own_other=float(any(x > 0 for x in oo)) if oo else 0.0,
                     has_ss=float(len(ss) > 0), sib_same_state=np.mean([x > 0 for x in ss]) if ss else 0.0,
                     has_so=float(len(so) > 0), sib_other_state=np.mean([x > 0 for x in so]) if so else 0.0))
U = pd.concat([U.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
U = U.merge(plants[["pnorm", "STATE_CODE", "is_major", "nplant", "air_class", "rcra_class"]], left_on="REG", right_index=True)
U["lnplant"] = np.log(U.nplant); U["any_pen"] = (U.pen > 0).astype(float); U["lpen"] = np.where(U.pen > 0, np.log(U.pen.clip(lower=1)), np.nan)
U["cls"] = np.select([U.program == "water", U.program == "air", U.program == "rcra"],
                     [U.is_major.astype(str), U.air_class.fillna("na").astype(str), U.rcra_class.fillna("na").astype(str)], "na")
U = U.dropna(subset=["own_same"])
U.to_parquet(f"{ROOT}/extension/e21_units.parquet")
desc = U.groupby("program").agg(units=("REG", "size"), plants=("REG", "nunique"), firms=("pnorm", "nunique"), share_penalty=("any_pen", "mean"),
                                median_penalty_if_any=("pen", lambda s: s[s > 0].median()), share_epa=("epa", "mean"), share_judicial=("jud", "mean"),
                                mean_sib_same=("sib_same_state", "mean"), mean_sib_other=("sib_other_state", "mean"), mean_own_same=("own_same", "mean"))
pd.set_option("display.width", 250); print(desc.round(3).to_string(), flush=True)

# ---- regressions
X = "own_same + own_other + in_other + sib_same_state + has_ss + sib_other_state + has_so + epa + jud + lnplant + C(cls) + C(yr) + C(STATE_CODE)"
out = []
for prog in ["water", "air", "rcra", "pooled"]:
    D0 = U if prog == "pooled" else U[U.program == prog]
    for y in ["any_pen", "lpen"]:
        D = D0.dropna(subset=[y]).copy()
        Xp = X if D.in_other.nunique() > 1 else X.replace("in_other + ", "")  # in_other is constant (=1) for air and RCRA plants: every primary plant is an NPDES plant
        fml = f"{y} ~ {Xp}" + (" + C(program)" if prog == "pooled" else "")
        try:
            m = smf.ols(fml, D).fit(cov_type="cluster", cov_kwds={"groups": D.pnorm.astype("category").cat.codes})
            ex = m.model.exog; print(prog, y, "design rank", np.linalg.matrix_rank(ex), "of", ex.shape[1], flush=True)
        except Exception as ex:
            print(prog, y, "failed:", ex); continue
        for v in ["own_same", "own_other", "sib_same_state", "sib_other_state"]:
            out.append(dict(program=prog, outcome=y, variable=v, coef=m.params[v], se=m.bse[v], p=m.pvalues[v], n=int(m.nobs), firms=D.pnorm.nunique(), mean_y=D[y].mean()))
R = pd.DataFrame(out)
sib = R[R.variable.str.startswith("sib") & (R.program != "pooled")].index
R["p_holm"] = np.nan; R["p_bh"] = np.nan
R.loc[sib, "p_holm"] = multipletests(R.loc[sib, "p"], method="holm")[1]; R.loc[sib, "p_bh"] = multipletests(R.loc[sib, "p"], method="fdr_bh")[1]
print(R.round(4).to_string(index=False))
R.to_csv(f"{OUT_T}/E21_penalty_history.csv", index=False); desc.to_csv(f"{OUT_T}/E21_units_descriptives.csv")
