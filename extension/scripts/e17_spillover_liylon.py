"""E17: intra-firm enforcement spillovers in the Li & Lyon (2026) design: violation in t on penalties in t-1 at the plant and at four sibling categories
(same/different industry x same/different state), plant and year fixed effects, firm-clustered SE. Air (ICIS-Air) and water (NPDES)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
nai = pd.read_csv(f"{ROOT}/data/interim/frs/FRS_NAICS_CODES.csv", dtype=str, usecols=["REGISTRY_ID", "NAICS_CODE"])
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet")
nai = nai[nai.REGISTRY_ID.isin(set(grp.index))]; nai["n6"] = nai.NAICS_CODE.str.strip().str[:6]; nai["n4"] = nai.NAICS_CODE.str.strip().str[:4]
n6 = nai[nai.n6.str.len() == 6].groupby("REGISTRY_ID").n6.apply(set); n4 = nai[nai.n4.str.len() == 4].groupby("REGISTRY_ID").n4.apply(set)

# penalties / formal actions by plant-year
fac = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FACILITIES.csv", dtype=str, usecols=["PGM_SYS_ID", "REGISTRY_ID"])
fa = pd.read_csv(f"{ROOT}/data/interim/air/ICIS-AIR_FORMAL_ACTIONS.csv", dtype=str, usecols=["PGM_SYS_ID", "SETTLEMENT_ENTERED_DATE", "PENALTY_AMOUNT"]).merge(fac, on="PGM_SYS_ID")
fa["yr"] = pd.to_numeric(fa.SETTLEMENT_ENTERED_DATE.str[-4:], errors="coerce"); fa["pen"] = pd.to_numeric(fa.PENALTY_AMOUNT, errors="coerce").fillna(0) > 0
air_act = {"penalty": fa[fa.pen].groupby(["REGISTRY_ID", "yr"]).size(), "any formal action": fa.groupby(["REGISTRY_ID", "yr"]).size()}
ids = pd.read_parquet(f"{ROOT}/extension/data_npdes_ids_primary.parquet")
fw = pd.read_parquet(f"{ROOT}/data/interim/formal_enf.parquet", columns=["NPDES_ID", "SETTLEMENT_ENTERED_DATE", "FED_PENALTY_ASSESSED_AMT", "STATE_LOCAL_PENALTY_AMT"]).merge(ids, on="NPDES_ID")
fw["yr"] = pd.to_numeric(fw.SETTLEMENT_ENTERED_DATE.str[-4:], errors="coerce")
fw["pen"] = (pd.to_numeric(fw.FED_PENALTY_ASSESSED_AMT, errors="coerce").fillna(0) + pd.to_numeric(fw.STATE_LOCAL_PENALTY_AMT, errors="coerce").fillna(0)) > 0
water_act = {"penalty": fw[fw.pen].groupby(["REGISTRY_ID", "yr"]).size(), "any formal action": fw.groupby(["REGISTRY_ID", "yr"]).size()}
print("air penalties at primary plants 2010-2025:", int((air_act["penalty"].reset_index().yr.between(2010, 2025)).sum()), "| water penalties:", int((water_act["penalty"].reset_index().yr.between(2010, 2025)).sum()))

def build_x(P, act, ind):
    """P: panel (REG, yr, pnorm, STATE_CODE). act: Series indexed (REG, yr) of action counts. Returns lagged self and 4 sibling-category indicators."""
    A = set(act.index); P = P.copy()
    P["self"] = [((r, y - 1) in A) * 1.0 for r, y in zip(P.REG, P.yr)]
    plants = P.drop_duplicates("REG").set_index("REG")[["pnorm", "STATE_CODE"]]
    byfirm = plants.groupby("pnorm").apply(lambda g: list(zip(g.index, g.STATE_CODE)))
    acted = pd.DataFrame(list(A), columns=["REG", "yr"]); acted = acted[acted.REG.isin(plants.index)]
    acted_by = acted.groupby("yr").REG.apply(set).to_dict()
    cats = {c: np.zeros(len(P)) for c in ["sib_In1St1", "sib_In1St2", "sib_In2St1", "sib_In2St2"]}
    for k, (r, y, f, s) in enumerate(zip(P.REG, P.yr, P.pnorm, P.STATE_CODE)):
        prev = acted_by.get(y - 1)
        if not prev: continue
        ir = ind.get(r, set())
        for (q, qs) in byfirm[f]:
            if q == r or q not in prev: continue
            same_ind = len(ir & ind.get(q, set())) > 0; same_st = qs == s
            cats[f"sib_In{1 if same_ind else 2}St{1 if same_st else 2}"][k] = 1.0
    for c, v in cats.items(): P[c] = v
    return P

rows = []
setups = [("air", "any_air", air_act), ("water", "any_e90", water_act), ("air outcome, water sibling actions", "any_air", water_act)]
for prog, V, acts in setups:
    if prog.startswith("air"):
        P = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet")
    else:
        P = load_primary()
    P = P[P.pnorm.map(P.groupby("pnorm").REG.nunique()) >= 2]; P = P[P.yr >= 2011].reset_index(drop=True)
    for actname in ["penalty", "any formal action"]:
        for indname, ind in [("NAICS-6", n6), ("NAICS-4", n4)]:
            if actname != "penalty" and indname != "NAICS-6": continue
            X = build_x(P, acts[actname], ind)
            if prog == "air outcome, water sibling actions":  # own-plant regressor from the air program, sibling actions from water
                X["self"] = build_x(P, air_act[actname], ind)["self"]
            xs = ["self", "sib_In1St1", "sib_In1St2", "sib_In2St1", "sib_In2St2"]
            D = twoway_demean(X, [V] + xs, "REG", "yr")
            f = sm.OLS(D[V], D[xs]).fit(cov_type="cluster", cov_kwds={"groups": X.pnorm.astype("category").cat.codes})
            r = dict(program=prog, action=actname, industry=indname, plant_years=len(X), plants=X.REG.nunique(), firms=X.pnorm.nunique(), outcome_mean=X[V].mean())
            for x in xs: r[x] = f.params[x]; r[x + "_se"] = f.bse[x]; r[x + "_share_treated"] = X[x].mean()
            rows.append(r); print(prog, actname, indname, {x: (round(f.params[x], 4), round(f.bse[x], 4)) for x in xs}, flush=True)
out = pd.DataFrame(rows); out.to_csv(f"{OUT_T}/E17_liylon_spillover.csv", index=False)
pd.set_option("display.width", 300); print(out.round(4).to_string(index=False))
