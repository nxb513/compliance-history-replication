"""Build plant-year RCRA panel (violations, SNC, evaluations) on the water panel's plant-years."""
import pandas as pd, numpy as np
ROOT = "."; R = f"{ROOT}/data/interim/rcra"
py = pd.read_parquet(f"{ROOT}/data/interim/one_bad_plant_panel.parquet"); pr = py[py.nplant >= 3].copy()
links = pd.read_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
idmap = links[links.PGM_SYS_ACRNM == "RCRAINFO"][["PGM_SYS_ID", "REGISTRY_ID"]].drop_duplicates().rename(columns={"PGM_SYS_ID": "ID_NUMBER"})
rd = lambda f, cols: pd.read_csv(f"{R}/{f}", dtype=str, encoding="latin-1", usecols=cols)
fac = rd("RCRA_FACILITIES.csv", ["ID_NUMBER", "FED_WASTE_GENERATOR", "OPERATING_TSDF"]).merge(idmap, on="ID_NUMBER")
fac["lqg"] = fac.FED_WASTE_GENERATOR.str.strip() == "1"
fac["tsdf"] = fac.OPERATING_TSDF.fillna("").str.strip().str.replace("-", "").str.len() > 0
gen = {"1": 3, "2": 2, "3": 1}
fac["g"] = fac.FED_WASTE_GENERATOR.str.strip().map(gen).fillna(0)
cls = fac.groupby("REGISTRY_ID").agg(g=("g", "max"), tsdf=("tsdf", "max"))
cls["rcra_class"] = np.where(cls.tsdf | (cls.g == 3), "LQG or TSDF", cls.g.map({2: "small quantity generator", 1: "very small quantity generator", 0: "not a current generator"}))
print("primary plants linked to RCRAInfo:", idmap.REGISTRY_ID.nunique(), "| class:", cls.rcra_class.value_counts().to_dict())
v = rd("RCRA_VIOLATIONS.csv", ["ID_NUMBER", "DATE_VIOLATION_DETERMINED", "VIOL_DETERMINED_BY_AGENCY"]).merge(idmap, on="ID_NUMBER")
v["yr"] = pd.to_numeric(v.DATE_VIOLATION_DETERMINED.str[-4:], errors="coerce"); v = v[v.yr.between(2010, 2025)]
print("RCRA violation records 2010-2025:", len(v), "| agency:", v.VIOL_DETERMINED_BY_AGENCY.str.strip().value_counts().to_dict())
va = v.groupby(["REGISTRY_ID", "yr"]).size().rename("n_rcra_viol").reset_index()
s = rd("RCRA_VIOSNC_HISTORY.csv", ["ID_NUMBER", "YRMONTH", "SNC_FLAG"]).merge(idmap, on="ID_NUMBER")
s["yr"] = pd.to_numeric(s.YRMONTH.str[:4], errors="coerce"); s = s[s.yr.between(2010, 2025) & (s.SNC_FLAG.str.strip() == "Y")]
sa = s.groupby(["REGISTRY_ID", "yr"]).size().gt(0).rename("any_rcra_snc").reset_index()
e = rd("RCRA_EVALUATIONS.csv", ["ID_NUMBER", "EVALUATION_START_DATE"]).merge(idmap, on="ID_NUMBER")
e["yr"] = pd.to_numeric(e.EVALUATION_START_DATE.str[-4:], errors="coerce"); e = e[e.yr.between(2010, 2025)]
ea = e.groupby(["REGISTRY_ID", "yr"]).size().rename("n_rcra_eval").reset_index()
ren = lambda d: d.rename(columns={"REGISTRY_ID": "REG"})
panel = pr[pr.REG.isin(set(idmap.REGISTRY_ID))].merge(ren(va), on=["REG", "yr"], how="left").merge(ren(sa), on=["REG", "yr"], how="left").merge(ren(ea), on=["REG", "yr"], how="left")
panel["n_rcra_viol"] = panel.n_rcra_viol.fillna(0); panel["any_rcra"] = (panel.n_rcra_viol > 0).astype(float)
panel["any_rcra_snc"] = panel.any_rcra_snc.fillna(False).astype(float); panel["n_rcra_eval"] = panel.n_rcra_eval.fillna(0)
panel["rcra_class"] = panel.REG.map(cls.rcra_class)
panel.to_parquet(f"{ROOT}/extension/rcra_panel.parquet")
print("rcra panel:", panel.shape, "plants", panel.REG.nunique(), "firms", panel.pnorm.nunique())
print(panel.groupby("rcra_class").agg(plants=("REG", "nunique"), any_rcra=("any_rcra", "mean"), any_snc=("any_rcra_snc", "mean"), evals=("n_rcra_eval", "mean")).round(3).to_string())
