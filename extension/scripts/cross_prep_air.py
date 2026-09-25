"""Build plant-year Clean Air Act panel (violations, HPVs, evaluations) on the water panel's plant-years."""
import pandas as pd, numpy as np
ROOT = "."; A = f"{ROOT}/data/interim/air"
py = pd.read_parquet(f"{ROOT}/data/interim/one_bad_plant_panel.parquet"); pr = py[py.nplant >= 3].copy()
plants = set(pr.REG.astype(str))
fac = pd.read_csv(f"{A}/ICIS-AIR_FACILITIES.csv", dtype=str, usecols=["PGM_SYS_ID", "REGISTRY_ID", "AIR_POLLUTANT_CLASS_CODE", "AIR_OPERATING_STATUS_CODE"])
fac = fac[fac.REGISTRY_ID.isin(plants)]
rank = {"MAJ": 3, "SMI": 2, "MIN": 1}
fac["cls_rank"] = fac.AIR_POLLUTANT_CLASS_CODE.map(rank).fillna(0)
cls = fac.groupby("REGISTRY_ID").cls_rank.max().map({3: "major", 2: "synthetic minor", 1: "minor", 0: "other/unknown"})
print("primary plants with ICIS-Air ID:", fac.REGISTRY_ID.nunique(), "| class:", cls.value_counts().to_dict())
idmap = fac[["PGM_SYS_ID", "REGISTRY_ID"]].drop_duplicates()
vh = pd.read_csv(f"{A}/ICIS-AIR_VIOLATION_HISTORY.csv", dtype=str, usecols=["PGM_SYS_ID", "ENF_RESPONSE_POLICY_CODE", "EARLIEST_FRV_DETERM_DATE", "HPV_DAYZERO_DATE", "AGENCY_TYPE_DESC"])
vh = vh.merge(idmap, on="PGM_SYS_ID")
vh["date"] = np.where(vh.ENF_RESPONSE_POLICY_CODE == "HPV", vh.HPV_DAYZERO_DATE.fillna(vh.EARLIEST_FRV_DETERM_DATE), vh.EARLIEST_FRV_DETERM_DATE.fillna(vh.HPV_DAYZERO_DATE))
vh["yr"] = pd.to_numeric(vh.date.str[-4:], errors="coerce")
vh = vh[vh.yr.between(2010, 2025)]
print("air violation records 2010-2025 at primary plants:", len(vh), "| by type:", vh.ENF_RESPONSE_POLICY_CODE.value_counts().to_dict(), "| agency:", vh.AGENCY_TYPE_DESC.value_counts().to_dict())
va = vh.groupby(["REGISTRY_ID", "yr"]).agg(n_air_viol=("PGM_SYS_ID", "size"), any_hpv=("ENF_RESPONSE_POLICY_CODE", lambda s: (s == "HPV").any())).reset_index()
ev = pd.read_csv(f"{A}/ICIS-AIR_FCES_PCES.csv", dtype=str, usecols=["PGM_SYS_ID", "ACTUAL_END_DATE", "COMP_MONITOR_TYPE_CODE"]).merge(idmap, on="PGM_SYS_ID")
ev["yr"] = pd.to_numeric(ev.ACTUAL_END_DATE.str[-4:], errors="coerce"); ev = ev[ev.yr.between(2010, 2025)]
ea = ev.groupby(["REGISTRY_ID", "yr"]).size().rename("n_air_eval").reset_index()
panel = pr[pr.REG.isin(set(fac.REGISTRY_ID))].merge(va.rename(columns={"REGISTRY_ID": "REG"}), on=["REG", "yr"], how="left").merge(ea.rename(columns={"REGISTRY_ID": "REG"}), on=["REG", "yr"], how="left")
panel["n_air_viol"] = panel.n_air_viol.fillna(0); panel["any_air"] = (panel.n_air_viol > 0).astype(float)
panel["any_hpv"] = panel.any_hpv.fillna(False).astype(float); panel["n_air_eval"] = panel.n_air_eval.fillna(0)
panel["air_class"] = panel.REG.map(cls)
panel.to_parquet(f"{ROOT}/extension/air_panel.parquet")
print("air panel:", panel.shape, "plants", panel.REG.nunique(), "firms", panel.pnorm.nunique())
print("plant-years with any air violation:", round(panel.any_air.mean(), 3), "| any HPV:", round(panel.any_hpv.mean(), 3))
print(panel.groupby("air_class").agg(plants=("REG", "nunique"), any_air=("any_air", "mean"), any_hpv=("any_hpv", "mean"), evals=("n_air_eval", "mean")).round(3).to_string())
