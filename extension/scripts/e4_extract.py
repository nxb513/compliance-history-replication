"""E4 step 1: extract E90 effluent violations 2010-2025 for primary-sample plants from the 15.6 GB file."""
import duckdb, pandas as pd, time, os
ROOT = "."
py = pd.read_parquet(f"{ROOT}/data/interim/one_bad_plant_panel.parquet")
plants = set(py[py.nplant >= 3].REG.astype(str))
m = pd.read_parquet(f"{ROOT}/data/interim/us_npdes_facility_master.parquet", columns=["NPDES_ID", "REGISTRY_ID"])
m = m[m.REGISTRY_ID.astype(str).isin(plants)].drop_duplicates()
m.to_parquet(f"{ROOT}/extension/data_npdes_ids_primary.parquet")
print("NPDES permits of primary plants:", m.NPDES_ID.nunique(), "for plants:", m.REGISTRY_ID.nunique())
con = duckdb.connect()
con.execute("PRAGMA threads=8")
con.register("ids", m[["NPDES_ID", "REGISTRY_ID"]])
t = time.time()
src = f"{ROOT}/data/interim/part2/NPDES_EFF_VIOLATIONS.csv"
q = f"""
COPY (
  SELECT v.NPDES_ID, ids.REGISTRY_ID, v.PERM_FEATURE_NMBR, v.PARAMETER_CODE, v.PARAMETER_DESC,
         v.MONITORING_PERIOD_END_DATE, v.EXCEEDENCE_PCT, v.STATISTICAL_BASE_SHORT_DESC, v.VALUE_TYPE_CODE,
         v.RNC_DETECTION_CODE
  FROM read_csv('{src}', all_varchar=true, header=true, quote='"', ignore_errors=true) v
  JOIN ids ON v.NPDES_ID = ids.NPDES_ID
  WHERE v.VIOLATION_CODE = 'E90'
    AND TRY_CAST(right(v.MONITORING_PERIOD_END_DATE, 4) AS INTEGER) BETWEEN 2010 AND 2025
) TO '{ROOT}/extension/e90_param_primary.parquet' (FORMAT PARQUET)
"""
con.execute(q)
print("extract done in", round(time.time() - t), "s")
d = pd.read_parquet(f"{ROOT}/extension/e90_param_primary.parquet")
print("rows:", len(d), "plants:", d.REGISTRY_ID.nunique(), "parameters:", d.PARAMETER_CODE.nunique())
