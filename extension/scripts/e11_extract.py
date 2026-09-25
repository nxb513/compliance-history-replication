"""Extract ECHO Exporter rows for primary-sample plants (by FRS REGISTRY_ID)."""
import duckdb, pandas as pd
ROOT = "."
grp = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet").reset_index()
con = duckdb.connect(); con.register("g", grp[["REG"]].astype(str))
cols = ["REGISTRY_ID","FAC_STATE","FAC_LAT","FAC_LONG","FAC_DERIVED_STCTY_FIPS","FAC_DERIVED_CB2010","FAC_DERIVED_HUC","FAC_PERCENT_MINORITY","FAC_POP_DEN",
        "FAC_NAA_FLAG","FAC_IMP_WATER_FLG","FAC_CHESAPEAKE_BAY_FLG","FAC_US_MEX_BORDER_FLG","FAC_INDIAN_CNTRY_FLG","FAC_MAJOR_FLAG","FAC_ACTIVE_FLAG",
        "FAC_INSPECTION_COUNT","FAC_INFORMAL_COUNT","FAC_FORMAL_ACTION_COUNT","FAC_TOTAL_PENALTIES","FAC_PENALTY_COUNT","FAC_QTRS_WITH_NC","FAC_PROGRAMS_WITH_SNC",
        "AIR_FLAG","NPDES_FLAG","RCRA_FLAG","TRI_FLAG","GHG_FLAG","CWA_INSPECTION_COUNT","CWA_INFORMAL_COUNT","CWA_FORMAL_ACTION_COUNT","CWA_PENALTIES","CWA_QTRS_WITH_NC",
        "CWA_SNC_FLAG","CWA_13QTRS_COMPL_HISTORY","CAA_EVALUATION_COUNT","CAA_FORMAL_ACTION_COUNT","CAA_PENALTIES","RCRA_INSPECTION_COUNT","RCRA_FORMAL_ACTION_COUNT","RCRA_PENALTIES",
        "GHG_CO2_RELEASES","TRI_ON_SITE_RELEASES","FEC_NUMBER_OF_CASES","FEC_TOTAL_PENALTIES","FAC_DATE_LAST_INSPECTION_EPA","FAC_DATE_LAST_INSPECTION_STATE"]
sel = ", ".join(f'e."{c}"' for c in cols)
df = con.execute(f"""select {sel} from read_csv('{ROOT}/data/interim/echo/ECHO_EXPORTER.csv', all_varchar=true, header=true, ignore_errors=true) e
                     join g on e.REGISTRY_ID = g.REG""").df()
df.to_parquet(f"{ROOT}/extension/echo_exporter_primary.parquet")
print(df.shape, "unique plants", df.REGISTRY_ID.nunique(), "of", len(grp))
for c in ["FAC_PERCENT_MINORITY","FAC_POP_DEN","FAC_NAA_FLAG","FAC_IMP_WATER_FLG","CWA_FORMAL_ACTION_COUNT","CWA_PENALTIES","FAC_DERIVED_CB2010"]:
    s = df[c]; print(c, "nonnull", s.notna().sum(), "| sample", s.dropna().head(5).tolist())
