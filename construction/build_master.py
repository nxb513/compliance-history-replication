import pandas as pd, numpy as np
OUT="./data/interim"
pm=pd.read_parquet(f"{OUT}/us_permit_master.parquet")
va=pd.read_parquet(f"{OUT}/us_violation_agg.parquet")
fe=pd.read_parquet(f"{OUT}/us_formal_enf_agg.parquet")
inf=pd.read_parquet(f"{OUT}/us_informal_agg.parquet")
insp=pd.read_parquet(f"{OUT}/us_insp_agg.parquet")
q=pd.read_parquet(f"{OUT}/us_qncr_agg.parquet")
na=pd.read_parquet(f"{OUT}/us_naics.parquet"); si=pd.read_parquet(f"{OUT}/us_sic.parquet")

m=pm.merge(va,on="NPDES_ID",how="left").merge(fe,on="NPDES_ID",how="left")\
    .merge(inf,on="NPDES_ID",how="left").merge(insp,on="NPDES_ID",how="left")\
    .merge(q,on="NPDES_ID",how="left").merge(na,on="NPDES_ID",how="left").merge(si,on="NPDES_ID",how="left")
# REGISTRY_ID from ICIS (inspections preferred, else informal); FRS added later
m["REGISTRY_ID"]=m.reg_insp.fillna(m.reg_inf)
# counts -> 0
for c in ["n_viol","n_e90","n_report","n_formal","n_informal","n_insp","n_qtr","n_years","qncr_e90","qncr_rep"]:
    m[c]=pd.to_numeric(m[c],errors="coerce").fillna(0)
for c in ["fed_penalty","state_penalty"]: m[c]=pd.to_numeric(m[c],errors="coerce")
# flags
m["is_potw"]=(m.FACILITY_TYPE_INDICATOR=="POTW").astype(int)
m["is_nonpotw"]=(m.FACILITY_TYPE_INDICATOR=="NON-POTW").astype(int)
m["is_federal"]=(m.FACILITY_TYPE_INDICATOR=="FEDERAL").astype(int)
m["is_major"]=(m.MAJOR_MINOR_STATUS_FLAG=="M").astype(int)
m["is_individual"]=(m.PERMIT_TYPE_CODE=="NPD").astype(int)
m["is_general"]=(m.PERMIT_TYPE_CODE=="GPC").astype(int)
m["is_active"]=m.PERMIT_STATUS_CODE.isin(["EFF","ADC"]).astype(int)
m["has_dmr"]=((m.n_qtr>0)|(m.n_viol>0)).astype(int)
m["has_e90"]=(m.n_e90>0).astype(int)
m["has_enf"]=((m.n_formal>0)|(m.n_informal>0)).astype(int)
m["has_penalty"]=((m.fed_penalty.fillna(0)>0)|(m.state_penalty.fillna(0)>0)).astype(int)
m["TRI_linked"]=np.nan  # filled after FRS
m["GHGRP_linked"]=np.nan
m.to_parquet(f"{OUT}/us_npdes_facility_master.parquet")
print("MASTER rows:",len(m),"cols:",m.shape[1])

def C(cond): return int(cond.sum())
print("\n===== BASIC COUNTS (national) =====")
print("total NPDES facilities:            ",len(m))
print("individual-NPDES (NPD):            ",C(m.is_individual==1))
print("NON-POTW individual-NPDES:         ",C((m.is_individual==1)&(m.is_nonpotw==1)))
print("  of which MAJOR:                  ",C((m.is_individual==1)&(m.is_nonpotw==1)&(m.is_major==1)))
print("with DMR/QNCR history:             ",C(m.has_dmr==1))
print("with any E90 violation:            ",C(m.has_e90==1))
print("with any enforcement action:       ",C(m.has_enf==1))
print("  formal:",C(m.n_formal>0),"| informal:",C(m.n_informal>0))
print("with monetary penalty:             ",C(m.has_penalty==1))
print("with >=1 inspection:               ",C(m.n_insp>0))
print("with REGISTRY_ID (ICIS-derived):   ",C(m.REGISTRY_ID.notna()))
print("\nby major/minor:"); print(m.MAJOR_MINOR_STATUS_FLAG.value_counts(dropna=False).head(6).to_string())
print("\nby permit type:"); print(m.PERMIT_TYPE_CODE.value_counts(dropna=False).head(6).to_string())
print("\nDMR history depth (QNCR years) among individual-NPD:")
ind=m[m.is_individual==1]
for yr in [1,3,5]: print(f"  >= {yr} yr: {C(ind.n_years>=yr)}")
print("\nNON-POTW individual-NPD DMR depth:")
ni=m[(m.is_individual==1)&(m.is_nonpotw==1)]
for yr in [1,3,5]: print(f"  >= {yr} yr: {C(ni.n_years>=yr)}")
print("  with E90:",C(ni.has_e90==1),"| with enforcement:",C(ni.has_enf==1),"| with penalty:",C(ni.has_penalty==1))
print("\ntop states (all NPDES):"); print(m.STATE_CODE.value_counts().head(8).to_string())

# JOIN A: facilities -> violations ; JOIN B: violations -> enforcement
ve=pd.read_csv("./data/interim/part1/NPDES_VIOLATION_ENFORCEMENTS.csv",usecols=["NPDES_VIOLATION_ID","ENF_IDENTIFIER"],dtype=str,encoding="latin-1")
print("\n===== JOIN DIAGNOSTICS (prelim) =====")
print("A. facilities with violations / total:",C(m.n_viol>0),"/",len(m),f"({C(m.n_viol>0)/len(m)*100:.1f}%)")
print("B. violation->enforcement link rows:",len(ve),"| distinct violation_ids:",ve.NPDES_VIOLATION_ID.nunique(),
      "| distinct enf_ids:",ve.ENF_IDENTIFIER.nunique(),"| dup violation_id rows:",len(ve)-ve.NPDES_VIOLATION_ID.nunique())
