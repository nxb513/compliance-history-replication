import pandas as pd, numpy as np, time
P1="./data/interim/part1"; OUT="./data/interim"; t0=time.time()
def rd(f,**k): return pd.read_csv(f"{P1}/{f}",dtype=str,encoding="latin-1",**k)

# permit master (facilities + latest permit version)
fac=rd("ICIS_FACILITIES.csv",usecols=["NPDES_ID","FACILITY_NAME","FACILITY_TYPE_CODE","LOCATION_ADDRESS","CITY","COUNTY_CODE","STATE_CODE","ZIP","GEOCODE_LATITUDE","GEOCODE_LONGITUDE"])
perm=rd("ICIS_PERMITS.csv",usecols=["EXTERNAL_PERMIT_NMBR","VERSION_NMBR","FACILITY_TYPE_INDICATOR","PERMIT_TYPE_CODE","MAJOR_MINOR_STATUS_FLAG","PERMIT_STATUS_CODE","ISSUING_AGENCY","TOTAL_DESIGN_FLOW_NMBR","ORIGINAL_ISSUE_DATE","EFFECTIVE_DATE","EXPIRATION_DATE","TERMINATION_DATE"])
perm["VERSION_NMBR"]=pd.to_numeric(perm.VERSION_NMBR,errors="coerce")
perm=perm.sort_values("VERSION_NMBR").drop_duplicates("EXTERNAL_PERMIT_NMBR",keep="last")
pm=fac.merge(perm,left_on="NPDES_ID",right_on="EXTERNAL_PERMIT_NMBR",how="left")
pm["POTW_FLAG"]=np.where(pm.FACILITY_TYPE_INDICATOR=="POTW",1,np.where(pm.FACILITY_TYPE_INDICATOR=="NON-POTW",0,np.nan))
pm.to_parquet(f"{OUT}/us_permit_master.parquet"); print("permit master:",len(pm),"| non-POTW:",(pm.FACILITY_TYPE_INDICATOR=='NON-POTW').sum(),"| NPD:",(pm.PERMIT_TYPE_CODE=='NPD').sum())

# formal enforcement agg
fe=rd("NPDES_FORMAL_ENFORCEMENT_ACTIONS.csv")
fe["fed"]=pd.to_numeric(fe.FED_PENALTY_ASSESSED_AMT,errors="coerce"); fe["sl"]=pd.to_numeric(fe.STATE_LOCAL_PENALTY_AMT,errors="coerce")
fe["d"]=pd.to_datetime(fe.SETTLEMENT_ENTERED_DATE,errors="coerce")
feA=fe.groupby("NPDES_ID").agg(n_formal=("ENF_IDENTIFIER","nunique"),enf_first=("d","min"),enf_last=("d","max"),
    fed_penalty=("fed","sum"),state_penalty=("sl","sum")).reset_index()
feA.to_parquet(f"{OUT}/us_formal_enf_agg.parquet"); print("facilities w/ formal enf:",len(feA))

# informal + inspections (also give REGISTRY_ID)
inf=rd("NPDES_INFORMAL_ENFORCEMENT_ACTIONS.csv",usecols=["NPDES_ID","REGISTRY_ID","ENF_IDENTIFIER"])
infA=inf.groupby("NPDES_ID").agg(n_informal=("ENF_IDENTIFIER","nunique"),reg_inf=("REGISTRY_ID","first")).reset_index()
insp=rd("NPDES_INSPECTIONS.csv",usecols=["NPDES_ID","REGISTRY_ID","ACTIVITY_ID"])
inspA=insp.groupby("NPDES_ID").agg(n_insp=("ACTIVITY_ID","nunique"),reg_insp=("REGISTRY_ID","first")).reset_index()
infA.to_parquet(f"{OUT}/us_informal_agg.parquet"); inspA.to_parquet(f"{OUT}/us_insp_agg.parquet")
print("facilities w/ informal:",len(infA),"| w/ inspections:",len(inspA))

# QNCR history (quarterly noncompliance counts)
q=rd("NPDES_QNCR_HISTORY.csv")
for c in ["NUME90Q","NUMD8090Q","NUMCVDT","NUMSVCD","NUMPSCH"]: q[c]=pd.to_numeric(q[c],errors="coerce").fillna(0)
q["yr"]=q.YEARQTR.str[:4]
qA=q.groupby("NPDES_ID").agg(n_qtr=("YEARQTR","nunique"),qtr_first=("YEARQTR","min"),qtr_last=("YEARQTR","max"),
    n_years=("yr","nunique"),qncr_e90=("NUME90Q","sum"),qncr_rep=("NUMD8090Q","sum")).reset_index()
qA.to_parquet(f"{OUT}/us_qncr_agg.parquet"); print("facilities in QNCR:",len(qA))

# NAICS / SIC primary
na=rd("NPDES_NAICS.csv")
naP=na[na.PRIMARY_INDICATOR_FLAG=="Y"].drop_duplicates("NPDES_ID")[["NPDES_ID","NAICS_CODE"]].rename(columns={"NAICS_CODE":"naics"})
si=rd("NPDES_SICS.csv"); siP=si[si.PRIMARY_INDICATOR_FLAG=="Y"].drop_duplicates("NPDES_ID")[["NPDES_ID","SIC_CODE"]].rename(columns={"SIC_CODE":"sic"})
naP.to_parquet(f"{OUT}/us_naics.parquet"); siP.to_parquet(f"{OUT}/us_sic.parquet")
print("naics primary:",len(naP),"| sic primary:",len(siP),"| t=",int(time.time()-t0),"s")
