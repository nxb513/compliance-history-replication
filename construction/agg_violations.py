import pandas as pd, numpy as np, time
P2="./data/interim/part2/NPDES_EFF_VIOLATIONS.csv"
t0=time.time(); parts=[]; n=0
uc=["NPDES_ID","VIOLATION_CODE","MONITORING_PERIOD_END_DATE","EXCEEDENCE_PCT"]
for i,ch in enumerate(pd.read_csv(P2,usecols=uc,dtype=str,encoding="latin-1",chunksize=3_000_000)):
    n+=len(ch)
    ch["d"]=pd.to_datetime(ch.MONITORING_PERIOD_END_DATE,errors="coerce")
    ch["e90"]=(ch.VIOLATION_CODE=="E90").astype("int32")
    ch["d80"]=(ch.VIOLATION_CODE=="D80").astype("int32")
    ch["d90"]=(ch.VIOLATION_CODE=="D90").astype("int32")
    g=ch.groupby("NPDES_ID").agg(n=("VIOLATION_CODE","size"),e90=("e90","sum"),d80=("d80","sum"),
                                 d90=("d90","sum"),dmin=("d","min"),dmax=("d","max"))
    parts.append(g.reset_index())
    if i%3==0: print(f"chunk {i} seen={n:,} t={time.time()-t0:.0f}s",flush=True)
allg=pd.concat(parts,ignore_index=True)
fin=allg.groupby("NPDES_ID").agg(n_viol=("n","sum"),n_e90=("e90","sum"),n_d80=("d80","sum"),
                                 n_d90=("d90","sum"),viol_first=("dmin","min"),viol_last=("dmax","max")).reset_index()
fin["n_report"]=fin.n_d80+fin.n_d90
fin.to_parquet("./data/interim/us_violation_agg.parquet")
print("DONE facilities with violations:",len(fin),"| total rows scanned:",f"{n:,}","| t=",int(time.time()-t0),"s")
print(fin[["n_viol","n_e90","n_report"]].sum().to_string())
