import pandas as pd, numpy as np, glob
OUT="./data/interim"
pl=pd.read_csv(f"{OUT}/frs/FRS_PROGRAM_LINKS.csv",usecols=["PGM_SYS_ACRNM","PGM_SYS_ID","REGISTRY_ID"],dtype=str,encoding="latin-1")
npd=pl[pl.PGM_SYS_ACRNM=="NPDES"].dropna(subset=["PGM_SYS_ID","REGISTRY_ID"])
tri_regs=set(pl.loc[pl.PGM_SYS_ACRNM=="TRIS","REGISTRY_ID"].dropna())
ghg_regs=set(pl.loc[pl.PGM_SYS_ACRNM.isin(["E-GGRT","GGRT","GHGRP"]),"REGISTRY_ID"].dropna())
print("NPDES program links:",len(npd),"| distinct NPDES_ID:",npd.PGM_SYS_ID.nunique(),
      "| TRIS regs:",len(tri_regs),"| GHGRP regs:",len(ghg_regs))
# NPDES_ID -> REGISTRY_ID (first)
npd1=npd.sort_values("REGISTRY_ID").drop_duplicates("PGM_SYS_ID").rename(columns={"PGM_SYS_ID":"NPDES_ID","REGISTRY_ID":"REG_FRS"})[["NPDES_ID","REG_FRS"]]

m=pd.read_parquet(f"{OUT}/us_npdes_facility_master.parquet")
m=m.merge(npd1,on="NPDES_ID",how="left")
# JOIN C
mc=m.REG_FRS.notna()
print(f"\nJOIN C facilities->FRS: matched {mc.sum():,}/{len(m):,} ({mc.mean()*100:.1f}%)")
# finalize REGISTRY_ID: prefer FRS, fallback ICIS-derived
m["REGISTRY_ID"]=m.REG_FRS.fillna(m.REGISTRY_ID)
m["TRI_linked"]=m.REG_FRS.isin(tri_regs).astype(int)  # via FRS registry
m["GHGRP_linked"]=m.REG_FRS.isin(ghg_regs).astype(int)
m.drop(columns=["REG_FRS"]).to_parquet(f"{OUT}/us_npdes_facility_master.parquet")
print("master updated. REGISTRY_ID coverage:",f"{m.REGISTRY_ID.notna().mean()*100:.1f}%")

# JOIN E validation via TRI files' own FRS ID
tri_frs=set()
for f in glob.glob(f"{OUT.replace('interim','raw')}/tri/TRI_*_US.csv"):
    d=pd.read_csv(f,usecols=["3. FRS ID"],dtype=str,encoding="latin-1")
    tri_frs.update(d["3. FRS ID"].dropna().unique())
print(f"\nJOIN E: distinct FRS IDs across TRI 2010-2024: {len(tri_frs):,}")
npdes_regs=set(m.REGISTRY_ID.dropna())
inter=tri_frs & npdes_regs
print(f"TRI FRS-IDs also an NPDES-facility REGISTRY_ID: {len(inter):,} ({len(inter)/len(tri_frs)*100:.1f}% of TRI facilities)")

# firm-level relevance: NPDES non-POTW industrial that are TRI-linked
print("\n=== TRI/GHGRP linkage of the target population ===")
print("all NPDES facilities TRI-linked:",int(m.TRI_linked.sum()),"| GHGRP-linked:",int(m.GHGRP_linked.sum()))
ni=m[(m.is_individual==1)&(m.is_nonpotw==1)]
print("non-POTW individual-NPD: TRI-linked",int(ni.TRI_linked.sum()),"| GHGRP-linked",int(ni.GHGRP_linked.sum()))
nie=ni[ni.has_e90==1]
print("non-POTW individual-NPD WITH E90: TRI-linked",int(nie.TRI_linked.sum()),"of",len(nie))
# facilities->firm via REGISTRY: how many NPDES per registry among target
reg_np=ni.groupby("REGISTRY_ID").NPDES_ID.nunique()
print("distinct REGISTRY (firms/sites) in non-POTW-NPD:",ni.REGISTRY_ID.nunique(),"| multi-permit registries:",int((reg_np>1).sum()))
