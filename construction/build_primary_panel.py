import pandas as pd, numpy as np, glob, re, warnings, time; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import statsmodels.formula.api as smf, statsmodels.api as sm
OUT="./data/interim"; FEC=f"{OUT}/fec"; P1=f"{OUT}/part1"; FIG="./figures"; TAB="./tables"
rng=np.random.default_rng(5)
def hdr(t): print("\n"+"="*66+"\n"+t+"\n"+"="*66)

# ========== BUILD plant x year panel ==========
tcols=["3. FRS ID","17. STANDARD PARENT CO NAME","15. PARENT CO NAME","16. PARENT CO DB NUM","21. FEDERAL FACILITY"]
tri=pd.concat([pd.read_csv(f,usecols=tcols,dtype=str,encoding="latin-1") for f in glob.glob("./data/raw/tri/TRI_*_US.csv")],ignore_index=True).dropna(subset=["3. FRS ID"])
tri["pname"]=tri["17. STANDARD PARENT CO NAME"].fillna(tri["15. PARENT CO NAME"]).astype(str).str.upper().str.strip()
def norm(s):
    s=re.sub(r"[^A-Z0-9 ]"," ",str(s)); s=re.sub(r"\b(INC|INCORPORATED|LLC|L L C|CORP|CORPORATION|CO|COMPANY|LP|L P|LTD|LIMITED|HOLDINGS|HOLDING|GROUP|THE|USA|US|PLC|SA|NV|AG|INTERNATIONAL|INTL|ENTERPRISES|INDUSTRIES)\b"," ",s); return re.sub(r"\s+"," ",s).strip()
tri["pnorm"]=tri.pname.map(norm)
tl=tri.sort_values("16. PARENT CO DB NUM").drop_duplicates("3. FRS ID",keep="last")
reg2pn=dict(zip(tl["3. FRS ID"],tl.pnorm)); reg2duns=dict(zip(tl["3. FRS ID"],tl["16. PARENT CO DB NUM"])); reg2fed=dict(zip(tl["3. FRS ID"],tl["21. FEDERAL FACILITY"].astype(str).str.upper()))
PUB=re.compile(r"(US DEP|DEPARTMENT OF|DEPT OF|ARMY|NAVY|AIR FORCE|US DOE|US EPA|CITY OF|COUNTY OF|TOWN OF|VILLAGE OF|STATE OF|COMMONWEALTH|UNIVERSITY|COLLEGE|REGENTS|SCHOOL|MUNICIPAL|SANITAT|WATER DISTRICT|SEWER|PUBLIC WORKS|METROPOLITAN|UTILITY DISTRICT|RECLAMATION|POSTAL|VETERANS|NASA|TENNESSEE VALLEY|BUREAU OF|GOVERNMENT|NATIONAL LABORATOR|PORT AUTHORITY|HOUSING AUTHORITY|TRANSIT|CORRECTION|PRISON|ARMY CORPS)")
m=pd.read_parquet(f"{OUT}/us_npdes_facility_master.parquet")
m["pnorm"]=m.REGISTRY_ID.map(reg2pn); m["duns"]=m.REGISTRY_ID.map(reg2duns); m["tri_fed"]=m.REGISTRY_ID.map(reg2fed)
m["is_private"]=((m.pnorm.notna())&(m.pnorm!="")&(~m.pnorm.fillna("").str.contains(PUB))&(m.is_potw==0)&(m.is_federal==0)&(~m.tri_fed.isin(["YES","Y"]))).astype(int)
p=m[m.is_private==1].dropna(subset=["REGISTRY_ID"]).copy(); p["naics2"]=p.naics.astype(str).str[:2]; p["is_major"]=pd.to_numeric(p.is_major,errors="coerce").fillna(0)
plant=p.sort_values("n_e90",ascending=False).drop_duplicates("REGISTRY_ID").set_index("REGISTRY_ID")[["pnorm","duns","naics2","is_major","STATE_CODE","n_insp","n_formal"]]
npdes2reg=dict(zip(m.NPDES_ID,m.REGISTRY_ID))
qn=pd.read_csv(f"{P1}/NPDES_QNCR_HISTORY.csv",dtype=str,encoding="latin-1")
qn["yr"]=pd.to_numeric(qn.YEARQTR.str[:4],errors="coerce"); qn["qq"]=pd.to_numeric(qn.YEARQTR.str[4:5],errors="coerce")
qn=qn[(qn.qq.between(1,4))&(qn.yr.between(2010,2025))]
for c in ["NUME90Q","NUMSVCD","NUMD8090Q"]: qn[c]=pd.to_numeric(qn[c],errors="coerce").fillna(0)
qn["REG"]=qn.NPDES_ID.map(npdes2reg); qn=qn.dropna(subset=["REG"]); qn=qn[qn.REG.isin(set(plant.index))]
py=qn.groupby(["REG","yr"]).agg(e90=("NUME90Q","sum"),ser=("NUMSVCD","sum"),rep=("NUMD8090Q","sum"),nq=("qq","size")).reset_index()
py=py.join(plant,on="REG")
# firm sizes / samples
fp=py.groupby("pnorm").REG.nunique(); fyrs=py.groupby("pnorm").yr.nunique()
py["nplant"]=py.pnorm.map(fp); py["fyrs"]=py.pnorm.map(fyrs)
py["any_e90"]=(py.e90>0).astype(float); py["le90"]=np.log1p(py.e90)
py.to_parquet(f"{OUT}/one_bad_plant_panel.parquet")
def dstrict(pn):
    d=py[py.pnorm.isin(pn)]; return d
SAMP={"A >=2 plants":py[py.nplant>=2],"B >=3 plants":py[py.nplant>=3],
      "C DUNS-confirmed":py[py.nplant>=3][py[py.nplant>=3].pnorm.map(py.groupby('pnorm').duns.nunique())==1],
      "D >=10 yrs":py[(py.nplant>=3)&(py.fyrs>=10)],"E majors":py[(py.nplant>=3)&(py.is_major==1)]}
print("plant-year rows:",len(py),"| plants:",py.REG.nunique(),"| firms:",py.pnorm.nunique())
for k,v in SAMP.items(): print(f"  sample {k}: firms={v.pnorm.nunique()} plants={v.REG.nunique()} rows={len(v)}")
PR=py[py.nplant>=3].copy()   # primary strict analysis frame

# ========== 3. NESTED VARIANCE PARTITION (SS shares) + bootstrap ==========
hdr("3. NESTED VARIANCE PARTITION firm / plant-within-firm / residual(time)")
def partition(df,col):
    y=df[col].astype(float).values; g=y.mean()
    fm=df.groupby("pnorm")[col].transform("mean").values
    pm=df.groupby("REG")[col].transform("mean").values
    ss_firm=((fm-g)**2).sum(); ss_plant=((pm-fm)**2).sum(); ss_res=((y-pm)**2).sum(); tot=ss_firm+ss_plant+ss_res
    return ss_firm/tot, ss_plant/tot, ss_res/tot
def boot_part(df,col,B=200):
    firms=df.pnorm.unique(); out=[]
    idxby=df.groupby("pnorm").indices
    for _ in range(B):
        pick=rng.choice(firms,len(firms),replace=True)
        rows=np.concatenate([idxby[f] for f in pick])
        out.append(partition(df.iloc[rows],col))
    a=np.array(out); return a.mean(0),np.percentile(a,[2.5,97.5],0)
rowsv=[]
for col,lab in [("e90","E90 count"),("any_e90","any-E90"),("le90","log(1+E90)"),("ser","serious count")]:
    f,pl_,r=partition(PR,col); mean,ci=boot_part(PR,col,B=150)
    print(f"  {lab:14} firm={f:.2f} plant={pl_:.2f} resid={r:.2f}  | plant/firm ratio={pl_/max(f,1e-9):.1f}  (plant CI[{ci[0][1]:.2f},{ci[1][1]:.2f}])")
    rowsv.append(dict(outcome=lab,firm=round(f,3),plant=round(pl_,3),resid=round(r,3),plant_ci_lo=round(ci[0][1],3),plant_ci_hi=round(ci[1][1],3),firm_ci_lo=round(ci[0][0],3),firm_ci_hi=round(ci[1][0],3)))
pd.DataFrame(rowsv).to_csv(f"{TAB}/formal_variance_decomposition.csv",index=False)
# MixedLM cross-check (log outcome, strict sample)
try:
    t0=time.time(); sub=PR.copy(); sub["pid"]=sub.REG.astype(str)
    md=smf.mixedlm("le90~C(yr)",sub,groups=sub["pnorm"],re_formula="1",vc_formula={"plant":"0+C(pid)"}).fit(method="lbfgs",maxiter=60)
    vf=md.cov_re.iloc[0,0]; vc=md.vcomp[0]; ve=md.scale; totv=vf+vc+ve
    print(f"  MixedLM(le90): firm var share={vf/totv:.2f} plant var share={vc/totv:.2f} resid={ve/totv:.2f}  (t={time.time()-t0:.0f}s)")
except Exception as e: print("  MixedLM failed/slow:",str(e)[:80])

# ========== 4. FIRM-SIZE GRADIENT ==========
hdr("4. FIRM-SIZE GRADIENT of firm vs plant variance (le90)")
grows=[]
for lo,hi,lab in [(3,4,"3-4"),(5,9,"5-9"),(10,19,"10-19"),(20,999,"20+")]:
    s=py[(py.nplant>=lo)&(py.nplant<=hi)]
    if s.pnorm.nunique()<15: continue
    f,pl_,r=partition(s,"le90"); grows.append(dict(size=lab,firms=s.pnorm.nunique(),firm=round(f,3),plant=round(pl_,3),resid=round(r,3)))
    print(f"  {lab:6} firms={s.pnorm.nunique():4} firm-share={f:.2f} plant-share={pl_:.2f}  -> {'PLANT' if pl_>f else 'FIRM'} dominant")
pd.DataFrame(grows).to_csv(f"{TAB}/firm_size_gradient.csv",index=False)
print("  gradient: firm-share rises with size ->",[ (g['size'],g['firm']) for g in grows])

# ========== 5. WORST-PLANT SHARE with bootstrap CI ==========
hdr("5. WORST-PLANT SHARE with bootstrap CIs")
def worst_shares(df):
    g=df.groupby(["pnorm","REG"]).e90.sum().reset_index()
    f=g.groupby("pnorm").e90.agg(["max","sum"]); f=f[f["sum"]>0]; return (f["max"]/f["sum"])
ws=worst_shares(PR)
def boot_med(x,B=300):
    x=np.asarray(x); return np.percentile([np.median(rng.choice(x,len(x),replace=True)) for _ in range(B)],[2.5,97.5])
ci=boot_med(ws.values); print(f"  overall median worst-share (K>=3) = {ws.median():.2f}  95% CI [{ci[0]:.2f},{ci[1]:.2f}]  (n={len(ws)})")

# ========== 6. RANDOM BENCHMARK with simulation interval + tail p ==========
hdr("6. RANDOM-ALLOCATION BENCHMARK (exposure-preserving) with 95% sim interval")
g=PR.groupby(["pnorm","REG"]).agg(e90=("e90","sum"),nq=("nq","sum")).reset_index()
def bench(sub,draws=300):
    firmv=sub.groupby("pnorm").agg(V=("e90","sum"),K=("REG","nunique"))
    firmv=firmv[(firmv.V>0)&(firmv.K>=3)]
    act=[]; sim_excess=[]
    for pn,r in firmv.iterrows():
        gg=sub[sub.pnorm==pn]; V=int(r.V); probs=gg.nq.values/gg.nq.sum()
        a=gg.e90.max()/V; al=rng.multinomial(V,probs,size=draws); sm=al.max(axis=1)/V
        act.append(a); sim_excess.append((a - sm.mean(), (a<=sm).mean()))
    return np.array(act),np.array(sim_excess)
rows=[]
for lo,hi,lab in [(3,4,"3-4"),(5,9,"5-9"),(10,19,"10-19"),(20,999,"20+")]:
    sub=g[g.pnorm.map(PR.groupby('pnorm').REG.nunique()).between(lo,hi)]
    if sub.pnorm.nunique()<15: continue
    a,se=bench(sub); tailp=se[:,1].mean()
    rows.append(dict(size=lab,firms=len(a),actual=round(np.mean(a),3),excess=round(np.mean(se[:,0]),3),share_p_lt05=round((se[:,1]<0.05).mean(),2)))
    print(f"  {lab:6} firms={len(a):4} actual={np.mean(a):.2f} mean-excess={np.mean(se[:,0]):+.2f}  share of firms with sim-p<0.05: {(se[:,1]<0.05).mean()*100:.0f}%")
pd.DataFrame(rows).to_csv(f"{TAB}/random_benchmark.csv",index=False)

# ========== 7-8. SURVIVAL of bad-plant status + transitions ==========
hdr("7-8. BAD-PLANT SURVIVAL (national top-k) + transition matrices")
py=py.sort_values(["REG","yr"])
def topk_flag(qthr):
    thr=py.groupby("yr").e90.transform(lambda s:s[s>0].quantile(qthr) if (s>0).any() else np.inf)
    return (py.e90>=thr)&(py.e90>0)
def km_survival(flag,maxh=6):
    d=py.assign(bad=flag.values); d["prev"]=d.groupby("REG").bad.shift(1).fillna(False)
    d["entry"]=d.bad&(~d.prev)
    surv=[1.0]; # discrete hazard: among currently-bad, P(exit next yr)
    for h in range(1,maxh):
        cur=d[d.bad]; nxt=d.assign(bad_next=d.groupby("REG").bad.shift(-1))
        pass
    # simpler: build spells
    spells=[]
    for reg,gg in d.groupby("REG"):
        b=gg.bad.values; i=0
        while i<len(b):
            if b[i]:
                j=i
                while j+1<len(b) and b[j+1]: j+=1
                cens=(j==len(b)-1)
                spells.append((j-i+1,cens)); i=j+1
            else: i+=1
    spells=pd.DataFrame(spells,columns=["dur","cens"])
    S=[]; n=len(spells)
    for h in range(1,maxh+1):
        atrisk=(spells.dur>=h).sum(); exit_=((spells.dur==h)&(~spells.cens)).sum()
        S.append(atrisk);
    # KM
    surv=1.0; km=[]
    for h in range(1,maxh+1):
        atrisk=(spells.dur>=h).sum(); exit_=((spells.dur==h)&(~spells.cens)).sum()
        haz=exit_/atrisk if atrisk>0 else 0; surv*= (1-haz); km.append(surv)
    return km,spells,len(spells)
survrows=[]; kms={}
for qthr,lab in [(0.90,"top10"),(0.95,"top5"),(0.99,"top1")]:
    flag=topk_flag(qthr); km,sp,ns=km_survival(flag); kms[lab]=km
    med=next((h for h,s in enumerate(km,1) if s<0.5),">6")
    survrows.append(dict(defn=lab,n_spells=ns,S1=round(km[0],2),S3=round(km[2],2),S5=round(km[4],2),median_dur=med))
    print(f"  {lab}: spells={ns} S(1)={km[0]:.2f} S(3)={km[2]:.2f} S(5)={km[4]:.2f} median dur={med}")
pd.DataFrame(survrows).to_csv(f"{TAB}/bad_plant_survival.csv",index=False)
# transition matrices 1/3/5 yr
bk=pd.cut(py.e90,[-1,0,2,9,29,1e9],labels=["clean","low","med","high","extreme"])
pyt=py.assign(bk=bk)
def transmat(h):
    d=pyt.copy(); d["bk_f"]=d.groupby("REG").bk.shift(-h); d["yr_f"]=d.groupby("REG").yr.shift(-h)
    d=d[d.yr_f==d.yr+h]; return pd.crosstab(d.bk,d.bk_f,normalize="index").round(2)
for h in [1,3,5]:
    tm=transmat(h); print(f"\n{h}-year transition (P):"); print(tm.to_string()); tm.to_csv(f"{TAB}/transition_matrix_{h}yr.csv")
tm1=transmat(1)

# ========== 9. STATE DEPENDENCE vs stable type ==========
hdr("9. STATE DEPENDENCE vs STABLE PLANT TYPE (any-E90)")
d=py.sort_values(["REG","yr"]).copy(); d["lag"]=d.groupby("REG").any_e90.shift(1)
d=d.dropna(subset=["lag"])
ar=smf.ols("any_e90~lag",d).fit(); print(f"  raw AR(1) coef on lagged any-E90 = {ar.params['lag']:.2f} (R2={ar.rsquared:.2f})")
# within-plant (plant demeaned) AR
d["lag_dm"]=d.lag-d.groupby("REG").lag.transform("mean"); d["y_dm"]=d.any_e90-d.groupby("REG").any_e90.transform("mean")
arw=smf.ols("y_dm~lag_dm",d).fit(); print(f"  within-plant AR coef = {arw.params['lag_dm']:.2f}  (persistence remaining after removing stable plant identity)")
# how much does plant identity alone predict?
pm=d.groupby("REG").any_e90.transform("mean"); r2_plant=1-((d.any_e90-pm)**2).sum()/((d.any_e90-d.any_e90.mean())**2).sum()
print(f"  R2 of stable plant identity (plant means) for any-E90 = {r2_plant:.2f}  -> {'plant TYPE dominates' if r2_plant>0.4 else 'mixed'}")

# ========== 10. SIBLING CORRELATION ==========
hdr("10. WITHIN-FIRM SIBLING CORRELATION (firms K>=3)")
# ICC from partition already; add P(sibling bad | focal bad) vs matched base
pp=PR.groupby("REG").agg(bad=("any_e90","max"),pnorm=("pnorm","first"),naics2=("naics2","first")).reset_index()
def sib_prob(df):
    num=den=0
    for pn,gg in df.groupby("pnorm"):
        b=gg.bad.values; k=len(b)
        for i in range(k):
            if b[i]==1:
                den+=1; num+= (b.sum()-1)/(k-1) if k>1 else 0
    return num/den if den else np.nan
psib=sib_prob(pp); base=pp.bad.mean()
print(f"  P(a sibling plant ever-bad | focal plant ever-bad) = {psib:.2f}  vs firm-level base rate {base:.2f}  lift {psib/base:.2f}x")
print(f"  ICC(firm) for any-E90 (from partition) = {rowsv[1]['firm']:.2f}  (high ICC=>firm-wide; low=>plant-specific)")
sibrows=[]
for lo,hi,lab in [(3,4,"3-4"),(5,9,"5-9"),(10,999,"10+")]:
    s=pp[pp.pnorm.map(PR.groupby('pnorm').REG.nunique()).between(lo,hi)]
    if s.pnorm.nunique()>=15: sibrows.append(dict(size=lab,psib=round(sib_prob(s),2),base=round(s.bad.mean(),2)))
pd.DataFrame(sibrows).to_csv(f"{TAB}/sibling_correlations.csv",index=False); print("  by size:",sibrows)

# ========== 11-12. INDUSTRY + MAJOR/MINOR ==========
hdr("11-12. INDUSTRY-SPECIFIC + MAJOR/MINOR")
irows=[]
for s in PR.naics2.value_counts().head(6).index:
    sub=PR[PR.naics2==s]; f,pl_,_=partition(sub,"le90"); wsh=worst_shares(sub).median()
    irows.append(dict(naics2=s,firms=sub.pnorm.nunique(),worst_share=round(wsh,2),firm_var=round(f,2),plant_var=round(pl_,2)))
    print(f"  NAICS {s}: worst={wsh:.2f} firm-var={f:.2f} plant-var={pl_:.2f} (firms={sub.pnorm.nunique()})")
pd.DataFrame(irows).to_csv(f"{TAB}/industry_results.csv",index=False)
for lab,sub in [("majors",PR[PR.is_major==1]),("non-majors",PR[PR.is_major==0])]:
    f,pl_,_=partition(sub,"le90"); print(f"  {lab}: worst={worst_shares(sub).median():.2f} firm-var={f:.2f} plant-var={pl_:.2f} (firms={sub.pnorm.nunique()})")

# ========== 13. LINKAGE ROBUSTNESS ==========
hdr("13. PARENT-LINKAGE ROBUSTNESS (median worst-share + plant-var)")
lrows=[]
for k,v in SAMP.items():
    if v.pnorm.nunique()<15: continue
    f,pl_,_=partition(v,"le90"); wsh=worst_shares(v).median()
    lrows.append(dict(sample=k,firms=v.pnorm.nunique(),worst_share=round(wsh,2),plant_var=round(pl_,2),firm_var=round(f,2)))
    print(f"  {k:18} worst={wsh:.2f} plant-var={pl_:.2f} firm-var={f:.2f}")
pd.DataFrame(lrows).to_csv(f"{TAB}/linkage_robustness.csv",index=False)

# ========== 14-15. ENFORCEMENT tail + escalation ==========
hdr("14-15. ENFORCEMENT tracks persistent bad plants? (secondary)")
persist=py.groupby("REG").agg(yrs_bad=("any_e90","sum"),e90=("e90","sum")).reset_index()
persist=persist.join(plant[["n_insp","n_formal","pnorm"]],on="REG")
persist["n_insp"]=pd.to_numeric(persist.n_insp,errors="coerce").fillna(0); persist["n_formal"]=pd.to_numeric(persist.n_formal,errors="coerce").fillna(0)
cf=pd.read_csv(f"{FEC}/CASE_FACILITIES.csv",usecols=["REGISTRY_ID"],dtype=str,encoding="latin-1"); fec=set(cf.REGISTRY_ID.dropna())
persist["in_fec"]=persist.REG.isin(fec)
persist["cat"]=pd.cut(persist.yrs_bad,[-1,0,2,5,99],labels=["never","1-2y","3-5y","6+y"])
et=persist.groupby("cat").agg(n=("REG","size"),insp=("n_insp","mean"),formal=("n_formal","mean"),fed=("in_fec","mean")).round(2)
print(et.to_string()); et.to_csv(f"{TAB}/enforcement_tail.csv")
print("  -> does enforcement rise with years-bad?",list(et.formal.round(2)))

# ========== FIGURES ==========
vdf=pd.DataFrame(rowsv); x=range(len(vdf))
plt.figure(figsize=(7,4)); plt.bar([i-.2 for i in x],vdf.firm,.4,label="firm (culture)",color="#72B7B2"); plt.bar([i+.2 for i in x],vdf.plant,.4,label="plant (one bad plant)",color="#4C78A8")
plt.xticks(list(x),vdf.outcome,fontsize=8,rotation=10); plt.ylabel("variance share"); plt.legend(); plt.title("Formal plant vs firm variance partition (K>=3, w/ 95% CI on plant)")
for i,rr in vdf.iterrows(): plt.plot([i+.2,i+.2],[rr.plant_ci_lo,rr.plant_ci_hi],color="k",lw=1)
plt.tight_layout(); plt.savefig(f"{FIG}/plant_vs_firm_variance_formal.png",dpi=120); plt.close()
gr=pd.DataFrame(grows); plt.figure(figsize=(7,4)); plt.plot(gr["size"],gr.firm,"o-",label="firm share"); plt.plot(gr["size"],gr.plant,"s-",label="plant share")
plt.xlabel("plants per firm"); plt.ylabel("variance share"); plt.legend(); plt.title("Firm-size gradient: plant vs firm variance"); plt.tight_layout(); plt.savefig(f"{FIG}/firm_size_gradient.png",dpi=120); plt.close()
rb=pd.DataFrame(rows); plt.figure(figsize=(7,4)); xx=range(len(rb)); plt.bar([i-.2 for i in xx],rb.actual,.4,label="actual",color="#4C78A8"); plt.bar([i+.2 for i in xx],rb.actual-rb.excess,.4,label="random(exposure)",color="#F58518")
plt.xticks(list(xx),rb["size"]); plt.ylabel("worst-plant share"); plt.legend(); plt.title("Actual vs simulated random concentration"); plt.tight_layout(); plt.savefig(f"{FIG}/actual_vs_random_concentration_simulation.png",dpi=120); plt.close()
plt.figure(figsize=(7,4))
for lab in ["top10","top5","top1"]: plt.plot(range(1,7),kms[lab],"o-",label=lab)
plt.ylim(0,1); plt.xlabel("years in bad status"); plt.ylabel("survival S(h)"); plt.legend(); plt.title("Bad-plant survival (national top-k by E90)"); plt.tight_layout(); plt.savefig(f"{FIG}/bad_plant_survival_curve.png",dpi=120); plt.close()
plt.figure(figsize=(5.5,4.5)); import numpy as _np; plt.imshow(tm1.values,cmap="Blues",vmin=0,vmax=1)
plt.xticks(range(5),tm1.columns,rotation=45,fontsize=8); plt.yticks(range(5),tm1.index,fontsize=8)
for i in range(5):
    for j in range(5): plt.text(j,i,f"{tm1.values[i,j]:.2f}",ha="center",va="center",fontsize=7)
plt.title("1-year transition matrix"); plt.tight_layout(); plt.savefig(f"{FIG}/transition_matrix_heatmap.png",dpi=120); plt.close()
sb=pd.DataFrame(sibrows); plt.figure(figsize=(6.5,4)); plt.bar([i-.2 for i in range(len(sb))],sb.psib,.4,label="P(sibling bad|focal bad)",color="#4C78A8"); plt.bar([i+.2 for i in range(len(sb))],sb.base,.4,label="base rate",color="#F58518")
plt.xticks(range(len(sb)),sb["size"]); plt.legend(fontsize=8); plt.ylabel("prob"); plt.title("Sibling correlation by firm size"); plt.tight_layout(); plt.savefig(f"{FIG}/sibling_correlation_by_firm_size.png",dpi=120); plt.close()
ir=pd.DataFrame(irows); plt.figure(figsize=(7,4)); plt.bar([i-.2 for i in range(len(ir))],ir.firm_var,.4,label="firm",color="#72B7B2"); plt.bar([i+.2 for i in range(len(ir))],ir.plant_var,.4,label="plant",color="#4C78A8")
plt.xticks(range(len(ir)),["N"+s for s in ir.naics2]); plt.legend(fontsize=8); plt.ylabel("variance share"); plt.title("Industry variance decomposition"); plt.tight_layout(); plt.savefig(f"{FIG}/industry_variance_decomposition.png",dpi=120); plt.close()
print("\nDONE. saved tables + 7 figures + panel.")
