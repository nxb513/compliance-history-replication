"""E22: Florida HPV-penalty shock (Blundell 2020) as a causal test. (a) air violations at FL plants (validation); (b) water violations at the same FL plants
(cross-program effect); (c) water and air violations at out-of-state sister plants of firms with a FL plant (firm channel)."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
YRS = range(2010, 2018); CTRL = ["AL", "SC", "GA"]
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet"); air = air[air.yr.between(2010, 2017)].copy()
pr = load_primary(); pr = pr[pr.yr.between(2010, 2017)].copy()

def event(D, y, tr, fe=("REG", "yr"), extra_fe=None):
    D = D.copy(); xs = []
    for t in YRS:
        if t == 2011: continue
        D[f"e{t}"] = ((D.yr == t) & (D[tr] == 1)).astype(float); xs.append(f"e{t}")
    D["post"] = ((D.yr >= 2013) & (D[tr] == 1)).astype(float); D["trans"] = ((D.yr == 2012) & (D[tr] == 1)).astype(float); D["pre10"] = ((D.yr == 2010) & (D[tr] == 1)).astype(float)
    g2 = extra_fe if extra_fe else fe[1]
    R = twoway_demean(D, [y] + xs + ["post", "trans", "pre10"], fe[0], g2)
    cl = dict(cov_type="cluster", cov_kwds={"groups": D.pnorm.astype("category").cat.codes})
    fe_ = sm.OLS(R[y], R[xs]).fit(**cl); fp = sm.OLS(R[y], R[["post", "trans", "pre10"]]).fit(**cl)
    return fe_, fp

out, ev = [], []
# (a) and (b): FL vs AL/SC/GA air-regulated plants
S = air[air.STATE_CODE.isin(["FL"] + CTRL)].copy(); S["tr"] = (S.STATE_CODE == "FL").astype(int)
for y, lab in [("any_air", "E22a air violation (validation)"), ("any_e90", "E22b water violation, same plants")]:
    fe_, fp = event(S, y, "tr")
    out.append(dict(test=lab, plants_treated=S[S.tr == 1].REG.nunique(), plants_control=S[S.tr == 0].REG.nunique(), mean_treated_pre=S[(S.tr == 1) & (S.yr <= 2011)][y].mean(),
                    post=fp.params.post, post_se=fp.bse.post, pre2010=fp.params.pre10, pre2010_se=fp.bse.pre10, trans2012=fp.params.trans))
    for k in fe_.params.index: ev.append(dict(test=lab, year=int(k[1:]), coef=fe_.params[k], se=fe_.bse[k]))
    # placebo states
    pl = []
    for s, n in air[~air.STATE_CODE.isin(["FL"] + CTRL)].groupby("STATE_CODE").REG.nunique().items():
        if n < 30: continue
        Q = air[air.STATE_CODE.isin([s] + CTRL)].copy(); Q["tr"] = (Q.STATE_CODE == s).astype(int)
        pl.append(event(Q, y, "tr")[1].params.post)
    pl = np.array(pl); out[-1].update(placebo_states=len(pl), placebo_share_below=float((pl <= fp.params.post).mean()), placebo_share_abs_ge=float((np.abs(pl) >= abs(fp.params.post)).mean()))
    print(out[-1], flush=True)
# (c) firm channel: non-FL plants; treated if firm had an air-regulated FL plant in 2010-2012
flf = set(air[(air.STATE_CODE == "FL") & (air.yr <= 2012)].pnorm)
ctrlf = set(air[(air.STATE_CODE.isin(CTRL)) & (air.yr <= 2012)].pnorm) - flf
for base, y, lab in [(pr, "any_e90", "E22c water violation, out-of-state sister plants"), (air, "any_air", "E22c air violation, out-of-state sister plants")]:
    W = base[base.STATE_CODE != "FL"].copy(); W["tr"] = W.pnorm.isin(flf).astype(int); W["sy"] = W.STATE_CODE.astype(str) + "_" + W.yr.astype(str)
    fe_, fp = event(W, y, "tr", extra_fe="sy")
    out.append(dict(test=lab, plants_treated=W[W.tr == 1].REG.nunique(), plants_control=W[W.tr == 0].REG.nunique(), firms_treated=W[W.tr == 1].pnorm.nunique(),
                    mean_treated_pre=W[(W.tr == 1) & (W.yr <= 2011)][y].mean(), post=fp.params.post, post_se=fp.bse.post, pre2010=fp.params.pre10, pre2010_se=fp.bse.pre10, trans2012=fp.params.trans))
    for k in fe_.params.index: ev.append(dict(test=lab, year=int(k[1:]), coef=fe_.params[k], se=fe_.bse[k]))
    print(out[-1], flush=True)
    # placebo: firms with an AL/SC/GA plant but no FL plant, plants outside FL/AL/SC/GA
    P = base[~base.STATE_CODE.isin(["FL"] + CTRL)].copy(); P = P[~P.pnorm.isin(flf)]; P["tr"] = P.pnorm.isin(ctrlf).astype(int); P["sy"] = P.STATE_CODE.astype(str) + "_" + P.yr.astype(str)
    fe_p, fpp = event(P, y, "tr", extra_fe="sy")
    out.append(dict(test=lab.replace("E22c", "E22c placebo (firms with AL/SC/GA plant, no FL plant)"), plants_treated=P[P.tr == 1].REG.nunique(), plants_control=P[P.tr == 0].REG.nunique(),
                    firms_treated=P[P.tr == 1].pnorm.nunique(), post=fpp.params.post, post_se=fpp.bse.post, pre2010=fpp.params.pre10, pre2010_se=fpp.bse.pre10))
    print(out[-1], flush=True)
out, ev = pd.DataFrame(out), pd.DataFrame(ev)
pd.set_option("display.width", 250); print(out.round(4).to_string(index=False)); print(ev.round(4).to_string(index=False))
out.to_csv(f"{OUT_T}/E22_florida_did.csv", index=False); ev.to_csv(f"{OUT_T}/E22_florida_event.csv", index=False)
