"""E23: major county floods (NOAA Storm Events) as exogenous plant-level shocks. Stacked event study for hit plants and for unaffected
out-of-state sister plants. Pre-specified in Addendum 7."""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.api as sm
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 10e6
# ---- storm events: county-year flood damage
fr = []
for f in sorted(glob.glob(f"{ROOT}/data/raw/noaa_storms/StormEvents_details-*.csv.gz")):
    d = pd.read_csv(f, dtype=str, usecols=["STATE_FIPS", "CZ_TYPE", "CZ_FIPS", "EVENT_TYPE", "YEAR", "DAMAGE_PROPERTY"])
    fr.append(d[(d.CZ_TYPE == "C") & d.EVENT_TYPE.isin(["Flood", "Flash Flood"])])
se = pd.concat(fr, ignore_index=True)
def dollars(s):
    s = str(s).strip().upper()
    if s in ("", "NAN"): return 0.0
    mult = {"K": 1e3, "M": 1e6, "B": 1e9}.get(s[-1], 1.0); num = s[:-1] if s[-1] in "KMB" else s
    try: return float(num) * mult
    except ValueError: return 0.0
se["dmg"] = se.DAMAGE_PROPERTY.map(dollars); se["fips"] = se.STATE_FIPS.str.zfill(2) + se.CZ_FIPS.str.zfill(3); se["yr"] = se.YEAR.astype(int)
cy = se.groupby(["fips", "yr"]).dmg.sum()
print("flood events (county-based):", len(se), "| county-years with damage >= threshold:", int((cy >= THR).sum()), flush=True)
# ---- plants: county
ex = pd.read_parquet(f"{ROOT}/extension/echo_exporter_primary.parquet").set_index("REGISTRY_ID")
fips = ex.FAC_DERIVED_STCTY_FIPS.astype(str).str.zfill(5)
pr = load_primary(); pr["fips"] = pr.REG.map(fips)
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet")[["REG", "yr", "any_air"]]; rc = pd.read_parquet(f"{ROOT}/extension/rcra_panel.parquet")[["REG", "yr", "any_rcra"]]
P = pr.merge(air, on=["REG", "yr"], how="left").merge(rc, on=["REG", "yr"], how="left")
dmg = lambda fp, y: cy.get((fp, y), 0.0)
plants = P.drop_duplicates("REG").set_index("REG")[["pnorm", "STATE_CODE", "fips"]]
major = {(f, y) for (f, y), v in cy.items() if v >= THR}; any1m = {(f, y) for (f, y), v in cy.items() if v >= 1e6}
hit_year = {}
for r, row in plants.iterrows():
    for e in range(2011, 2023):
        if (row.fips, e) in major and not any((row.fips, e - k) in major for k in (1, 2, 3)):
            hit_year[r] = e; break
hy = pd.Series(hit_year); print("hit plants:", len(hy), "| firms:", plants.loc[hy.index].pnorm.nunique(), "| by year:", hy.value_counts().sort_index().to_dict(), flush=True)
firm_hits = {}  # (firm, year) -> set of states of hit plants
for r, e in hy.items(): firm_hits.setdefault((plants.at[r, "pnorm"], e), set()).add(plants.at[r, "STATE_CODE"])
firm_hit_years = {}
for (f, e) in firm_hits: firm_hit_years.setdefault(f, set()).add(e)
clean_cty = lambda fp, e: not any((fp, e + k) in any1m for k in range(-3, 4))
stacks = []
for e in sorted(set(hy.values)):
    win = P[P.yr.between(e - 3, e + 3)].copy()
    hit = set(hy[hy == e].index)
    firms_e = {f for (f, yy) in firm_hits if yy == e}
    role = []
    for r, f, s, fp in zip(win.REG, win.pnorm, win.STATE_CODE, win.fips):
        if r in hit: role.append("hit plant")
        elif f in firms_e and s not in firm_hits[(f, e)] and clean_cty(fp, e) and r not in hy.index: role.append("unaffected sister plant")
        elif not (set(range(e - 3, e + 4)) & firm_hit_years.get(f, set())) and clean_cty(fp, e): role.append("control")
        else: role.append("drop")
    win["role"] = role; win = win[win.role != "drop"]; win["stk"] = e; win["k"] = win.yr - e; stacks.append(win)
S = pd.concat(stacks, ignore_index=True); S["ps"] = S.REG + "_" + S.stk.astype(str); S["ys"] = S.yr.astype(str) + "_" + S.stk.astype(str)
print("plant-stacks by role:", S.drop_duplicates("ps").role.value_counts().to_dict(), flush=True)
out, tests = [], []
for g in ["hit plant", "unaffected sister plant"]:
    for y in ["any_e90", "any_air", "any_rcra"]:
        D = S[S.role.isin([g, "control"])].dropna(subset=[y]).copy(); D["tr"] = (D.role == g).astype(float); xs = []
        for k in [-3, -2, 0, 1, 2, 3]: D[f"k{k}"] = ((D.k == k) & (D.tr == 1)).astype(float); xs.append(f"k{k}")
        R = twoway_demean(D, [y] + xs, "ps", "ys")
        f = sm.OLS(R[y], R[xs]).fit(cov_type="cluster", cov_kwds={"groups": D.pnorm.astype("category").cat.codes})
        pre = f.wald_test("k-3 = 0, k-2 = 0", scalar=True)
        w = np.array([1 / 3, 1 / 3, 1 / 3]); b = f.params[["k0", "k1", "k2"]].values; V = f.cov_params().loc[["k0", "k1", "k2"], ["k0", "k1", "k2"]].values
        post, post_se = float(w @ b), float(np.sqrt(w @ V @ w))
        for k in xs: out.append(dict(group=g, outcome=y, event_time=int(k[1:]), coef=f.params[k], se=f.bse[k]))
        tests.append(dict(threshold=THR, group=g, outcome=y, treated_plant_stacks=D[D.tr == 1].ps.nunique(), treated_firms=D[D.tr == 1].pnorm.nunique(), control_plant_stacks=D[D.tr == 0].ps.nunique(),
                          mean_at_minus1=D[(D.tr == 1) & (D.k == -1)][y].mean(), post_avg_0_2=post, post_se=post_se, pretrend_p=float(pre.pvalue)))
        print(tests[-1], flush=True)
out, tests = pd.DataFrame(out), pd.DataFrame(tests)
tag = f"{int(THR/1e6)}M"
out.to_csv(f"{OUT_T}/E23_flood_event_coefficients_{tag}.csv", index=False); tests.to_csv(f"{OUT_T}/E23_flood_tests_{tag}.csv", index=False)
pd.set_option("display.width", 250); print(tests.round(4).to_string(index=False))
