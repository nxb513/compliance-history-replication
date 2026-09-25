"""E15: is the local firm component about the state regulator or about physical proximity? All plant pairs (about 26 million),
sibling vs unrelated, same vs different state, by distance. Firm bootstrap via firm-pair moment matrices."""
import sys, os, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
import statsmodels.formula.api as smf
pa = pd.read_parquet(f"{ROOT}/extension/plant_attributes_primary.parquet")
ex = pd.read_parquet(f"{ROOT}/extension/echo_exporter_primary.parquet").set_index("REGISTRY_ID")
pl = pa.copy(); pl["lat"] = pd.to_numeric(ex.FAC_LAT, errors="coerce").reindex(pl.index); pl["lon"] = pd.to_numeric(ex.FAC_LONG, errors="coerce").reindex(pl.index)
pl = pl.dropna(subset=["lat", "lon"])
obs = "l_params + l_outfalls + l_flow + flow_missing + permit_age + age_missing + major + epa_issued + comp_SWI + comp_PRE + C(n3)"
pl["resid_y"] = smf.ols(f"y ~ {obs}", pl).fit().resid
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet").groupby("REG").any_air.mean()
pl["air"] = air.reindex(pl.index)
pl = pl.reset_index()
N = len(pl); firm = pl.pnorm.astype("category").cat.codes.values.astype(np.int64); F = int(firm.max()) + 1
st = pl.st.values; n3 = pl.n3.astype(str).values
lat = np.radians(pl.lat.values); lon = np.radians(pl.lon.values)
z = {"raw": ((pl.y - pl.y.mean()) / pl.y.std()).values, "resid": ((pl.resid_y - pl.resid_y.mean()) / pl.resid_y.std()).values}
has_air = pl.air.notna().values; za = np.where(has_air, (pl.air - pl.air.mean()) / pl.air.std(), 0.0); zw = z["raw"]
edges6 = np.array([50, 100, 200, 400, 800]); edges4 = np.array([100, 200, 400])  # 6 bins for point estimates; 4 groups for bootstrap
print("plants", N, "firms", F, flush=True)

# point-estimate accumulators: cell = sib(2) x same(2) x bin6(6) x sameNAICS3(2) -> 48 cells; moments for symmetric corr: n, S1=sum(a+b), S2=sum(a^2+b^2), SAB=sum(ab)
NC = 48
acc = {k: np.zeros((NC, 4)) for k in ["raw", "resid", "cross"]}
# bootstrap accumulators: cell4 = sib x same x grp4(4) = 16, per firm pair (fmin, fmax): moments n, S1, S2, SAB (float64 sums, 16*F*F each)
BC = 16
bacc = {k: np.zeros((4, BC * F * F), dtype=np.float64) for k in ["raw", "cross"]}
t0 = time.time(); chunk_I, chunk_J = [], []
def flush(I, J):
    I = np.concatenate(I); J = np.concatenate(J)
    dlat = lat[J] - lat[I]; dlon = lon[J] - lon[I]
    h = np.sin(dlat / 2) ** 2 + np.cos(lat[I]) * np.cos(lat[J]) * np.sin(dlon / 2) ** 2
    d = 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))
    sib = (firm[I] == firm[J]).astype(np.int64); same = (st[I] == st[J]).astype(np.int64); sn3 = (n3[I] == n3[J]).astype(np.int64)
    b6 = np.searchsorted(edges6, d, side="right"); g4 = np.searchsorted(edges4, d, side="right")
    cell = ((sib * 2 + same) * 6 + b6) * 2 + sn3
    for k in ["raw", "resid"]:
        a, b = z[k][I], z[k][J]
        for m, w in enumerate([np.ones_like(a), a + b, a * a + b * b, a * b]):
            acc[k][:, m] += np.bincount(cell, weights=w, minlength=NC)
    # cross-program (water x air), symmetric over both orders, only pairs where both plants have air data
    both = has_air[I] & has_air[J]
    if both.any():
        c2 = cell[both]; wi, wj, ai, aj = zw[I][both], zw[J][both], za[I][both], za[J][both]
        # symmetric: pairs (w_i, a_j) and (w_j, a_i); moments on X = water, Y = air stacked
        for m, w in enumerate([2 * np.ones_like(wi), wi + wj, ai + aj, wi * aj + wj * ai]):
            acc["cross"][:, m] += np.bincount(c2, weights=w, minlength=NC)
    fmin = np.minimum(firm[I], firm[J]); fmax = np.maximum(firm[I], firm[J]); bcell = (sib * 2 + same) * 4 + g4
    idx = bcell * F * F + fmin * F + fmax
    a, b = z["raw"][I], z["raw"][J]
    for m, w in enumerate([np.ones_like(a), a + b, a * a + b * b, a * b]):
        bacc["raw"][m] += np.bincount(idx, weights=w, minlength=BC * F * F)
    if both.any():
        idx2 = idx[both]
        for m, w in enumerate([2 * np.ones_like(wi), wi + wj, ai + aj, wi * aj + wj * ai]):
            bacc["cross"][m] += np.bincount(idx2, weights=w, minlength=BC * F * F)
for i in range(N - 1):
    j = np.arange(i + 1, N); chunk_I.append(np.full(len(j), i)); chunk_J.append(j)
    if sum(len(x) for x in chunk_I) > 1_500_000 or i == N - 2:
        flush(chunk_I, chunk_J); chunk_I, chunk_J = [], []
print("pairs done in", round(time.time() - t0), "s", flush=True)

def corr_sym(M):  # M = [n, S1, S2, SAB] for symmetric pairs
    n, s1, s2, sab = M; m = s1 / (2 * n); var = s2 / (2 * n) - m ** 2; cov = sab / n - m ** 2
    return cov / var
def corr_cross(M):  # M = [2n, SX, SY, SXY] with X water, Y air stacked over both orders; var from global standardization approx: use moments
    n2, sx, sy, sxy = M; mx = sx / n2; my = sy / n2; return (sxy / n2 - mx * my)
bins6 = ["0-50", "50-100", "100-200", "200-400", "400-800", "800+"]
rows = []
for k in ["raw", "resid", "cross"]:
    A = acc[k]
    for sib in [1, 0]:
        for same in [1, 0]:
            for b in range(6):
                for scope, sn in [("all pairs", [0, 1]), ("different NAICS-3", [0])]:
                    cells = [((sib * 2 + same) * 6 + b) * 2 + s for s in sn]; M = A[cells].sum(axis=0)
                    if M[0] < 30: val = np.nan
                    else: val = corr_sym(M) if k != "cross" else corr_cross(M)
                    rows.append(dict(outcome={"raw": "water, raw", "resid": "water, residualized", "cross": "water x air (covariance of standardized rates)"}[k],
                                     pair="sibling" if sib else "unrelated", state="same state" if same else "different state", distance_km=bins6[b], scope=scope,
                                     pairs=int(M[0] if k != "cross" else M[0] / 2), corr=val))
E = pd.DataFrame(rows); E.to_csv(f"{OUT_T}/E15_pair_correlations_by_distance.csv", index=False)
pd.set_option("display.width", 250)
print(E[E.scope == "all pairs"].pivot_table(index=["outcome", "pair", "state"], columns="distance_km", values="corr", sort=False)[bins6].round(3).to_string())
print(E[E.scope == "all pairs"].pivot_table(index=["outcome", "pair", "state"], columns="distance_km", values="pairs", sort=False)[bins6].to_string())

# ---- bootstrap for pre-specified contrasts (raw water outcome and cross-program)
def stats_from(Bm, kind):
    # Bm: (4, 16) aggregated moments -> dict of contrasts
    f = corr_sym if kind == "raw" else corr_cross
    def c(sib, same, groups):
        cells = [(sib * 2 + same) * 4 + g for g in groups]; return f(Bm[:, cells].sum(axis=1))
    out = {}
    for sib, lab in [(1, "sib"), (0, "unrel")]:
        out[f"T1 {lab}: same minus different state, <=200 km"] = c(sib, 1, [0, 1]) - c(sib, 0, [0, 1])
        out[f"T1 {lab}: same minus different state, <=100 km"] = c(sib, 1, [0]) - c(sib, 0, [0])
    for same, lab in [(1, "same state"), (0, "different state")]:
        for gs, gl in [([0], "<=100 km"), ([1, 2], "100-400 km"), ([3], ">400 km")]:
            out[f"T2 sibling minus unrelated, {lab}, {gl}"] = c(1, same, gs) - c(0, same, gs)
    out["T3 sibling premium (all states) <=100 km minus >400 km"] = (c(1, 1, [0]) + 0) and ((f(Bm[:, [(1*2+s)*4+0 for s in [0,1]]].sum(axis=1)) - f(Bm[:, [(0*2+s)*4+0 for s in [0,1]]].sum(axis=1)))
                                                                  - (f(Bm[:, [(1*2+s)*4+3 for s in [0,1]]].sum(axis=1)) - f(Bm[:, [(0*2+s)*4+3 for s in [0,1]]].sum(axis=1))))
    out["T3 sibling premium, same state, <=100 km minus 100-400 km"] = (c(1, 1, [0]) - c(0, 1, [0])) - (c(1, 1, [1, 2]) - c(0, 1, [1, 2]))
    out["T3 sibling premium, different state, <=100 km minus >400 km"] = (c(1, 0, [0]) - c(0, 0, [0])) - (c(1, 0, [3]) - c(0, 0, [3]))
    out["level: sibling same state <=100 km"] = c(1, 1, [0]); out["level: sibling different state <=100 km"] = c(1, 0, [0])
    out["level: unrelated same state <=100 km"] = c(0, 1, [0]); out["level: unrelated different state <=100 km"] = c(0, 0, [0])
    out["level: sibling same state >400 km"] = c(1, 1, [3]); out["level: sibling different state >400 km"] = c(1, 0, [3])
    return out
rng = np.random.default_rng(161); res = []
for kind in ["raw", "cross"]:
    Bfull = bacc[kind].reshape(4, BC, F, F)
    est = stats_from(Bfull.sum(axis=(2, 3)), kind)
    draws = {k: [] for k in est}
    for b in range(200):
        w = np.bincount(rng.integers(0, F, F), minlength=F).astype(float)
        W = np.outer(w, w); np.fill_diagonal(W, w)  # siblings (same firm) weighted by w_f, unrelated firm pairs by w_f*w_g
        Bm = np.einsum("mcfg,fg->mc", Bfull, W, optimize=True)
        for k, v in stats_from(Bm, kind).items(): draws[k].append(v)
    for k in est:
        res.append(dict(outcome="water, raw" if kind == "raw" else "water x air", contrast=k, estimate=est[k], ci_low=np.nanpercentile(draws[k], 2.5), ci_high=np.nanpercentile(draws[k], 97.5)))
    print(kind, "bootstrap done", flush=True)
R = pd.DataFrame(res); R.to_csv(f"{OUT_T}/E15_contrasts_bootstrap.csv", index=False); print(R.round(3).to_string(index=False))
