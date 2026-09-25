"""Shared helpers for the Phase 4 extension (descriptive; primary sample)."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
ROOT = "."
OUT_T = f"{ROOT}/extension/tables"
OUT_F = f"{ROOT}/extension/figures"
RNG = np.random.default_rng(20260924)

def load_primary():
    py = pd.read_parquet(f"{ROOT}/data/interim/one_bad_plant_panel.parquet")
    pr = py[py.nplant >= 3].copy()
    pr["any_rep"] = (pr.rep > 0).astype(float)
    pr["lrep"] = np.log1p(pr.rep)
    return pr.reset_index(drop=True)

def partition(df, col, levels=("pnorm", "REG")):
    """Nested sum-of-squares shares. levels ordered from highest to lowest (last = plant).
    Returns dict of shares for each level plus residual."""
    y = df[col].astype(float).values
    prev = np.full(len(y), y.mean())
    out = {}; ss = []
    for lv in levels:
        key = df[list(lv)] if isinstance(lv, (list, tuple)) else df[lv]
        if isinstance(lv, (list, tuple)):
            m = df.groupby(list(lv))[col].transform("mean").values
        else:
            m = df.groupby(lv)[col].transform("mean").values
        ss.append(((m - prev) ** 2).sum()); prev = m
    ss.append(((y - prev) ** 2).sum())
    tot = sum(ss)
    names = [lv if isinstance(lv, str) else "x".join(lv) for lv in levels] + ["resid"]
    return {n: s / tot for n, s in zip(names, ss)}

def firm_bootstrap(df, stat_fn, B=300, seed=None):
    """Resample parent firms with replacement; stat_fn(df_boot) -> dict or float."""
    rng = np.random.default_rng(seed) if seed is not None else RNG
    firms = df.pnorm.unique(); idx = df.groupby("pnorm").indices
    res = []
    for b in range(B):
        pick = rng.choice(firms, len(firms), replace=True)
        rows = np.concatenate([idx[f] for f in pick])
        d = df.iloc[rows].copy()
        # relabel duplicated firms/plants so resampled copies are distinct clusters
        rep_id = np.concatenate([np.full(len(idx[f]), k) for k, f in enumerate(pick)])
        d["pnorm"] = d["pnorm"].astype(str) + "_" + rep_id.astype(str)
        d["REG"] = d["REG"].astype(str) + "_" + rep_id.astype(str)
        res.append(stat_fn(d))
    return res

def ci(a, lo=2.5, hi=97.5):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return (np.percentile(a, lo), np.percentile(a, hi)) if len(a) else (np.nan, np.nan)

def reml_shares(df, col):
    import statsmodels.formula.api as smf
    d = df[[col, "pnorm", "REG", "yr"]].dropna().copy()
    d["pid"] = d.REG.astype(str)
    # 2026-09-25: lbfgs alone can stop at a lower optimum; use the best of three optimizers (see reml_best)
    r = reml_best(df, col)
    return {"firm": r["firm"], "plant": r["plant"], "resid": r["resid"], "converged": True, "method": r["method"]}

def reml_best(df, col):
    """REML firm/plant shares, best log-likelihood across lbfgs, Powell and Nelder-Mead (lbfgs alone can stop at a lower optimum)."""
    import statsmodels.formula.api as smf
    d = df[[col, "pnorm", "REG", "yr"]].dropna().copy(); d["pid"] = d.REG.astype(str)
    fits = []
    for meth in ["lbfgs", "powell", "nm"]:
        try:
            md = smf.mixedlm(f"{col} ~ C(yr)", d, groups=d["pnorm"], re_formula="1", vc_formula={"plant": "0 + C(pid)"}).fit(method=meth, maxiter=2000, reml=True)
            fits.append((md.llf, meth, md))
        except Exception:
            pass
    llf, meth, md = max(fits, key=lambda t: t[0])
    vf = float(md.cov_re.iloc[0, 0]); vp = float(md.vcomp[0]); ve = float(md.scale); t = vf + vp + ve
    return {"firm": vf / t, "plant": vp / t, "resid": ve / t, "method": meth, "llf": llf, "llf_spread": max(f[0] for f in fits) - min(f[0] for f in fits)}

def twoway_demean(d, cols, g1, g2, iters=50):
    """Residualize cols on two sets of fixed effects (alternating projections)."""
    out = d[cols].astype(float).copy()
    for _ in range(iters):
        for g in (g1, g2):
            out = out - out.groupby(d[g]).transform("mean")
    return out
