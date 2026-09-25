"""R4 (Addendum 9): REML plant/firm shares on the DUNS-confirmed subset (sample C of the original design) vs the primary sample."""
import sys, os, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from ext_common import *
P = load_primary()
dc = P.groupby("pnorm").duns.nunique(); strict = set(dc[dc == 1].index)
air = pd.read_parquet(f"{ROOT}/extension/air_panel.parquet"); rc = pd.read_parquet(f"{ROOT}/extension/rcra_panel.parquet")
out = []
for prog, D, y in [("water", P, "any_e90"), ("air", air, "any_air"), ("rcra", rc, "any_rcra")]:
    for samp, S in [("primary", D), ("DUNS-confirmed", D[D.pnorm.isin(strict)])]:
        r = reml_best(S, y)
        out.append(dict(program=prog, sample=samp, firms=S.pnorm.nunique(), plants=S.REG.nunique(), mean=S[y].mean(), firm=r["firm"], plant=r["plant"], resid=r["resid"], method=r["method"], llf_spread=r.get("llf_spread")))
        print(out[-1], flush=True)
R = pd.DataFrame(out); R.to_csv(f"{OUT_T}/R4_linkage_strict_reml.csv", index=False)
pd.set_option("display.width", 200); print(R.round(3).to_string(index=False))
