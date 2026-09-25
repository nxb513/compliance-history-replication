"""B5: Holm and Benjamini-Hochberg adjusted p-values for own-chronic (F1) and sibling (F2) coefficients in E9, R2, B2, and for the E14 sibling post-event tests (F3)."""
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
T = "./extension/tables"
rows = []
for _, x in pd.read_csv(f"{T}/E9_air_lpm.csv").iterrows(): rows.append(("E9", f"air, {x['sample']}, ever violated", x.own_chronic, x.own_se, x.sibling_chronic, x.sib_se))
for _, x in pd.read_csv(f"{T}/E9_rcra_lpm.csv").iterrows(): rows.append(("E9", f"rcra, {x['sample']}, ever violated", x.own_chronic, x.own_se, x.sibling_chronic, x.sib_se))
for _, x in pd.read_csv(f"{T}/R2_E9_chronic_thresholds.csv").iterrows(): rows.append(("R2", f"{x.program}, {x.outcome}, {x.definition}", x.own, x.own_se, x.sibling, x.sib_se))
for _, x in pd.read_csv(f"{T}/B2_detection_source.csv").iterrows(): rows.append(("B2", f"{x.program}, detected by {x.detected_by}", x.own, x.own_se, x.sibling, x.sib_se))
d = pd.DataFrame(rows, columns=["table", "test", "own", "own_se", "sibling", "sib_se"])
out = []
for fam, b, se in [("F1 own chronic", "own", "own_se"), ("F2 sibling of chronic", "sibling", "sib_se")]:
    p = 2 * (1 - norm.cdf(np.abs(d[b] / d[se])))
    out.append(pd.DataFrame(dict(family=fam, table=d.table, test=d.test, coef=d[b], se=d[se], p=p, p_holm=multipletests(p, method="holm")[1], p_bh=multipletests(p, method="fdr_bh")[1])))
e14 = pd.read_csv(f"{T}/E14_event_tests.csv"); e14 = e14[e14.group != "enforced plant"]
out.append(pd.DataFrame(dict(family="F3 E14 sibling post-event", table="E14", test=e14.group, coef=e14.mean_post, se=np.nan, p=e14.post_joint_p,
                             p_holm=multipletests(e14.post_joint_p, method="holm")[1], p_bh=multipletests(e14.post_joint_p, method="fdr_bh")[1])))
out = pd.concat(out, ignore_index=True)
pd.set_option("display.width", 250); print(out.round(4).to_string(index=False))
print(out.groupby("family").agg(tests=("p", "size"), raw_p_below_05=("p", lambda s: (s < 0.05).sum()), holm_below_05=("p_holm", lambda s: (s < 0.05).sum()), bh_below_05=("p_bh", lambda s: (s < 0.05).sum())))
out.to_csv(f"{T}/B5_multiple_testing.csv", index=False)
