"""Figures for the EPG manuscript. Black-and-white legible (no tints): groups are told apart by line style and marker."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
T = "./extension/tables"; O = "./figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False})

# Figure 1: pair correlations of long-run water violation rates by distance and state border (E15)
d = pd.read_csv(f"{T}/E15_pair_correlations_by_distance.csv"); d = d[(d.outcome == "water, raw") & (d.scope == "all pairs")]
bins = ["0-50", "50-100", "100-200", "200-400", "400-800", "800+"]; x = np.arange(len(bins))
fig, ax = plt.subplots(figsize=(6.3, 3.6))
spec = [("sibling", "same state", "Same firm, same state", "-", "o", "black"), ("sibling", "different state", "Same firm, different state", "--", "s", "black"),
        ("unrelated", "same state", "Different firms, same state", "-", "^", "0.45"), ("unrelated", "different state", "Different firms, different state", "--", "D", "0.45")]
for pair, st, lab, ls, mk, col in spec:
    s = d[(d.pair == pair) & (d.state == st)].set_index("distance_km").reindex(bins)
    y = s["corr"].where(s.pairs >= 100)
    ax.plot(x, y, ls=ls, marker=mk, color=col, mfc="white" if ls == "--" else col, lw=1.4, ms=5.5, label=lab)
ax.axhline(0, color="0.6", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels([b + " km" for b in bins]); ax.set_xlabel("Distance between the two plants"); ax.set_ylabel("Correlation of long-run violation rates")
ax.legend(frameon=False, fontsize=8, loc="upper right"); fig.tight_layout(); fig.savefig(f"{O}/Figure1.png", dpi=600); fig.savefig(f"{O}/Figure1.tif", dpi=600); plt.close(fig)

# Figure 2: out-of-sample AUC by source of history (B1)
b = pd.read_csv(f"{T}/B1_auc_single.csv")
outs = [("air", "any air violation", "Air: any violation"), ("air", "any high priority air violation", "Air: high priority violation"),
        ("rcra", "any RCRA violation", "Hazardous waste: any violation"), ("rcra", "any RCRA significant noncomplier", "Hazardous waste: significant noncomplier"),
        ("water", "any effluent violation", "Water: any effluent violation")]
src = [("own program history", "Plant's own record, same program", "o", "black", "black"), ("own water history", "Plant's own water record", "o", "black", "white"),
       ("same-state siblings' program history", "Sister plants in the same state, same program", "s", "0.45", "0.45"),
       ("other-state siblings' program history", "Sister plants in other states, same program", "s", "0.45", "white"),
       ("firm (all siblings') program history", "All sister plants (firm record), same program", "^", "0.45", "0.45")]
fig, ax = plt.subplots(figsize=(6.3, 3.9))
for i, (pg, oc, lab) in enumerate(outs):
    for j, (info, slab, mk, col, fc) in enumerate(src):
        info2 = info.replace("program history", "water history") if pg == "water" and "program" in info else info
        if pg == "water" and info == "own water history": continue
        r = b[(b.program == pg) & (b.outcome == oc) & (b.information == info2)]
        if len(r): ax.plot(r.auc.iloc[0], len(outs) - 1 - i + (j - 2) * 0.12, marker=mk, color=col, mfc=fc, mec=col, ms=6, ls="none", label=slab if i == 0 else None)
ax.axvline(0.5, color="0.6", lw=0.6, ls=":")
ax.set_yticks(range(len(outs))); ax.set_yticklabels([o[2] for o in outs][::-1]); ax.set_xlim(0.48, 0.86)
ax.set_xlabel("Out-of-sample AUC (history 2010 to 2017, outcome 2018 to 2025)")
ax.legend(frameon=False, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.45, -0.52), ncol=2)
fig.tight_layout(); fig.savefig(f"{O}/Figure2.png", dpi=600, bbox_inches="tight"); fig.savefig(f"{O}/Figure2.tif", dpi=600, bbox_inches="tight"); plt.close(fig)

# Figure 3: sister plants around a firm's first formal water enforcement action (E14)
e = pd.read_csv(f"{T}/E14_event_coefficients.csv")
fig, ax = plt.subplots(figsize=(6.3, 3.3))
for g, lab, off, mk, fc in [("siblings, same state", "Sister plants, same state", -0.09, "o", "black"), ("siblings, other state", "Sister plants, other states", 0.09, "s", "white")]:
    s = e[e.group == g].sort_values("event_time")
    ax.errorbar(s.event_time + off, s.coef, yerr=1.96 * s.se, fmt=mk, color="black", mfc=fc, ms=5.5, capsize=2.5, lw=1, label=lab)
ax.axhline(0, color="0.6", lw=0.6); ax.axvline(-0.5, color="0.6", lw=0.6, ls=":")
ax.set_xticks(range(-3, 4)); ax.set_xlabel("Years relative to the firm's first formal action at one of its plants")
ax.set_ylabel("Change in probability of an\neffluent violation (vs year -1)")
ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=2); fig.tight_layout()
fig.savefig(f"{O}/Figure3.png", dpi=600, bbox_inches="tight"); fig.savefig(f"{O}/Figure3.tif", dpi=600, bbox_inches="tight"); plt.close(fig)
print("figures written")
