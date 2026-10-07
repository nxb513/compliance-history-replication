"""Figures for the manuscript. Black and white: groups are told apart by line style and marker.
Sized for a 119 mm column, Arial lettering at 8 pt, saved as EPS and PDF (vector, fonts embedded) and PNG.
Fig. 1: pair correlations by distance and state border (E15). Fig. S1: sister plants around a first formal action (E14).
Fig. 2: out-of-sample AUC by source of history (B1)."""
import os
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
T = "./extension/tables" if os.path.isdir("./extension/tables") else "D:/eco_research/extension/tables"
O = os.path.dirname(os.path.abspath(__file__))
W = 119 / 25.4  # figure width in inches (119 mm)
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "legend.fontsize": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False, "ps.fonttype": 42, "pdf.fonttype": 42})


def save(fig, name):
    # shrink the canvas until the cropped figure (labels and legend included) is at most 119 mm wide
    for _ in range(5):
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        if bb.width <= W + 0.01: break
        fig.set_size_inches(fig.get_figwidth() - (bb.width - W), fig.get_figheight())
    for ext, kw in [("eps", {}), ("pdf", {}), ("png", {"dpi": 600})]:
        fig.savefig(f"{O}/{name}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)


# Fig. 1: pair correlations of long-run water violation rates by distance and state border (E15)
d = pd.read_csv(f"{T}/E15_pair_correlations_by_distance.csv"); d = d[(d.outcome == "water, raw") & (d.scope == "all pairs")]
bins = ["0-50", "50-100", "100-200", "200-400", "400-800", "800+"]; x = np.arange(len(bins))
fig, ax = plt.subplots(figsize=(W, 3.0))
spec = [("sibling", "same state", "Same firm, same state", "-", "o", "black"), ("sibling", "different state", "Same firm, different state", "--", "s", "black"),
        ("unrelated", "same state", "Different firms, same state", "-", "^", "0.45"), ("unrelated", "different state", "Different firms, different state", "--", "D", "0.45")]
for pair, st, lab, ls, mk, col in spec:
    s = d[(d.pair == pair) & (d.state == st)].set_index("distance_km").reindex(bins)
    y = s["corr"].where(s.pairs >= 100)
    ax.plot(x, y, ls=ls, marker=mk, color=col, mfc="white" if ls == "--" else col, lw=1.0, ms=4.5, label=lab)
ax.axhline(0, color="0.6", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(bins); ax.set_xlabel("Distance between the two plants (km)"); ax.set_ylabel("Correlation of long-run violation rates")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
save(fig, "Figure_1")

# Fig. 2: sister plants around a firm's first formal water enforcement action (E14)
e = pd.read_csv(f"{T}/E14_event_coefficients.csv")
fig, ax = plt.subplots(figsize=(W, 2.8))
for g, lab, off, mk, fc in [("siblings, same state", "Sister plants, same state", -0.09, "o", "black"), ("siblings, other state", "Sister plants, other states", 0.09, "s", "white")]:
    s = e[e.group == g].sort_values("event_time")
    ax.errorbar(s.event_time + off, s.coef, yerr=1.96 * s.se, fmt=mk, color="black", mfc=fc, ms=4.5, capsize=2, lw=0.8, label=lab)
ax.axhline(0, color="0.6", lw=0.5); ax.axvline(-0.5, color="0.6", lw=0.5, ls=":")
ax.set_xticks(range(-3, 4)); ax.set_xlabel("Years relative to the firm's first formal action at one of its plants")
ax.set_ylabel("Change in probability of an\neffluent violation (vs year -1)")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
save(fig, "Figure_S1")

# Fig. 3: out-of-sample AUC by source of history (B1)
b = pd.read_csv(f"{T}/B1_auc_single.csv")
outs = [("air", "any air violation", "Air,\nany violation"), ("air", "any high priority air violation", "Air,\nhigh priority violation"),
        ("rcra", "any RCRA violation", "Hazardous waste,\nany violation"), ("rcra", "any RCRA significant noncomplier", "Hazardous waste,\nsignificant noncomplier"),
        ("water", "any effluent violation", "Water,\nany effluent violation")]
src = [("own program history", "Plant's own record, same program", "o", "black", "black"), ("own water history", "Plant's own water record", "o", "black", "white"),
       ("same-state siblings' program history", "Sister plants in the same state, same program", "s", "0.45", "0.45"),
       ("other-state siblings' program history", "Sister plants in other states, same program", "s", "0.45", "white"),
       ("firm (all siblings') program history", "All sister plants (firm record), same program", "^", "0.45", "0.45")]
fig, ax = plt.subplots(figsize=(W, 3.6))
for i, (pg, oc, lab) in enumerate(outs):
    for j, (info, slab, mk, col, fc) in enumerate(src):
        info2 = info.replace("program history", "water history") if pg == "water" and "program" in info else info
        if pg == "water" and info == "own water history": continue
        r = b[(b.program == pg) & (b.outcome == oc) & (b.information == info2)]
        if len(r): ax.plot(r.auc.iloc[0], len(outs) - 1 - i + (j - 2) * 0.12, marker=mk, color=col, mfc=fc, mec=col, ms=5, ls="none", label=slab if i == 0 else None)
ax.axvline(0.5, color="0.6", lw=0.5, ls=":")
ax.set_yticks(range(len(outs))); ax.set_yticklabels([o[2] for o in outs][::-1]); ax.set_xlim(0.48, 0.86)
ax.set_xlabel("Out-of-sample AUC (history 2010 to 2017, outcome 2018 to 2025)")
fig.subplots_adjust(left=0.36, right=0.98, top=0.98, bottom=0.36)
fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=1)
save(fig, "Figure_2")
print("figures written to", O)
