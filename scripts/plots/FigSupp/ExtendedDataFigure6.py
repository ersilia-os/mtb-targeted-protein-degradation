"""
Extended Data Figure 6: on-target/off-target selectivity and docking-vs-Boltz-2 method
agreement across the final catalytic (script 99) and non-catalytic (script 100) hit selections.

Promotes the two exploratory scripts drafted in tmp/ (supp_fig_selectivity.py,
supp_fig_selection_funnel.py) into a single 2x2 manuscript-style supplementary figure.

(a) All 1,095 filtered hits (script 70): best on-target Mtb docking score (min across the 4 CAT
    + 8 NON-CAT curated pockets) vs. best off-target human docking score (389-pocket AF2
    counter-screen, scripts 90-97). The 244 prioritized compounds (script 101) are highlighted.
(b) Same axes, restricted to the 244 prioritized compounds, colored by selection origin
    (catalytic-only / non-catalytic-only / both, script 101's `origins` column).
(c-d) Docking-vs-Boltz-2 agreement within each hit_type tier (single/dual/multi) of the
    catalytic (c) and non-catalytic (d) selections: each (targets, compound_id) top-N slot is
    classified by which method(s) ranked it there.

Usage:
    python ExtendedDataFigure6.py
"""
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "..", "src"))

import matplotlib.pyplot as plt
import pandas as pd
import stylia
from stylia.config import get_fg_color

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

output_dir = os.path.join(root, "..", "..", "..", "output")
plots_dir = os.path.join(output_dir, "plots", "FigSupp", "ExtendedDataFigure6")
os.makedirs(plots_dir, exist_ok=True)

AUDIT_CSV = os.path.join(output_dir, "101_merge_selections", "audit_input.csv")
CAT_CSV = os.path.join(output_dir, "99_catalytic_selection", "catalytic_hits.csv")
NONCAT_CSV = os.path.join(output_dir, "100_noncatalytic_selection", "noncatalytic_hits.csv")

CURATED_DOCKING_COLS = [
    "docking_alaS_CAT", "docking_alaS_NONCAT",
    "docking_aspS_CAT", "docking_aspS_NONCAT",
    "docking_lysS_CAT", "docking_lysS_NONCAT",
    "docking_pheST_CAT", "docking_pheST_NONCAT",
]
ORIGIN_LABELS = {
    "catalytic": "Catalytic only",
    "non-catalytic": "Non-catalytic only",
    "catalytic|non-catalytic": "Both",
}
HIT_TYPE_ORDER = ["single", "dual", "multi"]
AGREEMENT_ORDER = ["docking only", "boltz2 only", "both methods"]
AGREEMENT_LABELS = {"docking only": "Docking only", "boltz2 only": "Boltz-2 only", "both methods": "Both methods"}

# ===========================================================================
# Data
# ===========================================================================
audit = pd.read_csv(AUDIT_CSV)
audit["best_on_target"] = audit[CURATED_DOCKING_COLS].min(axis=1)
audit["origin_label"] = audit["origins"].map(ORIGIN_LABELS)

selectivity_data = audit[["compound_id", "prioritized", "origin_label", "best_on_target", "human_best_af2"]]
selectivity_data_path = os.path.join(plots_dir, "ExtendedDataFigure6_selectivity_data.csv")
selectivity_data.to_csv(selectivity_data_path, index=False)
print(f"Saved {len(selectivity_data)} row(s) to {selectivity_data_path}")


def slot_agreement_counts(df):
    """Per hit_type, classify each (targets, compound_id) selection slot by which method(s)
    ranked it top-N, and return counts as a hit_type x agreement-category table."""
    methods_per_slot = df.groupby(["hit_type", "targets", "compound_id"])["method"].apply(set)

    def classify(methods):
        if methods == {"docking", "boltz2"}:
            return "both methods"
        return f"{next(iter(methods))} only"

    labels = methods_per_slot.apply(classify)
    counts = labels.reset_index(name="agreement").groupby(["hit_type", "agreement"]).size().unstack(fill_value=0)
    for col in AGREEMENT_ORDER:
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[AGREEMENT_ORDER]
    present = [h for h in HIT_TYPE_ORDER if h in counts.index]
    return counts.loc[present]


cat_counts = slot_agreement_counts(pd.read_csv(CAT_CSV))
noncat_counts = slot_agreement_counts(pd.read_csv(NONCAT_CSV))
funnel_data = pd.concat([cat_counts.assign(pocket_type="catalytic"), noncat_counts.assign(pocket_type="non-catalytic")])
funnel_data_path = os.path.join(plots_dir, "ExtendedDataFigure6_funnel_data.csv")
funnel_data.to_csv(funnel_data_path)
print(f"Saved {len(funnel_data)} row(s) to {funnel_data_path}")

# ===========================================================================
# Plot
# ===========================================================================
PANEL_LABEL_Y_PAD = 0.015


def add_panel_label(fig, ax, letter):
    renderer = fig.canvas.get_renderer()
    bbox = ax.get_tightbbox(renderer)
    x0 = fig.transFigure.inverted().transform((bbox.x0, 0))[0]
    top_y = fig.transFigure.inverted().transform((0, bbox.y1))[1]
    fig.text(x0, top_y + PANEL_LABEL_Y_PAD, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="bottom",
              transform=fig.transFigure)


def plot_diagonal(ax, series_1, series_2):
    lo = min(series_1.min(), series_2.min())
    hi = max(series_1.max(), series_2.max())
    nc = stylia.NamedColors()
    ax.plot([lo, hi], [lo, hi], linestyle="--", color=nc.silver, zorder=1)


def plot_all_hits(ax):
    nc = stylia.NamedColors()
    bg = audit[~audit["prioritized"]]
    fg = audit[audit["prioritized"]]
    ax.scatter(bg["best_on_target"], bg["human_best_af2"], color=nc.silver, alpha=0.4, label="Filtered hits (1,095)")
    ax.scatter(fg["best_on_target"], fg["human_best_af2"], color=nc.crimson, label="Prioritized (244)")
    plot_diagonal(ax, audit["best_on_target"], audit["human_best_af2"])
    ax.legend()
    stylia.label(ax, xlabel="Best on-target score (Mtb)", ylabel="Best off-target score (human)")


def plot_prioritized_by_origin(ax):
    fg = audit[audit["prioritized"]]
    pal = stylia.CategoricalPalette("npg")
    colors = pal.get(len(ORIGIN_LABELS))
    color_by_label = dict(zip(ORIGIN_LABELS.values(), colors))
    for label, color in color_by_label.items():
        sub = fg[fg["origin_label"] == label]
        ax.scatter(sub["best_on_target"], sub["human_best_af2"], color=color, label=label)
    plot_diagonal(ax, audit["best_on_target"], audit["human_best_af2"])
    ax.legend()
    stylia.label(ax, xlabel="Best on-target score (Mtb)", ylabel="Best off-target score (human)")


def plot_stacked_bars(ax, counts, title):
    nc = stylia.NamedColors()
    color_by_cat = {"docking only": nc.cobalt, "boltz2 only": nc.tangerine, "both methods": nc.crimson}
    x = range(len(counts.index))
    bottom = [0] * len(counts.index)
    for cat in AGREEMENT_ORDER:
        values = counts[cat].values
        ax.bar(x, values, bottom=bottom, color=color_by_cat[cat], label=AGREEMENT_LABELS[cat])
        bottom = [b + v for b, v in zip(bottom, values)]
    ax.set_xticks(list(x))
    ax.set_xticklabels([h.capitalize() for h in counts.index])
    ax.legend(loc="upper left")
    stylia.label(ax, xlabel="Hit type", ylabel="Selection slots (targets x top-N)", title=title)


def main():
    fig, axs = stylia.create_figure(2, 2, width=0.6, height=0.6)
    panel_axes = []

    ax_a = axs.next()
    plot_all_hits(ax_a)
    panel_axes.append(ax_a)

    ax_b = axs.next()
    plot_prioritized_by_origin(ax_b)
    panel_axes.append(ax_b)

    ax_c = axs.next()
    plot_stacked_bars(ax_c, cat_counts, "Catalytic")
    panel_axes.append(ax_c)

    ax_d = axs.next()
    plot_stacked_bars(ax_d, noncat_counts, "Non-catalytic")
    panel_axes.append(ax_d)

    # stylia.save_figure() always finishes with a bare plt.tight_layout() (default padding),
    # which would override any spacing set here - so the extra gap between panels is applied
    # with our own tight_layout(h_pad=.., w_pad=..) call, and the file is written directly with
    # the same savefig() call stylia.save_figure() uses internally, instead of going through it.
    fig.tight_layout(h_pad=4.0, w_pad=3.0)

    # Needs a real renderer for add_panel_label's tightbbox measurements.
    fig.canvas.draw()
    for ax, letter in zip(panel_axes, ["a", "b", "c", "d"]):
        add_panel_label(fig, ax, letter)

    pdf_path = os.path.join(plots_dir, "ExtendedDataFigure6.pdf")
    png_path = os.path.join(plots_dir, "ExtendedDataFigure6.png")
    for path in (pdf_path, png_path):
        plt.savefig(path, dpi=600, transparent=False, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
