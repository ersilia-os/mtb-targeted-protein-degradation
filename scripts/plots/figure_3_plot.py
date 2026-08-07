"""
Figure 3: multi-target docking hits from figure_3_calculations.py. A 2-column, 5-row grid; panels a,
b, g and h have real content so far, panels c-f are empty placeholders labeled for now, pending
content decisions.

Panel a follows the nested-overlaid-bars convention from
notebooks/46_docking_exploration_deliverables.ipynb: one bar per protein per cutoff, drawn loosest
cutoff (most compounds) first and strictest (fewest compounds) last, so each stricter cutoff's bar
sits fully inside the previous one instead of a standard grouped/stacked layout. The population is
figure_3_calculations.py's SELECTED_SET_CUTOFF multi-target hit set (compounds already selected as
hitting >= MIN_TARGETS genes at that cutoff), not the full screened library - so the y-axis ceiling
is that set's own size, and every bar height is a count out of that fixed total. Proteins sorted by
their compound count at cutoff -11, descending.

Panel b is a chord (Circos) diagram of pairwise protein hit overlap, at CIRCOS_CUTOFF, within that
same selected set - following notebooks/46_docking_exploration_IIIa_colors.ipynb's precedent
(pycirclize's Circos.chord_diagram, rasterized via FigureCanvasAgg into a stylia axes), colored with
figure_1's own gene->color mapping so protein colors stay consistent across every figure. Genes with
zero hits at the cutoff are dropped outright (a chord diagram sector needs a nonzero span); among
the rest, low-degree genes get a blank (padded-space) label instead of their real name - same
declutter heuristic as that notebook - since their sector is too small for real text without
overlapping its neighbors.

Panel g is a nested 2x4 grid (figure_3_plot.py's SHOWCASE_GRID_ROWS/SHOWCASE_GRID_COLS), filled only
in its last SHOWCASE_GRID_FILLED_COLS columns per row per request, of PyMOL renders of ONE showcase
compound (figure_3_calculations.py's SHOWCASE_COMPOUND_ID, rank 1 by n_targets/best-score in the
cutoff-11 multi-target hit set) docked into each hit gene's ACTUAL best-scoring pocket
(figure_3_calculations.py's compute_showcase_compound_pockets - the same per-gene min-score pocket
compute_gene_min_scores uses to call a hit, not merely the highest-P2Rank-probability pocket
notebooks/46_pocket_visualisation.ipynb simplifies to). Docked-pose SDFs were only ever retained for
14 of 276 pockets, so only the genes whose actual best pocket is among those 14 are shown here
(currently alaS, aspS, lysS, pheS - 4 of the compound's 10 hit genes, filling the grid's 2x2 block
exactly) - the rest have a docking score but no retained 3D pose, so they're omitted rather than
faked. Rendering follows scripts/47b_reference_pocket_visualization.py's established PyMOL recipe
(reference structure translucent cartoon, docked ligand pose extracted from docking.tar.gz and
colored with the project's standard orange-carbon convention), simplified to one ligand/no pocket
sphere/no PDB or AlphaFill overlays, and zoomed tightly on the ligand alone so all panels are framed
comparably. Each render carries three corner badges: gene+UniProt AC (colored via figure_1's
gene_to_color), an InterPro-derived CAT/Other/NA pocket-domain classification (real data from
output/pocket_detection_data_interpro.tsv, not fabricated), and the docking score.

Panel h is a rank-ordered area plot of P2Rank pocket probabilities across all 276 detected pockets
(output/pocket_detection_data.csv, via figure_3_calculations.py's compute_pocket_scores), sorted
probability descending - x-axis is that sort rank (1-276, arbitrary pocket identity), not a
biological quantity.

Usage:
    python figure_3_plot.py [--rerun]
"""
import argparse
import json
import os
import string
import sys
import tarfile
import tempfile

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

os.environ["QT_QPA_PLATFORM"] = "offscreen"

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
import stylia
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.lines import Line2D
from pycirclize import Circos
from pymol import cmd
from stylia.config import get_size
from stylia.figure.figure import stylize

from docking_utils import LIBRARIES

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_3")
os.makedirs(plots_dir, exist_ok=True)

N_ROWS = 3
N_COLS = 2
# Trailing rows appended after the N_ROWS x N_COLS grid, each a single full-width panel (spanning
# both columns) rather than split 2-up like the rows above - same total width as a and b combined.
N_EXTRA_ROWS = 2

# Row 0 (panel a/b) and the first extra row (panel g, the showcase image grid) are made taller per
# request - rows 1-(N_ROWS-1) (the still-empty placeholders) and panel h are left as-is. A share
# ratio alone can't do that (it only redistributes a FIXED total figure height between rows), so the
# overall figure height is bumped too, or the taller rows' content still gets squeezed/clipped
# against the figure's own edges regardless of how large its ratio is. Columns are 2.5:1 (a:b) rather
# than the grid's default even split - panel a's bar chart genuinely needs the width for 21 xtick
# labels, while panel b's roughly-square chord diagram just leaves a wide empty gutter in an equally-
# wide column. Panel g's PyMOL structure renders need more vertical room than a simple line/area plot
# (panel h) to read as recognizable pocket/ligand geometry rather than tiny blobs.
ROW_HEIGHT_RATIOS = [5.5] + [1] * (N_ROWS - 1) + [6, 1]
FIGURE_HEIGHT = 1.2
COLUMN_WIDTH_RATIOS = [1.2, 1]

# Loosest -> strictest, matching the draw order (background bar first, tightest cutoff on top) and
# figure_3_calculations.py's PANEL_A_CUTOFFS/SELECTED_SET_CUTOFF.
CUTOFFS = [-8, -9, -10, -11, -12]
SELECTED_SET_CUTOFF = -11
SORT_BY_CUTOFF = -11
YLIM_MAX = 370


def plot_protein_hit_counts(ax):
    counts = pd.read_csv(os.path.join(plots_dir, "figure_3_selected_set_protein_hit_counts.csv"))
    counts = counts.sort_values(f"count_cutoff{abs(SORT_BY_CUTOFF)}", ascending=False).reset_index(drop=True)

    palette = stylia.CategoricalPalette("npg").get(len(CUTOFFS))
    for x, row in counts.iterrows():
        for c, cutoff in enumerate(CUTOFFS):
            label = str(cutoff) if x == 0 else None
            ax.bar(x, row[f"count_cutoff{abs(cutoff)}"], color=palette[c], edgecolor="black",
                   linewidth=stylia.LINEWIDTH, zorder=2, label=label)

    ax.set_ylim(0, YLIM_MAX)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts["gene"], rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    # Outside the axes to the right, vertically centered - no xlim padding trick needed since it
    # doesn't compete with the plot area at all.
    ax.legend(title="Docking score", loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=stylia.FONTSIZE_SMALL)
    stylia.label(ax, xlabel="", ylabel="Number of compounds", abc="a")


CIRCOS_CUTOFF = -11

# Genes contributing less than this fraction of the diagram's total pairwise overlap would get a
# blank (padded-space, still correctly colored) sector label instead of their real name - per
# request, every gene keeps its real label now except CIRCOS_HIDE_LABELS below, so this is 0 (never
# triggers). Was 1/25, matching notebooks/46_docking_exploration_IIIa_colors.ipynb's threshold.
CIRCOS_DECLUTTER_FRACTION = 0

# Blanked regardless of the threshold above - every gene here has a tiny sector (count_cutoff11 of
# 1-15, vs. 45+ for the smallest of the labeled genes) and they all cluster together at the same
# point on the circle, so their real-text labels stack on top of each other into an unreadable
# smear rather than reading as separate sectors.
CIRCOS_HIDE_LABELS = {"trpS", "serS", "gltS", "leuS", "metS", "pheT", "proS", "thrS"}


def plot_circos_overlap(ax):
    hits_path = os.path.join(plots_dir, f"figure_3_multi_target_hits_cutoff{abs(CIRCOS_CUTOFF)}.csv")
    selected = pd.read_csv(hits_path)
    gene_cols = [c.removeprefix("score_") for c in selected.columns if c.startswith("score_")]

    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        gene_to_color = json.load(f)["gene_to_color"]

    hit_sets = {g: set(selected.loc[selected[f"score_{g}"] <= CIRCOS_CUTOFF, "compound_id"]) for g in gene_cols}
    genes = [g for g in gene_cols if hit_sets[g]]

    matrix = pd.DataFrame(0, index=genes, columns=genes, dtype=int)
    for g1 in genes:
        for g2 in genes:
            if g1 != g2:
                matrix.loc[g1, g2] = len(hit_sets[g1] & hit_sets[g2])

    node_strength = matrix.sum(axis=0) + matrix.sum(axis=1)
    order = node_strength.sort_values(ascending=False).index
    matrix = matrix.loc[order, order]
    node_strength = node_strength.loc[order]

    # Declutter: low-degree genes' sectors are too small for real text without the labels
    # overlapping each other, so give them a unique blank label instead (still correctly colored
    # via cmap, just not named) - same heuristic and threshold as the notebook precedent.
    cmap = {g: gene_to_color[g] for g in genes}
    threshold = (node_strength.sum() / 2) * CIRCOS_DECLUTTER_FRACTION
    rename_map = {}
    for i, g in enumerate(matrix.index, start=1):
        if node_strength[g] < threshold or g in CIRCOS_HIDE_LABELS:
            blank = " " * i
            rename_map[g] = blank
            cmap[blank] = gene_to_color[g]
    matrix = matrix.rename(index=rename_map, columns=rename_map)

    # Chord diagrams double-count a symmetric matrix's upper+lower triangles - zero the lower
    # triangle (diagonal is already 0, self-overlap is meaningless) so each pair's link is drawn once.
    matrix.values[np.tril_indices_from(matrix.values, k=0)] = 0

    circos = Circos.chord_diagram(
        matrix, space=2, cmap=cmap,
        label_kws=dict(size=30),
        link_kws=dict(ec="black", lw=1, alpha=1),
    )
    fig_circos = circos.plotfig(figsize=(8, 8), dpi=200)
    canvas = FigureCanvas(fig_circos)
    canvas.draw()
    w, h = fig_circos.canvas.get_width_height()
    img = np.frombuffer(canvas.tostring_argb(), dtype=np.uint8).reshape(h, w, 4)[:, :, [1, 2, 3, 0]]
    plt.close(fig_circos)

    # plotfig()'s own canvas leaves a lot of blank margin around the actual diagram - crop to the
    # non-white content bounding box (as figure_1_plot.py's autocrop_to_content does for the PyMOL
    # renders) so the diagram fills panel b instead of sitting small in a mostly-empty box.
    nonwhite = np.where((img[..., :3] < 250).any(axis=2))
    y0, y1 = nonwhite[0].min(), nonwhite[0].max()
    x0, x1 = nonwhite[1].min(), nonwhite[1].max()
    img = img[y0:y1 + 1, x0:x1 + 1]

    ax.imshow(img)
    ax.axis("off")
    stylia.label(ax, xlabel="", ylabel="", abc="b")


def plot_pocket_scores(ax):
    scores = pd.read_csv(os.path.join(plots_dir, "figure_3_pocket_scores.csv"))
    nc = stylia.NamedColors()
    ax.plot(scores["pocket_rank"], scores["pocket_probability"], color=nc.crimson, linewidth=stylia.LINEWIDTH, zorder=2)
    ax.fill_between(scores["pocket_rank"], scores["pocket_probability"], color=nc.crimson, alpha=0.3, zorder=1)
    ax.set_xlim(scores["pocket_rank"].min(), scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    stylia.label(ax, xlabel="Pocket number", ylabel="P2Rank probability", abc="h")


# Matches scripts/47b_reference_pocket_visualization.py's COLOR_LIGAND_DOCKED / notebooks/
# 46_pocket_visualisation.ipynb's convention - orange carbons for a docked (as opposed to
# experimental/AlphaFill/homolog) ligand pose, kept consistent across every PyMOL render in the repo.
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]
COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]

# Angstrom cube centered on the ligand alone (not the whole structure, unlike figure_1_plot.py's/
# 47b's renders) so every panel-g sub-image is framed at a comparable, tightly-zoomed scale
# regardless of each protein's own size - a framing choice, tune by eye if too tight or too sparse.
SHOWCASE_ZOOM_BOX = 22.0


def _color_ligand(selection, carbon_color, color_name):
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def _zoom_to_fixed_box(selection, box_size, box_name="_zoom_box"):
    (xmin, ymin, zmin), (xmax, ymax, zmax) = cmd.get_extent(selection)
    cx, cy, cz = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax)
    h = box_size / 2.0
    cmd.delete(box_name)
    corners = [
        (cx - h, cy - h, cz - h), (cx - h, cy - h, cz + h),
        (cx - h, cy + h, cz - h), (cx - h, cy + h, cz + h),
        (cx + h, cy - h, cz - h), (cx + h, cy - h, cz + h),
        (cx + h, cy + h, cz - h), (cx + h, cy + h, cz + h),
    ]
    for i, (x, y, z) in enumerate(corners):
        cmd.pseudoatom(f"{box_name}_{i}", pos=[x, y, z])
    cmd.group(box_name, f"{box_name}_*")
    cmd.zoom(box_name, buffer=0.0, complete=1)
    cmd.delete(box_name)
    cmd.delete(f"{box_name}_*")


def render_showcase_pocket(row, rerun=False):
    png_path = os.path.join(plots_dir, f"figure_3_showcase_{row['gene']}.png")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing render for {row['gene']}: {png_path}")
        return png_path

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    structure_path = os.path.join(root, "..", "..", "output", "aligned_relaxed_structures",
                                   row["uniprot_ac"], f"{row['structure_name']}.pdb")
    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color", COLOR_STRUCTURE)
    cmd.color("structure_color", "structure")
    cmd.show("cartoon", "structure")
    cmd.hide("lines", "structure")
    cmd.hide("sticks", "structure")
    # No pocket-sphere overlay (per request) - gene identity/color is already carried by panel g's
    # own corner badge, so the cartoon's only job here is to stay out of the ligand's way. Higher
    # transparency than the original 0.3 for exactly that reason.
    # cartoon_transparency, not the generic transparency setting - the latter governs surface
    # representation and is a no-op on a cartoon (confirmed by the first render still showing a
    # fully-opaque loop occluding the ligand despite "transparency" being set).
    cmd.set("cartoon_transparency", 0.6, "structure")

    tar_path = os.path.join(LIBRARIES["REAL"], row["pocket_name"], "docking.tar.gz")
    with tarfile.open(tar_path, "r|gz") as tf:
        for member in tf:
            if member.name != f"docking/{row['compound_id']}_out.sdf":
                continue
            data = tf.extractfile(member).read()
            with tempfile.NamedTemporaryFile(suffix=".sdf", mode="wb", delete=False) as tmp:
                tmp.write(data)
                tmp_path = tmp.name
            cmd.load(tmp_path, "ligand")
            os.remove(tmp_path)
            break
    cmd.show("sticks", "ligand")
    cmd.hide("lines", "ligand")
    _color_ligand("ligand", COLOR_LIGAND_DOCKED, "ligC_docked")

    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    # PyMOL's default transparency_mode doesn't alpha-blend cartoon geometry correctly during ray
    # tracing at this render's tight, ligand-centered zoom - a nearby loop rendered edge-on showed up
    # as a solid opaque black blob over the ligand despite cartoon_transparency being set. Mode 2
    # fixes the blending (confirmed by isolated testing) and also makes the transparency actually
    # visible instead of still reading as near-opaque.
    cmd.set("transparency_mode", 2)
    cmd.set("antialias", 2)
    _zoom_to_fixed_box("ligand", box_size=SHOWCASE_ZOOM_BOX)
    cmd.ray(1200, 1200)
    cmd.png(png_path, dpi=600)

    return png_path


# Panel g's nested grid: SHOWCASE_GRID_COLS wide, filled only in its last SHOWCASE_GRID_FILLED_COLS
# columns per row (per request - mimics a reference layout where the first columns are reserved for
# other content this repo doesn't have yet), SHOWCASE_GRID_ROWS tall. 4 genes currently have a
# retained pose, filling exactly a 2x2 block in the last 2 of 4 columns.
SHOWCASE_GRID_ROWS = 2
SHOWCASE_GRID_COLS = 4
SHOWCASE_GRID_FILLED_COLS = 2


def _corner_badge(sub_ax, dot_color, label, loc, anchor_point):
    # Reuses figure_1_plot.py's colored-circle Line2D-handle-as-legend trick (already used
    # project-wide for gene color coding), called multiple times per axes via add_artist() so each
    # corner badge survives the next call instead of replacing it (ax.legend() alone only keeps the
    # most recent one). anchor_point is in DATA (pixel) coordinates, not axes-fraction - with
    # adjustable="datalim" the axes BOX can be letterboxed around the actual image (blank padding on
    # the sides), so an axes-fraction anchor drifts off the visible image; pixel coordinates always
    # track the image's own corners regardless of that padding.
    handle = [Line2D([0], [0], marker="o", color="w", label=label,
                      markerfacecolor=dot_color, markeredgecolor="black", markeredgewidth=0.5,
                      markersize=stylia.MARKERSIZE)]
    legend = sub_ax.legend(handles=handle, loc=loc, bbox_to_anchor=anchor_point,
                            bbox_transform=sub_ax.transData, frameon=True,
                            framealpha=0.85, edgecolor="black", fontsize=stylia.FONTSIZE_SMALL,
                            handletextpad=0.3, borderpad=0.3, labelspacing=0)
    sub_ax.add_artist(legend)


def plot_showcase_strip(ax, rerun=False):
    # keep_default_na=False - otherwise pandas silently reads the literal string "NA" in
    # interpro_categories back as a real NaN (pandas' own default missing-value token list).
    pockets = pd.read_csv(os.path.join(plots_dir, "figure_3_showcase_compound_pockets.csv"), keep_default_na=False)
    pockets = pockets[pockets["has_pose"]].sort_values("gene").reset_index(drop=True)

    with open(os.path.join(root, "..", "..", "output", "plots", "figure_1", "color_mapping.json")) as f:
        gene_to_color = json.load(f)["gene_to_color"]
    nc = stylia.NamedColors()

    fig = ax.figure
    pos = ax.get_position()
    subgs = ax.get_subplotspec().subgridspec(SHOWCASE_GRID_ROWS, SHOWCASE_GRID_COLS, wspace=0.05, hspace=0.15)
    ax.remove()

    first_col = SHOWCASE_GRID_COLS - SHOWCASE_GRID_FILLED_COLS
    for i, row in pockets.iterrows():
        grid_row, grid_col = divmod(i, SHOWCASE_GRID_FILLED_COLS)
        sub_ax = stylize(fig.add_subplot(subgs[grid_row, first_col + grid_col]))
        img = mpimg.imread(render_showcase_pocket(row, rerun=rerun))

        # Inline autocrop-to-content, same convention as panel b's circos crop, but the threshold is
        # scaled for mpimg.imread's float32-in-[0,1] PNG data rather than panel b's uint8 ARGB buffer.
        nonwhite = np.where((img[..., :3] < 250 / 255).any(axis=2))
        y0, y1 = nonwhite[0].min(), nonwhite[0].max()
        x0, x1 = nonwhite[1].min(), nonwhite[1].max()
        img = img[y0:y1 + 1, x0:x1 + 1]

        sub_ax.imshow(img)
        sub_ax.axis("off")
        # Same fix as figure_1_plot.py's show_zoomed_image: without this, imshow's default aspect
        # handling shrinks/repositions the AXES BOX to match the image, leaving blank axes-fraction
        # padding around the actually-visible image - and the corner badges below anchor to that
        # full (padded) box via bbox_to_anchor, not the visible image, so they'd float off it.
        sub_ax.set_aspect("equal", adjustable="datalim")

        # imshow's default (origin="upper") data extent puts (0, 0) at the image's top-left corner
        # and (w, h) at its bottom-right - these are the pixel-space anchor points _corner_badge uses.
        img_h, img_w = img.shape[:2]

        _corner_badge(sub_ax, gene_to_color[row["gene"]], f"{row['gene']} ({row['uniprot_ac']})",
                      loc="upper left", anchor_point=(0, 0))

        # CAT/Other/NA - real InterPro-derived classification from figure_3_calculations.py's
        # compute_showcase_compound_pockets (output/pocket_detection_data_interpro.tsv), not
        # fabricated. Filled dot when this pocket has any curated annotation, open dot for NA.
        interpro_label = row["interpro_categories"]
        interpro_color = "white" if interpro_label == "NA" else nc.lime
        _corner_badge(sub_ax, interpro_color, interpro_label, loc="upper right", anchor_point=(img_w, 0))

        _corner_badge(sub_ax, nc.amber, f"Docking score: {row['score']:.3f}",
                      loc="lower left", anchor_point=(0, img_h))

        stylia.label(sub_ax, xlabel="", ylabel="", abc="g" if i == 0 else None)

    # va="top" (text drawn BELOW y, i.e. inside panel g's own row) rather than va="bottom" at the
    # same y - the row boundary sits right against the row above it, and va="bottom" would draw
    # upward into that neighboring row's own margin/labels.
    compound_id = pockets["compound_id"].iloc[0]
    fig.text((pos.x0 + pos.x1) / 2, pos.y1, f"Compound {compound_id} (EnamineREAL10B)",
              ha="center", va="top", fontsize=stylia.FONTSIZE_SMALL)


def main(rerun=False):
    # Bypasses stylia.create_figure() (a plain nrows x ncols grid, every cell the same span) since
    # the trailing N_EXTRA_ROWS need a single axes spanning both columns instead of two split cells
    # - built by hand the same way figure_1_plot.py's more complex layouts are, then each axes is
    # run through stylia's own stylize() so it still matches every other panel's look.
    size = get_size()
    fig = plt.figure(figsize=(size, FIGURE_HEIGHT * size))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(N_ROWS + N_EXTRA_ROWS, N_COLS,
                           height_ratios=ROW_HEIGHT_RATIOS, width_ratios=COLUMN_WIDTH_RATIOS)

    labels = iter(string.ascii_lowercase[:N_ROWS * N_COLS + N_EXTRA_ROWS])
    for row in range(N_ROWS):
        for col in range(N_COLS):
            label = next(labels)
            ax = stylize(fig.add_subplot(gs[row, col]))
            if label == "a":
                plot_protein_hit_counts(ax)
            elif label == "b":
                plot_circos_overlap(ax)
            else:
                stylia.label(ax, xlabel="", ylabel="", abc=label)
    # Images row placed before the P2Rank distribution row (per request: showcase images sit
    # between the a/b row and the pocket-score distribution) - labels are assigned by physical
    # position, so the showcase strip is now "g" and the distribution is now "h".
    for extra_row in range(N_EXTRA_ROWS):
        next(labels)
        ax = stylize(fig.add_subplot(gs[N_ROWS + extra_row, :]))
        if extra_row == 0:
            plot_showcase_strip(ax, rerun=rerun)
        else:
            plot_pocket_scores(ax)

    for ext in ("png", "pdf"):
        output_path = os.path.join(plots_dir, f"figure_3.{ext}")
        stylia.save_figure(output_path)
        print(f"Saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force PyMOL re-rendering of panel g's showcase-compound PNGs, "
                              "even if they already exist")
    args = parser.parse_args()
    main(rerun=args.rerun)
