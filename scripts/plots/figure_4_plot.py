"""
Figure 4, panel a: PyMOL screenshot of PDB 7K98, the crystal structure of M. tuberculosis PheRS
(pheS/pheT, the project's Tier-1 target) in its "preaminoacylation" complex with cognate tRNA(Phe)
and the phenylalanyl-adenylate mimic F-AMS (ligand code W5Y). This structure is already used
downstream in the pipeline as 7K98_pocket_6 ("pocket at dimerization interface",
output/selected_pockets.csv, scripts 48-50/65/71-75) - panel 4a introduces the real experimental
structure behind that pocket.

Source: data/structures/pdbe_database/P9WFU3/P9WFU3_updated/7k98_updated.cif. Citation (from the
deposited file's own JRNL header, not invented): Michalska K, Jedrzejczak R, Wower J, Chang C,
Baragana B, Gilbert IH, Forte B, Joachimiak A. "Mycobacterium tuberculosis Phe-tRNA synthetase:
structural insights into tRNA recognition and aminoacylation." Nucleic Acids Res. 2021;49(9):5351.
DOI: 10.1093/nar/gkab272. Deposited 2020-09-28, resolution 2.19 A.

Shows the full deposited complex: PheRS alpha subunit (pheS, chains A/D), PheRS beta subunit
(pheT, chains B/E, one copy - chain B - also shown as an opaque surface), tRNA(Phe) (chains C/F),
and the native co-crystallized ligand W5Y (F-AMS, fuchsia sticks - not a docked pose, so
intentionally NOT the project's usual orange COLOR_LIGAND_DOCKED). Solvent (waters) and metallic
atoms (the bound Mg2+ ions), along with Cl-/glycerol, are not shown. Every custom color (pheS
cobalt, pheT amber, tRNA turquoise, ligand fuchsia) is sourced from stylia.ArticleColors() rather
than a native PyMOL color name or a hand-picked RGB, so the palette stays tied to the same source
as every other stylia-driven plot in the project - and pheS/pheT get distinct colors rather than
the project's own near-identical gene_to_color hexes for this pair (adjacent on the alphabetical
spectral palette), since distinguishing the two subunits is this panel's whole point.

Also shown: the 3 curated pheS/pheT pockets from output/selected_pockets.csv, as large spheres
(add_pocket_spheres()), one distinct article-palette color per pocket (POCKET_COLORS). That file
actually lists 4 pheS/pheT rows, but
alphafold3_P9WFU3_model_2_pocket_2 ("falls at interface") and 7K98_pocket_6 ("pocket at
dimerization interface") describe the same physical pocket, independently detected on an AF3
model and on this real structure - shown once, via 7K98_pocket_6's own native centroid (already
in this structure's own coordinate frame, unlike the other two pockets, which are P2Rank pockets
detected on separate monomeric models and need aligning onto "structure" first - see
_add_aligned_pocket_sphere()). User-confirmed dedup.

Panel size comes from output/plots/figure_4/panel_layout.csv (columns: panel, x, delta_x, y,
delta_y, cm - x/y reserved for future external-layout positioning), same convention as
figure_3_plot.py's own per-panel-standalone-PDF setup.

Panel b: 4x3 grid, one snapshot per output/selected_pockets.csv pocket (structure surface + best
docked ligand pose - see plot_pocket_grid_panel/render_pocket_snapshot for the full recipe and
its own POCKET_VIEWS per-pocket hand-tuned camera views). Still being tuned one pocket at a time.

Panel c: two robustness checks on the docking scores, density-colored (see _density -
histogram-based, not a gaussian_kde, to stay fast at 10k+ points). Left: mean Uni-Dock score vs.
std across script 65's 5 replicates, over the ~2,923 aggregated hits (output/62_aggregate_hits)
x 12 pockets - are the docking scores themselves reproducible? Right: mean docking score vs.
Boltz-2-predicted IC50 (output/75_boltz2_collect_affinities, converted from its
log10(IC50 uM) affinity_pred_value), over the ~1,095 FILTERED hits (output/70_filtering) x 12
pockets - a different, later-stage/smaller compound population than the left subplot, since
Boltz-2 only ran on the final filtered set - does an orthogonal structural method agree with the
docking score? Both subplots exclude a handful of extreme positive-score "invalid pose" outliers
from the view only (score >= 0: 60/35,076 points left, 15/13,140 right) - never dropped from the
underlying data, only out of frame; counts printed to console each run.

Panel d: the 2 protein-level (pheST/aspS/lysS/alaS) UpSet plots already rendered by script 68
(save_protein_upset / _score_upsets) - CAT pockets at docking score <= -11
(output/68_plot_results/upset_score_11_CAT.png) and NON-CAT pockets at docking score <= -10
(output/68_plot_results/upset_score_10_NONCAT.png), the two thresholds picked out of script 68's
own DOCKING_SCORE_THRESHOLDS sweep ([-12,-11,-10,-9], run for "all"/CAT/NON-CAT variants). These 2
PNGs are composited side by side, not regenerated here - script 68 is a numbered pipeline script
(can't be imported, its filename isn't a valid module identifier), and its own renders are already
the right content. Keeping them in sync with any future re-run of script 68 (different pockets,
thresholds, or --annotate) just means rerunning `68_plot_results.py --only-upset` before this panel.

Usage:
    python figure_4_plot.py [--rerun] [--subpanels a,b,c,d]
"""
import argparse
import json
import os
import re
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

os.environ["QT_QPA_PLATFORM"] = "offscreen"

import tarfile
import tempfile

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
import stylia
from matplotlib.lines import Line2D
from pymol import cmd
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr, rankdata, spearmanr
from stylia.figure.figure import stylize

from docking_utils import LIBRARIES, load_scores

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

ROOT = os.path.join(root, "..", "..")
plots_dir = os.path.join(ROOT, "output", "plots", "figure_4")
renders_dir = os.path.join(plots_dir, "renders")
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(renders_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c", "d"]
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

# Panel c
DOCKING_RESULTS_DIR_65 = os.path.join(ROOT, "output", "65_aggregated_docking", "docking_results")
FILTERED_HITS_CSV = os.path.join(ROOT, "output", "70_filtering", "filtered_hits.csv")
AFFINITY_RESULTS_CSV = os.path.join(ROOT, "output", "75_boltz2_collect_affinities", "affinity_results.csv")

# Panel d - script 68's own UpSet renders (see module docstring), composited as-is, not regenerated.
UPSET_RESULTS_DIR_68 = os.path.join(ROOT, "output", "68_plot_results")
UPSET_CAT_PNG = os.path.join(UPSET_RESULTS_DIR_68, "upset_score_11_CAT.png")
UPSET_NONCAT_PNG = os.path.join(UPSET_RESULTS_DIR_68, "upset_score_10_NONCAT.png")

SOURCE_CIF = os.path.join(ROOT, "data", "structures", "pdbe_database", "P9WFU3", "P9WFU3_updated",
                           "7k98_updated.cif")

# Every custom color below is sourced from stylia's own article palette (not a native PyMOL
# color name or a raw hand-picked RGB), so this render's palette stays tied to the same source
# as every stylia-driven plot in the project.
ac = stylia.ArticleColors()
COLOR_PHES = ac.cobalt      # alpha subunit, chains A/D
COLOR_PHET = ac.amber       # beta subunit, chains B/E
COLOR_TRNA = ac.turquoise   # replaces the native PyMOL "palecyan" used elsewhere in the project
                             # (scripts/77_pocket_annotation/11_build_pymol_sessions.py's
                             # COLOR_EXPERIMENTAL_TRNA) with the closest article-palette color.
# Closest article-palette match to the project's usual "native/experimental ligand" magenta
# (scripts/47b_reference_pocket_visualization.py's COLOR_LIGAND_PDB), still visually distinct
# from the orange COLOR_LIGAND_DOCKED reserved project-wide for docked poses (W5Y is not a
# docked pose).
COLOR_LIGAND = ac.fuchsia
COLOR_MG = "green"
# One distinct article-palette color per pocket (none reused from COLOR_PHES/PHET/TRNA/LIGAND
# above), keyed by pocket_name.
POCKET_COLORS = {
    "swissmodel_P9WFU3_model_0_pocket_1": ac.crimson,      # pheS, CAT
    "7K98_pocket_6": ac.lime,                              # pheS, NON-CAT (dimerization interface)
    "alphafold2_P9WFU1_model_0_pocket_1": ac.orchid,       # pheT, NON-CAT
}

# output/selected_pockets.csv's 3 curated pheS/pheT pockets (see module docstring dedup note),
# rendered as large spheres. 7K98_pocket_6 is already native to this structure's own coordinate
# frame (output/48_detect_pocket_multimers/); the other two are P2Rank pockets detected on
# separate monomeric models (output/pocket_detection_data.csv) and need aligning onto "structure"
# first - offset=0 for both P9WFU3 and P9WFU1 vs 7K98 in
# processed/pocket_annotation/chain_manifest.json, so pocket residue numbers correspond directly,
# no renumbering needed (unlike some other structure pairs in that manifest).
SELECTED_POCKETS_CSV = os.path.join(ROOT, "output", "selected_pockets.csv")
POCKET_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
MULTIMER_POCKET_DATA_CSV = os.path.join(ROOT, "output", "48_detect_pocket_multimers",
                                         "pocket_detection_data.csv")
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures")
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets")
# Extended beyond pheS/pheT (panel a's own pockets) to all 5 genes in output/selected_pockets.csv,
# for panel b's pocket grid below.
GENE_TO_UNIPROT = {"pheS": "P9WFU3", "pheT": "P9WFU1", "aspS": "P9WFW3", "lysS": "P9WFU9", "alaS": "P9WFW7"}
GENE_TO_TARGET_CHAIN = {"pheS": "A", "pheT": "B"}  # which "structure" chain to align onto
# alphafold3_P9WFU3_model_2_pocket_2 ("falls at interface") and 7K98_pocket_6 ("pocket at
# dimerization interface") describe the same physical pocket, independently detected on an AF3
# model and on this real structure - shown once, via 7K98_pocket_6's own native coordinates
# (user-confirmed dedup).
SKIP_POCKET_NAMES = {"alphafold3_P9WFU3_model_2_pocket_2"}

POCKET_SPHERE_RADIUS = 8.0
# Slightly transparent - fully opaque (0.0) read as too flat/solid; this keeps the sphere
# reading as a clean marker without completely hiding what's behind it.
POCKET_SPHERE_TRANSPARENCY = 0.15


def load_panel_sizes(panels):
    """panel -> (width_in, height_in), from output/plots/figure_4/panel_layout.csv - same
    convention as figure_3_plot.py's own load_panel_sizes."""
    if not os.path.exists(panel_layout_path):
        raise FileNotFoundError(
            f"{panel_layout_path} not found. Expected columns: panel, x, delta_x, y, delta_y (cm), "
            "one row per panel: " + ", ".join(PANEL_LETTERS))
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: (df.loc[p, "delta_x"] / 2.54, df.loc[p, "delta_y"] / 2.54) for p in panels}


def save_panel(fig, letter, use_tight_layout=True, tight_pad=1.08, tight_w_pad=None, subplots_adjust=None):
    # No bbox_inches="tight" - saves at exactly panel_layout.csv's delta_x/delta_y, same
    # convention as figure_3_plot.py's save_panel. use_tight_layout=False (panel b's grid, and
    # panel c below) skips plt.tight_layout() - it silently fails ("Axes not compatible with
    # tight_layout" warning; confirmed on panel c's shared-y/log-scale axes too, not just panel
    # b's bbox_to_anchor corner badges) and produces IDENTICAL output regardless of tight_pad,
    # clipping axis-label text at the saved edge since nothing here uses bbox_inches="tight" to
    # grow the canvas to fit. subplots_adjust (dict of left/right/top/bottom/wspace fractions)
    # replaces that failed auto-layout with explicit, hand-tuned margins; falls back to the old
    # near-zero-margin default (fine for panel b's label-free grid) when not given.
    if use_tight_layout:
        plt.tight_layout(pad=tight_pad, w_pad=tight_w_pad)
    else:
        fig.subplots_adjust(**(subplots_adjust or dict(left=0.01, right=0.99, top=0.99, bottom=0.01)))
    output_path = os.path.join(plots_dir, f"Fig_4{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_4{letter}.pdf")


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
    """Crop a PyMOL render to its non-white bounding box (plus a small margin) - same helper as
    figure_1_plot.py's autocrop_to_content."""
    gray = img[..., :3].mean(axis=2)
    max_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
    mask = gray < background_frac * max_val
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return img

    h, w = img.shape[:2]
    y0, y1 = rows[0], rows[-1]
    x0, x1 = cols[0], cols[-1]
    pad_y = int((y1 - y0) * padding_frac)
    pad_x = int((x1 - x0) * padding_frac)
    y0 = max(y0 - pad_y, 0)
    y1 = min(y1 + pad_y, h - 1)
    x0 = max(x0 - pad_x, 0)
    x1 = min(x1 + pad_x, w - 1)
    return img[y0:y1 + 1, x0:x1 + 1]


def _color_ligand(selection, carbon_color, color_name):
    """Color-by-element with a custom carbon color - same helper as figure_3_plot.py's
    _color_ligand / 47b's color_ligand."""
    cmd.util.cbag(selection)
    cmd.set_color(color_name, carbon_color)
    cmd.color(color_name, f"{selection} and elem C")


def _parse_pocket_residues(resi_str):
    """"D_148 D_154 ..." -> (chain, [148, 154, ...]) - same chain_resn format used throughout
    the project's pocket_detection_data.csv files. Assumes a single chain (true for every pocket
    used here)."""
    tokens = resi_str.split()
    chain = tokens[0].split("_")[0]
    resnums = [int(tok.split("_")[1]) for tok in tokens]
    return chain, resnums


def _add_native_pocket_sphere(obj_name, centroid_xyz, color_name):
    """7K98_pocket_6's own centroid, already in "structure"'s coordinate frame - no alignment."""
    cmd.pseudoatom(obj_name, pos=list(centroid_xyz), vdw=POCKET_SPHERE_RADIUS)
    # A freshly created pseudoatom defaults to some small marker representation (nonbonded
    # cross and/or nb_spheres dot) that stays on even after enabling spheres - hiding "nonbonded"
    # alone didn't clear it, so hide everything unconditionally first instead of guessing which
    # rep it actually is.
    cmd.hide("everything", obj_name)
    cmd.show("spheres", obj_name)
    cmd.color(color_name, obj_name)
    cmd.set("sphere_transparency", POCKET_SPHERE_TRANSPARENCY, obj_name)


def _add_aligned_pocket_sphere(obj_name, mobile_path, mobile_chain, resnums, centroid_xyz, target_chain, color_name):
    """Aligns a monomeric model's pocket (via its own resi CA atoms) onto "structure"'s
    corresponding chain, then keeps only a sphere at that pocket's centroid - now transformed
    into "structure"'s own coordinate frame."""
    resi_str = "+".join(str(r) for r in resnums)
    mobile_obj = f"_mobile_{obj_name}"
    cmd.load(mobile_path, mobile_obj)
    cmd.pseudoatom(f"_centroid_{obj_name}", pos=list(centroid_xyz), vdw=POCKET_SPHERE_RADIUS)
    cmd.create(mobile_obj, f"{mobile_obj} or _centroid_{obj_name}")
    cmd.delete(f"_centroid_{obj_name}")

    rms, n_atoms = cmd.align(
        f"{mobile_obj} and chain {mobile_chain} and resi {resi_str} and name CA and polymer.protein",
        f"structure and chain {target_chain} and resi {resi_str} and name CA and polymer.protein",
    )[:2]
    print(f"    Aligned {obj_name} onto structure chain {target_chain}: RMSD={rms:.2f} A over {n_atoms} atoms")

    cmd.create(obj_name, f"{mobile_obj} and resn PSD")
    cmd.delete(mobile_obj)
    # Same default-marker issue as _add_native_pocket_sphere - hide everything before showing spheres.
    cmd.hide("everything", obj_name)
    cmd.show("spheres", obj_name)
    cmd.color(color_name, obj_name)
    cmd.set("sphere_transparency", POCKET_SPHERE_TRANSPARENCY, obj_name)


def add_pocket_spheres():
    """Adds the 3 curated pheS/pheT pockets from output/selected_pockets.csv (see module
    docstring / SKIP_POCKET_NAMES) as large spheres, one distinct POCKET_COLORS color each."""
    for pocket_name, color in POCKET_COLORS.items():
        cmd.set_color(f"pocket_color_{pocket_name}", list(color))

    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    selected = selected[selected["gene_name"].isin(GENE_TO_TARGET_CHAIN) & ~selected["pocket_name"].isin(SKIP_POCKET_NAMES)]

    multimer_data = pd.read_csv(MULTIMER_POCKET_DATA_CSV)
    monomer_data = pd.read_csv(POCKET_DATA_CSV)

    for _, row in selected.iterrows():
        gene, pocket_name = row["gene_name"], row["pocket_name"]
        obj_name = f"pocket_{pocket_name}"

        if "_model_" not in pocket_name:
            # Native multimer pocket (e.g. 7K98_pocket_6) - already in "structure"'s own frame.
            structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)
            prow = multimer_data[(multimer_data["File name"] == f"{structure_name}.pdb")
                                  & (multimer_data["Pocket number"] == int(pocket_number))].iloc[0]
            centroid = [float(v) for v in prow["Pocket centroid coordinate (x y z)"].split()]
            _add_native_pocket_sphere(obj_name, centroid, f"pocket_color_{pocket_name}")
            print(f"  {gene}: {pocket_name} - native centroid, no alignment needed")
        else:
            structure_name, pocket_number = pocket_name.rsplit("_pocket_", 1)
            uniprot_ac = GENE_TO_UNIPROT[gene]
            prow = monomer_data[(monomer_data["Uniprot AC"] == uniprot_ac)
                                 & (monomer_data["File name"] == f"{structure_name}.pdb")
                                 & (monomer_data["Pocket number"] == int(pocket_number))].iloc[0]
            centroid = [float(v) for v in prow["Pocket centroid coordinate (x y z)"].split()]
            mobile_chain, resnums = _parse_pocket_residues(prow["Pocket residues (chain_resn)"])
            mobile_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
            _add_aligned_pocket_sphere(obj_name, mobile_path, mobile_chain, resnums, centroid,
                                        GENE_TO_TARGET_CHAIN[gene], f"pocket_color_{pocket_name}")


PNG_PATH = os.path.join(renders_dir, "figure_4a_7K98.png")
PSE_PATH = os.path.join(renders_dir, "figure_4a_7K98.pse")


def render_7k98(rerun=False):
    if os.path.exists(PNG_PATH) and not rerun:
        print(f"Reusing existing render: {PNG_PATH}")
        return PNG_PATH

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(SOURCE_CIF, "structure")
    cmd.hide("everything", "structure")

    # Only one copy of each half of the biological assembly is shown - chain A (pheS), chain B
    # (pheT), chain C (tRNA) - not the second D/E/F copy: all 3 pocket spheres below are only
    # meaningful relative to this A/B copy anyway (7K98_pocket_6's own residues are A+B, and the
    # two aligned pockets were placed onto chain A / chain B specifically), so D/E/F would just
    # be a redundant, unannotated duplicate competing for space in this already-tight framing
    # (user-confirmed trade-off against showing the full symmetric tetramer).
    cmd.set_color("phes_color", list(COLOR_PHES))
    cmd.set_color("phet_color", list(COLOR_PHET))
    cmd.show("cartoon", "structure and chain A and polymer.protein")
    cmd.color("phes_color", "structure and chain A and polymer.protein")
    cmd.show("cartoon", "structure and chain B and polymer.protein")
    cmd.color("phet_color", "structure and chain B and polymer.protein")

    # tRNA(Phe), chain C only (see above).
    cmd.set_color("trna_color", list(COLOR_TRNA))
    cmd.show("cartoon", "structure and chain C and polymer.nucleic")
    cmd.color("trna_color", "structure and chain C and polymer.nucleic")

    # Native co-crystallized ligand F-AMS (W5Y), chain A's copy only (see above).
    ligand_sel = "structure and resn W5Y and chain A"
    cmd.show("sticks", ligand_sel)
    _color_ligand(ligand_sel, COLOR_LIGAND, "ligC_w5y")

    # Chain B (pheT/beta, one of its two copies) also gets an opaque surface on top of its
    # cartoon, to show the subunit's molecular envelope. Computed on a standalone copy
    # (cmd.create), not directly on "structure and chain B" - PyMOL's surface calculation
    # otherwise still "sees" the rest of "structure" (chain A/C and beyond) as neighboring atoms,
    # which cuts a flat, artificially trimmed face where the selection boundary sits rather than
    # a true molecular envelope. The copy inherits chain B's already-applied phet_color (colored
    # above, before this point).
    cmd.create("chainB_surface", "structure and chain B and polymer.protein")
    # cmd.create carries over whatever reps were already on for these atoms (cartoon, from
    # earlier) - hidden here so it isn't drawn twice (once via "structure", once via this copy).
    cmd.hide("cartoon", "chainB_surface")
    cmd.show("surface", "chainB_surface")
    cmd.set("transparency", 0, "chainB_surface")

    # Solvent (waters) and metallic atoms (the bound Mg2+ ions) are never shown - along with
    # Cl-/glycerol, none of these crystallization/cryo additives are part of the rendered view.
    cmd.hide("everything", "structure and (solvent or metals)")

    # 3 curated pheS/pheT pockets (output/selected_pockets.csv), as large crimson spheres.
    add_pocket_spheres()

    # Standard project render settings, matching figure_1_plot.py / figure_3_plot.py /
    # 47b_reference_pocket_visualization.py's identical block.
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.02)
    cmd.set("spec_reflect", 0)
    cmd.set("specular", 0)
    cmd.set("antialias", 2)

    # Hand-tuned view (via PyMOL's own get_view), matching figure_1_plot.py's MY_VIEWS
    # per-structure convention.
    # cmd.set_view((
    #     0.394189388, -0.309228092, 0.865439832,
    #     -0.092989191, 0.923435509, 0.372304469,
    #     -0.914309084, -0.227235228, 0.335254282,
    #     -0.000425093, -0.000246342, -259.175048828,
    #     -6.897897720, 7.815885544, -67.098320007,
    #     159.797851562, 358.588134766, -20.000000000,
    # ))

    cmd.set_view ((
        0.302491635,    0.855443597,   -0.420365244,\
        -0.657938540,    0.506502688,    0.557280958,\
        0.689641178,    0.107999124,    0.716047466,\
        -0.000425093,   -0.000246342, -203.064010620,\
        -6.897897720,    7.815885544,  -67.098320007,\
    103.686798096,  302.477264404,  -20.000000000 ))


    cmd.ray(2400, 2400)
    cmd.png(PNG_PATH, dpi=600)
    cmd.save(PSE_PATH)
    print(f"Saved render: {PNG_PATH}")
    print(f"Saved session: {PSE_PATH}")
    return PNG_PATH


# ===========================================================================
# Panel b (prototype): 4x3 grid, one snapshot per output/selected_pockets.csv pocket - structure
# surface + pocket centroid sphere + top-1 docked ligand pose (user-confirmed style, over a
# plain surface+sphere or a cartoon+residues alternative). Grid order matches
# selected_pockets.csv's own row order, filled row-major (row 1 = pheS's 3 pockets, row 2 =
# pheT + aspS's first 2, etc.) - not one row per gene, since row counts differ per gene (4/3/2/3)
# and wouldn't tile a plain 4x3 grid evenly. Built one pocket at a time (--preview-pocket) before
# assembling the full grid, mirroring how panel a's own recipe (view/colors/pockets) was tuned
# through several rounds of visual feedback rather than decided upfront.
#
# Each cell also gets a gene-colored circle badge (_corner_badge, top-left) + a CAT/NON-CAT
# label (bottom center, white alpha-blended box) - same conventions as figure_1_plot.py's
# render_structure_panel (colored-circle-as-legend-handle trick) and figure_3_plot.py's
# _corner_badge / docking-score text box, applied here to every cell (not just the ones that
# rendered - see plot_pocket_grid_panel). Colors come from figure_1's own gene_to_color
# (output/plots/figure_1/color_mapping.json), matching every other figure's gene coloring -
# not the target-group scheme this used to have before the user asked to match prior figures.
# ===========================================================================

GENE_COLOR_MAPPING_PATH = os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")
SITE_TYPE_LABELS = {"CAT": "Catalytic", "NON-CAT": "Non-catalytic"}


def _corner_badge(ax, dot_color, label, loc="upper left", anchor_point=(0, 1)):
    """Colored-circle-as-legend-handle badge - same trick as figure_1_plot.py's
    render_structure_panel / figure_3_plot.py's _corner_badge. anchor_point defaults to
    axes-fraction (ax.legend()'s own default bbox_transform), so it works whether or not this
    cell actually has an image."""
    # stylia only defines 3 fontsize/markersize tiers (SMALL/plain/BIG - no "MEDIUM") - FONTSIZE_BIG
    # (8pt) read fine but its matched MARKERSIZE_BIG (30, a 3x jump from plain MARKERSIZE=10) made
    # the circle balloon past the label. Settled on the plain middle tier for both - FONTSIZE
    # (6pt, "medium") is legible at this grid's small per-cell size without the BIG tier's circle
    # blowing up, and MARKERSIZE at its own plain value actually fits that text.
    handle = [Line2D([0], [0], marker="o", color="w", label=label,
                      markerfacecolor=dot_color, markeredgecolor="black", markeredgewidth=0.5,
                      markersize=stylia.MARKERSIZE * 0.6)]
    legend = ax.legend(handles=handle, loc=loc, bbox_to_anchor=anchor_point, frameon=True,
                        framealpha=0.85, edgecolor="black", fontsize=stylia.FONTSIZE,
                        handletextpad=0.35, borderpad=0.3, labelspacing=0)
    ax.add_artist(legend)


def _site_type_label(ax, site_type):
    """CAT/NON-CAT label, bottom center (with a small gap above the axes' own bottom border,
    not flush against it), white alpha-blended box - same text-box style as figure_3_plot.py's
    docking-score label (_draw_compound_pose_row). FONTSIZE (stylia's plain "medium" tier) to
    match _corner_badge's text size."""
    ax.text(0.5, 0.06, SITE_TYPE_LABELS[site_type], transform=ax.transAxes,
            ha="center", va="bottom", fontsize=stylia.FONTSIZE, color="black",
            bbox=dict(facecolor="white", edgecolor="black", alpha=0.85, boxstyle="square,pad=0.3"))

POCKET_GRID_N_ROWS = 4
POCKET_GRID_N_COLS = 3

# Same neutral structure-surface color used throughout the project's pocket-visualization
# renders (scripts 47b/51, notebooks 46) - unlike panel a, this panel's surfaces aren't
# gene-colored, since here color is reserved for the target-group badge.
COLOR_STRUCTURE_NEUTRAL = [0.7804, 0.8275, 0.8667]
# Orange - the project's standard docked-pose color, verbatim from scripts/47b/51/figure_3_plot.py.
COLOR_LIGAND_DOCKED = [0xF5 / 255, 0xA6 / 255, 0x3A / 255]

# Same fixed Angstrom cube (centered on the ligand) as figure_3_plot.py's SHOWCASE_ZOOM_BOX, so
# every pocket snapshot in the grid is framed at the same physical scale regardless of protein size.
# 16.0 (figure_3's own value) put the camera essentially inside the pocket cavity here - close
# enough that the far cavity wall rendered as a solid black backface and the ligand washed out
# under it - backed off to give the surface room to actually curve around the ligand. Used as a
# fallback (auto orient + fixed-box zoom) for any pocket without its own file in VIEWS_DIR.
POCKET_SNAPSHOT_ZOOM_BOX = 28.0

# One user-dropped <pocket_name>.txt per pocket, each holding a raw PyMOL "set_view (...)" block
# (literally pasted from PyMOL's own get_view, backslash line continuations and the "### cut
# above here ###" trailer and all - _load_pocket_view() below just regexes out the 18 numbers,
# so the exact paste formatting doesn't matter). Same per-structure hand-tuned-view idea as
# figure_1_plot.py's MY_VIEWS, but file-based instead of a hardcoded dict, since each of the 12
# pockets sits on its own structure/coordinate frame and views get dropped in one at a time.
VIEWS_DIR = os.path.join(plots_dir, "views")
os.makedirs(VIEWS_DIR, exist_ok=True)


def _load_pocket_view(pocket_name):
    """Returns an 18-float cmd.set_view() tuple parsed out of VIEWS_DIR/<pocket_name>.txt, or
    None if that file doesn't exist yet OR is still an empty placeholder (every pocket got one
    pre-created via "create one empty txt file per pocket", ready for a view to be dropped in
    later) - both cases fall back to auto-orient + fixed-box zoom. Only raises if the file has
    SOME content that doesn't parse to exactly 18 numbers, a real formatting mistake."""
    view_path = os.path.join(VIEWS_DIR, f"{pocket_name}.txt")
    if not os.path.isfile(view_path):
        return None
    with open(view_path) as f:
        text = f.read()
    if not text.strip():
        return None
    numbers = [float(n) for n in re.findall(r"[-+]?\d+\.\d+", text)]
    if len(numbers) != 18:
        raise ValueError(f"{view_path} should contain exactly 18 numbers (a PyMOL set_view "
                          f"block), found {len(numbers)}")
    return tuple(numbers)


# Multimeric (experimental PDB, e.g. 7K98_pocket_6) pipeline paths - scripts 48-50, same
# convention as scripts/51_selected_pockets_visualization.py's own MULTIMER_STRUCTURES_DIR /
# MULTIMER_DOCKING_DIR.
MULTIMER_STRUCTURES_DIR = os.path.join(ROOT, "output", "48_detect_pocket_multimers", "stripped_structures")
MULTIMER_DOCKING_DIR = os.path.join(ROOT, "output", "50_unidock_docking_multimers")


def _resolve_pocket(gene_name, pocket_name):
    """(structure_path, docking_dir) for a pocket_name, branching on the monomeric
    ("_model_" in pocket_name - AF2/AF3/SwissModel/Chai1) vs multimeric (e.g. 7K98_pocket_6)
    convention - same resolution logic as scripts/51_selected_pockets_visualization.py's
    resolve_pocket()."""
    if "_model_" in pocket_name:
        structure_name, _pocket_number = pocket_name.rsplit("_pocket_", 1)
        uniprot_ac = GENE_TO_UNIPROT[gene_name]
        structure_path = os.path.join(ALIGNED_DIR, uniprot_ac, f"{structure_name}.pdb")
        docking_dir = os.path.join(LIBRARIES["REAL"], pocket_name)
    else:
        pdb_code, _pocket_number = pocket_name.rsplit("_pocket_", 1)
        structure_path = os.path.join(MULTIMER_STRUCTURES_DIR, f"{pdb_code}.pdb")
        docking_dir = os.path.join(MULTIMER_DOCKING_DIR, pocket_name)
    return structure_path, docking_dir


def _extract_ligand_sdf(docking_dir, compound_id):
    """Extracts docking/<compound_id>_out.sdf from docking_dir/docking.tar.gz to a temp file,
    returning its path (caller must remove it) - same tarfile idiom as figure_3_plot.py's
    render_showcase_pocket / scripts/51's load_ligand_objects."""
    tar_path = os.path.join(docking_dir, "docking.tar.gz")
    with tarfile.open(tar_path, "r|gz") as tf:
        for member in tf:
            if member.name == f"docking/{compound_id}_out.sdf":
                data = tf.extractfile(member).read()
                break
        else:
            raise FileNotFoundError(f"{compound_id}_out.sdf not found in {tar_path}")
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".sdf")
    os.close(tmp_fd)
    with open(tmp_path, "wb") as f:
        f.write(data)
    return tmp_path


def render_pocket_snapshot(gene_name, site_type, pocket_name, rerun=False):
    """Renders one output/selected_pockets.csv pocket: structure surface (neutral gray) + the
    pocket's best-scoring REAL round-2 docked pose (orange sticks, COLOR_LIGAND_DOCKED - matches
    every other docked-pose render project-wide). No pocket-centroid marker (dropped,
    user-confirmed - see comment below). Both monomeric and multimeric pockets are handled (see
    _resolve_pocket). Also saves a .pse session alongside the PNG, so an exact camera view
    worked out interactively in PyMOL can be copied back in via cmd.set_view() - same workflow
    as panel a's own hand-tuned MY_VIEWS-style view."""
    png_path = os.path.join(renders_dir, f"figure_4b_{pocket_name}.png")
    pse_path = os.path.join(renders_dir, f"figure_4b_{pocket_name}.pse")
    if os.path.exists(png_path) and not rerun:
        print(f"Reusing existing render: {png_path}")
        return png_path

    # Pocket centroid sphere dropped (user-confirmed) - pocket_N.pdb is the full cluster of pocket
    # residue atoms, not a single point, so a sphere there rendered as a huge, nearly-opaque blob
    # burying the ligand rather than a discreet marker. Just surface + ligand for now, matching
    # the user-supplied mockup exactly.
    structure_path, docking_dir = _resolve_pocket(gene_name, pocket_name)
    scores = load_scores(os.path.join(docking_dir, "report.csv"))
    best_compound = scores.idxmin()
    print(f"  {pocket_name}: best compound {best_compound}, score {scores.min():.3f}")
    sdf_path = _extract_ligand_sdf(docking_dir, best_compound)

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    cmd.load(structure_path, "structure")
    cmd.set_color("structure_color_neutral", COLOR_STRUCTURE_NEUTRAL)
    cmd.color("structure_color_neutral", "structure")
    cmd.hide("everything", "structure")
    cmd.show("surface", "structure")
    # 0.3, matching scripts/51_selected_pockets_visualization.py's own surface+ligand
    # convention - 0.1 left the surface nearly opaque, muting the ligand underneath it to a
    # washed-out gray instead of its own orange.
    cmd.set("transparency", 0.3, "structure")

    cmd.load(sdf_path, "ligand")
    os.remove(sdf_path)
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
    cmd.set("antialias", 2)

    view = _load_pocket_view(pocket_name)
    if view is not None:
        cmd.set_view(view)
    else:
        cmd.orient("ligand")
        _zoom_to_fixed_box_b("ligand", box_size=POCKET_SNAPSHOT_ZOOM_BOX)

    cmd.ray(1200, 1200)
    cmd.png(png_path, dpi=600)
    cmd.save(pse_path)
    print(f"Saved render: {png_path}")
    print(f"Saved session: {pse_path}")
    return png_path


def _zoom_to_fixed_box_b(selection, box_size, box_name="_zoom_box_b"):
    """Same fixed-cube zoom helper as figure_3_plot.py's _zoom_to_fixed_box."""
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


# Default "zoom in a bit" applied in the compositing step (not by touching each pocket's own
# PyMOL camera/view file) - a plain center-crop keeping this fraction of the width and height,
# cropped evenly off all 4 sides. Used whenever a pocket has no override in POCKET_ZOOMS.
POCKET_SNAPSHOT_ZOOM_FRAC = 0.8

# Per-pocket center-crop fraction override, hardcoded (not file-based, unlike POCKET_VIEWS) -
# fill in as each of the 12 pockets gets fine-tuned; a pocket not listed here just uses
# POCKET_SNAPSHOT_ZOOM_FRAC.
POCKET_ZOOMS = {}
POCKET_ZOOMS["swissmodel_P9WFU3_model_0_pocket_1"] = 0.6
POCKET_ZOOMS["alphafold2_P9WFU1_model_0_pocket_1"] = 0.6


def _center_crop(img, frac):
    h, w = img.shape[:2]
    y0, y1 = int(h * (1 - frac) / 2), int(h * (1 + frac) / 2)
    x0, x1 = int(w * (1 - frac) / 2), int(w * (1 + frac) / 2)
    return img[y0:y1, x0:x1]


def plot_pocket_grid_panel(letter, size, rerun=False):
    """Panel b: POCKET_GRID_N_ROWS x POCKET_GRID_N_COLS grid, one cell per output/
    selected_pockets.csv row, filled row-major in that file's own order (see module docstring).
    Every cell gets a gene-colored circle badge (_corner_badge) + CAT/NON-CAT label
    (_site_type_label), even ones that fail to render, so the grid stays informative about
    which pocket a blank cell belongs to. Both monomeric and multimeric pockets are handled
    (render_pocket_snapshot / _resolve_pocket) - the recipe is still being tuned one pocket at a
    time, so any pocket that fails to render for some other reason is left as a blank (but
    still labeled) cell with a printed warning, rather than crashing the whole panel."""
    selected = pd.read_csv(SELECTED_POCKETS_CSV)
    with open(GENE_COLOR_MAPPING_PATH) as f:
        gene_to_color = json.load(f)["gene_to_color"]

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(POCKET_GRID_N_ROWS, POCKET_GRID_N_COLS, hspace=0.05, wspace=0.05)

    for i, row in selected.iterrows():
        ax = stylize(fig.add_subplot(gs[i // POCKET_GRID_N_COLS, i % POCKET_GRID_N_COLS]))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("black")
            spine.set_linewidth(stylia.LINEWIDTH)
        stylia.label(ax, xlabel="", ylabel="")

        gene_name, site_type, pocket_name = row["gene_name"], row["site_type"], row["pocket_name"]
        _corner_badge(ax, gene_to_color[gene_name], gene_name)
        _site_type_label(ax, site_type)

        try:
            png_path = render_pocket_snapshot(gene_name, site_type, pocket_name, rerun=rerun)
        except Exception as e:
            print(f"  Warning: failed to render {pocket_name} ({gene_name}, {site_type}): {e}. Leaving blank.")
            continue
        img = _center_crop(mpimg.imread(png_path), POCKET_ZOOMS.get(pocket_name, POCKET_SNAPSHOT_ZOOM_FRAC))
        ax.imshow(img, interpolation="none")
        ax.set_aspect("equal", adjustable="datalim")

    save_panel(fig, letter, use_tight_layout=False)


def plot_structure_panel(letter, size, rerun=False):
    png_path = render_7k98(rerun=rerun)

    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)

    img = autocrop_to_content(mpimg.imread(png_path))
    ax.imshow(img, interpolation="none")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    stylia.label(ax, xlabel="", ylabel="")

    # use_tight_layout=False (same near-zero-margin path as panel b) so the render fills the
    # whole panel box instead of tight_layout's conservative padding around a blank/off axis.
    save_panel(fig, letter, use_tight_layout=False)


def load_replicate_robustness():
    """Pool script 65's 12 per-pocket results.csv files (compound, score, score_std - mean/std
    across 5 Uni-Dock replicate runs) into one long (compound_id, pocket_name, score, score_std)
    table over the ~2,923 aggregated hits x 12 pockets."""
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    frames = []
    for pocket_name in pockets:
        df = pd.read_csv(os.path.join(DOCKING_RESULTS_DIR_65, pocket_name, "results.csv"),
                          usecols=["compound", "score", "score_std"])
        df = df.rename(columns={"compound": "compound_id"})
        df["pocket_name"] = pocket_name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)[["compound_id", "pocket_name", "score", "score_std"]]


def load_docking_vs_affinity():
    """Mean docking score (output/70_filtering/filtered_hits.csv, one column per pocket) vs.
    Boltz-2-predicted IC50 in nM (output/75_boltz2_collect_affinities/affinity_results.csv's
    affinity_pred_value, which is log10(IC50 in uM) - see that script's own docstring), over the
    ~1,095 filtered hits x 12 pockets, joined on (compound_id, pocket_name)."""
    pockets = pd.read_csv(SELECTED_POCKETS_CSV)["pocket_name"].tolist()
    scores = pd.read_csv(FILTERED_HITS_CSV, usecols=["compound_id"] + pockets)
    scores = scores.melt(id_vars="compound_id", value_vars=pockets, var_name="pocket_name", value_name="score")

    affinity = pd.read_csv(AFFINITY_RESULTS_CSV, usecols=["compound_id", "pocket_name", "affinity_pred_value",
                                                           "affinity_probability_binary"])
    merged = scores.merge(affinity, on=["compound_id", "pocket_name"], how="inner")
    merged["ic50_nM"] = (10 ** merged["affinity_pred_value"]) * 1000
    return merged


def _density(x, y, bins=100, smoothing_sigma=2.0):
    """Fast (O(n + bins^2), not O(n^2)) per-point density: 2D histogram, Gaussian-smoothed (in
    bin units) to read as a continuous field rather than blocky per-bin cells, + bin-index
    lookup - avoids a gaussian_kde, which doesn't scale to this many points."""
    x, y = np.asarray(x), np.asarray(y)
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    hist = gaussian_filter(hist, sigma=smoothing_sigma)
    xidx = np.clip(np.digitize(x, xedges) - 1, 0, hist.shape[0] - 1)
    yidx = np.clip(np.digitize(y, yedges) - 1, 0, hist.shape[1] - 1)
    return hist[xidx, yidx]


def _density_scatter(ax, x, y, low_cutoff=0.35, log_x=False, log_y=False, zorder=2):
    """log_x/log_y=True bins density on log10(x)/log10(y) instead of the raw value - needed
    whenever the caller also sets ax.set_xscale("log")/ax.set_yscale("log"), otherwise the
    linear-space histogram bins can't resolve density among data that spans multiple orders of
    magnitude (e.g. IC50 in nM): the bulk of the points collapse into a couple of bins near zero,
    so the coloring stops tracking what the log-scale scatter actually looks dense. x/y themselves
    stay in their original (linear) units for plotting - only the density calculation's own
    binning changes."""
    x, y = np.asarray(x), np.asarray(y)
    density = _density(np.log10(x) if log_x else x, np.log10(y) if log_y else y)

    # Sort ascending by density so the densest (darkest) points are drawn last, on top - otherwise
    # sparse light points plotted later in the array can sit over the dense core.
    order = np.argsort(density)
    x, y, density = x[order], y[order], density[order]

    # Percentile-rank density into [0, 1] (robust to skew, same spirit as FadingColormap's own
    # "uniform" quantile transform), then rescale into [low_cutoff, 1] so the sparsest points
    # still get a visible color instead of fading out to the colormap's near-white end.
    frac = (rankdata(density, method="average") - 1) / (len(density) - 1)
    frac = low_cutoff + (1 - low_cutoff) * frac

    cm = stylia.FadingColormap("turquoise")
    ax.scatter(x, y, color=cm.cmap(frac), linewidth=0, zorder=zorder)


def plot_robustness_affinity_panel(letter, size):
    # sharey=True: both subplots put "Mean docking score" on the y-axis, so they're directly
    # comparable at a glance and share one set of y-ticks/limits (driven by the left subplot's
    # explicit [-15, -4] setting below).
    fig, axs = plt.subplots(1, 2, figsize=size, sharey=True, gridspec_kw={"wspace": 0.1})
    for ax in axs:
        stylize(ax)

    # Both axes have identical height (same subplot row), so one shared axes-fraction y-value
    # here puts both x-axis labels at the exact same figure height - matplotlib's automatic
    # label placement instead follows each axis's own tick-label height, and the right subplot's
    # log-scale "10^n" ticks are taller than the left's plain decimals, which otherwise leaves
    # the two "xlabel"s misaligned even though the axes themselves already line up.
    X_LABEL_Y_C = -0.16

    # Left: std across 5 replicates (script 65) vs. mean docking score, zoomed to the
    # score in [-15, -4] / std in [-0.1, 2.5] window - excludes the extreme positive-score
    # "invalid pose" outliers from the VIEW only (never dropped from the underlying data).
    robustness = load_replicate_robustness()
    r_pearson, _ = pearsonr(robustness["score"], robustness["score_std"])
    r_spearman, _ = spearmanr(robustness["score"], robustness["score_std"])
    print(f"\nPanel c (left) - replicate robustness: {len(robustness):,} points "
          f"(full, unfiltered): Pearson r={r_pearson:.4f}, Spearman r={r_spearman:.4f}")

    left_view = robustness[robustness["score"].between(-15, -4) & robustness["score_std"].between(-0.1, 2.5)]
    print(f"  In view [-15,-4] x [-0.1,2.5]: {len(left_view):,} / {len(robustness):,} points "
          f"({len(robustness) - len(left_view):,} excluded from view)")

    ax = axs[0]
    _density_scatter(ax, left_view["score_std"], left_view["score"])
    ax.set_xlim([-0.1, 2.5])
    ax.set_ylim([-15, -4])
    stylia.label(ax, xlabel="Standard deviation (kcal/mol)",
                 ylabel="Mean docking score (kcal/mol)")
    ax.xaxis.set_label_coords(0.5, X_LABEL_Y_C)

    # Right: Boltz-2-predicted IC50 (nM) vs. mean docking score, over the ~1,095 filtered hits -
    # a different, later-stage compound population than the left subplot (see module docstring).
    # Sharing the y-axis with the left subplot means points with score outside [-15, -4] (a few,
    # up towards 0) fall out of view here too - view-only, same non-destructive zoom as the left.
    dock_vs_aff = load_docking_vs_affinity()
    r_pearson, _ = pearsonr(dock_vs_aff["score"], dock_vs_aff["ic50_nM"])
    r_spearman, _ = spearmanr(dock_vs_aff["score"], dock_vs_aff["ic50_nM"])
    print(f"Panel c (right) - docking vs. Boltz-2: {len(dock_vs_aff):,} points "
          f"(full, unfiltered): Pearson r={r_pearson:.4f}, Spearman r={r_spearman:.4f}")

    right_view = dock_vs_aff[dock_vs_aff["score"] < 0]
    print(f"  In view (score < 0): {len(right_view):,} / {len(dock_vs_aff):,} points "
          f"({len(dock_vs_aff) - len(right_view):,} excluded from view)")

    # Low-confidence Boltz-2 binder calls (affinity_probability_binary < 0.5 - not to be confused
    # with the affinity_pred_value/IC50 itself) are drawn first, in flat gray, so they sit behind
    # the density-colored high-confidence points rather than competing with them for attention.
    low_conf = right_view[right_view["affinity_probability_binary"] < 0.5]
    high_conf = right_view[right_view["affinity_probability_binary"] >= 0.5]
    print(f"  Low-confidence (prob < 0.5, grayed out): {len(low_conf):,} / {len(right_view):,}")

    ax = axs[1]
    ax.scatter(low_conf["ic50_nM"], low_conf["score"], color="lightgray", linewidth=0, zorder=1)
    _density_scatter(ax, high_conf["ic50_nM"], high_conf["score"], log_x=True)
    ax.set_xscale("log")
    stylia.label(ax, xlabel="Boltz-2 predicted IC50 (nM)", ylabel="")
    ax.xaxis.set_label_coords(0.5, X_LABEL_Y_C)
    ax.tick_params(axis="y", length=0)

    save_panel(fig, letter, use_tight_layout=False,
               subplots_adjust=dict(left=0.10, right=0.98, top=0.97, bottom=0.20, wspace=0.12))


def plot_upset_panel(letter, size):
    """Panel d: script 68's own CAT (score <= -11) and NON-CAT (score <= -10) protein-level
    UpSet plots (see module docstring), composited side by side - same render-externally-then-
    imshow pattern as panel a's plot_structure_panel (autocrop_to_content reused as-is)."""
    missing = [p for p in (UPSET_CAT_PNG, UPSET_NONCAT_PNG) if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            "Panel d needs script 68's own UpSet renders, missing: " + ", ".join(missing) +
            ". Run `python 68_plot_results.py --only-upset` first.")

    fig, axs = plt.subplots(1, 2, figsize=size)
    fig.patch.set_facecolor("white")
    for ax, png_path in zip(axs, (UPSET_CAT_PNG, UPSET_NONCAT_PNG)):
        stylize(ax)
        img = autocrop_to_content(mpimg.imread(png_path))
        ax.imshow(img, interpolation="none")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
    stylia.label(axs[0], xlabel="", ylabel="")

    save_panel(fig, letter)


def main(rerun=False, subpanels=None):
    if subpanels is None:
        subpanels = PANEL_LETTERS
    sizes = load_panel_sizes(subpanels)

    if "a" in subpanels:
        plot_structure_panel("a", sizes["a"], rerun=rerun)
    if "b" in subpanels:
        plot_pocket_grid_panel("b", sizes["b"], rerun=rerun)
    if "c" in subpanels:
        plot_robustness_affinity_panel("c", sizes["c"])
    if "d" in subpanels:
        plot_upset_panel("d", sizes["d"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force PyMOL re-rendering of every panel's PNG(s), even if they "
                              "already exist (default: False, reusing whatever's already rendered)")
    parser.add_argument("--subpanels", type=str, default=None,
                         help=f"Comma-separated subset of panels to generate, from "
                              f"{{{','.join(PANEL_LETTERS)}}}. Default: all.")
    args = parser.parse_args()

    if args.subpanels is None:
        subpanels = PANEL_LETTERS
    else:
        subpanels = [p.strip() for p in args.subpanels.split(",")]
        invalid = [p for p in subpanels if p not in PANEL_LETTERS]
        if invalid:
            parser.error(f"Unknown panel(s) {invalid}, must be from {PANEL_LETTERS}")
    main(rerun=args.rerun, subpanels=subpanels)
