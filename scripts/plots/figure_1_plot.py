"""
Figure 1: protein structures (panel a) + sequence identity / structural RMSD heatmaps
(panels b/c), from figure_1_calculations.py. Panel d is a blank placeholder (previously a
PocketVec cosine-distance heatmap, dropped per request).

Built on figure_3/figure_4's panel-per-file architecture instead of the original,
oversized monolithic composite this script replaced: each panel is its own standalone
figure, saved as its own PDF (Fig_1a.pdf ... Fig_1d.pdf) at an EXACT physical size read
from output/plots/figure_1/panel_layout.csv (columns: panel, x, delta_x, y, delta_y,
padding, in cm) - so stylia's fixed-point-size text renders at its true intended size on
the page, instead of shrinking along with an oversized canvas once placed at final print
width. Panel letters are baked onto each saved PDF (add_panel_label, from figure_3) and a
Nature two-column-width guideline check (MAX_WIDTH_IN, also from figure_3) warns if any
panel's delta_x exceeds stylia.SIZE.

Merging into one positioned master PDF (Fig_1_full.pdf, via pypdf) happens in this same
file (merge_panels(), called at the end of main()) rather than a separate script.

panel_layout.csv is a first-draft starting point, meant to be tuned iteratively by
rendering and looking, not solved analytically up front.

Usage:
    python figure_1_plot.py [--rerun] [--subpanels a,b,c,d] [--show-grids]
"""
import argparse
import glob
import os
import sys

os.environ["QT_QPA_PLATFORM"] = "offscreen"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import json
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde
from pypdf import PdfReader, PdfWriter, Transformation
import pymupdf
import pymol
from pymol import cmd
import stylia
from stylia.config import get_fg_color
from stylia.figure.figure import stylize

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

data_dir = os.path.join(root, "..", "..", "data")
output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)

PANEL_LETTERS = ["a", "b", "c", "d"]
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

# ===========================================================================
# Panel a data + helpers: protein structures + VI/PDB/ChEMBL/pocket circles
# (ported verbatim from the original monolithic-composite version of this script)
# ===========================================================================

pocket_detection_data = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))

with open(os.path.join(plots_dir, "color_mapping.json")) as f:
    mappings = json.load(f)
uniprot_to_gene = mappings["uniprot_to_gene"]
cmap_dict = mappings["gene_to_color"]
gene_to_vi = mappings["gene_to_vulnerability_index"]
VI_GENOME_MIN = mappings["vulnerability_index_genome_min"]
VI_GENOME_MAX = mappings["vulnerability_index_genome_max"]
gene_to_pdb_count = mappings["gene_to_pdb_count"]
gene_to_chembl_binders = mappings["gene_to_chembl_binders"]
gene_to_chembl_total = mappings["gene_to_chembl_total"]
gene_to_unique_pocket_count = mappings["gene_to_unique_pocket_count"]


def zoom_to_fixed_box(selection="structure", box_size=30.0, box_name="_zoom_box"):
    (xmin, ymin, zmin), (xmax, ymax, zmax) = cmd.get_extent(selection)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
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


PROTEINS = sorted(set(pocket_detection_data['Uniprot AC']))
PROTEINS = np.array(PROTEINS)[np.argsort([uniprot_to_gene[i] for i in PROTEINS])]

MY_VIEWS = {
    "P9WFW7": (0.312430620,    0.466541409,    0.827481985, -0.380785346,   -0.736536384,    0.559036851,  0.870279729,   -0.489755034,   -0.052467018,  0.000000000,    0.000000000, -466.590057373, -1.569999695,    4.530498505,   -0.682498932,367.863159180,  565.316955566,   20.000000000 ),
    "P9WFW5": (0.5257094502449036, 0.24445058405399323, -0.8147862553596497, -0.4274793863296509, 0.9040127396583557, -0.004594181664288044, 0.7354535460472107, 0.35071954131126404, 0.5797445774078369, 0.0, 0.0, -491.1474304199219, 1.5335006713867188, 0.069000244140625, -2.2099990844726562, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFW3": (-0.26671209931373596, 0.51364666223526, 0.815495491027832, -0.777970016002655, -0.6141894459724426, 0.13241274654865265, 0.5688825249671936, -0.5991148352622986, 0.563413143157959, 4.291534423828125e-06, 8.13603401184082e-06, -374.25726318359375, 5.0578107833862305, 5.926779747009277, 4.082901954650879, 270.3341064453125, 478.18035888671875, 20.0),
    "P9WFW1": (-0.4216419458389282, 0.8592755198478699, -0.2895902097225189, -0.27892550826072693, 0.1809750199317932, 0.9431059956550598, 0.8627971410751343, 0.47842830419540405, 0.16336770355701447, 0.0, 0.0, -491.1474304199219, 4.153499603271484, -1.8014984130859375, 3.6859970092773438, 387.224365234375, 595.0704956054688, 20.0),
    "P9WQA1": ( 0.684969127,    0.360620856,    0.633064270, -0.128735647,   -0.795332789,    0.592346072,0.717109144,   -0.487235487,   -0.498353481,0.000000000,    0.000000000, -466.590057373,-1.930500031,   -0.642501831,   -2.583999634,367.863159180,  565.316955566,   20.000000000),
    "P9WN61": (0.762495219707489, 0.0944238156080246, -0.6400659680366516, -0.6421377658843994, -0.010534999892115593, -0.7665179967880249, -0.07912039756774902, 0.9954747557640076, 0.0525999590754509, 0.0, 0.0, -491.1474304199219, -1.5144996643066406, -0.6930007934570312, 6.256999969482422, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFV9": (-0.7302274703979492, -0.11305446177721024, 0.673784077167511, 0.1144353523850441, 0.9520406723022461, 0.2837631106376648, -0.673551619052887, 0.2843160629272461, -0.6822696328163147, 0.0, 0.0, -491.1474304199219, -0.9449996948242188, 8.261497497558594, -2.5309982299804688, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFV7": (0.634964108467102, 0.6521565318107605, -0.4141402542591095, 0.548267126083374, -0.002736924681812525, 0.8362982869148254, 0.5442638397216797, -0.7580783367156982, -0.35929417610168457, 0.0, 0.0, -491.1474304199219, 2.8590011596679688, -1.5230026245117188, -2.671001434326172, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFV5": (0.941296398639679, 0.1540936380624771, 0.30035918951034546, -0.33244943618774414, 0.2685926556587219, 0.9040654301643372, 0.05863688141107559, -0.9508472084999084, 0.3040543496608734, 0.0, 0.0, -491.1474304199219, 0.2584991455078125, 2.062999725341797, -3.6795005798339844, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFV3": (-0.0221111997961998, 0.9986055493354797, 0.04795076698064804, -0.21365615725517273, 0.04213476553559303, -0.9760009050369263, -0.9766592383384705, -0.03182630613446236, 0.2124253511428833, 0.0, 0.0, -497.0749816894531, 0.4790000915527344, -6.293498992919922, 10.033000946044922, 393.1518859863281, 600.998046875, 20.0),
    "P9WFV1": (-0.609279215335846, 0.7763136625289917, 0.16159172356128693, -0.13472627103328705, -0.3021676242351532, 0.943687379360199, 0.7814245223999023, 0.5531991124153137, 0.2886964678764343, 4.6759843826293945e-05, -2.6777386665344238e-05, -494.4760437011719, 10.133469581604004, 6.652966499328613, -5.180924415588379, 390.5534973144531, 598.3994140625, 20.0),
    "P9WFU9": (  0.108211756,   -0.602419376,    0.790810108, 0.695626080,    0.614177227,    0.372677028,-0.710205436,    0.509780943,    0.485519290, 0.000000000,    0.000000000, -466.590057373, 0.877498627,   -3.341499329,    6.893997192,67.863159180,  565.316955566,   20.000000000 ),
    "P9WFU5": (0.3144778311252594, 0.6041325926780701, 0.7322072386741638, 0.2726598381996155, -0.7963241934776306, 0.5399298667907715, 0.9092639088630676, 0.02984808385372162, -0.4151495099067688, 0.0, 0.0, -491.1474304199219, -3.6910018920898438, -1.782501220703125, -0.1999969482421875, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFU3": (0.8389383554458618, -0.3773128390312195, -0.3921956419944763, -0.49929937720298767, -0.24691034853458405, -0.8305028676986694, 0.21652206778526306, 0.8925665020942688, -0.3955335021018982, 0.0, 0.0, -491.1474304199219, 15.207504272460938, -7.323997497558594, -2.7709999084472656, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFU1": (0.649043918,   -0.306319863,   -0.696354270, -0.756129503,   -0.159016922,   -0.634806931, 0.083722532, 0.938553512, -0.334826291,  0.000000000, 0.000000000, -466.590057373, -0.643997192,    4.859497070,   -4.352996826,367.863159180,  565.316955566,   20.000000000 ),
    "P9WFT9": (0.8664358854293823, -0.32233402132987976, 0.381300687789917, 0.4621317684650421, 0.8068275451660156, -0.3680553734302521, -0.18900713324546814, 0.4951080083847046, 0.8480231165885925, 0.0, 0.0, -491.1474304199219, 1.9814987182617188, -6.362499237060547, 1.6375007629394531, 387.224365234375, 595.0704956054688, 20.0),
    "P9WFT7": (0.406122983,    0.380642772,    0.830766201, 0.126066312,    0.877084732,   -0.463493139, -0.905079246,    0.292967290,    0.308219969,  0.000000000,    0.000000000, -371.862823486, -4.729499817,   -3.773498535,    5.320499420,273.135955811,  470.589752197,   20.000000000),
    "P9WFT5": "",
    "P9WFT3": (0.244085222,   -0.812733531,   -0.529044509,-0.957535088,   -0.115654342,   -0.264106989, 0.153462872,    0.571041286,   -0.806449711, 0.000000000,    0.000000000, -466.590057373, 3.377502441,    1.659999847,   -0.513000488,67.863159180,  565.316955566,   20.000000000 ),
    "P9WFT1": "",
    "P9WFS9": ""
}

ZOOM_DICT = {
    "argS": 1.5, "alaS": 1.4, "aspS": 1.4, "gatA": 1.6, "gatB": 1.4, "thrS": 0.85,
    "ileS": 1.3, "valS": 1.0, "leuS": 1.3, "gltS": 1.4, "proS": 1.35, "glyS": 1.7,
    "serS": 1.5, "cysS1": 1.5, "metS": 1.5, "lysS": 1.5, "hisS": 1.6, "trpS": 1.7,
    "pheS": 1.25, "pheT": 1.1, "tyrS": 1.6,
}


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
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


def show_zoomed_image(ax, img, zoom=1.0):
    h, w = img.shape[:2]
    if zoom > 1:
        crop_w = int(w / zoom)
        crop_h = int(h / zoom)
        cx, cy = w // 2, h // 2
        x0 = max(cx - crop_w // 2, 0)
        x1 = min(cx + crop_w // 2, w)
        y0 = max(cy - crop_h // 2, 0)
        y1 = min(cy + crop_h // 2, h)
        img = img[y0:y1, x0:x1]

    img = autocrop_to_content(img)

    if zoom < 1:
        h2, w2 = img.shape[:2]
        pad_h = int(h2 / zoom) - h2
        pad_w = int(w2 / zoom) - w2
        pad_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
        pad_widths = [(pad_h // 2, pad_h - pad_h // 2), (pad_w // 2, pad_w - pad_w // 2)]
        if img.ndim == 3:
            pad_widths.append((0, 0))
        img = np.pad(img, pad_widths, mode="constant", constant_values=pad_val)

    ax.imshow(img, interpolation="none", resample=False)
    ax.set_aspect("equal", adjustable="datalim")


POCKET_SPHERE_RADIUS = 12.0


def render_protein_structures(proteins=PROTEINS, rerun=False):
    for protein in proteins:
        gene = uniprot_to_gene[protein]
        png_path = os.path.join(plots_dir, f"figure_{protein}_{gene}.png")
        if os.path.exists(png_path) and not rerun:
            print(f"Reusing existing render for {gene} ({protein}): {png_path}")
            continue

        df = pocket_detection_data[pocket_detection_data['Uniprot AC'] == protein].sort_values('Pocket score', ascending=False)
        ref_st = df['File name'].tolist()[0]

        ref_pockets = df[df['File name'] == ref_st].sort_values('Pocket probability', ascending=False)
        all_centroids = [
            np.array([float(v) for v in row['Pocket centroid coordinate (x y z)'].split()])
            for _, row in ref_pockets.iterrows()
        ]
        top_centroids = [all_centroids[0]]
        min_separation = 2 * POCKET_SPHERE_RADIUS + 1
        for centroid in all_centroids[1:]:
            if np.linalg.norm(centroid - top_centroids[0]) > min_separation:
                top_centroids.append(centroid)
                break

        pymol.finish_launching(['pymol', '-cq'])
        cmd.reinitialize()

        cmd.load(os.path.join(output_dir, "aligned_relaxed_structures", protein, ref_st), "structure")
        cmd.set_color("structure_color", [0.90, 0.90, 0.90])
        cmd.color("structure_color", "structure")
        cmd.hide("everything", "structure")
        cmd.show("surface", "structure")
        cmd.set("transparency", 0, "structure")

        centroid_color_rgb = mcolors.to_rgb(cmap_dict[gene])
        cmd.set_color("pocket_centroid_color", list(centroid_color_rgb))
        for i, (cx, cy, cz) in enumerate(top_centroids):
            obj_name = f"pocket_centroid_{i}"
            cmd.pseudoatom(obj_name, pos=[float(cx), float(cy), float(cz)], vdw=POCKET_SPHERE_RADIUS)
            cmd.show("spheres", obj_name)
            cmd.color("pocket_centroid_color", obj_name)
            cmd.set("sphere_transparency", 0.1, obj_name)

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
        if MY_VIEWS[protein] != "":
            cmd.set_view(MY_VIEWS[protein])
        zoom_to_fixed_box("structure", box_size=95)
        cmd.ray(1200, 1200)
        cmd.png(png_path, dpi=600)

        for i, centroid in enumerate(all_centroids):
            if any(np.array_equal(centroid, tc) for tc in top_centroids):
                continue
            cx, cy, cz = centroid
            obj_name = f"pocket_all_{i}"
            cmd.pseudoatom(obj_name, pos=[float(cx), float(cy), float(cz)], vdw=POCKET_SPHERE_RADIUS)
            cmd.show("spheres", obj_name)
            cmd.color("pocket_centroid_color", obj_name)
            cmd.set("sphere_transparency", 0.1, obj_name)

        cmd.save(os.path.join(plots_dir, f"session_{protein}_{gene}.pse"))


N_COLS = 7


def apply_grid_frame(ax, show_grids):
    if show_grids:
        ax.set_axis_on()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("red")
            spine.set_linewidth(0.5)
    else:
        ax.set_axis_off()


def render_structure_panel(ax, protein, show_grids=False):
    apply_grid_frame(ax, show_grids)
    gene = uniprot_to_gene[protein]
    file = os.path.join(plots_dir, f"figure_{protein}_{gene}.png")
    if os.path.exists(file):
        img = mpimg.imread(file)
        zoom = ZOOM_DICT[gene]
        show_zoomed_image(ax, img, zoom=zoom)
    stylia.label(ax, xlabel="", ylabel="")


def render_gene_label_row(ax, protein, show_grids=False):
    """Gene-name marker + text, in its own thin dedicated row - kept separate from both
    the structure image above (whose content varies per protein, e.g. pocket spheres near
    the bottom edge) and the VI/PDB/ChEMBL/pocket circle cluster below, so this label can
    never collide with either regardless of per-protein structure geometry or circle size."""
    apply_grid_frame(ax, show_grids)
    gene = uniprot_to_gene[protein]
    color = cmap_dict[gene]
    handle = [Line2D([0], [0], marker='o', color='w', label=gene, markerfacecolor=color,
                      markeredgecolor='black', markeredgewidth=0.5, markersize=stylia.MARKERSIZE)]
    # loc="center" (not the row-boundary-straddling bbox_to_anchor trick this used before the
    # gene label got its own dedicated row) - centers the marker+text as one unit under the
    # structure image above, safe now that this row has nothing else in it to collide with.
    ax.legend(handles=handle, loc="center", frameon=False, fontsize=stylia.FONTSIZE, handletextpad=0.3)
    stylia.label(ax, xlabel="", ylabel="")


PLACEHOLDER_CIRCLE_POSITIONS = [(0.15, 0.8), (0.85, 0.8), (0.15, 0.2), (0.85, 0.2)]  # [VI, PDBs, ChEMBL binders, distinct pockets]
CIRCLE_SIZE_MIN = 1
CIRCLE_SIZE_MAX = 300
VI_TRNA_MIN = min(gene_to_vi.values())
PDB_TRNA_MAX = max(gene_to_pdb_count.values())
CHEMBL_BINDERS_TRNA_MAX = max(gene_to_chembl_binders.values())
POCKET_COUNT_TRNA_MAX = max(gene_to_unique_pocket_count.values())
EMPTY_CIRCLE_SIZE = 10


def readable_text_color(hex_color):
    r, g, b = mcolors.to_rgb(hex_color)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if luminance < 0.5 else "black"


def render_placeholder_panel(ax, protein, show_grids=False):
    apply_grid_frame(ax, show_grids)
    gene = uniprot_to_gene[protein]
    color = cmap_dict[gene]
    s_max = stylia.MARKERSIZE_BIG * 20

    vi_x, vi_y = PLACEHOLDER_CIRCLE_POSITIONS[0]
    vi = gene_to_vi[gene]
    vi_frac = (VI_GENOME_MAX - vi) / (VI_GENOME_MAX - VI_TRNA_MIN)
    vi_size_value = CIRCLE_SIZE_MIN + vi_frac * (CIRCLE_SIZE_MAX - CIRCLE_SIZE_MIN)
    s = s_max * (vi_size_value / CIRCLE_SIZE_MAX)
    ax.scatter([vi_x], [vi_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
    ax.text(vi_x, vi_y, f"VI\n{vi:.1f}", ha='center', va='center', color=readable_text_color(color), fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)

    pdb_x, pdb_y = PLACEHOLDER_CIRCLE_POSITIONS[1]
    pdb_count = gene_to_pdb_count[gene]
    pdb_label = "PDB" if pdb_count == 1 else "PDBs"
    if pdb_count > 0:
        pdb_frac = 0.5 + 0.5 * (pdb_count - 1) / (PDB_TRNA_MAX - 1)
        s = s_max * pdb_frac
        ax.scatter([pdb_x], [pdb_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
        ax.text(pdb_x, pdb_y, f"{pdb_count}\n{pdb_label}", ha='center', va='center', color=readable_text_color(color), fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)
    else:
        ax.scatter([pdb_x], [pdb_y], s=EMPTY_CIRCLE_SIZE, facecolor='none', edgecolor=color, linewidth=stylia.LINEWIDTH, zorder=1)

    chembl_x, chembl_y = PLACEHOLDER_CIRCLE_POSITIONS[2]
    if gene in gene_to_chembl_binders:
        n_binders = gene_to_chembl_binders[gene]
        binder_label = "ligand" if n_binders == 1 else "ligands"
        if n_binders > 0:
            chembl_frac = 0.5 + 0.5 * (n_binders / CHEMBL_BINDERS_TRNA_MAX)
            s = s_max * chembl_frac
            ax.scatter([chembl_x], [chembl_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
        chembl_text_color = readable_text_color(color) if n_binders > 0 else "black"
        ax.text(chembl_x, chembl_y, f"{n_binders}\n{binder_label}", ha='center', va='center', color=chembl_text_color, fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)
    else:
        ax.scatter([chembl_x], [chembl_y], s=EMPTY_CIRCLE_SIZE, facecolor='none', edgecolor=color, linewidth=stylia.LINEWIDTH, zorder=1)

    pocket_x, pocket_y = PLACEHOLDER_CIRCLE_POSITIONS[3]
    n_pockets = gene_to_unique_pocket_count[gene]
    pocket_label = "pock."
    if n_pockets > 0:
        pocket_frac = 0.5 + 0.5 * (n_pockets - 1) / (POCKET_COUNT_TRNA_MAX - 1)
        s = s_max * pocket_frac
        ax.scatter([pocket_x], [pocket_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
    pocket_text_color = readable_text_color(color) if n_pockets > 0 else "black"
    ax.text(pocket_x, pocket_y, f"{n_pockets}\n{pocket_label}", ha='center', va='center', color=pocket_text_color, fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)

    AARS_CLASS_LABELS = {
        "alaS": "Class II", "argS": "Class I", "aspS": "Class II", "cysS1": "Class I",
        "gltS": "Class I", "glyS": "Class II", "hisS": "Class II", "ileS": "Class I",
        "leuS": "Class I", "lysS": "Class II", "metS": "Class I", "pheS": "Class II",
        "pheT": "Class II", "proS": "Class II", "serS": "Class II", "thrS": "Class II",
        "trpS": "Class I", "tyrS": "Class I", "valS": "Class I",
    }
    CLASS_CIRCLE_SIZE_FRAC = 0.8
    if gene in AARS_CLASS_LABELS:
        ax.scatter([0.5], [0.5], s=s_max * CLASS_CIRCLE_SIZE_FRAC, color=color, edgecolor='black', linewidth=0, zorder=1)
        ax.text(0.5, 0.5, AARS_CLASS_LABELS[gene], ha='center', va='center',
                color=readable_text_color(color), fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)
    else:
        ax.scatter([0.5], [0.5], s=EMPTY_CIRCLE_SIZE, facecolor='none', edgecolor=color, linewidth=stylia.LINEWIDTH, zorder=1)

    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-0.3, 1.3)
    ax.set_autoscale_on(False)
    stylia.label(ax, xlabel="", ylabel="")


def cleanup_intermediate_files(store_pymol, remove_pngs):
    if not store_pymol:
        for f in glob.glob(os.path.join(plots_dir, "session_*.pse")):
            os.remove(f)
    if remove_pngs:
        for f in glob.glob(os.path.join(plots_dir, "figure_*.png")):
            os.remove(f)


# ===========================================================================
# Panel b/c data + helper: SeqId / structural RMSD heatmaps
# (ported from the original monolithic-composite version; PocketVec dropped, abc-title handling dropped
# since panel-lettering is now centralized in save_panel/add_panel_label)
# ===========================================================================

seqid_df = pd.read_csv(os.path.join(plots_dir, "sequence_identity_matrix.tsv"), sep="\t", index_col=0)
alphabetical_order = seqid_df.index.tolist()

matrix = seqid_df.loc[alphabetical_order, alphabetical_order].values
gene_labels = [uniprot_to_gene[uid] for uid in alphabetical_order]

struct_df = pd.read_csv(os.path.join(plots_dir, "structural_similarity_matrix.tsv"), sep="\t", index_col=0)
struct_matrix = struct_df.loc[alphabetical_order, alphabetical_order].values

CBAR_VMIN = 20
CBAR_VMAX = 50
STRUCT_VMIN = 0
STRUCT_VMAX = 20
STRUCT_GRADIENT_MAX = 8


def clipped_gradient_cmap(cmap, vmin, vmax, gradient_max, n=256):
    frac = (gradient_max - vmin) / (vmax - vmin)
    n_gradient = max(1, int(round(n * frac)))
    colors_gradient = cmap(np.linspace(0, 1, n_gradient))
    end_color = cmap(1.0)
    colors_flat = np.tile(end_color, (n - n_gradient, 1))
    return mcolors.ListedColormap(np.vstack([colors_gradient, colors_flat]))


def plot_comparison_heatmap(ax, matrix, labels, cmap, vmin, vmax, cbar_label, ticks, line_color,
                             dist_values=None, overlay_values=None, dist_vmax=None, box_vmax=None,
                             box_ytick_exclude=()):
    cmap = cmap.copy()
    cmap.set_bad(color="#BBBBBB")
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    # Interleaved labels (not all 21 genes on both axes) - even index -> x only, odd index -> y
    # only, so every gene appears exactly once across the two axes, per request.
    x_labels = [lab if i % 2 == 0 else "" for i, lab in enumerate(labels)]
    y_labels = [lab if i % 2 == 1 else "" for i, lab in enumerate(labels)]
    ax.set_xticklabels(x_labels, rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    ax.set_yticklabels(y_labels, fontsize=stylia.FONTSIZE_SMALL)
    ax.yaxis.tick_right()

    divider = make_axes_locatable(ax)
    n = matrix.shape[0]

    tax = divider.append_axes("top", size="26%", pad=0.05)
    global_only_matrix = np.where(np.triu(np.ones((n, n), dtype=bool)), matrix, matrix.T)
    box_data = [np.delete(global_only_matrix[:, i], i) for i in range(n)]
    box_data = [d[~np.isnan(d)] for d in box_data]
    bp = tax.boxplot(
        box_data, positions=range(n), widths=0.6, whis=(0, 100), showfliers=False, patch_artist=True,
        boxprops=dict(edgecolor="black", linewidth=stylia.LINEWIDTH),
        whiskerprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        capprops=dict(color="black", linewidth=stylia.LINEWIDTH),
        medianprops=dict(color="black", linewidth=stylia.LINEWIDTH),
    )
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(cmap_dict[labels[i]])
        box.set_alpha(0.8)
    tax.set_xlim(-0.5, n - 0.5)
    box_ymax = box_vmax if box_vmax is not None else vmax
    tax.set_ylim(vmin, box_ymax)
    tax.set_yticks([t for t in range(0, int(box_ymax) + 1, 10) if t >= vmin and t not in box_ytick_exclude])
    tax.tick_params(axis="x", labelbottom=False, bottom=False)
    tax.yaxis.tick_right()
    tax.tick_params(axis="y", labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    for spine in ("top", "left"):
        tax.spines[spine].set_visible(False)
    tax.set_axisbelow(True)
    tax.grid(visible=True, linewidth=stylia.LINEWIDTH * 0.5, color="#DDDDDD", alpha=0.6)
    stylia.label(tax, xlabel="", ylabel="")

    cax = divider.append_axes("left", size="5%", pad=0.23)
    cbar = ax.figure.colorbar(im, cax=cax, ticklocation="right")
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    cbar.outline.set_linewidth(stylia.LINEWIDTH)

    dax = divider.append_axes("left", size="26%", pad=0.03, sharey=None if dist_vmax is not None else cax)
    pairwise_values = dist_values if dist_values is not None else matrix[np.triu_indices(n, k=1)]
    pairwise_values = pairwise_values[~np.isnan(pairwise_values)]
    y_max = dist_vmax if dist_vmax is not None else vmax
    kde = gaussian_kde(pairwise_values)
    y = np.linspace(vmin, y_max, 200)
    density = kde(y)
    dax.plot(-density, y, color=line_color, linewidth=stylia.LINEWIDTH)
    gradient = np.linspace(vmin, vmax, 256).reshape(-1, 1)
    grad_im = dax.imshow(gradient, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto",
                          origin="lower", extent=[-density.max(), 0, vmin, vmax])
    clip_vertices = np.column_stack([np.concatenate([-density, [0, 0]]), np.concatenate([y, [y[-1], y[0]]])])
    grad_im.set_clip_path(Polygon(clip_vertices, closed=True, transform=dax.transData))

    max_density = density.max()
    if overlay_values is not None:
        overlay_values = overlay_values[~np.isnan(overlay_values)]
        overlay_density = gaussian_kde(overlay_values)(y)
        dax.plot(-overlay_density, y, color="black", linewidth=stylia.LINEWIDTH * 0.5, linestyle="--")
        max_density = max(max_density, overlay_density.max())

    dax.set_xlim(-max_density * 1.02, 0)
    dax.set_xticks([])
    if dist_vmax is not None:
        dax.set_ylim(vmin, y_max)
        dax.tick_params(axis="y", labelleft=True, left=True, labelsize=stylia.FONTSIZE_SMALL,
                         length=2, width=stylia.LINEWIDTH)
        dax.spines["left"].set_linewidth(stylia.LINEWIDTH)
    else:
        dax.tick_params(axis="y", labelleft=False, left=False)
    for spine in ("top", "right") + (() if dist_vmax is not None else ("left",)):
        dax.spines[spine].set_visible(False)
    dax.set_axisbelow(True)
    dax.grid(visible=True, linewidth=stylia.LINEWIDTH * 0.5, color="#DDDDDD", alpha=0.6)
    stylia.label(dax, xlabel="", ylabel=cbar_label)

    stylia.label(ax, xlabel="", ylabel="")


# ===========================================================================
# Panel d data + helper: P2Rank probability curve + domain-annotation strips
# (migrated from figure_3_plot.py's own panel c, along with
# figure_1_calculations.py's compute_pocket_scores - figure 3's panel c no
# longer needs this data, so it was moved here rather than duplicated)
# ===========================================================================

# Vertical gap between the P2Rank prob. row and the 4-domain-row block below it -
# independent of MOSAIC's own (default) row/column gaps.
ROW_BLOCK_HSPACE = 0.4

# Vertical gap among just the 4 domain rows themselves, via their own inner subgridspec.
DOMAIN_ROWS_HSPACE = 0.45

# Row heights for the 5 stacked rows: P2Rank prob. gets 4 parts, each of the 4 domain
# bands gets 1 part (8 parts total) - so the probability curve occupies exactly
# 4/8 = 1/2 of the total stack height, versus an equal 1/5 split.
ROW_BLOCK_HEIGHT_RATIOS = [4, 1, 1, 1, 1]

# Fixed absolute left margin, wide enough to fit the longest row label ("Anticodon
# binding") at its own fixed point size - an absolute width (converted to a fraction of
# this panel's own delta_x in build_panel_d), not a fixed fraction, so it stays correct
# regardless of how panel d's own width is later tuned in panel_layout.csv.
ROW_LABEL_LEFT_MARGIN_IN = 0.75

# Extra headroom (axes-fraction) added on top of each label's own anchor y, so the label
# sits with real empty space above the plot box instead of flush against it.
LABEL_GAP = 0.02

# Typical P2Rank pocket-probability cutoff, per request - source/cite before this value
# appears in a manuscript legend (e.g. the P2Rank paper or its own documentation). Only
# used to stop the fill in plot_pocket_scores_row below the cutoff; no other visual
# change (no reference line, no tick relabeling).
P2RANK_CUTOFF = 0.2


def plot_pocket_scores_row(ax):
    """Panel d's 1st row (of 5, alongside the 4 plot_domain_strip_row bands) -
    pocket_rank on the x-axis (not y), probability on y. Only fills above P2RANK_CUTOFF -
    sub-cutoff pockets are left unfilled."""
    scores = pd.read_csv(os.path.join(plots_dir, "figure_1_pocket_scores.csv"))
    nc = stylia.NamedColors()
    ax.plot(scores["pocket_rank"], scores["pocket_probability"], color=nc.orchid, linewidth=stylia.LINEWIDTH, zorder=2)
    above_cutoff = scores["pocket_probability"].where(scores["pocket_probability"] >= P2RANK_CUTOFF)
    ax.fill_between(scores["pocket_rank"], above_cutoff, color=nc.orchid, alpha=0.3, zorder=1)
    ax.set_xlim(scores["pocket_rank"].min(), scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    # Explicit tick at 0.5 (unlabeled) in addition to the 0/1 endpoints - a visible
    # midpoint reference without cluttering the axis with a "0.5" label.
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "", "1"])
    stylia.label(ax, xlabel="", ylabel="P2Rank prob.")
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")
    # ha alone doesn't control distance from the axis - only how the text aligns around its own
    # anchor point, and this row's real "0"/"1" tick labels throw off matplotlib's automatic
    # anchor placement (Axis._autolabelpos) when combined with rotation=0, pushing the anchor (and
    # so the whole label) too far left. set_label_coords pins the anchor explicitly, close to the
    # tick labels, and disables that automatic repositioning so it stays put on every redraw.
    ax.yaxis.set_label_coords(-0.022, 0.5, transform=ax.transAxes)


# The 4 domain bands: column in figure_1_pocket_scores.csv -> band title -> its own
# "present" color, against a shared white "absent" background. Catalytic uses a
# ligand-evidence confidence threshold (CATALYTIC_CONFIDENCE_MIN); the other 3 use plain
# InterPro label presence, mutually exclusive with Catalytic - see
# figure_1_calculations.py's compute_pocket_scores/CURATED_LABEL_COLUMNS.
DOMAIN_STRIP_COLUMNS = [
    ("is_catalytic", "Catalytic", "crimson"),
    ("is_trna_binding", "tRNA binding", "cobalt"),
    ("is_editing", "Editing", "amber"),
    ("is_anticodon_binding", "Anticodon binding", "lime"),
]

# Each domain-row "present" mark's width, in pocket_rank units (1 unit = 1 pocket's own
# rank spacing) - 2x a pocket's natural 1-unit spacing, for better legibility.
DOMAIN_ROW_BAR_WIDTH = 2.0


def plot_domain_strip_row(ax, column, title, color_name):
    """One row of panel d (4 stacked rows total, 1 per DOMAIN_STRIP_COLUMNS entry) - a
    thin colored cell per pocket (all 276, same pocket_rank x-order - sorted by P2Rank
    probability descending, rank 1 leftmost), in this band's own color where `column` is
    True, white otherwise - figure_1_calculations.py's compute_pocket_scores, itself from
    output/77_pocket_annotation/pocket_detection_interpro_updated.csv."""
    scores = pd.read_csv(os.path.join(plots_dir, "figure_1_pocket_scores.csv"))
    nc = stylia.NamedColors()
    present_color = getattr(nc, color_name)
    colors = [present_color if v else nc.white for v in scores[column]]
    ax.bar(scores["pocket_rank"], 1, width=DOMAIN_ROW_BAR_WIDTH, color=colors, edgecolor="none", zorder=2)
    ax.set_xlim(scores["pocket_rank"].min(), scores["pocket_rank"].max())
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    stylia.label(ax, xlabel="", ylabel=title)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.label.set_ha("right")
    ax.yaxis.label.set_va("center")


# ===========================================================================
# Panel-saving scaffolding (ported from figure_3_plot.py / figure_4_plot.py)
# ===========================================================================

PANEL_LABEL_MARGIN = 0.02


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page), fixed regardless of
    padding - from figure_3_plot.py."""
    fig.text(PANEL_LABEL_MARGIN, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
              transform=fig.transFigure)


MAX_WIDTH_IN = stylia.SIZE  # Nature two-column guideline (~7.09in/18cm) - from figure_3_plot.py


def load_panel_sizes(panels):
    if not os.path.exists(panel_layout_path):
        raise FileNotFoundError(
            f"{panel_layout_path} not found. Expected columns: panel, x, delta_x, y, delta_y (cm), "
            "one row per panel: " + ", ".join(PANEL_LETTERS))
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: (df.loc[p, "delta_x"] / 2.54, df.loc[p, "delta_y"] / 2.54) for p in panels}


def load_panel_padding(panels):
    df = pd.read_csv(panel_layout_path).set_index("panel")
    if "padding" not in df.columns:
        return {p: 0.0 for p in panels}
    missing = [p for p in panels if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")
    return {p: df.loc[p, "padding"] for p in panels}


def apply_padding(fig, padding):
    scale = 1 - padding
    for ax in fig.axes:
        pos = ax.get_position()
        x0 = 0.5 + (pos.x0 - 0.5) * scale
        x1 = 0.5 + (pos.x1 - 0.5) * scale
        y0 = 0.5 + (pos.y0 - 0.5) * scale
        y1 = 0.5 + (pos.y1 - 0.5) * scale
        ax.set_position([x0, y0, x1 - x0, y1 - y0])


def save_panel(fig, letter, use_tight_layout=True, tight_pad=1.08, tight_w_pad=None,
               subplots_adjust=None, padding=0.0):
    """Saves at exactly panel_layout.csv's delta_x/delta_y (no bbox_inches="tight").
    use_tight_layout=False + subplots_adjust is the escape hatch for panels where
    tight_layout() warns/fails (e.g. panel a's packed grid) - from figure_4_plot.py.
    fig.set_tight_layout(False) defeats stylia's own auto-relayout rcParam, which would
    otherwise silently override this layout at savefig time."""
    fig.set_tight_layout(False)
    if use_tight_layout:
        plt.tight_layout(pad=tight_pad, w_pad=tight_w_pad)
    else:
        fig.subplots_adjust(**(subplots_adjust or dict(left=0.01, right=0.99, top=0.99, bottom=0.01)))
    apply_padding(fig, padding)
    add_panel_label(fig, letter)
    output_path = os.path.join(plots_dir, f"Fig_1{letter}.pdf")
    plt.savefig(output_path, dpi=600, transparent=False)
    plt.close(fig)
    print(f"Saved Fig_1{letter}.pdf")


CM_TO_PT = 72 / 2.54


def merge_panels():
    """Pastes each Fig_1{letter}.pdf onto one positioned master PDF (Fig_1_full.pdf) per
    panel_layout.csv's x/y - pure translation (no rescaling, true vector via pypdf), since
    each panel already saved at its exact target size. Ported from figure_4_merge.py, kept
    in this same file per request rather than split into a separate merge script."""
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS)

    writer = PdfWriter()
    master_page = writer.add_blank_page(width=total_width_cm * CM_TO_PT, height=total_height_cm * CM_TO_PT)

    for p in PANEL_LETTERS:
        panel_path = os.path.join(plots_dir, f"Fig_1{p}.pdf")
        if not os.path.exists(panel_path):
            print(f"  Skipping {p}: {panel_path} not found - Fig_1_full.pdf will be missing it.")
            continue
        panel_page = PdfReader(panel_path).pages[0]
        x_top_cm, y_top_cm, delta_y_cm = df.loc[p, "x"], df.loc[p, "y"], df.loc[p, "delta_y"]
        x_pt = x_top_cm * CM_TO_PT
        y_bottom_pt = (total_height_cm - y_top_cm - delta_y_cm) * CM_TO_PT
        master_page.merge_transformed_page(panel_page, Transformation().translate(tx=x_pt, ty=y_bottom_pt))
        print(f"  Placed Fig_1{p}.pdf at x={x_top_cm}cm, y={y_top_cm}cm (top-left)")

    output_path = os.path.join(plots_dir, "Fig_1_full.pdf")
    with open(output_path, "wb") as f:
        writer.write(f)
    print(f"Saved merged master figure ({total_width_cm:.2f} x {total_height_cm:.2f} cm) to {output_path}")

    # Flattened PNG alongside the vector PDF (user request), rendered at the same dpi=600
    # save_panel's own PDFs already use for their embedded raster content - pymupdf renders the
    # already-merged page directly (no external Poppler binary, unlike pdftoppm/pdf2image).
    png_path = os.path.join(plots_dir, "Fig_1_full.png")
    pdf_doc = pymupdf.open(output_path)
    zoom = 600 / 72  # pymupdf's Pixmap render defaults to 72dpi
    pix = pdf_doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    pix.save(png_path)
    pdf_doc.close()
    print(f"Saved {png_path}")


# ===========================================================================
# Panel builders
# ===========================================================================

GENE_LABEL_ROW_RATIO = 0.22  # thin dedicated row for the gene-name marker+text, own comment above render_gene_label_row


def build_panel_a(size, padding, show_grids=False):
    n_bands = len(PROTEINS) // N_COLS
    height_ratios = [1.0, GENE_LABEL_ROW_RATIO, 1.0] * n_bands

    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(3 * n_bands, N_COLS, height_ratios=height_ratios, hspace=0.0, wspace=0.0)
    cells = iter((row, col) for row in range(3 * n_bands) for col in range(N_COLS))

    for band_start in range(0, len(PROTEINS), N_COLS):
        band = PROTEINS[band_start:band_start + N_COLS]
        for protein in band:
            row, col = next(cells)
            render_structure_panel(stylize(fig.add_subplot(gs[row, col])), protein, show_grids=show_grids)
        for protein in band:
            row, col = next(cells)
            render_gene_label_row(stylize(fig.add_subplot(gs[row, col])), protein, show_grids=show_grids)
        for protein in band:
            row, col = next(cells)
            render_placeholder_panel(stylize(fig.add_subplot(gs[row, col])), protein, show_grids=show_grids)

    save_panel(fig, "a", use_tight_layout=False, padding=padding)


def build_panel_b(size, padding):
    nc = stylia.NamedColors()
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)

    seqid_cmap = stylia.FadingColormap("crimson", transformation=None).cmap
    seqid_n = matrix.shape[0]
    seqid_local_values = np.clip(matrix[np.tril_indices(seqid_n, k=-1)], CBAR_VMIN, CBAR_VMAX)
    seqid_global_values = np.clip(matrix[np.triu_indices(seqid_n, k=1)], CBAR_VMIN, CBAR_VMAX)
    plot_comparison_heatmap(ax, matrix, gene_labels, seqid_cmap, CBAR_VMIN, CBAR_VMAX,
                             "Sequence identity (%)", ticks=np.arange(CBAR_VMIN, CBAR_VMAX + 1, 5),
                             line_color=nc.crimson, dist_values=seqid_local_values,
                             overlay_values=seqid_global_values, box_vmax=40)
    save_panel(fig, "b", padding=padding)


def build_panel_c(size, padding):
    nc = stylia.NamedColors()
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)

    struct_cmap_full = stylia.FadingColormap("cobalt", transformation=None).cmap.reversed()
    struct_cmap = clipped_gradient_cmap(struct_cmap_full, STRUCT_VMIN, STRUCT_VMAX, STRUCT_GRADIENT_MAX)
    struct_n = struct_matrix.shape[0]
    struct_local_values = struct_matrix[np.tril_indices(struct_n, k=-1)]
    struct_global_values = struct_matrix[np.triu_indices(struct_n, k=1)]
    plot_comparison_heatmap(ax, struct_matrix, gene_labels, struct_cmap, STRUCT_VMIN, STRUCT_VMAX,
                             "Structural RMSD (Å)", ticks=np.arange(STRUCT_VMIN, STRUCT_VMAX + 1, 5),
                             line_color=nc.cobalt, dist_values=struct_local_values,
                             overlay_values=struct_global_values)
    save_panel(fig, "c", padding=padding)


def build_panel_d(size, padding):
    """P2Rank probability curve + 4 domain bands, migrated from figure_3_plot.py's own
    panel c (previously a blank placeholder here, before that a PocketVec cosine-distance
    heatmap). Outer gridspec splits P2Rank prob. from the 4-domain-row block
    (ROW_BLOCK_HSPACE gap, ROW_BLOCK_HEIGHT_RATIOS split); an inner subgridspec then packs
    the 4 domain rows with their own tighter DOMAIN_ROWS_HSPACE gap."""
    fig = plt.figure(figsize=size)
    fig.patch.set_facecolor("white")
    outer_gs = fig.add_gridspec(2, 1, height_ratios=[ROW_BLOCK_HEIGHT_RATIOS[0], sum(ROW_BLOCK_HEIGHT_RATIOS[1:])],
                                 hspace=ROW_BLOCK_HSPACE)
    plot_pocket_scores_row(stylize(fig.add_subplot(outer_gs[0, 0])))
    domain_gs = outer_gs[1, 0].subgridspec(len(DOMAIN_STRIP_COLUMNS), 1, hspace=DOMAIN_ROWS_HSPACE)
    for i, (column, title, color_name) in enumerate(DOMAIN_STRIP_COLUMNS):
        ax = stylize(fig.add_subplot(domain_gs[i, 0]))
        plot_domain_strip_row(ax, column=column, title=title, color_name=color_name)
    # tight_layout() isn't compatible with this nested subgridspec (falls back to
    # matplotlib's stock default margins, clipping the longest row labels on the left) -
    # figure_1_plot.py's own save_panel() exposes an explicit escape hatch for exactly
    # this, so the margin is set directly here rather than via a manual post-hoc
    # fig.subplots_adjust() call after save_panel() runs (figure_3_plot.py's own
    # workaround, needed there only because its simpler save_panel() lacks this
    # parameter). Left margin is a FIXED ABSOLUTE width (ROW_LABEL_LEFT_MARGIN_IN),
    # converted to a fraction of this panel's own delta_x - stays correct regardless of
    # how panel d's own width is tuned in panel_layout.csv.
    save_panel(fig, "d", use_tight_layout=False,
               subplots_adjust=dict(left=ROW_LABEL_LEFT_MARGIN_IN / size[0], right=0.995),
               padding=padding)


def main(rerun=False, subpanels=None, show_grids=False, store_pymol=False, remove_pngs=False):
    if subpanels is None:
        subpanels = PANEL_LETTERS
    sizes = load_panel_sizes(subpanels)
    paddings = load_panel_padding(subpanels)

    for letter in subpanels:
        if sizes[letter][0] > MAX_WIDTH_IN:
            print(f"WARNING: panel '{letter}' delta_x is {sizes[letter][0] * 2.54:.2f}cm, exceeding "
                  f"the {MAX_WIDTH_IN * 2.54:.2f}cm Nature two-column guideline by "
                  f"{(sizes[letter][0] - MAX_WIDTH_IN) * 2.54:.2f}cm.")

    if "a" in subpanels:
        render_protein_structures(proteins=PROTEINS, rerun=rerun)
        build_panel_a(sizes["a"], paddings["a"], show_grids=show_grids)
        cleanup_intermediate_files(store_pymol=store_pymol, remove_pngs=remove_pngs)

    if "b" in subpanels:
        build_panel_b(sizes["b"], paddings["b"])

    if "c" in subpanels:
        build_panel_c(sizes["c"], paddings["c"])

    if "d" in subpanels:
        build_panel_d(sizes["d"], paddings["d"])

    merge_panels()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False,
                         help="Force PyMOL re-rendering of panel a's per-protein PNGs, even if they already exist")
    parser.add_argument("--remove", action="store_true", default=False,
                         help="Delete panel a's per-protein figure_*.png renders after building it (default: keep them)")
    parser.add_argument("--store-pymol", action="store_true", default=False,
                         help="Keep panel a's per-protein PyMOL .pse sessions instead of deleting them after use")
    parser.add_argument("--show-grids", action="store_true", default=False,
                         help="Draw panel a's GridSpec rectangles (red spines) to visualize the underlying grid layout")
    parser.add_argument("--subpanels", type=str, default=None,
                         help="Comma-separated subset of panels to (re)generate (e.g. 'b,c'), from "
                              f"{{{','.join(PANEL_LETTERS)}}}. Default: all. Fig_1_full.pdf is always "
                              "re-merged from whatever Fig_1{letter}.pdf files currently exist on disk.")
    args = parser.parse_args()
    if args.subpanels is None:
        subpanels = PANEL_LETTERS
    else:
        subpanels = [p.strip() for p in args.subpanels.split(",")]
        invalid = [p for p in subpanels if p not in PANEL_LETTERS]
        if invalid:
            parser.error(f"Unknown panel(s) {invalid}, must be from {PANEL_LETTERS}")
    main(rerun=args.rerun, subpanels=subpanels, show_grids=args.show_grids,
         store_pymol=args.store_pymol, remove_pngs=args.remove)
