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
import pymol
from pymol import cmd
import stylia
from stylia.figure.figure import stylize

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

data_dir = os.path.join(root, "..", "..", "data")
output_dir = os.path.join(root, "..", "..", "output")
plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_1")
os.makedirs(plots_dir, exist_ok=True)

# ===========================================================================
# Panel a: protein structures + VI/PDB/ChEMBL/pocket circles
# ===========================================================================

# Load pocket detection data
pocket_detection_data = pd.read_csv(os.path.join(output_dir, "pocket_detection_data.csv"))

# Uniprot->gene and gene->color mappings computed by figure_1_calculations.py
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
        (cx - h, cy - h, cz - h),
        (cx - h, cy - h, cz + h),
        (cx - h, cy + h, cz - h),
        (cx - h, cy + h, cz + h),
        (cx + h, cy - h, cz - h),
        (cx + h, cy - h, cz + h),
        (cx + h, cy + h, cz - h),
        (cx + h, cy + h, cz + h),
    ]

    for i, (x, y, z) in enumerate(corners):
        cmd.pseudoatom(f"{box_name}_{i}", pos=[x, y, z])

    cmd.group(box_name, f"{box_name}_*")
    cmd.zoom(box_name, buffer=0.0, complete=1)

    # cleanup
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
    "argS": 1.5,
    "alaS": 1.4,
    "aspS": 1.4,
    "gatA": 1.6,
    "gatB": 1.4,
    "thrS": 0.85,
    "ileS": 1.3,
    "valS": 1.0,
    "leuS": 1.3,
    "gltS": 1.4,
    "proS": 1.35,
    "glyS": 1.7,
    "serS": 1.5,
    "cysS1": 1.5,
    "metS": 1.5,
    "lysS": 1.5,
    "hisS": 1.6,
    "trpS": 1.7,
    "pheS": 1.25,
    "pheT": 1.1,
    "tyrS": 1.6,
}


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
    # PyMOL renders onto a fixed-size white canvas, so the protein blob is
    # usually surrounded by a lot of true whitespace beyond any hand-tuned zoom.
    # Crop tightly to the non-white bounding box (plus a small margin) so that
    # whitespace doesn't show up as a gap in the composite grid.
    # mpimg.imread returns float32 in [0, 1] for PNGs (not 0-255), so the
    # background threshold must scale with the array's own dtype/range.
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
        # Mirror image of the zoom>1 crop-in above, but applied after autocrop (padding
        # before autocrop would just get stripped straight back off): pad the
        # already-tightly-cropped content with whitespace so its canvas grows by 1/zoom,
        # shrinking how much of the fixed grid cell the protein itself fills.
        h2, w2 = img.shape[:2]
        pad_h = int(h2 / zoom) - h2
        pad_w = int(w2 / zoom) - w2
        pad_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
        pad_widths = [(pad_h // 2, pad_h - pad_h // 2), (pad_w // 2, pad_w - pad_w // 2)]
        if img.ndim == 3:
            pad_widths.append((0, 0))
        img = np.pad(img, pad_widths, mode="constant", constant_values=pad_val)

    ax.imshow(img, interpolation="none", resample=False)
    # adjustable="datalim" keeps the axes' own box fixed at its GridSpec-assigned
    # rectangle and pads the data view instead — otherwise matplotlib's default
    # ("box") shrinks/repositions the box itself for extreme-aspect images, which
    # throws off anything anchored via ax.transAxes (e.g. the gene-name label).
    ax.set_aspect("equal", adjustable="datalim")


POCKET_SPHERE_RADIUS = 12.0


def render_protein_structures(proteins=PROTEINS, rerun=False):
    for protein in proteins:
        gene = uniprot_to_gene[protein]
        png_path = os.path.join(plots_dir, f"figure_{protein}_{gene}.png")
        if os.path.exists(png_path) and not rerun:
            print(f"Reusing existing render for {gene} ({protein}): {png_path}")
            continue

        # Filter pockets
        df = pocket_detection_data[pocket_detection_data['Uniprot AC'] == protein].sort_values('Pocket score', ascending=False)

        # Select reference structure
        ref_st = df['File name'].tolist()[0]

        # Pocket centroids to display: highest-probability pocket among this reference
        # structure's own pockets (same coordinate frame as the loaded structure), plus
        # the next-highest-probability pocket far enough away not to overlap it (distance
        # > 2*POCKET_SPHERE_RADIUS + 1, i.e. the two spheres' radii plus a small buffer).
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

        # Create pymol session
        pymol.finish_launching(['pymol', '-cq'])
        cmd.reinitialize()

        # Load structure — surface only, no pockets
        cmd.load(os.path.join(output_dir, "aligned_relaxed_structures", protein, ref_st), "structure")
        cmd.set_color("structure_color", [0.90, 0.90, 0.90])
        cmd.color("structure_color", "structure")
        cmd.hide("everything", "structure")
        cmd.show("surface", "structure")
        cmd.set("transparency", 0, "structure")

        # Pocket centroid markers: top pocket(s), POCKET_SPHERE_RADIUS sphere each,
        # protein's own color, almost opaque. PyMOL's sphere_transparency runs 0
        # (opaque) -> 1 (invisible), so "almost opaque" is a LOW value here, not 0.9.
        centroid_color_rgb = mcolors.to_rgb(cmap_dict[gene])
        cmd.set_color("pocket_centroid_color", list(centroid_color_rgb))
        for i, (cx, cy, cz) in enumerate(top_centroids):
            obj_name = f"pocket_centroid_{i}"
            cmd.pseudoatom(obj_name, pos=[float(cx), float(cy), float(cz)], vdw=POCKET_SPHERE_RADIUS)
            cmd.show("spheres", obj_name)
            cmd.color("pocket_centroid_color", obj_name)
            cmd.set("sphere_transparency", 0.1, obj_name)

        # Fancy visualization
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

        # After saving the PNG, add every other detected pocket for this reference
        # structure to the session (the PNG only shows the top pocket(s) selected above).
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
    # --show-grids draws each subplot's actual GridSpec-assigned rectangle (spines
    # only, no ticks/labels) so the underlying grid layout can be inspected directly.
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
        color = cmap_dict[gene]
        # Matches notebook 46's approach: ax.legend() with a small bbox_to_anchor,
        # leaving more of the row's height to the structure image itself. Safe now
        # that show_zoomed_image uses adjustable="datalim" (axes box position no
        # longer shifts per-column based on image content).
        legend_handles = [Line2D([0], [0], marker='o', color='w', label=gene, markerfacecolor=color, markeredgecolor='black', markeredgewidth=0.5, markersize=stylia.MARKERSIZE)]
        ax.legend(handles=legend_handles, loc="upper center", frameon=False, bbox_to_anchor=(0.5, 0.05), fontsize=stylia.FONTSIZE, handletextpad=0.3)

    stylia.label(ax, xlabel="", ylabel="")


PLACEHOLDER_CIRCLE_POSITIONS = [(0.3, 0.8), (0.7, 0.8), (0.3, 0.2), (0.7, 0.2)]  # [VI, PDBs, ChEMBL binders, distinct pockets]

# Shared size scale for every circle in the panel: rendered marker size always maps
# linearly onto [CIRCLE_SIZE_MIN, CIRCLE_SIZE_MAX], whatever the underlying metric.
CIRCLE_SIZE_MIN = 1
CIRCLE_SIZE_MAX = 300

# VI: size CIRCLE_SIZE_MIN at the least vulnerable gene in the whole H37Rv proteome
# (VI_GENOME_MAX), size CIRCLE_SIZE_MAX at the most vulnerable gene among these 21
# tRNA synthetases (VI_TRNA_MIN, not genome-wide) — more negative VI = more vulnerable.
VI_TRNA_MIN = min(gene_to_vi.values())

# PDBs: size CIRCLE_SIZE_MIN at 0 associated experimental structures, size
# CIRCLE_SIZE_MAX at the most PDB-covered of these 21 tRNA synthetases.
PDB_TRNA_MAX = max(gene_to_pdb_count.values())

# ChEMBL binders: size CIRCLE_SIZE_MIN at 0 confirmed binders, size CIRCLE_SIZE_MAX at
# the most-binder-rich of the (only 3) genes with any ChEMBL target mapping at all.
CHEMBL_BINDERS_TRNA_MAX = max(gene_to_chembl_binders.values())

# Distinct pockets: size CIRCLE_SIZE_MIN at the fewest distinct pockets, size
# CIRCLE_SIZE_MAX at the most, among these 21 tRNA synthetases.
POCKET_COUNT_TRNA_MAX = max(gene_to_unique_pocket_count.values())

# Fixed placeholder size (points^2, absolute — not scaled by s_max) for "no data" circles
# (0 PDBs, no ChEMBL mapping): deliberately small and constant across every panel, so it
# reads as "nothing here" rather than as a data point on the row's own size scale.
EMPTY_CIRCLE_SIZE = 10


def readable_text_color(hex_color):
    # Standard YIQ perceptual luminance; below 0.5 the background reads as "dark"
    # enough that black text loses contrast, so switch to white.
    r, g, b = mcolors.to_rgb(hex_color)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if luminance < 0.5 else "black"


def render_placeholder_panel(ax, protein, show_grids=False):
    apply_grid_frame(ax, show_grids)

    gene = uniprot_to_gene[protein]
    color = cmap_dict[gene]

    # s_max is an empirically-validated points^2 area (not derived from
    # ax.transData): the axes box isn't finalized until stylia.save_figure()'s
    # internal tight_layout() call, which runs after this function, so any
    # live geometry read here (e.g. via transData) reflects stale, pre-layout
    # axes bounds. This value matches the fixed-size circles validated earlier
    # at this same layout to not overlap.
    s_max = stylia.MARKERSIZE_BIG * 20

    # Upper-left circle: Bosch et al. 2021 Vulnerability Index (VI).
    vi_x, vi_y = PLACEHOLDER_CIRCLE_POSITIONS[0]
    vi = gene_to_vi[gene]
    vi_frac = (VI_GENOME_MAX - vi) / (VI_GENOME_MAX - VI_TRNA_MIN)
    vi_size_value = CIRCLE_SIZE_MIN + vi_frac * (CIRCLE_SIZE_MAX - CIRCLE_SIZE_MIN)
    s = s_max * (vi_size_value / CIRCLE_SIZE_MAX)
    ax.scatter([vi_x], [vi_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
    ax.text(vi_x, vi_y, f"VI\n{vi:.1f}", ha='center', va='center', color=readable_text_color(color), fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)

    # Upper-right circle: number of associated experimental PDB structures. Zero PDBs
    # gets an empty (unfilled) placeholder circle and no text label, instead of a
    # scaled marker with a "0 PDBs" label.
    pdb_x, pdb_y = PLACEHOLDER_CIRCLE_POSITIONS[1]
    pdb_count = gene_to_pdb_count[gene]
    pdb_label = "PDB" if pdb_count == 1 else "PDBs"
    if pdb_count > 0:
        # 1 PDB = 50% of max size, then scales linearly up to 100% at PDB_TRNA_MAX
        # (mirrors the ChEMBL ligands / distinct pockets circles' floor-then-scale rule).
        pdb_frac = 0.5 + 0.5 * (pdb_count - 1) / (PDB_TRNA_MAX - 1)
        s = s_max * pdb_frac
        ax.scatter([pdb_x], [pdb_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
        ax.text(pdb_x, pdb_y, f"{pdb_count}\n{pdb_label}", ha='center', va='center', color=readable_text_color(color), fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)
    else:
        ax.scatter([pdb_x], [pdb_y], s=EMPTY_CIRCLE_SIZE, facecolor='none', edgecolor=color, linewidth=stylia.LINEWIDTH, zorder=1)

    # Bottom-left circle: ChEMBL ligands (IC50 <= 10 uM). Genes with no ChEMBL target
    # mapping at all get an empty (unfilled) placeholder circle and no text label,
    # instead of "0 ligands" (a real measured zero, drawn with no circle, stays as-is).
    chembl_x, chembl_y = PLACEHOLDER_CIRCLE_POSITIONS[2]
    if gene in gene_to_chembl_binders:
        n_binders = gene_to_chembl_binders[gene]
        binder_label = "ligand" if n_binders == 1 else "ligands"
        if n_binders > 0:
            # Any nonzero count starts at 50% of max size, then scales proportionally
            # up to 100% at CHEMBL_BINDERS_TRNA_MAX (unlike VI/PDBs, which scale from 0).
            chembl_frac = 0.5 + 0.5 * (n_binders / CHEMBL_BINDERS_TRNA_MAX)
            s = s_max * chembl_frac
            ax.scatter([chembl_x], [chembl_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
        chembl_text_color = readable_text_color(color) if n_binders > 0 else "black"
        ax.text(chembl_x, chembl_y, f"{n_binders}\n{binder_label}", ha='center', va='center', color=chembl_text_color, fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)
    else:
        ax.scatter([chembl_x], [chembl_y], s=EMPTY_CIRCLE_SIZE, facecolor='none', edgecolor=color, linewidth=stylia.LINEWIDTH, zorder=1)

    # Bottom-right circle: distinct pockets per protein, deduplicated across all of its
    # structures (see figure_1_calculations.py: greedy dedup by pocket centroid distance,
    # threshold matching notebooks/08_coherence_detected_pockets.ipynb).
    pocket_x, pocket_y = PLACEHOLDER_CIRCLE_POSITIONS[3]
    n_pockets = gene_to_unique_pocket_count[gene]
    pocket_label = "pock."
    if n_pockets > 0:
        # 1 pocket = 50% of max size, then scales linearly up to 100% at
        # POCKET_COUNT_TRNA_MAX (mirrors the ChEMBL binders circle's floor-then-scale rule).
        pocket_frac = 0.5 + 0.5 * (n_pockets - 1) / (POCKET_COUNT_TRNA_MAX - 1)
        s = s_max * pocket_frac
        ax.scatter([pocket_x], [pocket_y], s=s, color=color, edgecolor='black', linewidth=0, zorder=1)
    pocket_text_color = readable_text_color(color) if n_pockets > 0 else "black"
    ax.text(pocket_x, pocket_y, f"{n_pockets}\n{pocket_label}", ha='center', va='center', color=pocket_text_color, fontsize=stylia.FONTSIZE_SMALL, zorder=2, linespacing=1.6)

    # Set (and lock) fixed data coordinates AFTER all scatter/text calls — each
    # ax.scatter() call above triggers matplotlib's default autoscale-on-new-artist
    # behavior, which silently overrides limits set beforehand.
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(-0.3, 1.3)
    ax.set_autoscale_on(False)

    stylia.label(ax, xlabel="", ylabel="")


def cleanup_intermediate_files(store_pymol, remove_pngs):
    # session_*.pse and per-protein figure_*.png are written unconditionally by
    # render_protein_structures; drop them here if the caller asked to remove them.
    # (figure_1.png/.pdf, the composite, isn't matched by either pattern.)
    if not store_pymol:
        for f in glob.glob(os.path.join(plots_dir, "session_*.pse")):
            os.remove(f)
    if remove_pngs:
        for f in glob.glob(os.path.join(plots_dir, "figure_*.png")):
            os.remove(f)


# ===========================================================================
# Panels B/C/D: SeqId, structural RMSD, and PocketVec comparison heatmaps
# ===========================================================================

genes_sorted = sorted(uniprot_to_gene.values())

# scripts/plots/figure_1_calculations.py's output: pairwise sequence identity — upper
# triangle = global (Needleman-Wunsch, scripts/14_calculate_SeqId.py's SeqId_matrix.tsv),
# lower triangle = local (Smith-Waterman, that same script's LocalSeqId_matrix.tsv),
# diagonal = self-identity (~100%, meaningfully computed by both aligners already).
seqid_df = pd.read_csv(os.path.join(plots_dir, "sequence_identity_matrix.tsv"), sep="\t", index_col=0)
alphabetical_order = seqid_df.index.tolist()  # already alphabetical by gene name

matrix = seqid_df.loc[alphabetical_order, alphabetical_order].values
gene_labels = [uniprot_to_gene[uid] for uid in alphabetical_order]

# scripts/plots/figure_1_calculations.py's output: pairwise structural comparison —
# upper triangle = mean RMSD, lower triangle = min RMSD (both pooled across all
# structure-file combinations for that protein pair), diagonal = placeholder 0 (same-
# protein comparison not yet computed). Same alphabetical order as the SeqId matrix
# above, so rows/columns line up between the two panels.
struct_df = pd.read_csv(os.path.join(plots_dir, "structural_similarity_matrix.tsv"), sep="\t", index_col=0)
struct_matrix = struct_df.loc[alphabetical_order, alphabetical_order].values

# scripts/plots/figure_1_calculations.py's output: pocket-level PocketVec comparison,
# aggregated to protein level as min cosine distance between every deduplicated pocket
# of protein i and protein j (symmetric; diagonal = min distance between a protein's
# own distinct pockets). Same alphabetical order as the other panels.
pocketvec_df = pd.read_csv(os.path.join(plots_dir, "pocketvec_similarity_matrix.tsv"), sep="\t", index_col=0)
pocketvec_matrix = pocketvec_df.loc[alphabetical_order, alphabetical_order].values

CBAR_VMIN = 20
CBAR_VMAX = 40  # values above 40 (up to ~67% for local identity) saturate
STRUCT_VMIN = 0
STRUCT_VMAX = 20  # axis/colorbar/distribution extent
STRUCT_GRADIENT_MAX = 8  # local RMSD's true max (7.79) - full color gradient covers all of local, flat white beyond
POCKETVEC_VMIN = 0.12
POCKETVEC_VMAX = 0.22  # covers all but 1.4% of pairs (true max ~0.225)


def clipped_gradient_cmap(cmap, vmin, vmax, gradient_max, n=256):
    """Compress cmap's full gradient into [vmin, gradient_max]; hold its end color flat for (gradient_max, vmax]."""
    frac = (gradient_max - vmin) / (vmax - vmin)
    n_gradient = max(1, int(round(n * frac)))
    colors_gradient = cmap(np.linspace(0, 1, n_gradient))
    end_color = cmap(1.0)
    colors_flat = np.tile(end_color, (n - n_gradient, 1))
    return mcolors.ListedColormap(np.vstack([colors_gradient, colors_flat]))


def plot_comparison_heatmap(ax, matrix, labels, cmap, vmin, vmax, cbar_label, ticks, line_color, abc=None,
                             dist_values=None, overlay_values=None, dist_vmax=None):
    # Copy (not mutate) the caller's colormap - some pairs may have no reliable value
    # (e.g. panel c's 10%-coverage-filtered global RMSD, NaN for pairs where nothing
    # clears the cutoff) and should render as a distinct grey "no data" cell rather than
    # being confusable with a real value at the colormap's pale end.
    cmap = cmap.copy()
    cmap.set_bad(color="#BBBBBB")
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=stylia.FONTSIZE_SMALL)
    ax.set_yticklabels(labels, fontsize=stylia.FONTSIZE_SMALL)
    # Row labels move to the right side of the heatmap — the left side is now taken up
    # by the distribution + colorbar panels appended below.
    ax.yaxis.tick_right()

    divider = make_axes_locatable(ax)
    n = matrix.shape[0]

    # Top panel: one boxplot per protein, column-aligned with the heatmap, each showing
    # that protein's pairwise comparison values against the other proteins (self-
    # comparison excluded). whis=(0, 100) stretches whiskers to the actual min/max (no
    # outlier detection/fliers). No sharex here — boxplot()'s own tick-labeling clobbers
    # a shared axis's labels, so alignment is matched explicitly via xlim instead
    # (imshow's default extent for n categories is exactly [-0.5, n-0.5]).
    tax = divider.append_axes("top", size="20%", pad=0.05)
    # NaN-safe (e.g. panel c's 4 pairs with no coverage-filtered global value) - dropped
    # rather than plotted, so the boxplot summarizes only the pairs with real data.
    box_data = [np.delete(matrix[:, i], i) for i in range(n)]
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
    tax.set_ylim(vmin, vmax)
    # vmin excluded here — sitting right on the axis's own bottom border, it's redundant
    # and crowds the boxplot; the colorbar (below) still shows the full tick set.
    tax.set_yticks([t for t in ticks if t != vmin])
    tax.tick_params(axis="x", labelbottom=False, bottom=False)
    tax.tick_params(axis="y", labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    for spine in ("top", "right"):
        tax.spines[spine].set_visible(False)
    # stylia.create_figure() applies this same thin/light grid to ax automatically, but
    # axes appended manually via make_axes_locatable don't get it — set it explicitly so
    # tax/lax match the heatmap's own grid instead of falling back to matplotlib's bolder default.
    tax.set_axisbelow(True)
    tax.grid(visible=True, linewidth=stylia.LINEWIDTH * 0.5, color="#DDDDDD", alpha=0.6)
    stylia.label(tax, xlabel="", ylabel="")

    # Colorbar axis carved out of the same divider as the heatmap axis (rather than
    # fig.colorbar's default shrink=<1 sizing), so its height matches the heatmap's
    # height/ylim exactly instead of floating shorter in the middle of the panel.
    # Appended first (before the distribution panel) so it's the one adjacent to the
    # heatmap; ticks face right, into the open gap left by ax's labels moving right.
    # pad is wide enough to fit the widest tick labels used across all three panels
    # (PocketVec's "0.20"-style decimals), kept the same for all three for consistency.
    cax = divider.append_axes("left", size="5%", pad=0.28)
    cbar = ax.figure.colorbar(im, cax=cax, ticklocation="right")
    cbar.set_ticks(ticks)
    # Label lives on dax (further out) instead of here, so it reads as a single shared
    # label for the colorbar+distribution pair rather than being squeezed next to cax.
    # Explicit length/width — the default colorbar tick size renders disproportionately
    # thick/blocky at 600 dpi compared to the rest of the figure's line weight.
    cbar.ax.tick_params(labelsize=stylia.FONTSIZE_SMALL, length=2, width=stylia.LINEWIDTH)
    cbar.outline.set_linewidth(stylia.LINEWIDTH)

    # Distribution panel: single KDE of pairwise comparison values across all proteins
    # (upper-triangle values, i.e. one value per protein pair — excludes the diagonal
    # self-comparison). Appended after the colorbar, so it sits further left, outermost.
    # sharey=cax ties its value axis directly to the colorbar's own axis — no second,
    # redundant set of tick labels for the same quantity. Not used when dist_vmax is
    # given (the density panel then needs its own, wider range - see below).
    dax = divider.append_axes("left", size="20%", pad=0.03, sharey=None if dist_vmax is not None else cax)
    # dist_values overrides the default (upper-triangle) source for the main filled
    # curve - e.g. panel b passes the lower triangle (local identity) here instead.
    pairwise_values = dist_values if dist_values is not None else matrix[np.triu_indices(n, k=1)]
    pairwise_values = pairwise_values[~np.isnan(pairwise_values)]  # e.g. panel c's 4 no-data pairs
    # dist_vmax lets the density panel's own y-range extend past the heatmap/colorbar's
    # vmax - needed when the plotted quantity's full distribution doesn't fit the
    # heatmap's (deliberately narrower, for cell contrast) color range: without this, a
    # KDE evaluated only up to vmax would show a curve still rising at the edge rather
    # than its true peak-and-decline shape, and would look artificially truncated.
    y_max = dist_vmax if dist_vmax is not None else vmax
    kde = gaussian_kde(pairwise_values)
    y = np.linspace(vmin, y_max, 200)
    density = kde(y)
    # Baseline (zero density) sits against the colorbar side (dax's right edge); the
    # curve bulges away from it, toward the outer margin (dax's left edge).
    dax.plot(-density, y, color=line_color, linewidth=stylia.LINEWIDTH)
    # Fill under the curve with the same colormap/scale as the heatmap's colorbar
    # instead of a flat color, so the fill's color at a given height directly mirrors
    # what that value looks like on the heatmap/colorbar. Extent intentionally stays
    # vmin/vmax (not y_max) even when dist_vmax is set, so the fill only covers the
    # heatmap's own color-mapped range - the curve *line* above that (if any) is left
    # unfilled, honestly showing it's outside what the heatmap itself can represent.
    gradient = np.linspace(vmin, vmax, 256).reshape(-1, 1)
    grad_im = dax.imshow(gradient, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto",
                          origin="lower", extent=[-density.max(), 0, vmin, vmax])
    clip_vertices = np.column_stack([np.concatenate([-density, [0, 0]]), np.concatenate([y, [y[-1], y[0]]])])
    grad_im.set_clip_path(Polygon(clip_vertices, closed=True, transform=dax.transData))

    # Optional second distribution (e.g. global identity/RMSD, overlaid on the main
    # curve) as a thin dashed BLACK line on the same axis - no fill, just a lightweight
    # reference for comparison against the main curve.
    max_density = density.max()
    if overlay_values is not None:
        overlay_values = overlay_values[~np.isnan(overlay_values)]
        overlay_density = gaussian_kde(overlay_values)(y)
        dax.plot(-overlay_density, y, color="black", linewidth=stylia.LINEWIDTH * 0.5, linestyle="--")
        max_density = max(max_density, overlay_density.max())

    dax.set_xlim(-max_density * 1.02, 0)
    dax.set_xticks([])
    if dist_vmax is not None:
        # dax now has its own scale, different from the heatmap/colorbar's - needs its
        # own visible tick labels (previously relied entirely on cax's, right next to
        # it, since both shared one scale) so the reader can see the two panels differ.
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
    # abc lives here (upper-left of dax, the distribution panel) rather than tax/ax.
    # Extra pad (beyond stylia.label's default) lifts it further above dax's own top edge.
    stylia.label(dax, xlabel="", ylabel=cbar_label)
    if abc is not None:
        dax.set_title(abc, loc="left", fontweight="bold", fontsize=stylia.FONTSIZE_BIG, pad=22)

    stylia.label(ax, xlabel="", ylabel="")
    return dax


# ===========================================================================
# Master composite: structures (panel a, left) + heatmaps (panels b/c/d, right)
# ===========================================================================

ROW_HEIGHT_SCALE = 0.9  # shrinks unused vertical margin in width-constrained structure/circle rows

# Fractions of stylia.get_size(). Structures widened beyond its original standalone
# GRID_WIDTH=1.0 — at 1.0, each of the 7 structure columns is too narrow for the VI/PDB/
# ligand/pocket circles to sit side by side without touching (circle diameter ~0.38in
# vs ~0.18in of horizontal clearance at width=1.0). Heatmap stack keeps its original
# standalone proportions (width=0.5, height=1.5).
STRUCT_WIDTH_FRAC = 1.4
HEATMAP_WIDTH_FRAC = 0.5
HEATMAP_HEIGHT_FRAC = 1.5


def render_master_figure(proteins=PROTEINS, show_grids=False):
    n_bands = len(proteins) // N_COLS
    full_n_bands = len(PROTEINS) // N_COLS
    height_ratios = [1.0, 1.0] * n_bands
    struct_height_frac = sum(height_ratios) / (2 * full_n_bands) * ROW_HEIGHT_SCALE

    size = stylia.get_size()
    master_width_frac = STRUCT_WIDTH_FRAC + HEATMAP_WIDTH_FRAC
    master_height_frac = max(struct_height_frac, HEATMAP_HEIGHT_FRAC)

    fig = plt.figure(figsize=(master_width_frac * size, master_height_frac * size))
    fig.patch.set_facecolor("white")
    master_gs = fig.add_gridspec(1, 2, width_ratios=[STRUCT_WIDTH_FRAC, HEATMAP_WIDTH_FRAC], wspace=0.15)

    # --- Left: structures + circles grid (panel a) ---
    left_gs = master_gs[0, 0].subgridspec(2 * n_bands, N_COLS, height_ratios=height_ratios, hspace=0.0, wspace=0.0)
    first_structure_ax = None
    left_cells = iter((row, col) for row in range(2 * n_bands) for col in range(N_COLS))

    for band_start in range(0, len(proteins), N_COLS):
        band = proteins[band_start:band_start + N_COLS]

        for protein in band:
            row, col = next(left_cells)
            ax = stylize(fig.add_subplot(left_gs[row, col]))
            if first_structure_ax is None:
                first_structure_ax = ax
            render_structure_panel(ax, protein, show_grids=show_grids)

        for protein in band:
            row, col = next(left_cells)
            ax = stylize(fig.add_subplot(left_gs[row, col]))
            render_placeholder_panel(ax, protein, show_grids=show_grids)

    stylia.label(first_structure_ax, xlabel="", ylabel="")

    # --- Right: SeqId / structural RMSD / PocketVec heatmaps (panels b/c/d) ---
    right_gs = master_gs[0, 1].subgridspec(3, 1, hspace=0.0)
    nc = stylia.NamedColors()

    seqid_cmap = stylia.FadingColormap("crimson", transformation=None).cmap
    seqid_n = matrix.shape[0]
    # Main filled curve = local identity (lower triangle); thin dashed overlay = global
    # identity (upper triangle), for direct comparison on the same axis. Local identity
    # reaches up to ~67% (11% of pairs exceed CBAR_VMAX=40) - clipped so the KDE reflects
    # an honest pile-up at the boundary rather than silently truncating; global identity
    # already fits entirely within CBAR_VMAX, so its clip is a no-op.
    seqid_local_values = np.clip(matrix[np.tril_indices(seqid_n, k=-1)], CBAR_VMIN, CBAR_VMAX)
    seqid_global_values = np.clip(matrix[np.triu_indices(seqid_n, k=1)], CBAR_VMIN, CBAR_VMAX)
    b_dax = plot_comparison_heatmap(stylize(fig.add_subplot(right_gs[0, 0])), matrix, gene_labels, seqid_cmap,
                                     CBAR_VMIN, CBAR_VMAX, "Sequence identity (%)",
                                     ticks=np.arange(CBAR_VMIN, CBAR_VMAX + 1, 5), line_color=nc.crimson, abc="b",
                                     dist_values=seqid_local_values, overlay_values=seqid_global_values)

    # Different color family from panel b (cobalt vs crimson) so the two panels are
    # visually distinct at a glance; still reversed so dark reads as "more similar" here
    # too, even though for RMSD that means a LOW value. Gradient itself is compressed
    # into [0, STRUCT_GRADIENT_MAX] and held flat (white) beyond that, out to
    # STRUCT_VMAX - so the axis/colorbar reads to 20 while the color only distinguishes
    # differences within the tighter, more informative 0-10 Å range.
    struct_cmap_full = stylia.FadingColormap("cobalt", transformation=None).cmap.reversed()
    struct_cmap = clipped_gradient_cmap(struct_cmap_full, STRUCT_VMIN, STRUCT_VMAX, STRUCT_GRADIENT_MAX)
    struct_n = struct_matrix.shape[0]
    # Distribution panel: main filled curve = local RMSD (lower triangle), single shared
    # 0-20 scale with the heatmap/colorbar (no second axis). Only 2/206 (~1%) global
    # values exceed 20 (max ~22), so no clipping needed for the overlay KDE either.
    # Global RMSD (upper triangle) is overlaid as a thin dashed black line, same
    # convention as panel b.
    struct_local_values = struct_matrix[np.tril_indices(struct_n, k=-1)]
    struct_global_values = struct_matrix[np.triu_indices(struct_n, k=1)]
    plot_comparison_heatmap(stylize(fig.add_subplot(right_gs[1, 0])), struct_matrix, gene_labels, struct_cmap,
                             STRUCT_VMIN, STRUCT_VMAX, "Structural RMSD (Å)",
                             ticks=np.arange(STRUCT_VMIN, STRUCT_VMAX + 1, 5), line_color=nc.cobalt, abc="c",
                             dist_values=struct_local_values, overlay_values=struct_global_values)

    # Third distinct color family (lime). Reversed for the same reason as panel c — for
    # cosine distance, LOW value = more similar. Same tick set drives both the boxplot and
    # the colorbar (via the shared `ticks` param), so they can't drift out of sync.
    pocketvec_cmap = stylia.FadingColormap("lime", transformation=None).cmap.reversed()
    plot_comparison_heatmap(stylize(fig.add_subplot(right_gs[2, 0])), pocketvec_matrix, gene_labels, pocketvec_cmap,
                             POCKETVEC_VMIN, POCKETVEC_VMAX, "PocketVec cosine distance",
                             ticks=np.arange(POCKETVEC_VMIN, POCKETVEC_VMAX + 0.001, 0.02), line_color=nc.lime, abc="d")

    # "a" is placed via fig.text() at panel b's own measured title height (rather than
    # first_structure_ax.set_title(pad=...)) because panel b's "b" label sits on dax,
    # which starts BELOW the boxplot (tax) appended above it — the structures grid has no
    # such panel above its top row, so the two axes' own tops don't line up on their own.
    fig.canvas.draw()
    b_title_bbox = b_dax.title.get_window_extent()
    title_y_fig = fig.transFigure.inverted().transform((0, b_title_bbox.y0))[1]
    a_x_fig = first_structure_ax.get_position().x0 - 0.015
    fig.text(a_x_fig, title_y_fig, "a", fontweight="bold", fontsize=stylia.FONTSIZE_BIG,
             ha="left", va="bottom", transform=fig.transFigure)

    # Skips stylia.save_figure()'s wrapper (its hardcoded plt.tight_layout() call doesn't
    # handle nested GridSpecs predictably) — spacing here is fully controlled by the
    # subgridspecs' own hspace/wspace and each heatmap's make_axes_locatable pads.
    # pad_inches adds a small uniform margin around the tightly-cropped content (default
    # bbox_inches="tight" crop otherwise leaves ~zero space on any side, including left).
    for ext in ("png", "pdf"):
        output_path = os.path.join(plots_dir, f"figure_1.{ext}")
        plt.savefig(output_path, dpi=600, transparent=False, bbox_inches="tight", pad_inches=0.15)
        print(f"Saved to {output_path}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", default=False, help="Force PyMOL re-rendering for every protein, even if its figure_*.png already exists")
    parser.add_argument("--remove", action="store_true", default=False, help="Delete the per-protein figure_*.png renders after building the composite (default: keep them)")
    parser.add_argument("--store-pymol", action="store_true", default=False, help="Keep the per-protein PyMOL .pse sessions instead of deleting them after use")
    parser.add_argument("--show-grids", action="store_true", default=False, help="Draw each subplot's GridSpec rectangle (red spines) to visualize the underlying grid layout")
    args = parser.parse_args()

    render_protein_structures(proteins=PROTEINS, rerun=args.rerun)
    render_master_figure(show_grids=args.show_grids)
    cleanup_intermediate_files(store_pymol=args.store_pymol, remove_pngs=args.remove)
