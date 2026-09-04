"""
Figure 1: protein structures (panel a) + sequence identity / structural RMSD heatmaps
(panels b/c), from figure_1_calculations.py, + lysS's 5 pipeline structures superposed as a
cartoon-only PyMOL ensemble render (panel d, render_lysS_ensemble). Panel d's earlier content
(P2Rank pocket-probability curve + domain-annotation strips) has been split out into its own
standalone supplementary figure, scripts/plots/figure_supp_pocket.py.

Built on figure_4/figure_5's panel-per-file architecture instead of the original,
oversized monolithic composite this script replaced: each panel is its own standalone
figure, saved as its own PDF (Fig_1a.pdf ... Fig_1d.pdf) at an EXACT physical size read
from output/plots/figure_1/panel_layout.csv (columns: panel, x, delta_x, y, delta_y,
padding, in cm) - so stylia's fixed-point-size text renders at its true intended size on
the page, instead of shrinking along with an oversized canvas once placed at final print
width. Panel letters are baked onto each saved PDF (add_panel_label, from figure_4) and a
Nature two-column-width guideline check (MAX_WIDTH_IN, also from figure_4) warns if any
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
from PIL import Image
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


def autocrop_to_content_with_offset(img, padding_frac=0.05, background_frac=0.98):
    """autocrop_to_content, but also returns the (x_offset, y_offset) of the returned crop's
    top-left corner within `img`'s original pixel grid - needed to translate a pixel location
    computed against the original, uncropped render (e.g. a projected point) into the cropped
    image's own pixel coordinates."""
    gray = img[..., :3].mean(axis=2)
    max_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
    mask = gray < background_frac * max_val
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return img, 0, 0
    h, w = img.shape[:2]
    y0, y1 = rows[0], rows[-1]
    x0, x1 = cols[0], cols[-1]
    pad_y = int((y1 - y0) * padding_frac)
    pad_x = int((x1 - x0) * padding_frac)
    y0 = max(y0 - pad_y, 0)
    y1 = min(y1 + pad_y, h - 1)
    x0 = max(x0 - pad_x, 0)
    x1 = min(x1 + pad_x, w - 1)
    return img[y0:y1 + 1, x0:x1 + 1], x0, y0


def autocrop_to_content(img, padding_frac=0.05, background_frac=0.98):
    cropped, _, _ = autocrop_to_content_with_offset(img, padding_frac, background_frac)
    return cropped


def crop_for_display(img, zoom=1.0):
    """Same crop pipeline as show_zoomed_image (center-crop for zoom>1, then autocrop, then
    pad for zoom<1), but also returns the (x_offset, y_offset) of the returned crop's top-left
    corner within `img`'s original pixel grid. A pixel location computed against the original,
    uncropped image (e.g. a pocket centroid's projected screen position) converts into this
    crop's own data coordinates via (orig_x - x_offset, orig_y - y_offset) - valid because
    ax.imshow's default data coordinates are just the displayed array's own row/col indices."""
    h, w = img.shape[:2]
    x_off, y_off = 0, 0
    if zoom > 1:
        crop_w = int(w / zoom)
        crop_h = int(h / zoom)
        cx, cy = w // 2, h // 2
        x0 = max(cx - crop_w // 2, 0)
        x1 = min(cx + crop_w // 2, w)
        y0 = max(cy - crop_h // 2, 0)
        y1 = min(cy + crop_h // 2, h)
        img = img[y0:y1, x0:x1]
        x_off += x0
        y_off += y0

    img, cx0, cy0 = autocrop_to_content_with_offset(img)
    x_off += cx0
    y_off += cy0

    if zoom < 1:
        h2, w2 = img.shape[:2]
        pad_h = int(h2 / zoom) - h2
        pad_w = int(w2 / zoom) - w2
        pad_val = 1.0 if np.issubdtype(img.dtype, np.floating) else 255
        top, left = pad_h // 2, pad_w // 2
        pad_widths = [(top, pad_h - top), (left, pad_w - left)]
        if img.ndim == 3:
            pad_widths.append((0, 0))
        img = np.pad(img, pad_widths, mode="constant", constant_values=pad_val)
        x_off -= left
        y_off -= top

    return img, x_off, y_off


def show_zoomed_image(ax, img, zoom=1.0):
    cropped, _, _ = crop_for_display(img, zoom)
    ax.imshow(cropped, interpolation="none", resample=False)
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


# lysS's 5 pipeline structures (matches the "# structures" count in pocket_detection_data.csv
# for P9WFU9, not the full output/aligned_relaxed_structures/P9WFU9/ directory - that directory
# also holds alphafold3 models 2-4 and chai1 models 0/2-4, which were never carried into the
# pocket-detection pipeline). Already aligned to a shared coordinate frame (see
# figure_1_calculations.py's structural-RMSD section), so loading them as separate PyMOL objects
# in the same session superposes them with no extra cmd.align/super step needed.
LYSS_UNIPROT = "P9WFU9"
LYSS_STRUCTURES = [
    "alphafold2_P9WFU9_model_0.pdb",
    "alphafold3_P9WFU9_model_0.pdb",
    "alphafold3_P9WFU9_model_1.pdb",
    "chai1_P9WFU9_model_1.pdb",
    "swissmodel_P9WFU9_model_0.pdb",
]

# Fixed camera view for the lysS ensemble render, in the same cmd.get_view()-tuple format as
# panel a's MY_VIEWS - picked interactively by the user from the saved .pse session.
LYSS_VIEW = (
    0.219232261, -0.238682613, 0.946026027,
    0.975314498, 0.079742722, -0.205896914,
    -0.026288977, 0.967810452, 0.250276804,
    0.000000000, 0.000000000, -419.003662109,
    1.055000305, -3.290000916, 6.551498413,
    392.760864258, 445.246520996, 20.000000000,
)

# Opaque cartoon (0, not the transparency briefly tried before) - PyMOL's ray_trace_mode=1
# black silhouette/contour outlines are only drawn on fully opaque geometry, not transparent
# cartoon, and those outlines were requested (kept) over the transparency, per explicit choice.
CARTOON_TRANSPARENCY = 0.0

# Same same-pocket centroid-distance cutoff as figure_1_calculations.py's dedup pass (that
# script's own POCKET_DEDUP_DISTANCE_THRESHOLD) - duplicated here as a literal rather than
# imported, since figure_1_calculations.py is a top-level script, not an importable module.
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14

NB_SPHERES_SIZE = 4.0

# PyMOL's ray_trace_mode black outline is ~1-2 raw pixels wide regardless of ray_trace_gain
# (verified empirically - gain from 0.001 to 0.3 made no visible difference); the only lever
# that actually thins it is supersampling: ray-trace at RAY_SUPERSAMPLE x the final pixel
# size, then downsample with antialiasing (PIL LANCZOS) - the outline's fixed pixel width
# shrinks relative to the whole image once blended into surrounding white during that
# downsample, making it read as a thinner line.
RAY_RESOLUTION = 1200
RAY_SUPERSAMPLE = 2

# Grayscale RGB triples for the canonical pocket spheres, in order (per request: white, dark
# gray, black) - each still gets a black ray_trace_mode outline, so even the white/black
# spheres read as distinct shaded rings against the white page background.
POCKET_COLORS = [
    (1.00, 1.00, 1.00),  # white
    (0.35, 0.35, 0.35),  # dark gray
    (0.00, 0.00, 0.00),  # black
]

# Saturated, mutually-exclusive marker colors used ONLY for the label-position detection pass
# (never seen in the actual figure) - one per canonical pocket group, in POCKET_COLORS order.
POCKET_DETECT_COLORS = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]


def _lysS_canonical_pocket_groups():
    """Every pocket centroid detected for lysS across its 5 pipeline structures, each assigned
    to a canonical (deduplicated) pocket group. Canonical centroids come from the same greedy
    accept-if-farther-than-6.14A pass as figure_1_calculations.py's dedup (sorted by Pocket
    score descending); every other (rejected) centroid is then assigned to its NEAREST
    canonical centroid - valid because a rejected centroid is, by construction of that greedy
    pass, always within the threshold of at least one canonical centroid already accepted
    before it. Returns (centroids, group_indices, n_canonical)."""
    lys_pockets = pocket_detection_data[pocket_detection_data['Uniprot AC'] == LYSS_UNIPROT] \
        .sort_values('Pocket score', ascending=False)
    all_centroids = [
        np.array([float(v) for v in row['Pocket centroid coordinate (x y z)'].split()])
        for _, row in lys_pockets.iterrows()
    ]
    accepted_centroids = []
    for centroid in all_centroids:
        if all(np.linalg.norm(centroid - c) > POCKET_DEDUP_DISTANCE_THRESHOLD for c in accepted_centroids):
            accepted_centroids.append(centroid)
    group_indices = [
        int(np.argmin([np.linalg.norm(centroid - c) for c in accepted_centroids]))
        for centroid in all_centroids
    ]
    return all_centroids, group_indices, len(accepted_centroids)


def _flatten_to_opaque_white(png_path, target_resolution):
    """Composites a ray-traced PyMOL PNG (RGBA, background left transparent by cmd.ray/cmd.png
    regardless of cmd.bg_color - confirmed empirically) onto an OPAQUE white canvas, then
    downsamples to target_resolution x target_resolution with antialiasing (PIL LANCZOS - see
    RAY_SUPERSAMPLE's own comment for why), overwriting png_path. Compositing to opaque RGB
    BEFORE resizing is required, not optional: resizing the raw RGBA directly (tried first) left
    transparent regions with corrupted (0,0,0,0) - i.e. BLACK, not the intended white - pixels
    afterward (a Pillow resize-of-unpremultiplied-alpha artifact), which silently broke
    autocrop_to_content's brightness-based content detection downstream (it ignores alpha and
    only looks at RGB, so it saw near-black "content" almost everywhere, and effectively
    stopped cropping the large blank margin around the actual structure)."""
    with Image.open(png_path) as im:
        im = im.convert("RGBA")
        opaque_white = Image.new("RGBA", im.size, (255, 255, 255, 255))
        flattened = Image.alpha_composite(opaque_white, im).convert("RGB")
        flattened.resize((target_resolution, target_resolution), Image.LANCZOS).save(png_path)


def render_lysS_ensemble(rerun=False):
    """All 5 of lysS's pipeline structures superposed in one cartoon-only render (no surface),
    every structure colored in lysS's own color from color_mapping.json - shows their
    conformational spread now that they share one coordinate frame. All detected pocket
    centroids across those 5 structures are also shown as nb_spheres, colored by canonical
    (deduplicated) pocket group - see _lysS_canonical_pocket_groups. Also saves a .pse session
    (plots_dir/session_{LYSS_UNIPROT}_lysS_ensemble.pse) so the camera angle can be picked
    interactively in PyMOL, then hardcoded into LYSS_VIEW above via cmd.get_view().

    Returns (png_path, centroid_px): centroid_px maps each canonical pocket group index (int,
    same indexing as POCKET_COLORS) to that pocket's (x, y) pixel location in png_path's own
    RAW, uncropped pixel grid - found via a second, quick detection-only ray pass where the
    cartoon is hidden and each pocket group is temporarily recolored to a unique saturated
    marker (POCKET_DETECT_COLORS), never seen in the actual figure, then locating each marker
    color's pixel centroid. Needed so build_panel_c can draw a leader line from a "Pocket A/B/C"
    label to the right sphere, since PyMOL doesn't otherwise expose a screen-space projection
    API - crop_for_display's returned offsets convert this raw-grid location into whatever
    cropped/zoomed data-coordinate system that panel actually displays."""
    gene = "lysS"
    png_path = os.path.join(plots_dir, f"figure_{LYSS_UNIPROT}_{gene}_ensemble.png")
    centroid_px_path = os.path.join(plots_dir, f"figure_{LYSS_UNIPROT}_{gene}_ensemble_centroid_px.json")
    if os.path.exists(png_path) and os.path.exists(centroid_px_path) and not rerun:
        print(f"Reusing existing ensemble render for {gene} ({LYSS_UNIPROT}): {png_path}")
        with open(centroid_px_path) as f:
            centroid_px = {int(k): tuple(v) for k, v in json.load(f).items()}
        return png_path, centroid_px

    pymol.finish_launching(['pymol', '-cq'])
    cmd.reinitialize()

    color_rgb = mcolors.to_rgb(cmap_dict[gene])
    cmd.set_color("lysS_color", list(color_rgb))

    for i, fname in enumerate(LYSS_STRUCTURES):
        obj_name = f"structure_{i}"
        cmd.load(os.path.join(output_dir, "aligned_relaxed_structures", LYSS_UNIPROT, fname), obj_name)
        cmd.hide("everything", obj_name)
        cmd.show("cartoon", obj_name)
        cmd.color("lysS_color", obj_name)
        cmd.set("cartoon_transparency", CARTOON_TRANSPARENCY, obj_name)

    centroids, group_indices, n_canonical = _lysS_canonical_pocket_groups()
    if n_canonical > len(POCKET_COLORS):
        raise ValueError(f"lysS has {n_canonical} canonical pockets, more than the "
                          f"{len(POCKET_COLORS)} curated POCKET_COLORS - add more colors.")
    for i in range(n_canonical):
        cmd.set_color(f"lysS_pocket_{i}", list(POCKET_COLORS[i]))
    cmd.set("nb_spheres_size", NB_SPHERES_SIZE)
    for j, (centroid, group) in enumerate(zip(centroids, group_indices)):
        obj_name = f"pocket_centroid_{j}"
        cmd.pseudoatom(obj_name, pos=[float(v) for v in centroid])
        cmd.show("nb_spheres", obj_name)
        cmd.color(f"lysS_pocket_{group}", obj_name)
    print(f"  lysS: {len(centroids)} pocket centroids -> {n_canonical} canonical groups: {group_indices}")

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
    if LYSS_VIEW:
        cmd.set_view(LYSS_VIEW)
    zoom_to_fixed_box("all", box_size=95)
    session_path = os.path.join(plots_dir, f"session_{LYSS_UNIPROT}_{gene}_ensemble.pse")
    cmd.save(session_path)
    hires = RAY_RESOLUTION * RAY_SUPERSAMPLE

    cmd.ray(hires, hires)
    cmd.png(png_path, dpi=600)
    _flatten_to_opaque_white(png_path, RAY_RESOLUTION)

    # Detection-only pass (never saved as a visible figure asset): hide the cartoons and
    # recolor each pocket group to a unique saturated marker, so its pixel centroid can be
    # found unambiguously in the rendered raster - reuses the exact same camera/zoom already
    # locked in above, so positions line up with the real production render pixel-for-pixel.
    for i in range(len(LYSS_STRUCTURES)):
        cmd.hide("everything", f"structure_{i}")
    for i in range(n_canonical):
        cmd.set_color(f"lysS_pocket_detect_{i}", list(POCKET_DETECT_COLORS[i]))
    for j, group in enumerate(group_indices):
        cmd.color(f"lysS_pocket_detect_{group}", f"pocket_centroid_{j}")
    cmd.set("ray_trace_mode", 0)
    detect_path = os.path.join(plots_dir, f"_tmp_{LYSS_UNIPROT}_{gene}_detect.png")
    cmd.ray(hires, hires)
    cmd.png(detect_path, dpi=600)
    _flatten_to_opaque_white(detect_path, RAY_RESOLUTION)
    with Image.open(detect_path) as im:
        detect_arr = np.array(im.convert("RGB")).astype(float) / 255.0
    os.remove(detect_path)

    centroid_px = {}
    for i in range(n_canonical):
        dist = np.linalg.norm(detect_arr - np.array(POCKET_DETECT_COLORS[i]), axis=-1)
        ys, xs = np.where(dist < 0.3)
        if len(xs):
            centroid_px[i] = (float(xs.mean()), float(ys.mean()))
    with open(centroid_px_path, "w") as f:
        json.dump(centroid_px, f)

    cmd.delete("all")
    print(f"Saved PyMOL session: {session_path}")
    return png_path, centroid_px


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
# Panel-saving scaffolding (ported from figure_4_plot.py / figure_5_plot.py)
# ===========================================================================

PANEL_LABEL_MARGIN = 0.02


def add_panel_label(fig, letter):
    """Bold panel letter at the top-left of the FIGURE (page), fixed regardless of
    padding - from figure_4_plot.py."""
    fig.text(PANEL_LABEL_MARGIN, 1 - PANEL_LABEL_MARGIN, letter, fontweight="bold",
              fontsize=stylia.FONTSIZE_BIG, color=get_fg_color(), ha="left", va="top",
              transform=fig.transFigure)


MAX_WIDTH_IN = stylia.SIZE  # Nature two-column guideline (~7.09in/18cm) - from figure_4_plot.py


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
    tight_layout() warns/fails (e.g. panel a's packed grid) - from figure_5_plot.py.
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
    each panel already saved at its exact target size. Ported from figure_5_merge.py, kept
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


# Extra crop-in factor on top of autocrop_to_content, per request ("too much blank space") -
# tune by eye if the structure still isn't filling the panel, or clips at this value.
LYSS_PANEL_ZOOM = 1.3

# Canonical pocket group index (POCKET_COLORS order: white, dark gray, black) -> (label text,
# (x, y) label position in AXES FRACTION coords, horizontal/vertical text alignment). Per
# request: "Pocket A" (top-right) must point at the black sphere (group 2, not the original
# white/group 0), "Pocket B" sits mid-left (not bottom-left), and "Pocket C" took over group 0's
# (white) former bottom-right slot. A black leader line (ax.annotate's arrowprops) connects
# each corner/edge label to its actual sphere position regardless of these swaps.
LYSS_POCKET_LABELS = {
    2: ("Pocket A", (0.95, 0.95), "right", "top"),
    1: ("Pocket B", (0.05, 0.5), "left", "center"),
    0: ("Pocket C", (0.95, 0.05), "right", "bottom"),
}


def build_panel_c(size, padding):
    """Structural RMSD heatmap - moved here (from panel d) per request, swapping places with
    the lysS ensemble render (now panel d, build_panel_d)."""
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


def build_panel_d(size, padding, rerun=False):
    """lysS's 5 pipeline structures, superposed as cartoon-only renders (no surface) all in
    lysS's own color - see render_lysS_ensemble - cropped in further (LYSS_PANEL_ZOOM) to
    remove blank margin. Each canonical pocket's sphere is labeled "Pocket A/B/C (n/total)" -
    n being how many of lysS's raw pocket detections (across its 5 structures) that canonical
    pocket absorbed during dedup, out of the total pocket structures for lysS - at a fixed
    corner/edge with a leader line to its actual (detected) position - see
    render_lysS_ensemble's centroid_px return and LYSS_POCKET_LABELS. Moved here (from panel c)
    per request, swapping places with the structural RMSD heatmap (now panel c,
    build_panel_c)."""
    png_path, centroid_px = render_lysS_ensemble(rerun=rerun)
    _, group_indices, _ = _lysS_canonical_pocket_groups()
    group_counts = {i: group_indices.count(i) for i in set(group_indices)}

    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor("white")
    stylize(ax)
    ax.set_axis_off()

    img = mpimg.imread(png_path)
    cropped, x_off, y_off = crop_for_display(img, zoom=LYSS_PANEL_ZOOM)
    ax.imshow(cropped, interpolation="none", resample=False)
    ax.set_aspect("equal", adjustable="datalim")

    for i, (label_text, label_pos, ha, va) in LYSS_POCKET_LABELS.items():
        if i not in centroid_px:
            continue
        full_label = f"{label_text} ({group_counts.get(i, 0)})"
        raw_x, raw_y = centroid_px[i]
        data_xy = (raw_x - x_off, raw_y - y_off)
        # Label background matches that pocket's own sphere color (POCKET_COLORS), with text
        # color picked for contrast (readable_text_color) - white/dark gray/black spheres need
        # black/white/white label text respectively.
        facecolor = POCKET_COLORS[i]
        ax.annotate(full_label, xy=data_xy, xycoords="data", xytext=label_pos,
                    textcoords="axes fraction", ha=ha, va=va, fontsize=stylia.FONTSIZE,
                    color=readable_text_color(facecolor),
                    bbox=dict(boxstyle="round,pad=0.3", facecolor=facecolor, edgecolor="black",
                              linewidth=stylia.LINEWIDTH),
                    arrowprops=dict(arrowstyle="-", color="black", linewidth=stylia.LINEWIDTH))

    stylia.label(ax, xlabel="", ylabel="")
    save_panel(fig, "d", use_tight_layout=False, padding=padding)


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
        build_panel_d(sizes["d"], paddings["d"], rerun=rerun)

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
