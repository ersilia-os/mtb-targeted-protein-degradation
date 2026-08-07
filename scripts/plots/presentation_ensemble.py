#!/usr/bin/env python3
"""
PyMOL session + per-structure PNGs for all 12 aligned aspS (P9WFW3) structures, plus a circular
summary figure built from those PNGs with matplotlib (this is a structure-image grid, not a data
plot, so stylia doesn't apply here - see feedback_stylia_not_for_molecule_grids memory).

Each structure is loaded as a cartoon (no surface), merged with the pockets P2Rank detected AND
that survived script 08's quality gate (probability > 0.2 / Top-3, then per-residue pLDDT/QSQE
filter) - i.e. every pocket_*.pdb under output/detected_pockets/P9WFW3/<structure>/pockets/. Many
of these raw per-structure pocket instances are the same physical pocket detected independently
across models, so they're grouped by identity using the same 6.14A centroid-distance dedup as
scripts/plots/figure_1_calculations.py (score-descending greedy accept), and each distinct pocket
gets its own color from stylia.SpectralColormap("npg") - same palette convention figure_1 uses
for its gene colors. Structure stays the neutral gray from scripts/69_pymol_visualizations.py.

All outputs saved to the repo's tmp/ (ad hoc, not a pipeline output).

Usage:
    python presentation_ensemble.py
    python presentation_ensemble.py --skip-render     # reuse existing PNGs, just rebuild the circular figure
    python presentation_ensemble.py --include-name    # write each structure's name next to it in the circular figure
"""
import argparse
import glob
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stylia
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")

GENE = "aspS"
UNIPROT_AC = "P9WFW3"
ALIGNED_DIR = os.path.join(ROOT, "output", "aligned_relaxed_structures", UNIPROT_AC)
DETECTED_POCKETS_DIR = os.path.join(ROOT, "output", "detected_pockets", UNIPROT_AC)
COLOR_MAPPING_PATH = os.path.join(ROOT, "output", "plots", "figure_1", "color_mapping.json")
POCKET_DATA_CSV = os.path.join(ROOT, "output", "pocket_detection_data.csv")
OUT_PATH = os.path.join(ROOT, "tmp", f"{GENE}_{UNIPROT_AC}_aligned_structures.pse")

COLOR_STRUCTURE = [0.7804, 0.8275, 0.8667]
POCKET_RADIUS = 5.0  # Angstrom, set directly via vdw (not sphere_scale, which is a multiplier)
POCKET_SPHERE_TRANSPARENCY = 0.0
SURFACE_TRANSPARENCY = 0.6  # applied to the final merged object (cmd.create doesn't inherit it)
POCKET_DEDUP_DISTANCE_THRESHOLD = 6.14  # matches scripts/plots/figure_1_calculations.py

PNG_DIR = os.path.join(ROOT, "tmp", f"{GENE}_{UNIPROT_AC}_structure_pngs")
PNG_WIDTH_PX = 2400
PNG_HEIGHT_PX = 2400
PNG_DPI = 300
CENTER_PNG_PATH = os.path.join(PNG_DIR, "_all_structures_cartoon.png")

CIRCLE_FIGURE_PATH = os.path.join(ROOT, "tmp", f"{GENE}_{UNIPROT_AC}_circular_summary.png")
CIRCLE_RADIUS = 0.9
CIRCLE_IMAGE_ZOOM = 0.13
CENTER_IMAGE_ZOOM = 0.27


def hex_to_rgb_floats(hex_color):
    hex_color = hex_color.lstrip("#")
    return [int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def compute_pocket_groups():
    """Returns:
    - pocket_group: (structure_name, pocket_number) -> distinct-pocket group id, via the same
      greedy dedup as figure_1_calculations.py (score-descending, 6.14A centroid threshold).
    - n_groups: number of distinct pockets.
    - representative_keys: {(structure_name, pocket_number)} of the highest-probability instance
      within each group - the one shown as a sphere in the combined all-structures overview.
    - group_color_index: group id -> palette index, ordered by group size descending, so color
      assignment tracks how often a pocket recurs across structures rather than its dedup order."""
    df = pd.read_csv(POCKET_DATA_CSV)
    df = df[df["Uniprot AC"] == UNIPROT_AC].sort_values("Pocket score", ascending=False).reset_index(drop=True)

    accepted_centroids = []
    pocket_group = {}
    row_groups = []
    for _, row in df.iterrows():
        centroid = np.array([float(v) for v in row["Pocket centroid coordinate (x y z)"].split()])
        dists = [np.linalg.norm(centroid - c) for c in accepted_centroids]
        if not dists or min(dists) > POCKET_DEDUP_DISTANCE_THRESHOLD:
            accepted_centroids.append(centroid)
            group = len(accepted_centroids) - 1
        else:
            group = int(np.argmin(dists))
        structure_name = row["File name"].replace(".pdb", "")
        pocket_group[(structure_name, int(row["Pocket number"]))] = group
        row_groups.append(group)
    df["group"] = row_groups
    n_groups = len(accepted_centroids)

    representative_keys = set()
    for _, group_df in df.groupby("group"):
        best = group_df.loc[group_df["Pocket probability"].idxmax()]
        representative_keys.add((best["File name"].replace(".pdb", ""), int(best["Pocket number"])))

    group_sizes = df.groupby("group").size().to_dict()
    order_by_size = sorted(group_sizes, key=lambda g: -group_sizes[g])
    group_color_index = {g: i for i, g in enumerate(order_by_size)}

    return pocket_group, n_groups, representative_keys, group_color_index


def render_structures(structure_paths):
    import pymol
    from pymol import cmd

    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()

    pocket_group, n_groups, representative_keys, group_color_index = compute_pocket_groups()
    palette = stylia.SpectralColormap("npg").sample(n_groups)
    for i, rgb in enumerate(palette):
        cmd.set_color(f"pocket_group_{i}", list(rgb[:3]))

    cmd.set_color("structure_color", COLOR_STRUCTURE)

    obj_names = []
    for structure_path in structure_paths:
        structure_name = os.path.splitext(os.path.basename(structure_path))[0]

        cmd.load(structure_path, "structure")
        cmd.color("structure_color", "structure")
        cmd.show("cartoon", "structure")
        cmd.show("surface", "structure")
        cmd.hide("lines", "structure")
        cmd.hide("sticks", "structure")

        pocket_paths = sorted(glob.glob(os.path.join(DETECTED_POCKETS_DIR, structure_name, "pockets", "pocket_*.pdb")))
        for pocket_path in pocket_paths:
            pocket_obj = os.path.splitext(os.path.basename(pocket_path))[0]
            pocket_number = int(pocket_obj.rsplit("_", 1)[-1])
            group = pocket_group[(structure_name, pocket_number)]
            color_index = group_color_index[group]
            is_representative = (structure_name, pocket_number) in representative_keys
            cmd.load(pocket_path, pocket_obj)
            cmd.color(f"pocket_group_{color_index}", pocket_obj)
            cmd.alter(pocket_obj, f"vdw={POCKET_RADIUS}")
            cmd.alter(pocket_obj, f"b={1.0 if is_representative else 0.0}")
            cmd.rebuild(pocket_obj)
            cmd.show("spheres", pocket_obj)
            cmd.set("sphere_transparency", POCKET_SPHERE_TRANSPARENCY, pocket_obj)

        merge_selection = " or ".join(["structure"] + [os.path.splitext(os.path.basename(p))[0] for p in pocket_paths])
        cmd.create("_tmp_merged", merge_selection)
        cmd.delete("structure")
        for pocket_path in pocket_paths:
            cmd.delete(os.path.splitext(os.path.basename(pocket_path))[0])
        cmd.set_name("_tmp_merged", structure_name)
        cmd.set("transparency", SURFACE_TRANSPARENCY, structure_name)

        obj_names.append(structure_name)
        print(f"  Loaded: {structure_name} ({len(pocket_paths)} surviving pocket(s))")

    cmd.order(" ".join(obj_names), location="top", sort=0)

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

    cmd.set_view((
        0.310638905, 0.337588221, 0.888560593,
        -0.185486570, 0.938364446, -0.291664898,
        -0.932254851, -0.074213326, 0.354110122,
        0.000000000, 0.000000000, -532.742675781,
        1.185535431, 1.602592468, 2.187229156,
        471.395843506, 594.089599609, 20.000000000,
    ))

    os.makedirs(PNG_DIR, exist_ok=True)
    for obj_name in obj_names:
        for other in obj_names:
            cmd.disable(other)
        cmd.enable(obj_name)
        png_path = os.path.join(PNG_DIR, f"{obj_name}.png")
        cmd.png(png_path, width=PNG_WIDTH_PX, height=PNG_HEIGHT_PX, dpi=PNG_DPI, ray=1)
        print(f"  Saved PNG: {png_path}")

    # Combined overview (all 12 aligned structures + pockets at once, cartoon only - no
    # surface, since with 12 overlapping surfaces the shells would just occlude each other).
    # Only the highest-probability instance per distinct pocket is shown (b=1.0, tagged above)
    # to avoid up to 5 overlapping spheres (e.g. the chai1 cluster) cluttering the same site.
    for obj_name in obj_names:
        cmd.enable(obj_name)
        cmd.hide("surface", obj_name)
    cmd.hide("spheres", "hetatm and b < 0.5")
    cmd.png(CENTER_PNG_PATH, width=PNG_WIDTH_PX, height=PNG_HEIGHT_PX, dpi=PNG_DPI, ray=1)
    print(f"  Saved PNG: {CENTER_PNG_PATH}")
    cmd.show("spheres", "hetatm and b < 0.5")

    for obj_name in obj_names:
        cmd.disable(obj_name)
    cmd.enable(obj_names[0])

    cmd.save(OUT_PATH)
    print(f"Saved: {OUT_PATH}")


def build_circular_figure(obj_names, gene_color_hex, include_name=False):
    fig, ax = plt.subplots(figsize=(16, 16))

    n = len(obj_names)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    for obj_name, angle in zip(obj_names, angles):
        x, y = CIRCLE_RADIUS * np.cos(angle), CIRCLE_RADIUS * np.sin(angle)

        png_path = os.path.join(PNG_DIR, f"{obj_name}.png")
        img = plt.imread(png_path)
        imagebox = OffsetImage(img, zoom=CIRCLE_IMAGE_ZOOM)
        ax.add_artist(AnnotationBbox(imagebox, (x, y), frameon=False, zorder=2))

        if include_name:
            ax.text(x, y - 0.18, obj_name, ha="center", va="top", fontsize=11, zorder=2)

    center_img = plt.imread(CENTER_PNG_PATH)
    center_imagebox = OffsetImage(center_img, zoom=CENTER_IMAGE_ZOOM)
    ax.add_artist(AnnotationBbox(center_imagebox, (0, 0), frameon=False, zorder=2))

    # Center marker + label (gene_to_color dot + "aspS" text, matching figure_1/2's gene
    # legend recipe) left empty for now - re-add using gene_color_hex when needed.

    margin = CIRCLE_RADIUS * 1.6
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(CIRCLE_FIGURE_PATH, dpi=PNG_DPI, bbox_inches="tight")
    print(f"Saved: {CIRCLE_FIGURE_PATH}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skip-render", action="store_true",
                         help="Skip PyMOL rendering and reuse the PNGs already in tmp/; only rebuild the circular summary figure.")
    parser.add_argument("--include-name", action="store_true",
                         help="Write each structure's name next to it in the circular summary figure.")
    args = parser.parse_args()

    structure_paths = sorted(glob.glob(os.path.join(ALIGNED_DIR, "*.pdb")))
    if not structure_paths:
        raise SystemExit(f"No aligned structures found in {ALIGNED_DIR}")
    obj_names = [os.path.splitext(os.path.basename(p))[0] for p in structure_paths]

    with open(COLOR_MAPPING_PATH) as f:
        gene_color_hex = json.load(f)["gene_to_color"][GENE]

    if not args.skip_render:
        render_structures(structure_paths)

    build_circular_figure(obj_names, gene_color_hex, include_name=args.include_name)


if __name__ == "__main__":
    main()
