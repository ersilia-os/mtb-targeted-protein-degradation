"""
Merges figure_3's 6 standalone panel PDFs (Fig_3a.pdf ... Fig_3f.pdf, from figure_3_plot.py) into
one master PDF, each panel placed at its own (x, y) from output/plots/figure_3/panel_layout.csv.

x/y are top-left origin, y growing downward (page-layout convention) - converted below to PDF's
own bottom-left/y-up coordinate space for placement. Each panel PDF already saves at EXACTLY its
delta_x/delta_y size (see figure_3_plot.py's save_panel()), so placement is a pure translation, no
scaling - content stays true vector (no rasterization) via pypdf's page-merging.

Master page size is the tightest bounding box containing every placed panel (max over panels of
x + delta_x, and of y + delta_y) - no separate size to configure.

Usage:
    python figure_3_merge.py
"""
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "..", "src"))

import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation

plots_dir = os.path.join(root, "..", "..", "output", "plots", "figure_3")
panel_layout_path = os.path.join(plots_dir, "panel_layout.csv")

PANEL_LETTERS = ["a", "b", "c", "d", "e", "f"]
CM_TO_PT = 72 / 2.54


def main():
    if not os.path.exists(panel_layout_path):
        raise FileNotFoundError(f"{panel_layout_path} not found.")
    df = pd.read_csv(panel_layout_path).set_index("panel")
    missing = [p for p in PANEL_LETTERS if p not in df.index]
    if missing:
        raise ValueError(f"{panel_layout_path} is missing row(s) for panel(s): {missing}")

    total_width_cm = max(df.loc[p, "x"] + df.loc[p, "delta_x"] for p in PANEL_LETTERS)
    total_height_cm = max(df.loc[p, "y"] + df.loc[p, "delta_y"] for p in PANEL_LETTERS)

    writer = PdfWriter()
    master_page = writer.add_blank_page(width=total_width_cm * CM_TO_PT, height=total_height_cm * CM_TO_PT)

    for p in PANEL_LETTERS:
        panel_path = os.path.join(plots_dir, f"Fig_3{p}.pdf")
        panel_page = PdfReader(panel_path).pages[0]

        x_top_cm, y_top_cm, delta_y_cm = df.loc[p, "x"], df.loc[p, "y"], df.loc[p, "delta_y"]
        x_pt = x_top_cm * CM_TO_PT
        # top-left/y-down (x, y_top) -> this panel's bottom-left corner in the master's native
        # bottom-left/y-up PDF space.
        y_bottom_pt = (total_height_cm - y_top_cm - delta_y_cm) * CM_TO_PT

        master_page.merge_transformed_page(panel_page, Transformation().translate(tx=x_pt, ty=y_bottom_pt))
        print(f"Placed Fig_3{p}.pdf at x={x_top_cm}cm, y={y_top_cm}cm (top-left)")

    output_path = os.path.join(plots_dir, "Fig_3_full.pdf")
    with open(output_path, "wb") as f:
        writer.write(f)
    print(f"Saved merged master figure ({total_width_cm:.2f} x {total_height_cm:.2f} cm) to {output_path}")


if __name__ == "__main__":
    main()
