"""
Merges figure_2's two panels into a single vector figure, both panels staying true vector
(no rasterization at any point):
  a) plots/figure_2/scheme.svg      - hand-authored Inkscape schematic (180mm x 215mm)
  b) plots/figure_2/figure_2.svg    - physicochemical property distributions (figure_2_plot.py)

Stacked vertically (a on top, b below), panel widths matched, lowercase "a"/"b" labels added
in stylia's own panel-label style (bold, Arial, print-format font size - see
stylia/figure/figure.py's label()).

svgutils.compose.SVG() strips each source file's own <svg> wrapper and inlines its raw content,
so embedding is always 1 (child viewBox unit) = 1 (parent user unit) - no implicit rescaling.
The parent Figure below is declared directly in mm, so its user-unit space is literally
millimeters; scheme.svg's own viewBox is already mm-equivalent (its "180mm"/"215mm" declared
size matches its viewBox numbers 1:1), so it needs no scaling. figure_2.svg is matplotlib
output (fix_mpl=True treats its declared pt-labelled width as bare viewBox units), so it's
explicitly scaled to match the parent's mm-based width.

Usage:
    python figure_2_merge.py
"""
import os

from lxml import etree
from svgutils.compose import Element, Figure, SVG, Text, Unit

SVG_NS = "http://www.w3.org/2000/svg"


def white_background(width_mm, height_mm):
    """Opaque white backing rect - avoids relying on viewer/consumer background color."""
    rect = etree.Element(f"{{{SVG_NS}}}rect", {
        "x": "0", "y": "0", "width": str(width_mm), "height": str(height_mm), "fill": "white",
    })
    return Element(rect)

root = os.path.dirname(os.path.abspath(__file__))
plots_dir = os.path.join(root, "..", "..", "plots", "figure_2")

# Matches stylia's panel-label style for print format/article style: bold, Arial, black
# (stylia/figure/figure.py: label(abc=...) -> fontweight="bold", fontsize=get_fontsize_big()=8pt).
LABEL_SIZE_MM = 8 * 25.4 / 72  # 8pt converted to mm, since this canvas's user units are mm
LABEL_MARGIN_MM = 5  # from each panel's top-left corner
LABEL_KWARGS = dict(size=LABEL_SIZE_MM, weight="bold", font="Arial", color="black")

PANEL_GAP_MM = 6  # breathing room between panel a and b, so the "b" label clears the plots
DISTRIBUTIONS_SIZE_FACTOR = 0.85  # render panel b a bit smaller than panel a (full scheme width)


def main():
    scheme = SVG(os.path.join(plots_dir, "scheme.svg"))
    distributions = SVG(os.path.join(plots_dir, "figure_2.svg"), fix_mpl=True)

    # scheme's own viewBox is already mm-equivalent (declared "180mm"/"215mm" == viewBox 180/215),
    # so this round-trips to the same numbers - scheme itself needs no scaling.
    scheme_width_mm = Unit(scheme.width).to("mm").value
    scheme_height_mm = Unit(scheme.height).to("mm").value

    # distributions.width/.height are figure_2.svg's raw (post-fix_mpl) viewBox units, not mm -
    # scale so its width is a fraction of scheme's width in the shared mm-based parent canvas.
    scale = (scheme_width_mm * DISTRIBUTIONS_SIZE_FACTOR) / distributions.width
    distributions.scale(scale)
    distributions_width_mm = distributions.width * scale
    distributions_height_mm = distributions.height * scale

    panel_b_x = (scheme_width_mm - distributions_width_mm) / 2  # center the smaller panel
    panel_b_y = scheme_height_mm + PANEL_GAP_MM
    distributions.move(panel_b_x, panel_b_y)

    label_a = Text("a", LABEL_MARGIN_MM, LABEL_MARGIN_MM + LABEL_SIZE_MM, **LABEL_KWARGS)
    # Sits within the gap, well above panel_b_y where the actual plot content starts.
    label_b = Text("b", LABEL_MARGIN_MM, scheme_height_mm + LABEL_SIZE_MM, **LABEL_KWARGS)

    total_width_mm = scheme_width_mm
    total_height_mm = panel_b_y + distributions_height_mm

    fig = Figure(
        f"{total_width_mm}mm",
        f"{total_height_mm}mm",
        white_background(total_width_mm, total_height_mm),
        scheme,
        distributions,
        label_a,
        label_b,
    )
    output_path = os.path.join(plots_dir, "figure_2_merged.svg")
    fig.save(output_path)
    print(f"Panel a (scheme): {scheme_width_mm:.1f} x {scheme_height_mm:.1f} mm")
    print(f"Panel b (distributions): {scheme_width_mm:.1f} x {distributions_height_mm:.1f} mm (scaled by {scale:.4f})")
    print(f"Total: {total_width_mm:.1f} x {total_height_mm:.1f} mm")
    print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()
