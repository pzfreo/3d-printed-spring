"""Build a 3D WSO2 logo (brand symbol + wordmark) from the official SVG.

Exports a single combined solid plus a two-colour split (raised artwork and
backing plate as separate, co-located bodies) for multi-material printing.

Requires build123d:  pip install build123d
Usage:               python wso2_logo_3d.py [path/to/WSO2-Logo-Black.svg]
Outputs (next to this script):
    wso2_logo_3d.step                 combined single solid
    wso2_artwork.step / .stl          raised glyphs + symbol  (colour 1)
    wso2_plate.step   / .stl          rounded backing plate   (colour 2)
"""

import os
import sys
from build123d import (
    import_svg, Compound, scale, Pos, Box, Axis,
    extrude, fillet, export_step, export_stl,
)

TARGET_W = 120.0   # mm, final logo width
RELIEF   = 3.0     # mm, raised height of the artwork
PLATE_T  = 2.0     # mm, backing-plate thickness
MARGIN   = 6.0     # mm, border around the artwork
FILLET_R = 5.0     # mm, corner radius of the plate

HERE = os.path.dirname(os.path.abspath(__file__))


def build(svg_path):
    raw = Compound(list(import_svg(svg_path)))
    s = TARGET_W / raw.bounding_box().size.X

    # Scale to size and recentre on the origin
    art = scale(raw, by=s)
    c = art.bounding_box().center()
    art = Pos(-c.X, -c.Y, 0) * art

    # Raised artwork (colour 1): each closed region extruded up from z=0
    artwork = Compound([extrude(f, amount=RELIEF) for f in art.faces()])

    # Backing plate (colour 2): rounded box from z=-PLATE_T to z=0
    bb = artwork.bounding_box()
    plate = Pos(bb.center().X, bb.center().Y, -PLATE_T / 2) * Box(
        bb.size.X + 2 * MARGIN, bb.size.Y + 2 * MARGIN, PLATE_T
    )
    plate = fillet(plate.edges().filter_by(Axis.Z), radius=FILLET_R)

    return plate, artwork, plate + artwork


if __name__ == "__main__":
    svg = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "WSO2-Logo-Black.svg")
    plate, artwork, combined = build(svg)

    export_step(combined, os.path.join(HERE, "wso2_logo_3d.step"))
    export_step(artwork, os.path.join(HERE, "wso2_artwork.step"))
    export_stl(artwork, os.path.join(HERE, "wso2_artwork.stl"))
    export_step(plate, os.path.join(HERE, "wso2_plate.step"))
    export_stl(plate, os.path.join(HERE, "wso2_plate.stl"))

    print(f"combined: {combined.volume:.0f} mm^3, "
          f"size {[round(v, 1) for v in combined.bounding_box().size]} mm")
