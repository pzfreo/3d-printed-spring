# WSO2 logo (3D)

A 3D model of the WSO2 logo — the circular brand symbol plus the `wso2`
wordmark — generated from the official SVG with [build123d](https://github.com/gumyr/build123d).

The artwork is raised as a relief on a rounded backing plate, so all glyphs
stay as one connected, printable piece. Overall size **132 × 42 × 5 mm**.

## Files

| File | What it is |
| --- | --- |
| `WSO2-Logo-Black.svg` | Source artwork (official WSO2 black logo) |
| `wso2_logo_3d.py` | Build script |
| `wso2_logo_3d.step` | Combined single solid |
| `wso2_artwork.step` / `.stl` | Raised glyphs + symbol — **colour 1** |
| `wso2_plate.step` / `.stl` | Backing plate — **colour 2** |

## Two-colour printing

`wso2_artwork` and `wso2_plate` share the same coordinate origin, so importing
both into a slicer drops them into place automatically — assign one colour to
the artwork and another to the plate. The artwork sits on the `z = 0` plane on
top of the 2 mm plate.

## Rebuild / tweak

```bash
pip install build123d
python wso2_logo_3d.py
```

Adjust the constants at the top of `wso2_logo_3d.py`:

- `TARGET_W` — overall width (mm)
- `RELIEF` — raised height of the artwork (mm)
- `PLATE_T` — backing-plate thickness (mm)
- `MARGIN` — border around the artwork (mm)
- `FILLET_R` — plate corner radius (mm)

## Note

The WSO2 name and logo are trademarks of WSO2 LLC. This model is for personal
/ non-commercial use.
