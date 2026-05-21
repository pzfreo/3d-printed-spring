# 3d-printed-spring

An **orthoplanar spring with an isotropic resonant frequency**, designed in [build123d](https://github.com/gumyr/build123d) and validated with [Netgen / NGSolve](https://ngsolve.org/) finite-element modal analysis.

![Isometric render of the spring](docs/spring_iso.png)

## In plain English

Picture a flat disc the size of a £2 coin: a thick outer ring, a heavy little post in the middle, and three thin wiggly arms in the gap connecting them. The outer ring is what you'd clamp down, and the central post is a small proof mass. Flick the centre in any direction and it bounces back and forth on the arms like a mass on a spring. Different springs prefer to bounce in different directions — a typical flat spring is happy to bobble up and down but feels rigid edge-on. The goal here was a spring that feels **equally bouncy in all three directions** — up/down, left/right, front/back — so a tap produces essentially the same little hum whichever way you push the proof mass. That property is what "isotropic resonant frequency" means.

Making it work meant fighting against the natural anisotropy of a flat spring. Flat springs are floppy out of plane and stiff in plane; we wanted them equal. Two design moves get us there: the arms are made **taller than they are wide** (1.0 mm in plane × 2.7 mm thick), which softens sideways flexing relative to up-and-down, and each arm follows a **wiggly, snake-like path** so that pushing the centre sideways mostly bends the arms gently rather than stretching them lengthwise (stretching is what was making sideways feel too stiff). After tuning, an open-source physics simulator (Netgen + NGSolve, the same kind of finite-element analysis engineers use to predict how parts vibrate) confirmed all three bouncing directions land at roughly **33 Hz, within 2 % of each other**.

## Design parameters

| Parameter | Value | Notes |
|---|---|---|
| Outer diameter | 60 mm | clamp ring 5 mm wide (rim ID 50 mm) |
| Proof-mass hub diameter | 18 mm | with a 12 mm cylinder hanging below as mass |
| Plate / arm thickness (z) | 2.7 mm | |
| Arm width (in plane) | 1.0 mm | square-ish cross-section, slightly tall |
| Number of arms | 3 | 120° rotational symmetry → enforces f_x ≡ f_y |
| Arm sweep | 220° | with a sinusoidal 4-lobe radial wiggle, 4 mm amplitude |
| PLA mass | ~14.7 g | density 1.24 g/cm³ |
| First three resonances | **32.6 / 33.2 / 33.2 Hz** (z / x / y) | ratio 1.017 (1.7 %) |

See [`DESIGN.md`](DESIGN.md) for the engineering writeup — why orthoplanar springs are naturally anisotropic, the two knobs that fix it, the FEA workflow, and the parameter-sweep that landed on these numbers.

## Print assumptions

- **Material:** PLA, FDM (E ≈ 2.5 GPa assumed for FEA — anisotropic in real prints, so the *absolute* frequency will shift ±15 % from layer effects; the *isotropy ratio* is geometric and unaffected).
- **Nozzle / layer:** 0.4 mm nozzle, 0.2 mm layer height (3 perimeters cover the 1 mm arm width).
- **Orientation:** Flat — proof-mass post and clamp boss point *down* into the bed; print upside-down with supports under the rim/hub. This keeps the arms loaded in plane (perpendicular to layers) for the dominant axial mode.
- **Supports:** Only under the rim and hub bosses; the arms themselves are flat on the bed.

## Development

```
git clone https://github.com/pzfreo/3d-printed-spring.git
cd 3d-printed-spring
uv sync --group dev
```

The MCP server is preconfigured in `.mcp.json` — any Claude Code (or MCP-compatible client) opened in this repo launches the published [build123d-mcp](https://pypi.org/project/build123d-mcp/) server on first tool use.

Re-run the modal analysis (needs Python 3.12 because of the `cadquery-ocp` ↔ `netgen-occt` OCCT-version match described in build123d [issue #297](https://github.com/gumyr/build123d/issues/297)):

```
uv run --python 3.12 --with netgen-mesher --with ngsolve python -m spring.validate
```

Or sweep parameters:

```
uv run --python 3.12 --with netgen-mesher --with ngsolve python scripts/sweep_isotropy.py
```

## License

Apache-2.0.
