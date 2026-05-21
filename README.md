# 3d-printed-spring

Parametric 3D-printable spring designed in [build123d](https://github.com/gumyr/build123d).

## What this is for

> TBD — describe the spring's intended use (constraint range, application, what it replaces).

## Design parameters

> Fill these in as the design firms up.

| Parameter | Target | Notes |
|---|---|---|
| Spring rate (N/mm) | TBD | |
| Free length (mm) | TBD | |
| Outer diameter (mm) | TBD | |
| Wire / cross-section | TBD | round, square, ribbon, etc. |
| Number of active coils | TBD | |
| Pitch (mm) | TBD | |
| Direction | TBD | LH / RH |
| End condition | TBD | closed, open, ground, etc. |

## Print assumptions

> Pin these down before geometry — they constrain wire-profile minimum, layer adhesion margin, support strategy.

- **Printer:** TBD
- **Material:** TBD (TPU? PETG? PLA-CF? Resin?)
- **Nozzle / layer:** TBD
- **Orientation:** TBD (Z-axis vs flat)
- **Supports:** TBD

## Development

```
git clone https://github.com/pzfreo/3d-printed-spring.git
cd 3d-printed-spring
uv sync --group dev
```

The MCP server is preconfigured in `.mcp.json` — any Claude Code (or MCP-compatible client) opened in this repo launches the published [build123d-mcp](https://pypi.org/project/build123d-mcp/) server on first tool use.

## License

Apache-2.0.
