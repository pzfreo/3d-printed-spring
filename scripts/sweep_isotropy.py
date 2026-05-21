"""Parameter sweep to tune the orthoplanar spring to k_x = k_y = k_z.

Each row varies (arm_w, plate_t, arm_sweep_deg) and reports the first three
translational-mode frequencies plus the isotropy ratio f_xy / f_z.
"""
from __future__ import annotations
import sys
from spring.spring import SpringParams
from spring.validate import run, Material


def first_three_translational(results):
    """Pick first eigenmode dominated by each translation axis."""
    picked = {}
    for k, f, axis, vu, hub in results:
        # Require the hub to actually move — skip plate flexure modes
        # (hub displacement small) and rotational modes (hub displacement
        # exists but is much smaller than the largest pure-translation mode).
        if hub < 1.0:
            continue
        if axis not in picked:
            picked[axis] = (k, f, hub)
    return picked


def sweep():
    cases = [
        # arm_w, plate_t, sweep_deg, wiggle_amp, wiggle_lobes
        (1.0, 2.7, 220, 4.0, 4),     # was 1.019
        (1.0, 2.65, 220, 4.0, 4),
        (1.0, 2.6, 220, 4.0, 4),
        (1.0, 2.55, 220, 4.0, 4),
    ]
    print(f"{'w':>5} {'t':>5} {'sweep':>5} {'wA':>4} {'wL':>3}  "
          f"{'f_x':>7} {'f_y':>7} {'f_z':>7}  "
          f"{'f_xy/f_z':>8}")
    for w, t, sw, wa, wl in cases:
        p = SpringParams(arm_w=w, plate_t=t, arm_sweep_deg=sw,
                         wiggle_amp=wa, wiggle_lobes=wl)
        try:
            results = run(p, Material(), maxh_mm=1.5, n_modes=6)
        except Exception as e:
            print(f"{w:5.2f} {t:5.2f} {sw:5d} {wa:4.1f} {wl:3d}  ERROR: {e}",
                  file=sys.stderr)
            continue
        picked = first_three_translational(results)
        fx = picked.get("x", (None, float("nan"), 0))[1]
        fy = picked.get("y", (None, float("nan"), 0))[1]
        fz = picked.get("z", (None, float("nan"), 0))[1]
        ratio = ((fx + fy) / 2) / fz if fz else float("nan")
        print(f"{w:5.2f} {t:5.2f} {sw:5d} {wa:4.1f} {wl:3d}  "
              f"{fx:7.2f} {fy:7.2f} {fz:7.2f}  {ratio:8.3f}")


if __name__ == "__main__":
    sweep()
