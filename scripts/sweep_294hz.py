"""Search for a spring tuned to ~294 Hz (D4) with kx=ky=kz preserved.

Reaching 294 Hz at 60 mm OD requires opposite moves to the isotropy
fix — shorter, stiffer arms — so we add the part radius as a sweep
parameter and look for a smaller spring that keeps the wiggle topology
but at smaller absolute scale.
"""
from __future__ import annotations
import sys
from spring.spring import SpringParams
from spring.validate import run, Material


def first_three_translational(results):
    picked = {}
    for k, f, axis, vu, hub in results:
        if hub < 1.0:
            continue
        if axis not in picked:
            picked[axis] = (k, f, hub)
    return picked


def sweep():
    cases = [
        # OD, rim_in, hub, plate_t, arm_w, sw, wA, wL, hub_drop, rim_drop
        # Uniform scales of the 33 Hz design (60 → 40 → 30 → 25 mm OD)
        (40.0, 33.0, 12.0, 1.8, 0.7, 220, 2.7, 4, 8.0, 4.0),
        (30.0, 25.0,  9.0, 1.4, 0.5, 220, 2.0, 4, 6.0, 3.0),
        # Same 60 mm OD but tiny proof mass + short stiff arms
        (60.0, 50.0, 18.0, 2.7, 1.0,  90, 1.5, 2, 0.0, 6.0),
        (60.0, 50.0, 18.0, 2.7, 1.0, 120, 2.0, 3, 0.0, 6.0),
        # Mid scale: 40 mm OD, full-thickness arms, no proof mass
        (40.0, 33.0, 12.0, 2.7, 1.0, 220, 3.0, 4, 0.0, 4.0),
        (40.0, 33.0, 12.0, 2.7, 1.0, 220, 3.5, 4, 4.0, 4.0),
        # 30 mm OD with full-thickness arms (so wall is printable)
        (30.0, 24.0, 10.0, 2.7, 1.0, 220, 2.0, 4, 4.0, 3.0),
        (30.0, 24.0, 10.0, 2.7, 1.0, 220, 2.5, 4, 0.0, 3.0),
    ]
    hdr = ("OD", "rin", "hub", "t", "w", "sw", "wA", "wL", "hd")
    print(f"{hdr[0]:>5} {hdr[1]:>5} {hdr[2]:>5} {hdr[3]:>4} {hdr[4]:>4} "
          f"{hdr[5]:>4} {hdr[6]:>4} {hdr[7]:>3} {hdr[8]:>5}  "
          f"{'f_x':>7} {'f_y':>7} {'f_z':>7}  {'f_xy/f_z':>8} {'Δ%':>6}")
    for c in cases:
        OD, rin, hub, t, w, sw, wA, wL, hd, rd = c
        p = SpringParams(
            D_out=OD, D_rim_in=rin, D_hub=hub,
            plate_t=t, arm_w=w, arm_sweep_deg=sw,
            wiggle_amp=wA, wiggle_lobes=wL,
            hub_drop=hd, rim_drop=rd,
        )
        try:
            res = run(p, Material(), maxh_mm=min(1.0, t*0.6), n_modes=6)
        except Exception as e:
            print(f"{OD:5.1f} {rin:5.1f} {hub:5.1f} {t:4.1f} {w:4.1f} "
                  f"{sw:4d} {wA:4.1f} {wL:3d} {hd:5.1f}  ERROR: {e}",
                  file=sys.stderr)
            continue
        picked = first_three_translational(res)
        fx = picked.get("x", (None, float("nan"), 0))[1]
        fy = picked.get("y", (None, float("nan"), 0))[1]
        fz = picked.get("z", (None, float("nan"), 0))[1]
        ratio = ((fx + fy) / 2) / fz if fz else float("nan")
        fmean = (fx + fy + fz) / 3
        delta = (fmean - 294.0) / 294.0 * 100
        print(f"{OD:5.1f} {rin:5.1f} {hub:5.1f} {t:4.1f} {w:4.1f} "
              f"{sw:4d} {wA:4.1f} {wL:3d} {hd:5.1f}  "
              f"{fx:7.1f} {fy:7.1f} {fz:7.1f}  {ratio:8.3f} {delta:+6.1f}%")


if __name__ == "__main__":
    sweep()
