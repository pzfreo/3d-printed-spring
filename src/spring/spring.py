"""Parametric orthoplanar spring with central proof-mass hub.

The spring is one solid part. The outer rim is the clamped boundary and has a
downward boss for a stiff clamp. The hub is a central disc with a downward
cylindrical extension that serves as the proof mass. N planar serpentine arms
of square cross-section connect rim to hub in the same plane.

Square cross-section (w == t) gives k_z ~= k_xy per arm; N >= 3 rotational
symmetry forces k_x == k_y. That makes the first three translational
eigenfrequencies (x, y, z) approximately equal — i.e. "isotropic resonance".
"""

from dataclasses import dataclass
import numpy as np
from build123d import (
    BuildLine, Circle, Pos, Spline, trace, extrude,
)


@dataclass
class SpringParams:
    """Default parameters tuned to FEA-isotropic (f_xy/f_z ≈ 1.02) for PLA.

    With these defaults, NGSolve modal analysis reports the first three
    translational modes within ~2% of each other (≈33 Hz each).
    """
    D_out: float = 60.0       # outer diameter (mm)
    D_rim_in: float = 50.0    # inner edge of clamp rim (mm)
    D_hub: float = 18.0       # hub / proof-mass diameter (mm)
    plate_t: float = 2.7      # arm thickness in z (mm)
    arm_w: float = 1.0        # arm width in plane  (mm)
    N_arms: int = 3
    arm_sweep_deg: float = 220.0
    rim_drop: float = 6.0     # depth of clamp boss below plate (mm)
    rim_rise: float = 0.0     # height of clamp boss above plate (mm)
    hub_drop: float = 12.0    # depth of proof-mass below plate (mm)
    hub_rise: float = 0.0     # height of proof-mass above plate (mm); 0 = no upper boss
    # Sinusoidal radial wiggle in the arm path. Lengthens the arm path so
    # bending dominates the in-plane stiffness (axial-stretch contribution
    # is what makes a naive square-section orthoplanar spring lateral-stiff).
    # Tapered by sin(πt) so the arm still meets hub/rim radially.
    wiggle_amp: float = 4.0   # radial amplitude (mm)
    wiggle_lobes: int = 4     # number of lobes along the arm
    # Optional heat-set insert pockets at the top and bottom of the hub.
    # Set insert_depth > 0 to subtract two coaxial cylinders. The default
    # 8.0 mm hole fits a standard M6 heat-set insert (e.g. Ruthex M6).
    insert_d: float = 8.0
    insert_depth: float = 0.0  # 0 = no pockets


def big_ring_support_free_params(
    D_out: float = 80.0, wiggle_amp: float = 6.0
) -> "SpringParams":
    """Support-free long-ringing variant.

    All features above the plate plane — the plate sits flat on the bed and
    everything extrudes straight up, no overhangs, no supports.

    The hub is a single tall column above the plate (~21 mm). Two M6
    heat-set inserts press into the same column from opposite ends:
        - top pocket: drilled from the top face of the hub
        - bottom pocket: drilled from the plate's underside, going up into
          the base of the hub (the slicer just leaves a small first-layer
          hole at the hub centre)

    The clamp uses a rim_rise (above the plate) instead of a rim_drop, so
    the user clamps the rim from above or grips the rim from the side.

    Mild physics cost: the PLA proof mass is now asymmetric about the plate
    plane (CoM offset above by ~10 mm), so lateral pushes couple a tiny
    amount into rocking. At ~10 Hz this is well below the ~80 Hz rocking
    modes so the coupling is small. If you bolt heavier mass on the bottom
    than the top, the combined CoM moves back toward the plate plane.
    """
    return SpringParams(
        D_out=D_out,
        D_rim_in=D_out - 10.0,
        D_hub=14.0,
        plate_t=2.7,
        arm_w=1.0,
        N_arms=3,
        arm_sweep_deg=220.0,
        wiggle_amp=wiggle_amp,
        wiggle_lobes=4,
        rim_drop=0.0,           # nothing below the plate
        rim_rise=6.0,           # clamp boss now sits above
        hub_drop=0.0,
        hub_rise=18.0,          # single tall hub; both M6 inserts fit
        insert_d=8.0,
        insert_depth=8.0,
    )


def big_ring_params(D_out: float = 80.0, wiggle_amp: float = 6.0) -> "SpringParams":
    """Preset for a long-ringing spring with M6 heat-set inserts both sides.

    Slim, tall hub keeps the PLA proof mass small (~3 g) so user-added
    steel/brass masses (10-15 g each side) dominate. Larger outer rim
    (default 80 mm OD) gives a softer spring than the 33 Hz design,
    pushing the loaded resonant frequency below 10 Hz and the perceptible
    ring time toward 5 s at Q≈50. Same arm cross-section and wiggle topology
    as the isotropic 33 Hz design so f_xy/f_z stays near 1.

    Args:
        D_out: outer diameter in mm. Default 80 mm; raise toward 100 mm for
               an even softer spring (lower f, longer ring).
        wiggle_amp: radial wiggle amplitude in mm. Default 6 mm — scaled
               from the 33 Hz design (4 mm) by the larger radial gap.
    """
    # Scale only the outer/rim dimensions; the hub stays slim to keep PLA
    # proof mass small. Arm cross-section unchanged so isotropy tuning carries.
    return SpringParams(
        D_out=D_out,
        D_rim_in=D_out - 10.0,
        D_hub=14.0,           # 3 mm wall around an 8 mm M6 insert
        plate_t=2.7,
        arm_w=1.0,
        N_arms=3,
        arm_sweep_deg=220.0,
        wiggle_amp=wiggle_amp,
        wiggle_lobes=4,
        rim_drop=6.0,
        hub_drop=9.0,
        hub_rise=9.0,
        insert_d=8.0,
        insert_depth=8.0,
    )


def _arm_face(start_angle_rad: float, p: SpringParams, n_pts: int = 80):
    sweep_rad = np.radians(p.arm_sweep_deg)
    t_arr = np.linspace(0.0, 1.0, n_pts)
    # Smoothstep: zero slope at both ends → arm meets hub and rim radially.
    s = 3 * t_arr ** 2 - 2 * t_arr ** 3
    theta = start_angle_rad + sweep_rad * s
    # Small overlap into rim and hub so the boolean fuses cleanly.
    r_inner = p.D_hub / 2 - 0.2
    r_outer = p.D_rim_in / 2 + 0.2
    r = r_inner + (r_outer - r_inner) * t_arr
    if p.wiggle_amp > 0 and p.wiggle_lobes > 0:
        # Sinusoidal radial wiggle along the path, tapered at both ends so the
        # arm still meets hub/rim radially (sin(π·t) envelope).
        env = np.sin(np.pi * t_arr)
        r = r + p.wiggle_amp * env * np.sin(p.wiggle_lobes * np.pi * t_arr)
    pts = [(float(r[i] * np.cos(theta[i])), float(r[i] * np.sin(theta[i])))
           for i in range(n_pts)]
    with BuildLine() as ln:
        Spline(*pts)
    return trace(ln.line, line_width=p.arm_w)


def build_spring(p: SpringParams = SpringParams()):
    rim_face = Circle(p.D_out / 2) - Circle(p.D_rim_in / 2)
    hub_face = Circle(p.D_hub / 2)
    plane_face = rim_face + hub_face
    for k in range(p.N_arms):
        plane_face = plane_face + _arm_face(k * 2 * np.pi / p.N_arms, p)

    part = extrude(plane_face, amount=p.plate_t)
    if p.rim_drop > 0:
        part = part + Pos(0, 0, -p.rim_drop) * extrude(rim_face, amount=p.rim_drop)
    if p.rim_rise > 0:
        part = part + Pos(0, 0, p.plate_t) * extrude(rim_face, amount=p.rim_rise)
    if p.hub_drop > 0:
        part = part + Pos(0, 0, -p.hub_drop) * extrude(hub_face, amount=p.hub_drop)
    if p.hub_rise > 0:
        part = part + Pos(0, 0, p.plate_t) * extrude(hub_face, amount=p.hub_rise)

    if p.insert_depth > 0:
        insert_face = Circle(p.insert_d / 2)
        # Top pocket: drilled into the top face of the hub_rise (or plate if no rise)
        top_z = p.plate_t + p.hub_rise
        part = part - Pos(0, 0, top_z - p.insert_depth) * extrude(
            insert_face, amount=p.insert_depth
        )
        # Bottom pocket: drilled into the bottom face of the hub_drop (or plate if no drop)
        bot_z = -p.hub_drop
        part = part - Pos(0, 0, bot_z) * extrude(
            insert_face, amount=p.insert_depth
        )

    return part
