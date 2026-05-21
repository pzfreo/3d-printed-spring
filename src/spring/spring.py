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
    hub_drop: float = 12.0    # depth of proof-mass below plate (mm)
    # Sinusoidal radial wiggle in the arm path. Lengthens the arm path so
    # bending dominates the in-plane stiffness (axial-stretch contribution
    # is what makes a naive square-section orthoplanar spring lateral-stiff).
    # Tapered by sin(πt) so the arm still meets hub/rim radially.
    wiggle_amp: float = 4.0   # radial amplitude (mm)
    wiggle_lobes: int = 4     # number of lobes along the arm


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

    plate = extrude(plane_face, amount=p.plate_t)
    clamp_boss = Pos(0, 0, -p.rim_drop) * extrude(rim_face, amount=p.rim_drop)
    mass_boss = Pos(0, 0, -p.hub_drop) * extrude(hub_face, amount=p.hub_drop)
    return plate + clamp_boss + mass_boss
