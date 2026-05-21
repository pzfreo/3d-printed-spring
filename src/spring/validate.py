"""Modal analysis of the orthoplanar spring with netgen + NGSolve.

Workflow follows build123d issue #297: build the part in build123d, write a
BREP file, load it into netgen via OCCGeometry, tag faces by geometric
proximity, mesh, then run a generalised eigenvalue solve in NGSolve for the
linear-elastic modes of the spring with the outer rim clamped.

All units are SI (m, kg, s, Pa). The geometry produced by build_spring()
is in mm; we scale it by 1e-3 inside netgen before meshing.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import tempfile
import math
import numpy as np

import build123d as bd
from netgen.occ import OCCGeometry, gp_Pnt
import ngsolve as ngs
from ngsolve.eigenvalues import PINVIT

from .spring import SpringParams, build_spring


# PLA bulk properties (FDM-printed PLA, conservative typical values).
@dataclass
class Material:
    E: float = 2.5e9     # Young's modulus (Pa)
    nu: float = 0.36     # Poisson ratio
    rho: float = 1240.0  # density (kg/m^3)


def _scale_to_m(shape):
    """Return a netgen OCC shape scaled from mm to m."""
    return shape.Scale(gp_Pnt(0.0, 0.0, 0.0), 1.0e-3)


def _tag_clamp_face(shape, p):
    """Name the rim's outer annular face 'clamp'.

    For the original symmetric design (rim_drop > 0), the clamp face is the
    bottom annular face at z = -rim_drop. For the support-free design
    (rim_drop = 0, rim_rise > 0), it's the top annular face at z = plate_t +
    rim_rise. Both are annular planar faces whose centroid sits on the axis.
    """
    if p.rim_drop > 0:
        z_target = -p.rim_drop * 1e-3
    elif p.rim_rise > 0:
        z_target = (p.plate_t + p.rim_rise) * 1e-3
    else:
        raise RuntimeError(
            "no rim_drop or rim_rise — no rigid clamp surface to tag"
        )
    tagged = 0
    for f in shape.faces:
        c = f.center
        r = math.hypot(c.x, c.y)
        if abs(c.z - z_target) < 1e-4 and r < 1e-3:
            f.name = "clamp"
            tagged += 1
        else:
            f.name = "free"
    if tagged == 0:
        raise RuntimeError("did not find clamp face — geometric tagging failed")
    return tagged


def write_brep(spring_part: bd.Part, path: Path) -> Path:
    bd.export_brep(spring_part, str(path))
    return path


def load_geometry_si(brep_path: Path, p: SpringParams):
    """Load a BREP, scale to metres, tag the clamp face. Returns OCCGeometry."""
    raw = OCCGeometry(str(brep_path))
    shape_m = _scale_to_m(raw.shape)
    _tag_clamp_face(shape_m, p)
    return OCCGeometry(shape_m)


def mesh_spring(geo: OCCGeometry, maxh_mm: float = 1.2):
    """Generate a netgen mesh; maxh given in mm gets converted to m."""
    ngmesh = geo.GenerateMesh(maxh=maxh_mm * 1e-3)
    return ngs.Mesh(ngmesh)


def solve_modes(mesh: ngs.Mesh, mat: Material, n_modes: int = 8,
                order: int = 2, maxit: int = 60):
    mu = mat.E / (2 * (1 + mat.nu))
    lam = mat.E * mat.nu / ((1 + mat.nu) * (1 - 2 * mat.nu))

    fes = ngs.VectorH1(mesh, order=order, dirichlet="clamp")
    u, v = fes.TnT()

    def strain(w):
        return 0.5 * (ngs.grad(w) + ngs.grad(w).trans)

    def stress(w):
        e = strain(w)
        return 2 * mu * e + lam * ngs.Trace(e) * ngs.Id(3)

    a = ngs.BilinearForm(fes, symmetric=True)
    a += ngs.InnerProduct(stress(u), strain(v)) * ngs.dx
    a.Assemble()

    m = ngs.BilinearForm(fes, symmetric=True)
    m += mat.rho * (u * v) * ngs.dx
    m.Assemble()

    pre = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    lams, uvecs = PINVIT(a.mat, m.mat, pre=pre, num=n_modes,
                         maxit=maxit, printrates=False, GramSchmidt=True)
    return fes, uvecs, [float(x) for x in lams]


def classify_modes(fes, uvecs, lams, p: SpringParams):
    """Return per-mode (idx, f_Hz, dominant_axis, hub_unit_dir, |hub_disp|).

    For each eigenvector we copy it into a GridFunction and probe the
    displacement at the centre of the proof mass to decide which translation
    axis (x, y, z) it most resembles.
    """
    mesh = fes.mesh
    # Probe the centre of the solid PLA gap between the two insert pockets
    # (or just mid-plate if there are no pockets).
    if p.insert_depth > 0:
        top_of_bottom_pocket = -p.hub_drop + p.insert_depth
        bot_of_top_pocket = p.plate_t + p.hub_rise - p.insert_depth
        probe_z_mm = 0.5 * (top_of_bottom_pocket + bot_of_top_pocket)
    else:
        probe_z_mm = p.plate_t / 2.0
    probe = mesh(0.0, 0.0, probe_z_mm * 1e-3)
    gfu = ngs.GridFunction(fes)
    out = []
    for k, lam in enumerate(lams):
        f = math.sqrt(max(float(lam), 0.0)) / (2 * math.pi)
        gfu.vec.data = uvecs[k]
        vec = gfu(probe)
        v = np.array([float(vec[0]), float(vec[1]), float(vec[2])])
        norm = float(np.linalg.norm(v))
        v_unit = v / norm if norm > 1e-15 else v
        axis = ["x", "y", "z"][int(np.argmax(np.abs(v_unit)))]
        out.append((k + 1, f, axis, v_unit, norm))
    return out


def run(p: SpringParams = SpringParams(), mat: Material = Material(),
        maxh_mm: float = 1.2, n_modes: int = 8, brep_path: str | None = None):
    part = build_spring(p)
    if brep_path is None:
        brep_path = tempfile.NamedTemporaryFile(suffix=".brep", delete=False).name
    write_brep(part, Path(brep_path))
    geo = load_geometry_si(Path(brep_path), p)
    mesh = mesh_spring(geo, maxh_mm=maxh_mm)
    print(f"mesh: {mesh.ne} elements, {mesh.nv} vertices, "
          f"boundaries={mesh.GetBoundaries()}")
    fes, uvecs, lams = solve_modes(mesh, mat, n_modes=n_modes)
    print(f"DOFs: {fes.ndof}, eigenvalues found: {len(lams)}")
    results = classify_modes(fes, uvecs, lams, p)
    print()
    print(f"{'mode':>4}  {'f (Hz)':>9}  axis  hub-direction (x, y, z)        |hub|")
    for k, f, axis, vu, n in results:
        print(f"{k:>4}  {f:9.2f}   {axis}    "
              f"({vu[0]:+.3f}, {vu[1]:+.3f}, {vu[2]:+.3f})   {n:.2e}")
    return results


if __name__ == "__main__":
    run()
