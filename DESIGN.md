# Design notes — isotropic orthoplanar spring

This document explains *what* we built and *why* it works. The goal: a single
3D-printable PLA part that behaves as an isotropic resonator — a proof mass
suspended such that its three translational natural frequencies (kₓ, kᵧ, k_z)
are equal, so it rings at the same frequency regardless of which direction
you push it.

## 1. The target, precisely

For a small mass `m` suspended on a spring with stiffness `k` in some
direction, the natural frequency is

```
f = (1 / 2π) · √(k / m)
```

If a proof mass `m` is suspended on the same physical structure but feels
three different stiffnesses `kₓ`, `kᵧ`, `k_z` along the three axes, then it
rings at three different frequencies. "Isotropic resonance" means we want
**kₓ = kᵧ = k_z**, so all three natural frequencies coincide. This is the
useful property for suspensions, accelerometers, vibration isolators, and
gravity-wave-style payload mounts.

We restricted the topology to an **orthoplanar spring**: a single planar disc
of material where compliant arms in the disc's plane allow a central hub to
translate out of plane. This is the topology you get from a stamped washer
spring, and it's the natural fit for 3D printing flat.

## 2. Why orthoplanar springs are naturally anisotropic

Treat each arm as a long curved cantilever beam with rectangular cross
section: **width `w`** in the disc's plane (perpendicular to the arm's path)
and **thickness `t`** out of plane (the `z` direction).

A force applied at the hub end of one arm, with the rim end clamped, deforms
the arm in different ways depending on the load direction:

* **Out-of-plane load (z):** the only response is bending about the in-plane
  axis. Second moment of area is `Iₓᵧ = w·t³ / 12`. Bending compliance scales
  as `L³ / (E · Iₓᵧ)`, so the per-arm stiffness is

    `k_z ∝ E · w · t³ / L³`

* **In-plane load (x or y):** decomposes along the arm's path into an
  *axial* component (stretching the arm along its length) and a *transverse*
  component (bending the arm in-plane about the out-of-plane axis,
  `I_z = t · w³ / 12`).
  * Axial stretching stiffness: `E · A / L = E · w · t / L`
  * In-plane bending stiffness: `E · I_z / L³ = E · t · w³ / L³`

The two in-plane contributions are in series along the arc, so total
in-plane stiffness is somewhere between them — but for the small `L`
typical of a planar disc, **axial stretching dominates** because it scales
as `1/L`, not `1/L³`.

So the per-arm ratio is approximately

```
k_xy / k_z  ≈  (w / t)²   +   (axial-stretch contribution)
```

The first term is the bending-only ratio: a **square cross-section (w = t)
makes the bending contributions equal**, which was the starting intuition.
The second term is what gets you. For any reasonable arm length, axial
stretching is *much* stiffer than bending, so any in-plane component along
the arm's tangent direction adds substantial extra lateral stiffness.

With N ≥ 3 arms arranged with N-fold symmetry, kₓ = kᵧ falls out for free.
The hard target is **k_z = kₓ.**

## 3. Two design moves to kill the anisotropy

We have two knobs to lower `k_xy / k_z` toward 1:

**A. Make the cross-section taller than wide (`t > w`).** This drops the
bending part of the ratio because it goes as `(w/t)²`. With t = 2.7 mm and
w = 1.0 mm the bending ratio is `(1/2.7)² ≈ 0.14`, so out-of-plane is
*nominally stiffer*, leaving headroom for the axial-stretch contribution
to bring the lateral total back up.

**B. Lengthen the arm with a wiggle.** Axial stretching stiffness scales as
`1/L`, but in-plane bending stiffness scales as `1/L³`. So with the *same*
arm cross-section, a *longer* arm reduces stretching more slowly than it
reduces bending — but more importantly, a wiggly path forces the local
tangent direction to oscillate, so on average the in-plane load decomposes
roughly half into stretching and half into transverse bending, which
massively reduces the *effective* axial-stretch share.

Concretely, the arm path is

```
θ(t) = θ₀ + 220° · smoothstep(t)
r(t) = r_inner + (r_outer − r_inner)·t + A · sin(π·t) · sin(4π·t)
```

where `smoothstep(t) = 3t² − 2t³` ensures the arm meets hub and rim
radially (no kink), and the `sin(π·t)` envelope tapers the radial wiggle
amplitude `A = 4 mm` to zero at both ends. The result is a 4-lobe snake
that goes in-and-out four times as it sweeps 220° around.

These two moves, together, push `k_xy / k_z` from 1.54 down to 1.02.

## 4. FEA validation workflow (per build123d issue #297)

There's no closed-form solution for an N-armed wiggle-path orthoplanar
spring, so we used finite-element modal analysis to validate the tuning.

### Pipeline

1. **Geometry** in `build123d` → write to **BREP** (`bd.export_brep`).
2. **Load** into `netgen.occ.OCCGeometry(brep_path)`.
3. **Scale** the shape from millimetres to metres so all subsequent
   stiffness/mass quantities can be specified in plain SI units
   (`shape.Scale(gp_Pnt(0,0,0), 1e-3)`).
4. **Tag** the bottom annular face of the clamp boss with
   `face.name = "clamp"` by matching it geometrically:
   it is the only planar face at z = −rim_drop whose centroid sits on the
   axis. All other faces are named `"free"`.
5. **Mesh** with `geo.GenerateMesh(maxh=1.2e-3)` (max element size 1.2 mm).
6. **Linear elasticity** in NGSolve: assemble stiffness `K` and mass `M`
   over a vector H¹ space with Dirichlet boundary `"clamp"`:
   ```python
   fes  = ngs.VectorH1(mesh, order=2, dirichlet="clamp")
   u,v  = fes.TnT()
   def strain(w): return 0.5*(ngs.grad(w) + ngs.grad(w).trans)
   def stress(w):
       e = strain(w)
       return 2·μ·e  +  λ·Trace(e)·I₃
   a += InnerProduct(stress(u), strain(v)) * dx   # stiffness
   m += ρ · (u·v) * dx                            # mass
   ```
   with Lamé parameters `μ = E / 2(1+ν)`, `λ = E·ν / ((1+ν)(1−2ν))`.
   PLA constants: `E = 2.5 GPa`, `ν = 0.36`, `ρ = 1240 kg/m³`.
7. **Eigensolver** `PINVIT(a.mat, m.mat, pre=…, num=8)` (preconditioned
   inverse iteration with Cholesky preconditioning) returns eight smallest
   eigenvalues `λᵢ = ωᵢ²` and eigenvectors.
8. **Classify modes** by probing the eigenvector at the centre of the proof
   mass `(0, 0, −hub_drop/2)`: the dominant axis of the resulting
   displacement vector tells us which translation the mode corresponds to.

### The Python-3.12 caveat

`cadquery-ocp` (which build123d uses) and `netgen-occt` both link against
OpenCascade. On Python 3.9 they pin different OCCT versions and clash on
import. On Python 3.12 both ship OCCT 7.8.1, so they coexist in one
process — no subprocess hopping required. That's what the `uv run --python
3.12 …` incantation in the README is enforcing.

## 5. The tuning sweep

`scripts/sweep_isotropy.py` iterates over `(w, t, sweep°, wiggle_amp,
wiggle_lobes)` combinations, runs the full FEA for each, and reports
`f_x`, `f_y`, `f_z`, and the ratio `(f_x+f_y)/(2·f_z)`. Selected results:

| w (mm) | t (mm) | sweep° | wiggle A | lobes | f_x | f_y | f_z | f_xy/f_z |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2.0 | 2.0 | 220 | 0 | 0 | 78.1 | 78.1 | 50.7 | **1.539** |
| 1.0 | 2.0 | 220 | 0 | 0 | 32.3 | 32.3 | 26.1 | 1.238 |
| 0.8 | 2.0 | 300 | 0 | 0 | 16.0 | 16.0 | 13.7 | 1.164 |
| 1.0 | 2.5 | 220 | 4.0 | 4 | 32.1 | 32.1 | 30.7 | 1.044 |
| 1.0 | 2.5 | 220 | 4.0 | 5 | 34.7 | 34.7 | 35.5 | **0.977** |
| 1.0 | 2.7 | 220 | 4.0 | 4 | 33.2 | 33.2 | 32.6 | **1.019** |

Plain square-section S-arms ran at ratio 1.5 (lateral much too stiff).
Reducing `w` alone hit diminishing returns around 1.16. Adding the radial
wiggle was the move that broke through, with the (t = 2.7, wiggle 4×4) row
landing on the chosen design point at 1.7 % off isotropic.

## 6. Mesh-convergence sanity check

The chosen design point was re-meshed at two resolutions:

| maxh (mm) | elements | DOFs | f_z (Hz) | f_x (Hz) | f_y (Hz) |
|---:|---:|---:|---:|---:|---:|
| 1.5 | 14 545 | 83 475 | 32.59 | 33.19 | 33.20 |
| 1.2 | 22 717 | 127 032 | 32.63 | 33.17 | 33.18 |

Frequencies move by under 0.2 Hz between mesh sizes — the result is mesh
converged.

## 7. Higher modes & caveats

The first six FEA modes at the chosen design point:

| mode | f (Hz) | shape |
|:---:|:---:|---|
| 1 | 32.63 | translation in **z** |
| 2 | 33.17 | translation in **x** |
| 3 | 33.18 | translation in **y** |
| 4 | 74.3 | plate flexure (hub barely moves) |
| 5 | 103.5 | proof-mass rocking about y |
| 6 | 103.7 | proof-mass rocking about x |

Modes 1–3 are degenerate to within 1.7 %. Modes 5–6 (rocking about x/y) are
~3× higher and are *not* matched to modes 1–3 — they're a different physical
DOF (rotation, not translation). If your application also needs all three
rotational frequencies to coincide with the translational ones, this single
planar spring isn't enough; you'd need a multi-stage 3D suspension.

**Material:** the FEA uses an isotropic linear-elastic PLA at E = 2.5 GPa.
FDM-printed PLA is in fact anisotropic (along-layer ≈ 3.5 GPa, between
layers ≈ 2.0 GPa). The *isotropy ratio* between modes is geometric and so
is preserved; the *absolute frequencies* will shift ±15 % from material
property uncertainty alone.

**Print orientation matters.** Print the spring flat — disc lying on the
bed, bosses pointing down or up depending on your supports strategy — so the
arm's z direction (the 2.7 mm dimension) is across the printer's build
layers. The out-of-plane bending mode (mode 1) then stresses the arms
*along* the layers, which is the strong direction of FDM PLA. Printing
on edge would put layer-adhesion in the failure path of the dominant mode.

## 8. Long-ring variants

The 33 Hz reference design *rings* for a brief 0.3–1 s depending on the
clamp. For applications where you want a longer-ringing resonator — felt
as a slow visible wobble rather than a quick buzz — you want **low
frequency**, since the perceptible ring duration scales as `3 · Q / (π·f)`.

Empirically measured Q for the printed 33 Hz spring:

| Mount | Perceptible ring | Implied Q |
|---|---|---|
| Hand-held rim | ~1 s | ~35 |
| Rim pressed firmly against a fixed surface | ~2 s | ~70 |

PLA's intrinsic Q is the dominant loss. Material upgrades give modest
gains (PETG ~+30 %, polycarbonate ~+200 %), but the easiest knob is **f
itself**. The two new presets below preserve the isotropic arm topology
of section 3 but redesign the proof-mass + clamp geometry so the user
can attach external mass via M6 heat-set inserts, dragging f down into
the 5–15 Hz range.

### `big_ring_params()` — symmetric, needs supports

The 80 mm-OD scale-up. Slim, *tall* hub (14 mm diameter, 20.7 mm long)
centred on the plate plane: hub_rise = hub_drop = 9 mm. Two M6 pockets
(8.0 mm × 8.0 mm deep) at each end of the hub, drilled in opposite
directions. Slim hub keeps the PLA proof mass tiny (~3 g) so 20–30 g of
user-added steel dominates the modal mass.

FEA result with no external mass:

| mode | f (Hz) | axis |
|:--:|:--:|---|
| 1 | 28.4 | z |
| 2 | 29.0 | x |
| 3 | 29.0 | y |

Isotropy ratio: **f_xy / f_z = 1.020** (2 % off).

Projected with 25 g external mass: f ≈ 9.4 Hz, perceptible ring ≈ 5.1 s
(Q=50). At Q=70: ≈ 7.1 s.

**Print caveat:** because the hub extends both sides of the plate, the
plate has to bridge across either the hub-rise *or* the hub-drop in mid
air. Supports are required under one of the two hub bosses no matter
which orientation you choose.

### `big_ring_support_free_params()` — single-sided hub, prints flat

Everything above the plate plane. `hub_drop = 0`, `rim_drop = 0`. Both
the clamp boss (`rim_rise = 6 mm`) and the proof-mass column
(`hub_rise = 18 mm`) sit on top of the plate. The plate itself is the
bottom face of the part and lies flat on the bed during printing.

The two M6 pockets are in the same upward hub column, drilled from
opposite ends: the top pocket from the top face downward, the bottom
pocket from the *plate's underside* upward into the base of the hub
(the slicer simply leaves a small 8 mm hole in the first layer at the
hub centre — the pocket is just a vertical cylinder going up). 4.7 mm
of solid PLA isolates the two pockets.

This makes the proof mass asymmetric about the plate plane (its CoM
sits ~10 mm above the plate, not at it), which couples a small amount
of lateral push into rocking. In practice the rocking modes are at
~80 Hz — much higher than the loaded ~10 Hz translation — so the
coupling is weak.

FEA result with no external mass:

| mode | f (Hz) | axis |
|:--:|:--:|---|
| 1 | 26.9 | x |
| 2 | 26.9 | y |
| 3 | 28.2 | z |

Isotropy ratio: **f_xy / f_z = 0.954** (4.6 % off — slightly *softer*
in plane than out of plane, the opposite tilt from the symmetric
version).

Projected with 25 g external mass: f ≈ 8.9 Hz, perceptible ring ≈ 5.4 s
(Q=50). At Q=70: **7.5 s**. With 30 g: 5.8 s / 8.2 s respectively.

### Effective-mass approximation

The post-FEA scaling `f_loaded = f_PLA · √(m_hub_PLA / (m_hub_PLA + m_ext))`
treats the hub as a rigid proof mass. This is a rigid-body approximation:
it assumes the hub moves coherently while the arms flex, and that the
arms' own kinetic energy is negligible. For these designs the proof mass
dominates (the arms together weigh ~1 g), so the approximation should be
accurate to a few percent. If you want a stricter prediction, model the
external mass as a denser cylinder inside the FEA — the NGSolve setup
allows region-dependent mass density via a piecewise `CoefficientFunction`.

### Heat-set insert sizing

The 8.0 mm × 8.0 mm pockets are sized for **thermal install** with most
M6 heat-set inserts: the soldering iron heats the brass, the PLA melts
around the knurling, and the displaced material locks the insert in
place. Tolerances by brand:

| Brand / spec | Pocket vs. insert |
|---|---|
| Ruthex M6 (8.1 mm OD × 8.1 mm long) | 0.1 mm tight — fine via thermal install, insert sits ~0.1 mm proud |
| Voron-spec M6 (8.0 × 8.0) | Exact match |
| Generic brass knurled M6 (7.5–8.0 mm, 7–10 mm long) | Diameter fine; ≤8 mm long fits, longer needs depth bump |

For inserts ≥9 mm long, raise `insert_depth` to 10 mm and `hub_rise` to
20 mm to preserve the 4.7 mm solid-PLA gap between the two pockets.
