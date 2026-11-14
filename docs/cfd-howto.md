# Running the 3D CFD on the prepared AHMU1218142 geometries

This is a step-by-step guide for running the six prepared lumens through a 3D solver, and for checking
the results against the reduced-order model and the credibility plan. It covers ANSYS Fluent (the most
common university license), COMSOL (the CFD Module) and OpenFOAM (free).

The 3D model solves the Navier–Stokes equations in the patient lumen. Like the app, it computes the
hemodynamics of a **given** lumen. The occluded and recanalized lumens are prescribed (see
[cfd-variants-AHMU1218142.md](cfd-variants-AHMU1218142.md)), so the CFD does not simulate a stent,
balloon or retriever acting on clot or wall. That would need fluid–structure and contact mechanics,
which this project does not attempt.

## 0. What you have

The prepared geometries are in `data/cmha/prepared/AHMU1218142/variants/<variant>/`, where
`<variant>` is one of:
`healthy`, `occluded_icas70`, `occluded_thrombus95`, `recanalized_thrombectomy`,
`recanalized_angioplasty` or `recanalized_stent`. Each folder contains:

| File | Use |
|---|---|
| `*_ext_in10D_out5D_mm.stl` | Whole surface, one named solid per patch (`wall`, `inlet`, `outlet_1`, `outlet_2`), in mm (Fluent, COMSOL) |
| `*_ext_in10D_out5D_m.stl` | Same, in metres (OpenFOAM) |
| `*_{wall,inlet,outlet_1,outlet_2}_mm.stl` | One patch per file, for importers that ignore solid names |
| `cfd_setup.json` | Inlet centre, normal and radius; outlet sizes; lesion throat; measurement planes; starting outlet resistances; reduced-order predictions |

Regenerate them with `python3 scripts/cmha_variants.py AHMU1218142`. The data folder is git-ignored,
so it stays on your machine.

## 1. Two kinds of run

| Run | Inlet | Outlets | Variants | Answers |
|---|---|---|---|---|
| **A. Fixed flow** | Fully developed (parabolic) velocity, Q = 153 mL/min | Equal static pressure (0 Pa gauge; any value works for rigid, incompressible flow) | healthy, all three recanalized | ΔP and WSS at the same flow. Calibrates the outlet resistances for B. |
| **B. Territory (primary endpoint)** | Total pressure 85 mmHg | p = 10 mmHg + R·Q at each outlet, with R from run A | all six | Distal flow as a % of healthy, flow restored, flow split, ΔP, WSS |

Do **not** run the two occluded lumens at a fixed 153 mL/min. They would need 149 mmHg and
2 × 10⁵ mmHg of driving pressure, which is not physiological. Run them only with the territory boundary
condition, which is the boundary condition the app uses.

## 2. Fluid and flow regime (all solvers)

* Blood: Newtonian, ρ = 1060 kg/m³, μ = 0.0035 Pa·s (ν = 3.30 × 10⁻⁶ m²/s). To reproduce the
  dataset authors exactly, use μ = 0.00345 Pa·s; their paper states 0.004, but their case file uses 0.00345.
* Steady, laminar, rigid no-slip walls.
* Reynolds numbers expected from the reduced-order model:

  | Location | Re |
  |---|---|
  | Healthy M1 | 371 |
  | Stent residual throat | 539 |
  | ICAS 70% throat | 617 |
  | Thrombus 95% throat | 79 (almost no flow) |

  For the ICAS 70% case, steady residuals may stall or oscillate. A stenotic jet at Re ≈ 600 can be
  disturbed downstream (the app flags this). If that happens, run a transient laminar case with a time step
  of about 1 × 10⁻⁴ s and average over the last 0.5 s.
* Inlet profile for the 10 D variants. With the inlet centre **c**, inward normal **n** and radius R from
  `cfd_setup.json` (R = 1.3946 mm, A = 6.1102 mm²):
  u(r) = 2·V̄·(1 − r²/R²) along **n**, where r is the distance from the inlet axis. For Q = 153 mL/min,
  V̄ = 0.4173 m/s and u_max = 0.8347 m/s. Use the cap centre (`extension_end_center_mm`) as the profile
  centre.

## 3. ANSYS Fluent (2022 R1 or later)

**Mesh (Fluent Meshing, Watertight Geometry workflow)**

1. Import `*_mm.stl` with units mm. The four solids come in as face zones `wall`, `inlet`, `outlet_1`
   and `outlet_2`.
2. Add local sizing: a global maximum of 0.15 mm, a curvature/proximity minimum of 0.02 mm, and a body
   of influence around the lesion (planes in `cfd_setup.json`, ±2 mm). Aim for at least 10 cells across
   the throat:
   * ICAS 70%: throat 0.78 mm, so cells ≤ 0.08 mm
   * thrombus 95%: throat 0.13 mm, so cells ≤ 0.013 mm (locally)
3. Boundary layers on `wall`: 5 layers, first height 0.02 mm, growth 1.2.
4. Volume: poly-hexcore. Check that the minimum orthogonal quality is above 0.2.

**Solver**

5. Pressure-based, steady, laminar. Material `blood` with the properties above.
6. Inlet:
   * **Run A:** velocity inlet. Use the parabolic profile from section 2 as an expression, for example
     `2*0.4173[m/s]*(1 - (r_axis/1.3946[mm])^2)`, where `r_axis` is a named expression for the distance
     from the inlet axis.
   * **Run B:** pressure inlet, total gauge pressure `85[mmHg]`.
7. Outlets:
   * **Run A:** pressure outlet, 0 Pa.
   * **Run B:** pressure outlet with
     `p_out1 = 10[mmHg] + R1 * (-MassFlow(['outlet_1'])/1060[kg m^-3])`, and the same for `outlet_2`.
     Fluent reports outlet mass flow as negative, hence the minus sign.
   * Start from R1 = 5.03 × 10⁹ Pa·s/m³ and R2 = 1.68 × 10¹⁰ Pa·s/m³ (`territory_bc_initial`). If the
     explicit coupling oscillates, under-relax the outlet pressures, or ramp R from 0 over 200 iterations.
8. Calibrate R from run A on the healthy lumen. Read the area-averaged inlet pressure p_in and the outlet
   flows Q₁ and Q₂, with both outlets at 0. Then set Rᵢ = (85 mmHg − p_in − 10 mmHg)/Qᵢ, converted to SI.
   This makes the healthy territory run return exactly 153 mL/min with the healthy split. Use these
   recalibrated R values for every variant, unchanged.
9. Methods: Coupled, second-order pressure and momentum. Converge to scaled residuals ≤ 10⁻⁶, with flat
   monitors for inlet pressure, both outlet flows and the plane-averaged pressures.

**Reports to extract (each variant)**

* Distal flow Q₁ + Q₂ in mL/min, as a % of the healthy territory run, and the share Q₁/(Q₁ + Q₂).
* ΔP across the lesion: the area-averaged static pressure on the `upstream` and `downstream` planes
  (point and normal in `cfd_setup.json`), giving ΔP = p_up − p_down and Pd/Pa = p_down/p_up.
* Wall shear stress between those planes: area-weighted mean and maximum, and the wall area with
  WSS < 0.4 Pa (stagnation) and WSS > 10 Pa (elevated shear), matching the app's safety endpoints.
* Recirculation: the length downstream of the throat over which the axial velocity near the wall is
  negative.

## 4. COMSOL (Laminar Flow, CFD Module)

1. Import the `*_mm.stl` mesh. Then use **Create Geometry from Mesh**, using "Minimal" face detection
   so the caps stay separate faces. Identify the inlet and outlets with Ball selections at the cap
   centres in `cfd_setup.json`.
2. Laminar Flow, stationary. Blood as above.
3. Inlet:
   * **Run A:** Inlet → *Fully developed flow*, flow rate 2.55e-6 m³/s. COMSOL builds the parabolic
     profile itself.
   * **Run B:** Inlet → *Pressure* 85[mmHg] (total).
4. Outlets (run B): Outlet → *Pressure*, p₀ = `10[mmHg] + R1*intop1(u*nx+v*ny+w*nz)`, where `intop1`
   is an Integration coupling operator on that outlet face. COMSOL solves this fully coupled, so the
   resistance boundary condition is stable without relaxation.
5. Mesh: Physics-controlled "Finer", plus a Size node on the lesion faces (maximum 0.08 mm for ICAS 70%,
   0.015 mm for 95%) and boundary layers (5, stretching 1.2).
6. Results: Cut Planes at the measurement planes, giving average pressure; surface integrals of
   `spf.Tw` on the lesion wall for WSS.

## 5. OpenFOAM (free)

1. Use `*_m.stl`. Make a background mesh with `blockMesh` around the bounding box, then run
   `surfaceFeatureExtract` and `snappyHexMesh` with the patches named from the STL solids. Add 5
   surface layers and lesion refinement level 5–6.
2. Use `simpleFoam` (steady, incompressible, laminar; `transportProperties` ν = 3.302e-6).
3. Inlet:
   * **Run A:** `flowRateInletVelocity` with `volumetricFlowRate 2.55e-6` and `profile laminarBL`.
     Alternatively, use the 22 D variant from `cfd-geometry-AHMU1218142.md` with a plug inlet.
   * **Run B:** `totalPressure` with p₀ = 85 mmHg/ρ in kinematic units.
4. Outlets (run B): use a `codedFixedValue` patch setting p = (10 mmHg + R·Q)/ρ, where Q is the
   patch flux sum `sum(phi.boundaryField()[patchi])`.
5. Post-process with `surfaces` function objects on the measurement planes and `wallShearStress`.

## 6. Mesh independence (needed for credibility)

For each solver, run the healthy and ICAS 70% variants on three meshes with a refinement ratio of at
least 1.3. Compute the observed order and the GCI for distal flow, ΔP across the lesion and peak WSS,
using the same method as `python -m bvm verify` (Roache, safety factor 1.25). Report the GCI with every
3D number. This is credibility evidence category 8 for model B.

## 7. Comparing the 3D results with the reduced-order predictions

The 1D predictions are in [cfd-variants-AHMU1218142.md](cfd-variants-AHMU1218142.md) and in each
`cfd_setup.json`:

| Variant | 1D distal flow (territory) | 1D ΔP along M1 at that flow |
|---|---|---|
| healthy | 153.0 mL/min (100%) | 0.9 mmHg |
| occluded_icas70 | 75.3 mL/min (49%) | 38.5 mmHg |
| occluded_thrombus95 | 1.6 mL/min (1%) | 74.2 mmHg |
| recanalized_thrombectomy | 153.0 mL/min (100%) | 0.9 mmHg |
| recanalized_angioplasty | 146.4 mL/min (96%) | 4.1 mmHg |
| recanalized_stent | 148.5 mL/min (97%) | 3.1 mmHg |

For each variant, report the comparison error E = S_3D − S_1D, alongside the 3D GCI. Expect agreement for
the healthy and recanalized lumens. Expect the 1D model to under-predict the ICAS 70% loss, because the
Young–Tsai coefficient comes from idealized axisymmetric stenoses, not an irregular patient lumen. The
size of that gap is the model-form result.

Both are models, so this is code-to-code evidence (credibility categories 5 and 7, supplementary).
Validation still needs measurements: the blood-vessel-mimic flow loop in the
[credibility plan](credibility-assessment-plan.md).

## 8. Checklist before trusting a run

- [ ] The healthy run A inlet flow is 153.0 mL/min and the outlet split is near 0.77 (Murray estimate; the
      authors' CFD for this subject gave 0.769).
- [ ] The healthy run B reproduces 153 mL/min within 0.5% after R calibration.
- [ ] Residuals ≤ 10⁻⁶ and monitors flat, or for transient runs the averaged quantities change by < 1%
      over the last 20% of the averaging window.
- [ ] GCI < 5% for the reported endpoints on the finest mesh.
- [ ] Mass balance: Q_in − (Q₁ + Q₂) < 0.1% of Q_in.
- [ ] Viscosity and density match section 2 (and match what any comparison case used).
