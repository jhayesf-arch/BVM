# Credibility assessment plan

**Models:** (A) the BVM reduced-order MCA model (`bvm/`, `app/`) and (B) a planned patient-specific 3D CFD
model of a healthy MCA bifurcation.
**Framework:** FDA, *Assessing the Credibility of Computational Modeling and Simulation in Medical Device
Submissions*, final guidance, 2023 (docket FDA-2021-D-0980), which builds on ASME V&V 40-2018.
**Status:** draft plan. Evidence marked *done* is in the repository. Goals marked *proposed* need the
lab's agreement before evidence is gathered, because the guidance asks for goals to be set in advance.

This is academic research, not a regulatory submission. We use the guidance as an accepted structure
for deciding how much evidence the model needs and for reporting it. Nothing here makes a regulatory
claim.

---

## Step 1. Question of interest

> In a native middle cerebral artery, what percentage of the healthy reference distal MCA flow is
> restored after an intervention (angioplasty, stenting or thrombectomy) for M1 stenosis or occlusion
> of 50–95%?

Secondary questions use the same models: the pressure drop and Pd/Pa across the treated segment, how
flow divides between the distal (M2) territory and the lenticulostriate perforators, and the change in
wall shear stress (WSS). Safety screening covers low-shear (stagnation) length, elevated shear,
starved or covered perforators, and wall contact.

## Step 2. Context of use (COU)

| | Model A: reduced-order (this repo) | Model B: patient-specific 3D CFD (planned) |
|---|---|---|
| Role | Rapid screening of how geometry, wall stiffness, BC and device outcome change the endpoints, to choose blood vessel mimic (BVM) bench conditions and generate hypotheses | Resolve local flow features model A cannot: separation, recirculation, WSS distribution and outlet flow split at the bifurcation |
| What else informs the decision | BVM bench flow-loop measurements (primary); model B | BVM bench measurements; model A as a cross-check |
| Quantities of interest | Distal flow % of healthy; ΔP and Pd/Pa across lesion; perforator flow; mean and peak WSS | The same, plus WSS maps, recirculation length and the M2 flow split |
| Range | M1 diameter 1.5–5 mm; stenosis 0–95% (diameter); steady laminar flow, Re < 2000; wall E 0.01–10 MPa | One healthy control geometry (below) with a computational stenosis or thrombus and a recanalized lumen |
| Out of scope | Pulsatile effects, transition and turbulence, device–clot and device–wall mechanics, collaterals, autoregulation | Collaterals, autoregulation, device mechanics |

## Step 3. Model risk

**Model influence** is how much the answer rests on the model rather than on other evidence. It is
**low to medium**: the BVM bench measurements are the primary evidence, and the models choose
conditions and explain mechanisms.

**Decision consequence** is how much harm a wrong answer could cause. It is **low** for laboratory
research, where a wrong answer costs bench time but no patient is affected. It would become **high** if
results were ever used to support a device claim.

**Model risk: low to medium** for the current COU. Evidence below is planned to that level. Raising the
influence (for example, using the model *instead of* bench tests) or the consequence requires
reassessing this step and raising the credibility goals.

## Step 4. Credibility evidence, by the guidance's categories

| # | Category | Model A status | Artifact | Next action |
|---|---|---|---|---|
| 1 | Code verification | **Done** | `tests/test_verification.py`: exact Hagen–Poiseuille (1e-12); linear taper against the closed-form integral with observed order 2.0; rigid stenosis equal to the exact viscous integral plus Young–Tsai loss (1e-10); Laplace strain at the inlet; White correlation continuous at De = 11.6; territory BC equal to a series-resistance circuit; mass conservation with branches. The JS port matches Python to about 1e-7 to 1e-11. | Agreement between two ports of the same algorithm is weak evidence (they can share an error). It is kept as a regression check only. |
| 2 | Model calibration | **Done; weak by design** | Territory distal resistance set so the healthy reference carries the reference flow. Perforator resistances set to a literature-derived share. Young–Tsai K_t = 1.52 taken from Young & Tsai 1973. | Report as calibration, never as validation. |
| 3 | Bench test validation | **Planned: primary evidence** | Protocol below | Build stenosis test articles; run the flow loop |
| 4 | In vivo validation | Partial: published data only | MRI flows (Amin-Hanjani 2015); pressure wire (Shoemaker 2025) | Compare model Pd/Pa with published ICAS translesional pressure-ratio data at matched stenosis |
| 5 | Population-based validation | Partial | Healthy M1 geometry from 44 CTA controls (Song 2024) and 554 angiograms; flow from about 320 adults | Compare the *distribution* of model outputs over the 44 control geometries with population flow and Pd/Pa distributions |
| 6 | Emergent model behavior | Partial | The critical-stenosis curve emerges (flow holds until about 70% diameter stenosis, then collapses) without being imposed | The disturbed-flow flag is a threshold, not a prediction, and does **not** count. Model B should test transition onset, the guidance's own example. |
| 7 | Model plausibility | **Done** | Governing equations and assumptions in `README.md`; every anatomical input tagged measured, derived or assumed with a citation (`app/anatomy.js`) | Weak alone; keep current |
| 8 | Calculation verification / UQ | **Done** | `docs/verification-results.md` (`python -m bvm verify`): GCI under 0.4% for every endpoint at 401 stations; Monte Carlo propagation of input uncertainty with rank-correlation sensitivity | Rerun after any engine change. Add bench measurement uncertainty once known. |

**Model B geometry: done for AHMU1218142.** The surface is prepared and verified in
[cfd-geometry-AHMU1218142.md](cfd-geometry-AHMU1218142.md):
* watertight, manifold and consistently oriented, with the anatomy unchanged and volume within 0.2%;
* ports match the authors' Fluent zones to 0.07 mm, and the inflow is reproduced exactly;
* one input discrepancy found: the case file uses μ 0.00345 Pa·s, where the paper states 0.004;
* across 43 controls, Murray's D³ rule reproduces the authors' CFD flow split (r = 0.85,
  mean absolute error 0.05).

Occluded and recanalized lumens for the same subject are built and verified in
[cfd-variants-AHMU1218142.md](cfd-variants-AHMU1218142.md): every variant is watertight with no flipped
faces, is unchanged outside the lesion, and has a measured throat within 0.003 of the target. The
solver set-up, mesh-independence steps and the 1D-versus-3D comparison are in
[cfd-howto.md](cfd-howto.md).

Model B (3D CFD) needs its own category 1 and 8 evidence: solver verification, mesh independence by
GCI on at least three meshes, and outlet BC sensitivity. It then validates against the same bench data.

### Bench validation protocol (category 3, proposed)

* **Comparator:** BVM flow loop with test articles of straight M1-scale lumens (3.0 mm), printed
  cosine stenoses of 0, 50, 70 and 85% (diameter), 6 mm long, in substrates of measured modulus.
  A rigid set isolates the stenosis loss from wall compliance.
* **Measured quantities:** inlet flow (in-line flowmeter), pressure upstream and downstream of the
  stenosis (catheter-tip or side-port transducers at fixed distances), Pd/Pa, and, if PIV is
  available, centerline velocity and jet length.
* **Matched conditions:** steady flow at 60–180 mL/min. The blood analog has measured viscosity and
  density (Newtonian glycerol–water near 3.5 mPa·s and 1060 kg/m³), and temperature is logged. The
  model is run with the *measured* geometry, viscosity and flow, not the nominal ones.
* **Validation metric (ASME V&V 20 style):** comparison error E = S − D for each quantity, judged
  against the validation uncertainty u_val = √(u_num² + u_input² + u_D²). u_num comes from the GCI,
  u_input from propagating measured-input uncertainty, and u_D is the measurement uncertainty.
* **Applicability:** the test conditions span the COU (Re, stenosis severity, wall stiffness). Any
  gap, such as pulsatility or in vivo wall properties, is listed as a limitation.

## Step 5. Credibility factors and goals (proposed)

The gradation scales follow ASME V&V 40. Goals are set for **low-to-medium model risk**.

| Factor (V&V 40 unless noted) | Proposed goal for model A |
|---|---|
| Software quality assurance | Version control, automated test suite on every change (met) |
| Numerical code verification | Analytic solutions and observed order of accuracy (met) |
| Discretization error | GCI below 1% for every QoI at the default grid (met) |
| Numerical solver error | Fixed-point and root-finding tolerances at least 1e-6 below QoI changes of interest (met) |
| Use error | Inputs tagged with basis; worked scenarios reproducible from the URL hash (met) |
| Model form | Bench ΔP and flow within u_val for 0–70% stenosis; for 85%, document the bias rather than require agreement |
| Model inputs | Geometry and viscosity measured on each test article |
| Test samples | At least 3 articles per stenosis level |
| Test conditions | At least 3 flow rates spanning 60–180 mL/min |
| Equivalency of input parameters | Model run with the as-measured geometry, fluid properties and flow |
| Output comparison | Quantitative, with u_val, for ΔP, Pd/Pa and flow; WSS only qualitative (not measurable in the bench) |
| Relevance of QoIs to the COU | ΔP and flow are direct inputs to the primary endpoint |
| Relevance of validation to the COU | Bench range covers the COU's Re and stenosis range; pulsatility and collaterals are not covered |
| *Calibration goodness of fit* (custom) | Report the residual between the reference flow and its calibrated value (exact by construction) |
| *Population agreement* (custom) | Median of the model healthy-flow distribution within the published population IQR |

## Patient-specific geometry for model B

Source: Song et al. 2024, *Scientific Data*, CMHA dataset (Figshare 10.6084/m9.figshare.26965450, CC BY),
downloaded to `data/cmha/` (git-ignored; checksum verified). Screening of the 44 healthy-control STLs:

* All 44 are mm-scale manifold surfaces with exactly 3 open ports (one inlet, two M2 outlets). They must
  be capped and extended with straight inlet and outlet sections before CFD. Mean triangle edge is
  0.13–0.31 mm, so the surface should be refined before meshing.
* The segments are short (about 7–28 mm between port centres). The authors' protocol required only
  ≥ 6 mm inlet and ≥ 3 mm outlets, so flow extensions are needed.
* Each control also includes the CTA (NIfTI) and the authors' Fluent case (`fluent_healthy_artery_*.cas.gz`).
  Their case uses steady flow, μ 0.004 Pa·s, an age/sex-based inlet velocity and 100 mmHg at both outlets.
* Data quality notes: in 27 of 44 rows the "thicker branch" diameter is smaller than the "smaller
  branch" diameter. Outlet numbering does not follow branch size. The OSI column (1.3–4.7) is invalid,
  since OSI is bounded 0–0.5. The authors' CFD results are model outputs, so comparing against them is
  code-to-code comparison (supplementary only), not validation.
* First candidates have long segments and a near-median M1: **AHMU1218142** (M1 2.79 mm, branches
  2.62 and 1.75 mm, ports about 23 mm apart), AHMU1218126 and AHMU1218141. Note that `AHMU1218143`'s
  STL has a different file name (`healthy artery_AHMU1218143.stl`).
* Because the geometry is truncated at the bifurcation, a complete M1 occlusion in model B gives an
  upper bound on the deficit: collateral supply is absent, as it is in model A.

## Adequacy assessment (to complete after evidence)

For each factor: goal, achieved level, and whether the gap matters for the COU. Then a statement of
whether model A, model B or both are credible enough for the question of interest at the stated model
risk, or what further evidence is needed.

## Reporting

The guidance's Appendix 2 lists what a credibility report should contain. In this repository:

| Report item | Location |
|---|---|
| Question of interest, COU, model risk | This document, steps 1–3 |
| Model description, equations, assumptions, limitations | `README.md` (Model, Limitations) |
| Input sources and basis | `app/anatomy.js`; Anatomy tab in the app |
| Code verification | `tests/test_verification.py`, `docs/verification-results.md` |
| Calculation verification and UQ | `docs/verification-results.md` |
| Validation | Bench protocol above (results pending) |
| Credibility factors, goals, adequacy | Step 5 and the adequacy section |

## References

* FDA. *Assessing the Credibility of Computational Modeling and Simulation in Medical Device
  Submissions.* Final guidance, 2023. Docket FDA-2021-D-0980.
* ASME V&V 40-2018. *Assessing Credibility of Computational Modeling through Verification and
  Validation: Application to Medical Devices.*
* ASME V&V 20-2009. *Standard for Verification and Validation in Computational Fluid Dynamics and Heat
  Transfer.*
* Roache PJ. Perspective: a method for uniform reporting of grid refinement studies. *J Fluids Eng*
  1994;116:405–413.
* Song M et al. Intracranial aneurysm CTA images and 3D models dataset with clinical morphological and
  hemodynamic data. *Sci Data* 2024;11. doi:10.1038/s41597-024-04056-8.
* Young DF, Tsai FY. Flow characteristics in models of arterial stenoses. *J Biomech* 1973;6:395–410.
