# BVM: parametric 3D vessel model for ischemic disease

A computational companion to blood-vessel-mimic (BVM) work. It builds a 3D model of a native cerebral
artery segment (the middle cerebral artery M1 by default, plus four other stroke-prone sections with
literature geometry), solves steady hemodynamics on it, and exports CFD/print-ready geometry.
It is meant to support New Approach Methods (NAMs): you can sweep the same parameters you control
on the bench (lumen size, curvature, stenosis, substrate stiffness) and see how they change shear,
pressure, Reynolds number and distal perfusion.

## Interactive simulator (start here)

A browser app in `app/`. Move the sliders and the 3D vessel, flow animation and every readout update
instantly. There is no install and no build step: it needs only Python, to serve the files.

```bash
python3 -m http.server 8766 --directory app
```

Then open http://localhost:8766.

* **Pathology factors:** luminal diameter, tortuosity index, occlusion % (by diameter or area) and
  substrate stiffness (log scale from 10 kPa to 10 MPa, or a rigid wall).
* **3D viewport:** drag to orbit and scroll to zoom. The lumen can be colored by wall shear stress,
  pressure, velocity, Reynolds number or wall strain. Animated blood particles follow the local
  parabolic (Poiseuille) velocity profile in 40× slow motion, so they visibly speed up through a
  stenosis. A probe ring, moved with its slider or by clicking a chart, shows local values.
* **Reynolds number:** shown at the inlet and at the stenosis throat on a log-scale bar with the
  laminar (< 2000), transitional (2000–4000) and turbulent (≥ 4000) zones, plus a post-stenotic
  disturbance warning.
* **Flow restoration (top of the results):** in Brain territory mode, territory flow before and
  after the selected device, each as a % of the section's healthy flow, with normal, oligemia,
  penumbra-range and core-range bands. *Restored* is the share of the flow deficit recovered,
  (after − before) ÷ (healthy − before). This is volumetric flow in the modeled vessel, not an
  angiographic eTICI grade.
* **Charts:**
  * *Along the vessel* plots diameter, pressure, Re and WSS along the vessel.
  * *Sensitivity* sweeps any factor and plots a chosen output, for example the critical-stenosis
    flow curve.
  * *Sensitivity → Aspiration* compares RED 72 and ZOOM 71 head to head. Choose pump vacuum,
    catheter length or inner diameter as the x-axis, and aspiration flow, tip force index or catheter
    Re as the output. The estimate is Hagen–Poiseuille for a uniform bore, so flow scales with ID⁴:
    RED 72's extra 0.001″ gives about 5.8% more flow. Tubing, the clot and proximal ID steps are
    ignored, so it is a comparative index, not bench data.
  * The sensitivity sweep can also vary the *upstream pressure loss* and *taper*, and, with a device
    selected, plot *flow after device* or *flow restored* against any factor.
* **Scenarios:** grouped by anatomical section (below), plus the M1 parameter studies (50/70/85%
  stenosis, tortuous, narrow, stiffened, soft hydrogel).

### Anatomy tab: five stroke-prone sections

The left panel's **Anatomy** tab swaps the vessel for one of five sections. Each has one or more
geometry sources, and every value is labeled **measured** (from the cited study), **derived**
(computed from measured values, with the formula shown) or **assumed** (our choice, editable). Edited
values are highlighted against the source. The shape controls add **taper** (outlet ÷ inlet
diameter), a cosine **bulb** dilation and wall thickness. Slider ranges follow the section.

| Section | Geometry sources | Healthy flow | Section scenarios |
|---|---|---|---|
| MCA M1 | Lab baseline (3.0 mm × 20 mm, 120 mL/min); Żurada 2011 CTA (2.23 × 15.6 mm); 554 angiograms (2.39 × 15.6 mm) | 153 mL/min, the mean of L/R MCA in Amin-Hanjani 2015 (PC-MRA, ~320 adults) | embolic occlusion, 75% / 85% ICAS |
| Lateral lenticulostriate | 0.2 / 0.4 / 0.8 mm, within the 0.10–1.28 mm cast range (Marinković) | derived: 4.5 cm/s (7T flow MRI, Bouvy 2016) × lumen area ≈ 0.34 mL/min at 0.4 mm | lipohyalinosis (diffuse 50%), origin microatheroma 70%, branch occlusion |
| ICA bifurcation & carotid bulb | CCA 6.1 mm (ultrasound) tapering by the ICA/CCA ratio 0.81 (young) or 0.77 (older) from Thomas 2005; bulb +30% *assumed* | ICA 257 mL/min (Amin-Hanjani 2015) | bulb plaque 50 / 70 / 90% |
| Watershed distal cortical branch | 1.5 mm (angular) or 0.8 mm (temporopolar) cortical branch outer diameters; taper 0.6 *assumed* | derived by Murray's law from M1 | tortuous course, hypotension −25 / −40 mmHg, distal embolus |
| Basilar, proximal & mid | Vitošević 2019 CTA (3.17 × 32.4 mm), MRA n = 100 (3.06 × 31.3 mm), cadaver (3.05 × 25.6 mm) | 138 mL/min (124 at 61–80 y), Amin-Hanjani 2015 | mid-basilar ICAS 70%, thrombosis, dolichoectasia (> 4.5 mm, Smoker 1986) |

Notes and limits:

* **One flow path per section.** The ECA at the carotid bifurcation, the perforators off M1 and the
  basilar, and the collaterals are not modeled as branches.
* **Watershed inlet pressure.** It is assumed at 65 mmHg. A 2025 pressure-wire feasibility study
  (Shoemaker et al., n = 5, under anesthesia) measured 62 ± 29 mmHg in the MCA against 80 ± 18 mmHg
  radial. The **upstream pressure loss** input then models hypotension. It is applied after the
  territory is calibrated, and there is no autoregulation, so this is a worst case.
* **Stenosis % is relative to the local healthy diameter.** In the carotid bulb it is therefore
  ECST-like, not NASCET.
* **Bulb flow needs 3D CFD.** Flow separation and recirculation in the bulb are not captured by
  this model. Export the STL for 3D CFD.
* **Export:** STL (mm, with wall/inlet/outlet patches), CSV profiles, JSON, PNG, and a shareable link.
  The link encodes the current settings in the URL.

### Device panel (top of the 3D view)

* **Lesion switch.** Choose **Plaque** (chronic intracranial atherosclerotic stenosis, ICAS) or
  **Clot** (acute thrombotic occlusion).
* **Devices offered for each lesion:**

  | Lesion | Category | Devices |
  |---|---|---|
  | Clot | Stent retrievers | Trevo NXT (Stryker), Solitaire X (Medtronic), EMBOTRAP III (Cerenovus), Tigertriever (Rapid Medical, adjustable diameter) |
  | Clot | Aspiration catheters | RED 72 (Penumbra), ZOOM 71 (Imperative Care) |
  | Plaque | Angioplasty and stenting | Gateway PTA balloon and Wingspan stent (Stryker) |

* **Auto-sizing.** Sizes follow the published rules:
  * Gateway at ~80% of the normal vessel diameter
  * Wingspan 0.5–1.0 mm over the normal vessel, extending at least 3 mm past each end of the lesion
  * retrievers by their labeled vessel ranges
* **Views.** Switch between **Before**, **Device** (device in position) and **After**.
* **Details panel.** It shows:
  * the procedure inputs
  * IFU/sizing checks (✓ / ! / ✕)
  * device metrics: for aspiration, a Poiseuille flow estimate, the catheter Reynolds number and a tip-force index
  * a before→after table of flow, % of healthy, throat Re, ΔP, Pd/Pa and peak WSS
  * the patent numbers as published on each manufacturer's virtual patent-marking page, each linked to Google Patents
  * source links for every number
* **Patent status.** Each patent list is tagged as listed on the manufacturer's marking page, a
  related patent family only, or not found. The tags are honest about gaps: we found no marked
  patents for Solitaire X or RED 72, and EMBOTRAP is shown with its developer's patent family.
* **What is and isn't simulated.** The model does **not** simulate device–clot or device–wall
  mechanics. Outcomes are inputs you set or published means:
  * pass success and residual stenosis
  * balloon recoil (assumed 15%)
  * post-stent residual stenosis (31.9%, the HDE-study mean) and restenosis

  The simulator computes the hemodynamic consequence of the resulting lumen. A stented segment is
  treated as non-distensible, which shows the compliance mismatch at the stent edges.

`app/model.js` is a line-for-line port of the Python engine below. They agree to about 1e-7
(relative) on flow, Re, ΔP and WSS across mixed test cases, including tapered, bulbed and
perforator-scale vessels. Change both together.

## Python engine (batch runs and scripting)

The same model is available as the `bvm` Python package. Use it for scripted parameter studies and
for the test suite.

## Inputs

**Inputs you can vary**

| Input | Flag | Default | Notes |
|---|---|---|---|
| Luminal diameter | `--diameter-mm` | 3.0 | Unloaded lumen diameter (M1 is roughly 2.5–3.5 mm) |
| Tortuosity | `--tortuosity` | 1.0 | Tortuosity index = arc length / chord length; 1.0 is a straight vessel. Arc length stays fixed, so only curvature changes |
| Occlusion | `--occlusion-pct` | 0 | Cosine-shaped stenosis. `--occlusion-basis diameter` (NASCET-style, default) or `area` |
| Substrate stiffness | `--stiffness-kpa` | 1000 | Young's modulus of the wall/substrate. Use `inf` for a rigid wall |
| Taper | `--taper` | 1.0 | Outlet ÷ inlet unloaded diameter, linear along the arc length |
| Bulb | `--bulb-pct`, `--bulb-center`, `--bulb-length-mm` | 0 | Cosine-shaped local dilation, e.g. a carotid sinus |
| Upstream pressure loss | `--inlet-drop-mmhg` | 0 | Hypotension or proximal disease. Territory BC is calibrated without it, so it lowers perfusion |

**Fixed baseline values**

| Parameter | Value | Status |
|---|---|---|
| Dynamic viscosity μ | 0.0035 Pa·s | Newtonian CFD standard |
| Density ρ | 1060 kg/m³ | Standard |
| Kinematic viscosity ν | 3.30e-6 m²/s | Derived (μ/ρ) |
| Flow Q | 120 mL/min | Initial assumption (`--flow-ml-min`) |
| Mean velocity V | Calculated | V = Q / A(s) at every station |
| Reynolds number Re | Calculated | Re = ρVD/μ at every station |

## Outputs

* **Reynolds number and flow regime.** Re is reported at the inlet, at the stenosis throat and at its
  maximum. Each value is labeled with the straight-pipe criteria: *laminar* below 2000, *transitional*
  from 2000 to 4000, *turbulent* at 4000 and above. For stenoses of 50% area or more, the model also
  flags `post_stenotic_disturbance_possible` once the inlet Re reaches 300. Jets downstream of a
  stenosis can become disturbed well below 2000 (Varghese et al., *J Fluid Mech* 2007). The 300
  cutoff is a conservative screening choice, not a published threshold.
* Mean velocity, pressure, and the distal-to-proximal pressure ratio (Pd/Pa)
* Wall shear stress (WSS), including a low-WSS fraction below 0.4 Pa
* Dean number (curvature), Womersley number, and entrance length
* Wall strain and Moens–Korteweg pulse-wave speed
* With `--bc territory`: territory flow as a fraction of a healthy reference, plus an ischemia band
* Files:
  * `vessel_unloaded_mm.stl` and `vessel_loaded_mm.stl`: watertight STL with `wall`, `inlet` and
    `outlet` solids, in mm. Use them for SimVascular, snappyHexMesh or 3D-printing a mold
  * `vessel.vtk`: open in ParaView; carries per-point WSS, pressure and Re
  * `profiles.csv`: axial profiles along the vessel
  * `summary.json`: all summary results
  * PNG figures

## Usage

```bash
pip install -e '.[test]'
python -m bvm run                                   # baseline straight M1
python -m bvm run --occlusion-pct 70 --tortuosity 1.2 --stiffness-kpa 1000 --bc territory --out results/case1
python -m bvm sweep --param occlusion_pct --values 0:90:10 --bc territory
python -m bvm sweep --param stiffness_kpa --values 100,500,1000,5000,inf
pytest
```

```python
from bvm import Vessel, FlowConditions, solve
sol = solve(Vessel(diameter_mm=2.8, occlusion_pct=60), flow=FlowConditions(bc="territory"))
print(sol.summary()["reynolds"])
```

### Boundary conditions
* `fixed_flow` (default): Q is fixed at 120 mL/min whatever the occlusion. Use it to isolate local
  hemodynamics such as WSS and the throat Re.
* `territory`: the inlet pressure is fixed (85 mmHg) and the vessel drains into a lumped distal
  resistance that stands in for the MCA territory. That resistance is calibrated so the same vessel,
  straight and unoccluded, carries 120 mL/min. Occlusion and tortuosity then *reduce* flow, which is
  the ischemia-relevant case.

## Model

Each station along the centerline is solved as locally fully developed Poiseuille flow, with three
corrections:
* curvature (White 1929 Dean-number friction correlation)
* a Young–Tsai expansion loss (K_t = 1.52) downstream of the stenosis throat
* a thin-wall linear-elastic lumen, ε = P·r/(E·h), coupled to the local pressure. This is how
  substrate stiffness enters the model.

Baseline (straight, 3.0 mm, Q = 120 mL/min):
* rigid wall: V = 0.283 m/s, Re ≈ 257 (laminar)
* E = 1 MPa: the lumen dilates to about 3.17 mm, giving Re ≈ 243

## Limitations (read before citing numbers)
* **Reduced-order, steady, Newtonian.** This is not a 3D Navier–Stokes solution. It does not resolve
  recirculation, jet impingement or post-stenotic turbulence. Use the exported STL in a 3D CFD solver
  (OpenFOAM, SimVascular, FEniCS) to get local flow fields.
* **Entrance effects.** The entrance length (~39 mm at baseline) is longer than the 20 mm segment.
  The flow is still developing, so the inlet WSS is underestimated.
* **Wall law.** The linear thin-wall law is only valid below about 10% strain. Soft hydrogel
  substrates (tens of kPa) exceed this, and the model warns when they do. Those cases need
  hyperelastic FSI.
* **Territory model.** No autoregulation, collaterals or pulsatility are modeled. The ischemia bands
  map flow fraction onto approximate CBF thresholds (penumbra below ~40%, core below ~20% of normal).
  They are heuristic labels, not predictions of infarction.
