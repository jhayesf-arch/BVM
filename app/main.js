import { DEFAULTS, HIGH_WSS_PA, LOW_WSS_PA, MMHG, RE_LAMINAR_MAX, RE_TURBULENT_MIN, bloodFactors, hemoglobin, interp1, oxygenContent, solve, throatRatio, validate } from "./model.js?v=202611151830";
import { Viewer, cssGradient, SLOW_MOTION } from "./viewer.js?v=202611151830";
import { LineChart, fmt } from "./charts.js?v=202611151830";
import { ASPIRATION_DEFAULTS, DEVICES, aspirationEstimate, deviceById, defaultOptions, patentUrl, treat } from "./devices.js?v=202611151830";
import { REFS, SECTIONS, sectionBasis, sectionById, sectionParams } from "./anatomy.js?v=202611151830";
import { PatientViewer, solvePatient } from "./patient.js?v=202611151830";

const $ = (sel) => document.querySelector(sel);

// ───────────────────────── control definitions ─────────────────────────
const FACTORS = [
  { key: "diameterMm", label: "Luminal diameter", unit: "mm", min: 1.5, max: 5, step: 0.05,
    hint: "Inlet lumen diameter of the selected anatomical section. Narrower vessels raise shear and resistance." },
  { key: "tortuosity", label: "Tortuosity index", unit: "", min: 1, max: 1.6, step: 0.01,
    hint: "Arc length ÷ end-to-end distance. 1.00 is straight; the vessel length stays fixed." },
  { key: "occlusionPct", label: "Occlusion", unit: "%", min: 0, max: 95, step: 1,
    hint: "Smooth focal stenosis. Intracranial stenosis of 70% or more by diameter is classed as severe." },
  { key: "stiffnessKpa", label: "Substrate stiffness", unit: "kPa", min: 10, max: 10000, log: true,
    hint: "Young's modulus. Hydrogels ~1–100 kPa · PDMS ~0.5–3 MPa · native artery ~0.1–2 MPa (pressure-dependent)." },
];
const BC_CONTROLS = {
  territory: [
    { key: "flowMlMin", label: "Healthy reference flow", unit: "mL/min", min: 40, max: 250, step: 1 },
    { key: "inletMmHg", label: "Inlet pressure", unit: "mmHg", min: 40, max: 150, step: 1 },
    { key: "outflowMmHg", label: "Venous / ICP back-pressure", unit: "mmHg", min: 0, max: 30, step: 1 },
    { key: "inletDropMmHg", label: "Upstream pressure loss", unit: "mmHg", min: 0, max: 60, step: 1,
      hint: "Systemic hypotension or proximal disease. The territory was calibrated without it, so it lowers perfusion." },
  ],
  fixed_flow: [
    { key: "flowMlMin", label: "Volumetric flow Q", unit: "mL/min", min: 40, max: 250, step: 1 },
    { key: "inletMmHg", label: "Inlet pressure", unit: "mmHg", min: 40, max: 150, step: 1 },
  ],
};
const BC_HINT = {
  territory: "Inlet pressure is fixed and the vessel feeds a lumped downstream resistance, so occlusion reduces flow. This is the ischemia scenario.",
  fixed_flow: "Flow is held constant whatever the geometry. Use it to isolate local shear, Re and pressure loss.",
};
const ADVANCED = [
  { key: "stenosisLengthMm", label: "Stenosis length", unit: "mm", min: 1, max: 30, step: 0.5 },
  { key: "stenosisCenter", label: "Stenosis position", unit: "", min: 0.1, max: 0.9, step: 0.01, display: (v, p) => `${(v * p.lengthMm).toFixed(1)}`, unitFn: () => "mm" },
  { key: "waves", label: "Number of bends", unit: "", min: 0.5, max: 3, step: 0.5 },
  { key: "externalMmHg", label: "External pressure", unit: "mmHg", min: 0, max: 30, step: 1 },
  { key: "mu", label: "Blood viscosity μ at Hct 45%", unit: "Pa·s", min: 0.002, max: 0.006, step: 0.0001 },
  { key: "microDiamUm", label: "Microvascular diameter for Hct scaling", unit: "µm", min: 8, max: 100, step: 1 },
  { key: "cbfNormal", label: "Normal territory CBF", unit: "mL/100g/min", min: 30, max: 80, step: 1 },
  { key: "oefMax", label: "Maximum O₂ extraction fraction", unit: "", min: 0.6, max: 0.95, step: 0.01 },
  { key: "rho", label: "Blood density ρ", unit: "kg/m³", min: 1000, max: 1100, step: 1 },
];

// Blood composition and oxygen (bvm/oxygen.py)
const BLOOD = [
  { key: "hct", label: "Hematocrit", unit: "%", min: 15, max: 65, step: 1,
    hint: "Sets O₂ content (Hb = Hct × MCHC 33 g/dL) and, if coupled, the viscosity of the vessel and the microvascular bed. 45% is the reference." },
  { key: "sao2", label: "Arterial O₂ saturation", unit: "%", min: 70, max: 100, step: 1 },
  { key: "pao2", label: "Arterial PO₂", unit: "mmHg", min: 40, max: 150, step: 1 },
  { key: "cmro2", label: "O₂ demand (CMRO₂)", unit: "mL/100g/min", min: 1.5, max: 5, step: 0.1,
    hint: "Normal resting brain ≈ 3.3 (Kety & Schmidt 1948)." },
];

// Downstream tissue (territory BC only; bvm FlowConditions.collateral_pct / no_reflow_pct)
const TISSUE = [
  { key: "collateralPct", label: "Collateral supply", unit: "% of healthy flow", min: 0, max: 60, step: 1, territoryOnly: true,
    hint: "Leptomeningeal collaterals: the territory flow they would deliver with the vessel fully occluded. 0 = none (worst case). Roughly: poor < 10%, moderate 10–30%, good > 30% (assumed mapping of collateral grades)." },
  { key: "noReflowPct", label: "No-reflow (microvascular obstruction)", unit: "% of bed", min: 0, max: 90, step: 1, territoryOnly: true,
    hint: "Share of the distal capillary bed obstructed (microthrombi, swollen endothelium). Distal resistance ÷ (1 − share). Limits flow even after the artery is reopened; perforator beds are not affected." },
];

// Pulsatility (bvm/pulsatile.py); none of these change the mean results
const PULSE = [
  { key: "heartRate", label: "Heart rate", unit: "bpm", min: 40, max: 150, step: 1,
    hint: "Sets the harmonics of the cardiac cycle and the Womersley number." },
  { key: "pulsatilityIndex", label: "Healthy pulsatility index (PI)", unit: "", min: 0, max: 2, step: 0.01,
    hint: "Gosling PI = (Qmax − Qmin)/Qmean of the healthy vessel's inflow; fixes the inlet pressure pulse. TCD normal for the MCA ≈ 0.6–1.1. 0 = steady flow." },
  { key: "tauS", label: "Territory time constant τ = R·C", unit: "s", min: 0.05, max: 0.6, step: 0.01,
    hint: "Compliance of the downstream bed. ≈ 0.2 s in health (Kasprowicz et al.); longer τ damps pulsation more." },
];

// Anatomy-tab shape controls (healthy geometry of the section)
const SHAPE = [
  { key: "lengthMm", label: "Segment length", unit: "mm", min: 10, max: 40, step: 1 },
  { key: "taper", label: "Taper (outlet ÷ inlet diameter)", unit: "", min: 0.4, max: 1.2, step: 0.01,
    hint: "1.00 is a uniform tube. Below 1 the vessel narrows downstream." },
  { key: "bulbPct", label: "Bulb dilation", unit: "%", min: 0, max: 80, step: 1,
    hint: "Peak local widening, e.g. the carotid sinus. 0 = none." },
  { key: "bulbCenter", label: "Bulb position", unit: "", min: 0.1, max: 0.9, step: 0.01, display: (v, p) => `${(v * p.lengthMm).toFixed(1)}`, unitFn: () => "mm" },
  { key: "bulbLengthMm", label: "Bulb length", unit: "mm", min: 2, max: 30, step: 0.5 },
  { key: "wallMm", label: "Wall thickness", unit: "mm", min: 0.02, max: 1, step: 0.01 },
  { key: "nBranches", label: "Perforator branches", unit: "", min: 0, max: 12, step: 1,
    hint: "Side branches (e.g. lenticulostriates off M1). Each drains through its own resistance calibrated on the healthy vessel." },
  { key: "branchFlowPct", label: "Perforators' share of healthy flow", unit: "%", min: 0, max: 10, step: 0.1 },
  { key: "branchStart", label: "First branch at", unit: "", min: 0.05, max: 0.95, step: 0.01, display: (v, p) => `${(v * p.lengthMm).toFixed(1)}`, unitFn: () => "mm" },
  { key: "branchEnd", label: "Last branch at", unit: "", min: 0.05, max: 0.95, step: 0.01, display: (v, p) => `${(v * p.lengthMm).toFixed(1)}`, unitFn: () => "mm" },
];

const PRESETS = [
  { name: "Healthy M1 (baseline)", p: {} },
  { name: "Acute M1 clot · near-occlusion", p: { lesion: "thrombus", occlusionPct: 95, stenosisLengthMm: 8 } },
  { name: "Symptomatic ICAS · 75% (HDE mean)", p: { lesion: "icas", occlusionPct: 75 } },
  { name: "Moderate stenosis · 50%", p: { occlusionPct: 50 } },
  { name: "Severe stenosis · 70%", p: { occlusionPct: 70 } },
  { name: "Critical stenosis · 85%", p: { occlusionPct: 85 } },
  { name: "Tortuous vessel", p: { tortuosity: 1.3, waves: 2 } },
  { name: "Tortuous + 60% stenosis", p: { tortuosity: 1.3, waves: 2, occlusionPct: 60 } },
  { name: "Narrow vessel · 2.2 mm", p: { diameterMm: 2.2 } },
  { name: "Stiffened wall · 5 MPa", p: { stiffnessKpa: 5000 } },
  { name: "Soft hydrogel BVM · 50 kPa", p: { stiffnessKpa: 50 } },
  { name: "Rigid wall (classic CFD)", p: { rigid: true } },
];

// `floor` is the minimum top of the color scale, so a healthy vessel reads as mid-range
// rather than saturating the map.
const FIELDS = {
  wss: { label: "WSS (Pa)", cmap: "inferno", zero: true, floor: 5 },
  pressure: { label: "Pressure (mmHg)", cmap: "coolwarm", zero: false },
  velocity: { label: "Mean velocity (m/s)", cmap: "viridis", zero: true, floor: 0.5 },
  reynolds: { label: "Reynolds number", cmap: "viridis", zero: true, floor: 500 },
  strain: { label: "Wall strain (%)", cmap: "viridis", zero: true, floor: 10 },
  wallo2: { label: "Wall O₂ transfer index", cmap: "viridis", zero: true, floor: 1.2 },
};

const SENS_X = [...FACTORS, BC_CONTROLS.territory[3], ...TISSUE, PULSE[1], SHAPE[1], BLOOD[0], BLOOD[1]].map((f) => ({ key: f.key, label: `${f.label}${f.unit ? ` (${f.unit})` : ""}` }));
const sensDef = (key) => [...FACTORS, ...BC_CONTROLS.territory, ...TISSUE, ...PULSE, ...SHAPE, ...BLOOD].find((f) => f.key === key);
const SENS_Y = [
  { key: "flow_pct", label: "Flow, % of healthy", get: (s) => (s.flow_fraction ?? 1) * 100, bands: true },
  { key: "Q_mL_min", label: "Flow (mL/min)", get: (s) => s.Q_mL_min },
  { key: "Re_throat", label: "Reynolds number at throat", get: (s) => s.Re_throat, re: true },
  { key: "Re_inlet", label: "Reynolds number at inlet", get: (s) => s.Re_inlet, re: true },
  { key: "dP_mmHg", label: "Pressure drop (mmHg)", get: (s) => s.dP_mmHg },
  { key: "Pd_Pa", label: "Pd / Pa", get: (s) => s.Pd_Pa },
  { key: "wss_max", label: "Peak WSS (Pa)", get: (s) => s.wss_max },
  { key: "V_throat", label: "Throat velocity (m/s)", get: (s) => s.V_throat },
  { key: "post_pct", label: "Flow after device, % of healthy", post: true, bands: true, get: (s, b) => (b.flow_fraction ?? 1) * 100 },
  { key: "restored", label: "Flow restored, % of deficit", post: true, get: (s, b) => restoredPct(s, b) ?? NaN },
  { key: "dP_lesion", label: "Pressure drop across lesion (mmHg)", get: (s) => s.dP_lesion_mmHg ?? 0 },
  { key: "branch_pct", label: "Perforator flow, % of healthy", bands: true, get: (s) => (s.branch_fraction ?? NaN) * 100 },
  { key: "do2", label: "O₂ delivery to territory (mL O₂/min)", get: (s) => s.o2_distal.DO2_mL_min },
  { key: "cbf", label: "Territory CBF (mL/100 g/min)", get: (s) => s.o2_distal.CBF_mL_100g_min ?? NaN },
  { key: "oef", label: "Required O₂ extraction (territory)", oef: true, get: (s) => s.o2_distal.required_OEF ?? NaN },
  { key: "oef_perf", label: "Required O₂ extraction (perforator territory)", oef: true, get: (s) => s.o2_perforators?.required_OEF ?? NaN },
  { key: "pi_vessel", label: "Pulsatility index, vessel", get: (s) => s.pulse.PI_inlet },
  { key: "pi_terr", label: "Pulsatility index, territory", get: (s) => s.pulse.PI_territory ?? NaN },
  { key: "wss_sys", label: "Peak systolic WSS (Pa)", get: (s) => s.pulse.wss_max_systolic_Pa ?? NaN },
];

// share of the flow deficit recovered by the device; null when there is no deficit to recover
function restoredPct(a, b) {
  if (a.flow_fraction == null || b?.flow_fraction == null) return null;
  const deficit = 1 - a.flow_fraction;
  return deficit < 0.01 ? null : ((b.flow_fraction - a.flow_fraction) / deficit) * 100;
}

// ───────────────────────── state ─────────────────────────
let params = { ...DEFAULTS, ...readHash() };
let result = null;
const hashExtras = new URLSearchParams(location.hash.slice(1));
let deviceId = deviceById(hashExtras.get("device") ?? "") ? hashExtras.get("device") : null;
let deviceOpts = defaultOptions(deviceById(deviceId));
let view = ["before", "inplace", "after"].includes(hashExtras.get("view")) ? hashExtras.get("view") : "after";
let panelOpen = false;
let treatment = null, pre = null, post = null;
let lockedRange = null;
let probeFrac = 0.5;
const controls = new Map(); // key -> {sync()}

// ───────────────────────── control rendering ─────────────────────────
const toSlider = (f, v) => (f.log ? (Math.log10(v / f.min) / Math.log10(f.max / f.min)) * 1000 : v);
const fromSlider = (f, t) => (f.log ? f.min * (f.max / f.min) ** (t / 1000) : +t);
const decimals = (step) => Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
const show = (f, v) => (f.display ? f.display(v, params) : f.log ? (v >= 1000 ? (+v.toPrecision(3)).toString() : v.toFixed(0)) : v.toFixed(decimals(f.step)));

function makeControl(f, container, { primary = false, tag = "", hint = true } = {}) {
  const wrap = document.createElement("div");
  wrap.className = `factor ${primary ? "primary" : "compact"}`;
  const id = `in-${tag}${f.key}`;
  wrap.innerHTML = `
    <div class="factor-head">
      <label for="${id}"${f.hint && !hint ? ` title="${f.hint}"` : ""}>${f.label}</label>
      <span class="num"><input type="number" id="${id}-num" aria-label="${f.label} value" /><span class="unit">${f.unitFn ? f.unitFn() : f.unit}</span></span>
    </div>
    <input type="range" id="${id}" min="${f.log ? 0 : f.min}" max="${f.log ? 1000 : f.max}" step="${f.log ? 1 : f.step}" />
    ${f.hint && hint ? `<p class="hint">${f.hint}</p>` : ""}`;
  container.appendChild(wrap);
  const range = wrap.querySelector(`#${CSS.escape(id)}`), num = wrap.querySelector(`#${CSS.escape(id)}-num`);
  num.step = f.display ? "any" : f.log ? "any" : f.step;
  const sync = () => {
    range.value = toSlider(f, params[f.key]);
    num.value = show(f, params[f.key]);
    // patient mode always uses the territory BC, so its tissue sliders stay live
    const disabled = (f.key === "stiffnessKpa" && params.rigid) || (f.territoryOnly && params.bc === "fixed_flow" && !tag);
    range.disabled = num.disabled = disabled;
    wrap.style.opacity = disabled ? 0.45 : 1;
  };
  range.addEventListener("input", () => set(f.key, fromSlider(f, +range.value)));
  num.addEventListener("change", () => {
    let v = parseFloat(num.value);
    if (!isFinite(v)) return sync();
    if (["stenosisCenter", "bulbCenter", "branchStart", "branchEnd"].includes(f.key)) v = v / params.lengthMm;
    set(f.key, Math.min(Math.max(v, f.min), f.max));
  });
  controls.set(tag + f.key, { sync, wrap });
  return wrap;
}

function makeSegmented(container, options, key) {
  const seg = document.createElement("div");
  seg.className = "segmented small";
  seg.setAttribute("role", "radiogroup");
  seg.innerHTML = options.map(([v, l]) => `<button role="radio" data-value="${v}">${l}</button>`).join("");
  seg.addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) set(key, b.dataset.value); });
  container.appendChild(seg);
  const sync = () => seg.querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === params[key]));
  controls.set(`${key}-seg`, { sync });
}

function buildControls() {
  const host = $("#factors");
  for (const f of FACTORS) {
    const w = makeControl(f, host, { primary: true, hint: false });
    if (f.key === "occlusionPct") {
      const opts = document.createElement("div");
      opts.className = "inline-opts";
      w.insertBefore(opts, w.querySelector(".hint"));
      makeSegmented(opts, [["diameter", "% diameter"], ["area", "% area"]], "basis");
      opts.firstChild.style.flex = "1";
    }
    if (f.key === "stiffnessKpa") {
      const lab = document.createElement("label");
      lab.className = "check";
      lab.innerHTML = `<input type="checkbox" id="rigid" /> Rigid wall (ignore stiffness)`;
      w.appendChild(lab);
      const cb = lab.querySelector("input");
      cb.addEventListener("change", () => set("rigid", cb.checked));
      controls.set("rigid", { sync: () => (cb.checked = params.rigid) });
    }
  }
  const adv = $("#advanced");
  for (const f of ADVANCED) makeControl(f, adv);
  for (const f of SHAPE) makeControl(f, $("#shape"));
  for (const f of BLOOD) makeControl(f, $("#blood"), { hint: false });
  for (const f of TISSUE) makeControl(f, $("#tissue"), { hint: false });
  if ($("#pulse-controls")) for (const f of PULSE) makeControl(f, $("#pulse-controls"), { hint: false });
  if ($("#pt-tissue")) for (const f of TISSUE) makeControl(f, $("#pt-tissue"), { tag: "pt-", hint: false });
  const nnBoxes = ["#non-newtonian", "#pt-non-newtonian"].map($).filter(Boolean);
  for (const nn of nnBoxes) if (!nn.dataset.wired) { nn.addEventListener("change", () => set("nonNewtonian", nn.checked)); nn.dataset.wired = "1"; }
  controls.set("nonNewtonian", { sync: () => nnBoxes.forEach((nn) => (nn.checked = params.nonNewtonian)) });
  for (const f of BLOOD) makeControl(f, $("#pt-blood"), { tag: "pt-" });
  const phc = $("#pt-hct-couple");
  if (phc && !phc.dataset.wired) { phc.addEventListener("change", () => set("hctCouple", phc.checked)); phc.dataset.wired = "1"; }
  const hc = $("#hct-couple");
  if (!hc.dataset.wired) { hc.addEventListener("change", () => set("hctCouple", hc.checked)); hc.dataset.wired = "1"; }
  controls.set("hctCouple", { sync: () => {
    hc.checked = params.hctCouple;
    if (phc) phc.checked = params.hctCouple;
    const { muEff, distalScale } = bloodFactors(params);
    const derived = `Hb ${hemoglobin(params).toFixed(1)} g/dL · CaO₂ ${oxygenContent(params).toFixed(1)} mL O₂/dL · ` +
      `μ ${(muEff * 1e3).toFixed(2)} mPa·s · microvascular resistance ×${distalScale.toFixed(2)}`;
    $("#blood-derived").textContent = derived;
    const pd = $("#pt-blood-derived");
    if (pd) pd.textContent = derived.replace(/μ [\d.]+ mPa·s/, `μ in this M1 ${(bloodFactors({ ...params, diameterMm: ptMetas?.healthy?.profile.D_healthy_mm[0] ?? 2.79 }).muEff * 1e3).toFixed(2)} mPa·s`);
  } });
  renderBcControls();
}

function renderBcControls() {
  const host = $("#bc-controls");
  for (const key of ["flowMlMin", "inletMmHg", "outflowMmHg", "inletDropMmHg"]) controls.delete(key);
  host.innerHTML = "";
  for (const f of BC_CONTROLS[params.bc]) makeControl(f, host, { hint: false });
  $("#bc").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === params.bc));
  $("#bc-hint").textContent = BC_HINT[params.bc];
  syncControls();
}

function syncControls() { for (const c of controls.values()) c.sync(); }

function clampParams(q) {
  // keep the stenosis and bulb inside the vessel, and leave some driving pressure
  for (const [len, ctr] of [["stenosisLengthMm", "stenosisCenter"], ["bulbLengthMm", "bulbCenter"]]) {
    if (q[len] / 2 / q.lengthMm >= 0.5) q[len] = q.lengthMm * 0.9;
    const h = q[len] / 2 / q.lengthMm;
    q[ctr] = Math.min(Math.max(q[ctr], h), 1 - h);
  }
  if (q.basis === "diameter" && q.occlusionPct > 95) q.occlusionPct = 95;
  q.inletDropMmHg = Math.min(q.inletDropMmHg, q.inletMmHg - q.outflowMmHg - 2);
  if (q.branchStart > q.branchEnd) q.branchEnd = q.branchStart;
  q.nBranches = Math.round(q.nBranches);
  return q;
}

// Slider ranges follow the selected anatomical section (a perforator is 0.1–1.3 mm, a carotid 4–9 mm).
function applySectionRanges(sec) {
  const setR = (f, [min, max, step]) => Object.assign(f, { min, max, step });
  setR(FACTORS[0], sec.ranges.diameterMm);
  for (const list of Object.values(BC_CONTROLS)) setR(list[0], sec.ranges.flowMlMin);
  setR(SHAPE[0], sec.ranges.lengthMm);
}

function rebuildControls() {
  controls.clear();
  for (const id of ["#factors", "#advanced", "#shape", "#blood", "#pt-blood"]) $(id).innerHTML = "";
  buildControls();
}

function applySection(secId, srcId, extra = {}) {
  const sec = sectionById(secId);
  srcId = sec.sources.some((s) => s.id === srcId) ? srcId : sec.sources[0].id;
  params = clampParams({ ...DEFAULTS, bc: params.bc, lesion: params.lesion, ...sectionParams(sec, srcId), source: srcId, ...extra });
  const dev = deviceById(deviceId);
  if (dev && dev.lesion !== params.lesion) { deviceId = null; panelOpen = false; }
  applySectionRanges(sec);
  rebuildControls();
  lockedRange = null;
  schedule();
}

function set(key, value) {
  params = clampParams({ ...params, [key]: value });
  if (key === "bc") renderBcControls();
  $("#preset").value = "";
  syncControls();
  schedule();
  if (geomMode === "patient") renderPatient();
}

// ───────────────────────── viewer, charts ─────────────────────────
const viewer = new Viewer($("#viewport"));
const onProbe = (f) => { probeFrac = f; $("#probe").value = f; updateProbe(); };
const charts = {
  diam: new LineChart($("#c-diam"), { title: "Lumen diameter", unit: "mm", onProbe, zeroBased: true }),
  pres: new LineChart($("#c-pres"), { title: "Pressure", unit: "mmHg", onProbe }),
  re: new LineChart($("#c-re"), { title: "Reynolds number", unit: "", onProbe, zeroBased: true,
    thresholds: [{ y: RE_LAMINAR_MAX, label: "laminar limit 2000", color: "#f5b041" }] }),
  wss: new LineChart($("#c-wss"), { title: "Wall shear stress", unit: "Pa", onProbe, zeroBased: true,
    thresholds: [{ y: LOW_WSS_PA, label: "low-WSS 0.4 Pa", color: "#8a97a8" }] }),
};
const sensChart = new LineChart($("#c-sens"), { title: "", zeroBased: true });

// ───────────────────────── compute loop ─────────────────────────
let pending = false, firstFrame = true, sensTimer = null;
function schedule() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => { pending = false; compute(); });
  clearTimeout(sensTimer);
  sensTimer = setTimeout(computeSensitivity, 220);
}

function compute() {
  const err = $("#vp-error");
  const dev = deviceById(deviceId);
  try {
    pre = solve(params);
    treatment = dev ? treat(dev, params, deviceOpts) : null;
    post = treatment ? solve(treatment.post) : null;
    result = post && view === "after" ? post : pre;
    err.hidden = true;
  } catch (e) {
    err.textContent = e.message;
    err.hidden = false;
    return;
  }
  const field = $("#field").value, meta = FIELDS[field];
  const values = viewer.fieldValues(result, field);
  let range = fieldRange(values, meta);
  if ($("#lock-scale").checked) { lockedRange ??= { field, range }; if (lockedRange.field === field) range = lockedRange.range; else lockedRange = { field, range }; }
  viewer.update(result, field, meta.cmap, range, { reframe: firstFrame || geometryChanged() });
  firstFrame = false;
  viewer.setExtras(sceneExtras(dev));
  viewer.setProbe(probeFrac);
  renderDeviceBar(dev);
  renderDevicePanel(dev);
  renderAnatomy();
  updateRestoration(dev);
  renderEndpoints(dev);
  renderOxygen();
  renderPulse();
  syncSensOptions(dev);
  updateLegend(meta, range);
  const bl = $("#branch-legend");
  bl.hidden = !(params.nBranches > 0) || geomMode === "patient";   // the patient geometry has no perforators
  bl.title = "Side-branch stubs: perforating arteries leaving this vessel (lenticulostriates off M1, pontine perforators off the basilar). Color = branch flow as a share of its healthy flow.";
  $("#branch-legend-title").textContent = { m1: "Lenticulostriates", basilar: "Pontine perforators" }[params.section] ?? "Perforators";
  updateCharts();
  updateResults();
  updateProbe();
  writeHash();
}

let lastGeomKey = "";
function geometryChanged() {
  // reframe the camera only when the overall shape changes a lot (length/tortuosity/bends)
  const key = `${params.section}|${params.lengthMm}|${params.tortuosity > 1.0001}|${params.waves}`;
  const changed = key !== lastGeomKey;
  lastGeomKey = key;
  return changed;
}

function fieldRange(values, meta) {
  let lo = Math.min(...values), hi = Math.max(...values);
  if (meta.zero) lo = Math.min(0, lo);
  if (meta.floor) hi = Math.max(hi, meta.floor);
  if (hi - lo < 1e-9 * Math.max(1, Math.abs(hi))) { const m = Math.abs(hi) * 0.05 || 1; lo -= m; hi += m; }
  return [lo, hi];
}

function updateLegend(meta, [lo, hi]) {
  $("#legend-title").textContent = meta.label;
  $("#legend-bar").style.background = cssGradient(meta.cmap);
  $("#legend-ticks").innerHTML = [1, 0.75, 0.5, 0.25, 0].map((t) => `<span>${fmt(lo + t * (hi - lo))}</span>`).join("");
}

function updateCharts() {
  if (geomMode === "patient") return;
  const s = Array.from(result.geom.s, (v) => v * 1e3);
  const xr = [0, params.lengthMm];
  charts.diam.set([
    { x: s, y: Array.from(result.geom.radius, (r) => 2e3 * r), color: "#8a97a8", dashed: true, label: "unloaded" },
    { x: s, y: Array.from(result.radius, (r) => 2e3 * r), label: "pressurized" },
  ], { xLabel: "mm", probe: probeFrac, xRange: xr });
  charts.pres.set([{ x: s, y: Array.from(result.pressure, (p) => p / MMHG) }], { xLabel: "mm", probe: probeFrac, xRange: xr });
  charts.re.set([{ x: s, y: Array.from(result.reynolds) }], { xLabel: "mm", probe: probeFrac, xRange: xr });
  charts.wss.set([{ x: s, y: Array.from(result.wss) }], { xLabel: "mm", probe: probeFrac, xRange: xr });
}

function badge(el, cls, text) { el.className = `badge ${cls}`; el.textContent = text ?? cls.replace("-", " "); }
const logPct = (re) => Math.min(Math.max((Math.log10(Math.max(re, 10)) - 1) / 3, 0), 1) * 100;

function updateResults() {
  const S = result.summary, p = params, hasSten = p.occlusionPct > 0;
  $("#re-in").textContent = S.Re_inlet.toFixed(0);
  badge($("#re-in-badge"), S.regime_inlet);
  $("#re-throat").textContent = hasSten ? S.Re_throat.toFixed(0) : "—";
  if (hasSten) badge($("#re-throat-badge"), S.regime_throat); else badge($("#re-throat-badge"), "", "no stenosis");
  $("#re-mark-in").style.left = `${logPct(S.Re_inlet)}%`;
  $("#re-mark-thr").style.left = `${logPct(S.Re_throat)}%`;
  $("#re-mark-thr").style.display = hasSten ? "" : "none";
  const worst = S.regime_max;
  const pill = $("#regime-pill");
  pill.className = `pill ${worst}`;
  pill.textContent = `${worst[0].toUpperCase()}${worst.slice(1)} flow · Re max ${S.Re_max.toFixed(0)}`;
  const note = $("#re-note");
  if (S.disturbance) {
    note.className = "note alert";
    note.textContent = `Laminar by the pipe criterion, but the jet leaving a ${S.area_stenosis_pct.toFixed(0)}%-area stenosis at inlet Re ${S.Re_inlet.toFixed(0)} may break down into disturbed, transitional flow downstream. Confirm with 3D CFD.`;
  } else if (S.regime_max !== "laminar") {
    note.className = "note alert";
    note.textContent = `Re exceeds ${S.regime_max === "turbulent" ? RE_TURBULENT_MIN : RE_LAMINAR_MAX} in part of the vessel; this model's laminar losses underestimate pressure drop there.`;
  } else {
    note.className = "note";
    note.textContent = `Laminar everywhere (Re < 2000). Inertia dominates viscosity ${S.Re_max >= 100 ? "moderately" : "only weakly"}; the Dean number is ${S.Dean_max.toFixed(0)}${S.Dean_max > 11.6 ? ", so bends add secondary swirl and extra loss" : ""}.`;
  }

  const rows = [
    ["Mean velocity, inlet", `${fmt(S.V_inlet)} m/s`],
    ["Mean velocity, throat", hasSten ? `${fmt(S.V_throat)} m/s` : "—"],
    ["Pressurized diameter", `${fmt(S.D_inlet_mm)} mm`],
    ["Throat diameter", hasSten ? `${fmt(S.D_throat_mm)} mm` : "—"],
    ["Area stenosis", hasSten ? `${S.area_stenosis_pct.toFixed(0)}%` : "—"],
    ["Pressure drop ΔP", `${fmt(S.dP_mmHg)} mmHg`],
    ["Pd / Pa", S.Pd_Pa.toFixed(3), S.Pd_Pa < 0.8],
    ["WSS, inlet", `${fmt(S.wss_inlet)} Pa`],
    ["WSS, peak", `${fmt(S.wss_max)} Pa`, S.high_wss],
    ["Low-WSS (<0.4 Pa) length", `${(S.low_wss_fraction * 100).toFixed(0)}%`, S.low_wss_fraction > 0],
    ["Dean number, max", fmt(S.Dean_max)],
    ["Womersley number", fmt(S.Womersley)],
    ["Wall strain, max", `${(S.max_strain * 100).toFixed(1)}%`, S.max_strain > 0.1],
    ["Pulse-wave speed", S.wave_speed == null ? "rigid" : `${fmt(S.wave_speed)} m/s`],
  ];
  $("#metrics").innerHTML = rows.map(([k, v, flag]) => `<dt>${k}</dt><dd class="${flag ? "flag" : ""}">${v}</dd>`).join("");

  const card2 = $("#warnings-card");
  card2.hidden = !S.warnings.length;
  $("#warnings").innerHTML = S.warnings.map((w) => `<li>${w}</li>`).join("");
  $("#motion-pill").textContent = `· particles ${SLOW_MOTION}× slow-motion`;
}

const qFmt = (q) => (q < 1 ? q.toFixed(3) : q < 10 ? q.toFixed(2) : q.toFixed(1));
const BAND_COLOR = { normal: "var(--ok)", oligemia: "var(--warn)", "penumbra-range": "var(--orange)", "core-range": "var(--bad)" };

function updateRestoration(dev) {
  const card = $("#restore-card"), a = pre.summary, b = post?.summary;
  $("#flow-before").innerHTML = `${qFmt(a.Q_territory_mL_min)}<small> mL/min</small>`;
  $("#after-dev").textContent = dev ? `· ${dev.name}` : "";
  const meterA = $("#meter-after"), meterB = $("#meter-before");
  if (a.flow_fraction == null) {
    card.classList.add("disabled");
    $("#frac-before").textContent = "flow is prescribed";
    $("#flow-after").innerHTML = b ? `${qFmt(b.Q_territory_mL_min)}<small> mL/min</small>` : "—";
    $("#frac-after").textContent = b ? "flow is prescribed" : "no device";
    badge($("#restored"), "", "switch to Brain territory");
    badge($("#band"), "", "—");
    meterA.style.width = "0"; meterB.style.width = "0";
    $("#restore-note").textContent = "Restoration needs the Brain territory boundary condition. With fixed flow the device cannot change Q, only ΔP, shear and Re.";
    return;
  }
  card.classList.remove("disabled");
  const fa = a.flow_fraction * 100, refQ = a.reference_distal_mL_min;
  $("#frac-before").textContent = `${fa.toFixed(0)}% of healthy`;
  meterB.style.width = `${Math.min(fa, 100)}%`;
  if (!b) {
    $("#flow-after").innerHTML = "—";
    $("#frac-after").textContent = "pick a device";
    Object.assign(meterA.style, { width: `${Math.min(fa, 100)}%`, background: BAND_COLOR[a.ischemia_band] });
    const deficit = refQ - a.Q_territory_mL_min;
    badge($("#restored"), "", deficit > 0.01 * refQ ? `deficit ${qFmt(deficit)} mL/min` : "no flow deficit");
    badge($("#band"), a.ischemia_band);
    $("#restore-note").textContent = deficit > 0.01 * refQ
      ? "Choose a device in the viewport bar (Plaque → angioplasty/stent, Clot → thrombectomy) to model how much of this deficit it recovers."
      : "The vessel carries its healthy reference flow. Add an occlusion or upstream pressure loss to create a deficit.";
    return;
  }
  const fb = b.flow_fraction * 100, r = restoredPct(a, b);
  $("#flow-after").innerHTML = `${qFmt(b.Q_territory_mL_min)}<small> mL/min</small>`;
  $("#frac-after").textContent = `${fb.toFixed(0)}% of healthy`;
  Object.assign(meterA.style, { width: `${Math.min(fb, 100)}%`, background: BAND_COLOR[b.ischemia_band] });
  if (r == null) badge($("#restored"), "", "no flow deficit to restore");
  else if (r > 100.5) badge($("#restored"), "restored-good", "fully restored · flow above baseline");
  else badge($("#restored"), r >= 90 ? "restored-good" : r >= 50 ? "restored-part" : "restored-poor", `restored ${Math.max(r, 0).toFixed(0)}% of deficit`);
  badge($("#band"), b.ischemia_band, `after: ${b.ischemia_band.replace("-", " ")}`);
  const gain = b.Q_territory_mL_min - a.Q_territory_mL_min;
  const failed = !dev.balloon && !dev.stent && !deviceOpts.success;
  const causes = [post.params.occlusionPct > 0 && "the residual narrowing after treatment", params.inletDropMmHg > 0 && "the upstream pressure loss", params.noReflowPct > 0 && `no-reflow in ${params.noReflowPct}% of the microvascular bed`].filter(Boolean);
  const remaining = fb < 99 && !failed && causes.length ? ` The remaining ${(100 - fb).toFixed(0)}% gap comes from ${causes.length > 1 ? `${causes.slice(0, -1).join(", ")} and ${causes.at(-1)}` : causes[0]}.` : "";
  const qz = (q) => (q < 0.0005 ? "0" : qFmt(q));
  const collNote = params.collateralPct > 0 && params.bc === "territory" ? ` Collaterals supply ${qz(a.Q_collateral_mL_min)} mL/min before and ${qz(b.Q_collateral_mL_min)} after: they carry flow only while the vessel's outlet pressure is below normal.` : "";
  const hctNote = params.hct !== 45 ? ` "Healthy" is a healthy vessel at normal hematocrit (45%); at ${params.hct}% the blood's viscosity ${params.hct < 45 ? "is lower, so flow can exceed 100%" : "is higher, so flow stays below 100%"} even with a fully open lumen. See Oxygen delivery for the net effect.` : "";
  $("#restore-note").textContent = failed
    ? "Pass set as unsuccessful: the clot stays, so flow is unchanged."
    : `${gain >= 0 ? "+" : ""}${qFmt(gain)} mL/min after treatment.${remaining}${collNote} The bar's white line marks flow before treatment.${hctNote}`;
}

// ───────────────────────── pulsatility ─────────────────────────
const pulseChart = $("#c-pulse") ? new LineChart($("#c-pulse"), { title: "Flow over one cycle", unit: "mL/min", zeroBased: true }) : null;
function renderPulse() {
  if (!pulseChart) return;
  const a = pre.summary.pulse, b = post?.summary.pulse ?? null, terr = a.waveform_Q_territory_mL_min != null;
  const t = [...a.waveform_t, 1].map((x) => x * 60 / params.heartRate);
  const loop = (y) => [...y, y[0]];
  const series = [{ x: t, y: loop(a.waveform_Q_inlet_mL_min), color: "#8a97a8" }];
  if (terr) series.push({ x: t, y: loop(a.waveform_Q_territory_mL_min), color: "#8a97a8", dashed: true });
  if (b) {
    series.push({ x: t, y: loop(b.waveform_Q_inlet_mL_min), color: "#4cc2ff" });
    if (terr) series.push({ x: t, y: loop(b.waveform_Q_territory_mL_min), color: "#4cc2ff", dashed: true });
  }
  pulseChart.set(series, { xLabel: "s", xRange: [0, t[t.length - 1]] });
  const sw = (color, dashed) => `<i style="border-top:2px ${dashed ? "dashed" : "solid"} ${color}"></i>`;
  $("#pulse-legend").innerHTML = `<span>${sw("#8a97a8")}${sw("#4cc2ff")}into the vessel</span>` +
    (terr ? `<span>${sw("#8a97a8", 1)}${sw("#4cc2ff", 1)}to the territory</span>` : "") +
    `<span class="muted">grey before${b ? ", blue after" : ""}</span>`;
  const f = (v, d = 1) => (v == null || !isFinite(v) ? "—" : (+v.toFixed(d)).toLocaleString());
  const row = (k, x, y) => `<tr><td>${k}</td><td>${x}</td><td class="${b ? "" : "na"}">${b ? y : "—"}</td></tr>`;
  const range = (s, key) => `${f(s[`Q_${key}_systolic_mL_min`])} / ${f(s[`Q_${key}_diastolic_mL_min`])}`;
  let html = `<thead><tr><th></th><th>Before</th><th>After</th></tr></thead>` +
    row("Vessel inflow, systolic / diastolic, mL/min", range(a, "inlet"), b && range(b, "inlet")) +
    row("Pulsatility index, vessel", f(a.PI_inlet, 2), b && f(b.PI_inlet, 2));
  if (terr) html +=
    row("Territory flow, systolic / diastolic, mL/min", range(a, "territory"), b && range(b, "territory")) +
    row("Pulsatility index, territory", f(a.PI_territory, 2), b && f(b.PI_territory, 2));
  html += row("Peak systolic WSS, Pa", f(a.wss_max_systolic_Pa), b && f(b.wss_max_systolic_Pa)) +
    row("Lowest WSS at end diastole, Pa", f(a.wss_min_diastolic_Pa, 2), b && f(b.wss_min_diastolic_Pa, 2)) +
    row("Peak systolic Re", f(a.Re_max_systolic, 0), b && f(b.Re_max_systolic, 0)) +
    row("Systolic ΔP along vessel, mmHg", f(a.dP_vessel_systolic_mmHg, 1), b && f(b.dP_vessel_systolic_mmHg, 1));
  $("#pulse").innerHTML = html;
  const damp = terr && a.PI_inlet < 0.8 * params.pulsatilityIndex ? ` The lesion damps the vessel's pulsation (PI ${f(a.PI_inlet, 2)} vs ${f(params.pulsatilityIndex, 2)} healthy): the low-PI, slow-rising "tardus–parvus" pattern TCD shows distal to a tight stenosis.` : "";
  const rev = a.flow_reversal || b?.flow_reversal ? " Flow reverses in diastole, outside this model's quasi-steady range." : "";
  $("#pulse-note").textContent = (terr
    ? `Inlet pulse pressure ${f(a.inlet_pulse_pressure_mmHg)} mmHg, calibrated so the healthy vessel has PI ${f(params.pulsatilityIndex, 2)}. Peak systolic WSS is ${f(a.wss_max_systolic_Pa / pre.summary.wss_max, 2)}× the mean-flow value.`
    : `Fixed-flow mode prescribes the inflow waveform (PI ${f(params.pulsatilityIndex, 2)}); switch to Brain territory to see damping by the lesion and the downstream bed.`) + damp + rev;
}

// ───────────────────────── oxygen delivery ─────────────────────────
const OEF_TXT = { normal: "normal", compensated: "compensated (misery perfusion)", unmet: "demand unmet" };
function renderOxygen() {
  const a = pre.summary, b = post?.summary ?? null;
  const f = (v, d = 1) => (v == null || !isFinite(v) ? "—" : (+v.toFixed(d)).toLocaleString());
  const row = (k, x, y, cls = "") => `<tr><td>${k}</td><td>${x}</td><td class="${y == null ? "na" : cls}">${y ?? "—"}</td></tr>`;
  const oefCell = (o) => (o?.required_OEF == null ? "—" : `${f(o.required_OEF, 2)} <span class="badge oef-${o.band}">${OEF_TXT[o.band]}</span>`);
  const dir = (x, y, upGood) => (y == null || Math.abs(y - x) < 1e-9 ? "" : (y > x) === upGood ? "up" : "down");
  const ad = a.o2_distal, bd = b?.o2_distal, ap = a.o2_perforators, bp = b?.o2_perforators;
  let html = `<thead><tr><th></th><th>Before</th><th>After</th></tr></thead>` +
    row("Arterial O₂ content, mL O₂/dL", f(a.CaO2, 2), b && f(b.CaO2, 2)) +
    row("Delivered into the vessel, mL O₂/min", f(a.o2_inlet.DO2_mL_min, 2), b && f(b.o2_inlet.DO2_mL_min, 2), b && dir(a.o2_inlet.DO2_mL_min, b.o2_inlet.DO2_mL_min, true)) +
    row("Delivered to distal territory, mL O₂/min", f(ad.DO2_mL_min, 2), bd && f(bd.DO2_mL_min, 2), bd && dir(ad.DO2_mL_min, bd.DO2_mL_min, true)) +
    row("Territory CBF, mL/100 g/min", f(ad.CBF_mL_100g_min), bd && f(bd.CBF_mL_100g_min), bd && dir(ad.CBF_mL_100g_min, bd.CBF_mL_100g_min, true)) +
    row("Required O₂ extraction (distal)", oefCell(ad), bd ? oefCell(bd) : null);
  if (ap) html += row("Delivered to perforators, mL O₂/min", f(ap.DO2_mL_min, 3), bp && f(bp.DO2_mL_min, 3), bp && dir(ap.DO2_mL_min, bp.DO2_mL_min, true)) +
    row("Required O₂ extraction (deep)", oefCell(ap), bp ? oefCell(bp) : null);
  $("#oxygen").innerHTML = html;
  const tissue = ad.CBF_mL_100g_min == null
    ? "Fixed-flow mode prescribes the flow, so there is no territory to compare against: switch to Brain territory for CBF and required extraction."
    : (() => {
      // healthy delivery: reference flow at a normal hematocrit (45%) with the same SaO2/PaO2
      const healthyTissue = (params.cbfNormal * oxygenContent({ ...params, hct: 45 })) / 100;
      const pct = (o) => f((100 * o.DO2_tissue_mL_100g_min) / healthyTissue, 0);
      return `O₂ content is the same upstream and downstream (${f(a.CaO2, 1)} mL/dL): a blockage cuts delivery only through flow. ` +
        `The distal territory receives ${pct(ad)}% of healthy O₂ delivery${bd ? ` before and ${pct(bd)}% after treatment` : ""}. ` +
        `Demand ${params.cmro2} mL O₂/100 g/min; maximum extraction ${params.oefMax}.`;
    })();
  $("#oxygen-note").textContent = tissue + (params.hct !== 45 ? ` Hematocrit ${params.hct}%: viscosity ×${(a.mu_effective / params.mu).toFixed(2)} in the vessel, ×${a.distal_resistance_scale.toFixed(2)} downstream.` : "");
}

// ───────────────────────── secondary + safety endpoints ─────────────────────────
// Lesion span ∪ stent/balloon span (mm); measured at the same positions before and after.
function treatedSpan() {
  const L = params.lengthMm, sc = params.stenosisCenter * L, h = params.stenosisLengthMm / 2;
  let a = null, b = null;
  if (params.occlusionPct > 0) { a = sc - h; b = sc + h; }
  const o = treatment?.overlay;
  const span = o?.type === "stent" ? [o.startMm, o.endMm] : o?.type === "balloon" ? [o.centerMm - o.lengthMm / 2, o.centerMm + o.lengthMm / 2] : null;
  if (span) { a = a == null ? span[0] : Math.min(a, span[0]); b = b == null ? span[1] : Math.max(b, span[1]); }
  return a == null ? null : [Math.max(0, a), Math.min(L, b)];
}

function spanMetrics(res, [a, b]) {
  const sm = Array.from(res.geom.s, (v) => v * 1e3), w = res.wss;
  const pa = interp1(a, sm, res.pressure), pb = interp1(b, sm, res.pressure);
  let wInt = 0, len = 0, wMax = 0;
  for (let i = 1; i < sm.length; i++) {
    const x0 = Math.max(sm[i - 1], a), x1 = Math.min(sm[i], b);
    if (x1 > x0) { wInt += 0.5 * (w[i - 1] + w[i]) * (x1 - x0); len += x1 - x0; }
  }
  for (let i = 0; i < sm.length; i++) if (sm[i] >= a - 1e-9 && sm[i] <= b + 1e-9) wMax = Math.max(wMax, w[i]);
  return { dP: (pa - pb) / MMHG, pdpa: pb / pa, wssMean: len ? wInt / len : w[0], wssMax: wMax };
}

const areaTxt = (v) => (v > 99.5 && v < 100 ? v.toFixed(2) : v.toFixed(0));
const lengthWhere = (res, test) => (Array.from(res.wss).filter(test).length / res.wss.length) * params.lengthMm;

// thin-wall estimate of the pressure needed to stretch the healthy wall to the device's nominal diameter
function wallContact(dev) {
  if (!dev || !treatment || !(dev.stent || dev.category?.startsWith("Stent retriever"))) return null;
  const o = treatment.overlay;
  const centerMm = o.type === "stent" ? (o.startMm + o.endMm) / 2 : o.centerMm;
  const g = pre.geom, sm = Array.from(g.s, (v) => v * 1e3);
  const stretch = interp1(centerMm, sm, pre.radius) / interp1(centerMm, sm, g.radius);
  const dv = 2e3 * interp1(centerMm, sm, g.base) * stretch;  // loaded healthy diameter, mm
  const dn = dev.stent ? o.diameterMm : Math.min(deviceOpts.adjustDiameterMm ?? treatment.size.diameterMm, treatment.size.diameterMm);
  const eps = Math.max(0, (dn - dv) / dv);
  const pMmHg = params.rigid ? null : (params.stiffnessKpa * 1e3 * params.wallMm * 1e-3 * eps) / (dv / 2 * 1e-3) / MMHG;
  return { dv, dn, eps, pMmHg, retriever: !dev.stent };
}

function renderEndpoints(dev) {
  const a = pre.summary, b = post?.summary ?? null, span = treatedSpan() ?? [0, params.lengthMm];
  const segA = spanMetrics(pre, span), segB = post ? spanMetrics(post, span) : null;
  const n = (v, d) => (v == null || !isFinite(v) ? null : (+v.toFixed(d)).toLocaleString());
  const pct = (f) => (f == null ? "" : ` <span class="muted">(${(f * 100).toFixed(0)}%)</span>`);
  const row = (label, x, y, upGood = null, raw = null) => {
    let cls = "";
    if (raw && y != null && upGood != null && Math.abs(raw[1] - raw[0]) > 1e-6 * Math.max(1, Math.abs(raw[0]))) cls = (raw[1] > raw[0]) === upGood ? "up" : "down";
    return `<tr><td>${label}</td><td>${x ?? "—"}</td><td class="${y == null ? "na" : cls}">${y ?? "—"}</td></tr>`;
  };
  const group = (label) => `<tr class="group"><td colspan="3">${label}</td></tr>`;
  const head = `<thead><tr><th></th><th>Before</th><th>After</th></tr></thead>`;
  const hasBr = params.nBranches > 0;
  const spanTxt = treatedSpan() ? `${span[0].toFixed(1)}–${span[1].toFixed(1)} mm` : "whole segment, no lesion";

  let sec = head + group(`Treated segment · ${spanTxt}`) +
    row("Pressure drop ΔP, mmHg", n(segA.dP, 2), segB && n(segB.dP, 2), false, [segA.dP, segB?.dP]) +
    row("Pd/Pa across segment", n(segA.pdpa, 3), segB && n(segB.pdpa, 3), true, [segA.pdpa, segB?.pdpa]) +
    group("Downstream flow distribution, mL/min") +
    row("Distal territory", `${qFmt(a.Q_territory_mL_min)}${pct(a.flow_fraction)}`, b && `${qFmt(b.Q_territory_mL_min)}${pct(b.flow_fraction)}`, true, [a.Q_territory_mL_min, b?.Q_territory_mL_min]);
  if (a.Q_collateral_mL_min > 0 || b?.Q_collateral_mL_min > 0) sec +=
    row("· through the vessel outlet", qFmt(a.Q_distal_mL_min), b && qFmt(b.Q_distal_mL_min)) +
    row("· through collaterals", a.Q_collateral_mL_min < 0.0005 ? "0" : qFmt(a.Q_collateral_mL_min), b && (b.Q_collateral_mL_min < 0.0005 ? "0" : qFmt(b.Q_collateral_mL_min)));
  if (hasBr) sec +=
    row(`Perforators (${params.nBranches})`, `${qFmt(a.Q_branches_mL_min)}${pct(a.branch_fraction)}`, b && `${qFmt(b.Q_branches_mL_min)}${pct(b.branch_fraction)}`, true, [a.Q_branches_mL_min, b?.Q_branches_mL_min]) +
    row("Perforator share of inflow", `${(100 * a.Q_branches_mL_min / a.Q_mL_min).toFixed(2)}%`, b && `${(100 * b.Q_branches_mL_min / b.Q_mL_min).toFixed(2)}%`);
  const muTxt = (x) => { const lo = x.mu_apparent_min * 1e3, hi = x.mu_apparent_max * 1e3; return hi - lo < 0.005 ? lo.toFixed(2) : `${lo.toFixed(2)}–${hi.toFixed(2)}`; };
  sec += group(`Blood viscosity along the vessel, mPa·s${params.nonNewtonian ? " · shear-thinning" : ""}`) +
    row("Apparent viscosity, min–max", muTxt(a), b && muTxt(b)) +
    row("Highest vs high-shear value", `×${(a.mu_apparent_max / a.mu_effective).toFixed(2)}`, b && `×${(b.mu_apparent_max / b.mu_effective).toFixed(2)}`, false, [a.mu_apparent_max, b?.mu_apparent_max]);
  sec += group("Wall shear stress in segment, Pa") +
    row("Mean WSS", n(segA.wssMean, 2), segB && n(segB.wssMean, 2), null) +
    row("Peak WSS", n(segA.wssMax, 1), segB && n(segB.wssMax, 1), false, [segA.wssMax, segB?.wssMax]) +
    row("Change in mean WSS", "", segB && `${segB.wssMean >= segA.wssMean ? "+" : ""}${(100 * (segB.wssMean / segA.wssMean - 1)).toFixed(0)}%`);
  $("#secondary").innerHTML = sec;
  $("#secondary-note").textContent = a.flow_fraction == null
    ? "Fixed-flow mode: the inflow is prescribed, so the distribution shows how it splits, not how much reaches the territory."
    : hasBr ? "Distal = flow leaving the outlet (M2 territory). Perforators = deep territory (basal ganglia, internal capsule)." : "Add perforator branches in the Anatomy tab to see how flow divides between deep and distal territories.";
  $("#secondary-note").textContent += params.nonNewtonian
    ? " Shear-thinning blood thickens where shear is low (slow flow, wide lumen); the high-shear value is the viscosity set under Blood & oxygen."
    : " Blood is Newtonian, so viscosity is uniform; tick Shear-thinning blood under Blood & oxygen to see it vary.";

  // safety
  const lowA = lengthWhere(pre, (w) => w < LOW_WSS_PA), lowB = post ? lengthWhere(post, (w) => w < LOW_WSS_PA) : null;
  const hiA = lengthWhere(pre, (w) => w > HIGH_WSS_PA), hiB = post ? lengthWhere(post, (w) => w > HIGH_WSS_PA) : null;
  const starved = (s) => (s.branch_fractions ?? []).filter((f) => f < 0.5).length;
  let saf = head + group("Stagnation / low shear") +
    row(`Length with WSS < ${LOW_WSS_PA} Pa, mm`, n(lowA, 1), post && n(lowB, 1), false, [lowA, lowB]) +
    row("Minimum WSS, Pa", n(a.wss_min, 2), b && n(b.wss_min, 2), true, [a.wss_min, b?.wss_min]) +
    group("Elevated shear") +
    row("Peak WSS, Pa", n(a.wss_max, 1), b && n(b.wss_max, 1), false, [a.wss_max, b?.wss_max]) +
    row(`Length with WSS > ${HIGH_WSS_PA} Pa, mm`, n(hiA, 1), post && n(hiB, 1), false, [hiA, hiB]);
  if (hasBr && a.branch_fractions) saf +=
    group("Perforators") + row("Below 50% of healthy flow", `${starved(a)} of ${params.nBranches}`, b && `${starved(b)} of ${params.nBranches}`, false, [starved(a), b ? starved(b) : null]);
  const c = wallContact(dev);
  if (c) saf += group(`Wall contact · estimate${c.retriever ? " during deployment" : ""}`) +
    row("Oversizing vs local vessel", "", `${(c.eps * 100).toFixed(0)}% <span class="muted">(${c.dn.toFixed(2)} vs ${c.dv.toFixed(2)} mm)</span>`) +
    row("Pressure to hold that stretch", "", c.pMmHg == null ? "rigid wall" : `${c.pMmHg.toFixed(0)} mmHg <span class="muted">(${(c.pMmHg * MMHG / 1e3).toFixed(1)} kPa)</span>`);
  $("#safety").innerHTML = saf;

  // flags
  const flags = [];
  const areaPost = treatment ? 100 * (1 - throatRatio(treatment.post) ** 2) : null;
  if (params.occlusionPct > 0 && a.area_stenosis_pct >= 50)
    flags.push(["warn", `Before: a ${areaTxt(a.area_stenosis_pct)}%-area stenosis. Separation and recirculation (stagnation) are likely just downstream of the throat. Not resolved here; confirm with 3D CFD.`]);
  if (treatment && treatment.post.occlusionPct > 0 && areaPost >= 50)
    flags.push(["warn", `After: ${areaTxt(areaPost)}% area stenosis remains, so post-stenotic separation is still likely.`]);
  if (params.bulbPct >= 15) flags.push(["warn", `The ${params.bulbPct}% bulb dilation favors flow separation and low, oscillating shear on its outer wall, which this model cannot resolve.`]);
  if ((b ?? a).wss_max > HIGH_WSS_PA) flags.push(["warn", `Shear above ${HIGH_WSS_PA} Pa over ${n(b ? hiB : hiA, 1)} mm${b ? " after treatment" : ""}: a screening level for platelet activation in stenotic jets.`]);
  if (hasBr) {
    const s0 = params.stenosisCenter * params.lengthMm - params.stenosisLengthMm / 2, s1 = s0 + params.stenosisLengthMm;
    const origins = pre.geom.branchIdx.map((i) => pre.geom.s[i] * 1e3);
    const inPlaque = params.occlusionPct > 0 ? origins.filter((x) => x >= s0 && x <= s1).length : 0;
    if (inPlaque) flags.push(["warn", `${inPlaque} perforator origin${inPlaque > 1 ? "s lie" : " lies"} within the ${params.lesion === "thrombus" ? "clot" : "plaque"}. Origin narrowing itself is not modeled, only the pressure there.`]);
    if (dev?.stent && treatment) {
      const o = treatment.overlay, jailed = origins.filter((x) => x >= o.startMm && x <= o.endMm).length;
      if (jailed) flags.push(["bad", `The stent covers ${jailed} perforator origin${jailed > 1 ? "s" : ""}. Plaque shifted into a covered origin ("snowplow") can occlude it, a recognized cause of perforator stroke after intracranial stenting. That risk is not modeled.`]);
    }
  }
  if (c && c.eps > 0.1) flags.push(["warn", `Nominal oversizing of ${(c.eps * 100).toFixed(0)}% exceeds the linear wall law's ~10% range. Treat the contact pressure as an upper bound; self-expanding devices usually stop short of nominal size.`]);
  if (dev?.balloon) flags.push(["ok", "The balloon is sized below the vessel (~80%), so it loads the plaque rather than the healthy wall. Plaque mechanics and dissection risk are not modeled."]);
  if (!flags.length) flags.push(["ok", "No screening flags."]);
  $("#safety-flags").innerHTML = flags.map(([lvl, txt]) => `<li class="${lvl}">${esc(txt)}</li>`).join("");
}

function updateProbe() {
  if (!result) return;
  const n = result.radius.length, i = Math.round(probeFrac * (n - 1));
  const s = result.geom.s[i] * 1e3;
  $("#probe-s").textContent = `${s.toFixed(1)} mm`;
  $("#probe-head").textContent = `s = ${s.toFixed(1)} mm`;
  viewer.setProbe(probeFrac);
  for (const c of Object.values(charts)) c.setProbe(probeFrac);
  const rows = [
    ["Diameter", `${fmt(2e3 * result.radius[i])} mm`],
    ["Mean / peak velocity", `${fmt(result.velocity[i])} / ${fmt(2 * result.velocity[i])} m/s`],
    ["Reynolds number", `${result.reynolds[i].toFixed(0)}`],
    ["Pressure", `${fmt(result.pressure[i] / MMHG)} mmHg`],
    ["Wall shear stress", `${fmt(result.wss[i])} Pa`],
    ["Dean number", fmt(result.dean[i])],
  ];
  $("#probe-metrics").innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
}

// ───────────────────────── anatomy ─────────────────────────
const BASIS_ROWS = [
  ["diameterMm", "Inlet diameter", (p) => `${+p.diameterMm.toFixed(2)} mm`],
  ["taper", "Outlet diameter", (p) => `${+(p.diameterMm * p.taper).toFixed(2)} mm`],
  ["bulbPct", "Bulb peak", (p) => (p.bulbPct > 0 ? `+${p.bulbPct.toFixed(0)}%` : "none")],
  ["lengthMm", "Length", (p) => `${+p.lengthMm.toFixed(1)} mm`],
  ["tortuosity", "Tortuosity", (p) => p.tortuosity.toFixed(2)],
  ["wallMm", "Wall", (p) => `${+p.wallMm.toFixed(2)} mm`],
  ["flowMlMin", "Healthy flow", (p) => `${qFmt(p.flowMlMin)} mL/min`],
  ["inletMmHg", "Inlet pressure", (p) => `${p.inletMmHg.toFixed(0)} mmHg`],
  ["nBranches", "Perforators", (p) => (p.nBranches > 0 ? `${p.nBranches} · ${p.branchFlowPct}% of flow` : "none")],
];
let anatomyKey = "";

function buildAnatomy() {
  $("#sections").innerHTML = SECTIONS.map((s, i) =>
    `<button role="radio" data-sec="${s.id}"><span>${esc(s.name)}</span><small>${i + 1}</small></button>`).join("");
  $("#sections").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (b && b.dataset.sec !== params.section) applySection(b.dataset.sec);
  });
  $("#geom-source").addEventListener("change", (e) => applySection(params.section, e.target.value));
  $("#section-scenarios").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b) return;
    const sec = sectionById(params.section);
    applySection(sec.id, params.source, sec.scenarios[+b.dataset.i].p);
  });
  document.querySelectorAll(".side-tabs button").forEach((b) => b.addEventListener("click", () => {
    document.querySelectorAll(".side-tabs button").forEach((x) => { x.classList.toggle("active", x === b); x.setAttribute("aria-selected", x === b); });
    document.querySelectorAll(".side-body").forEach((x) => (x.hidden = x.dataset.sideBody !== b.dataset.side));
  }));
}

function renderAnatomy() {
  const sec = sectionById(params.section), srcId = params.source;
  $("#side-section").textContent = sec.short;
  $("#sections").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.sec === sec.id));
  if (anatomyKey !== `${sec.id}|${srcId}`) {
    anatomyKey = `${sec.id}|${srcId}`;
    $("#section-why").textContent = sec.why;
    $("#geom-source").innerHTML = sec.sources.map((s) => `<option value="${s.id}">${esc(s.label)}</option>`).join("");
    $("#section-scenarios").innerHTML = sec.scenarios.map((s, i) => `<button data-i="${i}">${esc(s.name)}</button>`).join("");
    const keys = new Set(), addRefs = (b) => [b[2] ?? []].flat().forEach((k) => keys.add(k));
    for (const src of sec.sources) Object.values(src.basis).forEach(addRefs);
    Object.values(sec.common.basis).forEach(addRefs);
    sec.scenarios.forEach((s) => Object.values(s.basis ?? {}).forEach(addRefs));
    $("#section-refs").innerHTML = [...keys].map((k) => `<li><a href="${REFS[k].url}" target="_blank" rel="noopener">${esc(REFS[k].label)}</a></li>`).join("");
    $("#section-note").textContent = {
      carotid: "Only the CCA → bulb → ICA path is modeled; the ECA branch, flow separation and recirculation in the bulb need 3D CFD (export the STL). Stenosis % is relative to the local bulb (ECST-like). The device catalog is intracranial, so its sizing checks will flag the carotid.",
      lsa: "Perforators are far smaller than any thrombectomy device; the catalog's sizing checks will flag them. Re is in single digits, so flow here is purely viscous.",
      watershed: "The territory has no autoregulation or leptomeningeal collaterals, so a pressure drop maps straight into lost flow. That makes this a worst case for the border zone.",
      m1: "M1 is the lab's reference vessel. Lenticulostriate perforators leave it but are not drawn as branches.",
      basilar: "Paramedian perforators are not modeled as branches; plaque over their origins is the clinical risk even when basilar flow is preserved.",
    }[sec.id];
  }
  $("#geom-source").value = srcId;
  const ref = sectionParams(sec, srcId), basis = sectionBasis(sec, srcId);
  $("#geom-basis").innerHTML = BASIS_ROWS.map(([k, label, show]) => {
    const b = basis[k];
    const edited = Math.abs((params[k] ?? 0) - (ref[k] ?? 0)) > 1e-9 || (k === "taper" && Math.abs(params.diameterMm - ref.diameterMm) > 1e-9);
    const refLink = [b?.[2] ?? []].flat().map((k, j, all) =>
      ` <a href="${REFS[k].url}" target="_blank" rel="noopener" title="${esc(REFS[k].label)}">ref${all.length > 1 ? j + 1 : ""}</a>`).join("");
    const why = edited ? `<span class="why">Edited from ${esc(BASIS_ROWS.find((r) => r[0] === k)[2](ref))}.</span>`
      : b ? `<span class="why"><span class="kind ${b[0]}">${b[0]}</span> ${esc(b[1])}${refLink}</span>` : "";
    return `<tr><td>${label}${why}</td><td class="v${edited ? " edited" : ""}">${show(params)}</td></tr>`;
  }).join("");
}

// ───────────────────────── devices ─────────────────────────
const GROUPS = [
  { lesion: "thrombus", label: "Acute clot · stent retrievers", test: (d) => d.category.startsWith("Stent retriever") },
  { lesion: "thrombus", label: "Acute clot · aspiration", test: (d) => d.aspiration },
  { lesion: "icas", label: "ICAS · angioplasty & stenting", test: (d) => d.balloon || d.stent },
];
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

function sceneExtras(dev) {
  return { ...lesionExtras(dev), branches: params.nBranches > 0 ? result.summary.branch_fractions ?? new Array(params.nBranches).fill(null) : null };
}

function lesionExtras(dev) {
  const lesionBefore = params.occlusionPct > 0 ? params.lesion : null;
  if (!dev || !treatment || view === "before") return { lesion: lesionBefore };
  if (view === "inplace") return { lesion: lesionBefore, overlay: treatment.overlay };
  // after: thrombectomy leaves only the underlying plaque; a failed pass leaves the clot
  const failed = !dev.balloon && !dev.stent && !deviceOpts.success;
  const lesionAfter = treatment.post.occlusionPct > 0 ? (failed ? params.lesion : "icas") : null;
  return { lesion: lesionAfter, overlay: dev.stent ? treatment.overlay : null };
}

function renderDeviceBar(dev) {
  $("#lesion").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === params.lesion));
  const sel = $("#device");
  const html = `<option value="">No device</option>` + GROUPS.filter((g) => g.lesion === params.lesion).map((g) =>
    `<optgroup label="${g.label}">${DEVICES.filter(g.test).map((d) => `<option value="${d.id}">${esc(d.name)} · ${esc(d.maker.split(" ")[0])}</option>`).join("")}</optgroup>`).join("");
  if (sel.dataset.lesion !== params.lesion) { sel.innerHTML = html; sel.dataset.lesion = params.lesion; }
  sel.value = deviceId ?? "";
  const size = $("#device-size"), viewSeg = $("#view"), tog = $("#device-toggle");
  size.hidden = viewSeg.hidden = tog.hidden = !dev;
  if (!dev) return;
  if (size.dataset.dev !== dev.id) {
    size.innerHTML = (dev.sizes.length > 1 ? `<option value="">Auto-size</option>` : "") +
      dev.sizes.map((s, i) => `<option value="${i}">${esc(s.label)}</option>`).join("");
    size.dataset.dev = dev.id;
  }
  size.value = deviceOpts.sizeIdx == null && dev.sizes.length > 1 ? "" : String(treatment?.sizeIdx ?? 0);
  if (deviceOpts.sizeIdx == null && dev.sizes.length > 1 && treatment)
    size.options[0].textContent = `Auto · ${dev.sizes[treatment.sizeIdx].label}`;
  viewSeg.querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === view));
  tog.textContent = panelOpen ? "Details ▴" : "Details ▾";
  tog.setAttribute("aria-expanded", panelOpen);
}

function optionInputs(dev) {
  const rows = [];
  const num = (key, label, hint, min, max, step) => rows.push(`<label class="opt"><span class="lab">${label}${hint ? `<small>${hint}</small>` : ""}</span>
      <input type="range" data-opt="${key}" min="${min}" max="${max}" step="${step}" value="${deviceOpts[key]}" />
      <input type="number" data-opt="${key}" min="${min}" max="${max}" step="${step}" value="${+(+deviceOpts[key]).toFixed(2)}" /></label>`);
  if (dev.balloon) {
    num("recoilPct", "Elastic recoil, %", "Assumed; no intracranial value in our sources", 0, 50, 1);
  } else if (dev.stent) {
    num("residualPct", "Residual stenosis after stenting, %", "Default = HDE study mean 31.9%", 0, 90, 0.1);
    num("isrPct", "In-stent restenosis at follow-up, %", "0 = none; ≥ 50% is usually called ISR", 0, 90, 1);
  } else {
    rows.push(`<label class="check"><input type="checkbox" data-opt="success" ${deviceOpts.success ? "checked" : ""} /> Successful pass (clot removed)</label>`);
    num("residualPct", "Underlying stenosis after clot removal, %", "0 = purely embolic clot; > 0 = clot on top of ICAS", 0, 90, 1);
    if (dev.aspiration) {
      num("vacuumKpa", "Pump vacuum, kPa below atmospheric", "Penumbra reports 98.9 kPa (−29.2 inHg) for RED 72 SILVER LABEL", 20, 101, 0.1);
      num("catheterLengthCm", "Catheter working length, cm", "Assumed; set from the IFU", 90, 160, 1);
    }
    if (dev.adjustable) {
      const sz = dev.sizes[treatment?.sizeIdx ?? 0];
      deviceOpts.adjustDiameterMm ??= sz.diameterMm;
      num("adjustDiameterMm", "Operator-set mesh diameter, mm", `Adjustable up to ${sz.diameterMm} mm`, 0.5, sz.diameterMm, 0.1);
    }
  }
  return rows.join("");
}

function compareTable() {
  const a = pre.summary, b = post.summary, rows = [
    ["Flow (mL/min)", a.Q_mL_min, b.Q_mL_min, 1, true],
    ...(a.flow_fraction != null ? [["Flow, % of healthy", a.flow_fraction * 100, b.flow_fraction * 100, 0, true]] : []),
    ["Throat diameter (mm)", a.D_throat_mm, b.D_throat_mm, 2, true],
    ["Reynolds number, throat", a.Re_throat, b.Re_throat, 0, null],
    ["Reynolds number, inlet", a.Re_inlet, b.Re_inlet, 0, null],
    ["Pressure drop ΔP (mmHg)", a.dP_mmHg, b.dP_mmHg, 2, false],
    ["Pd / Pa", a.Pd_Pa, b.Pd_Pa, 3, true],
    ["Peak WSS (Pa)", a.wss_max, b.wss_max, 1, false],
  ];
  const cell = (v, d) => (isFinite(v) ? (+v.toFixed(d)).toLocaleString() : "—");
  const body = rows.map(([k, x, y, d, upGood]) => {
    const cls = upGood == null || Math.abs(y - x) < 1e-9 * Math.max(1, Math.abs(x)) ? "" : (y > x) === upGood ? "up" : "down";
    return `<tr><td>${k}</td><td>${cell(x, d)}</td><td class="${cls}">${cell(y, d)}</td></tr>`;
  }).join("");
  const band = a.ischemia_band ? `<tr><td>Perfusion band</td><td>${a.ischemia_band.replace("-", " ")}</td><td>${b.ischemia_band.replace("-", " ")}</td></tr>` : "";
  const regime = `<tr><td>Flow regime (max)</td><td>${a.regime_max}</td><td>${b.regime_max}</td></tr>`;
  return `<table class="cmp"><thead><tr><th></th><th>Before</th><th>After</th></tr></thead><tbody>${body}${regime}${band}</tbody></table>`;
}

function renderDevicePanel(dev) {
  const panel = $("#device-panel");
  panel.hidden = !(dev && panelOpen && treatment);
  if (panel.hidden) return;
  const t = treatment, pt = dev.patents;
  const tagText = { marking: "on manufacturer's patent-marking page", family: "related patent family · not confirmed as product marking", "not-found": "no patents confirmed" }[pt.status];
  const scroll = panel.scrollTop;
  panel.innerHTML = `
    <button class="dp-close" aria-label="Close device details">×</button>
    <h3>${esc(dev.name)}</h3>
    <div class="maker">${esc(dev.maker)} · ${esc(dev.category)} · ${esc(t.size.label)}</div>
    <p class="summary">${esc(dev.summary)}</p>
    <p class="hint">${esc(dev.sizingNote)}</p>

    <h4>Procedure inputs</h4>
    ${optionInputs(dev)}

    <h4>Sizing &amp; IFU checks</h4>
    <ul class="checks">${t.checks.map((c) => `<li class="${c.level}">${esc(c.text)}</li>`).join("")}</ul>

    <h4>Device metrics</h4>
    <dl class="metrics">${t.metrics.map(([k, v, tip]) => `<dt${tip ? ` title="${esc(tip)}"` : ""}>${esc(k)}${tip ? " ⓘ" : ""}</dt><dd>${esc(v)}</dd>`).join("")}</dl>

    <h4>Hemodynamics before → after</h4>
    ${compareTable()}

    <h4>Patents <span class="tag ${pt.status}">${tagText}</span></h4>
    ${pt.nums.length ? `<div class="patents">${pt.nums.map((n) => `<a href="${patentUrl(n)}" target="_blank" rel="noopener">US ${n}</a>`).join("")}</div>` : ""}
    <p class="hint">${esc(pt.source.note ?? "")} Source: <a href="${pt.source.url}" target="_blank" rel="noopener">${esc(pt.source.label)}</a></p>

    <h4>Sources</h4>
    <ul class="srcs hint">${dev.sources.map((s) => `<li><a href="${s.url}" target="_blank" rel="noopener">${esc(s.label)}</a></li>`).join("")}</ul>
    <p class="hint">Outcome inputs (pass success, residual stenosis, recoil, restenosis) are set by you or taken from published means. The model does not simulate device–clot or device–wall mechanics; it computes the hemodynamic consequence of the resulting lumen.</p>`;
  panel.scrollTop = scroll;
  panel.querySelector(".dp-close").addEventListener("click", () => { panelOpen = false; schedule(); });
  panel.querySelectorAll("[data-opt]").forEach((el) => {
    const evt = el.type === "range" ? "input" : "change";
    el.addEventListener(evt, () => {
      const k = el.dataset.opt;
      deviceOpts = { ...deviceOpts, [k]: el.type === "checkbox" ? el.checked : parseFloat(el.value) };
      if (el.type === "range") { const twin = panel.querySelector(`input[type=number][data-opt="${k}"]`); if (twin) twin.value = el.value; }
      schedule();
    });
  });
}

function wireDevices() {
  $("#lesion").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b || b.dataset.value === params.lesion) return;
    const dev = deviceById(deviceId);
    if (dev && dev.lesion !== b.dataset.value) { deviceId = null; panelOpen = false; }
    set("lesion", b.dataset.value);
  });
  $("#device").addEventListener("change", (e) => {
    deviceId = e.target.value || null;
    deviceOpts = defaultOptions(deviceById(deviceId));
    panelOpen = !!deviceId;
    view = "after";
    schedule();
  });
  $("#device-size").addEventListener("change", (e) => {
    deviceOpts = { ...deviceOpts, sizeIdx: e.target.value === "" ? null : +e.target.value, adjustDiameterMm: null };
    schedule();
  });
  $("#view").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { view = b.dataset.value; schedule(); } });
  $("#device-toggle").addEventListener("click", () => { panelOpen = !panelOpen; schedule(); });
}

// ───────────────────────── sensitivity ─────────────────────────
// Aspiration comparison: Poiseuille estimate swept over one catheter/pump parameter.
const ASP_X = [
  { key: "vacuumKpa", label: "Pump vacuum (kPa below atm)", unit: "kPa", min: 20, max: 101 },
  { key: "catheterLengthCm", label: "Catheter working length (cm)", unit: "cm", min: 90, max: 160 },
  { key: "idIn", label: "Catheter inner diameter (in)", unit: "in", min: 0.035, max: 0.088 },
];
const ASP_Y = [
  { key: "flowMlMin", label: "Aspiration flow (mL/min)" },
  { key: "forceN", label: "Tip force index (N)" },
  { key: "re", label: "Reynolds number in catheter", re: true },
];
const ASP_COLORS = { red72: "#4cc2ff", zoom71: "#f5b041" };
let sensMode = "vessel";

let sensHasDevice = null;
function syncSensOptions(dev) {
  if (sensHasDevice === !!dev) return;
  sensHasDevice = !!dev;
  fillSensSelects(true);
}

function fillSensSelects(keep = false) {
  const xs = sensMode === "vessel" ? SENS_X : ASP_X;
  const ys = sensMode === "vessel" ? SENS_Y.filter((o) => !o.post || deviceById(deviceId)) : ASP_Y;
  const oldX = $("#sens-x").value, oldY = $("#sens-y").value;
  $("#sens-x").innerHTML = xs.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
  $("#sens-y").innerHTML = ys.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
  $("#sens-x").value = keep && xs.some((o) => o.key === oldX) ? oldX : sensMode === "vessel" ? "occlusionPct" : "vacuumKpa";
  if (keep && ys.some((o) => o.key === oldY)) $("#sens-y").value = oldY;
  $("#sens-mode").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === sensMode));
}

function buildSensitivity() {
  fillSensSelects();
  $("#sens-mode").addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b || b.dataset.value === sensMode) return;
    sensMode = b.dataset.value;
    fillSensSelects();
    computeSensitivity();
  });
  $("#sens-x").addEventListener("change", computeSensitivity);
  $("#sens-y").addEventListener("change", computeSensitivity);
}

function computeAspirationSensitivity() {
  const x = ASP_X.find((o) => o.key === $("#sens-x").value), y = ASP_Y.find((o) => o.key === $("#sens-y").value);
  const dev = deviceById(deviceId);
  const held = dev?.aspiration ? { vacuumKpa: deviceOpts.vacuumKpa, catheterLengthCm: deviceOpts.catheterLengthCm } : { ...ASPIRATION_DEFAULTS };
  const cats = DEVICES.filter((d) => d.aspiration);
  const N = 40, grid = Array.from({ length: N }, (_, k) => x.min + ((x.max - x.min) * k) / (N - 1));
  const est = (idMm, over = {}) => {
    const o = { ...held, ...over };
    return aspirationEstimate({ idMm, lengthCm: o.catheterLengthCm, vacuumKpa: o.vacuumKpa, mu: params.mu, rho: params.rho })[y.key];
  };
  let series, marker = null;
  if (x.key === "idIn") {
    series = [{ x: grid, y: grid.map((v) => est(v * 25.4)), color: "#8a97a8", label: "uniform-bore catheter" }];
    for (const c of cats) {
      const idIn = c.sizes[0].idMm / 25.4;
      series.push({ x: [idIn], y: [est(c.sizes[0].idMm)], color: ASP_COLORS[c.id], dots: true, dotRadius: 5.5, label: `${c.name} (${idIn.toFixed(3)}″)` });
    }
  } else {
    series = cats.map((c) => ({ x: grid, y: grid.map((v) => est(c.sizes[0].idMm, { [x.key]: v })), color: ASP_COLORS[c.id], label: c.name }));
    if (dev?.aspiration) marker = { x: held[x.key], y: est(dev.sizes[0].idMm) };
  }
  sensChart.opts.title = `${y.label} vs ${x.label.replace(/ \(.*\)/, "").toLowerCase()}`;
  sensChart.opts.thresholds = y.re ? [{ y: RE_LAMINAR_MAX, label: "laminar limit 2000", color: "#f5b041" }] : [];
  sensChart.opts.bands = [];
  sensChart.opts.legend = "top-left"; // curves rise to the right, so keep the legend clear of them
  sensChart.set(series, { xLabel: x.unit, marker });
  const heldTxt = [x.key !== "vacuumKpa" && `vacuum ${held.vacuumKpa} kPa`, x.key !== "catheterLengthCm" && `length ${held.catheterLengthCm} cm`].filter(Boolean).join(", ");
  const ratio = (cats[0].sizes[0].idMm / cats[1].sizes[0].idMm) ** 4;
  $("#sens-hint").textContent = `Hagen–Poiseuille for blood through a uniform-bore catheter (${heldTxt}${dev?.aspiration ? `, from the ${dev.name} panel` : ", defaults"}). ` +
    `Flow scales with ID⁴, so the 0.001″ ID difference gives RED 72 ${((ratio - 1) * 100).toFixed(1)}% more flow at equal vacuum and length. ` +
    `Tubing, the clot and proximal ID steps are ignored: a comparative index, not bench data.`;
}

function computeSensitivity() {
  if ($('[data-body="sensitivity"]').hidden) return;
  if (sensMode === "aspiration") return computeAspirationSensitivity();
  sensChart.opts.legend = null;
  $("#sens-hint").textContent = "All other inputs stay at their current values. The dot is the current configuration.";
  const xKey = $("#sens-x").value, y = SENS_Y.find((o) => o.key === $("#sens-y").value);
  const f = sensDef(xKey), dev = deviceById(deviceId);
  let hi = f.max;
  if (xKey === "occlusionPct") hi = params.basis === "diameter" ? 92 : 98;
  if (xKey === "inletDropMmHg") hi = params.inletMmHg - params.outflowMmHg - 5;
  const N = 28, xs = [], ys = [];
  const evalY = (p, n) => {
    const a = solve(p, n).summary;
    if (!y.post) return y.get(a);
    return y.get(a, solve(treat(dev, p, deviceOpts).post, n).summary);
  };
  for (let k = 0; k < N; k++) {
    const t = k / (N - 1);
    const v = f.log ? f.min * (hi / f.min) ** t : f.min + t * (hi - f.min);
    const p = { ...params, [xKey]: v, ...(xKey === "stiffnessKpa" ? { rigid: false } : {}) };
    if (validate(p).length) continue;
    try { const yy = evalY(p, 161); if (isFinite(yy)) { xs.push(f.log ? Math.log10(v) : v); ys.push(yy); } } catch { /* skip */ }
  }
  const curX = f.log ? Math.log10(params[xKey]) : params[xKey];
  const curY = pre ? (y.post ? (post ? y.get(pre.summary, post.summary) : null) : y.get(result.summary)) : null;
  sensChart.opts.title = `${y.label} vs ${f.label.toLowerCase()}${f.log ? " (log₁₀ kPa)" : ""}`;
  sensChart.opts.thresholds = y.re ? [{ y: RE_LAMINAR_MAX, label: "laminar limit 2000", color: "#f5b041" }]
    : y.oef ? [{ y: 0.45, label: "upper normal 0.45", color: "#3ecf8e" }, { y: params.oefMax, label: `maximum ${params.oefMax}`, color: "#ff5c5c" }] : [];
  sensChart.opts.bands = y.bands ? [
    { from: 0, to: 20, color: "rgba(255,92,92,.10)" }, { from: 20, to: 40, color: "rgba(255,138,61,.10)" },
    { from: 40, to: 80, color: "rgba(245,176,65,.08)" }, { from: 80, to: 200, color: "rgba(62,207,142,.06)" },
  ] : [];
  if (y.post) $("#sens-hint").textContent = `Each point applies ${dev.name} with the panel's current settings to the swept vessel. Points with no flow deficit are left out of "% of deficit".`;
  sensChart.set([{ x: xs, y: ys, dots: true }], {
    xLabel: f.log ? "log₁₀ kPa" : f.unit,
    marker: params.rigid && xKey === "stiffnessKpa" ? null : curY != null && isFinite(curY) ? { x: curX, y: curY } : null,
  });
}

// ───────────────────────── presets, hash, export ─────────────────────────
function buildPresets() {
  const sel = $("#preset");
  sel.innerHTML = `<option value="">Custom</option>` +
    SECTIONS.map((s) => `<optgroup label="${esc(s.name)}"><option value="s:${s.id}:-1">${esc(s.short)} · healthy</option>${
      s.scenarios.map((sc, i) => `<option value="s:${s.id}:${i}">${esc(s.short)} · ${esc(sc.name)}</option>`).join("")}</optgroup>`).join("") +
    `<optgroup label="MCA M1 · parameter studies">${PRESETS.map((p, i) => `<option value="p:${i}">${p.name}</option>`).join("")}</optgroup>`;
  sel.addEventListener("change", () => {
    if (sel.value === "") return;
    const [kind, a, b] = sel.value.split(":");
    if (kind === "s") {
      const sec = sectionById(a);
      applySection(sec.id, sec.sources[0].id, +b >= 0 ? sec.scenarios[+b].p : {});
    } else {
      applySection("m1", "lab", PRESETS[+a].p);
    }
    const v = sel.value; requestAnimationFrame(() => (sel.value = v));
  });
  $("#reset").addEventListener("click", () => { sel.value = "s:m1:-1"; sel.dispatchEvent(new Event("change")); });
}

function readHash() {
  try {
    const q = new URLSearchParams(location.hash.slice(1)), out = {};
    for (const [k, v] of q) {
      if (!(k in DEFAULTS) || k === "rigidSpan") continue;
      const d = DEFAULTS[k];
      out[k] = typeof d === "number" ? parseFloat(v) : typeof d === "boolean" ? v === "1" : v;
      if (typeof d === "number" && !isFinite(out[k])) delete out[k];
    }
    return out;
  } catch { return {}; }
}
function writeHash() {
  const q = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (k !== "rigidSpan" && v !== DEFAULTS[k]) q.set(k, typeof v === "boolean" ? (v ? "1" : "0") : String(+(+v).toPrecision(6) || v));
  if (deviceId) { q.set("device", deviceId); q.set("view", view); }
  if (geomMode === "patient") { q.set("geom", "patient"); q.set("lumen", ptVariant); q.set("ptcolor", ptColor); }
  history.replaceState(null, "", `${location.pathname}${location.search}${q.toString() ? `#${q}` : ""}`);
}

function download(name, data, type) {
  const a = document.createElement("a");
  a.href = typeof data === "string" && data.startsWith("data:") ? data : URL.createObjectURL(new Blob([data], { type }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function stl() {
  const nT = 64, g = result.geom, n = g.s.length, P = g.points, out = [];
  const vtx = (i, j) => {
    const r = result.radius[i], th = (2 * Math.PI * j) / nT, c = Math.cos(th), s = Math.sin(th), N = g.N[i], B = g.B[i];
    return [0, 1, 2].map((k) => ([P.x, P.y, P.z][k][i] + r * (c * N[k] + s * B[k])) * 1e3);
  };
  const ctr = (i) => [P.x[i] * 1e3, P.y[i] * 1e3, P.z[i] * 1e3];
  const facet = (a, b, c) => {
    const u = b.map((v, k) => v - a[k]), w = c.map((v, k) => v - a[k]);
    const nn = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
    const m = Math.hypot(...nn) || 1;
    out.push(` facet normal ${nn.map((v) => (v / m).toExponential(6)).join(" ")}\n  outer loop`,
      ...[a, b, c].map((p) => `   vertex ${p.map((v) => v.toExponential(6)).join(" ")}`), "  endloop\n endfacet");
  };
  out.push("solid wall");
  for (let i = 0; i < n - 1; i++) for (let j = 0; j < nT; j++) {
    const a = vtx(i, j), b = vtx(i, (j + 1) % nT), c = vtx(i + 1, j), d = vtx(i + 1, (j + 1) % nT);
    facet(a, b, d); facet(a, d, c);
  }
  out.push("endsolid wall", "solid inlet");
  for (let j = 0; j < nT; j++) facet(ctr(0), vtx(0, (j + 1) % nT), vtx(0, j));
  out.push("endsolid inlet", "solid outlet");
  for (let j = 0; j < nT; j++) facet(ctr(n - 1), vtx(n - 1, j), vtx(n - 1, (j + 1) % nT));
  out.push("endsolid outlet");
  return out.join("\n") + "\n";
}

function csv() {
  const g = result.geom, n = g.s.length;
  const cols = ["s_mm", "x_mm", "y_mm", "z_mm", "diameter_unloaded_mm", "diameter_loaded_mm", "curvature_1_per_mm",
    "pressure_mmHg", "velocity_m_s", "reynolds", "dean", "wss_Pa", "wall_strain"];
  const rows = [cols.join(",")];
  for (let i = 0; i < n; i++) rows.push([g.s[i] * 1e3, g.points.x[i] * 1e3, g.points.y[i] * 1e3, g.points.z[i] * 1e3,
    2e3 * g.radius[i], 2e3 * result.radius[i], g.curvature[i] * 1e-3, result.pressure[i] / MMHG, result.velocity[i],
    result.reynolds[i], result.dean[i], result.wss[i], result.strain[i]].map((v) => +v.toPrecision(6)).join(","));
  return rows.join("\n") + "\n";
}

function buildExport() {
  const btn = $("#export-btn"), menu = $("#export-menu");
  btn.addEventListener("click", (e) => { e.stopPropagation(); menu.hidden = !menu.hidden; btn.setAttribute("aria-expanded", !menu.hidden); });
  document.addEventListener("click", () => { menu.hidden = true; btn.setAttribute("aria-expanded", "false"); });
  menu.addEventListener("click", async (e) => {
    const kind = e.target.dataset.export;
    if (!kind || !result) return;
    if (geomMode === "patient" && (await exportPatient(kind))) return;
    const tag = `${params.section}_D${params.diameterMm}_TI${params.tortuosity}_occ${params.occlusionPct.toFixed(0)}_E${params.rigid ? "rigid" : Math.round(params.stiffnessKpa)}`;
    if (kind === "stl") download(`${tag}.stl`, stl(), "model/stl");
    if (kind === "csv") download(`${tag}_profiles.csv`, csv(), "text/csv");
    if (kind === "json") download(`${tag}.json`, JSON.stringify({
      inputs: params, view, results: result.summary,
      device: treatment ? { id: deviceId, name: deviceById(deviceId).name, size: treatment.size.label, options: deviceOpts,
        checks: treatment.checks, metrics: treatment.metrics, before: pre.summary, after: post.summary } : null,
      endpoints: { treated_span_mm: treatedSpan(), before: spanMetrics(pre, treatedSpan() ?? [0, params.lengthMm]),
        after: post ? spanMetrics(post, treatedSpan() ?? [0, params.lengthMm]) : null, wall_contact: wallContact(deviceById(deviceId)) },
    }, null, 2), "application/json");
    if (kind === "png") download(`${tag}.png`, viewer.snapshot());
    if (kind === "link") {
      try { await navigator.clipboard.writeText(location.href); btn.textContent = "Link copied ✓"; }
      catch { prompt("Copy this link:", location.href); }
      setTimeout(() => (btn.textContent = "Export ▾"), 1600);
    }
  });
}

// ───────────────────────── patient (CTA) geometry ─────────────────────────
const PATIENT_SID = "AHMU1218142";
const patient = new PatientViewer($("#viewport"));
let geomMode = hashExtras.get("geom") === "patient" ? "patient" : "parametric";
let ptVariant = hashExtras.get("lumen") ?? "healthy";
let ptColor = ["patch", "diameter", "wss"].includes(hashExtras.get("ptcolor")) ? hashExtras.get("ptcolor") : "patch";
let ptManifest = null, ptMeta = null, ptHealthy = null, ptMetas = null, ptResults = null;
const PT_SHORT = { healthy: "Healthy", occluded_icas70: "ICAS 70%", occluded_thrombus95: "Thrombus 95%",
  recanalized_thrombectomy: "Thrombectomy", recanalized_angioplasty: "Angioplasty", recanalized_stent: "Stent" };

async function setGeomMode(mode) {
  geomMode = mode;
  document.body.classList.toggle("patient-mode", mode === "patient");
  $("#geom-mode").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === mode));
  viewer.renderer.domElement.style.visibility = mode === "patient" ? "hidden" : "";
  patient.setActive(mode === "patient");
  if (mode === "patient") {
    document.querySelector('.tabs button[data-tab="profiles"]').click();
    await renderPatient();
  } else {
    document.querySelector(".legend").style.display = "";
    schedule();
  }
  writeHash();
}

async function renderPatient() {
  const err = $("#vp-error");
  try {
    ptManifest ??= await patient.manifest(PATIENT_SID);
    if (!ptManifest.variants.some((v) => v.key === ptVariant)) ptVariant = "healthy";
    const meta = await patient.load(PATIENT_SID, ptVariant);
    if (geomMode !== "patient") return;          // switched back while loading: leave the parametric view alone
    ptMeta = meta;
    if (!ptMetas) {
      const all = await Promise.all(ptManifest.variants.map((v) => fetch(`data/patient/${PATIENT_SID}/${v.key}.json`).then((r) => r.json())));
      ptMetas = Object.fromEntries(ptManifest.variants.map((v, i) => [v.key, all[i]]));
    }
    if (geomMode !== "patient") return;
    ptResults = Object.fromEntries(Object.entries(ptMetas).map(([k, m]) => [k, solvePatient(m, params)]));
    err.hidden = true;
  } catch (e) {
    err.textContent = `Could not load the patient geometry: ${e.message}`;
    err.hidden = false;
    return;
  }
  const list = $("#pt-variants");
  if (list && !list.children.length) {
    const groups = [["Baseline", ["healthy"]], ["Occluded", ["occluded_icas70", "occluded_thrombus95"]],
      ["Recanalized (prescribed outcome)", ["recanalized_thrombectomy", "recanalized_angioplasty", "recanalized_stent"]]];
    const by = Object.fromEntries(ptManifest.variants.map((v) => [v.key, v]));
    list.innerHTML = groups.map(([g, keys]) => `<div class="grp">${g}</div>` + keys.filter((k) => by[k]).map((k) => {
      const v = by[k], detail = v.label.includes("·") ? v.label.split("·").slice(1).join("·").trim() : v.label;
      return `<button role="radio" data-key="${k}"><b>${esc(PT_SHORT[k] ?? k)}</b><span class="flow" data-flow="${k}"></span>` +
        `<small>${esc(detail)}</small></button>`;
    }).join("")).join("");
  }
  list?.querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.key === ptVariant));
  list?.querySelectorAll("[data-flow]").forEach((el) => { el.textContent = `${ptResults[el.dataset.flow].flowPct.toFixed(0)}% flow`; });
  $("#pt-subject").textContent = ptManifest.subject;
  $("#pt-color").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === ptColor));
  // WSS scale fixed to the healthy lumen's range (floor 5 Pa) so lesions read as hot, not rescaled away
  let range = null;
  if (ptColor === "wss") {
    ptHealthy ??= await (await fetch(`data/patient/${PATIENT_SID}/healthy.json`)).json();
    if (geomMode !== "patient") return;
    range = [0, Math.max(5, 3 * Math.max(...ptHealthy.profile.wss_Pa.slice(1, -1)))];
  } else if (ptColor === "diameter") range = [0, Math.max(...ptMeta.profile.D_healthy_mm)];
  const used = patient.color(ptColor, range, ptResults[ptVariant].profile);
  const legend = document.querySelector(".legend");
  legend.style.display = used ? "" : "none";
  if (used) updateLegend({ label: ptColor === "wss" ? "WSS, 1D estimate (Pa)" : "Lumen diameter (mm)", cmap: ptColor === "wss" ? "inferno" : "viridis" }, used);
  $("#branch-legend").hidden = true;
  const key = [["#d04a3b", "M1 (lesion site)"], ["#b8665c", "bifurcation + M2"], ["#8a97a8", "flow extensions"],
    ["#4c8dff", "inlet"], ["#3ecf8e", "outlet 1 (larger M2)"], ["#f5b041", "outlet 2"]];
  $("#pt-key").innerHTML = (ptColor === "patch" ? key : key.slice(1))
    .map(([c, l]) => `<span><i style="background:${c}"></i>${l}</span>`).join("") +
    (ptColor === "patch" ? "" : `<span class="muted">M1 colored by ${ptColor === "wss" ? "Poiseuille WSS at the predicted flow" : "measured lumen diameter"}</span>`);
  $("#pt-attrib").innerHTML = `${esc(ptManifest.attribution)} <a href="https://doi.org/10.6084/m9.figshare.26965450" target="_blank" rel="noopener">Data</a> · ` +
    `<a href="https://github.com/jhayesf-arch/BVM/blob/main/docs/cfd-variants-${ptManifest.subject}.md" target="_blank" rel="noopener">How these were built</a>`;
  renderPatientTable();
  renderPatientOxygen();
  renderPatientCharts();
  writeHash();
}

function renderPatientTable() {
  const v = ptManifest.variants, by = Object.fromEntries(v.map((x) => [x.key, x])), R = ptResults;
  const restored = (k) => {
    const pre = R[k.includes("thrombectomy") ? "occluded_thrombus95" : "occluded_icas70"], post = R[k];
    if (!k.startsWith("recanalized") || !pre) return "—";
    const ref = R.healthy.flowPct, r = (100 * (post.flowPct - pre.flowPct)) / (ref - pre.flowPct);
    return r > 100.5 ? "full" : `${Math.max(r, 0).toFixed(0)}%`;
  };
  $("#pt-table").innerHTML = `<thead><tr><th></th><th>Distal flow</th><th>ΔP M1</th><th>Restored</th></tr></thead><tbody>` +
    v.map((x) => `<tr data-key="${x.key}" class="${x.key === ptVariant ? "sel" : ""}"><td>${PT_SHORT[x.key] ?? x.key}${x.checks_pass ? "" : " ⚠"}</td>` +
      `<td>${R[x.key].flowPct.toFixed(0)}%</td><td>${R[x.key].dpMmHg.toFixed(1)}</td><td>${restored(x.key)}</td></tr>`).join("") + `</tbody>`;
  const cur = by[ptVariant], r = R[ptVariant];
  $("#pt-note").textContent = `${cur.label}. ${cur.throat_diameter_mm ? `Measured throat ${cur.throat_diameter_mm.toFixed(2)} mm (${(100 * cur.throat_ratio).toFixed(1)}% of the local diameter). ` : ""}` +
    `Distal flow ${r.qMlMin.toFixed(1)} mL/min${r.qCollateralMlMin > 0.0005 ? ` (M1 ${r.qM1MlMin.toFixed(1)} + collaterals ${r.qCollateralMlMin.toFixed(1)})` : ""}. Flow % is relative to this healthy lumen at a normal hematocrit (45%). ` +
    `Device outcomes are prescribed lumens (app's device rules), not simulated device mechanics. Click a row or a lumen on the left to view it.`;
}

function renderPatientOxygen() {
  const v = ptManifest.variants, R = ptResults;
  const healthyTissue = (params.cbfNormal * oxygenContent({ ...params, hct: 45 })) / 100;   // normal blood, healthy lumen
  const f = (x, d = 1) => (x == null || !isFinite(x) ? "—" : (+x.toFixed(d)).toLocaleString());
  $("#pt-oxygen").innerHTML = `<thead><tr><th></th><th>O₂ to territory</th><th>CBF</th><th>Extraction needed</th></tr></thead><tbody>` +
    v.map((x) => {
      const o = R[x.key].oxygen;
      return `<tr data-key="${x.key}" class="${x.key === ptVariant ? "sel" : ""}"><td>${PT_SHORT[x.key] ?? x.key}</td>` +
        `<td>${f(o.DO2_mL_min, 1)} <span class="muted">(${f((100 * o.DO2_tissue_mL_100g_min) / healthyTissue, 0)}%)</span></td>` +
        `<td>${f(o.CBF_mL_100g_min, 0)}</td><td>${f(o.required_OEF, 2)} <span class="badge oef-${o.band}">${OEF_TXT[o.band].replace(" (misery perfusion)", "")}</span></td></tr>`;
    }).join("") + `</tbody>`;
  const r = R[ptVariant];
  $("#pt-oxygen-note").textContent =
    `mL O₂/min, with % of a healthy vessel's delivery at normal blood (45% hematocrit). CBF in mL/100 g/min. ` +
    `O₂ content is ${f(r.CaO2, 1)} mL/dL in every lumen, upstream and downstream, so lesions change delivery only through flow. ` +
    `Extraction needed = CMRO₂ ÷ (CBF·CaO₂/100): ≤ 0.45 normal; up to ${params.oefMax} compensated (misery perfusion); above that, demand unmet. ` +
    `${params.hct !== 45 ? `At ${params.hct}% hematocrit, viscosity in this M1 is ×${(r.muEff / 0.0035).toFixed(2)} and downstream resistance ×${r.distalScale.toFixed(2)}. ` : ""}` +
    `This geometry has no perforators (removed during segmentation), and collaterals and autoregulation are not modeled.`;
}

function renderPatientCharts() {
  const p = ptResults[ptVariant].profile, s = p.s_mm, xr = [0, s[s.length - 1]];
  charts.diam.set([
    { x: s, y: p.D_healthy_mm, color: "#8a97a8", dashed: true, label: "healthy" },
    { x: s, y: p.D_mm, label: PT_SHORT[ptVariant] ?? ptVariant },
  ], { xLabel: "mm along M1", xRange: xr });
  charts.pres.set([{ x: s, y: p.p_mmHg }], { xLabel: "mm along M1 (re: inlet)", xRange: xr });
  charts.re.set([{ x: s, y: p.re }], { xLabel: "mm along M1", xRange: xr });
  charts.wss.set([{ x: s, y: p.wss_Pa }], { xLabel: "mm along M1", xRange: xr });
}

async function exportPatient(kind) {
  const tag = `${ptManifest?.subject ?? PATIENT_SID}_${ptVariant}`;
  if (kind === "stl") { download(`${tag}_ext_in10D_out5D_mm.stl`, await (await fetch(`data/patient/${PATIENT_SID}/${ptVariant}.stl`)).blob(), "model/stl"); return true; }
  if (kind === "png") { download(`${tag}.png`, patient.snapshot()); return true; }
  if (kind === "csv") {
    const p = ptResults[ptVariant].profile, rows = ["s_mm,D_mm,D_healthy_mm,p_rel_inlet_mmHg,reynolds,wss_Pa"];
    p.s_mm.forEach((v, i) => rows.push([v, p.D_mm[i], p.D_healthy_mm[i], p.p_mmHg[i], p.re[i], p.wss_Pa[i]].join(",")));
    download(`${tag}_m1_profile.csv`, rows.join("\n") + "\n", "text/csv");
    return true;
  }
  if (kind === "json") {
    const r = ptResults[ptVariant];
    download(`${tag}.json`, JSON.stringify({ manifest: ptManifest, variant: ptVariant, blood: { hct: params.hct, hctCouple: params.hctCouple, sao2: params.sao2, pao2: params.pao2, cmro2: params.cmro2 },
      result: { distal_flow_mL_min: r.qMlMin, flow_pct: r.flowPct, dP_M1_mmHg: r.dpMmHg, CaO2: r.CaO2, oxygen: r.oxygen }, profile: r.profile }, null, 2), "application/json");
    return true;
  }
  return false;
}

// Null-safe: for ~10 min after a release GitHub Pages can serve a cached index.html with a newer main.js
// (or the reverse). A missing element must degrade that control, not stop the app from starting.
function wirePatient() {
  $("#geom-mode")?.addEventListener("click", (e) => { const b = e.target.closest("button"); if (b && b.dataset.value !== geomMode) setGeomMode(b.dataset.value); });
  $("#pt-variants")?.addEventListener("click", (e) => { const b = e.target.closest("button[data-key]"); if (b) { ptVariant = b.dataset.key; renderPatient(); } });
  $("#pt-color")?.addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { ptColor = b.dataset.value; renderPatient(); } });
  $("#pt-table")?.addEventListener("click", (e) => { const r = e.target.closest("tr[data-key]"); if (r) { ptVariant = r.dataset.key; renderPatient(); } });
  $("#pt-oxygen")?.addEventListener("click", (e) => { const r = e.target.closest("tr[data-key]"); if (r) { ptVariant = r.dataset.key; renderPatient(); } });
}

// ───────────────────────── wiring ─────────────────────────
function wireDisplay() {
  $("#field").addEventListener("change", () => { lockedRange = null; schedule(); });
  $("#lock-scale").addEventListener("change", () => { lockedRange = null; schedule(); });
  $("#opacity").addEventListener("input", (e) => viewer.setOpacity(+e.target.value));
  $("#particles").addEventListener("change", (e) => viewer.setParticles(e.target.checked));
  $("#probe").addEventListener("input", (e) => onProbe(+e.target.value));
  document.querySelectorAll(".tabs button").forEach((b) => b.addEventListener("click", () => {
    document.querySelectorAll(".tabs button").forEach((x) => x.classList.toggle("active", x === b));
    document.querySelectorAll(".tab-body").forEach((x) => (x.hidden = x.dataset.body !== b.dataset.tab));
    $("#tab-hint").textContent = b.dataset.tab === "profiles" ? "Click or drag a chart to move the probe" : "";
    if (b.dataset.tab === "sensitivity") computeSensitivity();
    else Object.values(charts).forEach((c) => c.draw());
  }));
}

applySectionRanges(sectionById(params.section));
params = clampParams(params);
buildControls();
buildAnatomy();
$("#bc").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) set("bc", b.dataset.value); });
buildPresets();
buildSensitivity();
buildExport();
wireDisplay();
wireDevices();
wirePatient();
syncControls();
compute();
setGeomMode(geomMode);
wireLanding();

// Landing page: the model choice. Shown on a fresh visit; a shared link (settings in the URL) opens the app directly.
// The top bar is unchanged; clicking its title returns here.
function wireLanding() {
  const landing = $("#landing"), lsel = $("#landing-preset");
  if (!landing || !lsel) return;
  lsel.innerHTML = $("#preset").innerHTML.replace(">Custom<", ">Current settings<");
  const open = () => { lsel.value = ""; landing.hidden = false; landing.querySelector(`button[data-landing="${geomMode}"]`)?.focus(); };
  const close = () => { landing.hidden = true; };
  const enter = (mode) => { close(); if (mode !== geomMode) setGeomMode(mode); };
  $("#landing-cards").addEventListener("click", (e) => { const b = e.target.closest("button[data-landing]"); if (b) enter(b.dataset.landing); });
  lsel.addEventListener("change", () => {
    if (lsel.value === "") return;
    const sel = $("#preset");
    sel.value = lsel.value;
    sel.dispatchEvent(new Event("change"));
    enter("parametric");
  });
  const brand = document.querySelector(".topbar .brand");
  if (brand) { brand.style.cursor = "pointer"; brand.title = "Back to the model choice"; brand.addEventListener("click", open); }
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && !landing.hidden) close(); });
  if (!location.hash) open();
}

// exposed for debugging and the Python parity check
window.bvm = { solve, stl, get params() { return params; }, get result() { return result; }, throatRatio };
