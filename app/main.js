import { DEFAULTS, LOW_WSS_PA, MMHG, RE_LAMINAR_MAX, RE_TURBULENT_MIN, solve, throatRatio, validate } from "./model.js";
import { Viewer, cssGradient, SLOW_MOTION } from "./viewer.js";
import { LineChart, fmt } from "./charts.js";

const $ = (sel) => document.querySelector(sel);

// ───────────────────────── control definitions ─────────────────────────
const FACTORS = [
  { key: "diameterMm", label: "Luminal diameter", unit: "mm", min: 1.5, max: 5, step: 0.05,
    hint: "Healthy M1 is roughly 2.5–3.5 mm. Narrower vessels raise shear and resistance." },
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
  { key: "stenosisLengthMm", label: "Stenosis length", unit: "mm", min: 2, max: 12, step: 0.5 },
  { key: "stenosisCenter", label: "Stenosis position", unit: "", min: 0.1, max: 0.9, step: 0.01, display: (v, p) => `${(v * p.lengthMm).toFixed(1)}`, unitFn: () => "mm" },
  { key: "waves", label: "Number of bends", unit: "", min: 0.5, max: 3, step: 0.5 },
  { key: "lengthMm", label: "Segment length", unit: "mm", min: 10, max: 40, step: 1 },
  { key: "wallMm", label: "Wall thickness", unit: "mm", min: 0.1, max: 1, step: 0.05 },
  { key: "externalMmHg", label: "External pressure", unit: "mmHg", min: 0, max: 30, step: 1 },
  { key: "mu", label: "Blood viscosity μ", unit: "Pa·s", min: 0.002, max: 0.006, step: 0.0001 },
  { key: "rho", label: "Blood density ρ", unit: "kg/m³", min: 1000, max: 1100, step: 1 },
  { key: "heartRate", label: "Heart rate (Womersley)", unit: "bpm", min: 40, max: 150, step: 1 },
];

const PRESETS = [
  { name: "Healthy M1 (baseline)", p: {} },
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
};

const SENS_X = FACTORS.map((f) => ({ key: f.key, label: `${f.label}${f.unit ? ` (${f.unit})` : ""}` }));
const SENS_Y = [
  { key: "flow_pct", label: "Flow, % of healthy", get: (s) => (s.flow_fraction ?? 1) * 100, bands: true },
  { key: "Q_mL_min", label: "Flow (mL/min)", get: (s) => s.Q_mL_min },
  { key: "Re_throat", label: "Reynolds number at throat", get: (s) => s.Re_throat, re: true },
  { key: "Re_inlet", label: "Reynolds number at inlet", get: (s) => s.Re_inlet, re: true },
  { key: "dP_mmHg", label: "Pressure drop (mmHg)", get: (s) => s.dP_mmHg },
  { key: "Pd_Pa", label: "Pd / Pa", get: (s) => s.Pd_Pa },
  { key: "wss_max", label: "Peak WSS (Pa)", get: (s) => s.wss_max },
  { key: "V_throat", label: "Throat velocity (m/s)", get: (s) => s.V_throat },
];

// ───────────────────────── state ─────────────────────────
let params = { ...DEFAULTS, ...readHash() };
let result = null;
let lockedRange = null;
let probeFrac = 0.5;
const controls = new Map(); // key -> {sync()}

// ───────────────────────── control rendering ─────────────────────────
const toSlider = (f, v) => (f.log ? (Math.log10(v / f.min) / Math.log10(f.max / f.min)) * 1000 : v);
const fromSlider = (f, t) => (f.log ? f.min * (f.max / f.min) ** (t / 1000) : +t);
const decimals = (step) => Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
const show = (f, v) => (f.display ? f.display(v, params) : f.log ? (v >= 1000 ? (+v.toPrecision(3)).toString() : v.toFixed(0)) : v.toFixed(decimals(f.step)));

function makeControl(f, container, { primary = false } = {}) {
  const wrap = document.createElement("div");
  wrap.className = `factor ${primary ? "primary" : "compact"}`;
  const id = `in-${f.key}`;
  wrap.innerHTML = `
    <div class="factor-head">
      <label for="${id}">${f.label}</label>
      <span class="num"><input type="number" id="${id}-num" aria-label="${f.label} value" /><span class="unit">${f.unitFn ? f.unitFn() : f.unit}</span></span>
    </div>
    <input type="range" id="${id}" min="${f.log ? 0 : f.min}" max="${f.log ? 1000 : f.max}" step="${f.log ? 1 : f.step}" />
    ${f.hint ? `<p class="hint">${f.hint}</p>` : ""}`;
  container.appendChild(wrap);
  const range = wrap.querySelector(`#${CSS.escape(id)}`), num = wrap.querySelector(`#${CSS.escape(id)}-num`);
  num.step = f.display ? "any" : f.log ? "any" : f.step;
  const sync = () => {
    range.value = toSlider(f, params[f.key]);
    num.value = show(f, params[f.key]);
    const disabled = f.key === "stiffnessKpa" && params.rigid;
    range.disabled = num.disabled = disabled;
    wrap.style.opacity = disabled ? 0.45 : 1;
  };
  range.addEventListener("input", () => set(f.key, fromSlider(f, +range.value)));
  num.addEventListener("change", () => {
    let v = parseFloat(num.value);
    if (!isFinite(v)) return sync();
    if (f.key === "stenosisCenter") v = v / params.lengthMm;
    set(f.key, Math.min(Math.max(v, f.min), f.max));
  });
  controls.set(f.key, { sync, wrap });
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
    const w = makeControl(f, host, { primary: true });
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
  $("#bc").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) set("bc", b.dataset.value); });
  renderBcControls();
}

function renderBcControls() {
  const host = $("#bc-controls");
  for (const key of ["flowMlMin", "inletMmHg", "outflowMmHg"]) controls.delete(key);
  host.innerHTML = "";
  for (const f of BC_CONTROLS[params.bc]) makeControl(f, host);
  $("#bc").querySelectorAll("button").forEach((b) => b.setAttribute("aria-checked", b.dataset.value === params.bc));
  $("#bc-hint").textContent = BC_HINT[params.bc];
  syncControls();
}

function syncControls() { for (const c of controls.values()) c.sync(); }

function set(key, value) {
  params = { ...params, [key]: value };
  // keep the stenosis inside the vessel
  const half = params.stenosisLengthMm / 2 / params.lengthMm;
  if (half >= 0.5) params.stenosisLengthMm = params.lengthMm * 0.9;
  const h = params.stenosisLengthMm / 2 / params.lengthMm;
  params.stenosisCenter = Math.min(Math.max(params.stenosisCenter, h), 1 - h);
  if (params.basis === "diameter" && params.occlusionPct > 95) params.occlusionPct = 95;
  if (key === "bc") renderBcControls();
  $("#preset").value = "";
  syncControls();
  schedule();
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
  try {
    result = solve(params);
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
  viewer.setProbe(probeFrac);
  updateLegend(meta, range);
  updateCharts();
  updateResults();
  updateProbe();
  writeHash();
}

let lastGeomKey = "";
function geometryChanged() {
  // reframe the camera only when the overall shape changes a lot (length/tortuosity/bends)
  const key = `${params.lengthMm}|${params.tortuosity > 1.0001}|${params.waves}`;
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

  // perfusion
  const card = $("#perfusion-card");
  $("#flow").innerHTML = `${S.Q_mL_min.toFixed(1)}<span class="muted" style="font-size:13px"> mL/min</span>`;
  if (S.flow_fraction != null) {
    card.classList.remove("disabled");
    const pct = S.flow_fraction * 100;
    $("#flow-frac").textContent = `${pct.toFixed(0)}% of healthy`;
    badge($("#band"), S.ischemia_band);
    const color = { normal: "var(--ok)", oligemia: "var(--warn)", "penumbra-range": "var(--orange)", "core-range": "var(--bad)" }[S.ischemia_band];
    Object.assign($("#flow-meter").style, { width: `${Math.min(pct, 100)}%`, background: color });
  } else {
    card.classList.add("disabled");
    $("#flow-frac").textContent = "flow is prescribed";
    badge($("#band"), "", "switch to Brain territory");
    Object.assign($("#flow-meter").style, { width: "100%", background: "#3a4658" });
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
  $("#motion-pill").textContent = `particles ${SLOW_MOTION}× slow-motion`;
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

// ───────────────────────── sensitivity ─────────────────────────
function buildSensitivity() {
  $("#sens-x").innerHTML = SENS_X.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
  $("#sens-y").innerHTML = SENS_Y.map((o) => `<option value="${o.key}">${o.label}</option>`).join("");
  $("#sens-x").value = "occlusionPct";
  $("#sens-x").addEventListener("change", computeSensitivity);
  $("#sens-y").addEventListener("change", computeSensitivity);
}

function computeSensitivity() {
  if ($('[data-body="sensitivity"]').hidden) return;
  const xKey = $("#sens-x").value, y = SENS_Y.find((o) => o.key === $("#sens-y").value);
  const f = FACTORS.find((o) => o.key === xKey);
  let hi = f.max;
  if (xKey === "occlusionPct") hi = params.basis === "diameter" ? 92 : 98;
  const N = 28, xs = [], ys = [];
  for (let k = 0; k < N; k++) {
    const t = k / (N - 1);
    const v = f.log ? f.min * (hi / f.min) ** t : f.min + t * (hi - f.min);
    const p = { ...params, [xKey]: v, ...(xKey === "stiffnessKpa" ? { rigid: false } : {}) };
    if (validate(p).length) continue;
    try { const r = solve(p, 161); xs.push(f.log ? Math.log10(v) : v); ys.push(y.get(r.summary)); } catch { /* skip */ }
  }
  const curX = f.log ? Math.log10(params[xKey]) : params[xKey];
  const curY = result ? y.get(result.summary) : null;
  sensChart.opts.title = `${y.label} vs ${f.label.toLowerCase()}${f.log ? " (log₁₀ kPa)" : ""}`;
  sensChart.opts.thresholds = y.re ? [{ y: RE_LAMINAR_MAX, label: "laminar limit 2000", color: "#f5b041" }] : [];
  sensChart.opts.bands = y.bands ? [
    { from: 0, to: 20, color: "rgba(255,92,92,.10)" }, { from: 20, to: 40, color: "rgba(255,138,61,.10)" },
    { from: 40, to: 80, color: "rgba(245,176,65,.08)" }, { from: 80, to: 200, color: "rgba(62,207,142,.06)" },
  ] : [];
  sensChart.set([{ x: xs, y: ys, dots: true }], {
    xLabel: f.log ? "log₁₀ kPa" : f.unit,
    marker: params.rigid && xKey === "stiffnessKpa" ? null : curY != null ? { x: curX, y: curY } : null,
  });
}

// ───────────────────────── presets, hash, export ─────────────────────────
function buildPresets() {
  const sel = $("#preset");
  sel.innerHTML = `<option value="">Custom</option>` + PRESETS.map((p, i) => `<option value="${i}">${p.name}</option>`).join("");
  sel.addEventListener("change", () => {
    if (sel.value === "") return;
    const keepBc = params.bc;
    params = { ...DEFAULTS, bc: keepBc, ...PRESETS[+sel.value].p };
    renderBcControls();
    schedule();
    const v = sel.value; requestAnimationFrame(() => (sel.value = v));
  });
  $("#reset").addEventListener("click", () => { sel.value = "0"; sel.dispatchEvent(new Event("change")); });
}

function readHash() {
  try {
    const q = new URLSearchParams(location.hash.slice(1)), out = {};
    for (const [k, v] of q) {
      if (!(k in DEFAULTS)) continue;
      const d = DEFAULTS[k];
      out[k] = typeof d === "number" ? parseFloat(v) : typeof d === "boolean" ? v === "1" : v;
      if (typeof d === "number" && !isFinite(out[k])) delete out[k];
    }
    return out;
  } catch { return {}; }
}
function writeHash() {
  const q = new URLSearchParams();
  for (const [k, v] of Object.entries(params)) if (v !== DEFAULTS[k]) q.set(k, typeof v === "boolean" ? (v ? "1" : "0") : String(+(+v).toPrecision(6) || v));
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
    const tag = `mca_D${params.diameterMm}_TI${params.tortuosity}_occ${params.occlusionPct.toFixed(0)}_E${params.rigid ? "rigid" : Math.round(params.stiffnessKpa)}`;
    if (kind === "stl") download(`${tag}.stl`, stl(), "model/stl");
    if (kind === "csv") download(`${tag}_profiles.csv`, csv(), "text/csv");
    if (kind === "json") download(`${tag}.json`, JSON.stringify({ inputs: params, results: result.summary }, null, 2), "application/json");
    if (kind === "png") download(`${tag}.png`, viewer.snapshot());
    if (kind === "link") {
      try { await navigator.clipboard.writeText(location.href); btn.textContent = "Link copied ✓"; }
      catch { prompt("Copy this link:", location.href); }
      setTimeout(() => (btn.textContent = "Export ▾"), 1600);
    }
  });
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

buildControls();
buildPresets();
buildSensitivity();
buildExport();
wireDisplay();
syncControls();
compute();

// exposed for debugging and the Python parity check
window.bvm = { solve, stl, get params() { return params; }, get result() { return result; }, throatRatio };
