// Reduced-order MCA vessel model — JavaScript port of bvm/ (Python is the reference
// implementation; keep the two in sync). Units: inputs in mm/kPa/mmHg/mL·min⁻¹, SI inside.

export const MMHG = 133.322;
export const MLMIN = 1e-6 / 60;
export const RE_LAMINAR_MAX = 2000;
export const RE_TURBULENT_MIN = 4000;
export const RE_STENOTIC_DISTURBANCE = 300;
export const YOUNG_TSAI_KT = 1.52;
export const LOW_WSS_PA = 0.4;
export const HIGH_WSS_PA = 10;
export const SMALL_STRAIN_LIMIT = 0.1;

export const DEFAULTS = Object.freeze({
  diameterMm: 3.0,
  tortuosity: 1.0,
  occlusionPct: 0,
  stiffnessKpa: 1000,
  rigid: false,
  lengthMm: 20,
  waves: 1,
  basis: "diameter",
  stenosisLengthMm: 6,
  stenosisCenter: 0.5,
  wallMm: 0.3,
  mu: 0.0035,
  rho: 1060,
  flowMlMin: 120,
  bc: "territory",
  inletMmHg: 85,
  outflowMmHg: 10,
  externalMmHg: 0,
  heartRate: 70,
  lesion: "icas",     // "icas" (atherosclerotic plaque) or "thrombus" (acute clot); affects rendering + device list
  rigidSpan: null,    // [startMm, endMm] held non-distensible (stented segment); null = none
  // anatomy (healthy shape; kept in the territory reference vessel)
  section: "m1",      // anatomical section id (app/anatomy.js)
  source: "lab",      // geometry source within the section
  taper: 1,           // outlet / inlet unloaded diameter, linear along arc length
  bulbPct: 0,         // peak local dilation, % of local diameter (carotid bulb)
  bulbCenter: 0.5,
  bulbLengthMm: 8,
  inletDropMmHg: 0,   // pressure lost upstream (hypotension, proximal disease); calibration ignores it
  // side branches (e.g. lenticulostriate perforators off M1), evenly spaced between the two fractions
  nBranches: 6,       // M1 default: lateral lenticulostriate perforators (see app/anatomy.js)
  branchStart: 0.3,
  branchEnd: 0.8,
  branchFlowPct: 1.9, // healthy share of inlet flow taken by all branches together
  // blood composition and oxygen (bvm/oxygen.py); mu above is the viscosity at hct = 45%
  hct: 45,            // hematocrit, %
  hctCouple: true,    // viscosity (vessel + microvascular bed) follows hematocrit, Pries et al. 1994
  microDiamUm: 25,    // effective resistance-vessel diameter downstream, µm (assumed)
  sao2: 98,           // arterial O2 saturation, %
  pao2: 95,           // mmHg
  mchc: 33,           // g/dL; Hb = Hct × MCHC
  cmro2: 3.3,         // mL O2/100 g/min
  cbfNormal: 50,      // mL/100 g/min at the healthy reference flow
  oefMax: 0.8,        // maximum O2 extraction fraction
  nonNewtonian: false, // shear-thinning blood in the vessel (Carreau, Cho & Kensey 1991)
  // downstream tissue (territory BC only; mirrors FlowConditions.collateral_pct / no_reflow_pct)
  collateralPct: 0,   // territory flow collaterals would supply with the vessel fully occluded, % of healthy
  noReflowPct: 0,     // % of the distal microvascular bed obstructed (no-reflow)
  // pulsatility (bvm/pulsatile.py): healthy inlet Gosling PI and the territory bed's time constant R·C
  pulsatilityIndex: 0.9,
  tauS: 0.2,
});

// Carreau blood (Cho & Kensey 1991) as a factor on the high-shear viscosity; mirrors bvm.hemodynamics
export const CARREAU = Object.freeze({ muRatio: 0.056 / 0.00345, lambda: 3.313, n: 0.3568 });
export const carreauFactor = (rate) => 1 + (CARREAU.muRatio - 1) * (1 + (CARREAU.lambda * rate) ** 2) ** ((CARREAU.n - 1) / 2);

export const HCT_REF = 45;

/** Relative apparent viscosity (to plasma), Pries et al. 1994 in vivo law (NIH Physiome JSim model 0157). */
export function priesRelativeViscosity(dUm, hctPct) {
  const d = dUm, h = hctPct / 100;
  const eta45 = 6 * Math.exp(-0.085 * d) + 3.2 - 2.44 * Math.exp(-0.06 * d ** 0.645);
  const w = 1 / (1 + 1e-11 * d ** 12);
  const c = (0.8 + Math.exp(-0.075 * d)) * (-1 + w) + w;
  const esl = (d / (d - 1.1)) ** 2;
  return (1 + (eta45 - 1) * (((1 - h) ** c - 1) / ((1 - 0.45) ** c - 1)) * esl) * esl;
}
export const viscosityRatio = (dUm, hctPct) => priesRelativeViscosity(dUm, hctPct) / priesRelativeViscosity(dUm, HCT_REF);

/** Effective viscosity in the vessel and microvascular resistance factor for a hematocrit. */
export function bloodFactors(p) {
  const hct = p.hct ?? HCT_REF;
  if (!(p.hctCouple ?? true) || hct === HCT_REF) return { muEff: p.mu, distalScale: 1 };
  return { muEff: p.mu * viscosityRatio(p.diameterMm * 1e3, hct), distalScale: viscosityRatio(p.microDiamUm ?? 25, hct) };
}

export const hemoglobin = (p) => ((p.hct ?? HCT_REF) * (p.mchc ?? 33)) / 100;
/** Arterial O2 content, mL O2/dL. */
export const oxygenContent = (p) => 1.34 * hemoglobin(p) * ((p.sao2 ?? 98) / 100) + 0.003 * (p.pao2 ?? 95);

export function oefBand(oef, oefMax) {
  if (oef <= 0.45) return "normal";
  if (oef <= oefMax) return "compensated";
  return "unmet";
}

export function territoryOxygen(qMlMin, qRefMlMin, p) {
  const cao2 = oxygenContent(p), out = { DO2_mL_min: (qMlMin * cao2) / 100 };
  if (qRefMlMin) {
    const cbf = ((p.cbfNormal ?? 50) * qMlMin) / qRefMlMin, tissue = (cbf * cao2) / 100;
    const oef = tissue > 0 ? (p.cmro2 ?? 3.3) / tissue : Infinity;
    Object.assign(out, { CBF_mL_100g_min: cbf, DO2_tissue_mL_100g_min: tissue, required_OEF: oef, band: oefBand(oef, p.oefMax ?? 0.8) });
  }
  return out;
}

export function flowRegime(re) {
  if (re < RE_LAMINAR_MAX) return "laminar";
  if (re < RE_TURBULENT_MIN) return "transitional";
  return "turbulent";
}

export function throatRatio(p) {
  const f = p.occlusionPct / 100;
  return p.basis === "diameter" ? 1 - f : Math.sqrt(1 - f);
}

export function validate(p) {
  const errs = [];
  if (!(p.diameterMm > 0)) errs.push("Diameter must be positive.");
  if (!(p.lengthMm > 0)) errs.push("Length must be positive.");
  if (p.tortuosity < 1) errs.push("Tortuosity index must be ≥ 1.");
  if (!(p.occlusionPct >= 0 && p.occlusionPct < 100)) errs.push("Occlusion must be in [0, 100).");
  if (!(p.stiffnessKpa > 0)) errs.push("Stiffness must be positive.");
  const half = p.stenosisLengthMm / 2, c = p.stenosisCenter * p.lengthMm;
  if (p.occlusionPct > 0 && (c - half < -1e-9 || c + half > p.lengthMm + 1e-9))
    errs.push("Stenosis extends past the vessel ends.");
  if (p.inletMmHg <= p.outflowMmHg) errs.push("Inlet pressure must exceed outflow pressure.");
  if (!(p.taper >= 0.2 && p.taper <= 2)) errs.push("Taper must be in [0.2, 2].");
  if (!(p.bulbPct >= 0 && p.bulbPct <= 150)) errs.push("Bulb dilation must be in [0, 150]%.");
  const bh = p.bulbLengthMm / 2, bc = p.bulbCenter * p.lengthMm;
  if (p.bulbPct > 0 && (bc - bh < -1e-9 || bc + bh > p.lengthMm + 1e-9)) errs.push("Bulb extends past the vessel ends.");
  if (!(p.nBranches >= 0) || !(p.branchStart > 0 && p.branchStart <= p.branchEnd && p.branchEnd < 1))
    errs.push("Branches need 0 < start ≤ end < 1.");
  if (!(p.branchFlowPct >= 0 && p.branchFlowPct < 100)) errs.push("Branch flow share must be in [0, 100)%.");
  if (!(p.hct > 0 && p.hct < 80)) errs.push("Hematocrit must be between 0 and 80%.");
  if (!((p.collateralPct ?? 0) >= 0 && (p.collateralPct ?? 0) < 100)) errs.push("Collateral supply must be in [0, 100)%.");
  if (!((p.noReflowPct ?? 0) >= 0 && (p.noReflowPct ?? 0) < 100)) errs.push("No-reflow must be in [0, 100)%.");
  if (!((p.pulsatilityIndex ?? 0.9) >= 0 && (p.pulsatilityIndex ?? 0.9) <= 2) || !((p.tauS ?? 0.2) > 0)) errs.push("Pulsatility index must be in [0, 2] and τ > 0.");
  if (p.inletDropMmHg < 0 || p.inletMmHg - p.inletDropMmHg <= p.outflowMmHg)
    errs.push("The upstream pressure drop leaves no driving pressure.");
  return errs;
}

// ---------- numerics (mirroring numpy/scipy behaviour) ----------
function linspace(a, b, n) {
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) out[i] = a + ((b - a) * i) / (n - 1);
  return out;
}

// numpy.gradient with uniform spacing, edge_order=1
function gradient(y, h) {
  const n = y.length, g = new Float64Array(n);
  g[0] = (y[1] - y[0]) / h;
  g[n - 1] = (y[n - 1] - y[n - 2]) / h;
  for (let i = 1; i < n - 1; i++) g[i] = (y[i + 1] - y[i - 1]) / (2 * h);
  return g;
}

function interp(x, xp, fp) {
  const out = new Float64Array(x.length);
  let j = 0;
  for (let i = 0; i < x.length; i++) {
    const xi = x[i];
    while (j < xp.length - 2 && xp[j + 1] < xi) j++;
    const t = (xi - xp[j]) / (xp[j + 1] - xp[j]);
    out[i] = fp[j] + Math.min(Math.max(t, 0), 1) * (fp[j + 1] - fp[j]);
  }
  return out;
}

function bisect(f, a, b, tol = 1e-13) {
  let fa = f(a);
  for (let k = 0; k < 200 && b - a > tol * Math.max(1, Math.abs(b)); k++) {
    const m = 0.5 * (a + b), fm = f(m);
    if ((fm > 0) === (fa > 0)) { a = m; fa = fm; } else b = m;
  }
  return 0.5 * (a + b);
}

const SINE_SAMPLES = 20000;
function sineArc(chord, amp, waves) {
  let L = 0, xp = 0, zp = 0;
  for (let i = 1; i < SINE_SAMPLES; i++) {
    const z = (chord * i) / (SINE_SAMPLES - 1);
    const x = amp * Math.sin((2 * Math.PI * waves * z) / chord);
    L += Math.hypot(x - xp, z - zp);
    xp = x; zp = z;
  }
  return L;
}

// ---------- geometry ----------
export function buildGeometry(p, n = 401) {
  const L = p.lengthMm * 1e-3;
  const s = linspace(0, L, n);
  const px = new Float64Array(n), py = new Float64Array(n), pz = new Float64Array(n);
  if (p.tortuosity <= 1 + 1e-9) {
    pz.set(s);
  } else {
    const chord = L / p.tortuosity;
    const amp = bisect((a) => sineArc(chord, a, p.waves) - L, 0, L);
    const zd = linspace(0, chord, SINE_SAMPLES), xd = new Float64Array(SINE_SAMPLES), sd = new Float64Array(SINE_SAMPLES);
    for (let i = 0; i < SINE_SAMPLES; i++) xd[i] = amp * Math.sin((2 * Math.PI * p.waves * zd[i]) / chord);
    for (let i = 1; i < SINE_SAMPLES; i++) sd[i] = sd[i - 1] + Math.hypot(xd[i] - xd[i - 1], zd[i] - zd[i - 1]);
    const scale = L / sd[SINE_SAMPLES - 1];
    for (let i = 0; i < SINE_SAMPLES; i++) sd[i] *= scale;
    px.set(interp(s, sd, xd));
    pz.set(interp(s, sd, zd));
  }
  const h = s[1] - s[0];
  const dx = gradient(px, h), dy = gradient(py, h), dz = gradient(pz, h);
  const T = [], N = [], B = [];
  for (let i = 0; i < n; i++) {
    const m = Math.hypot(dx[i], dy[i], dz[i]);
    T.push([dx[i] / m, dy[i] / m, dz[i] / m]);
  }
  const tx = gradient(T.map((t) => t[0]), h), ty = gradient(T.map((t) => t[1]), h), tz = gradient(T.map((t) => t[2]), h);
  const curvature = new Float64Array(n);
  for (let i = 0; i < n; i++) curvature[i] = Math.hypot(tx[i], ty[i], tz[i]);

  // rotation-minimizing frames (double reflection)
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const t0 = T[0];
  const ref = Math.abs(t0[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  let n0 = ref.map((r, k) => r - dot(ref, t0) * t0[k]);
  const n0m = Math.hypot(...n0);
  N.push(n0.map((v) => v / n0m));
  for (let i = 0; i < n - 1; i++) {
    const v1 = [px[i + 1] - px[i], py[i + 1] - py[i], pz[i + 1] - pz[i]], c1 = dot(v1, v1);
    const rL = N[i].map((v, k) => v - (2 / c1) * dot(v1, N[i]) * v1[k]);
    const tL = T[i].map((v, k) => v - (2 / c1) * dot(v1, T[i]) * v1[k]);
    const v2 = T[i + 1].map((v, k) => v - tL[k]), c2 = dot(v2, v2);
    let nn = c2 > 1e-30 ? rL.map((v, k) => v - (2 / c2) * dot(v2, rL) * v2[k]) : rL;
    const m = Math.hypot(...nn);
    N.push(nn.map((v) => v / m));
  }
  for (let i = 0; i < n; i++) {
    const t = T[i], q = N[i];
    B.push([t[1] * q[2] - t[2] * q[1], t[2] * q[0] - t[0] * q[2], t[0] * q[1] - t[1] * q[0]]);
  }

  // healthy shape (taper + bulb), then the cosine stenosis relative to the local diameter
  const r0 = (p.diameterMm / 2) * 1e-3, radius = new Float64Array(n), base = new Float64Array(n);
  const bump = (si, c, half) => { const xi = (si - c) / half; return Math.abs(xi) < 1 ? 0.5 * (1 + Math.cos(Math.PI * xi)) : 0; };
  const depth = p.occlusionPct > 0 ? 1 - throatRatio(p) : 0;
  for (let i = 0; i < n; i++) {
    let r = r0 * (1 + ((p.taper ?? 1) - 1) * (s[i] / L));
    if (p.bulbPct > 0) r *= 1 + (p.bulbPct / 100) * bump(s[i], p.bulbCenter * L, (p.bulbLengthMm / 2) * 1e-3);
    base[i] = r;
    radius[i] = depth ? r * (1 - depth * bump(s[i], p.stenosisCenter * L, (p.stenosisLengthMm / 2) * 1e-3)) : r;
  }
  let throat = 0;
  for (let i = 1; i < n; i++) {
    if (depth ? radius[i] / base[i] < radius[throat] / base[throat] : radius[i] < radius[throat]) throat = i;
  }
  const branchIdx = [];
  for (let k = 0; k < (p.nBranches ?? 0); k++) {
    const frac = p.branchStart + ((p.branchEnd - p.branchStart) * (k + 0.5)) / p.nBranches;
    branchIdx.push(Math.min(Math.max(Math.floor(frac * (n - 1) + 0.5), 1), n - 2));
  }
  return { s, points: { x: px, y: py, z: pz }, T, N, B, curvature, radius, base, throat, branchIdx };
}

export function whiteCurvatureFactor(de) {
  return de > 11.6 ? 1 / (1 - Math.pow(1 - Math.pow(11.6 / de, 0.45), 1 / 0.45)) : 1;
}

// Solve at inlet flow q. Side branches either take prescribed flows (opts.branchQ, used to
// calibrate) or drain through resistances opts.branchR to opts.pOut. Flow is piecewise constant
// between branch origins; pressure is marched from the inlet so each branch sees the local
// pressure, including any upstream stenosis loss. Mirrors bvm.hemodynamics._solve_at_flow.
export function solveAtFlow(geom, p, q, pIn, opts = {}) {
  const { s, curvature: kappa, radius: r0, branchIdx } = geom, n = s.length, nb = branchIdx.length;
  const { branchQ = null, branchR = null, pOut = 0 } = opts;
  const mu = p.mu, rho = p.rho, pExt = p.externalMmHg * MMHG;
  const eh = p.rigid ? Infinity : p.stiffnessKpa * 1e3 * p.wallMm * 1e-3;
  const iT = geom.throat, baseRatio = (geom.base[iT] / geom.base[0]) ** 2, occ = p.occlusionPct > 0;
  const sEnd = p.stenosisCenter * p.lengthMm * 1e-3 + (p.stenosisLengthMm / 2) * 1e-3;
  const ramp = new Float64Array(n);
  if (occ) for (let i = 0; i < n; i++) ramp[i] = Math.min(Math.max((s[i] - s[iT]) / Math.max(sEnd - s[iT], 1e-12), 0), 1);
  const at = new Map();
  branchIdx.forEach((idx, k) => at.set(idx, [...(at.get(idx) ?? []), k]));

  let r = Float64Array.from(r0);
  const vel = new Float64Array(n), re = new Float64Array(n), dean = new Float64Array(n), fc = new Float64Array(n), qn = new Float64Array(n);
  const pr = new Float64Array(n), strain = new Float64Array(n), rNew = new Float64Array(n), bq = new Float64Array(nb), muN = new Float64Array(n);
  const nonNewt = !!p.nonNewtonian;
  // under-relaxation ω = 0.5, halved (to 1/64) if the update stalls; mirrors bvm.hemodynamics
  let converged = false, qq = q, omega = 0.5, best = Infinity, stall = 0;
  const node = (i, qf) => {
    const m = nonNewt ? mu * carreauFactor((4 * qf) / (Math.PI * r[i] ** 3)) : mu;
    const v = qf / (Math.PI * r[i] * r[i]), rei = (rho * v * 2 * r[i]) / m, de = rei * Math.sqrt(r[i] * kappa[i]), f = whiteCurvatureFactor(de);
    return [v, rei, de, f, ((8 * m * qf) / (Math.PI * r[i] ** 4)) * f, m];
  };
  for (let it = 0; it < 4000; it++) {
    const a0 = Math.PI * r[0] ** 2 * baseRatio;
    let dp = 0, loss = 0, gPrev = 0, maxDiff = 0;
    qq = q;
    for (let i = 0; i < n; i++) {
      if (i > 0) dp += 0.5 * (gPrev + node(i, qq)[4]) * (s[i] - s[i - 1]);  // segment i-1 → i carries qq
      pr[i] = pIn - dp - loss * ramp[i];
      for (const k of at.get(i) ?? []) {
        bq[k] = branchQ ? branchQ[k] : Math.max(0, (pr[i] - pOut) / branchR[k]);
        qq -= bq[k];
      }
      [vel[i], re[i], dean[i], fc[i], gPrev, muN[i]] = node(i, qq);
      qn[i] = qq;
      if (occ && i === iT) loss = YOUNG_TSAI_KT * (rho / 2) * (qq / a0) ** 2 * (a0 / (Math.PI * r[iT] ** 2) - 1) ** 2;
    }
    for (let i = 0; i < n; i++) {
      strain[i] = p.rigid ? 0 : ((pr[i] - pExt) * r0[i]) / eh;
      if (p.rigidSpan && s[i] * 1e3 >= p.rigidSpan[0] - 1e-9 && s[i] * 1e3 <= p.rigidSpan[1] + 1e-9) strain[i] = 0;
      rNew[i] = r0[i] * (1 + Math.max(strain[i], -0.9));
      maxDiff = Math.max(maxDiff, Math.abs(rNew[i] - r[i]));
    }
    if (maxDiff < 1e-12) { converged = true; break; }
    if (maxDiff < 0.999 * best) { best = maxDiff; stall = 0; }
    else if (++stall >= 20) { omega = Math.max(omega / 2, 1 / 64); stall = 0; best = maxDiff; }
    for (let i = 0; i < n; i++) r[i] = (1 - omega) * r[i] + omega * rNew[i];
  }
  if (!converged) throw new Error("Wall–pressure coupling did not converge: the wall is very soft or the lumen is collapsing.");
  const wss = new Float64Array(n);
  for (let i = 0; i < n; i++) wss[i] = ((4 * muN[i] * qn[i]) / (Math.PI * r[i] ** 3)) * fc[i];
  return { q, qOut: qq, branchQ: bq, qNode: qn, radius: r, pressure: pr, velocity: vel, reynolds: re, dean, wss, strain, muNode: muN };
}

export function referenceParams(p) {
  return { ...p, occlusionPct: 0, tortuosity: 1, rigidSpan: null };
}

/** Solve a configuration. Returns geometry, per-station fields and a summary. */
export function solve(p, n = 401) {
  const errs = validate(p);
  if (errs.length) throw new Error(errs.join(" "));
  const geom = buildGeometry(p, n);
  const pCal = p.inletMmHg * MMHG, pIn = (p.inletMmHg - (p.inletDropMmHg ?? 0)) * MMHG, q0 = p.flowMlMin * MLMIN;
  const pOut = p.outflowMmHg * MMHG, nb = geom.branchIdx.length;
  const { muEff, distalScale } = bloodFactors(p), pEff = { ...p, mu: muEff };
  const refBq = new Float64Array(nb).fill(nb ? (q0 * p.branchFlowPct) / 100 / nb : 0);
  let branchR = new Float64Array(0), rDistal = null, ref = null, rColl = Infinity, pColl = 0, refP = null, refGeom = null;
  if (p.bc !== "fixed_flow" || nb) {
    refP = referenceParams(p); refGeom = buildGeometry(refP, n);
    ref = solveAtFlow(refGeom, refP, q0, pCal, { branchQ: refBq });
    branchR = Float64Array.from(refGeom.branchIdx, (idx, k) => ((ref.pressure[idx] - pOut) / refBq[k]) * distalScale);
    if (branchR.some((R) => !(R > 0)) || ref.pressure[n - 1] - pOut <= 0)
      throw new Error("The reference vessel alone drops more pressure than is available.");
    const rBed = ((ref.pressure[n - 1] - pOut) / ref.qOut) * distalScale;   // healthy, unobstructed bed
    rDistal = rBed / (1 - (p.noReflowPct ?? 0) / 100);
    // collaterals (mirrors bvm.hemodynamics.solve): alone they would carry collateralPct of the healthy
    // territory flow at full occlusion; supply at the healthy outlet pressure, so no flow in health
    const c = (p.collateralPct ?? 0) / 100;
    if (c > 0) rColl = rBed * (1 / c - 1);
    pColl = ref.pressure[n - 1] - (p.inletDropMmHg ?? 0) * MMHG;
  }
  let sol, referenceQ = null, qColl = 0;
  const collateral = (pEnd) => Math.max(0, (pColl - pEnd) / rColl);
  const opts = { branchR, pOut };
  if (p.bc === "fixed_flow") {
    sol = solveAtFlow(geom, pEff, q0, pIn, opts);
  } else {
    let qMax = (pIn - pOut) / rDistal;
    for (const R of branchR) qMax += (pIn - pOut) / R;
    const q = bisect((qq) => { const t = solveAtFlow(geom, pEff, qq, pIn, opts); return t.pressure[n - 1] - pOut - rDistal * (t.qOut + collateral(t.pressure[n - 1])); },
      qMax * 1e-6, qMax, qMax * 1e-12);
    sol = solveAtFlow(geom, pEff, q, pIn, opts);
    referenceQ = q0;
    qColl = collateral(sol.pressure[n - 1]);
  }
  const wssRef = (4 * p.mu * q0) / (Math.PI * ((p.diameterMm / 2) * 1e-3) ** 3);
  const result = { params: { ...p }, geom, ...sol, referenceQ, referenceQOut: ref?.qOut ?? null, referenceBranchQ: refBq,
    muEff, distalScale, wssRef, qColl };
  result.summary = summarize(result);
  result.summary.pulse = pulsatile({ p, pEff, geom, sol, qColl, q0, pIn, pCal, pOut, branchR, refBq, refP, refGeom, ref, rDistal, rColl });
  return result;
}

// ---------- pulsatility: small-signal harmonics about the mean state (mirrors bvm/pulsatile.py) ----------
const PULSE_N = 256, PULSE_OUT = 64, WAVE_AMPS = [1, 0.5, 0.25, 0.125], WAVE_PEAK = 0.15;
const cx = (re, im = 0) => ({ re, im });
const cadd = (a, b) => cx(a.re + b.re, a.im + b.im);
const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
const cinv = (a) => cdiv(cx(1), a);
const cscale = (a, k) => cx(a.re * k, a.im * k);

/** Complex harmonics c_k of the idealized inlet waveform w(t) = 1 + Re Σ c_k e^{jkωt}, scaled to a Gosling PI. */
export function waveformCoefficients(pi) {
  let gMax = -Infinity, gMin = Infinity;
  for (let i = 0; i < PULSE_N; i++) {
    let g = 0;
    WAVE_AMPS.forEach((a, k) => { g += a * Math.cos(2 * Math.PI * (k + 1) * (i / PULSE_N - WAVE_PEAK)); });
    gMax = Math.max(gMax, g); gMin = Math.min(gMin, g);
  }
  const sc = pi / (gMax - gMin);
  return WAVE_AMPS.map((a, k) => { const ph = -2 * Math.PI * (k + 1) * WAVE_PEAK; return cx(sc * a * Math.cos(ph), sc * a * Math.sin(ph)); });
}
export function synth(mean, coef, n) {
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let v = mean;
    coef.forEach((c, k) => { const th = 2 * Math.PI * (k + 1) * (i / n); v += c.re * Math.cos(th) - c.im * Math.sin(th); });
    out[i] = v;
  }
  return out;
}
function incrementalResistance(geom, p, q, pIn, branchQ) {
  const e = 1e-3 * q, n = geom.s.length;
  const lo = solveAtFlow(geom, p, q - e, pIn, { branchQ }).pressure[n - 1], hi = solveAtFlow(geom, p, q + e, pIn, { branchQ }).pressure[n - 1];
  return (lo - hi) / (2 * e);
}
function inertance(geom, radius, rho) {
  let L = 0;
  for (let i = 1; i < geom.s.length; i++) L += 0.5 * (1 / (Math.PI * radius[i] ** 2) + 1 / (Math.PI * radius[i - 1] ** 2)) * (geom.s[i] - geom.s[i - 1]);
  return rho * L;
}
function argmaxF(a) { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[k]) k = i; return k; }
function argminF(a) { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] < a[k]) k = i; return k; }

function pulsatile(c) {
  const { p, sol } = c, toMl = 1 / MLMIN, omega = (2 * Math.PI * p.heartRate) / 60;
  const coef = waveformCoefficients(p.pulsatilityIndex ?? 0.9), qBar = sol.q;
  const out = { heart_rate_bpm: p.heartRate, PI_reference: p.pulsatilityIndex ?? 0.9 };
  let qk, pk = null, qi, pi, waveT = null;
  if (p.bc !== "fixed_flow") {
    const n = c.geom.s.length, rD = c.rDistal, rDh = (c.ref.pressure[n - 1] - c.pOut) / c.ref.qOut, cap = (p.tauS ?? 0.2) / rDh;
    const rv = incrementalResistance(c.geom, c.pEff, qBar, c.pIn, sol.branchQ);
    const rvh = incrementalResistance(c.refGeom, c.refP, c.q0, c.pCal, c.refBq);
    const lv = inertance(c.geom, sol.radius, c.pEff.rho), lh = inertance(c.refGeom, c.ref.radius, c.refP.rho);
    const gc = c.qColl > 0 ? 1 / c.rColl : 0;
    qk = []; pk = []; const qtk = [];
    coef.forEach((ck, k) => {
      const w = (k + 1) * omega;
      const zinH = cadd(cx(rvh, w * lh), cinv(cx(1 / rDh, w * cap)));
      const pp = cmul(cscale(ck, c.q0), zinH);
      const y = cx(1 / rD + gc, w * cap);
      const q = cdiv(pp, cadd(cx(rv, w * lv), cinv(y)));
      pk.push(pp); qk.push(q); qtk.push(cscale(cdiv(q, y), 1 / rD));
    });
    const qTerr = sol.qOut + c.qColl;
    qi = synth(qBar, qk, PULSE_N); const qt = synth(qTerr, qtk, PULSE_N); pi = synth(c.pIn, pk, PULSE_N);
    const qtMax = max(qt), qtMin = min(qt), piMax = max(pi), piMin = min(pi);
    Object.assign(out, {
      tau_s: p.tauS ?? 0.2, incremental_resistance_ratio: rv / rvh,
      Q_territory_systolic_mL_min: qtMax * toMl, Q_territory_diastolic_mL_min: qtMin * toMl,
      PI_territory: qTerr > 0 ? (qtMax - qtMin) / qTerr : NaN,
      inlet_pulse_pressure_mmHg: (piMax - piMin) / MMHG,
    });
    waveT = Array.from(synth(qTerr, qtk, PULSE_OUT), (v) => v * toMl);
  } else {
    qk = coef.map((ck) => cscale(ck, qBar));
    qi = synth(qBar, qk, PULSE_N); pi = new Float64Array(PULSE_N).fill(c.pIn);
  }
  const qMax = max(qi), qMin = min(qi);
  Object.assign(out, {
    Q_inlet_systolic_mL_min: qMax * toMl, Q_inlet_diastolic_mL_min: qMin * toMl,
    PI_inlet: (qMax - qMin) / qBar, flow_reversal: qMin <= 0,
    waveform_t: Array.from({ length: PULSE_OUT }, (_, i) => i / PULSE_OUT),
    waveform_Q_inlet_mL_min: Array.from(synth(qBar, qk, PULSE_OUT), (v) => v * toMl),
    waveform_Q_territory_mL_min: waveT,
  });
  for (const [tag, i] of [["systolic", argmaxF(qi)], ["diastolic", argminF(qi)]]) {
    if (qi[i] <= 0) continue;
    const s = solveAtFlow(c.geom, c.pEff, qi[i], pi[i], { branchR: c.branchR, pOut: c.pOut });
    out[`wss_max_${tag}_Pa`] = max(s.wss); out[`wss_min_${tag}_Pa`] = min(s.wss);
    out[`Re_max_${tag}`] = max(s.reynolds);
    out[`dP_vessel_${tag}_mmHg`] = (s.pressure[0] - s.pressure[s.pressure.length - 1]) / MMHG;
  }
  return out;
}

function ischemiaBand(frac) {
  if (frac >= 0.8) return "normal";
  if (frac >= 0.4) return "oligemia";
  if (frac >= 0.2) return "penumbra-range";
  return "core-range";
}

/** Lévêque scaling of wall O2 mass transfer with wall shear rate, relative to the healthy inlet. */
export function wallO2Index(res) {
  const mu = res.muEff ?? res.params.mu, ref = res.wssRef / res.params.mu;
  return Array.from(res.wss, (w, i) => Math.cbrt(w / (res.muNode?.[i] ?? mu) / ref));
}

/** Linear interpolation of y(x) at x0 (x ascending), clamped like numpy.interp. */
export function interp1(x0, x, y) {
  if (x0 <= x[0]) return y[0];
  const n = x.length;
  if (x0 >= x[n - 1]) return y[n - 1];
  let j = 0;
  while (x[j + 1] < x0) j++;
  return y[j] + ((x0 - x[j]) / (x[j + 1] - x[j])) * (y[j + 1] - y[j]);
}

function argmin(a) { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] < a[k]) k = i; return k; }
function max(a) { let m = -Infinity; for (const v of a) m = Math.max(m, v); return m; }
function min(a) { let m = Infinity; for (const v of a) m = Math.min(m, v); return m; }

export function summarize(res) {
  const p = res.params, n = res.radius.length, iT = p.occlusionPct > 0 ? res.geom.throat : argmin(res.radius);
  const dIn = 2 * res.radius[0];
  const reIn = res.reynolds[0], reT = res.reynolds[iT], reMax = max(res.reynolds);
  const omega = (2 * Math.PI * p.heartRate) / 60, nu = p.mu / p.rho;
  const tr = throatRatio(p), areaSten = 100 * (1 - tr * tr);
  const entrance = 0.05 * reIn * dIn;
  const warnings = [];
  const maxStrain = max(res.strain);
  if (maxStrain > SMALL_STRAIN_LIMIT)
    warnings.push(`Wall strain reaches ${(maxStrain * 100).toFixed(0)}%, beyond the ~10% small-strain range of the linear wall law. Treat lumen dilation and Re here as qualitative.`);
  if (max(res.dean) > 2000) warnings.push("Dean number > 2000: outside the curvature correlation's range.");
  if (entrance > res.geom.s[n - 1])
    warnings.push(`Entrance length (${(entrance * 1e3).toFixed(0)} mm) is longer than the vessel, so flow is still developing and inlet WSS is underestimated.`);
  if (reMax >= RE_LAMINAR_MAX) warnings.push("Re exceeds 2000 somewhere: the laminar loss model is not valid there.");
  const out = {
    Q_mL_min: res.q / MLMIN,
    V_inlet: res.velocity[0], V_throat: res.velocity[iT], V_max: max(res.velocity),
    D_inlet_mm: dIn * 1e3, D_throat_mm: 2e3 * res.radius[iT],
    Re_inlet: reIn, Re_throat: reT, Re_max: reMax,
    regime_inlet: flowRegime(reIn), regime_throat: flowRegime(reT), regime_max: flowRegime(reMax),
    disturbance: p.occlusionPct > 0 && areaSten >= 50 && reIn >= RE_STENOTIC_DISTURBANCE,
    Dean_max: max(res.dean),
    Womersley: (dIn / 2) * Math.sqrt(omega / nu),
    entrance_mm: entrance * 1e3,
    P_in_mmHg: res.pressure[0] / MMHG, P_out_mmHg: res.pressure[n - 1] / MMHG,
    dP_mmHg: (res.pressure[0] - res.pressure[n - 1]) / MMHG,
    Pd_Pa: res.pressure[n - 1] / res.pressure[0],
    wss_inlet: res.wss[0], wss_max: max(res.wss), wss_min: min(res.wss),
    low_wss_fraction: res.wss.reduce((c, v) => c + (v < LOW_WSS_PA), 0) / n,
    high_wss: max(res.wss) > HIGH_WSS_PA,
    max_strain: maxStrain,
    wave_speed: p.rigid ? null : Math.sqrt((p.stiffnessKpa * 1e3 * p.wallMm * 1e-3) / (p.rho * dIn)),
    area_stenosis_pct: areaSten, diameter_stenosis_pct: 100 * (1 - tr),
    warnings,
  };
  out.Q_distal_mL_min = res.qOut / MLMIN;
  out.Q_collateral_mL_min = (res.qColl ?? 0) / MLMIN;
  out.Q_territory_mL_min = (res.qOut + (res.qColl ?? 0)) / MLMIN;
  out.Q_branches_mL_min = res.branchQ.reduce((a, b) => a + b, 0) / MLMIN;
  out.branch_mL_min = Array.from(res.branchQ, (b) => b / MLMIN);
  if (p.occlusionPct > 0) {
    const sc = p.stenosisCenter * p.lengthMm * 1e-3, half = (p.stenosisLengthMm / 2) * 1e-3;
    const pa = interp1(sc - half, res.geom.s, res.pressure), pb = interp1(sc + half, res.geom.s, res.pressure);
    out.dP_lesion_mmHg = (pa - pb) / MMHG;
    out.Pd_Pa_lesion = pb / pa;
  }
  // blood + oxygen (mirrors VesselSolution.summary()["blood"/"oxygen"])
  out.mu_effective = res.muEff ?? p.mu;
  out.distal_resistance_scale = res.distalScale ?? 1;
  out.mu_apparent_max = res.muNode ? max(res.muNode) : out.mu_effective;
  out.hemoglobin = hemoglobin(p);
  out.CaO2 = oxygenContent(p);
  const toMl = 1 / MLMIN, refB = res.referenceBranchQ?.reduce((a, b) => a + b, 0) ?? 0;
  out.o2_inlet = territoryOxygen(res.q * toMl, null, p);
  // tissue terms (CBF, required OEF) only exist with a territory: in fixed-flow mode the flow is prescribed
  out.o2_distal = territoryOxygen(out.Q_territory_mL_min, res.referenceQ != null ? res.referenceQOut * toMl : null, p);
  if (res.branchQ.length) out.o2_perforators = territoryOxygen(res.branchQ.reduce((a, b) => a + b, 0) * toMl, refB > 0 && res.referenceQ != null ? refB * toMl : null, p);
  if (res.wssRef) {
    const idx = wallO2Index(res);
    out.wall_o2_index_min = min(idx);
    out.wall_o2_index_below_0_8_fraction = idx.reduce((c, v) => c + (v < 0.8), 0) / idx.length;
  }
  if (res.referenceQ != null) {
    // distal (outlet) flow vs the healthy reference: the primary perfusion endpoint
    out.flow_fraction = (res.qOut + (res.qColl ?? 0)) / res.referenceQOut;
    out.ischemia_band = ischemiaBand(out.flow_fraction);
    out.reference_distal_mL_min = res.referenceQOut / MLMIN;
    const refB = res.referenceBranchQ.reduce((a, b) => a + b, 0);
    if (res.branchQ.length && refB > 0) {
      out.branch_fraction = res.branchQ.reduce((a, b) => a + b, 0) / refB;
      out.branch_fractions = Array.from(res.branchQ, (b, k) => b / res.referenceBranchQ[k]);
    }
  }
  return out;
}
