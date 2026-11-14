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
});

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

  // cosine stenosis
  const r0 = (p.diameterMm / 2) * 1e-3, radius = new Float64Array(n).fill(r0);
  if (p.occlusionPct > 0) {
    const sc = p.stenosisCenter * L, half = (p.stenosisLengthMm / 2) * 1e-3, depth = 1 - throatRatio(p);
    for (let i = 0; i < n; i++) {
      const xi = (s[i] - sc) / half;
      if (Math.abs(xi) < 1) radius[i] = r0 * (1 - depth * 0.5 * (1 + Math.cos(Math.PI * xi)));
    }
  }
  return { s, points: { x: px, y: py, z: pz }, T, N, B, curvature, radius };
}

export function whiteCurvatureFactor(de) {
  return de > 11.6 ? 1 / (1 - Math.pow(1 - Math.pow(11.6 / de, 0.45), 1 / 0.45)) : 1;
}

function solveAtFlow(geom, p, q, pIn) {
  const { s, curvature: kappa, radius: r0 } = geom, n = s.length;
  const mu = p.mu, rho = p.rho, pExt = p.externalMmHg * MMHG;
  const eh = p.rigid ? Infinity : p.stiffnessKpa * 1e3 * p.wallMm * 1e-3;
  let iT = 0;
  for (let i = 1; i < n; i++) if (r0[i] < r0[iT]) iT = i;
  const sEnd = p.stenosisCenter * p.lengthMm * 1e-3 + (p.stenosisLengthMm / 2) * 1e-3;

  let r = Float64Array.from(r0);
  const vel = new Float64Array(n), re = new Float64Array(n), dean = new Float64Array(n), fc = new Float64Array(n);
  const pr = new Float64Array(n), strain = new Float64Array(n), rNew = new Float64Array(n), g = new Float64Array(n);
  let converged = false;
  for (let it = 0; it < 200; it++) {
    for (let i = 0; i < n; i++) {
      const A = Math.PI * r[i] * r[i];
      vel[i] = q / A;
      re[i] = (rho * vel[i] * 2 * r[i]) / mu;
      dean[i] = re[i] * Math.sqrt(r[i] * kappa[i]);
      fc[i] = whiteCurvatureFactor(dean[i]);
      g[i] = ((8 * mu * q) / (Math.PI * r[i] ** 4)) * fc[i];
    }
    let loss = 0;
    if (p.occlusionPct > 0) {
      const a0 = Math.PI * r[0] ** 2, a1 = Math.PI * r[iT] ** 2;
      loss = YOUNG_TSAI_KT * (rho / 2) * (q / a0) ** 2 * (a0 / a1 - 1) ** 2;
    }
    let dp = 0, maxDiff = 0;
    for (let i = 0; i < n; i++) {
      if (i > 0) dp += 0.5 * (g[i] + g[i - 1]) * (s[i] - s[i - 1]);
      let ramp = 0;
      if (loss) ramp = Math.min(Math.max((s[i] - s[iT]) / Math.max(sEnd - s[iT], 1e-12), 0), 1);
      pr[i] = pIn - dp - loss * ramp;
      strain[i] = p.rigid ? 0 : ((pr[i] - pExt) * r0[i]) / eh;
      if (p.rigidSpan && s[i] * 1e3 >= p.rigidSpan[0] - 1e-9 && s[i] * 1e3 <= p.rigidSpan[1] + 1e-9) strain[i] = 0;
      rNew[i] = r0[i] * (1 + Math.max(strain[i], -0.9));
      maxDiff = Math.max(maxDiff, Math.abs(rNew[i] - r[i]));
    }
    if (maxDiff < 1e-12) { converged = true; break; }
    for (let i = 0; i < n; i++) r[i] = 0.5 * r[i] + 0.5 * rNew[i];
  }
  if (!converged) throw new Error("Wall–pressure coupling did not converge: the wall is too soft for this model.");
  const wss = new Float64Array(n);
  for (let i = 0; i < n; i++) wss[i] = ((4 * mu * q) / (Math.PI * r[i] ** 3)) * fc[i];
  return { q, radius: r, pressure: pr, velocity: vel, reynolds: re, dean, wss, strain };
}

function referenceParams(p) {
  return { ...p, occlusionPct: 0, tortuosity: 1, rigidSpan: null };
}

/** Solve a configuration. Returns geometry, per-station fields and a summary. */
export function solve(p, n = 401) {
  const errs = validate(p);
  if (errs.length) throw new Error(errs.join(" "));
  const geom = buildGeometry(p, n);
  const pIn = p.inletMmHg * MMHG, q0 = p.flowMlMin * MLMIN;
  let sol, referenceQ = null;
  if (p.bc === "fixed_flow") {
    sol = solveAtFlow(geom, p, q0, pIn);
  } else {
    const pOut = p.outflowMmHg * MMHG;
    const refP = referenceParams(p);
    const ref = solveAtFlow(buildGeometry(refP, n), refP, q0, pIn);
    const rDistal = (ref.pressure[n - 1] - pOut) / q0;
    if (rDistal <= 0) throw new Error("The reference vessel alone drops more pressure than is available.");
    const qMax = (pIn - pOut) / rDistal;
    const q = bisect((qq) => solveAtFlow(geom, p, qq, pIn).pressure[n - 1] - pOut - rDistal * qq, qMax * 1e-6, qMax, 1e-12);
    sol = solveAtFlow(geom, p, q, pIn);
    referenceQ = q0;
  }
  const result = { params: { ...p }, geom, ...sol, referenceQ };
  result.summary = summarize(result);
  return result;
}

function ischemiaBand(frac) {
  if (frac >= 0.8) return "normal";
  if (frac >= 0.4) return "oligemia";
  if (frac >= 0.2) return "penumbra-range";
  return "core-range";
}

function argmin(a) { let k = 0; for (let i = 1; i < a.length; i++) if (a[i] < a[k]) k = i; return k; }
function max(a) { let m = -Infinity; for (const v of a) m = Math.max(m, v); return m; }
function min(a) { let m = Infinity; for (const v of a) m = Math.min(m, v); return m; }

export function summarize(res) {
  const p = res.params, n = res.radius.length, iT = argmin(res.radius);
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
  if (res.referenceQ != null) {
    out.flow_fraction = res.q / res.referenceQ;
    out.ischemia_band = ischemiaBand(out.flow_fraction);
  }
  return out;
}
