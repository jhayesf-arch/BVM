// Browser self-test for the JS engine and device logic (verification evidence; see
// docs/credibility-assessment-plan.md). Open app/selftest.html locally or on the live site.
import * as M from "./model.js?v=202611171230";
import { ASPIRATION_DEFAULTS, DEVICES, aspirationEstimate, applyDevice, defaultOptions, treat } from "./devices.js?v=202611171230";
import * as C from "./cad.js?v=202611171230";

export function runSelfTest() {
  const results = [];
  const group = (name) => { const g = { name, cases: 0, fails: [] }; results.push(g); return g; };
  const near = (a, b, rel) => Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300);

  // ── engine: closed-form checks ──
  const eng = group("Engine: closed-form solutions");
  {
    eng.cases++;
    const r = M.solve({ ...M.DEFAULTS, rigid: true, nBranches: 0, bc: "fixed_flow" });
    const q = 120 * M.MLMIN, rad = 1.5e-3, dp = (8 * 0.0035 * q * 0.02) / (Math.PI * rad ** 4);
    if (!near(r.summary.dP_mmHg * M.MMHG, dp, 1e-12)) eng.fails.push(`Poiseuille ΔP ${r.summary.dP_mmHg * M.MMHG} vs ${dp}`);
    eng.cases++;
    const tw = (4 * 0.0035 * q) / (Math.PI * rad ** 3);
    if (!near(r.summary.wss_inlet, tw, 1e-12)) eng.fails.push(`Poiseuille WSS ${r.summary.wss_inlet} vs ${tw}`);
    eng.cases++;
    const t = M.solve({ ...M.DEFAULTS, rigid: true, nBranches: 0, bc: "territory", inletDropMmHg: 25 });
    if (!near(t.summary.flow_fraction, 50 / 75, 1e-6)) eng.fails.push(`territory pressure-ratio law ${t.summary.flow_fraction}`);
    eng.cases++;
    const b = M.solve({ ...M.DEFAULTS, nBranches: 8, branchFlowPct: 3, occlusionPct: 60 });
    const tot = b.qOut + b.branchQ.reduce((x, y) => x + y, 0);
    if (!near(tot, b.q, 1e-12)) eng.fails.push(`branch mass balance ${tot} vs ${b.q}`);
  }

  // ── blood & oxygen (mirrors tests/test_oxygen.py) ──
  const ox = group("Engine: hematocrit viscosity and oxygen delivery");
  {
    const chk = (ok, msg) => { ox.cases++; if (!ok) ox.fails.push(msg); };
    chk(Math.abs(M.priesRelativeViscosity(3000, 45) - 3.2) < 0.01, "Pries η(3 mm, 45%) ≠ 3.2");
    chk(M.priesRelativeViscosity(25, 45) < M.priesRelativeViscosity(3000, 45), "no Fåhræus–Lindqvist reduction at 25 µm");
    const f45 = M.bloodFactors(M.DEFAULTS);
    chk(f45.muEff === M.DEFAULTS.mu && f45.distalScale === 1, "Hct 45% must leave viscosity unchanged");
    const cao2 = 1.34 * (45 * 33 / 100) * 0.98 + 0.003 * 95;
    chk(near(M.oxygenContent(M.DEFAULTS), cao2, 1e-12), "CaO2 formula");
    const h = M.solve({ ...M.DEFAULTS }, 161).summary;
    chk(near(h.o2_distal.required_OEF, 3.3 / (50 * cao2 / 100), 1e-6), `baseline OEF ${h.o2_distal.required_OEF}`);
    const s = M.solve({ ...M.DEFAULTS, occlusionPct: 80 }, 161).summary;
    chk(near(s.o2_distal.DO2_mL_min / s.o2_inlet.DO2_mL_min, s.Q_distal_mL_min / s.Q_mL_min, 1e-12), "DO2 ratio ≠ flow ratio");
    const do2 = (hct) => M.solve({ ...M.DEFAULTS, hct }, 161).summary.o2_distal.DO2_mL_min;
    chk(do2(45) > do2(25) && do2(45) > do2(55), "no optimal hematocrit in a healthy vessel");
    const oef = (hct) => M.solve({ ...M.DEFAULTS, hct, occlusionPct: 70 }, 161).summary.o2_distal.required_OEF;
    chk(oef(25) > oef(35) && oef(35) > oef(45), "anemia should raise required extraction behind a stenosis");
    chk(M.solve({ ...M.DEFAULTS, occlusionPct: 95 }, 161).summary.o2_distal.band === "unmet", "95% occlusion should leave demand unmet");
  }

  // ── collaterals, no-reflow, shear-thinning blood (mirrors tests/test_collaterals_rheology.py) ──
  const tis = group("Engine: collaterals, no-reflow and non-Newtonian blood");
  {
    const chk = (ok, msg) => { tis.cases++; if (!ok) tis.fails.push(msg); };
    const S = (p) => M.solve({ ...M.DEFAULTS, ...p }, 161).summary;
    const h = S({ collateralPct: 40 });
    chk(h.Q_collateral_mL_min < 1e-6 && near(h.flow_fraction, 1, 1e-9), "collaterals must carry no flow in health");
    const occ = S({ occlusionPct: 99.5, stenosisLengthMm: 10, rigid: true, nBranches: 0, collateralPct: 30 });
    chk(Math.abs(occ.flow_fraction - 0.3) < 0.01, `collateral share at near-occlusion ${occ.flow_fraction}`);
    const g = [50, 80, 95].map((o) => S({ occlusionPct: o, collateralPct: 30 }).flow_fraction - S({ occlusionPct: o }).flow_fraction);
    chk(g[0] >= 0 && g[0] < g[1] && g[1] < g[2], "collateral gain should grow with stenosis");
    const t = S({ occlusionPct: 90, collateralPct: 25, noReflowPct: 30 });
    chk(near(t.Q_territory_mL_min, t.Q_distal_mL_min + t.Q_collateral_mL_min, 1e-12) && t.Q_collateral_mL_min > 0, "territory flow balance");
    const nr = S({ noReflowPct: 40 }).flow_fraction;
    chk(nr > 0.6 && nr < 0.62, `no-reflow 40% gives ${nr}`);
    chk(Math.abs(M.carreauFactor(0) - M.CARREAU.muRatio) < 1e-12 && M.carreauFactor(1e7) < 1.001, "Carreau limits");
    const n0 = M.solve({ ...M.DEFAULTS, nBranches: 0, bc: "fixed_flow" }, 161), n1 = M.solve({ ...M.DEFAULTS, nBranches: 0, bc: "fixed_flow", nonNewtonian: true }, 161);
    chk(Math.abs(n1.wss[0] / n0.wss[0] - M.carreauFactor((4 * n0.q) / (Math.PI * n0.radius[0] ** 3))) < 1e-3, "non-Newtonian WSS ratio");
    chk(near(S({ nonNewtonian: true }).flow_fraction, 1, 1e-9), "non-Newtonian healthy territory must calibrate to 100%");
    for (const o of [0, 70, 95]) for (const c of [0, 30]) for (const r of [0, 50]) {
      tis.cases++;
      try { S({ occlusionPct: o, collateralPct: c, noReflowPct: r, nonNewtonian: true, hct: 35 }); }
      catch (e) { tis.fails.push(`${o}% coll ${c} no-reflow ${r}: ${e.message}`); }
    }
  }

  // ── pulsatility (mirrors tests/test_pulsatile.py) ──
  const pg = group("Engine: pulsatile small-signal model");
  {
    const chk = (ok, msg) => { pg.cases++; if (!ok) pg.fails.push(msg); };
    const P = (p) => M.solve({ ...M.DEFAULTS, ...p }, 161).summary;
    for (const pi of [0, 0.6, 0.9, 1.4]) {
      const w = M.synth(1, M.waveformCoefficients(pi), 256);
      const mean = w.reduce((a, b) => a + b, 0) / w.length;
      chk(Math.abs(mean - 1) < 1e-12 && Math.abs(Math.max(...w) - Math.min(...w) - pi) < 1e-12, `waveform PI ${pi}`);
    }
    const h = P({});
    chk(near(h.pulse.PI_inlet, 0.9, 1e-9) && near(h.pulse.incremental_resistance_ratio, 1, 1e-9), "healthy reference must reproduce PI 0.9");
    const alt = P({ pulsatilityIndex: 1.3, tauS: 0.4, heartRate: 100 });
    chk(alt.Q_mL_min === h.Q_mL_min && alt.wss_max === h.wss_max, "pulsatility inputs must not change the mean results");
    chk(h.pulse.PI_territory > 0 && h.pulse.PI_territory < h.pulse.PI_inlet, "bed compliance should damp territory PI");
    const pis = [0, 50, 70, 90].map((o) => P({ occlusionPct: o }).pulse.PI_territory);
    chk(pis[0] >= pis[1] && pis[1] > pis[2] && pis[2] > pis[3], `stenosis should damp distal PI ${pis.map((x) => x.toFixed(3))}`);
    chk(P({ noReflowPct: 50 }).pulse.PI_inlet > 1.3 * h.pulse.PI_inlet, "no-reflow should raise vessel PI");
    const s = P({ occlusionPct: 60 });
    chk(s.pulse.wss_max_systolic_Pa > s.wss_max && s.wss_max > s.pulse.wss_max_diastolic_Pa, "systolic WSS > mean > diastolic");
    const z = P({ pulsatilityIndex: 0 });
    chk(Math.abs(z.pulse.PI_inlet) < 1e-12 && near(z.pulse.wss_max_systolic_Pa, z.wss_max, 1e-9), "PI 0 must give steady flow");
    const c = P({ occlusionPct: 80, collateralPct: 20 }), avg = (a) => a.reduce((x, y) => x + y, 0) / a.length;
    chk(near(avg(c.pulse.waveform_Q_inlet_mL_min), c.Q_mL_min, 1e-12) && near(avg(c.pulse.waveform_Q_territory_mL_min), c.Q_territory_mL_min, 1e-12), "waveforms must average to the mean flows");
  }

  // ── CAD import: centerline/lumen extraction on known shapes, and a parametric → STL → import round trip ──
  const cg = group("CAD import: centerline and lumen extraction");
  {
    const chk = (ok, msg) => { cg.cases++; if (!ok) cg.fails.push(msg); };
    const ext = (tri, o = {}) => C.extractLumen(C.weld(tri), { scale: 1e-3, ...o });
    const Dof = (e) => e.r.map((r) => 2e3 * r), L = (e) => e.s[e.s.length - 1] * 1e3;
    const poly = Math.sqrt((40 / (2 * Math.PI)) * Math.sin((2 * Math.PI) / 40));   // 40-gon area → equivalent diameter
    const st = (u) => [0, 0, 20 * u];
    const a = ext(C.tubeMesh(st, () => 1.5));
    chk(near(L(a), 20, 1e-6) && Dof(a).every((d) => near(d, 3 * poly, 1e-3)), `straight tube: L ${L(a)}, D ${Dof(a)[0]}`);
    chk(Math.abs(a.tortuosity - 1) < 1e-6 && Math.max(...a.x.map(Math.abs)) < 1e-9, "straight tube must have a straight centerline");
    const bump = (u) => (Math.abs(u - 0.5) < 0.15 ? 0.5 * (1 + Math.cos((Math.PI * (u - 0.5)) / 0.15)) : 0);
    const s70 = Dof(ext(C.tubeMesh(st, (u) => 1.5 * (1 - 0.7 * bump(u)), { nu: 300 })));
    chk(near(Math.min(...s70), 0.9 * poly, 5e-3), `70% stenosis throat ${Math.min(...s70)}`);
    const tp = ext(C.tubeMesh(st, (u) => 1.5 - 0.5 * u)), tr = ext(C.tubeMesh(st, (u) => 1.5 - 0.5 * u), { reverse: true });
    chk(Dof(tp)[0] > Dof(tp)[tp.r.length - 1] && Dof(tr)[0] < Dof(tr)[tr.r.length - 1], "inlet is the larger end unless reversed");
    const arc = ext(C.tubeMesh((u) => [10 - 10 * Math.cos((u * Math.PI) / 2), 0, 10 * Math.sin((u * Math.PI) / 2)], () => 1.5));
    chk(near(L(arc), 5 * Math.PI, 5e-3), `90° bend R 10 mm: length ${L(arc)} vs ${5 * Math.PI}`);
    const ub = ext(C.tubeMesh((u) => [6 - 6 * Math.cos(u * Math.PI), 0, 6 * Math.sin(u * Math.PI)], () => 1.2, { nu: 240 }));
    chk(near(L(ub), 6 * Math.PI, 5e-3), `U-bend (both caps in one plane): length ${L(ub)} vs ${6 * Math.PI}`);
    const hw = ext(C.tubeMesh(st, () => 1.5, { wall: 0.3 }));
    chk(near(hw.wallMm, 0.3 * poly, 0.02) && near(Dof(hw)[200], 3 * poly, 1e-3), `hollow wall: wall ${hw.wallMm}, lumen ${Dof(hw)[200]}`);
    const op = C.tubeMesh(st, () => 1.5), opn = ext(op.slice(0, op.length - 2 * 40 * 9));
    chk(opn.open && near(L(opn), 20, 1e-3), "open-ended surface (no caps)");
    const el = C.tubeMesh(st, () => 1, { nt: 120 });
    for (let i = 0; i < el.length; i += 3) el[i] *= 2;          // 2:1 elliptical lumen: rₕ⁴ = 2a³b³/(a²+b²)
    chk(near(ext(el).r[200] * 1e3, (16 / 5) ** 0.25, 2e-3), `elliptical lumen resistance radius ${ext(el).r[200] * 1e3}`);
    const stl = C.parseSTL(C.toBinarySTL(op));
    chk(stl.length === op.length && near(L(ext(stl)), 20, 1e-6), "binary STL round trip");
    chk(C.guessUnits(C.weld(C.tubeMesh((u) => [0, 0, 0.02 * u], () => 0.0015))) === "m", "metre-scale file guessed as metres");
    // round trip: parametric vessel → triangulated surface → import → same hemodynamics
    for (const extra of [{ occlusionPct: 70 }, { occlusionPct: 60, tortuosity: 1.2, waves: 1 }, { occlusionPct: 50, tortuosity: 1.1, waves: 2, taper: 0.8 }]) {
      const p = { ...M.DEFAULTS, ...extra, nBranches: 0, rigid: true }, r = M.solve(p), g = r.geom, n = g.s.length, nT = 64, tri = [];
      const v = (i, j) => { const th = (2 * Math.PI * j) / nT; return [0, 1, 2].map((k) => ([g.points.x, g.points.y, g.points.z][k][i] + r.radius[i] * (Math.cos(th) * g.N[i][k] + Math.sin(th) * g.B[i][k])) * 1e3); };
      const c = (i) => [g.points.x[i] * 1e3, g.points.y[i] * 1e3, g.points.z[i] * 1e3];
      for (let i = 0; i < n - 1; i++) for (let j = 0; j < nT; j++) tri.push(...v(i, j), ...v(i, (j + 1) % nT), ...v(i + 1, (j + 1) % nT), ...v(i, j), ...v(i + 1, (j + 1) % nT), ...v(i + 1, j));
      for (let j = 0; j < nT; j++) tri.push(...c(0), ...v(0, (j + 1) % nT), ...v(0, j), ...c(n - 1), ...v(n - 1, j), ...v(n - 1, (j + 1) % nT));
      const e = ext(Float32Array.from(tri));
      M.setCadGeometry(e);
      try {
        const b = M.solve({ ...p, geometry: "cad", lengthMm: L(e), diameterMm: 2e3 * e.r[0], taper: e.r[e.r.length - 1] / e.r[0], occlusionPct: 0 }).summary, a = r.summary;
        const tag = JSON.stringify(extra);
        chk(Math.abs(b.flow_fraction - a.flow_fraction) < 0.005, `${tag}: flow ${b.flow_fraction} vs ${a.flow_fraction}`);
        chk(near(b.dP_mmHg, a.dP_mmHg, 0.02), `${tag}: ΔP ${b.dP_mmHg} vs ${a.dP_mmHg}`);
        chk(near(b.wss_max, a.wss_max, 0.05), `${tag}: peak WSS ${b.wss_max} vs ${a.wss_max}`);
        chk(Math.abs(b.diameter_stenosis_pct - a.diameter_stenosis_pct) < 0.5, `${tag}: stenosis ${b.diameter_stenosis_pct}% vs ${a.diameter_stenosis_pct}%`);
      } finally { M.setCadGeometry(null); }
    }
    chk(!M.usesCad({ ...M.DEFAULTS, geometry: "cad" }), "cleared CAD geometry falls back to the parametric vessel");
  }

  // ── custom shape (shape editor) and devices acting on a lesion read from the geometry ──
  const sg = group("Custom shape and devices on drawn or imported lesions");
  {
    const chk = (ok, msg) => { sg.cases++; if (!ok) sg.fails.push(msg); };
    const P = { ...M.DEFAULTS, nBranches: 0 };
    const sl = M.solve({ ...P, occlusionPct: 60 }).summary, sh = M.solve({ ...P, shapeD: "0:1,0.35:1,0.5:0.4,0.65:1,1:1" }).summary;
    chk(Math.abs(sh.diameter_stenosis_pct - 60) < 0.1 && Math.abs(sh.flow_fraction - sl.flow_fraction) < 0.005,
      `drawn 60% narrowing ${sh.diameter_stenosis_pct}% / flow ${sh.flow_fraction} vs slider ${sl.flow_fraction}`);
    chk(JSON.stringify(M.solve({ ...P, shapeD: "0:1,1:1", shapeX: "0:0,1:0" }).summary) === JSON.stringify(M.solve(P).summary), "a flat shape equals no shape");
    for (const x of ["0:0,0.5:0.2,1:0", "0:0,0.25:0.12,0.75:-0.12,1:0"]) {
      const r = M.solve({ ...P, shapeX: x });
      chk(near(r.geom.s.at(-1), 0.02, 1e-12) && r.geom.s.at(-1) / r.geom.chord > 1.02 && !r.geom.lesion && r.summary.flow_fraction < 1,
        `drawn centerline ${x}: length ${r.geom.s.at(-1)}, tortuosity ${r.geom.s.at(-1) / r.geom.chord}`);
    }
    let threw = false;
    try { M.solve({ ...P, shapeX: "0:0,0.1:0.3,0.2:-0.3,0.3:0.3,0.4:-0.3,0.5:0.3,0.6:-0.3,1:0" }); } catch { threw = true; }
    chk(threw, "a centerline longer than the vessel must be rejected");
    const stent = DEVICES.find((d) => d.stent), ret = DEVICES.find((d) => d.category.startsWith("Stent retriever")), bal = DEVICES.find((d) => d.balloon);
    const p0 = { ...M.DEFAULTS, occlusionPct: 75 };
    chk(JSON.stringify(applyDevice(stent, p0, defaultOptions(stent)).post) === JSON.stringify(treat(stent, p0, defaultOptions(stent)).post), "slider lesion: applyDevice = treat");
    const ps = { ...P, shapeD: "0:1,0.35:1,0.5:0.25,0.65:1,1:1" };
    for (const dev of [stent, bal]) {
      const t = applyDevice(dev, ps, defaultOptions(dev)), b = M.solve(t.post).summary, a = M.solve(ps).summary;
      chk(Math.abs(t.lesionPre.occlusionPct - 75) < 0.1 && Math.abs(b.diameter_stenosis_pct - t.lesionPost.occlusionPct) < 0.1 && b.flow_fraction > a.flow_fraction,
        `${dev.name} on a drawn lesion: ${t.lesionPre.occlusionPct}% → ${b.diameter_stenosis_pct}% (planned ${t.lesionPost.occlusionPct}%)`);
    }
    const pt = { ...ps, lesion: "thrombus" };
    const ok = M.solve(applyDevice(ret, pt, { ...defaultOptions(ret), success: true }).post), bad = M.solve(applyDevice(ret, pt, { ...defaultOptions(ret), success: false }).post);
    chk(!ok.geom.lesion && near(ok.summary.flow_fraction, 1, 1e-6), `successful pass on a drawn clot leaves DS ${ok.summary.diameter_stenosis_pct}%`);
    chk(near(bad.summary.flow_fraction, M.solve(pt).summary.flow_fraction, 1e-12), "failed pass leaves the drawn clot");
    // imported design: the device treats the narrowing in the file
    const tri = C.tubeMesh((u) => [0, 0, 20 * u], (u) => 1.5 * (1 - 0.7 * (Math.abs(u - 0.5) < 0.15 ? 0.5 * (1 + Math.cos((Math.PI * (u - 0.5)) / 0.15)) : 0)), { nu: 300 });
    const e = C.extractLumen(C.weld(tri), { scale: 1e-3 });
    M.setCadGeometry(e);
    try {
      const pc = { ...P, geometry: "cad", lengthMm: 20, diameterMm: 2e3 * e.r[0], taper: 1 };
      const t = applyDevice(stent, pc, defaultOptions(stent)), b = M.solve(t.post).summary, a = M.solve(pc).summary;
      chk(Math.abs(t.lesionPre.occlusionPct - 70) < 1 && Math.abs(b.diameter_stenosis_pct - t.lesionPost.occlusionPct) < 0.2 && b.flow_fraction > a.flow_fraction + 0.2,
        `stent on an imported 70% narrowing: ${t.lesionPre.occlusionPct}% → ${b.diameter_stenosis_pct}%, flow ${a.flow_fraction} → ${b.flow_fraction}`);
    } finally { M.setCadGeometry(null); }
  }

  // ── engine: robustness across a parameter grid (both grids the app uses) ──
  const rob = group("Engine: converges across diameter × stenosis × length × branches × grid");
  for (const D of [1.5, 2, 2.5, 3, 4, 5]) for (const pre of [0, 50, 70, 85, 92, 95]) for (const L of [3, 4, 6, 10])
    for (const nb of [0, 6]) for (const n of [161, 401]) {
      rob.cases++;
      try { M.solve({ ...M.DEFAULTS, diameterMm: D, occlusionPct: pre, stenosisLengthMm: L, nBranches: nb }, n); }
      catch (e) { rob.fails.push(`D${D} ${pre}% L${L} branches${nb} n${n}: ${e.message}`); }
    }

  // ── devices: every device × size rule × outcome rule × flag ──
  const Ds = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5], pres = [50, 70, 85, 95], lens = [4, 6, 10, 16];
  for (const dev of DEVICES) {
    const g = group(`Device: ${dev.name}`);
    for (const D of Ds) for (const pre of pres) for (const len of lens) for (const success of (dev.balloon || dev.stent ? [true] : [true, false])) {
      g.cases++;
      const ctx = `D ${D} mm, ${pre}%, lesion ${len} mm${success ? "" : ", failed pass"}`;
      const fail = (m) => g.fails.push(`${ctx}: ${m}`);
      const p = { ...M.DEFAULTS, lesion: dev.lesion, diameterMm: D, occlusionPct: pre, stenosisLengthMm: len };
      const opts = { ...defaultOptions(dev), success };
      try {
        const t = treat(dev, p, opts);
        const errs = M.validate(t.post);
        if (errs.length) fail("post-treatment inputs invalid: " + errs.join(" "));
        const a = M.solve(p, 161).summary, b = M.solve(t.post, 161).summary;
        const postDS = t.post.occlusionPct;
        if (postDS > pre + 1e-9) fail(`stenosis worsened ${pre}% → ${postDS}%`);
        if (b.flow_fraction < a.flow_fraction - 1e-6) fail(`flow fell ${a.flow_fraction} → ${b.flow_fraction}`);
        if (!t.checks.every((c) => ["ok", "warn", "bad"].includes(c.level) && c.ok === (c.level === "ok"))) fail("check flag inconsistent with its level");
        const bad = (re) => t.checks.some((c) => c.level === "bad" && re.test(c.text));
        if (dev.balloon) {
          const bd = t.size.diameterMm, exp = Math.min(pre, 100 * (1 - Math.min(bd * (1 - opts.recoilPct / 100), D) / D));
          if (Math.abs(postDS - exp) > 1e-9) fail(`residual ${postDS} ≠ rule ${exp}`);
          const allowed = dev.sizes.filter((s) => s.diameterMm <= 0.8 * D + 1e-9);
          const want = allowed.length ? allowed[allowed.length - 1].diameterMm : dev.sizes[0].diameterMm;
          if (bd !== want) fail(`auto-size ${bd} mm ≠ largest ≤ 80% of vessel (${want})`);
        } else if (dev.stent) {
          if (Math.abs(postDS - Math.min(pre, 31.9)) > 1e-9) fail(`residual ${postDS} ≠ 31.9%`);
          const over = t.size.diameterMm - D;
          if (dev.sizes.some((s) => s.diameterMm - D >= 0.5 - 1e-9 && s.diameterMm - D <= 1 + 1e-9) && !(over >= 0.5 - 1e-9 && over <= 1 + 1e-9))
            fail(`stent ${t.size.diameterMm} mm not 0.5–1.0 mm over vessel though one exists`);
          if (dev.sizes.some((s) => s.diameterMm === t.size.diameterMm && s.lengthMm >= len + 6) && t.size.lengthMm < len + 6) fail("stent too short though a longer one exists");
          if (t.size.lengthMm - len < 6 - 1e-9 && !bad(/beyond each lesion end/)) fail("< 3 mm landing zone not flagged");
          if ((D < 2 || D > 4.5) && !bad(/indication/)) fail("vessel outside 2.0–4.5 mm not flagged");
        } else {
          if (!success && JSON.stringify(t.post) !== JSON.stringify(p)) fail("failed pass changed the vessel");
          if (success && Math.abs(b.flow_fraction - 1) > 1e-6) fail(`successful pass gives ${b.flow_fraction}, not 1`);
          if (dev.aspiration) {
            const s = t.size, e = aspirationEstimate({ idMm: s.idMm, lengthCm: ASPIRATION_DEFAULTS.catheterLengthCm, vacuumKpa: ASPIRATION_DEFAULTS.vacuumKpa, mu: p.mu, rho: p.rho });
            const r = s.idMm / 2e3, q = ((Math.PI * r ** 4 * ASPIRATION_DEFAULTS.vacuumKpa * 1e3) / (8 * p.mu * ASPIRATION_DEFAULTS.catheterLengthCm / 100)) * 6e7;
            if (!near(e.flowMlMin, q, 1e-12)) fail("aspiration flow ≠ Hagen–Poiseuille");
            if (s.odMm && s.odMm >= D && !bad(/OD/)) fail("catheter OD ≥ vessel not flagged");
          } else {
            const inRange = (x) => (x.vesselMinMm == null || D >= x.vesselMinMm) && (x.vesselMaxMm == null || D <= x.vesselMaxMm);
            if (dev.sizes.some((x) => inRange(x) && x.diameterMm >= D) && !(inRange(t.size) && t.size.diameterMm >= D)) fail(`auto-size ${t.size.label} outside its labeled range though one fits`);
            if (!inRange(t.size) && !bad(/labeled range/)) fail("out-of-range vessel not flagged");
          }
        }
        const wrong = treat(dev, { ...p, lesion: dev.lesion === "icas" ? "thrombus" : "icas" }, opts);
        if (!wrong.checks.some((c) => c.level === "bad" && /lesion is set/.test(c.text))) fail("wrong lesion type not flagged");
      } catch (e) { fail("threw: " + e.message); }
    }
  }
  return results;
}
