// Browser self-test for the JS engine and device logic (verification evidence; see
// docs/credibility-assessment-plan.md). Open app/selftest.html locally or on the live site.
import * as M from "./model.js?v=202611151712";
import { ASPIRATION_DEFAULTS, DEVICES, aspirationEstimate, defaultOptions, treat } from "./devices.js?v=202611151712";

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
