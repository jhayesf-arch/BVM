// Device catalog + treatment models for the MCA simulator.
//
// Every number below carries a source. Values we could not verify are marked `assumed: true`
// and surfaced as editable inputs in the UI rather than hidden constants. Patent numbers are
// listed exactly as published on each manufacturer's virtual patent-marking page; where the
// page did not list the product, we say so instead of guessing.

const GP = (n) => `https://patents.google.com/patent/US${n.replace(/,/g, "")}`;
const INCH = 25.4;

const SRC = {
  strykerPatents: { label: "Stryker neurovascular patent-marking page", url: "https://www.stryker.com/us/en/about/patents/neurovascular-intervention-patents.html" },
  wingspanDfu: { label: "Wingspan Stent System with Gateway PTA Balloon Catheter, US Directions for Use (Stryker)", url: "https://www.stryker.com/content/dam/stryker/neurovascular/products/wingspan-stent-system/resources/Wingspan-Gateway_DFU_US_Mar19.pdf" },
  bose2007: { label: "Bose et al., Stroke 2007;38:1531–1537 (Wingspan HDE study)", url: "https://doi.org/10.1161/STROKEAHA.106.477711" },
  sizingRules: { label: "Wingspan series: balloon at 80% of normal vessel, stent 0.5–1.0 mm over normal vessel", url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC8119103" },
  ptaGuide: { label: "Endovascular Today PTA balloon device guide (Gateway pressures, diameters)", url: "https://evtoday.com/device-guide/us/pta-balloons" },
  trevoIfu: { label: "Trevo NXT ProVue IFU, cell coverage / minimum vessel (Stryker)", url: "https://www.stryker.com/content/dam/stryker/neurovascular/products/trevo-nxt/resources/NV00045358-01A%20-%20Cell%20coverage-min%20vessel%20IFU.pdf" },
  trevo510k: { label: "FDA 510(k) K223305, Trevo NXT ProVue", url: "https://fda.innolitics.com/device/K223305" },
  solitaire: { label: "Medtronic Solitaire X product page (sizes, vessel ranges)", url: "https://www.medtronic.com/en-us/healthcare-professionals/products/neurological/neurovascular/stent-retrievers/solitaire-x-revascularization-device.html" },
  embotrap: { label: "J&J MedTech EMBOTRAP III product page", url: "https://www.jnjmedtech.com/en-US/product/embotrap-iii-revascularization-device" },
  tigerProtocol: { label: "TIGER trial protocol (Tigertriever sizing table)", url: "https://cdn.clinicaltrials.gov/large-docs/49/NCT03474549/Prot_SAP_000.pdf" },
  tiger510k: { label: "FDA 510(k) K230429, Tigertriever 21/17/13", url: "https://510k.innolitics.com/device/K230429" },
  rapidPatents: { label: "Rapid Medical patent-marking page", url: "https://www.rapid-medical.com/patents" },
  red510k: { label: "FDA 510(k) K211654, RED 72 (0.072 in ID, 0.085 in OD)", url: "https://fda.innolitics.com/device/K211654" },
  redSilver: { label: "Penumbra RED 72 SILVER LABEL (vacuum −29.2 inHg / −98.88 kPa)", url: "https://www.penumbrainc.com/red-72-sendit" },
  penumbraPatents: { label: "Penumbra patent-marking page", url: "https://www.penumbrainc.com/patents/" },
  zoom510k: { label: "FDA 510(k) K211476, ZOOM 71", url: "https://fda.innolitics.com/device/K211476" },
  imperativePatents: { label: "Imperative Care patent-marking page", url: "https://www.imperativecare.com/patents" },
  neuraviFamily: { label: "Google Patents: Neuravi Ltd clot-retrieval family (US9301769B2)", url: "https://patents.google.com/patent/US9301769B2/en" },
};

const patents = (nums, status, source) => ({ nums, status, source });

export const DEVICES = [
  // ───────────── acute thrombus: stent retrievers ─────────────
  {
    id: "trevo-nxt", name: "Trevo NXT ProVue", maker: "Stryker", category: "Stent retriever", lesion: "thrombus",
    summary: "Self-expanding retriever deployed across the clot through a microcatheter, then withdrawn with the clot integrated in its cells.",
    sizes: [
      { label: "3 × 32 mm", diameterMm: 3, lengthMm: 32, vesselMinMm: 1.4 },
      { label: "4 × 28 mm", diameterMm: 4, lengthMm: 28, vesselMinMm: 1.4 },
      { label: "4 × 41 mm", diameterMm: 4, lengthMm: 41, vesselMinMm: null },
      { label: "6 × 37 mm", diameterMm: 6, lengthMm: 37, vesselMinMm: 2.5 },
    ],
    sizingNote: "IFU: size the retriever to the vessel diameter at the deployment site; minimum vessel ID 1.4 mm (3×32, 4×28) and 2.5 mm (6×37). FDA K223305 labels the 3 and 4 mm models for 1.5–2.5 mm vessels.",
    patents: patents(["8,357,179", "8,529,596", "8,795,317", "8,795,345", "9,044,263", "9,072,537", "9,700,331"], "marking",
      { ...SRC.strykerPatents, note: "Listed for Trevo ProVue / XP ProVue Retrievers; Trevo NXT is not named on the first page of the listing." }),
    sources: [SRC.trevoIfu, SRC.trevo510k],
  },
  {
    id: "solitaire-x", name: "Solitaire X", maker: "Medtronic", category: "Stent retriever", lesion: "thrombus",
    summary: "Self-expanding retriever. The most-studied platform in the landmark thrombectomy trials, often used as a primary benchmark.",
    sizes: [
      { label: "3 × 20 mm", diameterMm: 3, lengthMm: 20, vesselMinMm: 1.0, vesselMaxMm: 3.0 },
      { label: "4 × 20 mm", diameterMm: 4, lengthMm: 20, vesselMinMm: 2.0, vesselMaxMm: 4.0 },
      { label: "4 × 40 mm", diameterMm: 4, lengthMm: 40, vesselMinMm: 2.0, vesselMaxMm: 4.0 },
      { label: "6 × 24 mm", diameterMm: 6, lengthMm: 24, vesselMinMm: 2.0, vesselMaxMm: 5.5 },
      { label: "6 × 40 mm", diameterMm: 6, lengthMm: 40, vesselMinMm: 2.0, vesselMaxMm: 5.5 },
    ],
    sizingNote: "Recommended vessel diameter is the smallest diameter at the thrombus site (Medtronic portfolio table). Lengths are usable lengths.",
    patents: patents([], "not-found", { label: "Medtronic virtual patent-marking list", url: "https://www.medtronic.com/patents", note: "Medtronic's marking list is a PDF we could not parse; no Solitaire patents confirmed." }),
    sources: [SRC.solitaire],
  },
  {
    id: "embotrap-3", name: "EMBOTRAP III", maker: "Cerenovus (J&J MedTech)", category: "Stent retriever", lesion: "thrombus",
    summary: "Dual-layer retriever: an outer cage with large inlet windows lets clot enter a reception space around an inner channel.",
    sizes: [
      { label: "5 × 22 mm", diameterMm: 5, lengthMm: 22, vesselMinMm: 1.5, vesselMaxMm: 5.0 },
      { label: "5 × 37 mm", diameterMm: 5, lengthMm: 37, vesselMinMm: 1.5, vesselMaxMm: 5.0 },
      { label: "6.5 × 45 mm", diameterMm: 6.5, lengthMm: 45, vesselMinMm: 1.5, vesselMaxMm: 6.5 },
    ],
    sizingNote: "22 and 37 mm models are indicated for 1.5–5.0 mm vessels; the 45 mm model for 1.5–6.5 mm.",
    patents: patents(["9,301,769", "10,588,649"], "family",
      { ...SRC.neuraviFamily, note: "Neuravi Ltd (EMBOTRAP's developer) clot-retrieval patent family. J&J's marking page returned 404, so these are not confirmed as the product's marking list." }),
    sources: [SRC.embotrap],
  },
  {
    id: "tigertriever", name: "Tigertriever", maker: "Rapid Medical", category: "Stent retriever (adjustable)", lesion: "thrombus",
    summary: "Operator-adjustable braided mesh: the physician sets the expanded diameter, and with it the radial force against the clot and wall.",
    adjustable: true,
    sizes: [
      { label: "Tigertriever 13", diameterMm: 2.5, lengthMm: 20.5, vesselMinMm: null, vesselMaxMm: 2.5 },
      { label: "Tigertriever 17", diameterMm: 3.0, lengthMm: 23, vesselMinMm: 0.5, vesselMaxMm: 3.0 },
      { label: "Tigertriever 21", diameterMm: 6.0, lengthMm: 32, vesselMinMm: 1.5, vesselMaxMm: 6.0 },
    ],
    sizingNote: "TIGER protocol table: 17 for 0.5–3.0 mm, 21 for 1.5–6.0 mm vessels; 13 expands to 2.5 mm OD (labeled range not found). Choose by the smallest vessel at the thrombus site.",
    patents: patents(["10,828,142", "11,083,473", "11,179,171", "11,813,156", "11,925,368", "11,980,380", "12,016,581"], "marking",
      { ...SRC.rapidPatents, note: "Listed as the Tigertriever patent family (not split by size)." }),
    sources: [SRC.tigerProtocol, SRC.tiger510k],
  },
  // ───────────── acute thrombus: aspiration ─────────────
  {
    id: "red72", name: "RED 72", maker: "Penumbra", category: "Aspiration catheter", lesion: "thrombus",
    summary: "Large-bore reperfusion catheter advanced to the clot face; pump vacuum ingests or corks the clot (contact aspiration).",
    aspiration: true,
    sizes: [{ label: "0.072″ ID", idMm: 0.072 * INCH, odMm: 0.085 * INCH }],
    sizingNote: "FDA K211654: 0.072 in min ID, 0.085 in max OD. Catheter length is not in our sources, so it is an editable assumption.",
    patents: patents([], "not-found", { ...SRC.penumbraPatents, note: "Penumbra's marking page does not list RED 72 by name (it lists a 'THNDRRED72KIT' Thunderbolt kit)." }),
    sources: [SRC.red510k, SRC.redSilver],
  },
  {
    id: "zoom71", name: "ZOOM 71", maker: "Imperative Care", category: "Aspiration catheter", lesion: "thrombus",
    summary: "Large-bore aspiration catheter with an asymmetric (beveled) tip; a direct comparator to RED 72 for aspiration efficiency.",
    aspiration: true,
    sizes: [{ label: "0.071″ ID", idMm: 0.071 * INCH, odMm: null }],
    sizingNote: "The 0.071 in ID is implied by the name and the ZOOM family range (0.035–0.071 in, K202182). OD and length are not in our sources.",
    patents: patents(["11,123,090", "11,147,949", "11,207,497", "11,311,303", "11,766,539", "12,042,160", "12,343,479", "12,350,443"], "marking",
      { ...SRC.imperativePatents, note: "Listed for ZOOM 71 (under 'Asymmetric Aspiration')." }),
    sources: [SRC.zoom510k],
  },
  // ───────────── chronic ICAS ─────────────
  {
    id: "gateway", name: "Gateway PTA Balloon", maker: "Stryker", category: "Angioplasty balloon", lesion: "icas",
    summary: "Submaximal angioplasty: an undersized balloon dilates the plaque, limiting dissection risk, alone or as predilation before Wingspan.",
    balloon: true,
    sizes: [1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0].map((d) => ({ label: `${d.toFixed(2)} mm`, diameterMm: d })),
    sizingNote: "Sized to ~80% of the normal vessel diameter. Diameters 1.5–4.0 mm; nominal 6 atm, rated burst 12 atm (14 atm for 2.25–3.0 mm). Intermediate size steps are assumed.",
    patents: patents(["6,280,457", "6,423,085", "6,921,410", "6,860,893", "7,070,607"], "marking",
      { ...SRC.strykerPatents, note: "Seen in a search-index excerpt of Stryker's page; not visible in the first page load (the list is paginated)." }),
    sources: [SRC.wingspanDfu, SRC.sizingRules, SRC.ptaGuide],
  },
  {
    id: "wingspan", name: "Wingspan Stent System", maker: "Stryker", category: "Intracranial stent (self-expanding nitinol)", lesion: "icas",
    summary: "Self-expanding nitinol stent placed after Gateway predilation for symptomatic 70–99% intracranial atherosclerotic stenosis (HDE).",
    stent: true,
    sizes: [2.5, 3.0, 3.5, 4.0, 4.5].flatMap((d) => [9, 15, 20].map((l) => ({ label: `${d.toFixed(1)} × ${l} mm`, diameterMm: d, lengthMm: l }))),
    sizingNote: "DFU: normal vessel 2.00–4.50 mm, lesion ≤ 14 mm, stent extends ≥ 3 mm past both ends of the lesion. Stent 0.5–1.0 mm larger than the normal vessel; Gateway predilation at 80%.",
    patents: patents(["8,092,508"], "marking", { ...SRC.strykerPatents, note: "Listed for Wingspan Stent System (claims cover the endoprosthesis delivery system)." }),
    sources: [SRC.wingspanDfu, SRC.bose2007, SRC.sizingRules],
  },
];

export const deviceById = (id) => DEVICES.find((d) => d.id === id);
export const patentUrl = GP;

/** Default procedure inputs for a device (all editable in the panel). */
export function defaultOptions(dev) {
  if (!dev) return {};
  if (dev.aspiration) return { success: true, residualPct: 0, vacuumKpa: 98.9, catheterLengthCm: 132 };
  if (dev.balloon) return { sizeIdx: null, balloonRatio: 0.8, recoilPct: 15 };
  if (dev.stent) return { sizeIdx: null, residualPct: 31.9, isrPct: 0 };
  return { success: true, residualPct: 0, sizeIdx: null, adjustDiameterMm: null };
}

const fmt = (v, d = 2) => (v == null || !isFinite(v) ? "—" : (+v.toFixed(d)).toString());

function autoSize(dev, p) {
  const D = p.diameterMm, lesion = p.stenosisLengthMm;
  if (dev.balloon) {
    const target = 0.8 * D;
    let k = 0;
    dev.sizes.forEach((s, i) => { if (s.diameterMm <= target + 1e-9) k = i; });
    return k;
  }
  if (dev.stent) {
    const needLen = lesion + 6;
    const ok = dev.sizes.map((s, i) => ({ s, i })).filter(({ s }) => s.diameterMm >= D + 0.5 - 1e-9 && s.diameterMm <= D + 1.0 + 1e-9);
    const pool = ok.length ? ok : dev.sizes.map((s, i) => ({ s, i })).filter(({ s }) => s.diameterMm >= D);
    const byLen = pool.filter(({ s }) => s.lengthMm >= needLen);
    const pick = (byLen.length ? byLen : pool).sort((a, b) => a.s.diameterMm - b.s.diameterMm || a.s.lengthMm - b.s.lengthMm)[0];
    return pick ? pick.i : dev.sizes.length - 1;
  }
  if (dev.aspiration) return 0;
  // retrievers: smallest device whose diameter ≥ vessel and whose labeled range contains it, long enough for the clot
  const cand = dev.sizes.map((s, i) => ({ s, i })).filter(({ s }) =>
    (s.vesselMinMm == null || D >= s.vesselMinMm) && (s.vesselMaxMm == null || D <= s.vesselMaxMm) && s.diameterMm >= D);
  const longEnough = cand.filter(({ s }) => s.lengthMm >= lesion);
  const pick = (longEnough.length ? longEnough : cand)[0];
  return pick ? pick.i : 0;
}

/**
 * Apply a device to the current (untreated) configuration.
 * Returns { size, sizeIdx, post (params for the treated vessel), overlay, checks[], metrics[] }.
 */
export function treat(dev, p, opts) {
  const D = p.diameterMm, L = p.lengthMm, lesionLen = p.stenosisLengthMm, sc = p.stenosisCenter * L;
  const lesionStart = sc - lesionLen / 2, lesionEnd = sc + lesionLen / 2;
  const preDS = p.basis === "diameter" ? p.occlusionPct : 100 * (1 - Math.sqrt(1 - p.occlusionPct / 100));
  const sizeIdx = opts.sizeIdx ?? autoSize(dev, p);
  const size = dev.sizes[Math.min(sizeIdx, dev.sizes.length - 1)];
  const checks = [], metrics = [];
  const check = (ok, text, level = "warn") => checks.push({ ok, text, level: ok ? "ok" : level });
  let post = { ...p, rigidSpan: null }, overlay = null;

  if (dev.lesion !== p.lesion)
    check(false, dev.lesion === "thrombus"
      ? "This is a thrombectomy device; the lesion is set to atherosclerotic plaque."
      : "This device treats atherosclerotic stenosis; the lesion is set to acute clot.", "bad");

  if (dev.aspiration) {
    const rId = size.idMm / 2 / 1e3, len = opts.catheterLengthCm / 100, dp = opts.vacuumKpa * 1e3;
    const q = (Math.PI * rId ** 4 * dp) / (8 * p.mu * len);
    const v = q / (Math.PI * rId * rId), re = (p.rho * v * 2 * rId) / p.mu;
    const force = Math.PI * rId * rId * dp;
    metrics.push(
      ["Aspiration flow, uniform-bore estimate", `${fmt(q * 6e7, 0)} mL/min`, "Hagen–Poiseuille for blood through the full catheter length at the set vacuum. It ignores tubing, the clot and any proximal ID step, so treat it as a comparative index rather than a bench value."],
      ["Reynolds number in catheter", `${fmt(re, 0)} (${re < 2000 ? "laminar" : "not laminar; estimate invalid"})`],
      ["Tip force index (tip area × vacuum)", `${fmt(force, 3)} N  ·  ${fmt(force / 9.81e-3, 1)} gf`, "Penumbra's 'thrombus removal force' metric."],
      ["Catheter ID / vessel diameter", fmt(size.idMm / D, 2)],
      ["Catheter OD / vessel diameter", size.odMm ? fmt(size.odMm / D, 2) : "OD not in sources"],
    );
    if (size.odMm) check(size.odMm < D, `Catheter OD ${fmt(size.odMm)} mm ${size.odMm < D ? "fits" : "does not fit"} the ${fmt(D)} mm vessel${size.odMm < D && size.odMm / D > 0.85 ? " (tight; good for corking)" : ""}.`, "bad");
    check(re < 2000, `Catheter flow is ${re < 2000 ? "laminar, so the Poiseuille estimate applies" : "not laminar, so the Poiseuille estimate is invalid"}.`);
    overlay = { type: "catheter", tipMm: Math.max(lesionStart - 0.3, 0.5), odMm: size.odMm ?? size.idMm + 0.3, idMm: size.idMm };
  } else if (dev.balloon) {
    const bd = size.diameterMm, ratio = bd / D;
    const lumen = Math.min(bd * (1 - opts.recoilPct / 100), D);
    const postDS = Math.min(preDS, 100 * (1 - lumen / D));
    post = { ...post, occlusionPct: postDS, basis: "diameter" };
    const bl = [9, 15, 20].find((x) => x >= lesionLen) ?? 20;
    check(ratio <= 1.0 + 1e-9, `Balloon ${fmt(bd)} mm is ${fmt(ratio * 100, 0)}% of the ${fmt(D)} mm normal vessel (guideline ≈ 80%; submaximal).`, ratio > 1 ? "bad" : "warn");
    if (ratio <= 1) checks[checks.length - 1].level = Math.abs(ratio - 0.8) <= 0.07 ? "ok" : "warn";
    check(D >= 2.0 && D <= 4.5, `Normal vessel ${fmt(D)} mm is ${D >= 2 && D <= 4.5 ? "within" : "outside"} the 2.00–4.50 mm Wingspan/Gateway indication.`);
    check(lesionLen <= 14, `Lesion length ${fmt(lesionLen, 1)} mm ${lesionLen <= 14 ? "≤" : ">"} 14 mm limit.`);
    metrics.push(
      ["Balloon / normal vessel", `${fmt(ratio * 100, 0)}%`],
      ["Post-dilation lumen", `${fmt(lumen)} mm after ${fmt(opts.recoilPct, 0)}% recoil (recoil is assumed)`],
      ["Diameter stenosis", `${fmt(preDS, 0)}% → ${fmt(postDS, 0)}%`],
      ["Inflation", "nominal 6 atm · rated burst 12–14 atm"],
    );
    overlay = { type: "balloon", centerMm: sc, lengthMm: bl, diameterMm: bd };
  } else if (dev.stent) {
    const sd = size.diameterMm, sl = size.lengthMm;
    const a = Math.max(sc - sl / 2, 0), b = Math.min(sc + sl / 2, L);
    const postDS = Math.min(preDS, Math.max(opts.residualPct, opts.isrPct));
    post = { ...post, occlusionPct: postDS, basis: "diameter", rigidSpan: [a, b] };
    const over = sd - D, margin = Math.min(lesionStart - (sc - sl / 2), (sc + sl / 2) - lesionEnd);
    check(D >= 2.0 && D <= 4.5, `Normal vessel ${fmt(D)} mm is ${D >= 2 && D <= 4.5 ? "within" : "outside"} the 2.00–4.50 mm indication.`, "bad");
    check(preDS >= 70, `Pre-treatment stenosis ${fmt(preDS, 0)}% (by diameter) ${preDS >= 70 ? "meets" : "is below"} the 70–99% indication.`);
    check(lesionLen <= 14, `Lesion length ${fmt(lesionLen, 1)} mm ${lesionLen <= 14 ? "≤" : ">"} 14 mm limit.`, "bad");
    check(over >= 0.5 - 1e-9 && over <= 1.0 + 1e-9, `Stent ${fmt(sd, 1)} mm is ${fmt(over, 2)} mm over the normal vessel (guideline 0.5–1.0 mm).`);
    check(margin >= 3 - 1e-9, `Stent covers ${fmt(margin, 1)} mm beyond each lesion end (DFU: ≥ 3 mm).`, "bad");
    check(sc - sl / 2 >= 0 && sc + sl / 2 <= L, `Stent ${sc - sl / 2 >= 0 && sc + sl / 2 <= L ? "fits within" : "runs past the end of"} the ${fmt(L, 0)} mm modeled segment.`);
    metrics.push(
      ["Diameter stenosis", `${fmt(preDS, 0)}% → ${fmt(postDS, 0)}%${opts.isrPct > opts.residualPct ? " (restenosis)" : ""}`],
      ["Stented (non-distensible) span", `${fmt(a, 1)}–${fmt(b, 1)} mm`],
      ["Clinical reference", "HDE study: 74.9% → 31.9% after stenting; 28% at 6 months"],
    );
    overlay = { type: "stent", startMm: a, endMm: b, diameterMm: sd };
  } else {
    // stent retriever
    const dd = dev.adjustable && opts.adjustDiameterMm ? Math.min(opts.adjustDiameterMm, size.diameterMm) : size.diameterMm;
    const inRange = (size.vesselMinMm == null || D >= size.vesselMinMm) && (size.vesselMaxMm == null || D <= size.vesselMaxMm);
    const lo = size.vesselMinMm, hi = size.vesselMaxMm;
    const rangeTxt = lo != null && hi != null ? `${lo}–${hi} mm`
      : lo != null ? `≥ ${lo} mm (no published maximum in our sources)`
      : hi != null ? `≤ ${hi} mm (no published minimum in our sources)` : "not in our sources";
    check(inRange, `Vessel ${fmt(D)} mm vs labeled range ${rangeTxt}.`, "bad");
    check(dd >= D, `Expanded diameter ${fmt(dd, 1)} mm ${dd >= D ? "≥" : "<"} vessel: ${dd >= D ? "wall apposition expected" : "under-sized, poor wall contact"}.`);
    check(size.lengthMm >= lesionLen, `Length ${fmt(size.lengthMm, 0)} mm ${size.lengthMm >= lesionLen ? "spans" : "is shorter than"} the ${fmt(lesionLen, 1)} mm clot.`);
    metrics.push(
      ["Device / vessel diameter", `${fmt(dd / D, 2)} (${dd >= D ? "+" : ""}${fmt((dd / D - 1) * 100, 0)}% oversize)`],
      ["Clot coverage margin", `${fmt(size.lengthMm - lesionLen, 1)} mm`],
    );
    overlay = { type: "retriever", centerMm: sc, lengthMm: Math.min(size.lengthMm, L - 0.5), diameterMm: Math.min(dd, D) };
  }

  if (!dev.balloon && !dev.stent) {
    // thrombectomy outcome: the clot is removed and any underlying stenosis remains
    post = opts.success ? { ...post, occlusionPct: Math.min(preDS, opts.residualPct), basis: "diameter" } : { ...p };
  }
  return { size, sizeIdx, post, overlay, checks, metrics };
}
