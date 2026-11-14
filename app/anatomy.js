// Anatomical sections most vulnerable to ischemic stroke, with geometry from published
// studies. Every value carries its basis: "measured" (taken from the cited study),
// "derived" (computed from measured values, formula given) or "assumed" (our choice, adjustable).
//
// The model is a single unbranched tube, so each section is represented by the one flow path
// that matters for it (e.g. CCA → bulb → ICA for the carotid bifurcation). Branches, collaterals
// and autoregulation are not modeled.

const pubmed = (q) => `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(q)}`;

export const REFS = {
  amin2015: { label: "Amin-Hanjani et al. 2015, J Cereb Blood Flow Metab 35:312 · PC-MRA flow in ~320 healthy adults",
    url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC4426749" },
  zurada2011: { label: "Żurada et al. 2011, Clin Anat · CTA of the M1 segment, n = 115",
    url: "https://doi.org/10.1002/ca.21051" },
  unair554: { label: "Neuroangiography patterns of the middle cerebral artery: 554 cerebral angiograms",
    url: "https://repository.unair.ac.id/96429/1/Neuroangiography%20patterns%20of%20the%20middle%20cerebral%20artery%20study%20of%20554%20cerebral%20angiography%20results.pdf" },
  mcaMicro2009: { label: "Microvascular anatomy of the MCA, 70 unfixed hemispheres (J Neurosurg 2009): cortical branch outer diameters",
    url: pubmed("microsurgical anatomy proximal segments middle cerebral artery 70 hemispheres") },
  marinkovic: { label: "Marinković et al. · lenticulostriate artery casts: 2–12 per side, 0.10–1.28 mm",
    url: pubmed("Marinkovic lenticulostriate arteries anatomy") },
  bouvy2016: { label: "Bouvy et al. 2016, NMR Biomed · 7T flow MRI of perforating arteries (basal ganglia 3.9–5.1 cm/s)",
    url: pubmed("Bouvy 7T perforating arteries blood flow velocity pulsatility") },
  thomas2005: { label: "Thomas et al. 2005, Stroke 36:2450 · carotid bifurcation geometry, 25 young vs 25 older adults",
    url: "https://doi.org/10.1161/01.STR.0000185679.62634.0a" },
  ccaUs: { label: "Healthy-adult carotid ultrasound: CCA lumen 6.1 ± 0.8 mm, ~1 cm below the bifurcation",
    url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC8311942" },
  krejza2006: { label: "Krejza et al. 2006, Stroke 37:1103 · CCA/ICA diameter vs sex and body size",
    url: pubmed("Krejza 2006 carotid artery diameter men women") },
  vitosevic2019: { label: "Vitošević et al. 2019, Turk Neurosurg · CTA basilar morphometry, n = 120",
    url: "https://turkishneurosurgery.org.tr/Articles/Details/2211" },
  regina: { label: "Intracranial vertebrobasilar morphometry on non-contrast MRA, n = 100 (Eur J Anat)",
    url: "https://eurjanat.com/articles/a-retrospective-study-of-intracranial-vertebrobasilar-system-morphometry-and-its-variations-using-non-contrast-magnetic-resonance-angiogram-mra/" },
  satapathy2018: { label: "Satapathy & Mohapatra 2018 · basilar artery in cadavers (Eastern India)",
    url: pubmed("Satapathy Mohapatra basilar artery morphometry") },
  smoker1986: { label: "Smoker et al. 1986, AJNR · dolichoectasia criterion: basilar diameter > 4.5 mm",
    url: pubmed("Smoker 1986 basilar artery dolichoectasia CT criteria") },
  shoemaker2025: { label: "Shoemaker et al. 2025, J Appl Physiol 138:693 · pressure wire under anesthesia (n = 5): radial 80 ± 18, ICA 70 ± 25, MCA 62 ± 29 mmHg",
    url: "https://doi.org/10.1152/japplphysiol.00825.2024" },
  murray: { label: "Murray's law (1926): flow ∝ diameter³ along a branching tree",
    url: pubmed("Murray 1926 physiological principle of minimum work vascular system") },
};

const M1_FLOW = 153; // mean of L 159 and R 146 mL/min (Amin-Hanjani 2015, all ages)
const M1_D = 2.39;   // mm, angiographic mean (554 angiograms)
const murray = (d) => +(M1_FLOW * (d / M1_D) ** 3).toFixed(1);
const lsaFlow = (dMm, vCmS) => +((vCmS / 100) * Math.PI * (dMm / 2e3) ** 2 * 6e7).toFixed(2);

// basis entry helpers: [kind, text, refKey?]
const M = (text, ref) => ["measured", text, ref];
const D = (text, ref) => ["derived", text, ref];
const A = (text, ref) => ["assumed", text, ref];

export const SECTIONS = [
  {
    id: "m1",
    name: "MCA · M1 segment",
    short: "MCA M1",
    why: "The most common site of ischemic stroke. It takes the largest share of ICA flow, so emboli from the heart or neck follow that stream and lodge in the proximal horizontal M1.",
    ranges: { diameterMm: [1.5, 5, 0.05], flowMlMin: [40, 250, 1], lengthMm: [10, 40, 1] },
    common: {
      p: { wallMm: 0.3, inletMmHg: 85, outflowMmHg: 10, taper: 1, bulbPct: 0, tortuosity: 1, waves: 1, stenosisLengthMm: 6, stenosisCenter: 0.5 },
      basis: {
        inletMmHg: A("Mean pressure at the M1 origin. Wire data show losses even proximally (MCA 62 ± 29 vs radial 80 ± 18 mmHg).", "shoemaker2025"),
        wallMm: A("Thin-wall estimate for a ~3 mm cerebral artery."),
      },
    },
    sources: [
      { id: "lab", label: "BVM lab baseline (your specification)",
        p: { diameterMm: 3.0, lengthMm: 20, flowMlMin: 120 },
        basis: { diameterMm: A("Lab specification"), lengthMm: A("Lab specification"), flowMlMin: A("Lab specification: Q ≈ 120 mL/min, a plausible initial assumption") } },
      { id: "zurada", label: "Żurada 2011 · CTA, n = 115",
        p: { diameterMm: 2.23, lengthMm: 15.6, flowMlMin: M1_FLOW },
        basis: { diameterMm: M("Mean internal diameter 2.23 mm", "zurada2011"), lengthMm: M("Mean length 15.62 mm", "zurada2011"),
          flowMlMin: M("MCA 159 ± 28 (L) / 146 ± 28 (R) mL/min → mean 153", "amin2015") } },
      { id: "angio554", label: "Angiography · n = 554",
        p: { diameterMm: 2.39, lengthMm: 15.6, flowMlMin: M1_FLOW },
        basis: { diameterMm: M("2.39 ± 0.49 mm", "unair554"), lengthMm: M("15.56 ± 7.75 mm", "unair554"),
          flowMlMin: M("MCA 159 / 146 mL/min (L/R) → 153", "amin2015") } },
      { id: "elderly", label: "Older adults (61–80 y) · flow",
        p: { diameterMm: 2.39, lengthMm: 15.6, flowMlMin: 140 },
        basis: { diameterMm: M("2.39 mm", "unair554"), lengthMm: M("15.56 mm", "unair554"),
          flowMlMin: M("MCA 146 (L) / 133 (R) mL/min at 61–80 y → 140", "amin2015") } },
    ],
    scenarios: [
      { name: "Embolic M1 occlusion · 95%", p: { lesion: "thrombus", occlusionPct: 95, stenosisLengthMm: 8 } },
      { name: "Symptomatic ICAS · 75%", p: { lesion: "icas", occlusionPct: 75 } },
      { name: "Severe ICAS · 85%", p: { lesion: "icas", occlusionPct: 85 } },
    ],
  },
  {
    id: "lsa",
    name: "Lateral lenticulostriate artery",
    short: "Lenticulostriate",
    why: "Small perforators leave M1 at near-right angles to supply the basal ganglia and internal capsule. They are end arteries with almost no collateral backup and see near-M1 pressure, so lipohyalinosis or a microatheroma at the origin causes lacunar infarcts.",
    ranges: { diameterMm: [0.1, 1.3, 0.01], flowMlMin: [0.02, 5, 0.01], lengthMm: [5, 30, 1] },
    common: {
      p: { inletMmHg: 85, outflowMmHg: 10, wallMm: 0.05, taper: 1, bulbPct: 0, tortuosity: 1, waves: 1.5, stenosisLengthMm: 4, stenosisCenter: 0.3 },
      basis: {
        inletMmHg: A("Arises directly from M1, so it sees near-M1 pressure. That is the 'high pressure on a tiny end artery' that drives lipohyalinosis."),
        wallMm: A("~50 µm wall for a 0.3–0.8 mm artery. Not from a cited measurement."),
        lengthMm: A("Course from the M1 origin into the basal ganglia; not taken from a cited measurement."),
      },
    },
    sources: [
      { id: "typical", label: "Typical lateral LSA · 0.4 mm",
        p: { diameterMm: 0.4, lengthMm: 20, flowMlMin: lsaFlow(0.4, 4.5) },
        basis: { diameterMm: M("Within the 0.10–1.28 mm range of cast studies; 0.4 mm is a representative pick", "marinkovic"),
          flowMlMin: D(`Mean velocity 4.5 cm/s × lumen area → ${lsaFlow(0.4, 4.5)} mL/min`, "bouvy2016") } },
      { id: "large", label: "Large lateral LSA · 0.8 mm",
        p: { diameterMm: 0.8, lengthMm: 20, flowMlMin: lsaFlow(0.8, 4.5) },
        basis: { diameterMm: M("Upper part of the 0.10–1.28 mm range", "marinkovic"),
          flowMlMin: D(`4.5 cm/s × lumen area → ${lsaFlow(0.8, 4.5)} mL/min`, "bouvy2016") } },
      { id: "small", label: "Small LSA · 0.2 mm",
        p: { diameterMm: 0.2, lengthMm: 15, flowMlMin: lsaFlow(0.2, 4.5) },
        basis: { diameterMm: M("Lower part of the 0.10–1.28 mm range", "marinkovic"),
          flowMlMin: D(`4.5 cm/s × lumen area → ${lsaFlow(0.2, 4.5)} mL/min`, "bouvy2016") } },
    ],
    scenarios: [
      { name: "Lipohyalinosis · diffuse 50% narrowing", p: { lesion: "icas", occlusionPct: 50, stenosisLengthMm: 16, stenosisCenter: 0.5 } },
      { name: "Microatheroma at the origin · 70%", p: { lesion: "icas", occlusionPct: 70, stenosisLengthMm: 3, stenosisCenter: 0.12 } },
      { name: "Branch occlusion · 90%", p: { lesion: "thrombus", occlusionPct: 90, stenosisLengthMm: 4, stenosisCenter: 0.3 } },
    ],
  },
  {
    id: "carotid",
    name: "ICA bifurcation & carotid bulb",
    short: "Carotid bulb",
    why: "Where the common carotid splits, flow separates in the widened bulb. The low, oscillating shear there makes the bulb (and the terminal ICA bifurcation) prime sites for atherosclerotic plaque and thrombus.",
    ranges: { diameterMm: [4, 9, 0.05], flowMlMin: [100, 500, 1], lengthMm: [20, 60, 1] },
    common: {
      p: { inletMmHg: 90, outflowMmHg: 10, wallMm: 0.6, tortuosity: 1, waves: 1, lengthMm: 40, bulbPct: 30, bulbCenter: 0.35, bulbLengthMm: 15,
        stenosisLengthMm: 10, stenosisCenter: 0.35, flowMlMin: 257 },
      basis: {
        flowMlMin: M("ICA 259 ± 50 (L) / 256 ± 52 (R) mL/min → 257. Only the ICA path is modeled; ECA flow leaves at the bifurcation.", "amin2015"),
        bulbPct: A("Bulb shape is illustrative: no cited normal bulb diameter. Adjust to a patient scan."),
        lengthMm: A("CCA → bulb → cervical ICA segment."),
        inletMmHg: A("Neck pressure, slightly above M1 inlet."),
        wallMm: A("~0.6 mm wall, on the order of carotid intima-media thickness."),
      },
    },
    sources: [
      { id: "young", label: "Thomas 2005 · young adults (24 ± 4 y)",
        p: { diameterMm: 6.1, taper: 0.81 },
        basis: { diameterMm: M("CCA lumen 6.1 ± 0.8 mm (ultrasound)", "ccaUs"), taper: M("ICA/CCA diameter ratio 0.81 ± 0.06", "thomas2005") } },
      { id: "older", label: "Thomas 2005 · older adults (63 ± 10 y)",
        p: { diameterMm: 6.1, taper: 0.77 },
        basis: { diameterMm: M("CCA lumen 6.1 ± 0.8 mm (ultrasound)", "ccaUs"), taper: M("ICA/CCA diameter ratio 0.77 ± 0.12", "thomas2005") } },
      { id: "male", label: "Larger CCA · 6.5 mm",
        p: { diameterMm: 6.5, taper: 0.81 },
        basis: { diameterMm: M("Krejza 2006 as cited: men 6.5 mm, women 6.1 mm", "krejza2006"), taper: M("ICA/CCA 0.81", "thomas2005") } },
    ],
    scenarios: [
      { name: "Bulb plaque · 50%", p: { lesion: "icas", occlusionPct: 50 } },
      { name: "Bulb plaque · 70% (severe)", p: { lesion: "icas", occlusionPct: 70 } },
      { name: "Critical bulb stenosis · 90%", p: { lesion: "icas", occlusionPct: 90, basis: "diameter" } },
    ],
  },
  {
    id: "watershed",
    name: "Watershed · distal cortical branch",
    short: "Watershed",
    why: "Border zones between the ACA, MCA and PCA territories are fed by the farthest, narrowest branches. They receive blood last, so they are the first to starve when systemic pressure falls or a proximal artery narrows.",
    ranges: { diameterMm: [0.5, 2.5, 0.05], flowMlMin: [1, 80, 0.1], lengthMm: [15, 60, 1] },
    common: {
      p: { inletMmHg: 65, outflowMmHg: 10, wallMm: 0.12, lengthMm: 40, taper: 0.6, bulbPct: 0, tortuosity: 1, waves: 3, stenosisLengthMm: 6, stenosisCenter: 0.5 },
      basis: {
        inletMmHg: A("Pressure already lost along proximal arteries. Wire data under anesthesia: MCA 62 ± 29 vs radial 80 ± 18 mmHg.", "shoemaker2025"),
        taper: A("Branch narrows toward the border zone. Pial arteries downstream are 50–300 µm."),
        lengthMm: A("Course from the sylvian fissure to the border zone."),
        wallMm: A("Thin cortical artery wall."),
      },
    },
    sources: [
      { id: "large", label: "Large cortical branch · 1.5 mm (angular)",
        p: { diameterMm: 1.5, flowMlMin: murray(1.5) },
        basis: { diameterMm: M("Angular artery outer diameter 1.5 mm (largest cortical branch)", "mcaMicro2009"),
          flowMlMin: D(`Murray scaling from M1 (153 mL/min at 2.39 mm): 153 × (1.5/2.39)³ = ${murray(1.5)} mL/min`, "murray") } },
      { id: "small", label: "Small cortical branch · 0.8 mm (temporopolar)",
        p: { diameterMm: 0.8, flowMlMin: murray(0.8) },
        basis: { diameterMm: M("Temporopolar artery outer diameter 0.8–0.9 mm (smallest cortical branch)", "mcaMicro2009"),
          flowMlMin: D(`Murray scaling: 153 × (0.8/2.39)³ = ${murray(0.8)} mL/min`, "murray") } },
    ],
    scenarios: [
      { name: "Tortuous course over the gyri · TI 1.2", p: { tortuosity: 1.2 } },
      { name: "Systemic hypotension · −25 mmHg", p: { inletDropMmHg: 25 } },
      { name: "Hypotension + proximal stenosis · −40 mmHg", p: { inletDropMmHg: 40 } },
      { name: "Distal branch embolus · 80%", p: { lesion: "thrombus", occlusionPct: 80, stenosisLengthMm: 5 } },
    ],
  },
  {
    id: "basilar",
    name: "Basilar artery · proximal & mid",
    short: "Basilar",
    why: "Formed by the two vertebral arteries, it feeds the brainstem. It is prone to progressive atherosclerosis, and plaque can cover the origins of its tiny paramedian perforators.",
    ranges: { diameterMm: [1.5, 6, 0.05], flowMlMin: [40, 300, 1], lengthMm: [15, 45, 1] },
    common: {
      p: { inletMmHg: 85, outflowMmHg: 10, wallMm: 0.3, taper: 1, bulbPct: 0, tortuosity: 1, waves: 1, stenosisLengthMm: 6, stenosisCenter: 0.4 },
      basis: {
        inletMmHg: A("Same as the M1 inlet assumption."),
        wallMm: A("Thin-wall estimate for a ~3 mm cerebral artery."),
      },
    },
    sources: [
      { id: "cta", label: "Vitošević 2019 · CTA, n = 120 (> 60 y)",
        p: { diameterMm: 3.17, lengthMm: 32.4, flowMlMin: 138 },
        basis: { diameterMm: M("3.17 mm in patients over 60 (2.87 mm under 60)", "vitosevic2019"),
          lengthMm: D("Mean of men 34.07 and women 30.79 mm", "vitosevic2019"), flowMlMin: M("BA 138 ± 41 mL/min", "amin2015") } },
      { id: "mra", label: "Non-contrast MRA · n = 100",
        p: { diameterMm: 3.06, lengthMm: 31.3, flowMlMin: 138 },
        basis: { diameterMm: M("3.01–3.12 mm", "regina"), lengthMm: M("31.31 mm", "regina"), flowMlMin: M("BA 138 ± 41 mL/min", "amin2015") } },
      { id: "cadaver", label: "Cadaver · Satapathy 2018",
        p: { diameterMm: 3.05, lengthMm: 25.6, flowMlMin: 124 },
        basis: { diameterMm: M("3.05 ± 0.41 mm", "satapathy2018"), lengthMm: M("25.58 ± 3.57 mm", "satapathy2018"),
          flowMlMin: M("BA 124 ± 31 mL/min at 61–80 y", "amin2015") } },
    ],
    scenarios: [
      { name: "Mid-basilar ICAS · 70%", p: { lesion: "icas", occlusionPct: 70, stenosisCenter: 0.5 } },
      { name: "Basilar thrombosis · 95%", p: { lesion: "thrombus", occlusionPct: 95, stenosisLengthMm: 10, stenosisCenter: 0.4 } },
      { name: "Dolichoectasia · 5 mm, tortuous", p: { diameterMm: 5.0, tortuosity: 1.3, waves: 1.5 },
        basis: { diameterMm: M("Diameter > 4.5 mm meets the dolichoectasia criterion", "smoker1986") } },
    ],
  },
];

export const sectionById = (id) => SECTIONS.find((s) => s.id === id) ?? SECTIONS[0];

/** Parameters for a section + source (geometry only; pathology inputs reset to healthy). */
export function sectionParams(sec, srcId) {
  const src = sec.sources.find((s) => s.id === srcId) ?? sec.sources[0];
  return {
    section: sec.id, occlusionPct: 0, inletDropMmHg: 0, rigidSpan: null,
    bulbCenter: 0.5, bulbLengthMm: 8, ...sec.common.p, ...src.p,
  };
}

/** Basis rows for display: key -> [kind, text, ref]. */
export function sectionBasis(sec, srcId) {
  const src = sec.sources.find((s) => s.id === srcId) ?? sec.sources[0];
  return { ...sec.common.basis, ...src.basis };
}
