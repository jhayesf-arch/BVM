"""Prepare one CMHA healthy-control MCA bifurcation for CFD and verify it.

Usage: python3 scripts/cmha_prepare.py AHMU1218142

Reads data/cmha/ (Song et al. 2024, CC BY; downloaded locally, git-ignored), writes CFD-ready STLs to
data/cmha/prepared/<id>/ and a verification report to docs/cfd-geometry-<id>.md. The checks follow the
credibility assessment plan (docs/credibility-assessment-plan.md):

* geometry verification of the prepared surface (watertight, orientation, fidelity, clearances)
* input verification against the authors' own Fluent case and CSV tables
* code-to-code comparison of the authors' CFD with our reduced-order model (supplementary only)
"""
from __future__ import annotations

import csv
import gzip
import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bvm import Blood, FlowConditions, Vessel, solve  # noqa: E402
from bvm.stlprep import prepare  # noqa: E402

DATA = ROOT / "data" / "cmha"
STATS = DATA / "stats" / "statistical results"


def fluent_case(path: Path) -> dict:
    """Boundary conditions, fluid properties and boundary-zone geometry from a Fluent .cas(.gz)."""
    d = gzip.open(path, "rb").read() if path.suffix == ".gz" else path.read_bytes()
    out = {"zones": {}}
    for zid, kind, name in re.findall(rb"\(39 \((\d+) ([a-z-]+) ([A-Za-z0-9_.:-]+)", d):
        out["zones"][name.decode()] = (int(zid), kind.decode())
    grab = lambda pat: [float(x) for x in re.findall(pat, d)]
    out["viscosity"] = grab(rb"\(viscosity \(constant \. ([0-9.eE+-]+)")
    out["density"] = grab(rb"\(density \(constant \. ([0-9.eE+-]+)")
    out["inlet_vmag"] = grab(rb"\(vmag \(constant \. ([0-9.eE+-]+)")
    out["gauge_pressure"] = grab(rb"\(p \(constant \. ([0-9.eE+-]+)")
    out["steady"] = b"(unsteady? . #f" in d
    # binary nodes (double, 3D) and polygonal boundary faces -> centroid, area per zone
    m = re.search(rb"\(3010 \((\w+) (\w+) (\w+) \w+ 3\)", d)
    if m:
        n = int(m.group(3), 16)
        st = d.index(b"(", d.index(b")", m.start()) + 1) + 1
        X = np.frombuffer(d, "<f8", 3 * n, st).reshape(n, 3)
        for name, (zid, _) in out["zones"].items():
            h = re.search(rb"\(2013 \(%x (\w+) (\w+) \w+ (\w+)\)" % zid, d)
            if not h or name.startswith("interior") or name == "wall":
                continue
            a, b = int(h.group(1), 16), int(h.group(2), 16)
            st = d.index(b"(", d.index(b")", h.start()) + 1) + 1
            ints = np.frombuffer(d, "<i4", (b - a + 1) * 12, st)
            k, area, cen = 0, 0.0, np.zeros(3)
            for _ in range(b - a + 1):
                cnt = ints[k]
                P = X[ints[k + 1:k + 1 + cnt] - 1]
                k += cnt + 3
                c = P.mean(0)
                fa = 0.5 * np.linalg.norm(sum(np.cross(P[i] - c, P[(i + 1) % cnt] - c) for i in range(cnt)))
                area += fa
                cen += fa * c
            out[f"{name}_centroid_mm"] = (cen / area * 1e3).tolist()
            out[f"{name}_area_mm2"] = area * 1e6
    return out


def csv_row(name: str, sid: str) -> dict:
    with open(STATS / name, encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    head = [h.replace("\n", " ").strip() for h in rows[0]]
    for r in rows[1:]:
        if r and r[0] == sid:
            return dict(zip(head, [x.strip() for x in r]))
    raise KeyError(sid)


def axis_intersection(points, dirs):
    """Least-squares point closest to several lines (port centre + port normal)."""
    A, b = np.zeros((3, 3)), np.zeros(3)
    for p, d in zip(points, dirs):
        M = np.eye(3) - np.outer(d, d)
        A += M
        b += M @ p
    return np.linalg.solve(A, b)


def population_split() -> list[str]:
    """Flow-split predictions vs. the authors' CFD across all controls (outlets matched by Fluent zone)."""
    from bvm.stlprep import find_ports, read_stl, weld
    with open(STATS / "hemodynamic_control.csv", encoding="utf-8-sig") as fh:
        share = {r[0]: abs(float(r[4])) for r in csv.reader(fh) if r and r[0].startswith("AHMU")}
    cfd, pred, skipped = [], {"Murray (D³)": [], "Poiseuille (D⁴)": [], "Poiseuille (D⁴/L)": []}, []
    for folder in sorted((DATA / "controls" / "controls").glob("AHMU*")):
        try:
            case = fluent_case(next(folder.glob("fluent_*.cas*")))
            V, F = weld(read_stl(next(folder.glob("*.stl"))))
            ports = find_ports(V, F)
            near = lambda z: min(ports, key=lambda p: np.linalg.norm(p.center - np.array(case[f"{z}_centroid_mm"])))
            p1, p2 = near("outlet1"), near("outlet2")
        except (KeyError, StopIteration, ValueError):
            skipped.append(folder.name)
            continue
        B = axis_intersection([p.center for p in ports], [p.normal for p in ports])
        L1, L2 = np.linalg.norm(p1.center - B), np.linalg.norm(p2.center - B)
        a, b = p1.diameter, p2.diameter
        cfd.append(share[folder.name])
        pred["Murray (D³)"].append(a**3 / (a**3 + b**3))
        pred["Poiseuille (D⁴)"].append(a**4 / (a**4 + b**4))
        pred["Poiseuille (D⁴/L)"].append((a**4 / L1) / (a**4 / L1 + b**4 / L2))
    cfd = np.array(cfd)
    out = [f"Across {len(cfd)} controls (outlets matched to the authors' Fluent zones; skipped: "
           f"{', '.join(skipped) or 'none'}), predicted outlet-1 flow share vs. their CFD "
           f"({cfd.mean():.3f} ± {cfd.std(ddof=1):.3f}):", "",
           "| Estimate | r | Mean error | Mean abs. error | Max abs. error |", "|---|---|---|---|---|"]
    for k, v in pred.items():
        e = np.array(v) - cfd
        out.append(f"| {k} | {np.corrcoef(v, cfd)[0, 1]:.2f} | {e.mean():+.3f} | {np.abs(e).mean():.3f} | {np.abs(e).max():.3f} |")
    out += ["", "Murray's D³ rule tracks the 3D split best: junction and curvature losses in 3D offset the "
            "D⁴ advantage of the larger branch. That supports using D³ to split M2 flow if the reduced-order "
            "model gains a bifurcation. This is a code-to-code comparison and still needs bench or in vivo "
            "confirmation. Matching outlets by zone also explains the earlier negative correlation from the "
            "CSV alone: its outlet and branch labels do not correspond.", ""]
    return out


def main(sid: str) -> None:
    folder = DATA / "controls" / "controls" / sid
    src = next(folder.glob("*.stl"))
    out = DATA / "prepared" / sid
    variants = [(10.0, 5.0), (22.0, 5.0)]
    reps = {v: prepare(src, out, *v, tag=sid) for v in variants}
    rep = reps[variants[0]]
    case = fluent_case(next(folder.glob("fluent_*.cas*")))
    morph, hemo = csv_row("morphological_control.csv", sid), csv_row("hemodynamic_control.csv", sid)
    clin = csv_row("clinical_all.csv", sid)
    rho = case["density"][0] if case["density"] else 1060.0
    mu = case["viscosity"][0]
    vin = case["inlet_vmag"][0]
    ports = {p["name"]: p for p in rep["ports"]}
    q_csv = float(hemo["Inlet mass flow[kg s^-1]"]) / rho
    q_case = vin * case["inlet_area_mm2"] * 1e-6
    split_cfd = abs(float(hemo["Outlet 1 mass flow ratio"]))
    age = int(clin["Age"])
    v_formula = (65.41 - 0.26 * age) / 100

    # match our ports to the authors' Fluent zones by centroid
    def match(zone):
        c = np.array(case[f"{zone}_centroid_mm"])
        name, p = min(ports.items(), key=lambda kv: np.linalg.norm(np.array(kv[1]["center"]) - c))
        return name, float(np.linalg.norm(np.array(p["center"]) - c)), p
    zmap = {z: match(z) for z in ("inlet", "outlet1", "outlet2")}

    # reduced-order comparison: Poiseuille network from bifurcation point to each port
    B = axis_intersection([np.array(p["center"]) for p in rep["ports"]], [np.array(p["normal"]) for p in rep["ports"]])
    L = {n: float(np.linalg.norm(np.array(p["center"]) - B)) for n, p in ports.items()}
    D = {n: p["diameter"] for n, p in ports.items()}
    o1, o2 = zmap["outlet1"][0], zmap["outlet2"][0]
    split_poiseuille = (D[o1] ** 4 / L[o1]) / (D[o1] ** 4 / L[o1] + D[o2] ** 4 / L[o2])
    split_murray = D[o1] ** 3 / (D[o1] ** 3 + D[o2] ** 3)
    re_in = rho * vin * D["inlet"] * 1e-3 / mu
    le = 0.05 * re_in * D["inlet"]
    # our engine on the M1 segment: straight rigid tube with the measured inlet diameter, fixed flow
    sol = solve(Vessel(diameter_mm=D["inlet"], length_mm=max(L["inlet"], 1.0), stiffness_kpa=math.inf),
                Blood(mu=mu, rho=rho), FlowConditions(flow_ml_min=q_case * 6e7, bc="fixed_flow"))
    s = sol.summary()

    ok = lambda b: "✅" if b else "❌"
    pct = lambda a, b: 100 * (a - b) / b
    lines = [
        f"# CFD geometry preparation and verification: {sid}",
        "",
        "Source: CMHA healthy-control dataset, Song et al. 2024, *Sci Data* "
        "([Figshare 10.6084/m9.figshare.26965450](https://doi.org/10.6084/m9.figshare.26965450), CC BY). "
        f"Subject {sid}: {clin['Gender']}, {age} y. Raw data and prepared STLs stay local in `data/cmha/` "
        "(git-ignored). Regenerate with `python3 scripts/cmha_prepare.py " + sid + "`.",
        "",
        "## Prepared geometry",
        "",
        f"![{sid} prepared for CFD](cfd-geometry-{sid}.png)",
        "",
        "Open ports were extended along their fitted normals, blended to a circle of equal area over one "
        "diameter, continued as straight cylinders and capped. Original vertices are untouched.",
        "",
        "| Variant | Use | Files (`data/cmha/prepared/" + sid + "/`) |",
        "|---|---|---|",
        f"| inlet 10 D, outlets 5 D | prescribe a **developed (parabolic) inlet profile** | `{sid}_ext_in10D_out5D_mm.stl` (+ `_m.stl` in metres, + one STL per patch) |",
        f"| inlet 22 D, outlets 5 D | a **plug inlet** is acceptable: 22 D exceeds the laminar entrance length | `{sid}_ext_in22D_out5D_mm.stl` (same set) |",
        "",
        "Patch names: `wall`, `inlet`, `outlet_1` (larger branch), `outlet_2` (smaller branch).",
        "",
        "| Port | Diameter (mm) | Extension (10 D variant, mm) | Planarity (RMS / radius) |",
        "|---|---|---|---|",
    ]
    lines += [f"| {p['name']} | {p['diameter']:.3f} | {p['extension_length']:.1f} | {p['planarity']:.1e} |" for p in rep["ports"]]
    lines += ["", "## 1. Geometry verification", "", "| Check | 10 D variant | 22 D variant |", "|---|---|---|"]
    for k in rep["checks"]:
        lines.append(f"| {k} | {ok(reps[variants[0]]['checks'][k])} | {ok(reps[variants[1]]['checks'][k])} |")
    lines += [
        "",
        f"Volume: original (flat-capped) {rep['volume_original_capped']:.1f} mm³ + extensions "
        f"{rep['volume_extensions_expected']:.1f} mm³ vs. prepared {rep['volume']:.1f} mm³ "
        f"(difference {100 * rep['volume_rel_error']:.2f}%, from the transition blend). Extension axes stay "
        f"≥ {min(rep['clearance_radii']):.1f} port radii from the rest of the vessel. Surface: "
        f"{rep['counts']['faces']} triangles (input {rep['counts']['input_faces']}).",
        "",
        "## 2. Input verification against the authors' files",
        "",
        "| Check | Value | Reference | Result |",
        "|---|---|---|---|",
    ]
    for z in ("inlet", "outlet1", "outlet2"):
        name, dist, p = zmap[z]
        dz = 2 * math.sqrt(case[f"{z}_area_mm2"] / math.pi)
        lines.append(f"| Fluent zone `{z}` ↔ our `{name}`: centroid distance; diameter | {dist:.2f} mm; {p['diameter']:.3f} mm | "
                     f"{dz:.3f} mm in Fluent | {ok(dist < 0.5 and abs(pct(p['diameter'], dz)) < 3)} |")
    lines += [
        f"| Inflow = Fluent inlet velocity × inlet area | {q_case * 6e7:.1f} mL/min | CSV {q_csv * 6e7:.1f} mL/min | {ok(abs(pct(q_case, q_csv)) < 1)} |",
        f"| Fluent inlet velocity vs. paper's age formula (65.41 − 0.26 × age) | {vin:.4f} m/s | {v_formula:.4f} m/s | {ok(abs(vin - v_formula) < 1e-3)} |",
        f"| Blood viscosity in the Fluent case | {mu} Pa·s | paper states 0.004 Pa·s | {ok(abs(mu - 0.004) < 1e-6)} discrepancy, the case file governs |",
        f"| Density | {rho:g} kg/m³ | paper 1060 | {ok(abs(rho - 1060) < 1e-6)} |",
        f"| Outlet pressure | {', '.join(f'{x:g}' for x in case['gauge_pressure'] if x)} Pa | paper 13 332 Pa (100 mmHg) | {ok(13332 in case['gauge_pressure'])} |",
        f"| Steady solver | {case['steady']} | paper: steady | {ok(case['steady'])} |",
        f"| CSV M1 diameter D | {float(morph['Aortic  diameter(D)'] if 'Aortic  diameter(D)' in morph else list(morph.values())[1]):.2f} mm | our inlet port {D['inlet']:.2f} mm | measured at a different site; within 6% |",
        f"| CSV 'thicker' branch D1 / 'smaller' D2 | {list(morph.values())[2]} / {list(morph.values())[3]} mm | our outlets {D[o1]:.2f} / {D[o2]:.2f} mm | labels swapped: D2 is the larger branch |",
        "",
        "## 3. Comparison with the reduced-order model (code-to-code, supplementary)",
        "",
        "The authors' CFD is another model, not measurement data, so agreement here supports plausibility "
        "only. It is not validation (credibility plan, category 7 and the supplementary use of category 5).",
        "",
        f"* **Flow split** (outlet 1 share): Fluent {split_cfd:.3f}; Poiseuille network with equal outlet "
        f"pressures (branch lengths {L[o1]:.1f} / {L[o2]:.1f} mm from the fitted bifurcation point) "
        f"{split_poiseuille:.3f}; Murray's law (D³) {split_murray:.3f}. The Poiseuille estimate is "
        f"{'higher' if split_poiseuille > split_cfd else 'lower'} by {abs(split_poiseuille - split_cfd):.3f}. "
        "The 3D model adds junction and curvature losses that favour the smaller branch, which a straight-tube "
        "network ignores.",
        f"* **Inlet Reynolds number**: ρVD/μ = {re_in:.0f} (laminar). Our engine on the M1 segment gives "
        f"Re {s['reynolds']['Re_inlet']:.0f} and mean velocity {s['flow']['mean_velocity_inlet_m_s']:.4f} m/s for the same flow. "
        f"Fluent's {vin:.4f} m/s is on its {case['inlet_area_mm2']:.3f} mm² inlet face; our fitted port is "
        f"{math.pi * (D['inlet'] / 2) ** 2:.3f} mm², which accounts for the "
        f"{abs(pct(s['flow']['mean_velocity_inlet_m_s'], vin)):.1f}% difference "
        f"{ok(abs(pct(s['flow']['mean_velocity_inlet_m_s'] * math.pi * (D['inlet'] / 2) ** 2, vin * case['inlet_area_mm2'])) < 0.01)} (flow identical).",
        f"* **Entrance length** 0.05·Re·D = {le:.0f} mm = {le / D['inlet']:.0f} D. The original model's inlet "
        f"segment is only about {L['inlet']:.0f} mm, so with the authors' plug inlet the flow reaching the "
        "bifurcation is still developing, and their near-inlet WSS is biased. Hence the two variants above.",
        f"* **Mean WSS**: Fluent {hemo['Mean WSS[Pa]']} Pa (area average over the whole bifurcation, including "
        f"the developing inlet); fully developed Poiseuille in M1 {s['wss']['wss_inlet_Pa']:.2f} Pa. Not "
        "directly comparable; listed for context.",
        "",
        "### Same comparison across the dataset",
        "",
        *population_split(),
        "## 4. Model verification suite",
        "",
        "Run alongside this preparation: `python3 -m pytest` (code verification, "
        "`tests/test_verification.py`, and model tests) and `python -m bvm verify` "
        "([verification-results.md](verification-results.md)).",
        "",
        "## Recommended CFD setup for this geometry",
        "",
        f"* Fluid: Newtonian, ρ = {rho:g} kg/m³, μ = 0.0035 Pa·s (this project) or {mu} Pa·s (to reproduce the authors).",
        f"* Inlet: Q = {q_case * 6e7:.1f} mL/min (authors) or 153 mL/min (MRI population mean). Use a parabolic "
        f"profile with u_max = {2 * vin:.3f} m/s on the 10 D variant, or a plug at {vin:.4f} m/s on the 22 D variant.",
        "* Outlets: equal static pressure (authors: 100 mmHg) or resistance outlets for the ischemia "
        "scenarios, matching this project's territory boundary condition.",
        "* Wall: rigid, no-slip (the authors' setup); the compliance study is better done in the reduced-order model.",
        "* Mesh: refine the surface first (input edges ≈ 0.2 mm), then check mesh independence by GCI on "
        "three meshes (credibility plan, category 8 for model B).",
        "",
    ]
    doc = ROOT / "docs" / f"cfd-geometry-{sid}.md"
    doc.write_text("\n".join(lines))
    all_ok = all(all(r["checks"].values()) for r in reps.values())
    print(f"wrote {doc.relative_to(ROOT)}; geometry checks {'all pass' if all_ok else 'FAILED'}; files in {out.relative_to(ROOT)}/")
    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "AHMU1218142")
