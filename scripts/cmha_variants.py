"""Occluded and recanalized versions of a CMHA healthy-control MCA, prepared for CFD and verified.

Usage: python3 scripts/cmha_variants.py AHMU1218142

The lesions are imposed computationally on the patient's own M1 lumen, using the same cosine
profile and the same device-outcome rules as the app (app/devices.js). The post-treatment lumens
are *prescribed* outcomes (published means or assumptions). No device–clot or device–wall
mechanics are simulated.

Writes CFD-ready STLs to data/cmha/prepared/<id>/variants/<variant>/ (git-ignored), a report to
docs/cfd-variants-<id>.md and a render to docs/cfd-variants-<id>.png.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from bvm import stlprep as sp  # noqa: E402
from bvm.hemodynamics import YOUNG_TSAI_KT  # noqa: E402
from cmha_prepare import axis_intersection  # noqa: E402

DATA = ROOT / "data" / "cmha"
WEB = ROOT / "app" / "data" / "patient"     # committed: compact geometry for the in-app viewer (CC BY)
MU, RHO = 0.0035, 1060.0           # this project's blood (the authors' case used μ = 0.00345)
MMHG = 133.322
GATEWAY_MM = [1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]   # app/devices.js
RECOIL, STENT_RESIDUAL = 0.15, 31.9


def m1_centerline(V, F):
    ports = sorted(sp.find_ports(V, F), key=lambda p: -p.diameter)
    B = axis_intersection([p.center for p in ports], [p.normal for p in ports])
    C, s = sp.resample(sp.march_centerline(V, ports[0], B), 0.25)
    A = sp.area_profile(V, F, C, s)
    D = 2 * np.sqrt(A / np.pi)
    # M1 ends where the section flares into the bifurcation (D > 1.4 × the median M1 diameter)
    med = np.median(D[2:len(D) // 2])
    end = next(i for i in range(4, len(D)) if D[i] > 1.4 * med)
    C, s = C[:end], s[:end]
    D = m1_diameters(V, F, C, s, ports[0])               # recompute on the trimmed centerline (end tangents)
    return C, s, D, ports


def m1_diameters(V, F, C, s, inlet):
    """Equivalent diameter along M1. Station 0 lies on the open port plane, where a slice is degenerate,
    so it takes the port's measured area (lesions never reach the inlet)."""
    A = sp.area_profile(V, F, C, s)
    A[0] = inlet.area
    return 2 * np.sqrt(A / np.pi)


def lesion_dp(s, D_healthy, D_lesion, q):
    """Reduced-order pressure drop along M1 (Poiseuille on the measured areas + Young–Tsai loss)."""
    r = D_lesion / 2e3
    visc = np.trapezoid(8 * MU * q / (np.pi * r**4), s * 1e-3)
    i = int(np.argmin(D_lesion / D_healthy))
    a0, a1 = np.pi * (D_healthy[i] / 2e3) ** 2, np.pi * (D_lesion[i] / 2e3) ** 2
    loss = YOUNG_TSAI_KT * RHO / 2 * (q / a0) ** 2 * (a0 / a1 - 1) ** 2 if a1 < 0.999 * a0 else 0.0
    return visc / MMHG, loss / MMHG, i


def territory_flow(s, Dh, Dl, q0=153 * 1e-6 / 60, p_in=85.0, p_v=10.0):
    """Distal flow when M1 drains into a resistance calibrated on the healthy lumen (app's territory BC)."""
    dp_h = sum(lesion_dp(s, Dh, Dh, q0)[:2])
    r_d = (p_in - p_v - dp_h) / q0
    lo, hi = q0 * 1e-6, q0 * 1.2
    for _ in range(80):
        q = 0.5 * (lo + hi)
        if sum(lesion_dp(s, Dh, Dl, q)[:2]) + r_d * q > p_in - p_v:
            hi = q
        else:
            lo = q
    return 0.5 * (lo + hi)


def main(sid: str) -> None:
    src = next((DATA / "controls" / "controls" / sid).glob("*.stl"))
    V0, F = sp.weld(sp.read_stl(src))
    C, s, D, ports = m1_centerline(V0, F)
    L1 = s[-1]
    center = 0.5 * L1
    i_c = int(np.argmin(abs(s - center)))
    d_local = float(np.median(D[max(0, i_c - 8):i_c + 9]))
    balloon = max(b for b in GATEWAY_MM if b <= 0.8 * d_local + 1e-9)
    angio_residual = 100 * (1 - min(balloon * (1 - RECOIL), d_local) / d_local)
    max_r = 1.6 * d_local / 2

    variants = [
        ("healthy", "Healthy baseline (original lumen)", None),
        ("occluded_icas70", "Occluded · ICAS plaque 70% (diameter), 6 mm", (70.0, 6.0)),
        ("occluded_thrombus95", "Occluded · embolic thrombus 95% (diameter), 8 mm", (95.0, 8.0)),
        ("recanalized_thrombectomy", "Recanalized · successful thrombectomy (0% residual)", (0.0, 8.0)),
        ("recanalized_angioplasty", f"Recanalized · Gateway {balloon:g} mm, 15% recoil → {angio_residual:.1f}% residual", (angio_residual, 6.0)),
        ("recanalized_stent", f"Recanalized · Wingspan, {STENT_RESIDUAL}% residual (HDE mean)", (STENT_RESIDUAL, 6.0)),
    ]
    q0 = 153e-6 / 60
    T = np.gradient(C, s, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True)

    def plane(si):
        i = int(np.argmin(abs(s - si)))
        return {"s_mm": float(s[i]), "point_mm": C[i].round(4).tolist(), "normal": T[i].round(5).tolist()}
    # starting outlet resistances for territory runs: healthy split by Murray's D³ (matched the authors'
    # CFD across 43 controls), venous 10 mmHg, M1 loss ≈ the healthy reduced-order ΔP
    o1, o2 = ports[1], ports[2]
    share1 = o1.diameter**3 / (o1.diameter**3 + o2.diameter**3)
    rows, renders = [], []
    for key, label, lesion in variants:
        if lesion and lesion[0] > 0:
            V, moved = sp.apply_stenosis(V0, C, s, center, lesion[1], lesion[0], max_radius=max_r)
        else:
            V, moved = V0.copy(), np.zeros(len(V0), bool)
        out = DATA / "prepared" / sid / "variants" / key
        out.mkdir(parents=True, exist_ok=True)
        ext = sp.extend(V, F, 10, 5)
        rep = sp.verify(V, F, ext)
        stem = f"{sid}_{key}_ext_in10D_out5D"
        sp.write_stl_patches(out / f"{stem}_mm.stl", ext.V, ext.F, ext.tags, ext.names)
        sp.write_stl_patches(out / f"{stem}_m.stl", ext.V, ext.F, ext.tags, ext.names, scale=1e-3)
        for k, name in enumerate(ext.names):
            sel = ext.tags == k
            sp.write_stl_patches(out / f"{stem}_{name}_mm.stl", ext.V, ext.F[sel], np.zeros(sel.sum(), int), [name])
        # lesion-specific verification
        Dl = m1_diameters(V, F, C, s, ports[0])
        n0, n1 = sp.face_normals(V0, F), sp.face_normals(V, F)
        flipped = int(np.sum(np.einsum("ij,ij->i", n0, n1) <= 0))
        outside_unchanged = bool(np.array_equal(V[~moved], V0[~moved]))
        target = None if not lesion else (1 - lesion[0] / 100)
        i_t = int(np.argmin(Dl / D))
        achieved = float(Dl[i_t] / D[i_t])
        qd = territory_flow(s, D, Dl)
        dpv, dpl, _ = lesion_dp(s, D, Dl, qd)            # at the flow the territory BC actually reaches
        dp_fixed = sum(lesion_dp(s, D, Dl, q0)[:2])
        span = (center - lesion[1] / 2, center + lesion[1] / 2) if lesion else (center - 3.0, center + 3.0)
        setup = {
            "subject": sid, "variant": key, "description": label, "units": "mm in *_mm.stl, m in *_m.stl",
            "patches": ext.names,
            "inlet": {"center_mm": ext.ports[0].center.round(4).tolist(), "inward_normal": (-ext.ports[0].normal).round(5).tolist(),
                      "radius_mm": ext.ports[0].diameter / 2,
                      "extension_end_center_mm": ext.end_centers[0].round(4).tolist(),
                      "note": "the inlet face is the cap at extension_end_center; its normal equals inward_normal"},
            "outlets": {n: {"diameter_mm": pp.diameter, "end_center_mm": e.round(4).tolist()}
                        for n, pp, e in zip(ext.names[2:], ext.ports[1:], ext.end_centers[1:])},
            "fluid": {"density_kg_m3": RHO, "viscosity_Pa_s": MU},
            "lesion": None if not lesion else {"percent_diameter": lesion[0], "length_mm": lesion[1],
                                                "throat_diameter_mm": float(Dl[i_t]), "throat_plane": plane(float(s[i_t]))},
            "measurement_planes": {"upstream": plane(span[0]), "downstream": plane(span[1]),
                                   "note": "area-average static pressure on these planes gives ΔP across the lesion (or the same segment when healthy)"},
            "territory_bc_initial": {"inlet_total_pressure_mmHg": 85, "venous_mmHg": 10, "healthy_inflow_mL_min": 153,
                                     "outlet_1_share": share1,
                                     "R_outlet_1_Pa_s_m3": (85 - 10 - 0.9) * MMHG / (q0 * share1),
                                     "R_outlet_2_Pa_s_m3": (85 - 10 - 0.9) * MMHG / (q0 * (1 - share1)),
                                     "note": "starting values; recalibrate from the healthy CFD run (docs/cfd-howto.md, step 3)"},
            "predictions_reduced_order": {},
        }
        export_web(sid, key, ext, len(F), C, s, D, Dl, qd, max_r, (center, lesion[1]) if lesion else None)
        rows.append(dict(key=key, label=label, rep=rep, setup=setup, out=out, flipped=flipped, outside=outside_unchanged,
                         moved=int(moved.sum()), target=target, achieved=achieved, s_t=float(s[i_t]),
                         d_t=float(Dl[i_t]), dpv=dpv, dpl=dpl, qd=qd, dp_fixed=dp_fixed, files=f"variants/{key}/{stem}_mm.stl"))
        renders.append((label, ext, len(F)))
        print(f"{key:26s} checks {'pass' if all(rep['checks'].values()) and not flipped and outside_unchanged else 'FAIL'}"
              f"  throat {achieved:.3f} (target {target if target is None else round(target, 3)})  "
              f"ΔP {dpv + dpl:.2f} mmHg at {qd * 6e7:.1f} mL/min  distal {100 * qd / q0:.1f}%")

    write_web_manifest(sid, rows, q0, L1, center, d_local)
    for r in rows:
        r["setup"]["predictions_reduced_order"] = {
            "distal_flow_mL_min_territory": r["qd"] * 6e7, "distal_flow_pct_of_healthy": 100 * r["qd"] / q0,
            "dP_M1_mmHg_at_that_flow": r["dpv"] + r["dpl"], "dP_M1_mmHg_if_153_forced": r["dp_fixed"]}
        (r["out"] / "cfd_setup.json").write_text(json.dumps(r["setup"], indent=2))
    # report
    ok = lambda b: "✅" if b else "❌"
    lines = [
        f"# Occluded and recanalized CFD geometries: {sid}",
        "",
        f"Built by `python3 scripts/cmha_variants.py {sid}` from the healthy control "
        "([cfd-geometry-" + sid + ".md](cfd-geometry-" + sid + ".md)). The STLs are in "
        f"`data/cmha/prepared/{sid}/variants/<variant>/` (git-ignored). Every variant has the same patches "
        "(`wall`, `inlet`, `outlet_1`, `outlet_2`) and extensions (inlet 10 D, outlets 5 D), so the runs differ only in the lesion.",
        "",
        f"![variants](cfd-variants-{sid}.png)",
        "",
        "## How the lesions are made",
        "",
        f"* **Where.** A centerline is marched through the patient's M1 from the inlet port. M1 ends where "
        f"the section flares into the bifurcation, giving {L1:.1f} mm along the curved centerline. Lesions "
        f"are centred mid-M1 (s = {center:.1f} mm), where the local healthy diameter is {d_local:.2f} mm.",
        "* **Shape.** The same cosine profile as the reduced-order model: each surface vertex near the lesion "
        "moves toward its nearest centerline point by 1 − depth·½(1 + cos πξ). The % is relative to the "
        "patient's *local* diameter, so the lumen keeps its own irregular shape, scaled down.",
        "* **What is prescribed, not simulated.** These are lesion and post-treatment *lumens* only, the "
        "same assumptions as the app's device panel:",
        "  * thrombus 95% over 8 mm (the app's embolic M1 preset)",
        "  * successful thrombectomy restores the original lumen (residual 0%)",
        f"  * Gateway sized to ≤ 80% of the vessel ({balloon:g} mm) with 15% assumed recoil → {angio_residual:.1f}% residual",
        f"  * Wingspan residual {STENT_RESIDUAL}% (HDE study mean)",
        "",
        "  Device struts, clot fragments, dissection and wall injury are not represented. The 3D CFD will "
        "show the hemodynamics of these lumens, not how a device produces them.",
        "* **Not fully occluded.** A 100% occlusion has no flow path, so CFD would return zero flow everywhere. "
        "Without collaterals, which this truncated geometry lacks, it is uninformative. 95% is the "
        "near-occlusion used instead.",
        "",
        "## Verification",
        "",
        "| Variant | Watertight, manifold, oriented | Volume check | Extensions clear | No flipped faces | Unchanged outside lesion | Throat / local D (target) |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        c = r["rep"]["checks"]
        topo = all(c[k] for k in ("watertight (no boundary edges)", "manifold (no edge in > 2 faces)",
                                  "consistent orientation", "outward normals (volume > 0)"))
        tgt = "—" if r["target"] is None or r["target"] == 1 else f"{r['achieved']:.3f} ({r['target']:.3f}) {ok(abs(r['achieved'] - r['target']) < 0.03)}"
        lines.append(f"| {r['key']} | {ok(topo)} | {ok(c['volume = original + extensions (±2%)'])} | "
                     f"{ok(c['extensions clear of the vessel'] and c['extensions clear of each other'])} | "
                     f"{ok(r['flipped'] == 0)} | {ok(r['outside'])} ({r['moved']} vertices moved) | {tgt} |")
    lines += [
        "",
        "Throat ratios are measured by slicing the deformed surface perpendicular to the centerline. They "
        "are not the requested values, so they test the whole construction. Small differences come from "
        "the patient lumen's irregular cross-section and the slice spacing (0.25 mm).",
        "",
        "## Reduced-order predictions to check the CFD against",
        "",
        "Computed from the *measured* area profile of each STL (Poiseuille along M1 plus the Young–Tsai "
        "expansion loss, μ = 0.0035 Pa·s). The territory column uses the app's boundary condition: inlet "
        "85 mmHg, venous 10 mmHg and a distal resistance calibrated so the healthy lumen carries "
        "153 mL/min. Branches and the bifurcation are lumped into that resistance.",
        "",
        "| Variant | Throat D (mm) | Distal flow, territory BC | ΔP along M1 at that flow (viscous + expansion), mmHg | ΔP if 153 mL/min were forced (fixed-flow run) | Restored vs. occluded |",
        "|---|---|---|---|---|---|",
    ]
    base = {r["key"]: r for r in rows}
    for r in rows:
        rest = ""
        if r["key"].startswith("recanalized"):
            pre = base["occluded_thrombus95" if "thrombectomy" in r["key"] else "occluded_icas70"]
            gain = (r["qd"] - pre["qd"]) / (q0 - pre["qd"])
            rest = f"{100 * gain:.0f}% of the deficit vs. {pre['key']}"
        forced = f"{r['dp_fixed']:.1f} mmHg" if r["dp_fixed"] < 75 else f"{r['dp_fixed']:.3g} mmHg: exceeds the 75 mmHg available, so do not run fixed-flow"
        dt = "—" if not r["moved"] else f"{r['d_t']:.2f}"
        lines.append(f"| {r['key']} | {dt} | {r['qd'] * 6e7:.1f} mL/min ({100 * r['qd'] / q0:.0f}%) | "
                     f"{r['dpv'] + r['dpl']:.2f} ({r['dpv']:.2f} + {r['dpl']:.2f}) | {forced} | {rest} |")
    lines += [
        "",
        "When the 3D results come in, compare ΔP across the lesion and the distal flow with this table. "
        "Agreement on the healthy and moderate cases, together with a quantified deviation at 95%, where "
        "3D jet and separation losses dominate, is the model-form evidence the credibility plan asks for. "
        "Because both are models, this is a code-to-code comparison; bench data from the vessel-mimic "
        "flow loop is still the validation.",
        "",
    ]
    (ROOT / "docs" / f"cfd-variants-{sid}.md").write_text("\n".join(lines))
    render(renders, ROOT / "docs" / f"cfd-variants-{sid}.png", sid)


def export_web(sid, key, ext, n_orig, C, s, D, Dl, qd, max_r, lesion_span=None):
    """Binary STL (faces grouped by patch) + per-face M1 station index, for the in-app viewer."""
    import struct
    from scipy.spatial import cKDTree
    out = WEB / sid
    out.mkdir(parents=True, exist_ok=True)
    order = np.argsort(ext.tags, kind="stable")          # wall (original faces first, then extensions), caps
    F = ext.F[order]
    tri = ext.V[F].astype("<f4")
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-30)
    rec = np.zeros(len(F), dtype=[("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")])
    rec["n"], rec["v"] = nrm, tri.reshape(-1, 9)
    head = f"BVM {sid} {key} | CMHA, Song et al. 2024 Sci Data, CC BY | mm".encode()[:80].ljust(80, b" ")
    (out / f"{key}.stl").write_bytes(head + struct.pack("<I", len(F)) + rec.tobytes())
    counts = {n: int(np.sum(ext.tags == k)) for k, n in enumerate(ext.names)}
    # M1 station for each original wall face (extensions and branches: -1)
    cen = ext.V[F[:counts["wall"]]].mean(1)
    dist, idx = cKDTree(C).query(cen)
    is_orig = np.zeros(counts["wall"], bool)
    is_orig[:n_orig] = True                               # stable sort kept original faces first
    inside = is_orig & (dist < max_r * 1.15) & (idx > 0) & (idx < len(s) - 1)
    face_s = np.where(inside, idx, -1).astype(int)
    r = Dl / 2e3
    # pressure along M1 relative to the inlet: Poiseuille on the measured areas + Young–Tsai loss ramped over
    # the diverging half of the lesion (same construction as the reduced-order engine)
    g = 8 * MU * qd / (np.pi * r**4)
    dp = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s) * 1e-3)])
    _, loss, i_t = lesion_dp(s, D, Dl, qd)
    if lesion_span and loss > 0:
        s_end = lesion_span[0] + lesion_span[1] / 2
        dp = dp + loss * MMHG * np.clip((s - s[i_t]) / max(s_end - s[i_t], 1e-9), 0, 1)
    # full precision: the in-app model re-solves from these, and a 0.13 mm throat is sensitive (ΔP ∝ r⁻⁴)
    prof = {"s_mm": s.tolist(), "D_mm": Dl.tolist(), "D_healthy_mm": D.tolist(),
            "wss_Pa": np.round(4 * MU * qd / (np.pi * r**3), 3).tolist(), "p_mmHg": np.round(-dp / MMHG, 4).tolist(),
            "re": np.round(4 * RHO * qd / (np.pi * MU * Dl * 1e-3), 2).tolist(), "q_mL_min": qd * 6e7}
    meta = {"counts": counts, "patch_order": ext.names, "n_original_wall": int(n_orig),
            "face_station": face_s.tolist(), "profile": prof,
            # what the in-app patient model needs to re-solve this lumen for other blood properties
            "lesion_span_mm": list(lesion_span) if lesion_span else None,
            "constants": {"mu_Pa_s": MU, "rho_kg_m3": RHO, "young_tsai_kt": YOUNG_TSAI_KT, "q_ref_mL_min": 153.0,
                          "p_in_mmHg": 85.0, "p_venous_mmHg": 10.0}}
    (out / f"{key}.json").write_text(json.dumps(meta, separators=(",", ":")))


def write_web_manifest(sid, rows, q0, L1, center, d_local):
    man = {
        "subject": sid,
        "attribution": "Geometry: CMHA dataset, Song M et al., Sci Data 2024 (doi:10.1038/s41597-024-04056-8), "
                       "Figshare 10.6084/m9.figshare.26965450, CC BY 4.0. Modified: flow extensions added; "
                       "lesions and post-treatment lumens imposed computationally (BVM project).",
        "m1_length_mm": L1, "lesion_center_mm": center, "local_diameter_mm": d_local,
        "units": "mm",
        "variants": [{"key": r["key"], "label": r["label"], "distal_flow_mL_min": r["qd"] * 6e7,
                      "distal_flow_pct": 100 * r["qd"] / q0, "dP_M1_mmHg": r["dpv"] + r["dpl"],
                      "throat_diameter_mm": r["d_t"] if r["moved"] else None,
                      "throat_ratio": r["achieved"] if r["moved"] else None,
                      "checks_pass": bool(all(r["rep"]["checks"].values()) and r["flipped"] == 0 and r["outside"])}
                     for r in rows],
    }
    (WEB / sid / "manifest.json").write_text(json.dumps(man, indent=1))


def render(items, path, sid):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(15, 8.6), facecolor="white")
    for k, (label, ext, n0) in enumerate(items):
        ax = fig.add_subplot(2, 3, k + 1, projection="3d")
        keep = np.arange(len(ext.F)) < n0             # show the patient lumen only (extensions hidden)
        tri = ext.V[ext.F[keep]]
        ax.add_collection3d(Poly3DCollection(tri, facecolors="#c0392b", edgecolor="none", alpha=0.95, shade=True))
        P = tri.reshape(-1, 3)
        lo, hi = P.min(0), P.max(0)
        mid, span = (lo + hi) / 2, (hi - lo).max() / 2
        for setter, m in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), mid):
            setter(m - span, m + span)
        ax.view_init(elev=28, azim=-120)
        ax.set_axis_off()
        ax.set_title(label, fontsize=10)
    fig.suptitle(f"{sid}: patient lumen with imposed lesions and prescribed post-treatment lumens "
                 "(flow extensions hidden)", fontsize=12)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=120)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "AHMU1218142")
