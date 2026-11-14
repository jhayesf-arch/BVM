"""Command line: ``python -m bvm run ...`` and ``python -m bvm sweep ...``."""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

from .export import write_all
from .hemodynamics import solve
from .params import Blood, FlowConditions, Vessel


def _add_common(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("vessel (modulated inputs)")
    g.add_argument("--diameter-mm", type=float, default=3.0)
    g.add_argument("--tortuosity", type=float, default=1.0, help="arc/chord index, 1.0 = straight")
    g.add_argument("--occlusion-pct", type=float, default=0.0)
    g.add_argument("--stiffness-kpa", type=float, default=1000.0, help="'inf' for a rigid wall")
    g = p.add_argument_group("vessel (other)")
    g.add_argument("--length-mm", type=float, default=20.0)
    g.add_argument("--tortuosity-waves", type=float, default=1.0)
    g.add_argument("--occlusion-basis", choices=["diameter", "area"], default="diameter")
    g.add_argument("--stenosis-length-mm", type=float, default=6.0)
    g.add_argument("--wall-thickness-mm", type=float, default=0.3)
    g.add_argument("--stenosis-center", type=float, default=0.5, help="fraction of arc length")
    g.add_argument("--taper", type=float, default=1.0, help="outlet/inlet diameter (1.0 = uniform)")
    g.add_argument("--bulb-pct", type=float, default=0.0, help="peak local dilation, e.g. carotid bulb")
    g.add_argument("--bulb-center", type=float, default=0.5)
    g.add_argument("--bulb-length-mm", type=float, default=8.0)
    g.add_argument("--n-branches", type=int, default=0, help="side branches, e.g. 6 lenticulostriates off M1")
    g.add_argument("--branch-start", type=float, default=0.3)
    g.add_argument("--branch-end", type=float, default=0.8)
    g = p.add_argument_group("blood and boundary conditions")
    g.add_argument("--mu", type=float, default=0.0035)
    g.add_argument("--rho", type=float, default=1060.0)
    g.add_argument("--flow-ml-min", type=float, default=120.0)
    g.add_argument("--bc", choices=["fixed_flow", "territory"], default="fixed_flow")
    g.add_argument("--inlet-pressure-mmhg", type=float, default=85.0)
    g.add_argument("--outflow-pressure-mmhg", type=float, default=10.0)
    g.add_argument("--external-pressure-mmhg", type=float, default=0.0)
    g.add_argument("--inlet-drop-mmhg", type=float, default=0.0,
                   help="upstream pressure loss (hypotension / proximal disease); lowers territory flow")
    g.add_argument("--branch-flow-pct", type=float, default=1.9, help="healthy share of inflow taken by all branches")


def _build(a: argparse.Namespace) -> tuple[Vessel, Blood, FlowConditions]:
    vessel = Vessel(diameter_mm=a.diameter_mm, length_mm=a.length_mm, tortuosity=a.tortuosity,
                    tortuosity_waves=a.tortuosity_waves, occlusion_pct=a.occlusion_pct,
                    occlusion_basis=a.occlusion_basis, stenosis_length_mm=a.stenosis_length_mm,
                    stiffness_kpa=a.stiffness_kpa, wall_thickness_mm=a.wall_thickness_mm,
                    stenosis_center=a.stenosis_center, taper=a.taper, bulb_pct=a.bulb_pct,
                    bulb_center=a.bulb_center, bulb_length_mm=a.bulb_length_mm, n_branches=a.n_branches,
                    branch_start=a.branch_start, branch_end=a.branch_end)
    flow = FlowConditions(flow_ml_min=a.flow_ml_min, inlet_pressure_mmhg=a.inlet_pressure_mmhg,
                          outflow_pressure_mmhg=a.outflow_pressure_mmhg,
                          external_pressure_mmhg=a.external_pressure_mmhg, bc=a.bc,
                          inlet_drop_mmhg=a.inlet_drop_mmhg, branch_flow_pct=a.branch_flow_pct)
    return vessel, Blood(mu=a.mu, rho=a.rho), flow


def _print_summary(s: dict) -> None:
    r, f, p, w = s["reynolds"], s["flow"], s["pressure"], s["wss"]
    print(f"Q = {f['Q_mL_min']:.1f} mL/min   V_inlet = {f['mean_velocity_inlet_m_s']:.3f} m/s   "
          f"V_throat = {f['mean_velocity_throat_m_s']:.3f} m/s")
    print(f"Re_inlet = {r['Re_inlet']:.0f} ({r['regime_inlet']})   "
          f"Re_throat = {r['Re_throat']:.0f} ({r['regime_throat']})   [{r['criteria']}]")
    if r["post_stenotic_disturbance_possible"]:
        print("  ! post-stenotic flow disturbance possible despite laminar Re (stenotic jets transition early)")
    print(f"Dean_max = {r['Dean_max']:.0f}   Womersley = {r['Womersley_inlet']:.2f}   "
          f"entrance length = {r['entrance_length_mm']:.0f} mm")
    print(f"ΔP = {p['dP_vessel_mmHg']:.2f} mmHg   Pd/Pa = {p['Pd_over_Pa']:.3f}")
    print(f"WSS inlet = {w['wss_inlet_Pa']:.2f} Pa   max = {w['wss_max_Pa']:.2f} Pa")
    if "perfusion" in s:
        pf = s["perfusion"]
        print(f"Territory flow = {pf['territory_flow_fraction']:.0%} of reference ({pf['ischemia_band']})")
    for msg in s["warnings"]:
        print(f"  warning: {msg}")


def _parse_values(text: str) -> list[float]:
    if ":" in text:
        start, stop, step = (float(x) for x in text.split(":"))
        return list(np.arange(start, stop + step / 2, step))
    return [float(x) for x in text.split(",")]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="bvm", description="Parametric MCA vessel model for ischemia studies")
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="solve one configuration and write geometry + results")
    _add_common(run)
    run.add_argument("--out", type=Path, default=Path("results/run"))
    run.add_argument("--no-plots", action="store_true")
    sw = sub.add_parser("sweep", help="vary one input and tabulate outputs")
    _add_common(sw)
    sw.add_argument("--param", required=True,
                    choices=["diameter_mm", "tortuosity", "occlusion_pct", "stiffness_kpa"])
    sw.add_argument("--values", required=True, help="start:stop:step or comma list, e.g. 0:90:10")
    sw.add_argument("--out", type=Path, default=Path("results/sweep"))
    a = ap.parse_args(argv)

    vessel, blood, flow = _build(a)
    a.out.mkdir(parents=True, exist_ok=True)

    if a.cmd == "run":
        sol = solve(vessel, blood, flow)
        paths = write_all(sol, a.out)
        if not a.no_plots:
            from .plotting import plot_3d, plot_profiles
            plot_3d(sol, a.out / "vessel_wss_3d.png", "wss")
            plot_profiles(sol, a.out / "profiles.png")
        _print_summary(sol.summary())
        print(f"\nwrote {', '.join(p.name for p in paths.values())} to {a.out}/")
        return

    rows = []
    for val in _parse_values(a.values):
        s = solve(replace(vessel, **{a.param: val}), blood, flow).summary()
        rows.append({a.param: val, "Q_mL_min": s["flow"]["Q_mL_min"],
                     "Re_inlet": s["reynolds"]["Re_inlet"], "regime_inlet": s["reynolds"]["regime_inlet"],
                     "Re_throat": s["reynolds"]["Re_throat"], "regime_throat": s["reynolds"]["regime_throat"],
                     "post_stenotic_disturbance_possible": s["reynolds"]["post_stenotic_disturbance_possible"],
                     "Dean_max": s["reynolds"]["Dean_max"], "dP_vessel_mmHg": s["pressure"]["dP_vessel_mmHg"],
                     "Pd_over_Pa": s["pressure"]["Pd_over_Pa"], "wss_max_Pa": s["wss"]["wss_max_Pa"],
                     "max_strain": s["wall"]["max_circumferential_strain"],
                     "territory_flow_fraction": s.get("perfusion", {}).get("territory_flow_fraction", math.nan)})
    with open(a.out / f"sweep_{a.param}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    from .plotting import plot_sweep
    plot_sweep(rows, a.param, a.out / f"sweep_{a.param}.png")
    hdr = f"{a.param:>14} {'Q mL/min':>9} {'Re_in':>6} {'Re_thr':>7} {'regime_thr':>12} {'ΔP mmHg':>8} {'WSSmax':>7}"
    print(hdr)
    for r in rows:
        print(f"{r[a.param]:>14.4g} {r['Q_mL_min']:>9.1f} {r['Re_inlet']:>6.0f} {r['Re_throat']:>7.0f} "
              f"{r['regime_throat']:>12} {r['dP_vessel_mmHg']:>8.2f} {r['wss_max_Pa']:>7.2f}")
    print(f"\nwrote sweep_{a.param}.csv and .png to {a.out}/")


if __name__ == "__main__":
    main()
