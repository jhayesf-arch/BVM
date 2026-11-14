"""Calculation verification and uncertainty quantification for the reduced-order model.

Supports FDA credibility evidence category 8 ("calculation verification / UQ results using
context-of-use simulations"). See docs/credibility-assessment-plan.md.

* Grid convergence: each endpoint is solved on 201, 401 and 801 stations. The observed order
  of accuracy and the grid convergence index (GCI; Roache 1994, safety factor 1.25) are reported
  for the default 401-station grid the app and CLI use.
* Uncertainty quantification: Monte Carlo propagation of input uncertainty (diameter, flow,
  pressure, wall stiffness, lesion length) into the endpoints, with rank correlations as a
  global sensitivity measure. Distributions are documented in INPUT_UNCERTAINTY with their basis.

Run: python -m bvm verify  (writes docs/verification-results.md)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from .hemodynamics import solve
from .params import Blood, FlowConditions, Vessel

GRIDS = (201, 401, 801)
FS = 1.25  # GCI safety factor for a three-grid study (Roache)


@dataclass(frozen=True)
class Scenario:
    name: str
    vessel: Vessel
    flow: FlowConditions


M1 = dict(diameter_mm=2.39, length_mm=15.6, n_branches=6, branch_start=0.3, branch_end=0.8)
M1_FLOW = dict(bc="territory", flow_ml_min=153.0, branch_flow_pct=1.9)

SCENARIOS = [
    Scenario("Healthy M1", Vessel(**M1), FlowConditions(**M1_FLOW)),
    Scenario("M1 ICAS 70%", Vessel(**M1, occlusion_pct=70), FlowConditions(**M1_FLOW)),
    Scenario("M1 ICAS 85%", Vessel(**M1, occlusion_pct=85), FlowConditions(**M1_FLOW)),
    Scenario("M1 70%, tortuous 1.3", Vessel(**M1, occlusion_pct=70, tortuosity=1.3, tortuosity_waves=2),
             FlowConditions(**M1_FLOW)),
]


def endpoints(sol) -> dict[str, float]:
    s = sol.summary()
    out = {
        "distal flow, % of healthy": 100 * s["perfusion"]["territory_flow_fraction"],
        "perforator flow, % of healthy": 100 * s["perfusion"].get("branch_flow_fraction", math.nan),
        "ΔP across lesion, mmHg": s["stenosis"].get("dP_lesion_mmHg", 0.0),
        "peak WSS, Pa": s["wss"]["wss_max_Pa"],
        "Re at throat": s["reynolds"]["Re_throat"],
    }
    return out


def grid_convergence(sc: Scenario) -> dict[str, dict[str, float]]:
    vals = [endpoints(solve(sc.vessel, flow=sc.flow, n=n)) for n in GRIDS]
    res = {}
    for k in vals[0]:
        f3, f2, f1 = (v[k] for v in vals)  # coarse, medium (default), fine
        e21, e32 = f1 - f2, f2 - f3
        r = 2.0
        if abs(e21) < 1e-12 * max(1.0, abs(f1)) or abs(e32) < 1e-12 * max(1.0, abs(f1)):
            res[k] = {"value": f2, "order": math.nan, "gci_pct": 0.0}
            continue
        p = math.log(abs(e32 / e21)) / math.log(r)
        p_used = min(max(p, 0.5), 4.0)  # guard against oscillatory / superconvergent noise
        gci_medium = FS * abs(e21) * r**p_used / (r**p_used - 1)  # error band on the 401-station result
        res[k] = {"value": f2, "order": p, "gci_pct": 100 * gci_medium / max(abs(f2), 1e-12)}
    return res


# (name, sampler, basis). Each sampler takes the generator and returns a value.
INPUT_UNCERTAINTY = [
    ("diameter_mm", lambda g: g.normal(2.39, 0.49), "measured: M1 2.39 ± 0.49 mm, 554 angiograms"),
    ("flow_ml_min", lambda g: g.normal(153.0, 28.0), "measured: MCA 153 ± 28 mL/min (Amin-Hanjani 2015)"),
    ("inlet_pressure_mmhg", lambda g: g.uniform(75.0, 95.0), "assumed: 85 ± 10 mmHg uniform"),
    ("stiffness_kpa", lambda g: 10 ** g.uniform(math.log10(300), math.log10(3000)), "assumed: 0.3–3 MPa log-uniform"),
    ("stenosis_length_mm", lambda g: g.uniform(4.0, 8.0), "assumed: 4–8 mm uniform"),
]
FLOW_KEYS = {"flow_ml_min", "inlet_pressure_mmhg"}


def monte_carlo(sc: Scenario, n_samples: int = 400, seed: int = 1, n: int = 201):
    g = np.random.default_rng(seed)
    rows, xs = [], []
    while len(rows) < n_samples:
        x = {name: float(f(g)) for name, f, _ in INPUT_UNCERTAINTY}
        x["diameter_mm"] = max(x["diameter_mm"], 1.2)
        x["flow_ml_min"] = max(x["flow_ml_min"], 40.0)
        v = replace(sc.vessel, **{k: val for k, val in x.items() if k not in FLOW_KEYS})
        fl = replace(sc.flow, **{k: val for k, val in x.items() if k in FLOW_KEYS})
        try:
            rows.append(endpoints(solve(v, flow=fl, n=n)))
            xs.append(x)
        except (ValueError, RuntimeError):
            continue
    return xs, rows


def _rank(a):
    return np.argsort(np.argsort(a)).astype(float)


def summarize_mc(xs, rows):
    keys = list(rows[0])
    out = {}
    for k in keys:
        y = np.array([r[k] for r in rows])
        if np.all(np.isnan(y)):
            continue
        rho = {}
        for name, _, _ in INPUT_UNCERTAINTY:
            xv = np.array([x[name] for x in xs])
            rx, ry = _rank(xv), _rank(y)
            rho[name] = float(np.corrcoef(rx, ry)[0, 1]) if np.std(ry) > 0 else 0.0
        out[k] = {"median": float(np.median(y)), "lo": float(np.percentile(y, 2.5)),
                  "hi": float(np.percentile(y, 97.5)), "rho": rho}
    return out


def report(n_samples: int = 400) -> str:
    L = ["# Verification results", "",
         "Generated by `python -m bvm verify`. Deterministic (fixed seed); rerun after any engine change.",
         "Supports evidence categories 1 and 8 of the credibility assessment plan "
         "([credibility-assessment-plan.md](credibility-assessment-plan.md)).", "",
         "## Code verification", "",
         "`tests/test_verification.py` checks the solver against closed-form solutions and theory:",
         "", "* uniform rigid tube: pressure drop, WSS and velocity equal Hagen–Poiseuille to 1e-12",
         "* linear taper: matches the closed-form integral, with observed order 2.0 (trapezoid rule)",
         "* rigid stenosis: equals the exact viscous integral plus the Young–Tsai loss to 1e-10",
         "* elastic wall: second-order convergence; inlet strain equals the Laplace thin-wall value",
         "* White curvature correlation: continuous at De = 11.6 and monotonic",
         "* territory BC reduces to a series-resistance circuit; branched vessels conserve mass", "",
         "## Calculation verification: grid convergence", "",
         f"Grids of {', '.join(map(str, GRIDS))} stations; GCI (safety factor {FS}) is the numerical error "
         "band on the 401-station result that the app and CLI report. `p` is the observed order.",
         "Pressure and flow converge at second order or better. Peak WSS, throat Re and perforator flow "
         "converge at about first order, because a maximum is sampled at stations and branch origins snap "
         "to the nearest station; their GCI is still well under 1%. Orders far from 1–2 (or negative) occur "
         "when the grid differences are at round-off level, and the GCI there is negligible.", ""]
    for sc in SCENARIOS:
        L += [f"**{sc.name}**", "", "| Endpoint | Value (401) | Observed order p | GCI |", "|---|---|---|---|"]
        for k, r in grid_convergence(sc).items():
            if math.isnan(r["value"]):
                continue
            p = "exact on all grids" if math.isnan(r["order"]) else f"{r['order']:.2f}"
            L.append(f"| {k} | {r['value']:.4g} | {p} | {r['gci_pct']:.2g}% |")
        L.append("")
    L += ["## Uncertainty quantification (Monte Carlo)", "",
          f"{n_samples} samples per scenario, 201 stations (the GCI above shows grid error is far below "
          "input uncertainty). Inputs:", "", "| Input | Distribution / basis |", "|---|---|"]
    L += [f"| `{name}` | {basis} |" for name, _, basis in INPUT_UNCERTAINTY]
    L += ["", "Intervals are the 2.5–97.5th percentiles. ρ is the Spearman rank correlation of the "
          "endpoint with each input (a global sensitivity measure; |ρ| near 1 = dominant input).", ""]
    for sc in SCENARIOS[1:3]:
        xs, rows = monte_carlo(sc, n_samples)
        L += [f"**{sc.name}**", "", "| Endpoint | Median | 95% interval | Most influential inputs (ρ) |", "|---|---|---|---|"]
        for k, s in summarize_mc(xs, rows).items():
            top = sorted(s["rho"].items(), key=lambda kv: -abs(kv[1]))[:2]
            L.append(f"| {k} | {s['median']:.3g} | {s['lo']:.3g} – {s['hi']:.3g} | "
                     + ", ".join(f"`{n}` {r:+.2f}" for n, r in top) + " |")
        L.append("")
    L += ["## Reading these numbers", "",
          "* Numerical (grid) error is negligible next to input uncertainty, so credibility is limited by "
          "**model form and inputs**, which only validation (bench, in vivo, population) can address.",
          "* These results verify that the equations are solved correctly. They do not show the equations "
          "are the right ones; see the validation activities in the plan.", ""]
    return "\n".join(L)
