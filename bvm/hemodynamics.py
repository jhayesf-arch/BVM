"""Reduced-order steady hemodynamics evaluated on the 3D vessel geometry.

Each centerline station is treated as locally fully-developed Poiseuille
flow, corrected for curvature (Dean flow) and with a Young–Tsai expansion
loss downstream of the stenosis throat. The wall deforms with the local
transmural pressure through a thin-wall linear-elastic law, which is how
substrate stiffness enters. This is a fast screening model; the exported
geometry is meant to be handed to a full 3D CFD/FSI solver for detail.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from .geometry import VesselGeometry, build_geometry
from .params import MMHG_TO_PA, ML_MIN_TO_M3_S, Blood, FlowConditions, Vessel

# Pipe-flow regime limits (steady flow in a straight, smooth tube).
RE_LAMINAR_MAX = 2000.0
RE_TURBULENT_MIN = 4000.0
# Stenotic jets break down far below the straight-pipe limit. DNS of a 75%-area
# stenosis (Varghese, Frankel & Fischer, J Fluid Mech 2007) found laminar flow at
# upstream Re = 500 for a perfectly axisymmetric throat, but transition downstream
# once a small eccentricity was added. So for >=50% area stenoses, an upstream Re
# of a few hundred is flagged as "post-stenotic disturbance possible" rather than
# trusted as laminar. The 300 cutoff is a conservative screening choice, not a
# published threshold; confirm with 3D CFD.
RE_STENOTIC_DISTURBANCE = 300.0
YOUNG_TSAI_KT = 1.52        # expansion-loss coefficient (Young & Tsai, 1973)
LOW_WSS_PA = 0.4            # below this: atheroprone / low-shear (Malek et al., 1999)
HIGH_WSS_PA = 10.0          # rough onset of platelet-activating shear in stenoses (screening flag)
SMALL_STRAIN_LIMIT = 0.10   # linear thin-wall law is not trustworthy beyond ~10% strain


def flow_regime(re: float) -> str:
    if re < RE_LAMINAR_MAX:
        return "laminar"
    if re < RE_TURBULENT_MIN:
        return "transitional"
    return "turbulent"


def white_curvature_factor(dean: np.ndarray) -> np.ndarray:
    """Curved/straight friction ratio, White (1929) laminar correlation (11.6 < De < 2000)."""
    dean = np.asarray(dean, dtype=float)
    out = np.ones_like(dean)
    m = dean > 11.6
    out[m] = 1.0 / (1.0 - (1.0 - (11.6 / dean[m]) ** 0.45) ** (1.0 / 0.45))
    return out


@dataclass
class VesselSolution:
    vessel: Vessel
    blood: Blood
    flow: FlowConditions
    geometry: VesselGeometry
    q: float                      # m³/s actually carried
    radius: np.ndarray            # loaded (pressurized) radius, m
    pressure: np.ndarray          # Pa, gauge
    velocity: np.ndarray          # mean velocity, m/s
    reynolds: np.ndarray
    dean: np.ndarray
    wss: np.ndarray               # Pa
    strain: np.ndarray            # circumferential strain of the wall
    reference_q: float | None = None   # territory BC: inlet flow of the healthy reference
    q_out: float = 0.0                  # m³/s leaving the outlet (distal territory)
    branch_q: np.ndarray = field(default_factory=lambda: np.zeros(0))   # m³/s into each side branch
    reference_q_out: float | None = None
    reference_branch_q: np.ndarray | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def s(self) -> np.ndarray:
        return self.geometry.s

    def summary(self) -> dict:
        b, v, fl = self.blood, self.vessel, self.flow
        i_t = self.geometry.throat_index if v.occlusion_pct > 0 else int(np.argmin(self.radius))
        d_in = 2 * self.radius[0]
        re_in, re_throat, re_max = float(self.reynolds[0]), float(self.reynolds[i_t]), float(self.reynolds.max())
        omega = 2 * math.pi * fl.heart_rate_bpm / 60.0
        p_in, p_out = self.pressure[0], self.pressure[-1]
        out = {
            "inputs": {
                "diameter_mm": v.diameter_mm,
                "tortuosity_index": v.tortuosity,
                "occlusion_pct": v.occlusion_pct,
                "occlusion_basis": v.occlusion_basis,
                "stiffness_kpa": v.stiffness_kpa if not v.is_rigid else "rigid",
                "length_mm": v.length_mm,
                "wall_thickness_mm": v.wall_thickness_mm,
                "taper": v.taper,
                "bulb_pct": v.bulb_pct,
                "boundary_condition": fl.bc,
                "inlet_drop_mmHg": fl.inlet_drop_mmhg,
            },
            "blood": {"mu_Pa_s": b.mu, "rho_kg_m3": b.rho, "nu_m2_s": b.nu},
            "flow": {
                "Q_mL_min": self.q / ML_MIN_TO_M3_S,
                "Q_distal_mL_min": self.q_out / ML_MIN_TO_M3_S,
                "Q_branches_mL_min": float(self.branch_q.sum()) / ML_MIN_TO_M3_S,
                "mean_velocity_inlet_m_s": float(self.velocity[0]),
                "mean_velocity_throat_m_s": float(self.velocity[i_t]),
                "loaded_diameter_inlet_mm": d_in * 1e3,
                "loaded_diameter_throat_mm": 2e3 * self.radius[i_t],
            },
            "reynolds": {
                "Re_inlet": re_in,
                "regime_inlet": flow_regime(re_in),
                "Re_throat": re_throat,
                "regime_throat": flow_regime(re_throat),
                "Re_max": re_max,
                "regime_max": flow_regime(re_max),
                "criteria": f"laminar < {RE_LAMINAR_MAX:.0f} <= transitional < {RE_TURBULENT_MIN:.0f} <= turbulent",
                "post_stenotic_disturbance_possible": bool(
                    v.occlusion_pct > 0 and v.area_stenosis_pct >= 50 and re_in >= RE_STENOTIC_DISTURBANCE
                ),
                "Dean_max": float(self.dean.max()),
                "Womersley_inlet": float(d_in / 2 * math.sqrt(omega / b.nu)),
                "entrance_length_mm": 0.05 * re_in * d_in * 1e3,
            },
            "pressure": {
                "P_inlet_mmHg": p_in / MMHG_TO_PA,
                "P_outlet_mmHg": p_out / MMHG_TO_PA,
                "dP_vessel_mmHg": (p_in - p_out) / MMHG_TO_PA,
                "Pd_over_Pa": float(p_out / p_in),
            },
            "wss": {
                "wss_inlet_Pa": float(self.wss[0]),
                "wss_max_Pa": float(self.wss.max()),
                "wss_min_Pa": float(self.wss.min()),
                "low_wss_fraction": float(np.mean(self.wss < LOW_WSS_PA)),
                "high_wss_flag": bool(self.wss.max() > HIGH_WSS_PA),
            },
            "wall": {
                "max_circumferential_strain": float(self.strain.max()),
                "moens_korteweg_wave_speed_m_s": (
                    None if v.is_rigid
                    else math.sqrt(v.stiffness_kpa * 1e3 * v.wall_thickness_mm * 1e-3 / (b.rho * d_in))
                ),
            },
            "stenosis": {
                "diameter_stenosis_pct": v.diameter_stenosis_pct,
                "area_stenosis_pct": v.area_stenosis_pct,
            },
            "warnings": list(self.warnings),
        }
        if v.occlusion_pct > 0:
            sc, half = v.stenosis_center * v.length_mm * 1e-3, v.stenosis_length_mm / 2 * 1e-3
            pa = float(np.interp(sc - half, self.s, self.pressure))
            pb = float(np.interp(sc + half, self.s, self.pressure))
            out["stenosis"]["dP_lesion_mmHg"] = (pa - pb) / MMHG_TO_PA
            out["stenosis"]["Pd_over_Pa_lesion"] = pb / pa
        if self.reference_q is not None:
            # distal (outlet) flow vs the healthy reference: the primary perfusion endpoint
            frac = self.q_out / self.reference_q_out
            out["perfusion"] = {
                "reference_Q_mL_min": self.reference_q / ML_MIN_TO_M3_S,
                "reference_Q_distal_mL_min": self.reference_q_out / ML_MIN_TO_M3_S,
                "territory_flow_fraction": frac,
                "ischemia_band": _ischemia_band(frac),
            }
            if len(self.branch_q) and self.reference_branch_q.sum() > 0:
                out["perfusion"]["branch_flow_fraction"] = float(self.branch_q.sum() / self.reference_branch_q.sum())
                out["perfusion"]["branch_flow_fractions"] = (self.branch_q / self.reference_branch_q).tolist()
        return out


def _ischemia_band(frac: float) -> str:
    # Approximate CBF thresholds relative to ~50 mL/100 g/min normal:
    # penumbra below ~20 (40%), infarct core below ~10–12 (~20%).
    # Heuristic only: no autoregulation or collaterals are modeled.
    if frac >= 0.8:
        return "normal"
    if frac >= 0.4:
        return "oligemia"
    if frac >= 0.2:
        return "penumbra-range"
    return "core-range"


def _solve_at_flow(geom: VesselGeometry, blood: Blood, fl: FlowConditions, q: float,
                   p_inlet: float, branch_q: np.ndarray | None = None,
                   branch_r: np.ndarray | None = None, p_out: float = 0.0) -> VesselSolution:
    """Solve at inlet flow ``q``. Side branches either take prescribed flows (``branch_q``,
    used to calibrate) or drain through resistances ``branch_r`` to ``p_out``.

    Flow is piecewise constant between branch origins; the pressure is marched from the inlet
    so each branch sees the local pressure, including any upstream stenosis loss.
    """
    v = geom.vessel
    s, kappa, r0 = geom.s, geom.centerline.curvature, geom.radius
    n = len(s)
    mu, rho = blood.mu, blood.rho
    p_ext = fl.external_pressure_mmhg * MMHG_TO_PA
    eh = math.inf if v.is_rigid else v.stiffness_kpa * 1e3 * v.wall_thickness_mm * 1e-3
    i_t = geom.throat_index
    base_ratio = (geom.base_radius[i_t] / geom.base_radius[0]) ** 2
    sc = v.stenosis_center * v.length_mm * 1e-3
    s_end = sc + v.stenosis_length_mm / 2 * 1e-3
    ramp = np.clip((s - s[i_t]) / max(s_end - s[i_t], 1e-12), 0.0, 1.0) if v.occlusion_pct > 0 else np.zeros(n)
    b_idx = geom.branch_index
    starts = [0] + sorted(set(int(i) for i in b_idx))
    ends = starts[1:] + [n - 1]

    r = r0.copy()
    for _ in range(200):
        area0 = np.pi * r[0] ** 2 * base_ratio
        p = np.empty(n); qn = np.empty(n); vel = np.empty(n); re = np.empty(n); dean = np.empty(n); fc = np.empty(n)
        bq = np.zeros(len(b_idx))
        qq, cum, loss = q, 0.0, 0.0
        for a, b in zip(starts, ends):
            if a > 0:  # side branches leaving at node a
                for k in np.flatnonzero(b_idx == a):
                    bq[k] = branch_q[k] if branch_q is not None else max(0.0, (p[a] - p_out) / branch_r[k])
                    qq -= bq[k]
            sl = slice(a, b + 1)
            rr = r[sl]
            v_ = qq / (np.pi * rr**2)
            re_ = rho * v_ * 2 * rr / mu
            de_ = re_ * np.sqrt(rr * kappa[sl])
            fc_ = white_curvature_factor(de_)
            g = 8 * mu * qq / (np.pi * rr**4) * fc_
            if v.occlusion_pct > 0 and a <= i_t < b + (b == n - 1):
                # Young–Tsai expansion loss with the flow through the throat, spread over the
                # diverging half of the stenosis; upstream area = healthy local area at the throat.
                loss = YOUNG_TSAI_KT * rho / 2 * (qq / area0) ** 2 * (area0 / (np.pi * r[i_t] ** 2) - 1.0) ** 2
            cumdp = cum + np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s[sl]))])
            p[sl] = p_inlet - cumdp - loss * ramp[sl]
            cum = cumdp[-1]
            m = slice(a, b + 1) if b == n - 1 else slice(a, b)
            k = b + 1 - a if b == n - 1 else b - a
            qn[m], vel[m], re[m], dean[m], fc[m] = qq, v_[:k], re_[:k], de_[:k], fc_[:k]
        strain = np.zeros_like(r) if v.is_rigid else (p - p_ext) * r0 / eh
        if v.rigid_span_mm is not None:
            a_, b_ = v.rigid_span_mm
            strain = np.where((s * 1e3 >= a_ - 1e-9) & (s * 1e3 <= b_ + 1e-9), 0.0, strain)
        r_new = r0 * (1.0 + np.maximum(strain, -0.9))
        if np.max(np.abs(r_new - r)) < 1e-12:
            break
        r = 0.5 * r + 0.5 * r_new
    else:
        raise RuntimeError("wall-pressure coupling did not converge (wall too soft?)")

    wss = 4 * mu * qn / (np.pi * r**3) * fc
    warnings = []
    if strain.max() > SMALL_STRAIN_LIMIT:
        warnings.append(
            f"max wall strain {strain.max():.0%} exceeds the {SMALL_STRAIN_LIMIT:.0%} small-strain limit; "
            "the linear thin-wall law is unreliable at this stiffness (use a hyperelastic FSI model)."
        )
    if dean.max() > 2000:
        warnings.append("Dean number > 2000, outside the White curvature correlation range.")
    le = 0.05 * re[0] * 2 * r[0]
    if le > s[-1]:
        warnings.append(
            f"entrance length ({le*1e3:.0f} mm) exceeds vessel length: flow is still developing, so the "
            "fully-developed Poiseuille WSS here underestimates true wall shear near the inlet."
        )
    if re.max() >= RE_LAMINAR_MAX:
        warnings.append("Re exceeds the laminar limit somewhere; the laminar loss model is not valid there.")
    sol = VesselSolution(v, blood, fl, geom, q, r, p, vel, re, dean, wss, strain, warnings=warnings)
    sol.q_out, sol.branch_q = float(qq), bq
    return sol


def _reference_vessel(v: Vessel) -> Vessel:
    return replace(v, occlusion_pct=0.0, tortuosity=1.0, rigid_span_mm=None)


def solve(vessel: Vessel, blood: Blood | None = None, flow: FlowConditions | None = None,
          reference: Vessel | None = None, n: int = 401) -> VesselSolution:
    """Solve the vessel. See FlowConditions for the two boundary-condition modes.

    For bc="territory" the distal resistance is calibrated on ``reference``
    (default: the same vessel, straight and unoccluded, keeping its taper and bulb)
    carrying flow_ml_min at the undropped inlet pressure. Side branches, if any, are
    calibrated on the same reference to take branch_flow_pct of that flow, in both modes.
    """
    blood = blood or Blood()
    flow = flow or FlowConditions()
    geom = build_geometry(vessel, n)
    p_cal = flow.inlet_pressure_mmhg * MMHG_TO_PA                 # calibration (healthy) inlet pressure
    p_in = (flow.inlet_pressure_mmhg - flow.inlet_drop_mmhg) * MMHG_TO_PA  # actual inlet pressure
    p_out = flow.outflow_pressure_mmhg * MMHG_TO_PA
    q0 = flow.flow_ml_min * ML_MIN_TO_M3_S

    nb = len(geom.branch_index)
    ref_bq = np.full(nb, q0 * flow.branch_flow_pct / 100 / nb) if nb else np.zeros(0)
    need_ref = flow.bc == "territory" or nb > 0
    if need_ref:
        ref_geom = build_geometry(reference or _reference_vessel(vessel), n)
        ref = _solve_at_flow(ref_geom, blood, flow, q0, p_cal, branch_q=ref_bq)
        p_branch = ref.pressure[ref_geom.branch_index] if nb else np.zeros(0)
        if np.any(p_branch - p_out <= 0) or ref.pressure[-1] - p_out <= 0:
            raise ValueError("reference vessel alone drops more pressure than is available; check pressures")
        with np.errstate(divide="ignore"):  # a branch given zero flow has infinite resistance
            branch_r = (p_branch - p_out) / ref_bq if nb else np.zeros(0)
        r_distal = (ref.pressure[-1] - p_out) / ref.q_out
    else:
        branch_r = np.zeros(0)

    if flow.bc == "fixed_flow":
        return _solve_at_flow(geom, blood, flow, q0, p_in, branch_r=branch_r, p_out=p_out)

    def residual(q: float) -> float:
        sol = _solve_at_flow(geom, blood, flow, q, p_in, branch_r=branch_r, p_out=p_out)
        return sol.pressure[-1] - p_out - r_distal * sol.q_out

    q_max = (p_in - p_out) * (1.0 / r_distal + float(np.sum(1.0 / branch_r)))
    q = brentq(residual, q_max * 1e-6, q_max, xtol=1e-14, rtol=1e-10)
    sol = _solve_at_flow(geom, blood, flow, q, p_in, branch_r=branch_r, p_out=p_out)
    sol.reference_q, sol.reference_q_out, sol.reference_branch_q = q0, ref.q_out, ref_bq
    return sol
