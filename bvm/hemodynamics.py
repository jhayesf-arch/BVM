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
    reference_q: float | None = None   # territory BC: flow in the healthy reference
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
        if self.reference_q is not None:
            frac = self.q / self.reference_q
            out["perfusion"] = {
                "reference_Q_mL_min": self.reference_q / ML_MIN_TO_M3_S,
                "territory_flow_fraction": frac,
                "ischemia_band": _ischemia_band(frac),
            }
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
                   p_inlet: float) -> VesselSolution:
    v = geom.vessel
    s, kappa, r0 = geom.s, geom.centerline.curvature, geom.radius
    mu, rho = blood.mu, blood.rho
    p_ext = fl.external_pressure_mmhg * MMHG_TO_PA
    eh = math.inf if v.is_rigid else v.stiffness_kpa * 1e3 * v.wall_thickness_mm * 1e-3
    i_t = geom.throat_index
    base_ratio = (geom.base_radius[i_t] / geom.base_radius[0]) ** 2
    sc = v.stenosis_center * v.length_mm * 1e-3
    s_end = sc + v.stenosis_length_mm / 2 * 1e-3

    r = r0.copy()
    for _ in range(200):
        area = np.pi * r**2
        vel = q / area
        re = rho * vel * 2 * r / mu
        dean = re * np.sqrt(r * kappa)
        fc = white_curvature_factor(dean)
        dp = cumulative_trapezoid(8 * mu * q / (np.pi * r**4) * fc, s, initial=0.0)
        if v.occlusion_pct > 0:
            # Young–Tsai expansion loss, spread over the diverging half of the stenosis.
            # upstream area = healthy local area at the throat (equals the inlet for a uniform tube)
            a0, a1 = area[0] * base_ratio, area[i_t]
            loss = YOUNG_TSAI_KT * rho / 2 * (q / a0) ** 2 * (a0 / a1 - 1.0) ** 2
            ramp = np.clip((s - s[i_t]) / max(s_end - s[i_t], 1e-12), 0.0, 1.0)
            dp = dp + loss * ramp
        p = p_inlet - dp
        strain = np.zeros_like(r) if v.is_rigid else (p - p_ext) * r0 / eh
        if v.rigid_span_mm is not None:
            a, b = v.rigid_span_mm
            strain = np.where((s * 1e3 >= a - 1e-9) & (s * 1e3 <= b + 1e-9), 0.0, strain)
        r_new = r0 * (1.0 + np.maximum(strain, -0.9))
        if np.max(np.abs(r_new - r)) < 1e-12:
            break
        r = 0.5 * r + 0.5 * r_new
    else:
        raise RuntimeError("wall-pressure coupling did not converge (wall too soft?)")

    wss = 4 * mu * q / (np.pi * r**3) * fc
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
    return VesselSolution(v, blood, fl, geom, q, r, p, vel, re, dean, wss, strain, warnings=warnings)


def _reference_vessel(v: Vessel) -> Vessel:
    return replace(v, occlusion_pct=0.0, tortuosity=1.0, rigid_span_mm=None)


def solve(vessel: Vessel, blood: Blood | None = None, flow: FlowConditions | None = None,
          reference: Vessel | None = None, n: int = 401) -> VesselSolution:
    """Solve the vessel. See FlowConditions for the two boundary-condition modes.

    For bc="territory" the distal resistance is calibrated on ``reference``
    (default: the same vessel, straight and unoccluded, keeping its taper and bulb)
    carrying flow_ml_min at the undropped inlet pressure.
    """
    blood = blood or Blood()
    flow = flow or FlowConditions()
    geom = build_geometry(vessel, n)
    p_cal = flow.inlet_pressure_mmhg * MMHG_TO_PA                 # calibration (healthy) inlet pressure
    p_in = (flow.inlet_pressure_mmhg - flow.inlet_drop_mmhg) * MMHG_TO_PA  # actual inlet pressure
    q0 = flow.flow_ml_min * ML_MIN_TO_M3_S

    if flow.bc == "fixed_flow":
        return _solve_at_flow(geom, blood, flow, q0, p_in)

    p_out = flow.outflow_pressure_mmhg * MMHG_TO_PA
    ref_geom = build_geometry(reference or _reference_vessel(vessel), n)
    ref = _solve_at_flow(ref_geom, blood, flow, q0, p_cal)
    r_distal = (ref.pressure[-1] - p_out) / q0
    if r_distal <= 0:
        raise ValueError("reference vessel alone drops more pressure than is available; check pressures")

    def residual(q: float) -> float:
        sol = _solve_at_flow(geom, blood, flow, q, p_in)
        return sol.pressure[-1] - p_out - r_distal * q

    q_max = (p_in - p_out) / r_distal
    q = brentq(residual, q_max * 1e-6, q_max, xtol=1e-14, rtol=1e-10)
    sol = _solve_at_flow(geom, blood, flow, q, p_in)
    sol.reference_q = q0
    return sol
