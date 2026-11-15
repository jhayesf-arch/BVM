"""Input parameters for the MCA vessel model.

All user-facing inputs are in clinical/lab units (mm, kPa, mmHg, mL/min);
everything is converted to SI internally.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

MMHG_TO_PA = 133.322
ML_MIN_TO_M3_S = 1e-6 / 60.0


@dataclass(frozen=True)
class Blood:
    """Newtonian blood. Defaults are the conventional CFD values."""

    mu: float = 0.0035   # dynamic viscosity at a hematocrit of 45%, Pa·s
    rho: float = 1060.0  # density, kg/m³
    hematocrit: float = 45.0               # %, sets O₂ content and (if coupled) viscosity
    viscosity_from_hematocrit: bool = True  # scale μ and the microvascular bed by Pries et al. 1994
    micro_diameter_um: float = 25.0         # effective diameter of the resistance vessels downstream (assumed)
    # Shear-thinning blood in the modeled vessel: Carreau law (Cho & Kensey 1991) scaled so its
    # high-shear limit equals mu (after any hematocrit coupling). See hemodynamics.carreau_factor.
    non_newtonian: bool = False

    @property
    def nu(self) -> float:
        """Kinematic viscosity, m²/s (derived: μ/ρ ≈ 3.30e-6)."""
        return self.mu / self.rho


@dataclass(frozen=True)
class Vessel:
    """Geometry and wall of a single vessel segment (default: MCA M1).

    The four modulated inputs are ``diameter_mm``, ``tortuosity``,
    ``occlusion_pct`` and ``stiffness_kpa``.
    """

    diameter_mm: float = 3.0        # unloaded luminal diameter (M1 is ~2.5–3.5 mm)
    length_mm: float = 20.0         # arc length along the centerline, kept fixed as tortuosity changes
    tortuosity: float = 1.0         # tortuosity index = arc length / chord length (1.0 = straight)
    tortuosity_waves: float = 1.0   # number of full sinusoidal bends when tortuosity > 1
    occlusion_pct: float = 0.0      # stenosis severity, 0 to <100 %
    occlusion_basis: str = "diameter"  # "diameter" (NASCET-style) or "area"
    stenosis_length_mm: float = 6.0    # axial length of the cosine-shaped narrowing
    stenosis_center: float = 0.5       # stenosis location as a fraction of arc length
    stiffness_kpa: float = 1000.0   # Young's modulus of the wall/substrate; math.inf = rigid
    wall_thickness_mm: float = 0.3
    # Arc-length span (mm) held non-distensible, e.g. a segment scaffolded by an implanted stent.
    rigid_span_mm: tuple[float, float] | None = None
    # Anatomy (healthy shape, kept in the territory reference vessel):
    taper: float = 1.0              # outlet / inlet unloaded diameter, linear along arc length
    bulb_pct: float = 0.0           # peak local dilation, % of the local diameter (e.g. carotid bulb)
    bulb_center: float = 0.5        # bulb location, fraction of arc length
    bulb_length_mm: float = 8.0     # axial length of the cosine-shaped dilation
    # Side branches (e.g. lenticulostriate perforators off M1), evenly spaced between
    # branch_start and branch_end (fractions of arc length). Each drains to the outflow
    # pressure through its own lumped resistance, calibrated on the healthy reference.
    n_branches: int = 0
    branch_start: float = 0.3
    branch_end: float = 0.8

    def __post_init__(self) -> None:
        if self.diameter_mm <= 0 or self.length_mm <= 0:
            raise ValueError("diameter_mm and length_mm must be positive")
        if self.tortuosity < 1.0:
            raise ValueError("tortuosity index must be >= 1.0 (1.0 = straight)")
        if not 0.0 <= self.occlusion_pct < 100.0:
            raise ValueError("occlusion_pct must be in [0, 100); 100% is a full occlusion with zero flow")
        if self.occlusion_basis not in ("diameter", "area"):
            raise ValueError("occlusion_basis must be 'diameter' or 'area'")
        if self.stiffness_kpa <= 0 or self.wall_thickness_mm <= 0:
            raise ValueError("stiffness_kpa and wall_thickness_mm must be positive")
        half = self.stenosis_length_mm / 2
        center = self.stenosis_center * self.length_mm
        if self.occlusion_pct > 0 and (center - half < 0 or center + half > self.length_mm):
            raise ValueError("stenosis extends past the ends of the vessel")
        if not 0.2 <= self.taper <= 2.0:
            raise ValueError("taper (outlet/inlet diameter) must be in [0.2, 2]")
        if not 0.0 <= self.bulb_pct <= 150.0:
            raise ValueError("bulb_pct must be in [0, 150]")
        bh, bc = self.bulb_length_mm / 2, self.bulb_center * self.length_mm
        if self.bulb_pct > 0 and (bc - bh < -1e-9 or bc + bh > self.length_mm + 1e-9):
            raise ValueError("bulb extends past the ends of the vessel")
        if self.n_branches < 0 or not 0.0 < self.branch_start <= self.branch_end < 1.0:
            raise ValueError("need n_branches >= 0 and 0 < branch_start <= branch_end < 1")

    @property
    def throat_radius_ratio(self) -> float:
        """Throat radius / reference radius."""
        p = self.occlusion_pct / 100.0
        return 1.0 - p if self.occlusion_basis == "diameter" else math.sqrt(1.0 - p)

    @property
    def area_stenosis_pct(self) -> float:
        return 100.0 * (1.0 - self.throat_radius_ratio**2)

    @property
    def diameter_stenosis_pct(self) -> float:
        return 100.0 * (1.0 - self.throat_radius_ratio)

    @property
    def is_rigid(self) -> bool:
        return math.isinf(self.stiffness_kpa)


@dataclass(frozen=True)
class FlowConditions:
    """Boundary conditions.

    bc="fixed_flow": the inlet flow is prescribed (Q = flow_ml_min no matter
        the occlusion). Good for isolating local hemodynamics (WSS, Re, ΔP).
    bc="territory": the inlet pressure is prescribed and the vessel drains
        into a lumped distal resistance (the MCA territory) calibrated so that
        a healthy reference vessel carries flow_ml_min. Occlusion then
        *reduces* flow, which is the ischemia-relevant case.
    """

    flow_ml_min: float = 120.0
    inlet_pressure_mmhg: float = 85.0     # mean pressure at the M1 inlet (assumption)
    outflow_pressure_mmhg: float = 10.0   # venous/ICP back-pressure of the territory (assumption)
    external_pressure_mmhg: float = 0.0   # pressure outside the wall (0 for a bench BVM, ~ICP in vivo)
    heart_rate_bpm: float = 70.0          # only used for the Womersley number
    bc: str = "fixed_flow"
    # Pressure lost upstream of the inlet (systemic hypotension, proximal disease). The territory
    # resistance is still calibrated at inlet_pressure_mmhg, so this lowers perfusion.
    inlet_drop_mmhg: float = 0.0
    # Share of the healthy inlet flow taken by all side branches together (Vessel.n_branches).
    branch_flow_pct: float = 1.9
    # Territory BC only. Leptomeningeal collaterals: territory flow they would supply with the vessel
    # fully occluded, as % of the healthy territory flow (0 = none, the worst case). They feed the
    # territory node from a pial supply at the healthy outlet pressure, so they carry no flow in health.
    collateral_pct: float = 0.0
    # Territory BC only. No-reflow (microvascular obstruction): % of the distal microvascular bed
    # obstructed; the distal resistance is divided by (1 - fraction). Perforator beds are not affected.
    no_reflow_pct: float = 0.0

    def __post_init__(self) -> None:
        if self.bc not in ("fixed_flow", "territory"):
            raise ValueError("bc must be 'fixed_flow' or 'territory'")
        if self.flow_ml_min <= 0:
            raise ValueError("flow_ml_min must be positive")
        if self.inlet_pressure_mmhg <= self.outflow_pressure_mmhg:
            raise ValueError("inlet pressure must exceed outflow pressure")
        if not 0.0 <= self.branch_flow_pct < 100.0:
            raise ValueError("branch_flow_pct must be in [0, 100)")
        if self.inlet_drop_mmhg < 0 or self.inlet_pressure_mmhg - self.inlet_drop_mmhg <= self.outflow_pressure_mmhg:
            raise ValueError("inlet_drop_mmhg must be >= 0 and leave the inlet above the outflow pressure")
        if not 0.0 <= self.collateral_pct < 100.0 or not 0.0 <= self.no_reflow_pct < 100.0:
            raise ValueError("collateral_pct and no_reflow_pct must be in [0, 100)")
