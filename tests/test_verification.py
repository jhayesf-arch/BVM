"""Code verification (FDA credibility evidence category 1): the solver against exact
solutions and expected convergence order. See docs/credibility-assessment-plan.md.

Each test checks the implementation against mathematics that is independent of it:
closed-form integrals, conservation laws, or the theoretical order of the discretization.
"""
import math

import numpy as np
import pytest

from bvm import FlowConditions, Vessel, solve
from bvm.geometry import radius_profile
from bvm.hemodynamics import YOUNG_TSAI_KT, white_curvature_factor
from bvm.params import ML_MIN_TO_M3_S, MMHG_TO_PA

RIGID = math.inf
MU, RHO = 0.0035, 1060.0
Q = 120 * ML_MIN_TO_M3_S


def dp_mmhg(sol):
    return (sol.pressure[0] - sol.pressure[-1]) / MMHG_TO_PA


@pytest.mark.parametrize("n", [51, 201, 801])
def test_poiseuille_uniform_tube_is_exact(n):
    sol = solve(Vessel(stiffness_kpa=RIGID), n=n)
    r, L = 1.5e-3, 0.02
    assert dp_mmhg(sol) * MMHG_TO_PA == pytest.approx(8 * MU * Q * L / (math.pi * r**4), rel=1e-12)
    assert sol.wss == pytest.approx(np.full(n, 4 * MU * Q / (math.pi * r**3)), rel=1e-12)
    assert sol.velocity == pytest.approx(np.full(n, Q / (math.pi * r**2)), rel=1e-12)


def test_tapered_tube_matches_closed_form_at_second_order():
    """∫ dx / r(x)^4 for a linear taper has a closed form; trapezoid error must fall as h²."""
    r0, L, t = 1.5e-3, 0.02, 0.6
    exact = 8 * MU * Q / math.pi * L / (3 * r0**4 * (t - 1)) * (1 - t**-3)
    errs = []
    for n in (101, 201, 401, 801):
        sol = solve(Vessel(stiffness_kpa=RIGID, taper=t), n=n)
        errs.append(abs(dp_mmhg(sol) * MMHG_TO_PA - exact) / exact)
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(3)]
    assert all(1.9 < p < 2.1 for p in orders), orders
    assert errs[2] < 1e-5  # default grid (401 stations)


def test_stenosis_pressure_drop_equals_viscous_integral_plus_young_tsai():
    v = Vessel(occlusion_pct=70, stiffness_kpa=RIGID)
    s = np.linspace(0, 0.02, 200_001)
    visc = np.trapezoid(8 * MU * Q / (np.pi * radius_profile(v, s) ** 4), s)
    a0, a1 = np.pi * 1.5e-3**2, np.pi * (1.5e-3 * 0.3) ** 2
    loss = YOUNG_TSAI_KT * RHO / 2 * (Q / a0) ** 2 * (a0 / a1 - 1) ** 2
    assert dp_mmhg(solve(v)) * MMHG_TO_PA == pytest.approx(visc + loss, rel=1e-10)


def test_elastic_wall_converges_at_second_order():
    vals = [dp_mmhg(solve(Vessel(occlusion_pct=70), n=n)) for n in (101, 201, 401, 801, 1601)]
    errs = [abs(v - vals[-1]) for v in vals[:-1]]
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(2)]
    assert all(1.8 < p < 2.6 for p in orders), orders


def test_thin_wall_strain_at_inlet_is_laplace():
    v = Vessel(stiffness_kpa=500.0, wall_thickness_mm=0.3)
    sol = solve(v)
    p_in = 85 * MMHG_TO_PA
    assert sol.strain[0] == pytest.approx(p_in * 1.5e-3 / (500e3 * 0.3e-3), rel=1e-9)


def test_white_correlation_limits_and_continuity():
    assert white_curvature_factor(np.array([0.0, 5.0, 11.6]))[...] == pytest.approx([1, 1, 1])
    just_above = white_curvature_factor(np.array([11.6 * (1 + 1e-6)]))[0]
    assert just_above == pytest.approx(1.0, abs=1e-6), "no jump at the De = 11.6 switch"
    de = np.linspace(12, 2000, 200)
    fc = white_curvature_factor(de)
    assert np.all(np.diff(fc) > 0), "friction ratio must rise monotonically with Dean number"


def test_territory_bc_is_a_series_resistance_circuit():
    """Rigid straight vessel: Q = ΔP / (R_vessel + R_distal) exactly, for any occlusion-free inlet drop."""
    fl = FlowConditions(bc="territory", inlet_pressure_mmhg=85, outflow_pressure_mmhg=10, inlet_drop_mmhg=20)
    sol = solve(Vessel(stiffness_kpa=RIGID), flow=fl)
    r_v = 8 * MU * 0.02 / (math.pi * 1.5e-3**4)
    r_d = (75 * MMHG_TO_PA - r_v * Q) / Q
    assert sol.q == pytest.approx(55 * MMHG_TO_PA / (r_v + r_d), rel=1e-8)


def test_branched_vessel_conserves_mass_and_pressure_continuity():
    fl = FlowConditions(bc="territory", branch_flow_pct=3.0)
    sol = solve(Vessel(n_branches=8, branch_start=0.1, branch_end=0.9, occlusion_pct=60), flow=fl)
    assert sol.q == pytest.approx(sol.q_out + sol.branch_q.sum(), rel=1e-12)
    assert np.all(np.diff(sol.pressure) <= 1e-9), "pressure must fall monotonically downstream"
