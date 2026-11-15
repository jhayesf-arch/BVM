"""Collateral supply, no-reflow and shear-thinning (Carreau) blood."""
import math

import numpy as np
import pytest

from bvm import Blood, FlowConditions, Vessel, solve
from bvm.hemodynamics import CARREAU_MU_RATIO, carreau_factor


def terr(**k):
    return FlowConditions(bc="territory", **k)


def frac(v, b=None, **k):
    return solve(v, b or Blood(), terr(**k)).summary()["perfusion"]["territory_flow_fraction"]


def test_defaults_leave_results_unchanged():
    v = Vessel(occlusion_pct=80, n_branches=6)
    a = solve(v, flow=terr()).summary()
    b = solve(v, Blood(non_newtonian=False), terr(collateral_pct=0, no_reflow_pct=0)).summary()
    assert a == b and a["flow"]["Q_collateral_mL_min"] == 0.0


def test_collaterals_carry_no_flow_in_the_healthy_vessel():
    s = solve(Vessel(n_branches=6), flow=terr(collateral_pct=40)).summary()
    assert s["flow"]["Q_collateral_mL_min"] < 1e-6      # root-finder tolerance only
    assert s["perfusion"]["territory_flow_fraction"] == pytest.approx(1.0, rel=1e-9)


def test_collateral_fraction_is_the_share_supplied_at_full_occlusion():
    # rigid, no branches: as the stenosis approaches occlusion the territory flow tends to collateral_pct
    v = Vessel(occlusion_pct=99.5, stiffness_kpa=math.inf, stenosis_length_mm=10)
    f0, f30 = frac(v), frac(v, collateral_pct=30)
    assert f0 < 0.01 and f30 == pytest.approx(0.30, abs=0.01)


def test_collaterals_raise_flow_more_when_the_stenosis_is_worse():
    gain = [frac(Vessel(occlusion_pct=o, n_branches=6), collateral_pct=30) - frac(Vessel(occlusion_pct=o, n_branches=6))
            for o in (50, 80, 95)]
    assert 0 <= gain[0] < gain[1] < gain[2]


def test_territory_node_balances_flow():
    s = solve(Vessel(occlusion_pct=90, n_branches=6), flow=terr(collateral_pct=25, no_reflow_pct=30))
    sm = s.summary()["flow"]
    assert sm["Q_territory_mL_min"] == pytest.approx(sm["Q_distal_mL_min"] + sm["Q_collateral_mL_min"], rel=1e-12)
    assert sm["Q_collateral_mL_min"] > 0


def test_no_reflow_caps_flow_even_with_a_fully_open_vessel():
    v = Vessel(n_branches=6)
    f = frac(v, no_reflow_pct=40)
    assert 0.6 < f < 0.62                    # the open vessel's own resistance takes a little of the drop
    assert frac(v, no_reflow_pct=40) < frac(v, no_reflow_pct=20) < 1


def test_no_reflow_lowers_oxygen_delivery_with_flow():
    s = solve(Vessel(n_branches=6), flow=terr(no_reflow_pct=50)).summary()
    assert s["oxygen"]["distal"]["CBF_mL_100g_min"] == pytest.approx(50 * s["perfusion"]["territory_flow_fraction"], rel=1e-9)


def test_carreau_factor_limits():
    assert carreau_factor(0.0) == pytest.approx(CARREAU_MU_RATIO)
    assert 1.0 < carreau_factor(1e7) < 1.001
    g = np.array([1, 10, 100, 1000])
    assert np.all(np.diff(carreau_factor(g)) < 0)          # shear-thinning
    assert 1.05 < float(carreau_factor(754)) < 1.15        # healthy M1 wall shear rate


def test_non_newtonian_raises_viscosity_and_wss_at_the_same_flow():
    n, b = solve(Vessel()), solve(Vessel(), Blood(non_newtonian=True))
    ratio = b.wss[0] / n.wss[0]
    assert ratio == pytest.approx(float(carreau_factor(4 * n.q / (np.pi * n.radius[0] ** 3))), rel=1e-3)
    assert np.all(b.mu_apparent >= 0.0035)


def test_non_newtonian_viscosity_is_highest_where_shear_is_lowest():
    s = solve(Vessel(occlusion_pct=70, n_branches=6), Blood(non_newtonian=True), terr())
    rate = s.wss / s.mu_apparent
    assert s.mu_apparent[np.argmin(rate)] == pytest.approx(s.mu_apparent.max())


def test_non_newtonian_healthy_territory_still_calibrates_to_full_flow():
    assert frac(Vessel(n_branches=6), Blood(non_newtonian=True)) == pytest.approx(1.0, rel=1e-9)
