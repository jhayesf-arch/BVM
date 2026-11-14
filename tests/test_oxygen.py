"""Hematocrit-coupled viscosity and oxygen delivery (bvm/oxygen.py)."""
import math

import pytest

from bvm import Blood, FlowConditions, Vessel, solve
from bvm.oxygen import Oxygen, pries_relative_viscosity, viscosity_ratio

TERR = FlowConditions(bc="territory")


def test_pries_law_values_and_fahraeus_lindqvist():
    assert pries_relative_viscosity(3000, 45) == pytest.approx(3.2, abs=0.01)   # large-vessel limit
    # plasma: 1 apart from the in vivo endothelial-layer factor (D/(D-1.1))^2 ≈ 1.0007 at 3 mm
    assert pries_relative_viscosity(3000, 0) == pytest.approx((3000 / 2998.9) ** 2, rel=1e-9)
    assert pries_relative_viscosity(25, 45) < pries_relative_viscosity(3000, 45)  # Fåhræus–Lindqvist
    hs = [20, 30, 40, 45, 50, 60]
    vals = [pries_relative_viscosity(3000, h) for h in hs]
    assert all(b > a for a, b in zip(vals, vals[1:]))
    assert viscosity_ratio(3000, 45) == 1.0 and viscosity_ratio(25, 45) == 1.0


def test_normal_hematocrit_leaves_every_result_unchanged():
    v = Vessel(occlusion_pct=70, n_branches=6)
    a = solve(v, Blood(hematocrit=45), TERR)
    b = solve(v, Blood(hematocrit=45, viscosity_from_hematocrit=False), TERR)
    assert a.q == b.q and a.mu_effective == 0.0035 and a.distal_scale == 1.0


def test_oxygen_content_and_baseline_extraction():
    s = solve(Vessel(), flow=TERR).summary()
    cao2 = 1.34 * (45 * 33 / 100) * 0.98 + 0.003 * 95
    assert s["blood"]["CaO2_mL_dL"] == pytest.approx(cao2)
    d = s["oxygen"]["distal"]
    assert d["CBF_mL_100g_min"] == pytest.approx(50.0, rel=1e-6)
    assert d["required_OEF"] == pytest.approx(3.3 / (50 * cao2 / 100), rel=1e-6)
    assert d["band"] == "normal"


def test_delivery_downstream_tracks_flow_because_content_is_unchanged():
    s = solve(Vessel(occlusion_pct=80, n_branches=6), flow=TERR).summary()
    o, f = s["oxygen"], s["flow"]
    assert o["distal"]["DO2_mL_min"] / o["inlet"]["DO2_mL_min"] == pytest.approx(f["Q_distal_mL_min"] / f["Q_mL_min"], rel=1e-12)


def test_optimal_hematocrit_in_a_healthy_vessel():
    do2 = {h: solve(Vessel(n_branches=6), Blood(hematocrit=h), TERR).summary()["oxygen"]["distal"]["DO2_mL_min"]
           for h in (25, 35, 45, 55)}
    assert do2[45] > do2[25] and do2[45] > do2[55] and do2[35] > do2[25]


def test_anemia_raises_required_extraction_behind_a_stenosis():
    oef = [solve(Vessel(occlusion_pct=70, n_branches=6), Blood(hematocrit=h), TERR).summary()["oxygen"]["distal"]["required_OEF"]
           for h in (25, 35, 45)]
    assert oef[0] > oef[1] > oef[2]


def test_oef_bands_and_unmet_demand_for_near_occlusion():
    s = solve(Vessel(occlusion_pct=95, n_branches=6), flow=TERR).summary()["oxygen"]["distal"]
    assert s["band"] == "unmet" and s["required_OEF"] > Oxygen().oef_max


def test_wall_o2_index_is_one_at_the_healthy_inlet():
    s = solve(Vessel(stiffness_kpa=math.inf), flow=TERR).summary()["oxygen"]
    assert s["wall_o2_index_min"] == pytest.approx(1.0, rel=1e-9)
