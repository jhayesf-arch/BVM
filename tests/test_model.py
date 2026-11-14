import math

import numpy as np
import pytest

from bvm import Blood, FlowConditions, Vessel, build_geometry, build_mesh, flow_regime, solve
from bvm.params import ML_MIN_TO_M3_S, MMHG_TO_PA

RIGID = math.inf


def test_blood_kinematic_viscosity():
    assert Blood().nu == pytest.approx(3.30e-6, rel=1e-2)


def test_baseline_velocity_and_reynolds_rigid():
    sol = solve(Vessel(stiffness_kpa=RIGID))
    q = 120 * ML_MIN_TO_M3_S
    area = math.pi * 0.0015**2
    v = q / area
    assert sol.velocity[0] == pytest.approx(v)
    assert sol.reynolds[0] == pytest.approx(1060 * v * 0.003 / 0.0035)  # ≈ 257
    assert sol.summary()["reynolds"]["regime_inlet"] == "laminar"


def test_poiseuille_pressure_drop_straight_rigid():
    sol = solve(Vessel(stiffness_kpa=RIGID))
    q, r, L = 120 * ML_MIN_TO_M3_S, 0.0015, 0.020
    expected = 8 * 0.0035 * L * q / (math.pi * r**4)
    assert sol.pressure[0] - sol.pressure[-1] == pytest.approx(expected, rel=1e-6)
    assert sol.wss[0] == pytest.approx(4 * 0.0035 * q / (math.pi * r**3))


@pytest.mark.parametrize("ti", [1.05, 1.2, 1.5])
def test_tortuosity_index_achieved(ti):
    cl = build_geometry(Vessel(tortuosity=ti)).centerline
    arc = np.sum(np.linalg.norm(np.diff(cl.points, axis=0), axis=1))
    chord = np.linalg.norm(cl.points[-1] - cl.points[0])
    assert arc == pytest.approx(0.020, rel=1e-3)
    assert arc / chord == pytest.approx(ti, rel=2e-3)


@pytest.mark.parametrize("pct,basis,ratio", [(50, "diameter", 0.5), (75, "area", 0.5), (70, "diameter", 0.3)])
def test_throat_radius(pct, basis, ratio):
    g = build_geometry(Vessel(occlusion_pct=pct, occlusion_basis=basis))
    assert g.radius.min() == pytest.approx(0.0015 * ratio, rel=1e-3)
    assert g.radius[0] == pytest.approx(0.0015)


@pytest.mark.parametrize("vessel", [Vessel(), Vessel(tortuosity=1.3, occlusion_pct=60)])
def test_mesh_is_closed_and_outward(vessel):
    g = build_geometry(vessel)
    mesh = build_mesh(g, n_theta=64)
    edges = np.sort(np.vstack([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]]), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(counts == 2)  # watertight, manifold
    expected = np.trapezoid(np.pi * g.radius**2, g.s)
    assert mesh.signed_volume() == pytest.approx(expected, rel=5e-3)  # positive ⇒ outward normals


def test_softer_wall_dilates_and_lowers_reynolds():
    stiff = solve(Vessel(stiffness_kpa=5000))
    soft = solve(Vessel(stiffness_kpa=500))
    assert soft.radius[0] > stiff.radius[0]
    assert soft.reynolds[0] < stiff.reynolds[0]  # same Q, larger D ⇒ lower Re (Re ∝ 1/D)
    p = 85 * MMHG_TO_PA
    assert stiff.strain[0] == pytest.approx(p * 0.0015 / (5000e3 * 0.0003), rel=1e-3)


def test_soft_substrate_warns():
    assert any("small-strain" in w for w in solve(Vessel(stiffness_kpa=50)).warnings)


def test_stenosis_raises_throat_reynolds_and_wss():
    base = solve(Vessel(stiffness_kpa=RIGID))
    sten = solve(Vessel(occlusion_pct=50, stiffness_kpa=RIGID))
    i = np.argmin(sten.radius)
    assert sten.reynolds[i] == pytest.approx(2 * base.reynolds[0], rel=1e-3)  # Re ∝ 1/D at fixed Q
    assert sten.wss.max() == pytest.approx(8 * base.wss[0], rel=1e-3)        # WSS ∝ 1/r³
    # Baseline Re_inlet ≈ 257 sits under the 300 screening cutoff; higher flow crosses it.
    assert not sten.summary()["reynolds"]["post_stenotic_disturbance_possible"]
    fast = solve(Vessel(occlusion_pct=50, stiffness_kpa=RIGID), flow=FlowConditions(flow_ml_min=200))
    assert fast.summary()["reynolds"]["post_stenotic_disturbance_possible"]


def test_territory_bc_calibrates_and_occlusion_reduces_flow():
    fl = FlowConditions(bc="territory")
    healthy = solve(Vessel(), flow=fl)
    assert healthy.q / ML_MIN_TO_M3_S == pytest.approx(120, rel=1e-6)
    fracs = [solve(Vessel(occlusion_pct=p), flow=fl).summary()["perfusion"]["territory_flow_fraction"]
             for p in (30, 60, 80, 90)]
    assert all(a > b for a, b in zip(fracs, fracs[1:]))
    assert fracs[0] > 0.95 and fracs[-1] < 0.5


def test_flow_regime_labels():
    assert flow_regime(500) == "laminar"
    assert flow_regime(3000) == "transitional"
    assert flow_regime(5000) == "turbulent"


@pytest.mark.parametrize("kw", [dict(occlusion_pct=100), dict(tortuosity=0.9), dict(diameter_mm=0),
                                dict(occlusion_pct=50, stenosis_center=0.05)])
def test_invalid_inputs(kw):
    with pytest.raises(ValueError):
        Vessel(**kw)
