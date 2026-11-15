"""Pulsatile small-signal model (bvm/pulsatile.py)."""
import numpy as np
import pytest

from bvm import Blood, FlowConditions, Vessel, solve
from bvm.pulsatile import N_SAMPLES, synth, waveform_coefficients


def pul(v=None, **k):
    return solve(v or Vessel(n_branches=6), flow=FlowConditions(bc="territory", **k)).summary()


def test_waveform_has_unit_mean_and_the_requested_pi():
    for pi in (0.0, 0.6, 0.9, 1.4):
        w = synth(1.0, waveform_coefficients(pi), N_SAMPLES)
        assert w.mean() == pytest.approx(1.0, abs=1e-12)
        assert w.max() - w.min() == pytest.approx(pi, abs=1e-12)


def test_mean_results_do_not_depend_on_pulsatility_inputs():
    a, b = pul(), pul(pulsatility_index=1.3, territory_tau_s=0.4, heart_rate_bpm=100)
    for key in ("flow", "pressure", "wss", "perfusion"):
        assert a[key] == b[key]


def test_healthy_reference_reproduces_the_chosen_pi():
    for pi in (0.6, 0.9, 1.2):
        p = pul(pulsatility_index=pi)["pulsatile"]
        assert p["PI_inlet"] == pytest.approx(pi, rel=1e-9)
        assert p["incremental_resistance_ratio"] == pytest.approx(1.0, rel=1e-9)


def test_bed_compliance_damps_territory_pulsatility():
    p = pul()["pulsatile"]
    assert 0 < p["PI_territory"] < p["PI_inlet"]
    longer = pul(territory_tau_s=0.4)["pulsatile"]
    assert longer["PI_territory"] < p["PI_territory"]


def test_zero_pi_gives_steady_flow():
    p = pul(pulsatility_index=0.0)["pulsatile"]
    s = pul(pulsatility_index=0.0)
    assert p["PI_inlet"] == pytest.approx(0, abs=1e-12) and p["inlet_pulse_pressure_mmHg"] == pytest.approx(0, abs=1e-9)
    assert p["wss_max_systolic_Pa"] == pytest.approx(s["wss"]["wss_max_Pa"], rel=1e-9)


def test_stenosis_damps_distal_pulsatility_tardus_parvus():
    pis = [pul(Vessel(occlusion_pct=o, n_branches=6))["pulsatile"]["PI_territory"] for o in (0, 50, 70, 90)]
    assert pis[0] >= pis[1] > pis[2] > pis[3]


def test_no_reflow_raises_vessel_pulsatility():
    assert pul(no_reflow_pct=50)["pulsatile"]["PI_inlet"] > pul()["pulsatile"]["PI_inlet"] * 1.3


def test_systolic_fields_exceed_the_mean():
    s = pul(Vessel(occlusion_pct=60, n_branches=6))
    p = s["pulsatile"]
    assert p["Q_inlet_systolic_mL_min"] > s["flow"]["Q_mL_min"] > p["Q_inlet_diastolic_mL_min"]
    assert p["wss_max_systolic_Pa"] > s["wss"]["wss_max_Pa"] > p["wss_max_diastolic_Pa"]
    assert p["Re_max_systolic"] > s["reynolds"]["Re_max"]


def test_waveforms_average_to_the_mean_flows():
    s = pul(Vessel(occlusion_pct=80, n_branches=6), collateral_pct=20)
    p = s["pulsatile"]
    assert np.mean(p["waveform_Q_inlet_mL_min"]) == pytest.approx(s["flow"]["Q_mL_min"], rel=1e-12)
    assert np.mean(p["waveform_Q_territory_mL_min"]) == pytest.approx(s["flow"]["Q_territory_mL_min"], rel=1e-12)


def test_fixed_flow_prescribes_the_waveform():
    p = solve(Vessel(), Blood(), FlowConditions()).summary()["pulsatile"]
    assert p["PI_inlet"] == pytest.approx(0.9, rel=1e-12) and p["waveform_Q_territory_mL_min"] is None
    assert p["Q_inlet_systolic_mL_min"] == pytest.approx(120 * (1 + max(synth(0, waveform_coefficients(0.9), N_SAMPLES))), rel=1e-12)
