"""Pulsatile flow about the steady (mean) solution: small-signal, frequency-domain lumped model.

The steady solver gives the mean flow. Pulsation is added as harmonics of the heart rate through a
linearized circuit around that mean state:

    inlet pressure pulse ──[ vessel: R_v + jωL ]── territory node ──┬── R_d (microvascular bed) ── venous
                                                                    ├── C   (bed compliance)
                                                                    └── R_c (collaterals, if they carry flow)

* R_v is the vessel's *incremental* resistance dΔP/dQ at the mean flow (finite difference of the full
  nonlinear solver, so it includes the stenosis expansion loss, which grows with Q², and curvature).
  L = ρ∫ds/A is the fluid inertance. Side-branch flows are held at their means (≈2% of the flow).
* The bed compliance comes from its time constant, C = τ / R_d,healthy, τ ≈ 0.2 s (Kasprowicz et al.,
  "time constant of the cerebral arterial bed", ~0.15–0.25 s in health). It is a property of the bed, so
  no-reflow (higher R_d) lengthens the time constant.
* Calibration: in the healthy reference vessel the inlet flow waveform is w(t)·Q_ref, an idealized shape
  (four harmonics peaking at 15% of the cycle; assumed, not a measured waveform) scaled to the chosen
  Gosling pulsatility index PI = (Q_max − Q_min)/Q_mean (TCD normal for the MCA ≈ 0.6–1.1). That fixes
  the inlet pressure pulse, which is then applied unchanged to the diseased or treated vessel.
* Peak-systolic and end-diastolic fields re-run the full nonlinear solver (quasi-steady) at the
  instantaneous inlet flow and pressure. With a Womersley number of 2–3 in M1, quasi-steady WSS is a
  reasonable estimate of the peak; the phase lag and near-wall overshoot need 3D pulsatile CFD.

This shows the classic distal damping behind a tight stenosis (low PI, "tardus–parvus") and how far a
device restores it, plus systolic WSS / Re that the mean flow understates.
"""
from __future__ import annotations

import math

import numpy as np

N_SAMPLES = 256          # samples per cycle for max/min/PI
N_OUT = 64               # samples per cycle returned for plotting
WAVE_AMPS = (1.0, 0.5, 0.25, 0.125)
WAVE_PEAK = 0.15         # fraction of the cycle at which the idealized waveform peaks


def waveform_coefficients(pi: float) -> np.ndarray:
    """Complex harmonic coefficients c_k (k = 1..4) of w(t) = 1 + Re Σ c_k e^{jkωt}, scaled to PI."""
    t = np.arange(N_SAMPLES) / N_SAMPLES
    g = sum(a * np.cos(2 * math.pi * (k + 1) * (t - WAVE_PEAK)) for k, a in enumerate(WAVE_AMPS))
    s = pi / (g.max() - g.min())
    return np.array([s * a * np.exp(-2j * math.pi * (k + 1) * WAVE_PEAK) for k, a in enumerate(WAVE_AMPS)])


def synth(mean: float, coef: np.ndarray, n: int) -> np.ndarray:
    t = np.arange(n) / n
    k = np.arange(1, len(coef) + 1)
    return mean + np.real(np.exp(2j * math.pi * np.outer(t, k)) @ coef)


def incremental_resistance(solve_at_flow, geom, blood, flow, q, p_in, branch_q) -> float:
    """−∂P_outlet/∂Q at fixed inlet pressure, branch flows held at branch_q."""
    e = 1e-3 * q
    lo = solve_at_flow(geom, blood, flow, q - e, p_in, branch_q=branch_q).pressure[-1]
    hi = solve_at_flow(geom, blood, flow, q + e, p_in, branch_q=branch_q).pressure[-1]
    return (lo - hi) / (2 * e)


def inertance(geom, radius, rho) -> float:
    s = geom.s
    inv_a = 1.0 / (np.pi * radius**2)
    return float(rho * np.sum(0.5 * (inv_a[1:] + inv_a[:-1]) * np.diff(s)))


def pulsatile_summary(ctx: dict) -> dict:
    """ctx is assembled by hemodynamics.solve(); see the module docstring for the circuit."""
    from .hemodynamics import _solve_at_flow   # local import: hemodynamics imports this module

    fl, sol = ctx["flow"], ctx["sol"]
    to_ml = 6e7
    omega = 2 * math.pi * fl.heart_rate_bpm / 60.0
    coef = waveform_coefficients(fl.pulsatility_index)
    q_bar = sol.q
    out = {"heart_rate_bpm": fl.heart_rate_bpm, "PI_reference": fl.pulsatility_index}

    if fl.bc == "territory":
        r_d, r_dh = ctx["r_distal"], ctx["r_bed_healthy"]
        cap = fl.territory_tau_s / r_dh
        rv = incremental_resistance(_solve_at_flow, ctx["geom"], ctx["blood_eff"], fl, q_bar, ctx["p_in"], sol.branch_q)
        rv_h = incremental_resistance(_solve_at_flow, ctx["ref_geom"], ctx["blood_ref"], fl, ctx["q0"], ctx["p_cal"], ctx["ref_bq"])
        l_v, l_h = inertance(ctx["geom"], sol.radius, ctx["blood_eff"].rho), inertance(ctx["ref_geom"], ctx["ref"].radius, ctx["blood_ref"].rho)
        g_c = 1.0 / ctx["r_coll"] if sol.q_collateral > 0 else 0.0
        q_k, p_k, qt_k = [], [], []
        for k, c in enumerate(coef, start=1):
            w = k * omega
            zin_h = rv_h + 1j * w * l_h + 1.0 / (1.0 / r_dh + 1j * w * cap)
            p = ctx["q0"] * c * zin_h                     # inlet pressure pulse, fixed by the healthy calibration
            y = 1.0 / r_d + 1j * w * cap + g_c
            q = p / (rv + 1j * w * l_v + 1.0 / y)
            p_k.append(p); q_k.append(q); qt_k.append(q / y / r_d)
        q_k, p_k, qt_k = map(np.array, (q_k, p_k, qt_k))
        q_terr = sol.q_out + sol.q_collateral
        qi, qt, pi_ = synth(q_bar, q_k, N_SAMPLES), synth(q_terr, qt_k, N_SAMPLES), synth(ctx["p_in"], p_k, N_SAMPLES)
        out.update({
            "tau_s": fl.territory_tau_s, "incremental_resistance_ratio": rv / rv_h,
            "Q_territory_systolic_mL_min": float(qt.max()) * to_ml, "Q_territory_diastolic_mL_min": float(qt.min()) * to_ml,
            "PI_territory": float((qt.max() - qt.min()) / q_terr) if q_terr > 0 else math.nan,
            "inlet_pulse_pressure_mmHg": float((pi_.max() - pi_.min()) / 133.322),
        })
        wave_t = synth(q_terr, qt_k, N_OUT)
    else:
        q_k = q_bar * coef
        qi, pi_ = synth(q_bar, q_k, N_SAMPLES), np.full(N_SAMPLES, ctx["p_in"])
        wave_t = None

    out.update({
        "Q_inlet_systolic_mL_min": float(qi.max()) * to_ml, "Q_inlet_diastolic_mL_min": float(qi.min()) * to_ml,
        "PI_inlet": float((qi.max() - qi.min()) / q_bar),
        "flow_reversal": bool(qi.min() <= 0),
        "waveform_t": (np.arange(N_OUT) / N_OUT).tolist(),
        "waveform_Q_inlet_mL_min": (synth(q_bar, q_k, N_OUT) * to_ml).tolist(),
        "waveform_Q_territory_mL_min": None if wave_t is None else (wave_t * to_ml).tolist(),
    })
    # quasi-steady fields at peak systole and end diastole (full nonlinear solver)
    for tag, i in (("systolic", int(np.argmax(qi))), ("diastolic", int(np.argmin(qi)))):
        if qi[i] <= 0:
            continue
        s = _solve_at_flow(ctx["geom"], ctx["blood_eff"], fl, float(qi[i]), float(pi_[i]),
                           branch_r=ctx["branch_r"], p_out=ctx["p_out"])
        out[f"wss_max_{tag}_Pa"] = float(s.wss.max())
        out[f"wss_min_{tag}_Pa"] = float(s.wss.min())
        out[f"Re_max_{tag}"] = float(s.reynolds.max())
        out[f"dP_vessel_{tag}_mmHg"] = float((s.pressure[0] - s.pressure[-1]) / 133.322)
    return out
