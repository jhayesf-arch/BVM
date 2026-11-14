"""Blood composition and oxygen delivery.

Hematocrit → viscosity
    Relative apparent viscosity follows the in vivo law of Pries et al. (Circ Res 1994;75:904), as written
    in the NIH Physiome JSim model 0157 ("InVivoViscosity"). Only the *ratio* η(D, Hct)/η(D, 45%) is used,
    so at a hematocrit of 45% every viscosity equals the project's specified value (μ = 3.5 mPa·s) and all
    earlier results are unchanged. The ratio is applied twice:
      * in the modeled vessel, at its inlet diameter (large vessels: η depends on Hct only);
      * in the downstream microvascular bed, at an effective diameter (default 25 µm, assumed), scaling the
        territory and perforator resistances. Most resistance sits there, so this term sets how much flow
        changes with hematocrit. The Fåhræus–Lindqvist effect makes it weaker than in large vessels.
    No autoregulation: in vivo, arterioles also dilate when O₂ content falls, which this does not capture.

Oxygen delivery (Fick principle)
    CaO₂ = 1.34·Hb·SaO₂ + 0.003·PaO₂ (mL O₂/dL), Hb = Hct·MCHC.
    DO₂ = Q·CaO₂. Content is the same upstream and downstream of a blockage (conducting arteries consume
    negligible O₂), so only flow changes delivery.
    Territory CBF = CBF_normal × (flow / healthy reference flow), with the healthy reference calibrated at
    a normal hematocrit (45%), so anemia can raise CBF above normal.
    Required OEF = CMRO₂ / (CBF·CaO₂/100). Normal ≈ 0.3–0.45. Up to OEF_max (≈ 0.8) the tissue
    compensates ("misery perfusion", Baron 1981; Powers 1991 stage II). Above it, demand cannot be met.

Wall O₂ transfer index
    For high-Péclet transport (O₂ in arteries) the concentration boundary layer is thin and the local
    mass-transfer coefficient scales as k ∝ γ̇_w^{1/3} (Lévêque). The index (γ̇_w/γ̇_ref)^{1/3}, relative
    to the healthy inlet, marks where shear-driven O₂ supply to the vessel wall falls. A 1D model cannot
    resolve recirculation zones, where 3D CFD shows the largest deficits (Moore & Ethier 1997; Tarbell 2003).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

HCT_REF = 45.0   # %; the hematocrit at which Blood.mu is specified


def pries_relative_viscosity(d_um: float, hct_pct: float) -> float:
    """Relative apparent viscosity (to plasma), Pries et al. 1994 in vivo law. d in µm, hct in %."""
    d, h = float(d_um), hct_pct / 100.0
    eta45 = 6.0 * math.exp(-0.085 * d) + 3.2 - 2.44 * math.exp(-0.06 * d**0.645)
    w = 1.0 / (1.0 + 1e-11 * d**12)
    c = (0.8 + math.exp(-0.075 * d)) * (-1.0 + w) + w
    esl = (d / (d - 1.1)) ** 2
    return (1.0 + (eta45 - 1.0) * (((1.0 - h) ** c - 1.0) / ((1.0 - 0.45) ** c - 1.0)) * esl) * esl


def viscosity_ratio(d_um: float, hct_pct: float) -> float:
    """η(D, Hct) / η(D, 45%)."""
    return pries_relative_viscosity(d_um, hct_pct) / pries_relative_viscosity(d_um, HCT_REF)


@dataclass(frozen=True)
class Oxygen:
    """Inputs for oxygen content and tissue demand (defaults: healthy adult at rest)."""

    sao2: float = 0.98            # arterial saturation, fraction
    pao2_mmhg: float = 95.0
    mchc_g_dl: float = 33.0       # Hb = Hct × MCHC
    cmro2: float = 3.3            # mL O₂/100 g/min (Kety & Schmidt 1948)
    cbf_normal: float = 50.0      # mL/100 g/min at the healthy reference flow
    oef_max: float = 0.8          # maximum extraction fraction (≈ 0.8–0.9 in PET studies)

    def hemoglobin(self, hct_pct: float) -> float:
        """g/dL."""
        return hct_pct * self.mchc_g_dl / 100.0

    def content(self, hct_pct: float) -> float:
        """Arterial O₂ content, mL O₂/dL."""
        return 1.34 * self.hemoglobin(hct_pct) * self.sao2 + 0.003 * self.pao2_mmhg


def oef_band(oef: float, oef_max: float) -> str:
    if oef <= 0.45:
        return "normal"
    if oef <= oef_max:
        return "compensated"     # misery perfusion: viable only by raising extraction
    return "unmet"


def territory_oxygen(q_ml_min: float, q_ref_ml_min: float | None, hct_pct: float, oxy: Oxygen) -> dict:
    cao2 = oxy.content(hct_pct)
    out = {"DO2_mL_min": q_ml_min * cao2 / 100.0}
    if q_ref_ml_min:
        cbf = oxy.cbf_normal * q_ml_min / q_ref_ml_min
        tissue_do2 = cbf * cao2 / 100.0              # mL O₂/100 g/min
        oef = oxy.cmro2 / tissue_do2 if tissue_do2 > 0 else math.inf
        out.update({"CBF_mL_100g_min": cbf, "DO2_tissue_mL_100g_min": tissue_do2,
                    "required_OEF": oef, "band": oef_band(oef, oxy.oef_max)})
    return out


def wall_o2_index(wss: np.ndarray, mu: float, wss_ref: float, mu_ref: float) -> np.ndarray:
    """Lévêque scaling of wall mass transfer with wall shear rate, relative to the healthy inlet."""
    return np.cbrt((np.asarray(wss) / mu) / (wss_ref / mu_ref))
