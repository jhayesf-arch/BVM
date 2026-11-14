"""Matplotlib figures: 3D lumen colored by a field, axial profiles, and sweeps."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from .geometry import build_mesh  # noqa: E402
from .hemodynamics import RE_LAMINAR_MAX, VesselSolution  # noqa: E402
from .params import MMHG_TO_PA  # noqa: E402


def plot_3d(sol: VesselSolution, path: Path, field: str = "wss") -> None:
    values = {"wss": (sol.wss, "WSS (Pa)"), "pressure": (sol.pressure / MMHG_TO_PA, "Pressure (mmHg)"),
              "reynolds": (sol.reynolds, "Reynolds number")}[field]
    mesh = build_mesh(sol.geometry, sol.radius, n_theta=36)
    wall = mesh.faces[mesh.face_tag == 0]
    tris = mesh.vertices[wall] * 1e3
    face_val = values[0][mesh.vertex_ring[wall]].mean(axis=1)
    norm = plt.Normalize(face_val.min(), face_val.max() if face_val.max() > face_val.min() else face_val.min() + 1)
    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(tris, facecolors=plt.cm.viridis(norm(face_val)), linewidth=0))
    pts = mesh.vertices * 1e3
    mid, span = (pts.max(0) + pts.min(0)) / 2, (pts.max(0) - pts.min(0)).max() / 2
    for setter, c in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), mid):
        setter(c - span, c + span)
    ax.set_xlabel("x (mm)"), ax.set_ylabel("y (mm)"), ax.set_zlabel("z (mm)")
    ax.view_init(elev=15, azim=-60)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax, shrink=0.6, label=values[1])
    v = sol.vessel
    ax.set_title(f"D={v.diameter_mm} mm, TI={v.tortuosity}, occlusion={v.occlusion_pct}% ({v.occlusion_basis}), "
                 f"E={v.stiffness_kpa} kPa")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_profiles(sol: VesselSolution, path: Path) -> None:
    s = sol.s * 1e3
    fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    axs[0, 0].plot(s, 2e3 * sol.geometry.radius, "--", label="unloaded")
    axs[0, 0].plot(s, 2e3 * sol.radius, label="pressurized")
    axs[0, 0].set_ylabel("Lumen diameter (mm)"), axs[0, 0].legend()
    axs[0, 1].plot(s, sol.pressure / MMHG_TO_PA)
    axs[0, 1].set_ylabel("Pressure (mmHg)")
    axs[1, 0].plot(s, sol.reynolds)
    axs[1, 0].axhline(RE_LAMINAR_MAX, color="r", ls=":", label="laminar limit (2000)")
    axs[1, 0].set_ylabel("Reynolds number"), axs[1, 0].legend()
    axs[1, 1].plot(s, sol.wss)
    axs[1, 1].axhline(0.4, color="gray", ls=":", label="low-WSS 0.4 Pa")
    axs[1, 1].set_ylabel("WSS (Pa)"), axs[1, 1].legend()
    for ax in axs[1]:
        ax.set_xlabel("Arc length (mm)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_sweep(rows: list[dict], param: str, path: Path) -> None:
    x = np.array([r[param] for r in rows])
    keys = [("Re_throat", "Re at throat"), ("wss_max_Pa", "Max WSS (Pa)"),
            ("dP_vessel_mmHg", "ΔP across vessel (mmHg)"), ("Q_mL_min", "Flow (mL/min)")]
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.5))
    for ax, (k, label) in zip(axs, keys):
        ax.plot(x, [r[k] for r in rows], "o-")
        ax.set_xlabel(param), ax.set_ylabel(label)
    axs[0].axhline(RE_LAMINAR_MAX, color="r", ls=":")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
