"""Parametric 3D vessel geometry: centerline, radius profile, and surface mesh."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .params import Vessel

TAGS = ("wall", "inlet", "outlet")


@dataclass
class Centerline:
    s: np.ndarray          # arc length, m (N,)
    points: np.ndarray     # m (N, 3)
    tangent: np.ndarray    # unit (N, 3)
    normal: np.ndarray     # rotation-minimizing frame (N, 3)
    binormal: np.ndarray   # (N, 3)
    curvature: np.ndarray  # 1/m (N,)


@dataclass
class VesselGeometry:
    vessel: Vessel
    centerline: Centerline
    radius: np.ndarray       # unloaded lumen radius, m (N,)
    base_radius: np.ndarray  # healthy (pre-stenosis) radius incl. taper and bulb, m (N,)

    @property
    def s(self) -> np.ndarray:
        return self.centerline.s

    @property
    def throat_index(self) -> int:
        """Station of the stenosis throat (deepest narrowing relative to the healthy shape)."""
        if self.vessel.occlusion_pct > 0:
            return int(np.argmin(self.radius / self.base_radius))
        return int(np.argmin(self.radius))


def _sinusoid_arc_length(chord: float, amp: float, waves: float, n: int = 20000) -> float:
    z = np.linspace(0.0, chord, n)
    x = amp * np.sin(2 * np.pi * waves * z / chord)
    return float(np.sum(np.hypot(np.diff(x), np.diff(z))))


def _rotation_minimizing_frames(points: np.ndarray, tangent: np.ndarray) -> np.ndarray:
    """Double-reflection method (Wang et al., 2008)."""
    normal = np.zeros_like(points)
    t0 = tangent[0]
    ref = np.array([1.0, 0.0, 0.0]) if abs(t0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    n0 = ref - ref.dot(t0) * t0
    normal[0] = n0 / np.linalg.norm(n0)
    for i in range(len(points) - 1):
        v1 = points[i + 1] - points[i]
        c1 = v1 @ v1
        r_l = normal[i] - (2 / c1) * (v1 @ normal[i]) * v1
        t_l = tangent[i] - (2 / c1) * (v1 @ tangent[i]) * v1
        v2 = tangent[i + 1] - t_l
        c2 = v2 @ v2
        n_next = r_l - (2 / c2) * (v2 @ r_l) * v2 if c2 > 1e-30 else r_l
        normal[i + 1] = n_next / np.linalg.norm(n_next)
    return normal


def build_centerline(v: Vessel, n: int = 401) -> Centerline:
    """Centerline with fixed arc length; tortuosity shortens the chord.

    Tortuosity > 1 is a planar sinusoid whose amplitude is solved so that
    arc length / chord length equals the requested tortuosity index.
    """
    L = v.length_mm * 1e-3
    s = np.linspace(0.0, L, n)
    if v.tortuosity <= 1.0 + 1e-9:
        points = np.column_stack([np.zeros(n), np.zeros(n), s])
    else:
        chord = L / v.tortuosity
        amp = brentq(lambda a: _sinusoid_arc_length(chord, a, v.tortuosity_waves) - L, 0.0, L)
        z_d = np.linspace(0.0, chord, 20000)
        x_d = amp * np.sin(2 * np.pi * v.tortuosity_waves * z_d / chord)
        s_d = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x_d), np.diff(z_d)))])
        s_d *= L / s_d[-1]  # remove residual discretization error
        points = np.column_stack([np.interp(s, s_d, x_d), np.zeros(n), np.interp(s, s_d, z_d)])

    d = np.gradient(points, s, axis=0)
    tangent = d / np.linalg.norm(d, axis=1, keepdims=True)
    curvature = np.linalg.norm(np.gradient(tangent, s, axis=0), axis=1)
    normal = _rotation_minimizing_frames(points, tangent)
    binormal = np.cross(tangent, normal)
    return Centerline(s, points, tangent, normal, binormal, curvature)


def _cosine_bump(s: np.ndarray, center: float, half: float) -> np.ndarray:
    xi = (s - center) / half
    return np.where(np.abs(xi) < 1.0, 0.5 * (1.0 + np.cos(np.pi * np.clip(xi, -1, 1))), 0.0)


def base_radius_profile(v: Vessel, s: np.ndarray) -> np.ndarray:
    """Healthy lumen radius: inlet diameter, linear taper, optional cosine bulb (no stenosis)."""
    L = v.length_mm * 1e-3
    r = v.diameter_mm / 2 * 1e-3 * (1.0 + (v.taper - 1.0) * s / L)
    if v.bulb_pct > 0:
        r = r * (1.0 + v.bulb_pct / 100 * _cosine_bump(s, v.bulb_center * L, v.bulb_length_mm / 2 * 1e-3))
    return r


def radius_profile(v: Vessel, s: np.ndarray) -> np.ndarray:
    """Unloaded lumen radius: healthy shape with a smooth cosine-shaped stenosis.

    The stenosis percentage is relative to the local healthy diameter (in a bulb this
    is ECST-like, i.e. against the local bulb, not the distal ICA).
    """
    r = base_radius_profile(v, s)
    if v.occlusion_pct > 0:
        depth = 1.0 - v.throat_radius_ratio
        bump = _cosine_bump(s, v.stenosis_center * v.length_mm * 1e-3, v.stenosis_length_mm / 2 * 1e-3)
        r = r * (1.0 - depth * bump)
    return r


def build_geometry(v: Vessel, n: int = 401) -> VesselGeometry:
    cl = build_centerline(v, n)
    return VesselGeometry(v, cl, radius_profile(v, cl.s), base_radius_profile(v, cl.s))


@dataclass
class TubeMesh:
    vertices: np.ndarray   # (V, 3), m
    faces: np.ndarray      # (F, 3) int, outward-oriented
    face_tag: np.ndarray   # (F,) int index into TAGS
    vertex_ring: np.ndarray  # (V,) centerline index of each vertex (for mapping fields)

    def signed_volume(self) -> float:
        a, b, c = (self.vertices[self.faces[:, k]] for k in range(3))
        return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def build_mesh(geom: VesselGeometry, radius: np.ndarray | None = None, n_theta: int = 48) -> TubeMesh:
    """Closed, watertight triangulated lumen surface with wall/inlet/outlet patches."""
    cl = geom.centerline
    r = geom.radius if radius is None else radius
    ns = len(cl.s)
    th = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    rings = cl.points[:, None, :] + r[:, None, None] * (
        np.cos(th)[None, :, None] * cl.normal[:, None, :]
        + np.sin(th)[None, :, None] * cl.binormal[:, None, :]
    )
    verts = np.vstack([rings.reshape(-1, 3), cl.points[0], cl.points[-1]])
    c_in, c_out = ns * n_theta, ns * n_theta + 1

    i, j = np.meshgrid(np.arange(ns - 1), np.arange(n_theta), indexing="ij")
    a = i * n_theta + j
    b = i * n_theta + (j + 1) % n_theta
    c = (i + 1) * n_theta + j
    d = (i + 1) * n_theta + (j + 1) % n_theta
    wall = np.vstack([np.column_stack([a.ravel(), b.ravel(), d.ravel()]),
                      np.column_stack([a.ravel(), d.ravel(), c.ravel()])])
    jj = np.arange(n_theta)
    inlet = np.column_stack([np.full(n_theta, c_in), (jj + 1) % n_theta, jj])
    last = (ns - 1) * n_theta
    outlet = np.column_stack([np.full(n_theta, c_out), last + jj, last + (jj + 1) % n_theta])

    faces = np.vstack([wall, inlet, outlet])
    tags = np.concatenate([np.zeros(len(wall), int), np.ones(n_theta, int), np.full(n_theta, 2)])
    ring_idx = np.concatenate([np.repeat(np.arange(ns), n_theta), [0, ns - 1]])
    return TubeMesh(verts, faces, tags, ring_idx)
