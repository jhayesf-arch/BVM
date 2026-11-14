"""Prepare an open vessel STL for CFD: find the ports, add flow extensions, cap and name patches.

Typical input is a segmented lumen surface trimmed at its inlet and outlets (e.g. the CMHA healthy
MCA bifurcations, Song et al. 2024). For each open port the tool:

1. orders the boundary loop and fits its plane (SVD), with the normal pointing out of the vessel;
2. sweeps the loop outward along that normal, blending it into a circle of the same area over one
   diameter (smooth transition), then continues as a straight cylinder to the requested length;
3. caps the end with a flat disc, named ``inlet`` or ``outlet_k``.

Original vertices are never moved, so the segmented anatomy is preserved exactly. The largest port
is taken as the inlet unless ``inlet_index`` is given; outlets are numbered by decreasing diameter.

Only numpy and scipy are needed. Units are those of the input (the CMHA STLs are in mm).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


# ───────────────────────────── I/O ─────────────────────────────
def read_stl(path) -> np.ndarray:
    """Triangles as an (F, 3, 3) array (binary or ASCII STL)."""
    data = Path(path).read_bytes()
    if len(data) >= 84:
        n = struct.unpack("<I", data[80:84])[0]
        if 84 + 50 * n == len(data):
            rec = np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")])
            return np.frombuffer(data[84:], dtype=rec, count=n)["v"].reshape(-1, 3, 3).astype(float)
    tris, cur = [], []
    for line in data.decode("latin-1").splitlines():
        p = line.split()
        if p[:1] == ["vertex"]:
            cur.append([float(x) for x in p[1:4]])
            if len(cur) == 3:
                tris.append(cur)
                cur = []
    return np.array(tris, float)


def weld(tris: np.ndarray, tol: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """Merge coincident vertices. Returns (vertices, faces) with degenerate faces removed."""
    pts = tris.reshape(-1, 3)
    key = np.round(pts / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    V = pts[first]
    F = inv.reshape(-1, 3)
    keep = (F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])
    return V, F[keep]


def write_stl_patches(path, V: np.ndarray, F: np.ndarray, tags: np.ndarray, names: list[str], scale: float = 1.0):
    """ASCII STL with one named solid per patch (Fluent Meshing / snappyHexMesh read the names)."""
    out = []
    for k, name in enumerate(names):
        out.append(f"solid {name}")
        for f in F[tags == k]:
            a, b, c = V[f] * scale
            nrm = np.cross(b - a, c - a)
            nrm /= np.linalg.norm(nrm) or 1.0
            out.append(f" facet normal {nrm[0]:.6e} {nrm[1]:.6e} {nrm[2]:.6e}\n  outer loop")
            out += [f"   vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}" for p in (a, b, c)]
            out.append("  endloop\n endfacet")
        out.append(f"endsolid {name}")
    Path(path).write_text("\n".join(out) + "\n")


# ───────────────────────────── topology ─────────────────────────────
def edge_stats(F: np.ndarray) -> dict:
    """Counts of boundary / non-manifold / inconsistently oriented edges."""
    d = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    und = np.sort(d, axis=1)
    _, cnt = np.unique(und, axis=0, return_counts=True)
    _, dcnt = np.unique(d, axis=0, return_counts=True)
    return {"boundary": int(np.sum(cnt == 1)), "nonmanifold": int(np.sum(cnt > 2)),
            "misoriented": int(np.sum(dcnt > 1))}


def boundary_loops(F: np.ndarray) -> list[np.ndarray]:
    """Ordered boundary loops, each following the direction the edge has in its face."""
    d = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    present = {(int(a), int(b)) for a, b in d}
    nxt = {a: b for a, b in present if (b, a) not in present}
    loops, seen = [], set()
    for s in list(nxt):
        if s in seen:
            continue
        loop, x = [], s
        while x not in seen:
            seen.add(x)
            loop.append(x)
            x = nxt[x]
        loops.append(np.array(loop))
    return loops


def signed_volume(V: np.ndarray, F: np.ndarray) -> float:
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# ───────────────────────────── ports ─────────────────────────────
@dataclass
class Port:
    loop: np.ndarray          # vertex indices in boundary order
    center: np.ndarray
    normal: np.ndarray        # unit, pointing out of the vessel
    u: np.ndarray             # in-plane basis
    v: np.ndarray
    area: float
    planarity: float          # RMS out-of-plane distance / equivalent radius

    @property
    def diameter(self) -> float:
        return 2 * np.sqrt(self.area / np.pi)


def analyze_port(V: np.ndarray, loop: np.ndarray, interior: np.ndarray) -> Port:
    P = V[loop]
    c = P.mean(0)
    _, _, vt = np.linalg.svd(P - c, full_matrices=False)
    n = vt[2]
    # outward: away from the vessel material near the port
    near = interior[np.linalg.norm(interior - c, axis=1) < 3 * np.sqrt(np.ptp(P, axis=0).max() ** 2)]
    if len(near) and np.dot(n, c - near.mean(0)) < 0:
        n = -n
    u = vt[0]
    v = np.cross(n, u)
    xy = np.column_stack([(P - c) @ u, (P - c) @ v])
    area = 0.5 * abs(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))
    r_eq = np.sqrt(area / np.pi)
    planarity = float(np.sqrt(np.mean(((P - c) @ n) ** 2)) / r_eq)
    return Port(loop, c, n, u, v, float(area), planarity)


def find_ports(V: np.ndarray, F: np.ndarray) -> list[Port]:
    loops = boundary_loops(F)
    on_loop = np.zeros(len(V), bool)
    for lp in loops:
        on_loop[lp] = True
    interior = V[~on_loop]
    return [analyze_port(V, lp, interior) for lp in loops]


# ───────────────────────────── extensions ─────────────────────────────
def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3 - 2 * t)


@dataclass
class Extended:
    V: np.ndarray
    F: np.ndarray
    tags: np.ndarray                 # 0 = wall, 1 = inlet, 2.. = outlets
    names: list[str]
    ports: list[Port]                # in tag order (inlet first)
    lengths: list[float]             # extension length per port
    end_centers: list[np.ndarray] = field(default_factory=list)


def extend(V: np.ndarray, F: np.ndarray, inlet_diameters: float = 10.0, outlet_diameters: float = 5.0,
           inlet_index: int | None = None, transition_diameters: float = 1.0, spacing: float | None = None) -> Extended:
    """Add straight flow extensions (lengths in port diameters) and flat caps to every open port."""
    ports = find_ports(V, F)
    if len(ports) < 2:
        raise ValueError(f"expected an open vessel with at least 2 ports, found {len(ports)}")
    order = sorted(range(len(ports)), key=lambda i: -ports[i].diameter)
    if inlet_index is not None:
        order.remove(inlet_index)
        order.insert(0, inlet_index)
    ports = [ports[i] for i in order]
    if spacing is None:  # mean edge length of the input surface
        e = np.concatenate([V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 1]]])
        spacing = float(np.mean(np.linalg.norm(e, axis=1)))

    Vs, Fs, Ts = [V], [F], [np.zeros(len(F), int)]
    nv = len(V)
    lengths, ends = [], []
    for k, p in enumerate(ports):
        D = p.diameter
        L = (inlet_diameters if k == 0 else outlet_diameters) * D
        Lt = min(transition_diameters * D, L)
        R = D / 2
        P = V[p.loop]
        rel = P - p.center
        proj = rel - np.outer(rel @ p.normal, p.normal)              # loop projected onto its plane
        th = np.arctan2(rel @ p.v, rel @ p.u)
        circ = R * (np.outer(np.cos(th), p.u) + np.outer(np.sin(th), p.v))
        # axial stations: fine through the transition, then the same spacing on the straight part
        s = np.unique(np.concatenate([np.linspace(0, Lt, max(3, int(np.ceil(Lt / spacing)) + 1)),
                                      np.linspace(Lt, L, max(2, int(np.ceil((L - Lt) / (2 * spacing))) + 1))]))
        rings = [p.loop]
        for sk in s[1:]:
            w = _smoothstep(sk / Lt)
            ring = p.center + sk * p.normal + (1 - w) * proj + w * circ
            Vs.append(ring)
            rings.append(np.arange(nv, nv + len(ring)))
            nv += len(ring)
        m = len(p.loop)
        faces = []
        for r0, r1 in zip(rings[:-1], rings[1:]):
            for i in range(m):
                a, b, c, d = r0[i], r0[(i + 1) % m], r1[i], r1[(i + 1) % m]
                faces += [(b, a, c), (b, c, d)]   # reverses the boundary edge a→b, so orientation matches
        Fs.append(np.array(faces))
        Ts.append(np.zeros(len(faces), int))
        end = p.center + L * p.normal
        Vs.append(end[None])
        ci = nv
        nv += 1
        last = rings[-1]
        cap = [(last[(i + 1) % m], last[i], ci) for i in range(m)]
        Fs.append(np.array(cap))
        Ts.append(np.full(len(cap), k + 1))
        lengths.append(float(L))
        ends.append(end)
    names = ["wall", "inlet"] + [f"outlet_{i}" for i in range(1, len(ports))]
    return Extended(np.vstack(Vs), np.vstack(Fs), np.concatenate(Ts), names, ports, lengths, ends)


# ───────────────────────────── checks ─────────────────────────────
def verify(original_V: np.ndarray, original_F: np.ndarray, ext: Extended) -> dict:
    """Geometry verification of the prepared surface. Returns metrics and pass/fail per check."""
    st = edge_stats(ext.F)
    vol = signed_volume(ext.V, ext.F)
    # volume of the original with flat caps (fan from each port centre) for a consistency check
    capF, capV = [original_F], [original_V]
    nv = len(original_V)
    for p in ext.ports:
        capV.append(p.center[None])
        lp = p.loop
        capF.append(np.array([(lp[(i + 1) % len(lp)], lp[i], nv) for i in range(len(lp))]))
        nv += 1
    vol_orig = signed_volume(np.vstack(capV), np.vstack(capF))
    vol_ext_expected = sum(np.pi * (p.diameter / 2) ** 2 * L for p, L in zip(ext.ports, ext.lengths))
    vol_err = (vol - vol_orig - vol_ext_expected) / vol
    # extensions must stay clear of the rest of the vessel and of each other
    tree = cKDTree(original_V)
    clearance = []
    n0 = len(original_V)
    for p, L in zip(ext.ports, ext.lengths):
        axis = p.center + np.outer(np.linspace(p.diameter * 1.5, L, 40), p.normal)
        dist, _ = tree.query(axis)
        clearance.append(float(dist.min() / (p.diameter / 2)))   # in port radii; > 1 means no overlap
    axis_gap = []
    for i in range(len(ext.ports)):
        for j in range(i + 1, len(ext.ports)):
            a = ext.ports[i].center + np.outer(np.linspace(0, ext.lengths[i], 60), ext.ports[i].normal)
            b = ext.ports[j].center + np.outer(np.linspace(0, ext.lengths[j], 60), ext.ports[j].normal)
            gap = np.min(np.linalg.norm(a[:, None] - b[None], axis=2)) - (ext.ports[i].diameter + ext.ports[j].diameter) / 2
            axis_gap.append(float(gap))
    checks = {
        "watertight (no boundary edges)": st["boundary"] == 0,
        "manifold (no edge in > 2 faces)": st["nonmanifold"] == 0,
        "consistent orientation": st["misoriented"] == 0,
        "outward normals (volume > 0)": vol > 0,
        "original vertices unchanged": bool(np.array_equal(ext.V[:n0], original_V)),
        "volume = original + extensions (±2%)": abs(vol_err) < 0.02,
        "extensions clear of the vessel": all(c > 1.0 for c in clearance),
        "extensions clear of each other": all(g > 0 for g in axis_gap),
    }
    return {"edges": st, "volume": vol, "volume_original_capped": vol_orig,
            "volume_extensions_expected": vol_ext_expected, "volume_rel_error": float(vol_err),
            "clearance_radii": clearance, "extension_axis_gaps": axis_gap, "checks": checks}


# ───────────────────────────── one-call preparation ─────────────────────────────
def prepare(src, out_dir, inlet_diameters: float = 10.0, outlet_diameters: float = 5.0,
            inlet_index: int | None = None, tag: str | None = None) -> dict:
    """Extend, verify and write a CFD-ready surface. Returns the verification report dict."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = tag or Path(src).stem
    V, F = weld(read_stl(src))
    ext = extend(V, F, inlet_diameters, outlet_diameters, inlet_index)
    rep = verify(V, F, ext)
    stem = f"{tag}_ext_in{inlet_diameters:g}D_out{outlet_diameters:g}D"
    write_stl_patches(out / f"{stem}_mm.stl", ext.V, ext.F, ext.tags, ext.names)
    write_stl_patches(out / f"{stem}_m.stl", ext.V, ext.F, ext.tags, ext.names, scale=1e-3)
    for k, name in enumerate(ext.names):
        sel = ext.tags == k
        write_stl_patches(out / f"{stem}_{name}_mm.stl", ext.V, ext.F[sel], np.zeros(sel.sum(), int), [name])
    rep["ports"] = [{"name": n, "center": p.center.tolist(), "normal": p.normal.tolist(), "diameter": p.diameter,
                     "area": p.area, "planarity": p.planarity, "extension_length": L, "end_center": e.tolist()}
                    for n, p, L, e in zip(ext.names[1:], ext.ports, ext.lengths, ext.end_centers)]
    rep["files"] = sorted(str(p.name) for p in out.glob(f"{stem}*.stl"))
    rep["counts"] = {"vertices": len(ext.V), "faces": len(ext.F), "input_faces": len(F)}
    return rep


# ───────────────────────────── centerline, slicing, lesions ─────────────────────────────
def march_centerline(V: np.ndarray, port: Port, stop: np.ndarray, step: float | None = None,
                     smooth: int = 5) -> np.ndarray:
    """Centerline from a port into the vessel until it comes within one port radius of ``stop``.

    Each step moves along the current direction and re-centres on the centroid of the surface
    vertices in a thin slab perpendicular to it (restricted to 1.5 radii, so branches are ignored).
    """
    R = port.diameter / 2
    h = step or R / 4
    p, d = port.center.copy(), -port.normal
    pts = [p.copy()]
    for _ in range(2000):
        q = p + h * d
        rel = V - q
        along = rel @ d
        radial = np.linalg.norm(rel - np.outer(along, d), axis=1)
        sel = (np.abs(along) < h) & (radial < 1.5 * R)
        if sel.sum() < 6:
            break
        c = V[sel].mean(0)
        c = c - ((c - q) @ d) * d                 # keep the step length along d
        d_new = c - p
        d = d_new / np.linalg.norm(d_new)
        p = c
        pts.append(p.copy())
        if np.linalg.norm(p - stop) < R:
            break
    P = np.array(pts)
    if smooth > 1 and len(P) > smooth:
        k = np.ones(smooth) / smooth
        P[smooth // 2:-(smooth // 2)] = np.column_stack([np.convolve(P[:, i], k, mode="valid") for i in range(3)])
    return P


def resample(P: np.ndarray, ds: float) -> tuple[np.ndarray, np.ndarray]:
    """Uniform arc-length resampling. Returns (points, arc length)."""
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    t = np.arange(0, s[-1] + 1e-12, ds)
    return np.column_stack([np.interp(t, s, P[:, i]) for i in range(3)]), t


def slice_area(V: np.ndarray, F: np.ndarray, point: np.ndarray, normal: np.ndarray) -> float:
    """Lumen cross-section area cut by a plane: area of the closed section loop nearest ``point``."""
    dist = (V - point) @ normal
    sgn = dist > 0
    cross = F[(sgn[F[:, 0]] != sgn[F[:, 1]]) | (sgn[F[:, 1]] != sgn[F[:, 2]])]
    seg = []
    for f in cross:
        ends = []
        for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            if sgn[a] != sgn[b]:
                t = dist[a] / (dist[a] - dist[b])
                ends.append(((min(a, b), max(a, b)), V[a] + t * (V[b] - V[a])))
        if len(ends) == 2:
            seg.append(ends)
    if not seg:
        return 0.0
    adj, pos = {}, {}
    for (k1, p1), (k2, p2) in seg:
        adj.setdefault(k1, []).append(k2)
        adj.setdefault(k2, []).append(k1)
        pos[k1], pos[k2] = p1, p2
    u = np.cross(normal, [1.0, 0, 0] if abs(normal[0]) < 0.9 else [0, 1.0, 0])
    u /= np.linalg.norm(u)
    w = np.cross(normal, u)
    best, seen = None, set()
    for start in adj:
        if start in seen:
            continue
        loop, prev, cur = [], None, start
        while cur not in seen:
            seen.add(cur)
            loop.append(pos[cur])
            nxt = [k for k in adj[cur] if k != prev and k not in seen]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
        L = np.array(loop)
        if len(L) < 3:
            continue
        d0 = np.linalg.norm(L.mean(0) - point)
        if best is None or d0 < best[0]:
            xy = np.column_stack([(L - point) @ u, (L - point) @ w])
            area = 0.5 * abs(np.dot(xy[:, 0], np.roll(xy[:, 1], -1)) - np.dot(xy[:, 1], np.roll(xy[:, 0], -1)))
            best = (d0, area)
    return float(best[1]) if best else 0.0


def area_profile(V, F, C, s) -> np.ndarray:
    T = np.gradient(C, s, axis=0)
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    return np.array([slice_area(V, F, c, t) for c, t in zip(C, T)])


def apply_stenosis(V: np.ndarray, C: np.ndarray, s: np.ndarray, center: float, length: float,
                   pct: float, basis: str = "diameter", max_radius: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Cosine-shaped concentric narrowing of the lumen about a centerline (same profile as the 1D model).

    Vertices are pulled toward their nearest centerline point by factor 1 − depth·bump(s). Only vertices
    within ``max_radius`` of the centerline and whose foot point lies inside the lesion are moved.
    Returns (new vertices, mask of moved vertices).
    """
    ratio = 1 - pct / 100 if basis == "diameter" else np.sqrt(1 - pct / 100)
    depth = 1 - ratio
    tree = cKDTree(C)
    dist, idx = tree.query(V)
    sv = s[idx]
    xi = (sv - center) / (length / 2)
    bump = np.where(np.abs(xi) < 1, 0.5 * (1 + np.cos(np.pi * np.clip(xi, -1, 1))), 0.0)
    if max_radius is not None:
        bump = np.where(dist < max_radius, bump, 0.0)
    moved = bump > 0
    Vn = V.copy()
    Vn[moved] = C[idx[moved]] + (1 - depth * bump[moved, None]) * (V[moved] - C[idx[moved]])
    return Vn, moved


def face_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
