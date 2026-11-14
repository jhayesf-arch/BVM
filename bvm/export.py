"""Writers for STL (CFD/3D printing), VTK (ParaView), CSV and JSON."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .geometry import TAGS, TubeMesh, build_mesh
from .hemodynamics import VesselSolution
from .params import MMHG_TO_PA


def write_stl(mesh: TubeMesh, path: Path, scale: float = 1e3) -> None:
    """ASCII STL with one named solid per patch (wall/inlet/outlet).

    Default scale writes millimetres. For OpenFOAM (metres) pass scale=1.
    Named solids become separate patches in snappyHexMesh / cfMesh.
    """
    v = mesh.vertices * scale
    lines = []
    for k, tag in enumerate(TAGS):
        lines.append(f"solid {tag}")
        for f in mesh.faces[mesh.face_tag == k]:
            a, b, c = v[f]
            n = np.cross(b - a, c - a)
            n /= np.linalg.norm(n) or 1.0
            lines.append(f" facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop")
            for p in (a, b, c):
                lines.append(f"   vertex {p[0]:.6e} {p[1]:.6e} {p[2]:.6e}")
            lines.append("  endloop\n endfacet")
        lines.append(f"endsolid {tag}")
    Path(path).write_text("\n".join(lines) + "\n")


def write_vtk(mesh: TubeMesh, path: Path, point_data: dict[str, np.ndarray], scale: float = 1e3) -> None:
    """Legacy VTK polydata; ``point_data`` arrays are per centerline station."""
    v = mesh.vertices * scale
    out = ["# vtk DataFile Version 3.0", "BVM vessel", "ASCII", "DATASET POLYDATA",
           f"POINTS {len(v)} float"]
    out += [f"{p[0]:.6e} {p[1]:.6e} {p[2]:.6e}" for p in v]
    out.append(f"POLYGONS {len(mesh.faces)} {4 * len(mesh.faces)}")
    out += [f"3 {a} {b} {c}" for a, b, c in mesh.faces]
    out.append(f"CELL_DATA {len(mesh.faces)}")
    out += ["SCALARS patch int 1", "LOOKUP_TABLE default"] + [str(t) for t in mesh.face_tag]
    out.append(f"POINT_DATA {len(v)}")
    for name, arr in point_data.items():
        out += [f"SCALARS {name} float 1", "LOOKUP_TABLE default"]
        out += [f"{x:.6e}" for x in np.asarray(arr)[mesh.vertex_ring]]
    Path(path).write_text("\n".join(out) + "\n")


def profile_table(sol: VesselSolution) -> dict[str, np.ndarray]:
    cl = sol.geometry.centerline
    return {
        "s_mm": sol.s * 1e3,
        "x_mm": cl.points[:, 0] * 1e3,
        "y_mm": cl.points[:, 1] * 1e3,
        "z_mm": cl.points[:, 2] * 1e3,
        "diameter_unloaded_mm": 2e3 * sol.geometry.radius,
        "diameter_loaded_mm": 2e3 * sol.radius,
        "curvature_1_per_mm": cl.curvature * 1e-3,
        "pressure_mmHg": sol.pressure / MMHG_TO_PA,
        "velocity_m_s": sol.velocity,
        "reynolds": sol.reynolds,
        "dean": sol.dean,
        "wss_Pa": sol.wss,
        "wall_strain": sol.strain,
    }


def write_csv(table: dict[str, np.ndarray], path: Path) -> None:
    cols = list(table)
    data = np.column_stack([table[c] for c in cols])
    np.savetxt(path, data, delimiter=",", header=",".join(cols), comments="", fmt="%.6g")


def write_all(sol: VesselSolution, out_dir: Path, n_theta: int = 48) -> dict[str, Path]:
    """Write loaded geometry (STL/VTK), axial profiles (CSV) and summary (JSON)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mesh_loaded = build_mesh(sol.geometry, sol.radius, n_theta)
    mesh_unloaded = build_mesh(sol.geometry, None, n_theta)
    table = profile_table(sol)
    paths = {
        "stl_unloaded": out_dir / "vessel_unloaded_mm.stl",
        "stl_loaded": out_dir / "vessel_loaded_mm.stl",
        "vtk": out_dir / "vessel.vtk",
        "csv": out_dir / "profiles.csv",
        "json": out_dir / "summary.json",
    }
    write_stl(mesh_unloaded, paths["stl_unloaded"])
    write_stl(mesh_loaded, paths["stl_loaded"])
    write_vtk(mesh_loaded, paths["vtk"], {k: v for k, v in table.items() if not k.endswith("_mm") or k.startswith("diameter")})
    write_csv(table, paths["csv"])
    paths["json"].write_text(json.dumps(sol.summary(), indent=2, default=float))
    return paths
