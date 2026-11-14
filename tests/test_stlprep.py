"""Flow-extension tool: verified on synthetic open vessels with known answers."""
import numpy as np
import pytest

from bvm import Vessel, build_geometry, build_mesh
from bvm import stlprep as sp


def open_tube(**kw):
    """Our own tube mesh with the inlet/outlet caps removed: an open vessel with 2 circular ports."""
    m = build_mesh(build_geometry(Vessel(**kw)), n_theta=48)
    keep = m.face_tag == 0
    tris = m.vertices[m.faces[keep]] * 1e3   # mm
    return sp.weld(tris)


@pytest.mark.parametrize("kw", [{}, {"tortuosity": 1.2, "tortuosity_waves": 1.5}, {"occlusion_pct": 50, "taper": 0.8}])
def test_extension_of_open_tube_is_watertight_and_preserves_anatomy(kw):
    V, F = open_tube(**kw)
    assert sp.edge_stats(F)["boundary"] > 0
    ext = sp.extend(V, F, inlet_diameters=10, outlet_diameters=5)
    rep = sp.verify(V, F, ext)
    assert all(rep["checks"].values()), rep["checks"]
    assert ext.names == ["wall", "inlet", "outlet_1"]


def test_ports_found_with_correct_size_direction_and_length():
    V, F = open_tube(diameter_mm=3.0, length_mm=20.0, taper=0.8)
    ports = sp.find_ports(V, F)
    d = sorted(p.diameter for p in ports)
    # 48-gon area is slightly below the circle's
    assert d[1] == pytest.approx(3.0, rel=0.005) and d[0] == pytest.approx(2.4, rel=0.005)
    ext = sp.extend(V, F, inlet_diameters=10, outlet_diameters=5)
    inlet, outlet = ext.ports
    assert inlet.normal == pytest.approx([0, 0, -1], abs=1e-6)     # tube runs along +z: inlet faces -z
    assert outlet.normal == pytest.approx([0, 0, 1], abs=1e-6)
    assert ext.lengths == pytest.approx([10 * inlet.diameter, 5 * outlet.diameter])
    assert ext.end_centers[0][2] == pytest.approx(-10 * inlet.diameter, rel=1e-6)


def test_inlet_override_and_patch_counts(tmp_path):
    V, F = open_tube(taper=0.8)
    ext = sp.extend(V, F, inlet_index=[i for i, p in enumerate(sp.find_ports(V, F)) if p.diameter < 2.5][0])
    assert ext.ports[0].diameter < ext.ports[1].diameter
    sp.write_stl_patches(tmp_path / "x.stl", ext.V, ext.F, ext.tags, ext.names)
    text = (tmp_path / "x.stl").read_text()
    assert [ln.split()[1] for ln in text.splitlines() if ln.startswith("solid ")] == ["wall", "inlet", "outlet_1"]
    back = sp.weld(sp.read_stl(tmp_path / "x.stl"), tol=1e-5)
    assert sp.edge_stats(back[1]) == {"boundary": 0, "nonmanifold": 0, "misoriented": 0}
