"""Parametric 3D blood vessel model of the middle cerebral artery (M1) for ischemia studies."""
from .geometry import build_geometry, build_mesh
from .hemodynamics import flow_regime, solve
from .params import Blood, FlowConditions, Vessel

__all__ = ["Blood", "FlowConditions", "Vessel", "build_geometry", "build_mesh", "flow_regime", "solve"]
