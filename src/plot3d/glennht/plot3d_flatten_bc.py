# plot3d_flatten_bc.py
"""plot3d-flatten deck boundary-condition dataclasses and YAML fragment
writer.

These mirror the ``boundary_conditions:`` list schema consumed by a
plot3d-flatten deck run configuration, covering scenarios such as:

- a single mixing-plane rotor/stator pair (``blades_full_ring`` per BC);
- a multi-row connected solve (surface-id convention ``100*row + {1..5}``,
  chained ``mixing_plane_partner`` seams across rows).

Only fields explicitly set (not ``None``) are emitted, so a plain wall BC
stays a 3-4 line mapping instead of a wall of ``null``s.

All dataclasses are keyword-only (``@dataclass(kw_only=True)``) so field
order in the class body can follow the schema's natural reading order
without fighting Python's "non-default argument follows default argument"
rule -- every field here has a name and every call site is expected to use
keyword arguments (matching how these run configurations are written by
hand).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import yaml

__all__ = [
    "Plot3DFlattenInletBC",
    "Plot3DFlattenOutletBC",
    "Plot3DFlattenWallBC",
    "write_boundary_conditions_yaml",
]


class _BCDumper(yaml.SafeDumper):
    """Dumper used only by this module, so the list-style override below
    can't leak into unrelated ``yaml.safe_dump()`` calls elsewhere in the
    process (PyYAML representers are registered per-``Dumper``-subclass,
    so subclassing rather than patching ``yaml.SafeDumper`` directly keeps
    this local)."""


def _represent_list(dumper: yaml.Dumper, data: list):
    """Inline (``[1, 2]``) for a list of plain scalars -- e.g. ``surfaces:``
    or a ``rotation.per_block[].blocks:`` entry -- block-style (one item per
    line) for anything containing a nested mapping or sequence, where an
    inline dump would be unreadable."""
    flow = all(isinstance(x, (int, float, str, bool, type(None))) for x in data)
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=flow)


_BCDumper.add_representer(list, _represent_list)


@dataclass(kw_only=True)
class Plot3DFlattenInletBC:
    """``type: inlet`` boundary condition.

    Args:
        name (str): BC name, e.g. ``"m_inlet"`` or ``"stator35_in"``.
        surfaces (List[int]): Surface id(s) this BC applies to.
        total_pressure (float): Total pressure, Pa.
        total_temperature (float): Total temperature, K.
        type (str): Always ``"inlet"``.
        flow_angle_deg (float, optional): Flow angle, degrees.
        turbulence_intensity (float, optional): Fractional turbulence
            intensity, e.g. ``0.02``.
        turbulence_length_scale (float, optional): Turbulence length
            scale, m.
        blades_full_ring (int, optional): Blade count for a full annulus
            (drives the mixing-plane pitch-change scaling).
        inlet_subtype (str, optional): e.g. ``"mixing_plane"``.
        mixing_plane_partner (str, optional): Name of the paired outlet
            BC on the upstream row.
        extra (Dict[str, Any], optional): Escape hatch for deck-specific
            fields not modelled above (e.g. ``mach`` in
            ``eee_4row.yaml``'s downstream mixing-plane inlets). Merged
            into the emitted mapping verbatim, after the named fields.
    """

    name: str
    surfaces: List[int]
    total_pressure: float
    total_temperature: float
    type: str = "inlet"
    flow_angle_deg: Optional[float] = None
    turbulence_intensity: Optional[float] = None
    turbulence_length_scale: Optional[float] = None
    blades_full_ring: Optional[int] = None
    inlet_subtype: Optional[str] = None
    mixing_plane_partner: Optional[str] = None
    extra: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "name": self.name,
            "type": self.type,
            "surfaces": list(self.surfaces),
            "total_pressure": self.total_pressure,
            "total_temperature": self.total_temperature,
        }
        if self.flow_angle_deg is not None:
            d["flow_angle_deg"] = self.flow_angle_deg
        if self.turbulence_intensity is not None:
            d["turbulence_intensity"] = self.turbulence_intensity
        if self.turbulence_length_scale is not None:
            d["turbulence_length_scale"] = self.turbulence_length_scale
        if self.blades_full_ring is not None:
            d["blades_full_ring"] = self.blades_full_ring
        if self.inlet_subtype is not None:
            d["inlet_subtype"] = self.inlet_subtype
        if self.mixing_plane_partner is not None:
            d["mixing_plane_partner"] = self.mixing_plane_partner
        if self.extra:
            d.update(self.extra)
        return d


@dataclass(kw_only=True)
class Plot3DFlattenOutletBC:
    """``type: outlet`` boundary condition.

    Args:
        name (str): BC name, e.g. ``"m_outlet"`` or ``"rotor35_out"``.
        surfaces (List[int]): Surface id(s) this BC applies to.
        back_pressure (float): Static back pressure, Pa.
        type (str): Always ``"outlet"``.
        extrapolation_order (int, optional): Extrapolation order at the
            outlet (e.g. ``0`` for the machine exit).
        subtype (str, optional): Outlet subtype, if any.
        mixing_plane_partner (str, optional): Name of the paired inlet BC
            on the downstream row.
        blades_full_ring (int, optional): Blade count for a full annulus.
        extra (Dict[str, Any], optional): See :class:`Plot3DFlattenInletBC`.
    """

    name: str
    surfaces: List[int]
    back_pressure: float
    type: str = "outlet"
    extrapolation_order: Optional[int] = None
    subtype: Optional[str] = None
    mixing_plane_partner: Optional[str] = None
    blades_full_ring: Optional[int] = None
    extra: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "name": self.name,
            "type": self.type,
            "surfaces": list(self.surfaces),
            "back_pressure": self.back_pressure,
        }
        if self.extrapolation_order is not None:
            d["extrapolation_order"] = self.extrapolation_order
        if self.subtype is not None:
            d["subtype"] = self.subtype
        if self.mixing_plane_partner is not None:
            d["mixing_plane_partner"] = self.mixing_plane_partner
        if self.blades_full_ring is not None:
            d["blades_full_ring"] = self.blades_full_ring
        if self.extra:
            d.update(self.extra)
        return d


@dataclass(kw_only=True)
class Plot3DFlattenWallBC:
    """``type: wall`` boundary condition.

    Args:
        name (str): BC name, e.g. ``"rotor35_blade"``.
        surfaces (List[int]): Surface id(s) this BC applies to.
        type (str): Always ``"wall"``.
        thermal (str): ``"adiabatic"`` (default) or an isothermal /
            heat-flux mode selected by ``wall_temperature`` /
            ``wall_heat_flux``.
        wall_temperature (float, optional): Fixed wall temperature, K.
        wall_heat_flux (float, optional): Fixed wall heat flux, W/m^2.
        rotating (bool, optional): Whether the wall spins with the row's
            frame (``rotation.per_block`` omega).
        wall_rotation_rate (float, optional): Explicit wall rotation rate
            override, rad/s, if different from the block's frame omega.
        extra (Dict[str, Any], optional): See :class:`Plot3DFlattenInletBC`.
    """

    name: str
    surfaces: List[int]
    type: str = "wall"
    thermal: str = "adiabatic"
    wall_temperature: Optional[float] = None
    wall_heat_flux: Optional[float] = None
    rotating: Optional[bool] = None
    wall_rotation_rate: Optional[float] = None
    extra: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "name": self.name,
            "type": self.type,
            "surfaces": list(self.surfaces),
            "thermal": self.thermal,
        }
        if self.wall_temperature is not None:
            d["wall_temperature"] = self.wall_temperature
        if self.wall_heat_flux is not None:
            d["wall_heat_flux"] = self.wall_heat_flux
        if self.rotating is not None:
            d["rotating"] = self.rotating
        if self.wall_rotation_rate is not None:
            d["wall_rotation_rate"] = self.wall_rotation_rate
        if self.extra:
            d.update(self.extra)
        return d


Plot3DFlattenBC = Any  # Plot3DFlattenInletBC | Plot3DFlattenOutletBC | Plot3DFlattenWallBC, kept loose for isinstance-free duck typing


def _boundary_conditions_yaml_text(
    bcs: List[Plot3DFlattenBC],
    *,
    rotation: Optional[Dict[str, Any]] = None,
) -> str:
    """Build the run-yaml fragment text for a list of GPU BC dataclasses.

    Factored out of :func:`write_boundary_conditions_yaml` so callers that
    need the YAML text in memory (e.g. to embed it as a string alongside
    other serialized data) don't have to round-trip through a temporary
    file just to get the same bytes :func:`write_boundary_conditions_yaml`
    would have written.

    Args:
        bcs (List[Plot3DFlattenInletBC | Plot3DFlattenOutletBC | Plot3DFlattenWallBC]): The boundary
            conditions to emit, in order.
        rotation (Dict[str, Any], optional): See
            :func:`write_boundary_conditions_yaml`.

    Returns:
        (str): The YAML text.
    """
    payload: Dict[str, Any] = {}
    if rotation is not None:
        payload["rotation"] = rotation
    payload["boundary_conditions"] = [bc.to_dict() for bc in bcs]

    return yaml.dump(payload, Dumper=_BCDumper, sort_keys=False, default_flow_style=False)


def write_boundary_conditions_yaml(
    bcs: List[Plot3DFlattenBC],
    filename: str,
    *,
    rotation: Optional[Dict[str, Any]] = None,
) -> str:
    """Serialize a list of GPU BC dataclasses to a run-yaml fragment.

    Args:
        bcs (List[Plot3DFlattenInletBC | Plot3DFlattenOutletBC | Plot3DFlattenWallBC]): The boundary
            conditions to emit, in order.
        filename (str): Output path, e.g. ``"stator_boundary_conditions.yaml"``.
        rotation (Dict[str, Any], optional): If given, also emit a
            top-level ``rotation:`` block (e.g. ``{"omega_x": 0.0,
            "per_block": [{"blocks": [0, 1], "omega_x": -1792.7}]}``),
            matching the ``rotation:`` block used in plot3d-flatten deck
            run configurations.

    Returns:
        (str): The YAML text that was written (useful for tests/inspection).
    """
    text = _boundary_conditions_yaml_text(bcs, rotation=rotation)

    with open(filename, "w") as fh:
        fh.write(text)

    return text
