# glennht_export_functions.py
from __future__ import annotations
from collections import defaultdict
import math
import json
import pathlib
from dataclasses import asdict, is_dataclass
from typing import Any, Iterable, Tuple, Dict, List
import os 

from .class_definitions import (
    Job, BCGroup, GIF, VolumeZone,
)

# ============================================================
# Physics + units helpers
# ============================================================
_AIR_MOLW_DEFAULT = 28.964                 # kg/kmol
_R_UNIV = 8314.4126                        # J/(kmol*K)
_SUTH = dict(mu0=1.716e-5, T0=273.15, S=110.4)  # air

_PRESSURE_TO_PA = {
    "pa": 1.0,
    "kpa": 1e3,
    "mpa": 1e6,
    "bar": 1e5,
    "mbar": 1e2,
    "atm": 101325.0,
    "psi": 6894.757293168,
}

def to_pa(value: float | None, unit: str | None) -> float | None:
    if value is None:
        return None
    u = (unit or "Pa").strip().lower()
    factor = _PRESSURE_TO_PA.get(u)
    if factor is None:
        raise ValueError(
            f"Unknown pressure unit '{unit}'. "
            f"Use one of: {', '.join(sorted(_PRESSURE_TO_PA))}"
        )
    return float(value) * factor

def ideal_R(molw: float) -> float:
    return _R_UNIV / molw

def mach_from_p0_over_p(p0_over_p: float, gamma: float) -> float:
    if p0_over_p is None or p0_over_p <= 1.0:
        return 0.0
    exp = (gamma - 1.0) / gamma
    term = p0_over_p ** exp - 1.0
    return math.sqrt(max(0.0, 2.0 * term / (gamma - 1.0)))

def T_from_T0(T0: float, M: float, gamma: float) -> float:
    return T0 / (1.0 + 0.5 * (gamma - 1.0) * M * M)

def mu_suth(T: float) -> float:
    mu0, T0, S = _SUTH["mu0"], _SUTH["T0"], _SUTH["S"]
    return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)

def a_sound(T: float, gamma: float, R: float) -> float:
    return math.sqrt(gamma * R * T)

def cp_from_gamma_R(gamma: float, R: float) -> float:
    return gamma * R / (gamma - 1.0)

def pr_and_k(mu: float, cp: float, pr_given: float | None, k_given: float | None) -> tuple[float, float]:
    # Pr = mu * cp / k
    if pr_given is not None and k_given is None:
        return pr_given, (mu * cp) / pr_given
    if pr_given is None and k_given not in (None, 0):
        return (mu * cp) / k_given, k_given
    pr_default = 0.706
    return pr_default, (mu * cp) / pr_default

def rpm_to_omegab(rpm: float | None) -> float:
    if rpm in (None, 0):
        return 0.0
    return 2.0 * math.pi * (rpm / 60.0)

# ============================================================
# Small utilities
# ============================================================
def ensure_extension(path: str | pathlib.Path, ext: str) -> str:
    p = pathlib.Path(path)
    return str(p if p.suffix.lower() == ext.lower() else p.with_suffix(ext))

def _asdict_soft(x):
    if is_dataclass(x):
        return asdict(x)
    if isinstance(x, dict):
        return x
    return {k: v for k, v in x.__dict__.items() if not k.startswith("_")}

# ============================================================
# Formatting + namelist helpers
# ============================================================
def _iter_fields(obj: Any) -> Iterable[Tuple[str, Any]]:
    if is_dataclass(obj):
        for k, v in asdict(obj).items():
            if v is not None:
                yield k, v
    elif isinstance(obj, dict):
        for k, v in obj.items():
            if v is not None:
                yield k, v
    else:
        for k, v in obj.__dict__.items():
            if not k.startswith("_") and v is not None:
                yield k, v

def _fmt_bool(v: bool) -> str:
    # For the detailed blocks, your style uses .TRUE./.FALSE.
    return ".TRUE." if v else ".FALSE."

def _fmt_value(v: Any) -> str:
    from enum import Enum, IntEnum
    if isinstance(v, bool):
        return _fmt_bool(v)
    if isinstance(v, (IntEnum, Enum)):
        return str(int(v))
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, (list, tuple)):
        return ",".join(_fmt_value(x) for x in v if x is not None)
    if is_dataclass(v) or isinstance(v, dict):
        inner = ",".join(f"{k}={_fmt_value(val)}" for k, val in _iter_fields(v))
        return f"({inner})"
    return f"'{str(v)}'"

def _export_namelist_block(header: str, obj: Any, *, exclude_names=()) -> str:
    """
    Create a Fortran namelist block like:
      &HEADER
      key1=..., key2=...
      &END
    Skips None, excludes any field in exclude_names and *_unit.
    """
    pairs = []
    for k, v in _iter_fields(obj):
        if k in exclude_names or k.endswith("_unit"):
            continue
        pairs.append(f"{k}={_fmt_value(v)}")
    inner = ", ".join(pairs)
    return f" &{header}\n{inner}\n &END\n"

def _write_bsurf_spec(w, bc):
    line = (
        f" &BSurf_Spec\n"
        f"BSurfID={bc.SurfaceID}, BCType={int(bc.BCType)}, BSurfName='{bc.Name}'"
    )
    if getattr(bc, "IsPostProcessing", False):
        line += ", BRefCond=T"
    if getattr(bc, "IsCalculateMassFlow", False) or getattr(bc, "ToggleProcessSurface", False):
        line += ", BCalc=T"
    w.write(line + "\n &END\n\n")

# ============================================================
# GIF + Volume Zone writers (dict inputs supported)
# ============================================================
def _write_gif_from_dict(w, gdict: Dict[str, Any]) -> None:
    sid1 = int(gdict.get("a", 0))
    sid2 = int(gdict.get("b", 0))
    name1 = f"surface {sid1}"
    name2 = f"surface {sid2}"
    bctype = 4  # GIF
    w.write(
        f" &BSurf_Spec\nBSurfID={sid1}, BCType={bctype}, BSurfName='{name1}'\n &END\n\n"
    )
    w.write(
        f" &BSurf_Spec\nBSurfID={sid2}, BCType={bctype}, BSurfName='{name2}'\n &END\n\n"
    )
    w.write(f" &GIF_Spec\nSurfID_1={sid1}, SurfID2={sid2}\n &END\n\n")

def _write_vzconditions(w, vz: Dict[str, Any]) -> None:
    """
    Convert your dict style:
      {"block_index": 1, "zone_type": "fluid"|"solid", "contiguous_id": 1}
    to the GHT namelist you showed (defaults baked in).
    """
    vzid = int(vz.get("contiguous_id", 0))
    ztype = str(vz.get("zone_type", "fluid")).strip().lower()
    vztype = 1 if ztype == "fluid" else 2

    if vztype == 1:
        # Fluid default
        w.write(
            " &VZConditions\n"
            f"VZid={vzid}, VZtype=1, OmegaVZ=0., VZMaterialName=Air,\n"
            "Fluid_Tref_prop=0., Fluid_k_Tref=285., Fluid_amu_Tref=285., Fluid_expnt=.7,UseDryAir=.TRUE.,\n"
            "!Fluid_cp=1002., Fluid_Pr=.7, Fluid_MW=28.964\n"
            " &END\n\n"
        )
    else:
        # Solid default
        w.write(
            " &VZConditions\n"
            f"VZid={vzid}, VZtype=2, OmegaVZ=0., VZMaterialName=CMC,\n"
            "Solid_Tref_prop=285., Solid_rho_Tref=2707. , Solid_condN_Tref=6.5, Solid_condT_Tref=6.5, Solid_condA_Tref=6.5,\n"
            "Solid_Csp_Tref=896.\n"
            " &END\n\n"
        )

# ============================================================
# Reference population (physics-driven, builds RefCond from Full)
# ============================================================
def populate_reference_from_inputs(
    job: Job,
    *,
    inlet_total_P0_Pa: float | None,
    reference_static_P_Pa: float | None,
    inlet_T0_K: float | None,
    refLen_m: float | None,
    MolW: float | None = None,
    gamma: float | None = None,
    Pr: float | None = None,
    k_override_WmK: float | None = None,
    rpm: float | None = None,
    rho_solid: float | None = None,
):
    """Fill job.ReferenceCondFull from physics and then mirror into job.ReferenceCond."""
    from glennht_export_classes import ReferenceCond  # local import to avoid circular typing

    rcfull = job.ReferenceCondFull
    MolW  = MolW  if MolW  not in (None, 0) else _AIR_MOLW_DEFAULT
    gamma = gamma if gamma not in (None, 0) else 1.4
    R     = ideal_R(MolW)
    T0    = inlet_T0_K if inlet_T0_K not in (None, 0) else (rcfull.refT0 or 300.0)

    # totals / length
    rcfull.refP0 = float(inlet_total_P0_Pa) if inlet_total_P0_Pa not in (None, 0) else (rcfull.refP0 or 101325.0)
    rcfull.reflen = float(refLen_m) if refLen_m not in (None, 0) else (rcfull.reflen or 1.0)

    # Mach from p0/p
    if reference_static_P_Pa not in (None, 0) and rcfull.refP0 not in (None, 0):
        p0_over_p = rcfull.refP0 / reference_static_P_Pa
        M = mach_from_p0_over_p(p0_over_p, gamma)
    else:
        M = 0.0

    # static temp, speed of sound, velocity
    T_static = T_from_T0(T0, M, gamma)
    a = a_sound(T_static, gamma, R)
    refVel = M * a

    # viscosity, density
    mu = mu_suth(T_static)
    Pref_static = reference_static_P_Pa if reference_static_P_Pa not in (None, 0) else rcfull.refP0
    rho = Pref_static / (R * T_static)

    # cp, k, Pr
    cp = cp_from_gamma_R(gamma, R)
    Pr_val, k_val = pr_and_k(mu, cp, Pr, k_override_WmK)

    # Reynolds
    Re = (rho * refVel * rcfull.reflen / mu) if (refVel and mu and rcfull.reflen) else 0.0

    # fill Full
    rcfull.refT0   = T0
    rcfull.refrho0 = rho
    rcfull.refVel  = refVel
    rcfull.refvisc = mu
    rcfull.refcond = k_val
    rcfull.refCp   = cp
    rcfull.MolW    = MolW
    rcfull.RgasUnv = _R_UNIV
    rcfull.Rgas    = R
    rcfull.gamma   = gamma
    rcfull.Re      = Re
    rcfull.Pr      = Pr_val
    rcfull.ndVisc  = 1.0
    rcfull.ndCond  = 1.0
    rcfull.Omegab  = rpm_to_omegab(rpm)
    rcfull.ReScalingFactor = 1.0
    rcfull.rho_solid = rho_solid
    rcfull.cond_solid = 20.0
    rcfull.csp_solid  = 896.0

    # compact RefCond derived from Full (stored back on the job)
    rc = ReferenceCond(
        useDimensionalVariables=False,
        refLen=rcfull.reflen,
        refP0=rcfull.refP0,
        refT0=rcfull.refT0,
        refRho0=rcfull.refrho0,
        refVel=rcfull.refVel,
        refVisc=rcfull.refvisc,
        refCond=rcfull.refcond,
        refCp=rcfull.refCp,
        MolW=rcfull.MolW,
        RgasUnv=rcfull.RgasUnv,
        Rgas=rcfull.Rgas,
        gamma=rcfull.gamma,
        Re=rcfull.Re,
        Pr=rcfull.Pr,
        ndVisc=1.0,
        ndCond=1.0,
        Omegab=rcfull.Omegab,
        ReScalingFactor=1.0,
        rho_solid=rcfull.rho_solid,
        cond_solid=rcfull.cond_solid,
        Csp_solid=rcfull.csp_solid,
    )
    job.ReferenceCond = rc

# ============================================================
# Boundary conditions export (.bcs) — robust defaults
# ============================================================
def export_to_boundary_condition(
    file_path_to_write: str,
    job_settings: Job,
    bc_group: BCGroup,
    gif_pairs: List[GIF] | List[Dict[str,Any]],
    volume_zones: List[VolumeZone] | List[Dict[str,Any]],
):
    file_path_to_write = ensure_extension(file_path_to_write, '.bcs')
    path = pathlib.Path(file_path_to_write)
    json_path = path.with_suffix(".json")

    # Sidecar JSON for debugging
    json_payload = {
        "inlets": [_asdict_soft(x) for x in bc_group.Inlets],
        "outlets": [_asdict_soft(x) for x in bc_group.Outlets],
        "slips": [_asdict_soft(x) for x in bc_group.SymmetricSlips],
        "walls": [_asdict_soft(x) for x in bc_group.Walls],
        "volume_zones": [x for x in volume_zones],
        "gif_pairs": [_asdict_soft(x) for x in gif_pairs],
        "job_settings": _asdict_soft(job_settings),
    }
    json_path.write_text(json.dumps(json_payload, indent=4))

    # ---- Robust reference defaults (use highest-pressure inlet if available)
    ref = job_settings.ReferenceCondFull

    def _reference_inlet(inlets: List[Any]) -> Tuple[Any | None, float | None]:
        best = None
        best_pa: float | None = None
        for inlet in inlets:
            p0 = getattr(inlet, "P0_const", None)
            if p0 is None:
                continue
            phys_pa = to_pa(p0, getattr(inlet, "P0_const_unit", "Pa"))
            if phys_pa is None:
                continue
            if best_pa is None or phys_pa > best_pa:
                best = inlet
                best_pa = phys_pa
        if best is None and inlets:
            return inlets[0], None
        return best, best_pa

    ref_inlet, ref_inlet_pa = _reference_inlet(bc_group.Inlets)

    # refLen: default 1.0 if missing
    if getattr(ref, "reflen", None) in (None, 0):
        ref.reflen = 1.0  # type: ignore

    # refP0: from highest-pressure inlet (convert to Pa) if missing
    if getattr(ref, "refP0", None) in (None, 0):
        if ref_inlet_pa not in (None, 0):
            ref.refP0 = ref_inlet_pa  # type: ignore

    # refT0: from first inlet if missing
    if getattr(ref, "refT0", None) in (None, 0):
        if ref_inlet and getattr(ref_inlet, "T0_const", None) is not None:
            ref.refT0 = ref_inlet.T0_const  # type: ignore

    def _dedupe_by_bc_id(objs: Iterable[Any]) -> List[Any]:
        seen: set[int] = set()
        unique: List[Any] = []
        for obj in objs:
            sid = obj.get("id") if isinstance(obj, dict) else getattr(obj, "SurfaceID", None)
            if sid is None or sid in seen:
                continue
            seen.add(sid)
            unique.append(obj)
        return unique

    # ---- Write .bcs file
    with path.open("w", encoding="utf-8") as w:
        # INLETS (normalize to refP0, refT0, refLen)
        for inlet in _dedupe_by_bc_id(bc_group.Inlets):
            if getattr(inlet, "P0_const", None) is not None:
                phys_pa = to_pa(inlet.P0_const, getattr(inlet, "P0_const_unit", "Pa"))
                if phys_pa is not None and ref.refP0 not in (None, 0):
                    inlet.P0_const = phys_pa / ref.refP0
            if getattr(inlet, "T0_const", None) is not None and ref.refT0 not in (None, 0):
                inlet.T0_const = inlet.T0_const / ref.refT0  # type: ignore
            if getattr(inlet, "twall_hub", None) is not None and ref.refT0 not in (None, 0):
                inlet.twall_hub = inlet.twall_hub / ref.refT0  # type: ignore
            if getattr(inlet, "twall_case", None) is not None and ref.refT0 not in (None, 0):
                inlet.twall_case = inlet.twall_case / ref.refT0  # type: ignore
            if getattr(inlet, "Ts_const", None) is not None and ref.reflen not in (None, 0):
                inlet.Ts_const = inlet.Ts_const / ref.reflen  # type: ignore
            _write_bsurf_spec(w, inlet)

        # OUTLETS (normalize back-pressure by refP0)
        for outlet in _dedupe_by_bc_id(bc_group.Outlets):
            if getattr(outlet, "Pback_const", None) is not None and ref.refP0 not in (None, 0):
                phys_pa = to_pa(outlet.Pback_const, getattr(outlet, "Pback_const_unit", "Pa"))
                if phys_pa is not None:
                    outlet.Pback_const = phys_pa / ref.refP0
            _write_bsurf_spec(w, outlet)

        # SLIPS / WALLS
        for slip in _dedupe_by_bc_id(bc_group.SymmetricSlips):
            _write_bsurf_spec(w, slip)
        for wall in _dedupe_by_bc_id(bc_group.Walls):
            _write_bsurf_spec(w, wall)

        # GIFS (dicts or dataclasses)
        for pair in gif_pairs:
            if isinstance(pair, dict):
                _write_gif_from_dict(w, pair)
            else:
                name1 = getattr(pair, "Name1", f"surface {pair.GIFSurface1}")
                name2 = getattr(pair, "Name2", f"surface {pair.GIFSurface2}")
                w.write(
                    f" &BSurf_Spec\nBSurfID={pair.GIFSurface1}, BCType={int(pair.BCType)}, "
                    f"BSurfName='{name1}'\n &END\n\n"
                )
                w.write(
                    f" &BSurf_Spec\nBSurfID={pair.GIFSurface2}, BCType={int(pair.BCType)}, "
                    f"BSurfName='{name2}'\n &END\n\n"
                )
                w.write(f" &GIF_Spec\nSurfID_1={pair.GIFSurface1}, SurfID2={pair.GIFSurface2}\n &END\n\n")

        # VZConditions (dict templates)
        volume_zone_unique = {d["contiguous_index"]: d for d in volume_zones}.values()
        for vz in volume_zone_unique:
            if isinstance(vz, dict):
                _write_vzconditions(w, vz)
            else:
                w.write(_export_namelist_block("VZConditions", vz)); w.write("\n")

        # keep the *first* object for each unique obj.subtype
        def first_by_subtype(objs: Iterable[T], subtype_attr: str, *, include_none=False) -> List[T]:
            seen = set()
            out: List[T] = []
            for o in objs:
                # works for both objects and dicts
                st = getattr(o, subtype_attr, None) if not isinstance(o, dict) else o.get(subtype_attr)
                if st is None and not include_none:
                    continue
                if st not in seen:
                    seen.add(st)
                    out.append(o)
            return out
        
        inlet_bc = first_by_subtype(bc_group.Inlets,'inlet_subType')
        outlet_bc = first_by_subtype(bc_group.Outlets,'outlet_subType')
        slip_bc = first_by_subtype(bc_group.SymmetricSlips,'slip_subType')
        wall_bc = first_by_subtype(bc_group.Walls,'wall_subType')
        
        # Detailed BC blocks (skip meta + *_unit)
        exclude = {"Name", "SurfaceID", "BCType"}
        for inlet in inlet_bc:
            w.write(_export_namelist_block("INLET_BC", inlet, exclude_names=exclude)); w.write("\n")
        for outlet in outlet_bc:
            w.write(_export_namelist_block("OUTLET_BC", outlet, exclude_names=exclude)); w.write("\n")
        for slip in slip_bc:
            w.write(_export_namelist_block("SLIP_BC", slip, exclude_names=exclude)); w.write("\n")
        for wall in wall_bc:
            w.write(_export_namelist_block("WALL_BC", wall, exclude_names=exclude)); w.write("\n")

# ============================================================
# Job export (styled like your sample)
# ============================================================
def export_to_job_file(
    job: Job,
    file_path_to_write: str | pathlib.Path,
    *,
    title: str | None = None,
    exec_serial: str = "GlennHT.serial",
    exec_mpi: str | None = None,           # optional: also print mpi line
):
    path = pathlib.Path(file_path_to_write)
    with path.open("w", encoding="utf-8") as w:
        if title:
            w.write(" &Title\n")
            w.write(f' TheTitle="{title}"\n')
            w.write(" &end\n\n")

        # Derive compact RefCond from Full if not already present
        if job.ReferenceCond is None and job.ReferenceCondFull is not None:
            # Attempt a simple populate if Full has enough fields; otherwise leave as-is.
            try:
                populate_reference_from_inputs(
                    job,
                    inlet_total_P0_Pa = job.ReferenceCondFull.refP0,
                    reference_static_P_Pa = None,
                    inlet_T0_K = job.ReferenceCondFull.refT0,
                    refLen_m = job.ReferenceCondFull.reflen,
                )
            except Exception:
                # Best effort—job.ReferenceCond may remain None if insufficient data.
                pass

        w.write(_export_namelist_block("JobFiles", job.JobFiles)); w.write("\n")
        w.write(_export_namelist_block("JobControl", job.JobControl)); w.write("\n")
        w.write(_export_namelist_block("TurbModelInput", job.TurbModelInput)); w.write("\n")
        w.write(_export_namelist_block("Plot3DParameters", job.Plot3DParameters)); w.write("\n")
        w.write(_export_namelist_block("InitialCond", job.InitialCond)); w.write("\n")
        w.write(_export_namelist_block("TimeStpControl", job.TimeStpControl)); w.write("\n")
        w.write(_export_namelist_block("SPDSchemeControl", job.SPDSchemeControl)); w.write("\n")
        w.write(_export_namelist_block("RKSchemeControl", job.RKSchemeControl)); w.write("\n")
        w.write(_export_namelist_block("MGSchemeControl", job.MGSchemeControl)); w.write("\n")
        w.write(_export_namelist_block("GasPropertiesInput", job.GasPropertiesInput)); w.write("\n")

        # Prefer Full (authoritative); also print compact if present
        if job.ReferenceCond is not None:
            w.write(_export_namelist_block("ReferenceCond", job.ReferenceCond)); w.write("\n")
        if job.ReferenceCondFull is not None:
            w.write(_export_namelist_block("ReferenceCondFull", job.ReferenceCondFull)); w.write("\n")

        # Footer like your example
        w.write(f'execFILE="{exec_serial}"\n')
        if exec_mpi:
            w.write(f'execFILE="{exec_mpi}"\n')

def summarize_contiguous(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Given a list of dicts with keys: 'block_index', 'zone_type', 'contiguous_id',
    return:
      - num_unique_contiguous_ids: int
      - zone_types_by_id: {contiguous_id: [unique zone types]}
      - ids_with_multiple_zone_types: [ids that have >1 zone type]
    """
    id_to_zone_types = defaultdict(set)

    for r in records:
        # will raise KeyError if required keys are missing (helpful fail)
        cid = r["contiguous_index"]
        zt = r["zone_type"]
        id_to_zone_types[cid].add(zt)

    zone_types_by_id = {cid: sorted(zts) for cid, zts in id_to_zone_types.items()}
    ids_with_multiple = [cid for cid, zts in zone_types_by_id.items() if len(zts) > 1]

    return {
        "num_unique_contiguous_indices": len(zone_types_by_id),
        "zone_types_by_id": zone_types_by_id,
        "ids_with_multiple_zone_types": sorted(ids_with_multiple),
    }


def export_to_glennht_conn(matches:List[Dict[str, Dict[int, str]]],outer_faces:List[Dict[str,int]],filename:str, 
                           gif_pairs:List[List[Dict[str, int]]],gif_faces:List[List[Dict[str, int]]],
                               volume_zones:List[Dict[str,Any]]):
    """Exports the connectivity to GlennHT format 

    Args:
        matches (Dict[str,Dict[int,str]]): Any matching faces between blocks 
        outer_faces (Dict[str,int]): Non matching faces of all blocks or surfaces to consider 
    """
    lines = list()
    
    # Print face matches
    blocks = ['block1','block2']
    nMatches = len(matches)
    lines.append(f'{nMatches}\n') # Print number of matches 
    for match in matches:                        
        for block in blocks:
            block_indx = match[block]['block_index']+1 # type: ignore # block1 and block2 are arbitrary names, the key is the block index 
            block_IMIN = match[block]['IMIN']+1 # type: ignore
            block_JMIN = match[block]['JMIN']+1 # type: ignore
            block_KMIN = match[block]['KMIN']+1 # type: ignore

            block_IMAX = match[block]['IMAX']+1 # type: ignore
            block_JMAX = match[block]['JMAX']+1 # type: ignore
            block_KMAX = match[block]['KMAX']+1 # type: ignore

            lines.append(f"{block_indx:3d}\t{block_IMIN:5d} {block_JMIN:5d} {block_KMIN:5d}\t{block_IMAX:5d} {block_JMAX:5d} {block_KMAX:5d}\n")
    # Print Surfaces 
    # Get total number of surfaces
    outer_faces.extend(gif_faces)
    outer_faces = sorted(outer_faces, key=lambda d: d["id"])

    lines.append(f"{len(outer_faces)}\n")
    for surface in outer_faces:
        block_indx = surface['block_index']+1 # type: ignore
        IMIN = surface['IMIN']+1 # type: ignore
        JMIN = surface['JMIN']+1 # type: ignore
        KMIN = surface['KMIN']+1 # type: ignore
        
        IMAX = surface['IMAX']+1 # type: ignore
        JMAX = surface['JMAX']+1 # type: ignore
        KMAX = surface['KMAX']+1 # type: ignore
        id = surface['id']       # type: ignore
        lines.append(f"{block_indx:3d}\t{IMIN:5d} {JMIN:5d} {KMIN:5d}\t{IMAX:5d} {JMAX:5d} {KMAX:5d}\t{id:4d}\n")
    
    # Print GIFs
    lines.append(f"{len(gif_pairs)}\n")
    for pairs in gif_pairs:
        lines.append(f"{pairs['a']} {pairs['b']} -2 1\n") # type: ignore
    
    summary = summarize_contiguous(volume_zones)
    # Print volume zones
    
    lines.append(f"{summary['num_unique_contiguous_indices']}\n")
    # Print what types of zones are there 
    for k,v in summary['zone_types_by_id'].items():
        lines.append(f"{k} ")
    lines.append("\n")
    # Print Zone Groups
    columns_to_print = 10
    for i,v in enumerate(volume_zones):
        if i % columns_to_print==0:
            lines.append(f"{v["contiguous_index"]}\n")
        else:
            lines.append(f"{v["contiguous_index"]} ")
    
    filename = ensure_extension(filename,'.ght_conn')
    with open(f'{filename}','w') as fp:
        fp.writelines(lines)

def ensure_extension(filename, default_ext=".txt"):
    base, ext = os.path.splitext(filename)
    if not ext:  # no extension present
        return filename + default_ext
    return filename
# ============================================================
# ---- Quick self-test / example usage -----------------------
# ============================================================
if __name__ == "__main__":
    # Minimal smoke test (requires your glennht_export_classes.py definitions)
    from glennht_export_classes import (
        Job, InletBC, OutletBC, BoundaryConditionType, BCGroup,
    )

    job = Job()

    # Populate physics from simple inputs (e.g., 60 bar total at inlet, 2 bar static outlet)
    populate_reference_from_inputs(
        job,
        inlet_total_P0_Pa = 60.0 * 1e5,
        reference_static_P_Pa = 2.0 * 1e5,
        inlet_T0_K = 300.0,
        refLen_m = 0.0254,
        rpm=0.0,
        rho_solid=2707.0,
    )

    # BCs with physical inputs; exporter normalizes to reference
    inlet = InletBC(
        BCType=BoundaryConditionType.Inlet, SurfaceID=10, Name="Inlet-60bar",
        P0_const=60.0, P0_const_unit="bar", T0_const=300.0
    )
    outlet = OutletBC(
        BCType=BoundaryConditionType.Outlet, SurfaceID=20, Name="Outlet-2bar",
        Pback_const=2.0, Pback_const_unit="bar"
    )
    bcg = BCGroup(Inlets=[inlet], Outlets=[outlet], SymmetricSlips=[], Walls=[])

    gifs = [
        {"id1": 101, "id2": 202, "pty": 12},
    ]
    volume_zones = [
        {"block_index": 1, "zone_type": "fluid", "contiguous_id": 1},
        {"block_index": 2, "zone_type": "solid", "contiguous_id": 1},
    ]

    export_to_boundary_condition("boundary_conditions.bcs", job, bcg, gifs=gifs, volume_zones=volume_zones)
    export_to_job_file(job, "jobfile", title="SingleBlade")
