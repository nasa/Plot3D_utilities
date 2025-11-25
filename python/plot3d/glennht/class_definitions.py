# glennht_export_classes.py
from __future__ import annotations
from dataclasses import dataclass, field
from enum import IntEnum, Enum
from typing import Any, Dict, List, Optional

# ----------------------------
# Enums (mirror your C#)
# ----------------------------
class BoundaryConditionType(IntEnum):
    Inlet = 1
    Outlet = 2
    SymmetryOrSlip = 3
    Wall = 4
    GIF = 1000

class InletBC_Subtype(IntEnum):
    Normal = 0
    AngleSpecified = 1
    AngleAndProfileSpecified = 2

class InletBC_Direction(IntEnum):
    Something = 1
    Annular = 2
    Cascade = 3

class OutletSubtype(IntEnum):
    UniformPressure = 0
    AveragePressure = 1
    UnsteadyFlow = 2

class SymmetricSlipSubtype(IntEnum):
    Symmetry = 0
    Slip = 1

class WallSubtype(IntEnum):
    SpecifiedWallHeatFlux = 0
    SpecifiedWallTemperature = 1
    BCWall_Subtype_Conjugate = 3

class GIFCoordinate(IntEnum):
    Cartesian = 0
    Polar = 1

class GIFType(IntEnum):
    Conjugate = 1
    StraightInterpolation = 2

class GIFOrder(IntEnum):
    ZeroOrder = 0
    LinearOrder = 1
    CubicOrder = 3

class TbModelType(IntEnum):
    NO_TURB_MODEL = 0
    K_OMEGA_TURB_MODEL = 1
    K_EPS_TURB_MODEL = 2
    ARSM_TURB_MODEL = 3
    RSM_TURB_MODEL = 4
    SST_TURB_MODEL = 5
    LES_TURB_MODEL = 7
    WL_KOMEGA_TURB_MODEL = 11
    K_OMEGA_GAMMA_TURB_MODEL = 12

# ----------------------------
# Base + concrete BCs
# ----------------------------
@dataclass
class BoundaryCondition:
    BCType: BoundaryConditionType
    SurfaceID: int
    Name: str = ""
    IsPostProcessing: bool = False
    IsCalculateMassFlow: bool = False
    ToggleProcessSurface: bool = False

@dataclass
class InletBC(BoundaryCondition):
    inlet_subType: Optional[InletBC_Subtype] = InletBC_Subtype.Normal
    inlet_ref_Mach_Nr: Optional[float] = 0.2

    T0_const: Optional[float] = None      # K (will be normalized by refT0)
    P0_const: Optional[float] = None      # physical input + unit → normalized
    P0_const_unit: str = "Pa"             # NOT exported (used only for conversion)
    Tu_const: Optional[float] = None
    Ts_const: Optional[float] = None
    ang1_const: Optional[float] = None
    bet1_const: Optional[float] = None

    annular_inlet: bool = False
    deltah: Optional[float] = None
    deltat: Optional[float] = None
    twall_hub: Optional[float] = None
    twall_case: Optional[float] = None

    have_inlet_prof: bool = False
    filen_inlet_profile: Optional[str] = None
    direction: Optional[InletBC_Direction] = None

@dataclass
class OutletBC(BoundaryCondition):
    outlet_subType: Optional[OutletSubtype] = None
    extrapolation_order: Optional[int] = None
    Pback_extrapolate_profile: bool = False
    Pback_const: Optional[float] = None
    Pback_const_unit: str = "Pa"          # NOT exported
    have_Pback_prof: bool = False
    annular_outlet: bool = False
    approx_Mach_out: Optional[float] = None
    filen_Pback_prof: Optional[str] = None

    mult_for_full_ring: Optional[int] = None
    p_over_pt_ratio_limit: Optional[float] = None

@dataclass
class SymmetricSlipBC(BoundaryCondition):
    slip_subType: Optional[SymmetricSlipSubtype] = None
    slip_omega: Optional[float] = None

@dataclass
class WallBC(BoundaryCondition):
    wall_subType: Optional[int] = None
    Twall_const: Optional[float] = None
    have_Twall_prof: bool = False
    filen_Twall_prof: Optional[str] = None

    Qwall_const: Optional[float] = 0.0
    have_Qwall_prof: bool = False
    filen_Qwall_prof: Optional[str] = None

    BEM_coupled_surf: bool = False
    Nr_wall_segments: Optional[int] = None
    segment_Omega: Optional[float] = None
    segment_xMin: Optional[float] = None

@dataclass
class GIF(BoundaryCondition):
    GIFSurface1: int = 0
    GIFSurface2: int = 0
    Name1: str = ""
    Name2: str = ""
    Coordinates: GIFCoordinate = GIFCoordinate.Cartesian
    Type: GIFType = GIFType.Conjugate
    Order: GIFOrder = GIFOrder.LinearOrder

    def __post_init__(self):
        self.BCType = BoundaryConditionType.GIF

@dataclass
class BCGroup:
    Inlets: List[InletBC] = field(default_factory=list)
    Outlets: List[OutletBC] = field(default_factory=list)
    SymmetricSlips: List[SymmetricSlipBC] = field(default_factory=list)
    Walls: List[WallBC] = field(default_factory=list)

# ----------------------------
# Job structures (mirror C#)
# ----------------------------
@dataclass
class JobFiles:
    DcmpFILE: Optional[str] = "ddcmp.dat"
    ConnFILE: Optional[str] = "connectivity.ght_conn"
    BCSpecFILE: Optional[str] = "boundary_conditions.bcs"
    GridFile: Optional[str] = "mesh.xyz"
    GridFileFormat: Optional[str] = "formatted"
    Plot3DFileFormat: Optional[str] = "formatted"
    SolnInFile: Optional[str] = "In.soln"
    SolnInFileFormat: Optional[str] = "unformatted"
    SolnOutFile: Optional[str] = "Out.soln"
    SolnOutFileFormat: Optional[str] = "unformatted"
    residFILE: Optional[str] = "his.subs"
    residFILE2: Optional[str] = "his.nosubs"

@dataclass
class JobControl:
    mRunLevel: Optional[int] = 0
    LUNout: Optional[int] = 6
    RestartSoln: bool = False
    SaveSoln: bool = True
    SaveTransientSoln: bool = False
    VerboseScreenOutput: bool = True

@dataclass
class TurbModelInput:
    TbModelType: TbModelType = TbModelType.K_OMEGA_TURB_MODEL
    PRNS_ResolutionParameter: float = 1.0

@dataclass
class Plot3DParameters:
    Plot3DParameterSet: str = "Standard"

@dataclass
class InitialCond:
    P0: Optional[float] = 1.0
    T0: Optional[float] = 1.0
    Minit: Optional[float] = 0.0
    alfa: Optional[float] = 0.0
    beta: Optional[float] = 0.0
    Tu: Optional[float] = 0.05
    Ts: Optional[float] = 0.05
    T0_solid: Optional[float] = 1.0
    annular_init: Optional[float] = 0.0

@dataclass
class TimeStpControl:
    UnsteadyFlow: bool = False
    FullyImplicitDiscr: bool = False
    EulerBackward: bool = False
    CranckNicholson: bool = False
    BlendedTime: bool = True
    TransientPlot3dFiles: bool = False
    have_previous_step: bool = False
    UseLowMPrecond: bool = False
    pcMinRefMach: Optional[float] = 1.0e-3
    dissRefMach: Optional[float] = 1.0e-3
    Implicitness: Optional[float] = 1.0
    dt_unst: Optional[float] = 1.0e-8
    CFLn: Optional[float] = 0.25
    CFLr: Optional[float] = 0.125
    cst: Optional[float] = 3.5
    cst_solid: Optional[float] = 7.0
    convergTolerance: Optional[float] = 5e-5
    nTimeSteps: Optional[int] = 50
    maxPseudoSteps: Optional[int] = 50
    maxPseudoSteps_solid: Optional[int] = 100
    fully_coupled_solid: bool = False
    ReinitializeTime: bool = False
    ResetTimeTo: Optional[float] = 0.0
    nTransBegin: Optional[int] = 0
    nfiles: Optional[int] = 10
    ninterval: Optional[int] = 10
    nHL: Optional[int] = 1
    time_avg_start: Optional[float] = None
    time_avg_end: Optional[float] = None
    restart_timeAvg: bool = False

@dataclass
class SPDSchemeControl:
    # GlennHT-style defaults (with shorthand preserved)
    NS_Central: str = "4*T"
    TB2_Upwind1: str = "4*T"
    NS_Upwind2: str = "4*F"

    ScalrCoeff_ArtDiss: bool = True
    useSecDiffArtDiss: bool = True
    useFrthDiffArtDiss: bool = True

    rk2: str = "4*0.12500"
    rk4: str = "4*0.032"

    # Keep other optional fields from your previous version
    NS_Upwind1: bool = False
    use_AUSM_Chima: bool = False
    use_AUSM_Liou_hTot: bool = False
    TB2_Central: bool = False
    TBRSM_Central: bool = False
    TBRSM_Upwind1: bool = True
    constArtDiss: bool = False
    scalarArtDiss: bool = True
    MatrxCoeff_ArtDiss: bool = False
    secDiffArtDiss: bool = True
    matrixArtDiss: bool = False
    frthDiffArtDiss: bool = True
    MachCutOff: float = 0.1
    ivanAlbada: int = 1

@dataclass
class RKSchemeControl:
    nStages: int = 4
    RKCoeff: str = "0.25,0.3333333,0.5,1.,6*0"
    compute_pdiff_in_stage: str = "T,T,T,T,6*F"
    compute_adiss_in_stage: str = "T,T,T,T,6*F"
    export_import_after_stage: str = "T,T,T,T,6*F"
    use_implicit_residual_smoothing: str = ".T."
    irs_neqs: Optional[int] = 1
    irs_use_GS: bool = True
    n_GS_iterations: Optional[int] = 3
    n_GlobalSweeps: Optional[int] = 1

@dataclass
class MGSchemeControl:
    FinestLevel: Optional[int] = 0
    CoarsestLevel: Optional[int] = 0
    pre_mg_sweeps: Optional[int] = 1
    mg_sweeps: Optional[int] = 0
    post_mg_sweeps: Optional[int] = 0
    SVFinestLevel: Optional[int] = 0
    SVCoarsestLevel: Optional[int] = 0
    SVpre_mg_sweeps: Optional[int] = 1
    SVmg_sweeps: Optional[int] = 0
    SVpost_mg_sweeps: Optional[int] = 0

@dataclass
class GasPropertiesInput:
    UseDryAir: bool = True
    Use_const_Cp: bool = True
    Use_const_trProp: bool = False
    Use_RefT: bool = True
    RefT_Properties: Optional[float] = -1.0E+99
    const_cp: Optional[float] = -1.0E+99
    const_visc: Optional[float] = -1.0E+99
    const_kth: Optional[float] = -1.0E+99
    Use_specialGas: bool = False
    SpecialGasMW: Optional[float] = -1.0E+99

@dataclass
class ReferenceCond:
    # Derived (not user input) — kept for export
    useDimensionalVariables: bool = False
    refLen: Optional[float] = 1.0
    refP0: Optional[float] = 101325.0  # Pa (total)
    refT0: Optional[float] = 300.0
    refRho0: Optional[float] = 1.1765823
    refVel: Optional[float] = 293.4588
    refVisc: Optional[float] = 1.84e-5
    refCond: Optional[float] = 0.02636
    refCp: Optional[float] = 1004.5784
    MolW: Optional[float] = 28.964
    RgasUnv: Optional[float] = 8314.4126
    Rgas: Optional[float] = 287.06023
    gamma: Optional[float] = 1.4
    Re: Optional[float] = 1.8765e7
    Pr: Optional[float] = 0.706
    ndVisc: Optional[float] = 1.0
    ndCond: Optional[float] = 1.0
    Omegab: Optional[float] = 0.0
    ReScalingFactor: Optional[float] = 1.0
    rho_solid: Optional[float] = None
    cond_solid: Optional[float] = None
    Csp_solid: Optional[float] = None

@dataclass
class ReferenceCondFull:
    # User-/physics-derived, authoritative
    reflen: Optional[float] = None
    refP0: Optional[float] = None
    refT0: Optional[float] = None
    refrho0: Optional[float] = None
    refVel: Optional[float] = None
    refvisc: Optional[float] = None
    refcond: Optional[float] = None
    refCp: Optional[float] = None
    MolW: Optional[float] = None
    RgasUnv: Optional[float] = None
    Rgas: Optional[float] = None
    gamma: Optional[float] = None
    Re: Optional[float] = None
    Pr: Optional[float] = None
    ndVisc: Optional[float] = None
    ndCond: Optional[float] = None
    Omegab: Optional[float] = None
    ReScalingFactor: Optional[float] = None
    rho_solid: Optional[float] = None
    cond_solid: Optional[float] = None
    csp_solid: Optional[float] = None

@dataclass
class Job:
    JobFiles: JobFiles = field(default_factory=JobFiles)
    JobControl: JobControl = field(default_factory=JobControl)
    TurbModelInput: TurbModelInput = field(default_factory=TurbModelInput)
    Plot3DParameters: Plot3DParameters = field(default_factory=Plot3DParameters)
    InitialCond: InitialCond = field(default_factory=InitialCond)
    TimeStpControl: TimeStpControl = field(default_factory=TimeStpControl)
    SPDSchemeControl: SPDSchemeControl = field(default_factory=SPDSchemeControl)
    RKSchemeControl: RKSchemeControl = field(default_factory=RKSchemeControl)
    MGSchemeControl: MGSchemeControl = field(default_factory=MGSchemeControl)
    GasPropertiesInput: GasPropertiesInput = field(default_factory=GasPropertiesInput)

    # Only carry the FULL reference; compact RefCond will be derived by the exporter
    ReferenceCondFull: ReferenceCondFull = field(default_factory=ReferenceCondFull)

    # Filled/derived later by exporter (kept optional to avoid being “input”)
    ReferenceCond: Optional[ReferenceCond] = None

# VolumeZone placeholder (dict-like)
VolumeZone = Dict[str, Any]
