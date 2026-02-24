# glennht_export_classes.py
"""Dataclass and enum definitions for GlennHT boundary condition and job configuration.

This module mirrors the C# class hierarchy used by GlennHT and provides
Python dataclasses for boundary conditions, job control parameters, scheme
settings, gas properties, and reference conditions.  These definitions are
consumed by the export functions in ``export_functions.py``.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from enum import IntEnum, Enum
from typing import Any, Dict, List, Optional

# ----------------------------
# Enums (mirror your C#)
# ----------------------------
class BoundaryConditionType(IntEnum):
    """Top-level boundary condition category."""
    Inlet = 1
    Outlet = 2
    SymmetryOrSlip = 3
    Wall = 4
    GIF = 1000

class InletBC_Subtype(IntEnum):
    """Inlet boundary condition sub-classification."""
    Normal = 1
    AngleSpecified = 2
    AngleAndProfileSpecified = 3

class InletBC_Direction(IntEnum):
    """Flow direction convention for inlet boundary conditions."""
    Something = 1
    Annular = 2
    Cascade = 3

class OutletSubtype(IntEnum):
    """Outlet boundary condition sub-classification."""
    UniformPressure = 0
    AveragePressure = 1
    UnsteadyFlow = 2

class SymmetricSlipSubtype(IntEnum):
    """Sub-type for symmetry or slip-wall boundary conditions."""
    Symmetry = 0
    Slip = 1

class WallSubtype(IntEnum):
    """Wall thermal boundary condition sub-type."""
    SpecifiedWallHeatFlux = 0
    SpecifiedWallTemperature = 1
    BCWall_Subtype_Conjugate = 3

class GIFCoordinate(IntEnum):
    """Coordinate system used for General Interface (GIF) interpolation."""
    Cartesian = 0
    Polar = 1

class GIFType(IntEnum):
    """Interpolation strategy for a General Interface patch."""
    Conjugate = 1
    StraightInterpolation = 2

class GIFOrder(IntEnum):
    """Polynomial order used for GIF interpolation."""
    ZeroOrder = 0
    LinearOrder = 1
    CubicOrder = 3

class TbModelType(IntEnum):
    """Turbulence model selector."""
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
    """Base class for all GlennHT boundary conditions.

    Attributes:
        BCType: High-level boundary condition category.
        SurfaceID: Integer surface identifier used in the connectivity file.
        Name: Human-readable label written as a comment in the ``.bcs`` file.
        IsPostProcessing: When True, marks the surface as a post-processing
            reference surface (``BRefCond=T`` in the namelist).
        IsCalculateMassFlow: When True, requests mass-flow calculation on this
            surface (``BCalc=T``).
        ToggleProcessSurface: Alternative flag that also sets ``BCalc=T``.
    """
    BCType: BoundaryConditionType
    SurfaceID: int
    Name: str = ""
    IsPostProcessing: bool = False
    IsCalculateMassFlow: bool = False
    ToggleProcessSurface: bool = False

@dataclass
class InletBC(BoundaryCondition):
    """Inlet boundary condition with optional profile and annular geometry.

    Physical pressure and temperature values are supplied in user-specified
    units; the exporter normalizes them by the reference conditions before
    writing the namelist.

    Attributes:
        inlet_subType: Selects how the flow angle is specified.
        surfID_inlet: Surface ID echoed inside the ``INLET_BC`` namelist block.
        inlet_ref_Mach_Nr: Reference Mach number at the inlet.
        T0_const: Total temperature in Kelvin; normalized by ``refT0`` on export.
        P0_const: Total pressure in the units given by ``P0_const_unit``;
            normalized by ``refP0`` on export.
        P0_const_unit: Unit string for ``P0_const`` (e.g. ``"Pa"``, ``"bar"``).
            Not written to the namelist.
        Tu_const: Turbulence intensity.
        Ts_const: Turbulent length scale; normalized by ``reflen`` on export.
        ang1_const: Flow angle (primary).
        bet1_const: Flow angle (secondary).
        annular_inlet: Enable annular inlet geometry.
        deltah: Hub delta height for annular geometry.
        deltat: Tip delta height for annular geometry.
        twall_hub: Hub wall temperature; normalized by ``refT0`` on export.
        twall_case: Casing wall temperature; normalized by ``refT0`` on export.
        have_inlet_prof: Whether a profile file is provided.
        filen_inlet_profile: Path to the inlet profile file.
        direction: Optional flow direction convention.
    """
    inlet_subType: InletBC_Subtype = InletBC_Subtype.Normal
    surfID_inlet:int = 1
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
    """Outlet boundary condition with back-pressure specification.

    Attributes:
        outlet_subType: Selects how the back-pressure is applied.
        surfID_outlet: Surface ID echoed inside the ``OUTLET_BC`` namelist.
        extrapolation_order: Order of extrapolation at the outlet patch.
        Pback_extrapolate_profile: Use a radial profile for back-pressure
            extrapolation when True.
        Pback_const: Back-pressure in the units given by ``Pback_const_unit``;
            normalized by ``refP0`` on export.
        Pback_const_unit: Unit string for ``Pback_const``.  Not written to the
            namelist.
        have_Pback_prof: Whether a back-pressure profile file is provided.
        annular_outlet: Enable annular outlet geometry.
        approx_Mach_out: Approximate exit Mach number (used for averaging).
        filen_Pback_prof: Path to the back-pressure profile file.
        mult_for_full_ring: Multiplier mapping a sector to the full annulus.
        p_over_pt_ratio_limit: Lower bound on static-to-total pressure ratio.
    """
    outlet_subType: OutletSubtype = OutletSubtype.UniformPressure
    surfID_outlet:int = 2
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
    """Symmetry or slip-wall boundary condition.

    Attributes:
        surfID_symmetricSlip: Surface ID echoed inside the ``SLIP_BC`` namelist.
        slip_subType: Distinguishes a true symmetry plane from a slip wall.
        slip_omega: Rotational speed (rad/s) for a moving slip wall.
    """
    surfID_symmetricSlip:int = 3
    slip_subType: SymmetricSlipSubtype = SymmetricSlipSubtype.Slip
    slip_omega: Optional[float] = None

@dataclass
class WallBC(BoundaryCondition):
    """Viscous wall boundary condition with optional heat-flux or temperature profile.

    Attributes:
        wall_subType: Thermal treatment applied at the wall.
        surfID_wall: Surface ID echoed inside the ``WALL_BC`` namelist.
        Twall_const: Constant wall temperature (K).
        have_Twall_prof: Whether a wall-temperature profile file is provided.
        filen_Twall_prof: Path to the wall-temperature profile file.
        Qwall_const: Constant wall heat flux (W/m²); defaults to adiabatic (0).
        have_Qwall_prof: Whether a wall heat-flux profile file is provided.
        filen_Qwall_prof: Path to the wall heat-flux profile file.
        BEM_coupled_surf: Enable BEM conjugate coupling on this surface.
        Nr_wall_segments: Number of wall segments for BEM coupling.
        segment_Omega: Rotational speed (rad/s) of the wall segment.
        segment_xMin: Axial start position of the rotating segment.
    """
    wall_subType: WallSubtype = WallSubtype.SpecifiedWallHeatFlux
    surfID_wall:int = 3
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
    """General Interface (GIF) patch pair linking two non-conformal surfaces.

    Attributes:
        GIFSurface1: Surface ID of the first GIF patch.
        GIFSurface2: Surface ID of the second GIF patch.
        Name1: Label for the first GIF surface.
        Name2: Label for the second GIF surface.
        Coordinates: Coordinate system used during interpolation.
        Type: Interpolation type applied across the interface.
        Order: Polynomial order of the interpolation.
    """
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
    """Container grouping all boundary condition lists for a single case.

    Attributes:
        Inlets: All inlet boundary conditions.
        Outlets: All outlet boundary conditions.
        SymmetricSlips: All symmetry and slip-wall boundary conditions.
        Walls: All viscous wall boundary conditions.
    """
    Inlets: List[InletBC] = field(default_factory=list)
    Outlets: List[OutletBC] = field(default_factory=list)
    SymmetricSlips: List[SymmetricSlipBC] = field(default_factory=list)
    Walls: List[WallBC] = field(default_factory=list)

# ----------------------------
# Job structures (mirror C#)
# ----------------------------
@dataclass
class JobFiles:
    """File path settings for all GlennHT input and output files.

    Attributes:
        DcmpFILE: Domain decomposition file.
        ConnFILE: Block connectivity file (``.ght_conn``).
        BCSpecFILE: Boundary condition specification file (``.bcs``).
        GridFile: Plot3D mesh file.
        GridFileFormat: Format of the grid file (``"formatted"`` or
            ``"unformatted"``).
        Plot3DFileFormat: Format of Plot3D solution files.
        SolnInFile: Restart solution file read at startup.
        SolnInFileFormat: Format of the restart file.
        SolnOutFile: Solution file written at the end of the run.
        SolnOutFileFormat: Format of the output solution file.
        residFILE: Residual history file with sub-iterations.
        residFILE2: Residual history file without sub-iterations.
    """
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
    """High-level solver execution controls.

    Attributes:
        mRunLevel: Solver run level (0 = normal solve).
        LUNout: Fortran logical unit number for screen output.
        RestartSoln: Read a restart solution when True.
        SaveSoln: Write the final solution when True.
        SaveTransientSoln: Save transient solution snapshots when True.
        VerboseScreenOutput: Enable verbose logging to screen when True.
    """
    mRunLevel: Optional[int] = 0
    LUNout: Optional[int] = 6
    RestartSoln: bool = False
    SaveSoln: bool = True
    SaveTransientSoln: bool = False
    VerboseScreenOutput: bool = True

@dataclass
class TurbModelInput:
    """Turbulence model selection and resolution parameters.

    Attributes:
        TbModelType: Turbulence model to activate.
        PRNS_ResolutionParameter: Near-wall resolution parameter used by some
            models.
    """
    TbModelType: TbModelType = TbModelType.K_OMEGA_TURB_MODEL
    PRNS_ResolutionParameter: float = 1.0

@dataclass
class Plot3DParameters:
    """Plot3D output parameter set selection.

    Attributes:
        Plot3DParameterSet: Named set controlling which flow variables are
            written to the Plot3D solution file.
    """
    Plot3DParameterSet: str = "Standard"

@dataclass
class InitialCond:
    """Non-dimensional initial flow conditions for the entire domain.

    All values are non-dimensional (normalized by the reference quantities
    defined in ``ReferenceCond``).

    Attributes:
        P0: Initial total pressure ratio.
        T0: Initial total temperature ratio.
        Minit: Initial Mach number.
        alfa: Initial flow angle (primary).
        beta: Initial flow angle (secondary).
        Tu: Initial turbulence intensity.
        Ts: Initial turbulent length scale.
        T0_solid: Initial solid-region temperature ratio.
        annular_init: Annular geometry initialization parameter.
    """
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
    """Time-stepping and convergence controls.

    Attributes:
        UnsteadyFlow: Enable time-accurate (unsteady) mode.
        FullyImplicitDiscr: Use fully implicit temporal discretization.
        EulerBackward: Use first-order backward Euler time integration.
        CranckNicholson: Use Crank-Nicolson (second-order) time integration.
        BlendedTime: Blend explicit and implicit time contributions.
        TransientPlot3dFiles: Write Plot3D files at each physical time step.
        have_previous_step: Indicates a previous time-step solution exists.
        UseLowMPrecond: Enable low-Mach preconditioning.
        pcMinRefMach: Minimum reference Mach for preconditioning.
        dissRefMach: Reference Mach number for artificial dissipation scaling.
        Implicitness: Blending factor between explicit (0) and implicit (1).
        dt_unst: Physical time step size (seconds) for unsteady runs.
        CFLn: CFL number for normal time advancement.
        CFLr: CFL number for residual smoothing.
        cst: Courant-number-like constant for the fluid solver.
        cst_solid: Courant-number-like constant for the solid solver.
        convergTolerance: Convergence criterion on the L2 residual.
        nTimeSteps: Total number of physical time steps.
        maxPseudoSteps: Maximum pseudo (inner) iterations per time step for
            fluid.
        maxPseudoSteps_solid: Maximum pseudo iterations per time step for
            solid.
        fully_coupled_solid: Solve fluid and solid simultaneously when True.
        ReinitializeTime: Reset the physical time counter when True.
        ResetTimeTo: Physical time value to reset to when ``ReinitializeTime``
            is True.
        nTransBegin: Physical time step index at which to start writing
            transient files.
        nfiles: Number of transient Plot3D snapshots to write.
        ninterval: Interval (in time steps) between transient snapshots.
        nHL: Number of harmonic levels for harmonic balance runs.
        time_avg_start: Start time for time-averaging.
        time_avg_end: End time for time-averaging.
        restart_timeAvg: Restart time-averaging from zero when True.
    """
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
    """Spatial discretization (SPD) scheme and artificial dissipation controls.

    Attributes:
        NS_Central: Central-difference stencil flags for Navier-Stokes (e.g.
            ``"4*T"``).
        TB2_Upwind1: First-order upwind flags for two-equation turbulence.
        NS_Upwind2: Second-order upwind flags for Navier-Stokes.
        ScalrCoeff_ArtDiss: Use scalar coefficients for artificial dissipation.
        useSecDiffArtDiss: Include second-difference artificial dissipation.
        useFrthDiffArtDiss: Include fourth-difference artificial dissipation.
        rk2: Second-stage Runge-Kutta artificial dissipation coefficients.
        rk4: Fourth-stage Runge-Kutta artificial dissipation coefficients.
        NS_Upwind1: Enable first-order upwind for Navier-Stokes.
        use_AUSM_Chima: Use Chima's AUSM flux scheme.
        use_AUSM_Liou_hTot: Use Liou's AUSM+ scheme with total enthalpy.
        TB2_Central: Central-difference for two-equation turbulence transport.
        TBRSM_Central: Central-difference for RSM turbulence transport.
        TBRSM_Upwind1: First-order upwind for RSM turbulence transport.
        constArtDiss: Use constant artificial dissipation coefficients.
        scalarArtDiss: Use scalar (non-matrix) artificial dissipation.
        MatrxCoeff_ArtDiss: Use matrix coefficients for artificial dissipation.
        secDiffArtDiss: Apply second-difference dissipation term.
        matrixArtDiss: Apply matrix-valued artificial dissipation.
        frthDiffArtDiss: Apply fourth-difference dissipation term.
        MachCutOff: Mach number below which dissipation is clipped.
        ivanAlbada: Van Albada limiter flag (1 = enabled).
    """
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
    """Runge-Kutta time-marching scheme and implicit residual smoothing settings.

    Attributes:
        nStages: Number of Runge-Kutta stages.
        RKCoeff: Comma-separated stage coefficients.
        compute_pdiff_in_stage: Per-stage flags for pressure-difference update.
        compute_adiss_in_stage: Per-stage flags for artificial dissipation
            update.
        export_import_after_stage: Per-stage flags for halo exchange.
        use_implicit_residual_smoothing: Enable implicit residual smoothing
            (``".T."`` or ``".F."``).
        irs_neqs: Number of equations smoothed implicitly.
        irs_use_GS: Use Gauss-Seidel sweeps for implicit smoothing.
        n_GS_iterations: Number of Gauss-Seidel iterations per stage.
        n_GlobalSweeps: Number of global smoothing sweeps.
    """
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
    """Multigrid scheme settings for fluid and solid solvers.

    Attributes:
        FinestLevel: Index of the finest grid level (0 = original mesh).
        CoarsestLevel: Index of the coarsest grid level.
        pre_mg_sweeps: Pre-smoothing sweeps on each level.
        mg_sweeps: Number of multigrid V-cycle sweeps.
        post_mg_sweeps: Post-smoothing sweeps on each level.
        SVFinestLevel: Finest level for the scalar variable (turbulence) solver.
        SVCoarsestLevel: Coarsest level for the scalar variable solver.
        SVpre_mg_sweeps: Pre-smoothing sweeps for the scalar variable solver.
        SVmg_sweeps: Multigrid sweeps for the scalar variable solver.
        SVpost_mg_sweeps: Post-smoothing sweeps for the scalar variable solver.
    """
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
    """Gas model and thermophysical property selections.

    Attributes:
        UseDryAir: Use the built-in dry-air property model when True.
        Use_const_Cp: Treat specific heat as constant when True.
        Use_const_trProp: Use constant transport properties (mu, k) when True.
        Use_RefT: Evaluate properties at a reference temperature when True.
        RefT_Properties: Reference temperature (K) for property evaluation.
        const_cp: Constant specific heat at constant pressure (J/kg/K).
        const_visc: Constant dynamic viscosity (Pa·s).
        const_kth: Constant thermal conductivity (W/m/K).
        Use_specialGas: Enable a user-defined special gas model when True.
        SpecialGasMW: Molecular weight (kg/kmol) of the special gas.
    """
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
    """Compact non-dimensional reference condition block written to the job file.

    This dataclass is derived by the exporter from ``ReferenceCondFull`` and
    is not intended to be populated directly by the user.

    Attributes:
        useDimensionalVariables: Solve in dimensional units when True.
        refLen: Reference length (m).
        refP0: Reference total pressure (Pa).
        refT0: Reference total temperature (K).
        refRho0: Reference density (kg/m³).
        refVel: Reference velocity (m/s).
        refVisc: Reference dynamic viscosity (Pa·s).
        refCond: Reference thermal conductivity (W/m/K).
        refCp: Reference specific heat at constant pressure (J/kg/K).
        MolW: Molecular weight (kg/kmol).
        RgasUnv: Universal gas constant (J/kmol/K).
        Rgas: Specific gas constant (J/kg/K).
        gamma: Ratio of specific heats.
        Re: Reference Reynolds number.
        Pr: Prandtl number.
        ndVisc: Non-dimensional viscosity (typically 1.0).
        ndCond: Non-dimensional thermal conductivity (typically 1.0).
        Omegab: Blade rotational speed (rad/s).
        ReScalingFactor: Factor applied to scale the Reynolds number.
        rho_solid: Solid material density (kg/m³).
        cond_solid: Solid thermal conductivity (W/m/K).
        Csp_solid: Solid specific heat (J/kg/K).
    """
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
    """Authoritative, physics-derived reference condition set.

    Populated by ``populate_reference_from_inputs`` and used as the source of
    truth when constructing the compact ``ReferenceCond`` block.

    Attributes:
        reflen: Reference length (m).
        refP0: Reference total pressure (Pa).
        refT0: Reference total temperature (K).
        refrho0: Reference density (kg/m³).
        refVel: Reference velocity (m/s).
        refvisc: Reference dynamic viscosity (Pa·s).
        refcond: Reference thermal conductivity (W/m/K).
        refCp: Reference specific heat at constant pressure (J/kg/K).
        MolW: Molecular weight (kg/kmol).
        RgasUnv: Universal gas constant (J/kmol/K).
        Rgas: Specific gas constant (J/kg/K).
        gamma: Ratio of specific heats.
        Re: Reference Reynolds number.
        Pr: Prandtl number.
        ndVisc: Non-dimensional viscosity (typically 1.0).
        ndCond: Non-dimensional thermal conductivity (typically 1.0).
        Omegab: Blade rotational speed (rad/s).
        ReScalingFactor: Factor applied to scale the Reynolds number.
        rho_solid: Solid material density (kg/m³).
        cond_solid: Solid thermal conductivity (W/m/K).
        csp_solid: Solid specific heat (J/kg/K).
    """
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
    """Top-level container for all GlennHT job configuration sections.

    Each field corresponds to a named Fortran namelist block that is written
    by :func:`~plot3d.glennht.export_functions.export_to_job_file`.

    Attributes:
        JobFiles: Input/output file paths.
        JobControl: Solver execution controls.
        TurbModelInput: Turbulence model selection.
        Plot3DParameters: Plot3D output parameter set.
        InitialCond: Non-dimensional initial conditions.
        TimeStpControl: Time-stepping and convergence parameters.
        SPDSchemeControl: Spatial discretization and dissipation settings.
        RKSchemeControl: Runge-Kutta scheme parameters.
        MGSchemeControl: Multigrid scheme parameters.
        GasPropertiesInput: Gas model and property selections.
        ReferenceCondFull: Authoritative physics-derived reference quantities;
            populated by :func:`~plot3d.glennht.export_functions.populate_reference_from_inputs`.
        ReferenceCond: Compact reference condition block derived from
            ``ReferenceCondFull`` by the exporter; ``None`` until populated.
    """
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

    # Filled/derived later by exporter (kept optional to avoid being "input")
    ReferenceCond: Optional[ReferenceCond] = None

# VolumeZone placeholder (dict-like)
VolumeZone = Dict[str, Any]
