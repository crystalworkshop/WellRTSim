# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

WellRTSim is a MATLAB-based wellbore reactive transport simulator that models coupled flow, heat transfer, and aqueous chemistry in geothermal wells. It integrates with PHREEQC (via the PhreeqcMatlab wrapper) for geochemical equilibrium and kinetic mineral reactions.

**Key technologies:**
- MATLAB (primary language)
- PHREEQC via PhreeqcMatlab external dependency for chemistry
- HDF5 for results storage
- IAPWS-IF97 for water/steam properties

## Setting Up and Running

### Prerequisites

1. Install MATLAB with required toolboxes (core functionality, HDF5 support)
2. Clone or obtain PhreeqcMatlab from https://github.com/simulkade/PhreeqcMatlab
3. PHREEQC database (phreeqc.dat) must be installed on your system

### Running a Simulation

1. In MATLAB, run `PhreeqcMatlab/startup.m` to initialize the PHREEQC interface and add it to the path
2. Edit `main.m` to set `SimDir` to the desired simulation case directory (e.g., `'Simulation/Krafla/KJ-9/'`)
3. Verify the chemistry setup in `SimDir/chemistry.md` (database path, aqueous species, kinetic minerals)
4. Run `main.m` - the simulation will:
   - Initialize steady-state conditions (pressure, enthalpy, velocity profile)
   - Set up graphics UI with well schematic, live transient plots, and wellhead diagnostics
   - Execute time-stepping loop with optional chemical transport
   - Write results to HDF5 files in `SimDir`

### Viewing Results

Add the `viewer/` folder to the path and call `result_view` (optionally with a case directory or `.h5` file, e.g. `result_view('Simulation/Krafla/KJ-9/')`) to launch the interactive results viewer for completed runs. `result_view.m` is only the entry point; the UI lives in the `rv_*` helpers (`rv_launchResultView`, `rv_createResultsViewer`, `rv_loadResultsH5`, `rv_readProfileSnapshot`, `rv_setupResultsViewerAxes`, `rv_updateResultsViewerPlots`). Additional documentation is in `user_manual.pdf`.

### Tests, Build, and Linting

There is no automated test, build, or lint harness — this is a research simulator run interactively from MATLAB. "Running" a change means executing `main.m` (or `result_view`) and inspecting console output, the live plots, and the HDF5 results. Do not look for a test runner or CI.

Simulation outputs are not version-controlled: `.gitignore` excludes `Simulation/**/*.h5`, `Simulation/**/*.mat`, `Simulation/**/*.phr`, the external `PhreeqcMatlab/` checkout, and `AGENTS.md`. The committed style guidance in `AGENTS.md` is mirrored in the "Code Style and Philosophy" section below.

## Code Style and Philosophy

From AGENTS.md: WellRTSim does not use defensive MATLAB code by default.

- Assume valid project inputs; fail loudly on violations
- Avoid excessive guard clauses and silent fallbacks
- Optimize for readability and scientific workflow over generic library robustness
- Add checks only for known invariants that break in this codebase

## Architecture

### State Structure

The central data structure is `state` - a MATLAB struct passed through all major functions. It contains:

- **Grid/geometry:** `n` (cell count), `x` (depths), `dx`, `Lp` (total depth), `Dp` (diameters)
- **Thermodynamic state:** `Y` (3×n matrix of [pressure; enthalpy; velocity]), `T` (temperature profile)
- **Time:** `tt` (current time), `dt` (timestep), `tfin` (final time)
- **Chemistry:** `chem` (struct with PHREEQC handle, species data, composition matrices), `C` (component mass fractions)
- **Graphics:** Figure handles, axes, UI controls (run/pause buttons, status label)
- **Feedzone/well:** Well casing geometry, deviation data, feedzone locations and PI values
- **Results:** HDF5 file handle, output buffers, save intervals

### Major Modules

**Initialization (`functions/init/`):**
- `initializeState.m` - Parse params.md and chemistry.md, set up grids, load well geometry
- `initializeSteadyState.m` - Compute steady-state P-H-U profile via top-down or bottom-up shooting
- `initializeGraphics.m` - Create figure with tabbed layout (well schematic, transient plots, wellhead diagnostics)
- `initResultsH5.m` - Open HDF5 file and create group structure

**Hydrodynamics (`functions/Hydrodynamics/`):**
- `OneStep.m` - Advance one timestep: calls `OneIter` in a nonlinear iteration loop, adjusts `dt` based on iteration count
- `OneIter.m` - Single Newton iteration: computes RHS, assembles Jacobian, solves linear system
- `RHS_v2.m` - Evaluates 3×n nonlinear residual (momentum, energy, mass conservation) at node i
- `TDMA.m` - Tridiagonal matrix algorithm solver

The state vector `Y = [P; H; U]` (pressure, enthalpy, velocity) is solved via implicit Euler with Newton linearization. Boundary conditions (top/bottom) and feedzone mass/energy injections are incorporated into RHS and Jacobian.

**Chemistry (`functions/Chemistry/`):**
- `initChemistryV2.m` - Parse chemistry setup, create PHREEQC handle, extract molar masses and formula weights
- `chemistryStepV2.m` - Main chemistry timestep: transport and PHREEQC equilibrium/kinetics
- `runPhreeqcStepV2.m` - Render PHREEQC input script from template, run per-cell equilibrium and kinetics
- `parseChemistrySetup.m` - Parse chemistry.md configuration file
- `IAPWS_IF97.m` - Industrial-formulation water properties (large static data file)
- `calculatePhaseProperties.m` - Evaluate density, enthalpy, viscosity for liquid and vapor phases

Chemistry is optional (`calc_chem == 1`). When enabled, it runs after hydrodynamic initialization reaches `stat_chem` time. Chemical components are tracked as mass fractions and transported via advection-dispersion. PHREEQC is called per cell to compute equilibrium activities and kinetic mineral reaction rates (Calcite, Quartz, Anhydrite are common examples).

**Graphics (`functions/Graphics/`):**
- `createWellPlots.m` - Left panel: well schematic (casing, perforations); right panel: well trajectory
- `plotResultsOnAxes.m` - Large function that plots pressure, enthalpy, velocity, temperature, or chemistry species vs. depth
- `setupTransientAxes.m`, `setupWellheadAxes.m` - Configure live plot axes
- `updateWellheadPlots.m`, `updateWellheadDiagnostics.m` - Update live wellhead time-series and diagnostics display

**Wellspec/Geometry (`functions/Wellspec/`, `functions/IO/`):**
- `read_geom.m` - Load CASTINGDETAIL.csv, DEVIATIONDETAIL.csv; create interpolants for diameter, roughness, offset, angle
- `appendProfileH5.m`, `appendWellheadH5.m`, `appendChemistryH5.m` - Write spatial profiles and time-series to HDF5

### Simulation Case Structure

Each case is a directory (e.g., `Simulation/Krafla/KJ-9/`) containing:

- **params.md** - Key-value parameters (grid size, time settings, boundary conditions, physical properties)
- **chemistry.md** - Chemistry configuration (aqueous species concentrations, gases, kinetic minerals, PHREEQC database path)
- **chemistry.pht** - PHREEQC template file (SOLUTION, EQUILIBRIUM_PHASES, KINETICS blocks with placeholders)
- **CASTINGDETAIL.csv** - Casing ID, OD, top/bottom depth, type (open-hole, slotted-liner, etc.)
- **DEVIATIONDETAIL.csv** - Well inclination and lateral offset vs. measured depth
- **Feedzones.csv** - Lateral feed entry locations and flow rates (optional)
- **InitTemperature.csv** - Initial temperature profile (optional; otherwise computed from geothermal gradient)
- **Pressure.csv**, **Temperature.csv** - Boundary condition timeseries (optional)
- **Results HDF5 files** - Output (generated during simulation)

Parameters are parsed as `key = value` pairs; units are specified in the parameter file (e.g., `P_unit = bar` converts input pressures to Pa).

## Important Implementation Details

### Pressure Units

Internal SI (Pa). User input pressure in params.md is specified as `P_unit` (bar or MPa). The parser scales:
- `state.P_top` and `state.P_bot` are immediately converted to Pa
- `state.pressureUnitScale` stores the conversion factor for display

### Time Integration

- Time stepping is adaptive: `dt` increases if iterations < 8, decreases if > maxiter
- Residual tolerance is `epsQ` (default 1e-5)
- Maximum iterations per step is `maxiter` (default 10)
- Chemistry typically runs semi-implicitly: transport is split from PHREEQC equilibrium/kinetics

### Steady-State Initialization

Two modes (set `IC_switch` in params.md):
1. **Top-down:** Start from wellhead P and H, integrate downward using PI and feedzone data
2. **Bottom-up:** Start from bottom P and T, shoot upward to match target flow rate (default and recommended)

Steady state solves the 1D two-phase flow equations (momentum + energy) in implicit form.

### PHREEQC Integration

PHREEQC is called per-cell (not well-averaged). Each cell's aqueous composition is tracked as mass fractions of primary components (Ca, C(4), Si, Na, K, S(6), S(-2), Cl). PHREEQC templates (chemistry.pht) define:
- Initial solution composition (from initial mass fractions)
- Equilibrium mineral phases (precipitation/dissolution constraints)
- Kinetic minerals with rate laws (e.g., surface area / reactive transport time)

Gas partition coefficients (CO2, H2S) can be computed analytically (KD correlations) or via PHREEQC.

### UI Controls

During simulation run:
- **Run button** - Start or resume
- **Pause button** - Pause and show intermediate state
- **Cancel button** - Stop and return to final state
- **Status label** - Current simulation time and completion %
- **Tabs** - Switch between well schematic, transient plots, and wellhead diagnostics

Live plots update every `pltf` steps (e.g., `pltf = 10` → plot every 10 timesteps).

## Common Tasks

### Modify Simulation Parameters

Edit the `params.md` file in your case directory. Key parameters:
- `n` - Number of grid cells
- `P_bot`, `T_bot` - Initial bottom conditions
- `P_top` - Outlet pressure (if `iBC_top == 3`)
- `Q_init` - Target initial flow rate
- `dt`, `dt_max` - Timestep control
- `tfin` - Final simulation time
- `calc_chem` - Enable chemistry (1) or disable (0)

### Add a New Chemistry Component

1. Add the component name to the `aqueous_names` list in chemistry.md
2. Specify initial concentration (ppm) and suffix (e.g., "as HCO3+")
3. Update the PHREEQC template (chemistry.pht) to include the species in the SOLUTION block
4. Ensure the PHREEQC database (phreeqc.dat) defines the species

The parser (`parseChemistrySetup.m`) will auto-extract molar masses from the database.

### Debug Initialization or Early Timesteps

Look at:
- `main.m` console output: steady-state P, H, Q values and initial plots
- HDF5 file: `appendProfileH5` writes initial state at t=0
- PHREEQC logs: If PHREEQC fails, check database path in chemistry.md and ensure PhreeqcMatlab is initialized

### Extract Results Post-Simulation

Use `viewer/result_view.m` or load HDF5 directly:
```matlab
data = h5read('Simulation/Krafla/KJ-9/Krafla_001.h5', '/Profiles/P');
time = h5read('Simulation/Krafla/KJ-9/Krafla_001.h5', '/Time');
```

HDF5 structure: `/Profiles/` (spatial), `/Wellhead/` (time-series), `/Chemistry/` (optional chemical data).

## Known Limitations and Quirks

- PHREEQC must be started externally (via PhreeqcMatlab/startup.m) before running main.m
- Graphics UI is tied to main.m; batch simulations should comment out graphics initialization
- Mineral kinetics currently use simple shrinking-core or first-order rate laws; extend `runPhreeqcStepV2.m` for other models
- Well inclination affects gravity terms; vertical wells are handled, but check gravity calculation in `initializeGravityCache` for highly deviated wells
- Chemistry is single-phase liquid only (gas phases optional but rarely used)

