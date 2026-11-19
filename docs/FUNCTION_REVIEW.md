# Function review and cleanup decisions

This review compares all 32 MATLAB files in the original `src/` tree at
commit `fb6a0fb0eb9b9e34e454277b64e76acbd88800dd` with the current source.
It also covers the new routing, output and direct SWMM helpers.

A function can be useful without being called by `ofm_urban`. Standalone
preprocessing utilities and optional model components are retained. Removal
requires a duplicate implementation, an obsolete implementation superseded
by a retained function, or an unfinished/unused file with no model call path.
A syntax or dependency check does not establish numerical validity.

## Roof drainage correction

The first solver rewrite removed the original roof/ground outlet branch.
That changed a model feature and was a mistake. The branch is restored in
`diffusive_face.m`, shared by `explicit.m` and `overland.m`.

* The original `DEM == 999` building marker remains supported.
* `SHAPE == 999` also identifies a building, allowing actual roof elevations.
* Adjacent roof cells use surface routing. A roof/ground face uses the
  original outlet law: `Q = 0.6 * 0.5 * sqrt(2 * 9.81) * h^(3/2)` in m3/s,
  where `h` is roof water depth above depression storage.
* Roof water is transferred into the adjacent ground cell. The transfer
  removes and adds the same volume. Total withdrawals are limited by the
  available roof water, including when several faces drain one roof cell.
* A dry roof cannot receive ground water through the simplified outlet.
* This preserves the original adjacent-cell downpipe representation. The
  source does not contain a separate map of individual downpipe locations,
  pipe diameters or roof-to-remote-outlet connections.
* The sequential scheme also retains its original steep-slope discharge
  expression. The adaptive scheme uses Manning routing on ordinary faces.

The legacy `DEM=999` convention is ambiguous with real terrain at exactly
999 metres. Direct calls to `diffusive_face` can provide explicit roof masks.
The OFM driver retains the legacy marker for existing configurations.

Tests check the original outlet equation, reversed face orientation, dry
roofs, roof storage, four outlets, actual roof elevations with a building
mask, rainfall reaching the ground, and volume balance in both schemes.

## Original MATLAB files

| Original file | Decision | Purpose, retained behavior, and check |
| --- | --- | --- |
| `ofm_urban.m` | Keep | Main simulation and optional model dispatch. Short ground, roof, boundary, hydrology, channel and sewer runs are tested. |
| `explicit.m` | Keep as canonical solver | Adaptive diffusive routing with shared transfers, sources, depression storage and roof outlets. Tested for conservation, symmetry and nonnegative depths. |
| `explicit2.m` | Merge into `explicit.m` | This was the original driver's selected adaptive implementation. Its intended routing features are implemented in the canonical solver. Internal calls and tests now use `explicit`. |
| `explicit3.m` | Remove duplicate | Another adaptive solver using `overland2`; no active driver selection. No separate model option is lost. |
| `overland.m` | Keep | Sequential scheme selected by `TSTP2D=1`. Predicted depths, roof outlets, steep-slope law and availability limits are retained. |
| `overland2.m` | Merge into `diffusive_face.m` | Original explicit-scheme face calculation. Land/roof conveyance and roof-to-ground drainage are represented in the shared function. |
| `overland3.m` | Remove unused alternative | An unselected predicted-depth variant. Both selected routing schemes and the roof outlet behavior remain available. |
| `cdepth.m` | Keep | Channel storage-to-depth conversion; undefined sinuosity variable corrected. Below/above-bank conversion tests pass. |
| `cflow.m` | Keep | Signed channel discharge and available-water cap. Dry and reverse-flow checks pass. |
| `chnroutin.m` | Keep | Channel transfers, junction handling and observation locations. Channel storage is tested; complex networks need further benchmark validation. |
| `qoutlet.m` | Keep | Surface outlet discharge limited by available water for fractional and long timesteps. Tested. |
| `qnode.m` | Keep | Free inlet discharge and storage cap. Tested. |
| `qmanhole.m` | Keep | Head-dependent surface-to-sewer flow and storage cap. Tested individually and with the actual SWMM engine. |
| `infiltration.m` | Keep | Green-Ampt infiltration and transition to ponding, with supply and pervious-area limits. Ponded equation, dry/impervious cases and coupled balance are tested. |
| `infiltration_old.m` | Remove unused alternative | Older recovery-zone formulation with a different parameter/state interface. Neither the original selected driver nor the current driver calls it. It is not treated as an identical equation to `infiltration.m`. |
| `runoffsimp.m` | Keep | Fast/slow/baseflow reservoirs and percolation. Timestep units and bounded withdrawals are tested with storage balance. |
| `intercsnow.m` | Keep | Canopy interception, snowfall, snowmelt and evaporation. Temperature conversion, zero-LAI/calm-wind handling, roughness geometry and shared ET demand corrected. Canopy/snow and complete-driver water balances pass. |
| `intercept_OLD.m` | Remove unused predecessor | Older canopy function with a different forcing interface. Its declared function name was `intercept`; the active driver uses `intercsnow`. The selected canopy/snow option is retained. |
| `evapotrans.m` | Keep | Energy-balance ET. Coordinate order, humidity fallback, east-positive solar time, nighttime radiation and stability convergence corrected. Day/night, calm wind and both sensible-heat signs are tested. |
| `coordinate.m` | Keep | UTM-to-latitude/longitude conversion used by evapotranspiration. Central-meridian check added. |
| `common.m` | Keep | Physical constants used by the hydrology options. This is an intentionally retained script. |
| `modswfip.m` | Keep | Validated scalar or multilayer soil inputs, explicit hydraulic-model selection and initial states. Whitespace/comments supported without evaluating expressions. |
| `soilground.m` | Keep | Layered SVAT with all 18 original process helpers. Undefined multilayer variables, mixed hydraulic models, rate units, storage bounds and groundwater/bypass accounting corrected. Retention, signed flow and water-balance tests pass. |
| `laivariant.m` | Restore and keep | Standalone seasonal LAI utility. Absence from the runtime call graph does not make it redundant. Winter/summer and output-shape checks pass. |
| `processing.m` | Restore and keep | Public batch-resampling utility. The main simulation samples native records with `sample_forcing` to save memory. Batch totals are tested. |
| `rainfest.m` | Keep | Gauge selection and gridded rainfall conversion. Gauge tie handling and rainfall units are tested. |
| `drainage.m` | Keep, replace backend | Persistent direct SWMM session. Hydraulic state, identifiers, inflows, flood volumes, synchronization and error cleanup are tested with EPA SWMM 5.2.4. |
| `drainage_old.m` | Remove unused predecessor | Unselected file/hotstart coupling using the old wrapper. Superseded by the retained direct session. |
| `swmm_old.m` | Remove unused wrapper | Old `swmm` class in a differently named file. The current simulation does not use this class; the retained native interface supplies the selected sewer coupling. |
| `visualize.m` | Keep | Optional inundation plotting. Requires Mapping Toolbox. Plot rendering is not covered by the numerical regression tests. |
| `parcompute.m` | Remove unused experiment | Demonstrates `parfor` and GPU arrays; no hydrological model call path or completed model feature. |
| `meshgen_notfinish.m` | Remove unfinished experiment | References undefined mesh input variables and is not called by the model. It is not a working mesh-generation option. |

Removed alternatives remain recoverable from the original Git commit above.
No additional original model functions are removed by this review.

## Generated files and superseded native bindings

| Original file | Decision and reason |
| --- | --- |
| `PCUW1D.mdf` | Remove generated synchronization input. The original network inputs remain under `Examples/`; direct coupling prepares one startup input per run. |
| `swmm5.dll` | Replace the old Windows-only library with the locally compiled, pinned EPA engine in `swmm_live`. |
| `read1D.mexw64` | Remove the old binary-output metadata reader; exchange results now come directly from engine memory. |
| `swmm5.h` | Replace the old interface declarations with EPA's pinned header under `third_party/swmm/include/`. The engine header has not been removed from the project. |

## New source files

| File | Purpose and check |
| --- | --- |
| `diffusive_face.m` | Shared face hydraulics, restored roof outlets, and conductance. Direct equation and conservation checks pass. |
| `sample_forcing.m` | Sample one native record interval. Totals, interpolation and gridded input checks pass. Its local `take` helper checks record bounds. |
| `read_ofm_series.m` | Read selected cell series from step snapshots. Times and values are tested. |
| `prepare_swmm_input.m` | Parse topology by identifiers and prepare the startup input. Supplied topology and reordered identifier tests pass. Its local `records` helper selects section data. |
| `build_swmm.m` | Compile the pinned engine and adapter. Build tested with MATLAB R2024b on Apple silicon. |
| `build_ofm.m` | Compile the C++ adaptive solver for local CPU threads. Build tested with MATLAB R2024b on Apple silicon; no OpenMP runtime required. |
| `ofm_route` / `native/ofm_route.cpp` | Optional threaded backend for the canonical adaptive solver. MATLAB/native flux and depth comparisons pass across thread counts, roof outlets, boundaries, masks and one-cell-wide grids. |
| `hydrology_cells.m` | One cell-process implementation shared by serial and worker execution. Model combinations retain the same water accounting. |
| `hydrology_executor.m` | Partition active cells, retain immutable worker parameters and select serial/parallel execution. Saved numerical states match serial results. |
| `build_gpu.m` | Build the Apple silicon Metal backend with a supported local C++ compiler and Apple frameworks. |
| `ofm_metal` / `native/ofm_metal.mm` / `native/ofm_metal_kernels.metal` | Selectable GPU adaptive routing. Compensated storage/transfer arithmetic and a CPU budget check retain conservative updates. CPU/GPU depth, roof, boundary, saved-state and direct-SWMM checks pass at documented GPU tolerances. |
| `swmm_live` / `native/swmm_live.c` | Compiled direct interface. Actual-engine tests pass for clocks, flooding, volume balance, no synchronization-file changes and session recovery. |
| `Contents.m` | MATLAB source index; documentation rather than a model routine. |

The intermediate `write_swmm_interval.m` introduced during the first repair
is also removed: direct coupling no longer writes interval exchange files.

## Soil-water local functions

All 18 local functions in `soilground.m` are retained:

| Functions | Role and review status |
| --- | --- |
| `soilpar`, `soilvar` | Use one hydraulic model consistently, valid layer capacities and per-layer state. CH field capacity and VG analytical storage are tested. |
| `inflow`, `byflfr`, `infpar` | Layer inflow, bypass and infiltration distribution. The parent integrates bounded paired transfers directly; `inflow` remains a reference flux-assembly helper. Bypass and capacity excess return to OFM. |
| `canopy`, `plntres`, `rootgrowth`, `tbylayer` | Handle zero snow density without snow, absent roots, growth start times, nonnegative uptake and finite allocation loops. Actual transpiration remains supply limited. |
| `dslop`, `vert`, `gwater` | Use actual adjacent-layer values, signed paired transfers and exact groundwater recession. Hydrostatic equilibrium, reverse flow and combined storage balance are tested. |
| `fk`, `fpsim` | Consistent CH or Mualem-VG conductivity/retention. Forward/inverse tests include wet/dry endpoints; zero water does not reappear on a later call. |
| `interp`, `iter`, `srfpar`, `swge` | Retained reference helpers. Corrected interpolation bounds, iteration signatures/storage and gradient limits, and missing evaporation constants. They are not selected by the present driver or claimed to have independent physical validation. |

## Optional hydrology corrections

The previously listed ET coordinate, vapour-pressure and convergence bugs,
undefined multilayer variables, incompatible helper calls and mixed soil
hydraulic parameterizations have been corrected. Soil state now uses a
layer dimension rather than assuming one layer. Layer rates are timestep
means in mm/day; ET demands enter in m/s. Soil evaporation and root uptake
use current unmet demand, and actual ET is counted once. Bypass, downslope
runoff and groundwater discharge return to the current surface step; deep
seepage is tracked separately.

The original example soil fractions of 5 were invalid. The example now
contains admissible illustrative values and is marked uncalibrated.
[Hydrology input and accounting notes](HYDROLOGY.md) define both supported
hydraulic models, units, initial states, process connections and assumptions.

## Validation limits

Passing these checks does not establish full physical validation. Observed
ET, snow and soil-profile data are still needed for calibration and field
validation. The energy-balance interface assumes a nearby standard meridian
and clear-sky emissivity at night because it has no timezone or measured-LW
input. Complex channel networks need independent benchmarks.

The restored roof outlet is the original simplified representation.
Individual downpipe capacities, locations and backwater effects require
additional input definitions and are not inferred from building edges.
These are stated model limits, not features removed during cleanup.

## Checks

`run_regression_tests` covers both routing schemes, the original roof outlet,
infiltration/runoff, channels, forcing and streamed output.
`run_hydrology_tests` covers ET, canopy/snow, both soil hydraulic models,
multiple layers and complete-driver water accounting for all supported
ET/canopy/soil combinations. `run_swmm_tests` checks direct exchange with the
compiled engine. These suites pass with MATLAB R2024b on Apple silicon.
`run_routing_parallel_tests` compares the MATLAB and threaded C++ solver,
and `run_parallel_tests` compares every saved numerical field between
serial and parallel runs. The tests cover worker hydrology and native routing.
`run_gpu_tests` compares CPU/Metal routing and complete runs, including
multilayer soil and CPU-worker hydrology. The actual SWMM suite also passes
with `Device='gpu'`. GPU precision differs from CPU FP64 and is documented
in [GPU setup and checks](GPU.md).
MATLAB code analysis reports no syntax errors; remaining notices include
legacy APIs, unused compatibility arguments and retained reference helpers.
Dependency review distinguishes runtime functions from public utilities and
the build script. These checks do not replace calibration or independent
physical benchmarks.
