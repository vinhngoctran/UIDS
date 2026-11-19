# Overland Flow-Urban Flood Model (OFM-Urban)

* Source codes: https://github.com/vinhngoctran/UIDS
* Developers: Vinh Ngoc Tran and Jongho Kim
* Contact address: School of Civil and Environmental Engineering, University of Ulsan, 93 Daehak-ro, Nam-gu, Ulsan (44610), Korea. Email: kjongho@ulsan.ac.kr; vinhtn@umich.edu
* Programing language: MATLAB 2021b


* OFM-Urban is an open-source Matlab-based software package for the % simulation of hydrological processes in complex watershed systems.

### First setup on a new computer (Windows, Linux, or macOS)

**Compile the SWMM interface once on each new computer before running a
configuration with `MODDRA=1`.** The compiled file is ignored by Git, so a
fresh clone needs a local build. After a successful build, reuse it for
subsequent simulations on that computer. Compilation is not required before
every run. Configurations with `MODDRA=0` can run without this build.

For threaded adaptive OFM routing, also compile `ofm_route` once using
`build_ofm` and a supported **C++** compiler. This works without Parallel
Computing Toolbox or a separate OpenMP installation. Cell hydrology can
additionally use MATLAB workers when Parallel Computing Toolbox is available.

On Apple silicon, `build_gpu` builds an optional Metal GPU routing backend.
Select it with `Device='gpu'`. It does not require CUDA or Parallel Computing
Toolbox. See the [GPU setup, precision and checks](docs/GPU.md) before use.

The interface is a MATLAB MEX file compiled for the local operating system
and CPU architecture. A Windows binary cannot be reused on Linux or macOS.
On a Mac, build using the MATLAB installation and architecture that will run
the model. See [MathWorks MEX platform compatibility](https://www.mathworks.com/help/matlab/matlab_external/platform-compatibility.html).

#### 1. Copy the complete repository and install a C/C++ compiler

Keep `src/`, `native/`, `third_party/swmm/`, `tests/` and the model inputs
together. The repository includes the SWMM C source. No separate SWMM
installation or CMake build is needed. Use a writable checkout because the
build writes its compiled module into `src/`.

Install a C compiler supported by your MATLAB release. Check the MathWorks
compiler lists for [Windows](https://www.mathworks.com/support/requirements/supported-compilers.html),
[Linux](https://www.mathworks.com/support/requirements/supported-compilers-linux.html)
or [macOS](https://www.mathworks.com/support/requirements/supported-compilers-mac.html)
for the required compiler version.

| System | Compiler setup before opening MATLAB |
| --- | --- |
| Windows PC | Install a supported MinGW-w64 compiler through MATLAB Add-Ons (search for **MinGW**), or install a supported Microsoft Visual Studio C/C++ toolchain. |
| Linux | Install a supported GCC version and the system C development tools using your distribution's package manager. |
| macOS | Install the supported Apple Xcode tools. For Command Line Tools, run `xcode-select --install` in Terminal. Use full Xcode if required by your MATLAB release's compiler list. |

For the Windows Add-Ons installation steps, see the
[MathWorks MinGW guide](https://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html).

#### 2. Select the compiler and build in MATLAB

Open MATLAB and set **Current Folder** to the repository root, which contains
`Main_OFM_Urban.m`, `src/`, `native/` and `third_party/`. Run these commands in
the MATLAB Command Window, not the operating system terminal:

```matlab
addpath(fullfile(pwd, 'src'));
mex -setup C
build_swmm;
% Threaded adaptive OFM solver:
mex -setup C++
build_ofm;
% Optional GPU routing on Apple silicon:
% build_gpu;
```

Select a supported **C** compiler if MATLAB offers more than one. MATLAB
retains the compiler selection for later builds. See
[MathWorks compiler selection](https://www.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html).

The build prints the path to `src/swmm_live.<mexext>`. The extension depends
on the local MATLAB platform; use `mexext` to check it. Do not rename a
compiled file to make it look like a binary for a different platform.

#### 3. Check the build before the first coupled run

Still in MATLAB at the repository root:

```matlab
expected = fullfile(pwd, 'src', ['swmm_live.' mexext]);
assert(isfile(expected), 'The SWMM interface has not been built here.');
assert(exist('swmm_live', 'file') == 3, 'The compiled interface is not on the MATLAB path.');
which swmm_live

addpath(fullfile(pwd, 'tests'));
run_regression_tests;
run_hydrology_tests;
run_swmm_tests;
run_routing_parallel_tests;
% Requires Parallel Computing Toolbox:
run_parallel_tests;
```

`which swmm_live` should point to this checkout's `src/` folder. The tests
should finish with **All OFM regression checks passed.**, **All ET, canopy
and layered soil checks passed.**, and **All direct SWMM and coupled exchange
checks passed.** The SWMM checks run the actual
compiled engine and exercise exchange in both directions.

The interface has been built and tested with MATLAB R2024b on an Apple
silicon Mac. Windows and Linux builds are provided through the same source
and build script; they have not yet been tested in this repository's current
validation run.

#### 4. Run the model

For later MATLAB sessions on the same computer, add `src/` to the path and
run your configuration. Select an existing `.inf` file whose input paths
match the new computer. Set `MODDRA=1` to enable the direct SWMM coupling.

```matlab
addpath(fullfile(pwd, 'src'));
Filename = 'your_model.inf'; % Replace with your existing configuration file.
Results = ofm_urban(Filename);
```

`Main_OFM_Urban.m` is another entry point; edit its `Filename` before running
it. Its default Test 01 configuration has sewer coupling disabled.

#### When to rebuild and how to recover a failed build

| Situation | Action |
| --- | --- |
| New computer, fresh clone, different operating system or CPU architecture | Run the setup and build steps on the destination computer. |
| Changes to `native/swmm_live.c` or the SWMM engine source | Rebuild before running the changed code. |
| Changes to `native/ofm_route.cpp` | Run `clear ofm_route; build_ofm` with the local C++ compiler. |
| MATLAB or compiler upgrade | Rebuild for the new environment before using the coupled model. |
| Repeated simulations with the same source and environment | Reuse the existing compiled module. |
| MATLAB cannot find a supported compiler | Install a compiler from the matching MathWorks list, restart MATLAB, then run `mex -setup C` and `build_swmm`. |
| `OFM:SWMMBuild` or a missing `swmm_live` message | Check Current Folder, add this checkout's `src/` to the path, and run `build_swmm`. |
| An invalid or incompatible MEX file message | Rebuild locally with the MATLAB installation that will run the model. |
| The compiled file is loaded, locked, or cannot be overwritten | Stop the simulation and restart MATLAB, then build before starting another run. Also check that `src/` is writable. |

### How to run

This task can be carried out through the following steps:

1. Download the source code using https://github.com/vinhngoctran/UIDS

2. Prepare the input files:
* Climate forcings (e.g. precipitation)
* DEM (digital elevation map)
* Land use/land cover
* Building footprint map
* Vegetation/LAI
* Albedo soil/vegetation
* Stormwater network (manhole, pipe, outfall)

3. Configure the Model Simulation

All information regarding input data and model setup is presented in the *.inf file (e.g., GN2020.inf). The required information to be provided includes:

* FOLDER: Source input folder
* TIMEOPT: Time step options, including start time and end time of the flood event, computational time step for OFM and SFM models, time step to save results
* MODOPT: Option to select simulation models. The current developed version includes models such as Channel routing model, drainage model, overland flow model, Conceptual-based Canopy interception, snowmelt, and evapotranspiration model, Energy balance-based Evapotranspiration model, Infiltration-runoff model, and soil water fluxes with vegetation model
* OUTLET: Parameters for the (channel) outlet of the watershed
* INFILTRATION: Parameters for the infiltration model
* RUNOFF: Parameters for the runoff model
* RAINFALL: Rainfall data source
* ROUGHNESS: Roughness (Manning) parameters
* DRAINAGE: Stormwater network data
* SOURCES: In-sources and out-sources within the watershed
* RESULTS: Name to save output files

Details of the format for each input file can be referenced in the examples within the Example Folder.

4. Run one configuration from the repository root. For `MODDRA=1`, first
complete the [new-computer setup](#first-setup-on-a-new-computer-windows-linux-or-macos).

```matlab
addpath('src');
Results = ofm_urban('Test_01.inf');
% To display the final maximum-depth map (requires Mapping Toolbox):
% Results = ofm_urban('Test_01.inf', struct('Plot', true));
```

`Main_OFM_Urban.m` runs one selected example and plots two point series.
Change `Filename` to select a different configuration.

### Source layout

`src/Contents.m` lists the model entry points, routing routines, optional
hydrology modules, and SWMM dependencies. The maintained adaptive solver is
`explicit.m`; `overland.m` provides the sequential scheme. Both use
`diffusive_face.m` for Manning fluxes.

The [function review](docs/FUNCTION_REVIEW.md) records a decision for every
original MATLAB source file. Distinct optional models and public utilities
are retained. `laivariant.m` and `processing.m` remain available for
preprocessing, even though the main simulation does not call them.

The selected solvers are consolidated in `explicit.m`, `overland.m` and
`diffusive_face.m`. Calls to the removed `explicit2`/`explicit3` aliases
should use `explicit`. The original roof-to-ground outlet behavior is
preserved in both routing schemes. Unused predecessors, unfinished
experiments and generated run files are listed with their removal reasons
in the review. The review also lists validation gaps in optional modules.

### Parallel execution

From the repository root, after running `build_ofm` once:

```matlab
addpath('src');
options = struct('Parallel','auto','NumWorkers',4);
Results = ofm_urban('your_model.inf',options);
disp(Results.Execution);
```

The adaptive solver (`TSTP2D=2`) has a compiled C++ implementation that
partitions shared-face calculations and cell updates between CPU threads.
Each face is calculated once. The cells gather those transfers after all
faces are finished, preserving the MATLAB solver's transfer order, roof
outlets, storage limits and boundary conditions. The native threads do not
call MATLAB APIs. Their count is capped by available hardware threads.

ET, canopy/snow, infiltration/runoff and soil calculations use `parfor`
blocks of independent active cells. Serial and worker execution call the
same `hydrology_cells` kernel. Static parameters stay in a
`parallel.pool.Constant` on workers; only current block state and forcing
are exchanged at each step. The default pool has up to four process workers.
Existing pools are reused and left open. These choices follow the MathWorks
guidance on [worker environments](https://www.mathworks.com/help/parallel-computing/choose-between-thread-based-and-process-based-environments.html)
and [reusing worker data](https://www.mathworks.com/help/parallel-computing/improve-parfor-performance.html).

| Option | Default and behavior |
| --- | --- |
| `Parallel` | `'auto'`: use native routing from 4096 active cells when built, and worker hydrology from `ParallelMinCells` when the toolbox is available. `true` or `'on'` bypasses size thresholds. `false` or `'off'` selects serial execution. |
| `NumWorkers` | Four when creating a pool or running native routing. An existing pool is reused; an explicit value limits how many of its workers participate. |
| `ParallelMinCells` | 8192 active cells for automatic worker hydrology. Smaller/light runs often spend more time scheduling workers than calculating. |
| `Routing` | `'auto'`, `'matlab'`, `'native'`, or `'gpu'`. Explicit native/GPU selection requires its local build and `TSTP2D=2`. `Parallel=false, Routing='native'` runs the native CPU solver with one thread. |
| `Device` | `'cpu'` by default. `'gpu'` selects the Apple Metal backend while optional hydrology continues on the CPU. |
| `PoolType` | `'Processes'` when starting a pool; `'Threads'` is also supported by the tested hydrology path. Existing pool type is retained. |

Automatic mode falls back to MATLAB routing when the compiled module is
missing, and to serial hydrology if the toolbox or pool is unavailable.
Explicit parallel hydrology requests report an error if they cannot run.
Check `Results.Execution.RoutingBackend`, `RoutingWorkers`,
`ParallelHydrology`, `PoolType` and `Workers` to see what actually ran.

The sequential scheme (`TSTP2D=1`), channel routing and SWMM exchange retain
their ordered updates. The SWMM engine itself remains on its existing build
and thread settings. Pool workers never call the SWMM MEX. Timesteps run in
order, and each step is still saved before the next one, without accumulating
an output queue or time history in RAM. Process workers add memory overhead;
reduce `NumWorkers` if RAM is limited. Gauge assignments are cached at
startup, and rainfall, boundary and output bookkeeping use grid operations.

For a strict serial comparison:

```matlab
serial = ofm_urban('your_model.inf',struct('Parallel',false,'Routing','matlab'));
```

`run_routing_parallel_tests` checks MATLAB/native equations, roof drainage,
wet/dry behavior, boundaries, conservation, one-cell-wide grids and thread
shutdown/recovery. `run_parallel_tests` compares every numeric result field
and surface-change grid in every saved step between serial and parallel runs.
Use `run_parallel_tests(128)` for a larger timing fixture. Pool startup, disk
speed, grid size and enabled models affect whole-run speed; routing-only
timings are not whole-model speedup claims.
Measured routing and whole-model fixture timings are recorded in the
[performance notes](docs/PERFORMANCE.md), together with reproduction commands.

### GPU execution on Apple silicon

Build once and select the GPU from MATLAB at the repository root:

```matlab
addpath('src');
mex -setup C++
build_gpu;
options = struct('Device','gpu','Parallel','auto','NumWorkers',4);
Results = ofm_urban('your_model.inf',options);
disp(Results.Execution);
```

This sends adaptive OFM routing (`TSTP2D=2`) to the Mac's Metal GPU. CPU
hydrology workers, channels, direct SWMM exchange and step-file saving keep
their existing roles. Set `Device='cpu'` to use CPU routing again. An explicit
GPU request reports an error if its platform/build is unavailable.

GPU storage and transfers use compensated FP32 pairs, while hydraulic powers
use FP32. They are not full FP64 arithmetic. Every GPU routing call checks
the water budget on the CPU, and tests compare GPU results with the CPU
solver. The [GPU guide](docs/GPU.md) records the precision rules, checks and
measured routing-only timings. On the tested M3 Ultra, the 1024 x 1024 routing
fixture was about 9.6 times faster than four CPU threads; whole-model gains
depend on hydrology, synchronization and disk costs.

### ET, canopy/snow and layered soil water

The optional hydrology implementation has been corrected, including the
ET coordinate order, humidity conversion and stability iterations;
canopy/snow temperature and water accounting; and the multilayer soil
variables, hydraulic curves and storage bounds. Combined model options
share one ET demand and count actual evaporation once. Bypass, downslope
runoff and groundwater discharge return to OFM within the same step.

`MODSWF=1` requires `MODINF=1` and either `MODINT=1` or `MODETS=1`. The soil
reader accepts scalar or per-layer parameters, initial soil potentials and
an explicit hydraulic-model selection: `imodel=0` for modified
Clapp-Hornberger, or `imodel=1` for Mualem-van Genuchten. Soil output arrays
use a third dimension for layers, without accumulating a time history.

See [hydrology inputs, units and water accounting](docs/HYDROLOGY.md) before
enabling these options. In particular, the weather column named `RH` is
specific humidity in **kg/kg**, and the example soil values are illustrative,
not calibrated site parameters.

### Buildings and roof drainage

Existing configurations using `DEM=999` to identify buildings retain the
original UIDS roof outlet. `SHAPE=999` can also identify roof cells when
actual roof elevations are available. Rainfall and stored roof water enter
the surface water balance. Roof-to-roof transfers use surface routing, while
each adjacent roof/ground face drains through the original expression:

```text
h = max(roof_water_depth - roof_depression_storage, 0)
Q = 0.6 * 0.5 * sqrt(2 * 9.81) * h^(3/2)  [m^3/s]
```

The same volume is removed from the roof and added to the neighbouring
ground cell. Combined withdrawals cannot drain more roof water than is
available. This is the original simplified downpipe/outlet representation
at building edges. The model does not infer individual downpipe locations,
diameters or connections to remote cells. Ground water does not flow onto
a dry roof through this outlet. Set roof impervious fractions and depression
storage in the input maps to represent the intended roof properties.

The legacy `DEM=999` marker is still reserved for roofs in the OFM driver.
Other high terrain elevations use ordinary land routing. See the function
review for the convention's limits and the roof regression checks.

### Direct OFM/SWMM coupling

For `MODDRA=1`, complete the
[new-computer setup](#first-setup-on-a-new-computer-windows-linux-or-macos)
before the first run. `build_swmm` compiles the included EPA SWMM 5.2.4
source and the OFM adapter into `src/swmm_live.<mexext>`. Later simulations
reuse that module. No separate DLL, binary-output reader, Python process or
SWMM installation is needed.

A single SWMM session stays open for the entire OFM run. At every surface
step, OFM passes the junction inflow vector directly to the engine. SWMM
advances its adaptive hydraulic steps to the exact OFM time and returns
node depths, link flows, outfall flows and flooded volumes as MATLAB arrays.
The adapter integrates flooding over every internal hydraulic step. OFM
adds that volume back to the mapped surface cells in the same surface step.
No inflow time-series files, hotstarts or binary-result reads are used for
synchronization. EPA binary time-series saving is disabled.

`DELTA_T` now controls the exchange interval. `DELTA_T1D` is the maximum
SWMM hydraulic step, not a separate file-exchange interval; the engine's
steps are also capped at each OFM endpoint. Fractional surface steps of at
least 0.001 seconds are supported. Use CMS units, whole-second start/end
times and one active, non-boundary surface cell per coupled junction.
Node and link values map by identifier rather than input-section order.

The prepared network is written and opened **once** under
`Results/<RES_GEN>/run_<id>/drainage/`, where the startup diagnostic report
also belongs. Relative external input paths are resolved against the original
SWMM input folder. Existing FLOW inflows on coupled junctions are replaced
by the direct OFM inflows. Other network forcing remains active. SWMM
surface ponding is disabled at coupled junctions because OFM stores that
water. The old Windows DLL, `read1D` MEX and interval-file writer have been
removed.

The coupling uses explicit exchange: inlet flow uses the available surface
water and sewer head at the start of the surface step. It does not iterate
the two solvers to convergence within a step. Check sensitivity to `DELTA_T`
when pressure or surface levels change rapidly. Result snapshots continue
to be saved after each computational step.

Additional snapshot fields are `SewerTimeSeconds`, `NodeFloodVolume` [m3
returned this step], `SewerInVolume`, `SewerReturnVolume`, `SewerOutVolume`
[cumulative m3 by junction/outfall], and `SewerVolume` [current network m3].
`node_surch` is the current step's average flooding rate in m3/s. It is not
replayed into the next step.

The engine is pinned to [EPA SWMM 5.2.4](https://github.com/USEPA/Stormwater-Management-Model/tree/v5.2.4).
Its source and local date-rounding correction are documented in
`third_party/swmm/README.md`; the direct interface is described in
`native/README.md`.

### Output written at every step

The model writes into a new `Results/<RES_GEN>/run_<id>/` directory for each
run. `Results.OutputDirectory` gives that directory. The output is:

* `metadata.mat`: grid, mask, coordinates, start/end times, timestep and planned step count.
* `step_000000001.mat`, `step_000000002.mat`, etc.: independent snapshots, saved after each completed computational step. Each contains `Results`, `Waterdepth` and `WaterVec` for that step only.
* `summary.mat`: final snapshot, maximum inundation, grid data, outlet peaks and volume totals.
* `MaxInundation.txt`: final maximum inundation grid.

`Results.Exp_flo` is now a **2D grid**, and the function returns the final
snapshot. It no longer returns a 3D depth history. Each snapshot records its
`Step`, `TimeSeconds`, `DateTime` and `StepsWritten`. Optional model fields
contain the current state. `ChannelDepth` stores channel depths when enabled.
`WaterVec` is the net surface depth change rate in m/s, not flow velocity.
`Exp_rain` is rainfall intensity in m/s and `Exp_vol` is cumulative infiltration
in metres. No-data cells in the `Exp_*` grids use -9999.

Snapshots are written every `DELTA_T`, including the final step. `NT_REP` no
longer controls saving. `NT_WAR` still controls console reporting, in seconds.
The solver keeps current grids and running maxima in memory. Rainfall,
weather and boundary records stay at their native time intervals and are
sampled as needed. Output histories and sewer inflow histories do not grow
in memory. Disk use and file count grow with the number of steps.

Load a snapshot or extract selected cells without assembling a 3D grid:

```matlab
step = load(fullfile(Results.OutputDirectory, 'step_000000010.mat'));
depth = step.Results.Exp_flo;

% [row column] pairs, with optional computational step numbers:
[t, h] = read_ofm_series(Results.OutputDirectory, [7 43; 7 63], ...
                        10:10:Results.StepsWritten);
plot(t/3600, h);
```

### Numerical changes and limits

* `TSTP2D=2` uses shared face fluxes, Manning conveyance driven by water surface slope, adaptive diffusion and total-outflow limits, and paired volume transfers. Face depth accounts for the higher bed so a dry obstruction cannot be crossed before overtopping. Sources are applied once before the routing substeps. Closed-domain routing conserves volume.
* `TSTP2D=1` retains the legacy sequential update with availability limits. Its old “implicit” label did not describe a global implicit equation solve. Use `TSTP2D=2` for the adaptive solver.
* Stage boundaries are clamped to nonnegative depths and held at their prescribed step-end stage during routing. Free boundaries retain zero depth. These boundaries can supply or remove water.
* Building cells retain the original roof-to-ground outlet law. Both routing schemes conserve the transferred volume and limit total roof drainage by available water. The sequential scheme also retains the original steep-slope outlet expression.
* Green-Ampt infiltration now uses the integrated ponded equation, accounts for the transition to ponding, and limits infiltration by available water and pervious area. `bf` is pervious-area cumulative infiltration in cm; `inf_rate` is cell-average m/s. Surface storage stays in `oflow`, so infiltration does not create a second surface store in `stor`.
* Runoff reservoirs integrate rates over `DELTA_T` and cap each withdrawal by storage. Channel transfers are committed within the current step; shared nodes are updated once. Inlet and outlet flow limits apply to fractional and large surface timesteps.
* SWMM receives direct inflow arrays every surface step and retains its hydraulic state. Flooded volumes return in the same step, and SWMM errors stop the run and close the session.

* The soil module now uses consistent hydraulic curves across initialization,
  conductivity and potential updates. Layer withdrawals and transfers respect
  available water and capacity. Groundwater recession integrates recharge
  within the step; bypass is tracked separately from groundwater recharge.

The checks below verify numerical behavior and water accounting, not field
calibration or full benchmark validation. The hydrology guide documents the
ET solar-clock and nighttime-radiation assumptions.
Diffusive-wave routing omits momentum effects, so Test 03 is outside its
intended scope.

The equation references are the US Army Corps of Engineers descriptions of
[diffusive-wave flow](https://www.hec.usace.army.mil/confluence/rasdocs/hecras/latest/technical-reference/hydraulic-equations/diffusive-wave-equation)
and [Green-Ampt infiltration](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/canopy-surface-infiltration-and-runoff-volume/infiltration/green-and-ampt-loss-model).

### Regression checks

From MATLAB at the repository root:

```matlab
addpath('tests');
run_regression_tests;
run_hydrology_tests;
% After build_swmm, exercise the actual sewer engine and coupled model:
run_swmm_tests;
```

The checks cover volume conservation, wetting/drying, depression storage,
water at rest, rainfall sources, roof drainage, native forcing intervals, infiltration,
reservoir balance, outlet limits, stage/free boundaries, zero-LAI canopy,
channel storage, SWMM input preparation, and snapshot sizes and times.
The direct SWMM checks exercise fractional synchronization times, identifier
mapping, flood return volumes, clock alignment, hydraulic volume balance,
error cleanup, and reopening sessions. They verify that synchronization
leaves startup file sizes and timestamps unchanged. A full coupled fixture
checks the surface exchange balance at every saved step. The general checks
also run a short section of the supplied Test 01 benchmark. Temporary
fixtures and result directories are removed after the checks.

The hydrology checks cover day/night ET, calm wind, vapour-pressure fallback,
canopy/snow conservation, CH and VG retention curves, multilayer flow and
storage bounds, and complete-driver water balances for each supported
ET/canopy/soil combination.

## Acknowledgments
Vinh Ngoc Tran, vinhtn@ornl.gov
