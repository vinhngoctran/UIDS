# Overland Flow-Urban Flood Model (OFM-Urban)

* Source codes: https://github.com/vinhngoctran/UIDS
* Developers: Vinh Ngoc Tran and Jongho Kim
* Contact address: School of Civil and Environmental Engineering, University of Ulsan, 93 Daehak-ro, Nam-gu, Ulsan (44610), Korea. Email: kjongho@ulsan.ac.kr; vinhtn@umich.edu
* Programing language: MATLAB 2021b


* OFM-Urban is an open-source Matlab-based software package for the % simulation of hydrological processes in complex watershed systems.

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

4. Run one configuration from the repository root:

```matlab
addpath('src');
Results = ofm_urban('Test_01.inf');
% To display the final maximum-depth map (requires Mapping Toolbox):
% Results = ofm_urban('Test_01.inf', struct('Plot', true));
```

`Main_OFM_Urban.m` runs one selected example and plots two point series.
Change `Filename` to select a different configuration.

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
* Elevations at or above 999 metres use the same Manning law as other elevations. The former special roof/weir branch based on the literal elevation 999 has been removed. Existing roof elevations still act as raised terrain; a separate roof drainage law would need explicit roof parameters.
* Green-Ampt infiltration now uses the integrated ponded equation, accounts for the transition to ponding, and limits infiltration by available water and pervious area. `bf` is pervious-area cumulative infiltration in cm; `inf_rate` is cell-average m/s. Surface storage stays in `oflow`, so infiltration does not create a second surface store in `stor`.
* Runoff reservoirs integrate rates over `DELTA_T` and cap each withdrawal by storage. Channel transfers are committed within the current step; shared nodes are updated once. Inlet and outlet flow limits apply to fractional and large surface timesteps.
* SWMM receives one interval-average inflow vector. Each coupling interval uses separate input/output hotstart files. SWMM errors now stop the run. Use CMS flow units and a `DELTA_T1D` that is a whole number of seconds and a whole multiple of `DELTA_T`.

The supplied `swmm5.dll` and `read1D.mexw64` are Windows binaries. Running
sewer coupling on another platform requires compatible SWMM and `read1D`
binaries. The sewer engine and the existing SVAT soil-water module have not
been validated by these changes. The SVAT input reader still describes a
single soil layer. The changes below are numerical regression checks, not
calibration or full benchmark validation. Diffusive-wave routing omits
momentum effects, so Test 03 is outside its intended scope.

The equation references are the US Army Corps of Engineers descriptions of
[diffusive-wave flow](https://www.hec.usace.army.mil/confluence/rasdocs/hecras/latest/technical-reference/hydraulic-equations/diffusive-wave-equation)
and [Green-Ampt infiltration](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/canopy-surface-infiltration-and-runoff-volume/infiltration/green-and-ampt-loss-model).

### Regression checks

From MATLAB at the repository root:

```matlab
addpath('tests');
run_regression_tests;
```

The checks cover volume conservation, wetting/drying, depression storage,
water at rest, rainfall sources, native forcing intervals, infiltration,
reservoir balance, outlet limits, stage/free boundaries, zero-LAI canopy,
channel storage, SWMM input construction, and snapshot sizes and times.
They also run a short section of the supplied Test 01 benchmark. Temporary
fixtures and result directories are removed after the checks.

## Acknowledgments
Vinh Ngoc Tran, vinhtn@umich.edu
