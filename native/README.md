# Direct SWMM interface

`swmm_live.c` embeds the pinned EPA SWMM engine in a MATLAB MEX module.
Build it with `addpath('src'); build_swmm`. Source files are compiled from
`third_party/swmm`; the generated module is ignored by Git.

The MATLAB `drainage` handle class owns the session and closes it on deletion.
The OFM driver also installs an `onCleanup` handler for failed runs. One
session can be active at a time because the EPA engine uses global state.

Native commands:

* `state = swmm_live('open', inputFile, reportFile, nodeIDs, linkIDs, outfallIDs)`
  opens the prepared network once, maps identifiers to indices, disables
  report/time-series saving, and calls `swmm_start(0)` once. IDs are cell
  arrays of character vectors. Values use SI units with FLOW_UNITS CMS.
* `state = swmm_live('advance', targetSeconds, inflow)` writes the inflow
  vector through `swmm_setValue(swmm_NODE_LATFLOW, ...)`, then calls
  `swmm_step` until the absolute target time. No file functions occur in
  the exchange adapter. SWMM's binary time-series saving stays disabled.
* `swmm_live('close')` ends and closes the engine once and releases maps.
  Closing an already closed session is safe. Engine errors also close it.

EPA's public `swmm_stride` accepts only integer seconds. The adapter includes
`swmm5.c` in its translation unit so it can set the private `RoutingDuration`
deadline for each exchange. This clips the last internal hydraulic step to
the synchronization time while preserving adaptive step selection. OFM
uses absolute targets (`t*DELTA_T`) to avoid accumulating clock error.

`NodeFloodVolume` is the sum of nonnegative SWMM node overflow rates times
their actual internal hydraulic step lengths, converted to m3. It is not an
endpoint rate multiplied by the OFM timestep. `OutVolume` is integrated over
the same steps. Coupled node ponding areas are set to zero before starting
because the surface model stores flooded water. All output arrays are copied
into MATLAB-owned memory, and each exchange retains only the current state.

SWMM still opens its network and initializes its output/report machinery at
startup, and final cleanup can write engine diagnostics. These files do not
carry exchange data. Normal OFM step snapshots remain independent of the
coupling interface. Networks configured with additional rainfall/runoff or
external routing files can still use those native SWMM inputs and outputs.

The direct interface uses explicit exchange, without a nonlinear iteration
between the two models. Smaller OFM timesteps give more frequent pressure
feedback. The pinned engine has one documented date-rounding correction in
`third_party/swmm/project.c`.
