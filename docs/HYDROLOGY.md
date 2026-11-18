# Optional hydrology: inputs, state and water balance

The optional ET, canopy/snow and soil modules remain separate model options.
The corrections repair their numerical implementation and their connection
to OFM. They do not supply calibrated parameters for a particular site.
All 18 original local soil-process functions are retained.

## Selecting the models

| Configuration flag | Behavior |
| --- | --- |
| `MODETS=1` | Estimate potential ET from the surface energy balance. Actual surface evaporation is limited by surface water and rainfall. |
| `MODINT=1` | Partition precipitation into canopy, ground and snow; drain canopy storage; melt snow; calculate evaporation using Penman-Monteith. |
| Both flags enabled | Use the energy-balance ET demand for the canopy and ground together. Canopy evaporation is deducted before the ground receives the remaining demand. |
| `MODSWF=1` | Route infiltration through a layered soil profile and groundwater reservoir. Requires `MODINF=1` and at least one of `MODETS` or `MODINT`. |

Without `MODSWF`, the existing conceptual runoff reservoirs receive
infiltration. With `MODSWF`, those reservoirs are replaced by the soil
profile, not run a second time on the same water.

After surface evaporation, unmet ground ET demand is multiplied by the
pervious fraction. Vegetation cover `1-exp(-0.5*LAI)` partitions that demand
between root uptake and top-layer soil evaporation. This is an explicit
model assumption. Each withdrawal is limited by accessible soil water;
the root resistance calculation also limits transpiration.

## Weather and energy-balance ET

The weather-file column `RH` means **specific humidity in kg/kg** in this
code. Do not put relative humidity percentages into it. Temperature is in
degrees Celsius, radiation in W/m2, wind speed in m/s, and atmospheric and
vapour pressure in kPa.

`evapotrans` receives `[TA SR WD VP RH PA]`. A `VP` value of `-9999` requests
conversion from specific humidity and atmospheric pressure:

```text
epsilon = molecular_mass_water / molecular_mass_dry_air
VP = RH * PA / (epsilon + (1-epsilon)*RH)
```

The conversion follows the specific-humidity definition documented by
[FAO WaPOR](https://www.fao.org/aquastat/py-wapor/et_look_rsts/meteo.html).
An explicitly supplied `VP` must be nonnegative. A 0.01 kPa numerical floor
protects logarithms at the dry endpoint.

`intercsnow` receives `[TA RH SR LR WD]`, including measured incoming
longwave radiation. It retains its constant-air-density humidity conversion.
The two forcing vectors have different column orders; the driver selects
the correct columns by their weather-file names.

The energy-balance formulation retains the original
[HRMET](https://github.com/samzipper/HRMET) radiation and sensible-heat
approach. The implementation now uses latitude/longitude in their correct
order and east-positive longitude, includes seconds in the solar clock,
and solves atmospheric stability with a bounded root search. Calm wind uses
a 0.1 m/s floor. Equal air and surface temperatures give zero sensible heat.
Negative latent-heat residuals produce zero ET; condensation is not modeled.

Two existing input limitations need attention when preparing a run:

* The solar calculation assumes the nearest 15-degree standard meridian.
  `G_daylightSavings` in `src/common.m` supplies the additional clock offset.
  Check that setting against the timestamps of the weather data.
* This ET interface has no measured-longwave input. At night, shortwave
  observations cannot determine cloud cover, so it uses clear-sky emissivity.
  The standalone canopy model uses the supplied longwave observations.

The first output of `evapotrans` is **potential ET**, despite the historical
internal name `act_evap`. The second output is supply-limited surface
evaporation. Both are m/s. An optional third output reports radiation,
sensible heat, coordinates, vapour pressure and stability iterations.

## Soil input file

`SFLX` in the `.inf` configuration selects the soil parameter file.
`modswfip` accepts whitespace-separated records and `%`, `#` or `!` comments.
Record names are case-insensitive. Duplicate records and invalid numeric
values cause errors. The reader does not evaluate MATLAB expressions.

`nlayer` selects the profile size. The fields `thick`, `stonef`, `frelden`,
`tini`, `qffc`, `qfpar`, `psim_init`, and each `par1` through `par10` accept
either one value for every layer or `nlayer` values. For example:

```text
nlayer 2
thick 300 400          % mm
stonef 0.10 0.15       % fractions
psim_init -33 -100     % kPa
par1 0.45 0.40         % saturated volumetric water content
```

These records supplement the other required records shown in
`Examples/Input/SOILFLUX.txt`; the fragment above is not a complete input.
Layer thickness must be positive and stone fractions must be in `[0,1)`.
Soil matric potentials are nonpositive kPa; `psicr` is a leaf potential in
MPa. Legacy positive suction magnitudes for `psicr` and CH `par4` are
normalized to negative values.

### Hydraulic parameterization

One model selection applies consistently to retention, its inverse and
hydraulic conductivity. The conventions follow the separate CH and
Mualem-van Genuchten branches in the
[LWFBrook90 reference source](https://github.com/pschmidtwalter/LWFBrook90R/blob/master/src/md_brook90.f95).

| Record | `imodel=0`: modified Clapp-Hornberger, default | `imodel=1`: Mualem-van Genuchten |
| --- | --- | --- |
| `par1` | Saturated water content, fraction | Saturated water content, fraction |
| `par2` | Field-capacity water content, fraction | Derived from field-capacity conductivity |
| `par3` | Field-capacity conductivity, mm/day | Field-capacity conductivity, mm/day, between 0 and `par6` |
| `par4` | Field-capacity matric potential, kPa | Derived field-capacity matric potential |
| `par5` | Derived field-capacity saturation fraction | Derived field-capacity effective saturation |
| `par6` | Derived saturated conductivity, mm/day | Saturated conductivity, mm/day |
| `par7` | Derived near-saturation coefficient, kPa | Alpha, 1/m, positive |
| `par8` | Derived near-saturation coefficient, dimensionless | n, greater than 1 |
| `par9` | Positive CH retention exponent b | Mualem tortuosity exponent, commonly 0.5 |
| `par10` | Saturation fraction at start of near-saturation branch | Residual volumetric water content, fraction |

All ten records remain in the file for compatibility. Derived values are
recomputed in the solver. For CH, field-capacity water content must be below
saturated water content, and `par10` must lie between field-capacity
saturation and 1. For VG, residual content must be below saturated content.
The illustrative example file uses CH. Its former water contents of 5 were
invalid fractions and have been replaced with admissible **uncalibrated**
values. Use measured or calibrated site parameters for model studies.

Other process controls retain their existing meanings:

* `ilayer` and `infexp` distribute infiltration with depth. `infexp=0`
  sends all infiltration to the top layer; `infexp=1` distributes it by
  thickness over the first `ilayer` layers.
* `qffc` and `qfpar` control wetness-dependent bypass. `bypar=0` disables
  that fraction, but water exceeding soil capacity still returns to OFM.
* `dslope` is in radians and `length` in metres for downslope drainage.
  `drain` multiplies bottom-layer gravity drainage.
* `gsc` is the groundwater recession coefficient in 1/day. `gsp` is the
  fraction of groundwater discharge lost as deep seepage. At `gsc=0`,
  new bottom drainage passes through without adding groundwater storage.
* The canopy/root parameters in the example file control root growth,
  root distribution and plant resistance. `qlayer` is retained for the
  reference surface-runoff helper; it does not activate another runoff model.

| Optional record | Default | Unit or meaning |
| --- | --- | --- |
| `imodel` | 0 | Hydraulic model, 0 or 1 |
| `psim_init` | -500 | Initial matric potential, kPa, scalar or per layer |
| `gwat_init` | 0 | Initial groundwater, mm |
| `max_step` | 300 | Maximum soil substep, seconds |
| `dpsimx` | 0.001 | Interface potential differences below this kPa value are treated as zero |
| `bypar` | 1 | Enable wetness-dependent bypass, 0 or 1 |

The soil integrator also limits substeps using available water and outgoing
fluxes. Signed transfers remove and add identical volumes between layers.
The lower layer supplies groundwater recharge; bypass is not added to
groundwater. Groundwater recession uses the exact linear-reservoir solution
for constant recharge within each substep. Check sensitivity to `max_step`
and the OFM timestep for a new parameter set.

## Outputs and accounting

Each saved step contains current state, not growing histories.

| `Results` field | Shape and unit |
| --- | --- |
| `psim` | Row x column x layer, kPa |
| `vrfli`, `infli`, `byfli`, `dsfli`, `ntfli`, `trani` | Row x column x layer, timestep-mean mm/day |
| `swat`, `gwat` | Row x column, total profile/groundwater storage in mm |
| `gwfl`, `seep` | Row x column, timestep-mean mm/day |
| `potential_et`, `act_evap`, `soil_evap`, `transpiration` | Row x column, m/s; soil fields require `MODSWF` |
| `total_evap`, `soil_et_total`, `seep_total` | Row x column, cumulative m |
| `soil_balance` | Row x column, soil-plus-groundwater balance residual for this step, mm |

MATLAB omits a trailing singleton dimension for one-layer profiles.
Positive `vrfli(i)` flows down from layer i; negative values flow upward.
`ntfli` is the net layer storage change rate. `act_evap` and `total_evap`
include canopy, surface and actual soil/root losses, each counted once.
`soil_et_total` is a component of `total_evap`, not an additional loss.

Bypass runoff, downslope runoff and groundwater discharge return to the
same OFM surface cell within the current step. Deep seepage leaves the
modeled system. `v_infil` records gross surface infiltration; some of that
water may return through bypass, so it is not the change in soil storage.
The soil balance includes infiltration, soil and groundwater storage,
evaporation, transpiration, bypass, downslope discharge and deep seepage.

## Checks and scope

Run `run_hydrology_tests` after adding `tests` to the MATLAB path. The suite
checks day/night and calm-wind ET, zero sensible heat, vapour-pressure
fallback, canopy/snow balance, dry and saturated soil bounds, both hydraulic
curves, unequal layers, signed flow, hydrostatic equilibrium, and the full
driver with each supported ET/canopy/soil combination. Complete driver
fixtures close the combined water balance, including initial soil water.

These checks establish execution, numerical bounds and water accounting.
They are not a field calibration or validation against observed ET,
snowmelt or soil-moisture profiles. Complex channel networks also need
independent benchmark data. Roof drainage retains the original adjacent
cell outlet representation described in the README; individual downpipe
locations, capacities and backwater behavior require additional inputs.
