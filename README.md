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

4. Run Simulation
The model can be run using the MATLAB GUI (with visualization function). The steps to run include:
* Run MATLAB
* Open the 'Main_OFM_Urban.m' code
* Specify the configuration file for model setup (*.inf file) in line 29
* Run the code

## Acknowledgments
Vinh Ngoc Tran, vinhtn@umich.edu
