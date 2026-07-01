# Generating Static Data

Here we describe how to generate the time-invariant input data that describes the domain state (the "static data"). Example static data for running a 1-day simulation can be found under [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data)


## Cutting a domain from existing static data

The ESMF-based script [Regrid_script.sh](../helpers/domains/Regrid_script.sh) takes a 2D lat/lon grid stored in a netcdf file as input and regrids some existing static data to it.
For all practical purposes, the only static data _needed_ is elevation and land use type.
Thus, this script uses a topographic data and land use data file, both of them containing lat/lon variables, and regrids their data to the target grid provided as input.
The topographic and land use data should ideally be at a higher resolution than the target lat/lon grid.
Please read the script header for more information about running it.

## Terrain-shaded radiation parameters

To generate a HICAR domain file containing all variables needed for running with terrain-shaded radiation, the python script gen_HICAR_dom.py, located in helpers/domains/, can be used.

HICAR relies on pre-computed static data to speed up some of it’s online calculations. To generate a HICAR domain file, an existing netCDF file with lat, lon, and a DEM is required; landuse categories and a land mask are optional (the land mask is derived from landuse when present, and the domain otherwise defaults to all-land). By default the input lat/lon variables are read as **lat**/**lon** and the terrain (DEM) variable as **topo**, but these names are configurable via the `lat_var`, `lon_var`, and `topo_var` arguments near the top of `helpers/domains/gen_HICAR_dom.py`. Additionally, a larger extent DEM of the same resolution is needed to generate parameters for terrain-shading of radiation. I.e., if you have a 50m resolution domain, a larger DEM with an extent ~20km beyond the boundaries of the target domain is also needed.

Once you have these two NetCDF files, you can use a python script to generate the rest of the variables used by HICAR.

First, install the conda environment located in the HICAR_dom.yml file found in helpers/domains

```bash
conda env create -f HICAR_dom.yml
```
Once this environment is installed, activate it with:
```bash
conda activate HICAR_dom
```

Now you will need to install HORAYZON, a python package to efficiently calculate the horizon line matrix and sky view factor (Steger et al., 2022). To do so, type:

```bash
git clone https://github.com/ChristianSteger/HORAYZON.git
cd HORAYZON
python -m pip install .
```

Your python environment should now be complete. To generate the domain file, open
`helpers/domains/gen_HICAR_dom.py` and edit the paths for the target domain,
radiation domain, and output domain. `HICAR_Domain.py` and `ProjHelpers.py`, both
contained in the `helpers/domains` directory, must be in the same directory as
`gen_HICAR_dom.py`.

Now run:

```bash
cd helpers/domains
python gen_HICAR_dom.py
```

## netCDF data type

Please ensure that the netCDF library which HICAR was compiled with supports parallel I/O for the netCDF file type which you are using. Parallel netCDF built on HDF5 only supports the netCDF-4 data type, while PnetCDF can only read classic netCDF data types.
