# Generating Static Data

HICAR relies on pre-computed static data to speed up some of it’s online calculations. To generate a HICAR domain file, an existing NetCDF file with lat, lon, DEM, landuse categories, and a land mask is needed. The lat and lon variables must be named **lat** and **lon**, and the terrain variable must be named **topo**. Additionally, a larger extent DEM of the same resolution is needed to generate parameters for terrain-shading of radiation. I.e., if you have a 50m resolution domain, a larger DEM with an extent ~20km beyond the boundaries of the target domain is also needed.

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

Your python environment should now be complete. To generate the domain file, open gen_HICAR_dom.py and edit the paths for the target domain, radiation domain, and output domain. HICAR_Domain.py and ProjHelpers.py, both contained in the helpers/domains directory, must be in the same directory as gen_HICAR_dom.py.

Now run:

```bash
python gen_HICAR_dom.py
```

## netCDF data type

Please ensure that the netCDF library which HICAR was compiled with supports parallel I/O for the netCDF file type which you are using. Parallel netCDF built on HDF5 only supports the netCDF-4 data type, while PnetCDF can only read classic netCDF data types.
