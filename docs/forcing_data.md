
# Forcing Data

For reading in forcing data, HICAR accepts a .txt file specified by the namelist variable 'forcing_file_list' in the 'forcing' section of the namelist. This .txt file should be a list of full file paths to the forcing data surrounded by full quotes, with one file path per line, i.e.:

```text
    "/your/path/to/forcing/file1.nc"
    "/your/path/to/forcing/file2.nc"
    "/your/path/to/forcing/file3.nc"
    "/your/path/to/forcing/file4.nc"
    "/your/path/to/forcing/file5.nc"
    ...
```

Such a forcing file list can be generated using the helper script in `helpers/filelist_script.sh`. Example:

```bash
./filelist_script.sh "/your/path/to/forcing/file*.nc" file_list_ouput.txt
```

As seen in the example, forcing files must be netCDF files. Please ensure that the netCDF library which HICAR was compiled with supports parallel I/O for the netCDF file type which you are using. Parallel netCDF built on HDF5 only supports the netCDF-4 data type, while PnetCDF can only read classic netCDF data types.

## Variables

Information about possible forcing variables can be obtained by asking the HICAR executable about a specific forcing variable. This can be done as:

`./HICAR -v pvar`

in the case of the pressure variable, giving the output:

```text
 ---------------------------------------------------------------------------------------------------
 Namelist Variable:           pvar
     Namelist Group:.........|Forcing
     Description:............|Name of the pressure variable in forcing file
     Units:..................|Pa
     Dimensions:.............|[T, Z, Y, X]
```

A list of all possible forcing variables can be found by generating the default namelist (see [namelist_options.md](namelist_options.md)), then looking under the "forcing" group.

A partial list of the available forcing variables is below, with optional variables in square brackets:

```Text
qvvar         = Water Vapor mixing ratio                  (kg/kg)
tvar          = Air Temperature                           (K with an optional offset)
pvar          = Pressure                                  (Pa)
uvar          = East-West wind                            (m/s)
vvar          = North-South wind                          (m/s)
[wvar]        = Vertical wind                             (m/s)
[hgtvar]      = Terrain Height                            (m)
zvar          = 3D model level heights                    (m)
latvar        = Latitude on mass (P/T/etc.) grid          (degrees)
lonvar        = Longitude on mass (P/T/etc.) grid         (degrees)
[ulat]        = Latitude on the EW-wind grid              (degrees)
[ulon]        = Longitude on the EW-wind grid             (degrees)
[vlat]        = Latitude on the NS-wind grid              (degrees)
[vlon]        = Longitude on the NS-wind grid             (degrees)
[sst_var]     = Sea Surface Temperature                   (K)
[swdown_var]  = Shortwave down at the surface             (W/m^2)
[lwdown_var]  = Longwave down at the surface              (W/m^2)
[shvar]       = Sensible heat flux from the surface       (W/m^2)
[lhvar]       = Latent heat flux from the surface         (W/m^2)
[pbl_var]     = Specified height of PBL                   (m)
```

## Variable data structure

Forcing variables must conform to a known dimension ordering, as follows:

4-D variables (spatial):
[time, Z, Y, X]

3-D variables (spatial):
[Z, Y, X]

2-D variables (spatial):
[Y, X]

Input fields in 3 spatial dimensions may either be 4D (with time), or 3D. If they are 3D, there must be a single time step per input file. Latitude and longitude may also be provided as 1D data. Information about how input variables have been interpreted by the model will be output towards the beginning of the model output, where applicable.

Spatial 2D or 3D forcing variables should always cover the full domain. HICAR cannot read in boundary-only forcing files, and expects that a spatially 3D forcing variable exists in 3 spatial dimensions.

Information about the expected ordering and meaning of a variable's dimensions can be obtained by asking the HICAR executable about a specific forcing variable, such as:

`./HICAR -v pvar`

in the case of the pressure variable.
