#!/bin/bash
# Example: full-physics, real-data alpine downscaling run (COSMO-style forcing
# over the Gaudergrat 250 m domain). Mirrors the physics + dynamics suite used
# for hectometric HICAR simulations.
#
# Only sets entries that differ from the model defaults; generate_examples.sh
# then strips everything left at default. Each line uses the same `sed` pattern
# as the test-case generators (tests/Test_Cases/input/nml_gen_scripts).
out_file=$1
s() { sed -i'.bak' "$1" "$out_file"; }

# --- run window / IO -------------------------------------------------------
# leading space anchors start_date so it does not also match re*start_date*
s "s/ start_date = ''/ start_date = '2017-02-14 00:00:00'/g"
s "s/end_date = ''/end_date = '2017-02-15 00:00:00'/g"
s "s/outputinterval = 3600/outputinterval = 3600/g"
s "s/inputinterval = 3600/inputinterval = 3600/g"
s "s/output_vars = ''/output_vars = 'qv','temperature','precipitation','swe','snow_height','hfls','hfss','u','v','w'/g"

# --- domain / forcing files ------------------------------------------------
s "s/init_conditions_file = ''/init_conditions_file = '..\/domains\/Gaudergrat_250m.nc'/g"
s "s/forcing_file_list = ''/forcing_file_list = 'forcing_file_list.txt'/g"
s 's/dx = 0.0/dx = 250.0/g'
s 's/nz = 500/nz = 40/g'

# --- physics suite ---------------------------------------------------------
s "s/ mp = 'none'/ mp = 'morrison'/g"
s "s/ pbl = 'none'/ pbl = 'ysu'/g"
s "s/ lsm = 'none'/ lsm = 'noahmp'/g"
s "s/ sfc = 'none'/ sfc = 'revmm5'/g"
s "s/ water = 'none'/ water = 'simple'/g"
s "s/ rad = 'none'/ rad = 'rrtmg'/g"

# --- dynamics --------------------------------------------------------------
s 's/RK3 = .False./RK3 = .True./g'
s 's/flux_corr = 0/flux_corr = 1/g'
s 's/h_order = 1/h_order = 3/g'
s 's/v_order = 1/v_order = 3/g'
s 's/cfl_reduction_factor = 0.9/cfl_reduction_factor = 1.6/g'
s 's/terrain_shading = .False./terrain_shading = .True./g'
s 's/Sx = .False./Sx = .True./g'

# --- forcing variable names (COSMO-style) ----------------------------------
s "s/ time_var = ''/ time_var = 'time'/g"
s "s/ pvar = ''/ pvar = 'P'/g"
s "s/ tvar = ''/ tvar = 'T'/g"
s "s/ qvvar = ''/ qvvar = 'QV'/g"
s "s/ uvar = ''/ uvar = 'U'/g"
s "s/ vvar = ''/ vvar = 'V'/g"
s "s/ hgtvar = ''/ hgtvar = 'HSURF'/g"
s "s/ zvar = ''/ zvar = 'HFL'/g"
s "s/ latvar = ''/ latvar = 'lat_1'/g"
s "s/ lonvar = ''/ lonvar = 'lon_1'/g"
s 's/qv_is_spec_humidity = .False./qv_is_spec_humidity = .True./g'
s 's/t_is_potential = .True./t_is_potential = .False./g'

# --- domain variable names -------------------------------------------------
s "s/ lat_hi = ''/ lat_hi = 'lat'/g"
s "s/ lon_hi = ''/ lon_hi = 'lon'/g"
s "s/ hgt_hi = ''/ hgt_hi = 'topo'/g"
s "s/ landvar = ''/ landvar = 'landmask'/g"
s "s/ vegtype_var = ''/ vegtype_var = 'landuse'/g"
s "s/LU_Categories = 'MODIFIED_IGBP_MODIS_NOAH'/LU_Categories = 'USGS'/g"