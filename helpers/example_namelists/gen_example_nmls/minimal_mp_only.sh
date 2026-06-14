#!/bin/bash
# Example: a minimal run — dynamics + microphysics only, no land-surface, PBL,
# radiation or surface-layer schemes (all left at their 'none' defaults). Useful
# as a lightweight smoke configuration or a starting point to add schemes one at
# a time.
#
# Only sets entries that differ from the model defaults; generate_examples.sh
# strips the rest.
out_file=$1
s() { sed -i'.bak' "$1" "$out_file"; }

# --- run window / IO -------------------------------------------------------
# leading space anchors start_date so it does not also match re*start_date*
s "s/ start_date = ''/ start_date = '2017-02-14 00:00:00'/g"
s "s/end_date = ''/end_date = '2017-02-14 06:00:00'/g"
s "s/output_vars = ''/output_vars = 'qv','temperature','precipitation','u','v','w'/g"

# --- domain / forcing files ------------------------------------------------
s "s/init_conditions_file = ''/init_conditions_file = '..\/domains\/Gaudergrat_250m.nc'/g"
s "s/forcing_file_list = ''/forcing_file_list = 'forcing_file_list.txt'/g"
s 's/dx = 0.0/dx = 250.0/g'
s 's/nz = 500/nz = 40/g'

# --- physics: microphysics only (pbl/lsm/rad/sfc/water stay 'none') --------
s "s/ mp = 'none'/ mp = 'morrison'/g"

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
