#!/bin/bash
# Example: a minimal run — dynamics + microphysics only, no land-surface, PBL,
# radiation or surface-layer schemes (all left at their 'none' defaults). Useful
# as a lightweight smoke configuration or a starting point to add schemes one at
# a time.
#
# Each line sets ONE namelist variable BY NAME in its group, independent of the
# model defaults; generate_examples.sh strips everything left at its default.
# See alpine_realdata.sh for notes on set_var.
out_file=$1
set_var() { "${PYTHON:-python3}" "$(dirname "$0")/../set_nml_var.py" "$out_file" "$1" "$2" --group "$3" || exit 1; }

# --- run window / IO -------------------------------------------------------
set_var start_date  "'2017-02-14 00:00:00'" general
set_var end_date    "'2017-02-14 06:00:00'" general
set_var output_vars "'qv','temperature','precipitation','u','v','w'" output

# --- domain / forcing files ------------------------------------------------
set_var init_conditions_file "'../domains/Gaudergrat_250m.nc'" domain
set_var forcing_file_list    "'file_list_TestCase.txt'" forcing
set_var dx 250.0 domain
set_var nz 40 domain

# --- physics: microphysics only (pbl/lsm/rad/sfc/water stay 'none') --------
set_var mp "'morrison'" physics

# --- forcing variable names (COSMO-style) ----------------------------------
set_var time_var "'time'" forcing
set_var pvar  "'P'" forcing
set_var tvar  "'T'" forcing
set_var qvvar "'QV'" forcing
set_var uvar  "'U'" forcing
set_var vvar  "'V'" forcing
set_var hgtvar "'HSURF'" forcing
set_var zvar   "'HFL'" forcing
set_var latvar "'lat_1'" forcing
set_var lonvar "'lon_1'" forcing
set_var qv_is_spec_humidity .True. forcing
set_var t_is_potential .False. forcing

# --- domain variable names -------------------------------------------------
set_var lat_hi  "'lat'" domain
set_var lon_hi  "'lon'" domain
set_var hgt_hi  "'topo'" domain
set_var landvar "'landmask'" domain
