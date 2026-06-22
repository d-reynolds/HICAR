#!/bin/bash
# Example: full-physics, real-data alpine downscaling run with a nested domain
# (COSMO-style forcing over the Gaudergrat 250 m domain, plus a 150 m child).
# Mirrors the physics + dynamics suite used for hectometric HICAR simulations.
#
# Each line sets ONE namelist variable BY NAME in its group, independent of the
# model defaults; generate_examples.sh strips everything left at its default.
# See alpine_realdata.sh for notes on set_var. Per-nest values are comma lists.
out_file=$1
set_var() { "${PYTHON:-python3}" "$(dirname "$0")/../set_nml_var.py" "$out_file" "$1" "$2" --group "$3" || exit 1; }

# --- general run settings --------------------------------------------------
set_var nests 2 general
set_var parent_nest "0, 1" general

# --- run window / IO -------------------------------------------------------
set_var start_date     "'2017-02-14 00:00:00'" general
set_var end_date       "'2017-02-15 00:00:00'" general
set_var outputinterval 3600 output
set_var inputinterval  3600 forcing
set_var output_vars    "'all'" output

# --- domain / forcing files ------------------------------------------------
set_var init_conditions_file "'../domains/Gaudergrat_250m.nc','../domains/Gaudergrat_nested_250m.nc'" domain
set_var forcing_file_list    "'file_list_TestCase.txt'" forcing
set_var dx "250.0, 150.0" domain
set_var nz 40 domain

# --- physics suite ---------------------------------------------------------
set_var mp    "'morrison'" physics
set_var pbl   "'ysu'" physics
set_var lsm   "'noahmp'" physics
set_var sfc   "'revmm5'" physics
set_var water "'simple'" physics
set_var rad   "'rrtmg'" physics

# --- dynamics --------------------------------------------------------------
set_var RK3 .True. time_parameters
set_var flux_corr 1 adv_parameters
set_var h_order "5, 3" adv_parameters
set_var v_order 3 adv_parameters
set_var cfl_reduction_factor 1.6 time_parameters
set_var terrain_shading .True. rad_parameters
set_var Sx .True. wind
set_var smooth_wind_distance 500.0 wind

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
set_var vegtype_var "'landuse'" domain
set_var svf_var "'svf'" domain
set_var hlm_var "'hlm'" domain
set_var LU_Categories "'USGS'" lsm_parameters
