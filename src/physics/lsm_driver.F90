!>----------------------------------------------------------
!! This module provides a wrapper to call various land surface models
!!
!! It sets up variables specific to the LSM to be used including both
!! history variables not currently stored in the domain level data
!! structure, and runtime parameters
!!
!! The main entry point to the code is lsm(domain,options,dt,model_time)
!!
!! <pre>
!! Call tree graph :
!!  lsm_init->[ allocate_noah_data,
!!              external initialization routines]
!!  lsm->[  sat_mr,
!!          calc_exchange_coefficient,
!!          external LSM routines]
!!
!! High level routine descriptions / purpose
!!   lsm_init           - allocates module data and initializes physics package
!!   lsm                - sets up and calls main physics package
!!  calc_exchange_coefficient - calculates surface exchange coefficient (for Noah)
!!  allocate_noah_data  - allocate module level data for Noah LSM
!!  apply_fluxes        - apply LSM fluxes (e.g. sensible and latent heat fluxes) to atmosphere
!!
!! Inputs: domain, options, dt, model_time
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!!      model_time      = time since beginning date (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module land_surface
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
    use module_sf_noahmpdrv, only : noahmplsm, noahmp_init
    use module_water_simple, only : water_simple
    use module_water_lake,   only : lake, lakeini, nlevsoil, nlevsnow, nlevlake
    use mod_atm_utilities,   only : sat_mr, calc_Richardson_nr, calc_solar_elevation
    use time_object,         only : Time_type
    use icar_constants,      only : kVARS, kLSM_NOAH, kLSM_NOAHMP, kSM_FSM
    use options_interface,   only : options_t
    use domain_interface,    only : domain_t
    use ieee_arithmetic
    use data_structures
    use mod_wrf_constants,   only : gravity, KARMAN, cp, R_d, XLV, rcp, STBOLT, epsilon
#ifdef FSM
    use module_sf_FSMdrv,   only : sm_FSM_init, sm_FSM
#endif
    implicit none

    private
    public :: lsm_init, lsm, lsm_var_request, lsm_apply_fluxes

    ! Noah LSM required variables.  Some of these should be stored in domain, but tested here for now
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

    ! LOTS of variables required by Noah, placed here "temporarily", this may be where some stay.
    ! Keeping them at the module level prevents having to allocate/deallocate every call
    ! also avoids adding LOTS of LSM variables to the main domain datastrucurt
    real,allocatable, dimension(:,:)    :: SMSTAV,SFCRUNOFF,UDRUNOFF, SW,                           &
                                           SNOWC, ACSNOW, ACSNOM, SNOALB,           &
                                           QGH, GSW, ALBEDO, ALBBCK, Z0, qz0,              &
                                           EMBCK, QSFC, CPM, SR, CHS, CHS2, CQS2,                   &
                                           VEGFRAC, SHDMIN,SHDMAX,SNOTIME,SNOPCX,&
                                           POTEVP, NOAHRES,z_atm,Ri,chklowq,  &
                                           flx4_2d,fvb_2d,fbur_2d,fgsn_2d, &
                                           nmp_snow, nmp_snowh

    logical :: MYJ, FRPCPN,ua_phys,RDLAI2D,USEMONALB
    real,allocatable, dimension(:,:,:)  :: SMCREL
    real,allocatable, dimension(:,:)    :: windspd
    real,allocatable, dimension(:)      :: Zs,DZs
    real :: XICE_THRESHOLD
    integer :: ITIMESTEP, update_interval, cur_vegmonth

!     real, parameter :: kappa=0.4 ! this should be karman from data_structure
    real, parameter :: freezing_threshold=273.15
    real, parameter :: SMALL_PRESSURE=0.1 ! note: 0.1Pa is very small 1e-10 wouldn't affect a single-precision float
    real, parameter :: SMALL_QV=1e-10
    real, parameter :: MAX_EXCHANGE_C = 0.5
    real, parameter :: MIN_EXCHANGE_C = 0.004

    character(len=kMAX_NAME_LENGTH) :: MMINLU
    logical :: FNDSOILW,FNDSNOWH,RDMAXALB
    integer :: num_soil_layers,ISURBAN,ISICE,ISWATER, ISLAKE
    integer :: exchange_term
    real*8  :: last_model_time
    real*8  :: Delta_t !! MJ added to detect the time for outputting 

    !Noah-MP specific
    real    :: NMP_SOILTSTEP
    integer :: IDVEG,IOPT_CRS,IOPT_BTR,IOPT_RUN,IOPT_SFC,IOPT_FRZ,IOPT_INF,IOPT_RAD,IOPT_ALB,IOPT_SNF,IOPT_TBOT, IOPT_TDRN, IOPT_NMPOUTPUT
    integer :: IOPT_STC, IOPT_GLA, IOPT_RSF, IOPT_SOIL, IOPT_PEDO, IOPT_CROP, IOPT_IRR, IOPT_IRRM, IZ0TLND, SF_URBAN_PHYSICS
    character(len=kMAX_NAME_LENGTH) :: landuse_name

    ! MJ added for FSM
    real,allocatable, dimension(:,:) :: current_snow, current_rain, current_precipitation, snowlsm, rainlsm

    ! Lake model: (allocated on lake init)
    real, allocatable, dimension (:,:)      ::      TH2 !, savedtke12d  lakedepth2d,
    integer :: lakeflag, lake_depth_flag, use_lakedepth, lake_count
    LOGICAL, allocatable, DIMENSION( :,: ) :: lake_or_not

contains


    subroutine lsm_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%landsurface == kLSM_NOAH .or. options%physics%landsurface == kLSM_BASIC) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%soil_water_content_liq,                    &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg, kVARS%albedo,                          &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%land_emissivity])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height, kVARS%QFX,  kVARS%land_emissivity, kVARS%soil_water_content_liq,            &  ! BK 2020/10/26
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0])!, kVARS%veg_type])    ! BK uncommented 2021/03/20
                         ! kVARS%soil_type, kVARS%land_mask, kVARS%vegetation_fraction]
        endif

        if (options%physics%landsurface == kLSM_NOAHMP) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%shortwave_direct, kVARS%shortwave_diffuse, kVARS%albedo,                                 &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%canopy_vapor_pressure, kVARS%canopy_temperature,                      &
                         kVARS%veg_leaf_temperature, kVARS%coeff_momentum_drag, kVARS%hpbl,                             &
                         kVARS%canopy_fwet, kVARS%snow_water_eq_prev, kVARS%water_table_depth, kVARS%water_aquifer,     &
                         kVARS%mass_leaf, kVARS%mass_root, kVARS%mass_stem, kVARS%mass_wood, kVARS%soil_carbon_fast,    &
                         kVARS%soil_carbon_stable, kVARS%eq_soil_moisture, kVARS%smc_watertable_deep, kVARS%recharge,   &
                         kVARS%recharge_deep, kVARS%storage_lake, kVARS%storage_gw, kVARS%mass_ag_grain,                &
                         kVARS%growing_degree_days, kVARS%plant_growth_stage, kVARS%temperature_2m_veg,                 &
                         kVARS%temperature_2m_bare, kVARS%mixing_ratio_2m_veg, kVARS%mixing_ratio_2m_bare,              &
                         kVARS%surface_rad_temperature, kVARS%net_ecosystem_exchange, kVARS%gross_primary_prod,         &
                         kVARS%net_primary_prod, kVARS%runoff_surface, kVARS%runoff_subsurface,                         &
                         kVARS%evap_canopy, kVARS%evap_soil_surface, kVARS%rad_absorbed_total, kVARS%rad_net_longwave,  &
                         kVARS%apar, kVARS%photosynthesis_total, kVARS%rad_absorbed_veg, kVARS%rad_absorbed_bare,       &
                         kVARS%stomatal_resist_total, kVARS%stomatal_resist_sun, kVARS%stomatal_resist_shade,           &
                         kVARS%lai, kVARS%sai, kVARS%snow_albedo_prev, kVARS%snow_age_factor, kVARS%canopy_water_ice,   &
                         kVARS%canopy_water_liquid, kVARS%vegetation_fraction_max, kVARS%crop_category,                 &
                         kVARS%date_planting, kVARS%date_harvest, kVARS%growing_season_gdd, kVARS%transpiration_rate,   &
                         kVARS%frac_within_gap, kVARS%frac_between_gap, kVARS%ground_temperature_canopy,                &
                         kVARS%ground_temperature_bare, kVARS%ch_veg, kVARS%ch_veg_2m, kVARS%ch_bare, kVARS%ch_bare_2m, &
                         kVARS%ch_under_canopy, kVARS%ch_leaf, kVARS%sensible_heat_veg, kVARS%sensible_heat_bare,       &
                         kVARS%sensible_heat_canopy, kVARS%evap_heat_veg, kVARS%evap_heat_bare, kVARS%evap_heat_canopy, &
                         kVARS%transpiration_heat, kVARS%ground_heat_veg, kVARS%ground_heat_bare, kVARS%snow_nlayers,   &
                         kVARS%net_longwave_veg, kVARS%net_longwave_bare, kVARS%net_longwave_canopy,                    &
                         kVARS%irr_frac_total, kVARS%irr_frac_sprinkler, kVARS%irr_frac_micro, kVARS%irr_frac_flood,    &
                         kVARS%irr_eventno_sprinkler, kVARS%irr_eventno_micro, kVARS%irr_eventno_flood,                 &
                         kVARS%irr_alloc_sprinkler, kVARS%irr_alloc_micro, kVARS%irr_alloc_flood, kVARS%irr_amt_flood,  &
                         kVARS%irr_evap_loss_sprinkler, kVARS%irr_amt_sprinkler, kVARS%irr_amt_micro,                   &
                         kVARS%evap_heat_sprinkler, kVARS%snowfall_ground, kVARS%rainfall_ground, kVARS%crop_type,      &
                         kVARS%ground_surf_temperature, kVARS%snow_temperature, kVARS%snow_layer_depth,                 &
                         kVARS%snow_layer_ice, kVARS%snow_layer_liquid_water, kVARS%soil_texture_1, kVARS%gecros_state, &
                         kVARS%soil_texture_2, kVARS%soil_texture_3, kVARS%soil_texture_4, kVARS%soil_sand_and_clay,    &
                         kVARS%vegetation_fraction_out, kVARS%latitude, kVARS%longitude, kVARS%cosine_zenith_angle,     &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%soil_water_content_liq, kVARS%xice,        &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%land_emissivity])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature,       &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,  kVARS%hpbl, kVARS%land_emissivity,  &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent, kVARS%QFX,                    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height, kVARS%canopy_water_ice, kVARS%canopy_vapor_pressure, kVARS%canopy_temperature,    &  ! BK 2020/10/26
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0,                                                  &
                         kVARS%ground_surf_temperature, kVARS%snow_nlayers, kVARS%veg_leaf_temperature,                 &
                         kVARS%canopy_water_liquid, kVARS%coeff_momentum_drag, kVARS%chs, kVARS%canopy_fwet,            &
                         kVARS%mass_leaf, kVARS%mass_root, kVARS%mass_stem, kVARS%mass_wood, kVARS%snow_water_eq_prev,  &
                         kVARS%snow_albedo_prev, kVARS%snow_temperature, kVARS%snow_layer_depth,  kVARS%snow_layer_ice, &
                         kVARS%snow_layer_liquid_water,kVARS%snowfall_ground, kVARS%rainfall_ground, kVARS%storage_lake,&
                         kVARS%storage_gw, kVARS%water_table_depth, kVARS%water_aquifer, kVARS%soil_carbon_fast,        &
                         kVARS%soil_carbon_stable, kVARS%lai, kVARS%sai, kVARS%soil_water_content_liq, kVARS%xice,      &
                         kVARS%mass_ag_grain, kVARS%growing_degree_days])!, kVARS%veg_type])    ! BK uncommented 2021/03/20
                         ! kVARS%soil_type, kVARS%land_mask, kVARS%vegetation_fraction]
        endif



        if (options%physics%snowmodel == kSM_FSM) then
            call options%alloc_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg, kVARS%slope_angle,                     &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity,                           &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%Nsnow,   &
                         kVARS%shd, kVARS%meltflux_out_tstep, kVARS%Sliq_out, &
                         kVARS%windspd_10m, kVARS%dSWE_salt, kVARS%dSWE_susp, kVARS%dSWE_subl, kVARS%dSWE_slide])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%exch_vars([kVARS%Ds, kVARS%Nsnow, kVARS%fsnow, kVARS%Sice, kVARS%Sliq, kVARS%snow_temperature])
             
             call options%restart_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,  kVARS%snowfall, kVARS%albedo, kVARS%QFX, kVARS%land_emissivity,            &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0,                                                  &
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%Nsnow  ])
        endif

       if (options%physics%watersurface > 0) then
            call options%alloc_vars( &
                         [kVARS%sst, kVARS%ustar, kVARS%surface_pressure, kVARS%water_vapor,            &
                         kVARS%temperature, kVARS%sensible_heat, kVARS%latent_heat, kVARS%land_mask,    &
                         kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2,                                            &
                         kVARS%humidity_2m, kVARS%temperature_2m, kVARS%skin_temperature, kVARS%u_10m, kVARS%v_10m])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%sst, kVARS%potential_temperature, kVARS%water_vapor, kVARS%skin_temperature,        &
                         kVARS%surface_pressure, kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,  &
                         kVARS%QFX,                                                                                 &
                         kVARS%humidity_2m, kVARS%temperature_2m])
        endif

        if (options%physics%watersurface == kWATER_LAKE ) then
            call options%alloc_vars( &
            [kVARS%lake_depth,kVARS%veg_type,kVARS%soil_type, kVARS%land_mask,kVARS%terrain,                        &
            kVARS%temperature,kVARS%pressure_interface, kVARS%dz_interface, kVARS%shortwave,  kVARS%longwave,       &
            kVARS%water_vapor, kVARS%latitude, kVARS%longitude, kVARS%sensible_heat, kVARS%latent_heat,             &
            kVARS%ground_heat_flux, kVARS%snow_water_equivalent, kVARS%t_lake3d, kVARS%dz_lake3d,                   &
            kVARS%t_soisno3d, kVARS%h2osoi_ice3d, kVARS%h2osoi_liq3d, kVARS%h2osoi_vol3d, kVARS%z3d,                &
            kVARS%dz3d, kVARS%watsat3d, kVARS%csol3d, kVARS%tkmg3d, kVARS%lakemask, kVARS%zi3d,                     &
            kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity, kVARS%xice,                        &
            kVARS%tksatu3d, kVARS%tkdry3d, kVARS%snl2d, kVARS%t_grnd2d,  kVARS%savedtke12d, kVARS%lakedepth2d,      & !  kVARS%snowdp2d, kVARS%h2osno2d,
            kVARS%lake_icefrac3d, kVARS%z_lake3d,kVARS%water_vapor, kVARS%potential_temperature     ])

            ! call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

            call options%restart_vars( &
            [kVARS%lake_depth,kVARS%veg_type,kVARS%soil_type, kVARS%land_mask,kVARS%terrain,                        &
            kVARS%temperature,kVARS%pressure_interface, kVARS%dz_interface, kVARS%shortwave,  kVARS%longwave,       &
            kVARS%water_vapor, kVARS%latitude, kVARS%longitude, kVARS%sensible_heat, kVARS%latent_heat,             &
            kVARS%ground_heat_flux, kVARS%snow_water_equivalent, kVARS%t_lake3d, kVARS%dz_lake3d,                   &
            kVARS%t_soisno3d, kVARS%h2osoi_ice3d, kVARS%h2osoi_liq3d, kVARS%h2osoi_vol3d, kVARS%z3d,                &
            kVARS%dz3d, kVARS%watsat3d, kVARS%csol3d, kVARS%tkmg3d, kVARS%lakemask, kVARS%zi3d,                     &
            kVARS%QFX, kVARS%land_emissivity, kVARS%xice,                                                           &
            kVARS%tksatu3d, kVARS%tkdry3d, kVARS%snl2d, kVARS%t_grnd2d, kVARS%savedtke12d, kVARS%lakedepth2d,       & !kVARS%snowdp2d, kVARS%h2osno2d,
            kVARS%lake_icefrac3d, kVARS%z_lake3d,kVARS%water_vapor, kVARS%potential_temperature ])
        endif



    end subroutine lsm_var_request

    subroutine calc_exchange_coefficient(wind,tskin,z0,airt,exchange_C)
        implicit none
        real, dimension(ims:ime,jms:jme),intent(in) :: wind, tskin, z0
        real, dimension(ims:ime,kms:kme,jms:jme),intent(in) :: airt
        real,dimension(ims:ime,jms:jme),intent(inout) :: exchange_C
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term, base_exchange_term
        
        lnz_atm_term = log((z_atm+z0)/z0)
        base_exchange_term=(75*karman**2 * sqrt((z_atm+z0)/z0)) / (lnz_atm_term**2)
        lnz_atm_term=(karman/lnz_atm_term)**2

        ! Richardson number
        exchange_C = 0

        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2+epsilon)
        ! Ri now is a function in atm_utlilities:
        ! calc_Richardson_nr(Ri, airt, tskin, z_atm, wind)

        ! "--------------------------------------------------"
        !  "Surface Richardson number"
        where(Ri<0)  exchange_C = lnz_atm_term * (1.0-(15.0*Ri)/(1.0+(base_exchange_term * sqrt((-1.0)*Ri))))
        where(Ri>=0) exchange_C = lnz_atm_term * 1.0/((1.0+15.0*Ri)*sqrt(1.0+5.0*Ri))

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_exchange_coefficient


! eqn A11 in Appendix A.2 of Chen et al 1997 (see below for reference)
    subroutine F2_formula(F2, z0, Ri)
        real, dimension(ims:ime,jms:jme), intent(inout) :: F2
        real, dimension(ims:ime,jms:jme), intent(in)    :: z0, Ri
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term
        
        lnz_atm_term = log((z_atm+z0)/z0)

        ! for the stable case from Mahrt (1987)
        where(Ri>=0) F2=exp(-Ri)
        ! for the unstable case from Holtslag and Beljaars (1989)
        where(Ri<0)  F2=(1-(15*Ri)/(1+((70.5*karman**2 * sqrt(-Ri * z_atm/z0))/(lnz_atm_term**2))) )

    end subroutine F2_formula
        
    
!From Appendix A.2 in Chen et al 1997
! Impact of Atmospheric Surface-layer Parameterizations in the new Land-surface Scheme of the Ncep Mesoscale ETA Model
! Boundary-Layer Meteorology 85:391-421
    subroutine calc_mahrt_holtslag_exchange_coefficient(wind,tskin,airt,znt,exchange_C)
        implicit none
        real, dimension(:,:),intent(in) :: wind,tskin,znt
        real, dimension(:,:,:),intent(in) :: airt
        real,dimension(:,:),intent(inout) :: exchange_C
        
        real, dimension(ims:ime,jms:jme)  :: lnz_atm_term, base_exchange_term
        
        lnz_atm_term = log((z_atm+znt)/znt)
        
        ! Richardson number
        Ri = gravity/airt(:,1,:) * (airt(:,1,:)-tskin)*z_atm/(wind**2+epsilon)

        call F2_formula(base_exchange_term,znt,Ri)

        exchange_C = karman**2 * base_exchange_term / lnz_atm_term**2

        where(exchange_C > MAX_EXCHANGE_C) exchange_C=MAX_EXCHANGE_C
        where(exchange_C < MIN_EXCHANGE_C) exchange_C=MIN_EXCHANGE_C
    end subroutine calc_mahrt_holtslag_exchange_coefficient


    subroutine surface_diagnostics(HFX, QFX, TSK, QSFC, CHS2, CQS2,T2, Q2, PSFC, &
                                    VEGFRAC, veg_type, land_mask, T2veg, T2bare, Q2veg, Q2bare)
        ! taken almost directly / exactly from WRF's module_sf_sfcdiags.F
        !-- HFX           net upward heat flux at the surface (W/m^2)
        !-- QFX           net upward moisture flux at the surface (kg/m^2/s)
        !-- TSK           surface temperature (K)
        !-- qsfc          specific humidity at lower boundary (kg/kg)
        implicit none
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  T2veg, T2bare, Q2veg, Q2bare
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  VEGFRAC
        INTEGER, DIMENSION(ims:ime, jms:jme ), INTENT(IN) ::  land_mask, veg_type
        integer :: i,j
        real :: rho

        ! !$omp parallel default(shared), private(i,j,rho)
        ! !$omp do
        do j=jts,jte
            do i=its,ite

                ! if ((domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(i,j)/=13).and.(domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(i,j)/=15).and.(domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(i,j)/=16).and.(domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(i,j)/=21)) then
                ! over glacier, urban and barren, noahmp veg 2m T is 0 or -9999e35
                if ((T2veg(i,j) > 200).and.(land_mask(i,j)==kLC_LAND)) then
                    T2(i,j) = VEGFRAC(i,j) * T2veg(i,j) &
                        + (1-VEGFRAC(i,j)) * T2bare(i,j)
                    Q2(i,j) = VEGFRAC(i,j) * Q2veg(i,j) &
                        + (1-VEGFRAC(i,j)) * Q2bare(i,j)
                else
                    ! over glacier we don't want to use the bare ground temperature though
                    if ((veg_type(i,j)/=ISICE)           &  ! was /=15  (15=snow/ice in MODIFIED_IGBP_MODIS_NOAH)
                        .and.(veg_type(i,j)/=ISLAKE)     &  ! was /=21  (ISLAKE  = options%lsm%lake_category)   # 17 is water, 21 is lakes (in MODIFIED_IGBP_MODIS_NOAH ) MODIFY FOR GENERIC LU TYPES!!
                        .and.(land_mask(i,j)==kLC_LAND)) then
                        T2(i,j) = T2bare(i,j)
                        Q2(i,j) = Q2bare(i,j)
                    else
                        RHO = PSFC(I,J)/(R_d * TSK(I,J))

                        if(CQS2(I,J).lt.1.E-3) then
                           Q2(I,J) = QSFC(I,J)
                        else
                           Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
                        endif

                        if(CHS2(I,J).lt.1.E-3) then
                           T2(I,J) = TSK(I,J)
                        else
                           T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP*CHS2(I,J))
                        endif
                    endif
                endif
                if (Q2(i,j) < SMALL_QV) Q2(i,j) = SMALL_QV

                ! TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
            enddo
        enddo
        ! !$omp end do
        ! !$omp end parallel
    end subroutine surface_diagnostics
    

    subroutine surface_diagnostics_FSM(HFX, QFX, TSK, QSFC, CHS2, CQS2,T2, Q2, PSFC)
        ! taken almost directly / exactly from WRF's module_sf_sfcdiags.F
        !-- HFX           net upward heat flux at the surface (W/m^2)
        !-- QFX           net upward moisture flux at the surface (kg/m^2/s)
        !-- TSK           surface temperature (K)
        !-- qsfc          specific humidity at lower boundary (kg/kg)

        implicit none
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  HFX, QFX, TSK, QSFC
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(INOUT) ::  Q2, T2
        REAL, DIMENSION(ims:ime, jms:jme ), INTENT(IN)    ::  PSFC, CHS2, CQS2
        integer :: i,j
        real :: rho
        real :: CP_=1005.

        !$omp parallel default(shared), private(i,j,rho)
        !$omp do
        do j=jts,jte
            do i=its,ite
                RHO = PSFC(I,J)/(R_d * TSK(I,J))
                !if(CQS2(I,J).lt.1.E-3) then
                !   Q2(I,J) = QSFC(I,J)
                !else
                   Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
                !endif
                !if(CHS2(I,J).lt.1.E-3) then
                !   T2(I,J) = TSK(I,J)
                !else
                   T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP_*CHS2(I,J))
                !endif
                if (Q2(i,j) < SMALL_QV) Q2(i,j) = SMALL_QV
                ! TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
                 !if ( isnan(HFX(I,J)) .or. abs(TSK(I,J)-T2(I,J))>30 ) write(*,*) "img-H222",i,j,this_image(), HFX(I,J), TSK(I,J), T2(I,J), CHS2(I,J) 
                 !if (STD_OUT_PE) write(*,*),"img-H222",i,j,this_image(), RHO, TSK(I,J), T2(I,J),PSFC(I,J)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
         !write(*,*),"RHO*CP",this_image(), RHO,CP,rd
    end subroutine surface_diagnostics_FSM
        

    subroutine lsm_apply_fluxes(domain,options,dt)
        ! add sensible and latent heat fluxes to the first atm level
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t), intent(in)   :: options
        real, intent(in) :: dt
        integer :: i,j,k
        integer, SAVE :: nz = 0
        real :: dTemp, lhdQV

        ! PBL scheme should handle the distribution of sensible and latent heat fluxes. If we are
        ! running the LSM without a PBL scheme, as may be done for High-resolution runs, then 
        ! run apply fluxes to still apply heat fluxes calculated by LSM
        if ( (options%physics%landsurface>0 .or. options%physics%watersurface>0 ) .and. (options%physics%boundarylayer==0)) then
            associate(density       => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d,             &
                    sensible_heat => domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,           &
                    latent_heat   => domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d,             &
                    dz            => domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,        &
                    pii           => domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,               &
                    th            => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d, &
                    qv            => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d          &
                )

            do j = jts, jte
            do k = kts, kts + nz
            do i = its, ite

                ! convert sensible heat flux to a temperature delta term
                ! (J/(s*m^2) * s / (J/(kg*K)) => kg*K/m^2) ... /((kg/m^3) * m) => K
                dTemp = (sensible_heat(i,j) * dt/cp)  &
                        / (density(i,k,j))
                ! add temperature delta converted back to potential temperature
                th(i,k,j) = th(i,k,j) + (dTemp / pii(i,k,j))

                ! convert latent heat flux to a mixing ratio tendancy term
                ! (J/(s*m^2) * s / (J/kg) => kg/m^2) ... / (kg/m^3 * m) => kg/kg
                lhdQV = (latent_heat(i,j) / XLV * dt) &
                        / (density(i,k,j))
                ! add water vapor in kg/kg
                qv(i,k,j) = qv(i,k,j) + lhdQV

            end do ! i
            end do ! k
            end do ! j

            ! write(*,*) MINVAL(lhdQV), MAXVAL(lhdQV), 'kg/kg (min/max) added to QV at', domain%sim_time%hour

            ! enforce some minimum water vapor content... just in case
            where(qv < SMALL_QV) qv = SMALL_QV
            end associate
        endif
    end subroutine lsm_apply_fluxes

    subroutine allocate_noah_data(num_soil_layers)
        implicit none
        integer, intent(in) :: num_soil_layers
        integer :: i

        ITIMESTEP=1
        
        if (allocated(SMSTAV)) deallocate(SMSTAV)
        if (allocated(SFCRUNOFF)) deallocate(SFCRUNOFF)
        if (allocated(UDRUNOFF)) deallocate(UDRUNOFF)
        if (allocated(SNOWC)) deallocate(SNOWC)
        if (allocated(ACSNOW)) deallocate(ACSNOW)
        if (allocated(ACSNOM)) deallocate(ACSNOM)
        if (allocated(SNOALB)) deallocate(SNOALB)
        if (allocated(QGH)) deallocate(QGH)
        if (allocated(GSW)) deallocate(GSW)
        if (allocated(chklowq)) deallocate(chklowq) 
        if (allocated(ALBEDO)) deallocate(ALBEDO)
        if (allocated(ALBBCK)) deallocate(ALBBCK)
        if (allocated(EMBCK)) deallocate(EMBCK)
        if (allocated(CPM)) deallocate(CPM)
        if (allocated(SR)) deallocate(SR)
        if (allocated(SHDMIN)) deallocate(SHDMIN)
        if (allocated(SHDMAX)) deallocate(SHDMAX)
        if (allocated(SNOTIME)) deallocate(SNOTIME)
        if (allocated(SNOPCX)) deallocate(SNOPCX)
        if (allocated(POTEVP)) deallocate(POTEVP)
        if (allocated(SMCREL)) deallocate(SMCREL)
        if (allocated(NOAHRES)) deallocate(NOAHRES)
        if (allocated(fbur_2d)) deallocate(fbur_2d)
        if (allocated(fgsn_2d)) deallocate(fgsn_2d)
        if (allocated(flx4_2d)) deallocate(flx4_2d)
        if (allocated(fvb_2d)) deallocate(fvb_2d)
        if (allocated(VEGFRAC)) deallocate(VEGFRAC)
        if (allocated(nmp_snow)) deallocate(nmp_snow)
        if (allocated(nmp_snowh)) deallocate(nmp_snowh)

        allocate(SMSTAV(ims:ime,jms:jme),source=0.5)!average soil moisture available for transp (between SMCWLT and SMCMAX)
        allocate(SFCRUNOFF(ims:ime,jms:jme), source=0.0)
        allocate(UDRUNOFF(ims:ime,jms:jme), source=0.0)
        allocate(SNOWC(ims:ime,jms:jme), source=0.0)
        allocate(ACSNOW(ims:ime,jms:jme), source=0.0)
        allocate(ACSNOM(ims:ime,jms:jme), source=0.0)
        allocate(SNOALB(ims:ime,jms:jme), source=0.8)

        allocate(QGH(ims:ime,jms:jme), source=0.02) ! saturated mixing ratio at ~20C
        allocate(GSW(ims:ime,jms:jme), source=0.0)
        allocate(chklowq(ims:ime,jms:jme), source=0.0)
        allocate(ALBEDO(ims:ime,jms:jme), source=0.17)
        allocate(ALBBCK(ims:ime,jms:jme), source=0.17)
        allocate(EMBCK(ims:ime,jms:jme), source=0.99)
        allocate(CPM(ims:ime,jms:jme), source=0.0)
        allocate(SR(ims:ime,jms:jme), source=0.0)

        allocate(SHDMIN(ims:ime,jms:jme), source=0.0)
        allocate(SHDMAX(ims:ime,jms:jme), source=100.0)
        allocate(SNOTIME(ims:ime,jms:jme), source=0.0)
        allocate(SNOPCX(ims:ime,jms:jme), source=0.0)
        allocate(POTEVP(ims:ime,jms:jme), source=0.0)
        allocate(SMCREL(ims:ime,num_soil_layers,jms:jme), source=0.0)
        allocate(NOAHRES(ims:ime,jms:jme), source=0.0)
        allocate(VEGFRAC(ims:ime,jms:jme), source=50.0)

        allocate(fbur_2d(ims:ime,jms:jme), source=0.0)
        allocate(fgsn_2d(ims:ime,jms:jme), source=0.0)
        allocate(flx4_2d(ims:ime,jms:jme), source=0.0)
        allocate(fvb_2d(ims:ime,jms:jme), source=0.0)

        allocate(nmp_snow(ims:ime,jms:jme))
        allocate(nmp_snowh(ims:ime,jms:jme))

        XICE_THRESHOLD = 1
        RDLAI2D = .false. !TLE check back on this one
        USEMONALB = .false.
        MYJ = .false.
        FRPCPN = .false. ! set this to true and calculate snow ratio to use microphysics based snow/rain partitioning
        ua_phys = .false.

        if (.not.(allocated(Zs))) allocate(Zs(num_soil_layers))
        if (.not.(allocated(DZs))) allocate(DZs(num_soil_layers))
        !DZs = [0.1,0.3,0.6,1.0]
        !if (options%physics%landsurface==kSM_FSM) then !! MJ added to adapt with FSM for the soil layer thickness
            DZs = [0.1,0.2,0.4,0.8]
        !endif
        Zs(1) = DZs(1)/2
        do i = 2,num_soil_layers
            Zs(i) = Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
        end do

    end subroutine allocate_noah_data

    subroutine lsm_init(domain,options, context_chng)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: context_chng
        integer :: i,j
        logical :: context_change
        real,allocatable,dimension(:,:) :: chstarxy

        if (options%physics%landsurface == 0) return   !! So we cannot (currently) run without lsm but with water.

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing LSM"

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "    max soil_deep_temperature on init: ", maxval(domain%vars_2d(domain%var_indx(kVARS%soil_deep_temperature)%v)%data_2d)
        if (STD_OUT_PE .and. .not.context_change) write(*,*) "    max skin_temperature on init: ", maxval(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d)

        exchange_term = 1

        ! module level variables for easy access... need to think about tiling to permit halo processing separately.
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        if (allocated(SW)) deallocate(SW)
        if (allocated(Z0)) deallocate(Z0)
        if (allocated(QZ0)) deallocate(QZ0)
        if (allocated(QSFC)) deallocate(QSFC)
        if (allocated(Ri)) deallocate(Ri)
        if (allocated(z_atm)) deallocate(z_atm)
        if (allocated(current_precipitation)) deallocate(current_precipitation)
        if (allocated(windspd)) deallocate(windspd)

        allocate(SW(ims:ime,jms:jme))
        allocate(Z0(ims:ime,jms:jme), source=domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d)
        allocate(QZ0(ims:ime,jms:jme), source=0.0)
        allocate(QSFC(ims:ime,jms:jme), source=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(:,kms,:))
        allocate(Ri(ims:ime,jms:jme), source=0.0)
        allocate(z_atm(ims:ime,jms:jme), source=50.0)
        allocate(current_precipitation(ims:ime,jms:jme), source=0.0)
        allocate(windspd(ims:ime,jms:jme), source=1.0)

        ! NOTE, these fields have probably not been initialized yet...
        ! windspd = sqrt(domain%u10**2+domain%v10**2)

        if (options%physics%landsurface > kLSM_BASIC) then
            if (options%physics%microphysics == 0) then
                write(*,*) "Land Surface models need microphysics to run"
                stop "Land Surface models need microphysics to run"
            endif

            ! MJ added:
            if (allocated(current_snow)) deallocate(current_snow)
            if (allocated(current_rain)) deallocate(current_rain)
            allocate(current_snow(ims:ime,jms:jme), source=0.0) ! MJ added 
            allocate(current_rain(ims:ime,jms:jme), source=0.0) ! MJ added 
            if (allocated(snowlsm)) deallocate(snowlsm)
            if (allocated(rainlsm)) deallocate(rainlsm)
            allocate(snowlsm(ims:ime,jms:jme), source=0.0) ! MJ added 
            allocate(rainlsm(ims:ime,jms:jme), source=0.0) ! MJ added 

        endif
        
        if (options%physics%landsurface==kLSM_NOAH .or. options%physics%landsurface==kLSM_NOAHMP .or. options%physics%snowmodel==kSM_FSM) then
            num_soil_layers=options%lsm%num_soil_layers ! to .nml?
            call allocate_noah_data(num_soil_layers)
        endif
        ! initial guesses (not needed?)
        domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d = domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d(:,kms,:)
        domain%vars_2d(domain%var_indx(kVARS%humidity_2m)%v)%data_2d = domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(:,kms,:)
        
        if(options%domain%snowh_var /="" .or. options%domain%swe_var /="") then 
            FNDSNOWH= .True.
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Find snow height in file i.s.o. calculating them from SWE: FNDSNOWH=", FNDSNOWH
        else
            FNDSNOWH=.False. ! calculate SNOWH from SNOW
        endif

        ! Noah Land Surface Model
        if (options%physics%landsurface==kLSM_NOAH) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Noah LSM"

            !These are still needed for NoahLSM, since it assumes that the exchange coefficients are multiplied by windspeed
            if (allocated(CHS)) deallocate(CHS)
            if (allocated(CHS2)) deallocate(CHS2)
            if (allocated(CQS2)) deallocate(CQS2)
            allocate(CHS(ims:ime,jms:jme))
            allocate(CHS2(ims:ime,jms:jme))
            allocate(CQS2(ims:ime,jms:jme))

            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm%urban_category
            ISICE   = options%lsm%ice_category
            ISWATER = options%lsm%water_category
            MMINLU  = options%lsm%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm%lake_category


            if (options%lsm%monthly_albedo) then
                if (.not.options%lsm%monthly_vegfrac) Then
                    print*, "ERROR, monthly albedo requires monthly vegfrac"
                    error stop
                endif
                ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :)
            else
                ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, 1, :)
            endif
            if (options%lsm%monthly_vegfrac) then
                VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, domain%sim_time%month, :)
            else
                VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, 1, :)
            endif
            cur_vegmonth = domain%sim_time%month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d
            ! prevents init from failing when processing water points that may have "soil_t"=0
            where(domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d<200) domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d=200
            where(domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d<0.0001) domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d=0.0001

            call LSM_NOAH_INIT( VEGFRAC,                             &
                                domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,& !SNOW, &  BK 18/03/2021
                                SNOWC,                               &
                                domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,          & !SNOWH, &   BK 18/03/2021
                                domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d,         &
                                domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,     & !-- SMSTAV      Soil moisture availability for evapotranspiration ( fraction between SMCWLT and SMCMXA)
                                domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d,   &
                                SFCRUNOFF,                           &
                                UDRUNOFF,                            &
                                ACSNOW,                              &
                                ACSNOM,                              &
                                domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di,                     &
                                domain%vars_2d(domain%var_indx(kVARS%soil_type)%v)%data_2di,                    &
                                domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,     &
                                domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d,   &
                                domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d,&
                                ZS,                                  &
                                DZS,                                 &
                                MMINLU,                              &
                                SNOALB,                              &
                                FNDSOILW,                            &
                                FNDSNOWH,                            &
                                RDMAXALB,                            &
                                num_soil_layers,                     &
                                .False.,                             & ! nlayers, is_restart (can't yet)
                                .True. ,                             & ! allowed_to_read (e.g. soilparm.tbl)
                                ids,ide, jds,jde, kds,kde,           &
                                ims,ime, jms,jme, kms,kme,           &
                                its,ite, jts,jte, kts,kte  )

            domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d
            domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d=0.01
            ! where(domain%veg_type==ISWATER) domain%land_mask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent
            where((domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di==ISWATER) .OR. (domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di==ISLAKE)) domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di=kLC_WATER  ! include lakes.
        endif

        ! Noah-MP Land Surface Model
        if (options%physics%landsurface==kLSM_NOAHMP) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Noah-MP LSM"

            FNDSOILW=.False. ! calculate SOILW (this parameter is ignored in LSM_NOAH_INIT)
            RDMAXALB=.False.

            ISURBAN = options%lsm%urban_category
            ISICE   = options%lsm%ice_category
            ISWATER = options%lsm%water_category
            MMINLU  = options%lsm%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm%lake_category

            if (options%lsm%monthly_albedo) then
                if (.not.options%lsm%monthly_vegfrac) Then
                    print*, "ERROR, monthly albedo requires monthly vegfrac"
                    error stop
                endif
                ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :)
            else
                ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, 1, :)
            endif
            if (options%lsm%monthly_vegfrac) then
                VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, domain%sim_time%month, :)
            else
                VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, 1, :)
            endif
            cur_vegmonth = domain%sim_time%month

            ! save the canopy water in a temporary variable in case this is a restart run because lsm_init resets it to 0
            domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d
            ! prevents init from failing when processing water points that may have "soil_t"=0
            where(domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d<200) domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d=200
            where(domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d<0.0001) domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d=0.0001

            IDVEG = options%lsm%nmp_dveg
            IOPT_CRS = options%lsm%nmp_opt_crs
            IOPT_SFC = options%lsm%nmp_opt_sfc
            IOPT_BTR = options%lsm%nmp_opt_btr
            IOPT_RUN = options%lsm%nmp_opt_run
            IOPT_INF = options%lsm%nmp_opt_inf
            IOPT_FRZ = options%lsm%nmp_opt_frz
            IOPT_INF = options%lsm%nmp_opt_inf
            IOPT_RAD = options%lsm%nmp_opt_rad
            IOPT_ALB = options%lsm%nmp_opt_alb
            IOPT_SNF = options%lsm%nmp_opt_snf
            IOPT_TBOT = options%lsm%nmp_opt_tbot
            IOPT_STC = options%lsm%nmp_opt_stc
            IOPT_GLA = options%lsm%nmp_opt_gla
            IOPT_RSF = options%lsm%nmp_opt_rsf
            IOPT_SOIL = options%lsm%nmp_opt_soil
            IOPT_PEDO = options%lsm%nmp_opt_pedo
            IOPT_CROP = options%lsm%nmp_opt_crop
            IOPT_IRR = options%lsm%nmp_opt_irr
            IOPT_IRRM = options%lsm%nmp_opt_irrm
            IOPT_TDRN = options%lsm%nmp_opt_tdrn
            IOPT_NMPOUTPUT = options%lsm%noahmp_output
            IZ0TLND = options%sfc%iz0tlnd
            NMP_SOILTSTEP = options%lsm%nmp_soiltstep
            SF_URBAN_PHYSICS = options%lsm%sf_urban_phys


            if (options%physics%snowmodel==kSM_FSM) IOPT_SNF = 4 !This will turn off precipitation partitioning in NoahMP, letting us remove snow from NoahMP


            !allocate dummy variable that doesn't do anything
            allocate(chstarxy(ims:ime,jms:jme), source=0.0)
            call NOAHMP_INIT ( MMINLU,                                  &
                                domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,   &
                                domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,             &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d,            &
                                domain%vars_2d(domain%var_indx(kVARS%soil_type)%v)%data_2di,                       &
                                domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di,                        &
                                domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,                &
                                domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,        &
                                domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d,      &
                                domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d , &
                                DZS , FNDSOILW , FNDSNOWH ,             &
                                domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,        &
                                domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di,           &
                                domain%vars_2d(domain%var_indx(kVARS%veg_leaf_temperature)%v)%data_2d,    &
                                domain%vars_2d(domain%var_indx(kVARS%ground_surf_temperature)%v)%data_2d, &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_water_ice)%v)%data_2d,        &
                                domain%vars_2d(domain%var_indx(kVARS%soil_deep_temperature)%v)%data_2d,   &
                                domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,                            &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_water_liquid)%v)%data_2d,     &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_vapor_pressure)%v)%data_2d,   &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_temperature)%v)%data_2d,      &
                                domain%vars_2d(domain%var_indx(kVARS%coeff_momentum_drag)%v)%data_2d,     &
                                domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d,                     &
                                domain%vars_2d(domain%var_indx(kVARS%canopy_fwet)%v)%data_2d,             &
                                domain%vars_2d(domain%var_indx(kVARS%snow_water_eq_prev)%v)%data_2d,      &
                                domain%vars_2d(domain%var_indx(kVARS%snow_albedo_prev)%v)%data_2d,        &
                                domain%vars_2d(domain%var_indx(kVARS%snowfall_ground)%v)%data_2d,         &
                                domain%vars_2d(domain%var_indx(kVARS%rainfall_ground)%v)%data_2d,         &
                                domain%vars_2d(domain%var_indx(kVARS%storage_lake)%v)%data_2d,            &
                                domain%vars_2d(domain%var_indx(kVARS%water_table_depth)%v)%data_2d,       &
                                domain%vars_2d(domain%var_indx(kVARS%water_aquifer)%v)%data_2d,           &
                                domain%vars_2d(domain%var_indx(kVARS%storage_gw)%v)%data_2d,              &
                                domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(:,1:3,:),        &
                                domain%vars_3d(domain%var_indx(kVARS%snow_layer_depth)%v)%data_3d(:,1:7,:),        &
                                domain%vars_3d(domain%var_indx(kVARS%snow_layer_ice)%v)%data_3d(:,1:3,:),          &
                                domain%vars_3d(domain%var_indx(kVARS%snow_layer_liquid_water)%v)%data_3d(:,1:3,:), &
                                domain%vars_2d(domain%var_indx(kVARS%mass_leaf)%v)%data_2d,               &
                                domain%vars_2d(domain%var_indx(kVARS%mass_root)%v)%data_2d,               &
                                domain%vars_2d(domain%var_indx(kVARS%mass_stem)%v)%data_2d,               &
                                domain%vars_2d(domain%var_indx(kVARS%mass_wood)%v)%data_2d,               &
                                domain%vars_2d(domain%var_indx(kVARS%soil_carbon_stable)%v)%data_2d,      &
                                domain%vars_2d(domain%var_indx(kVARS%soil_carbon_fast)%v)%data_2d,        &
                                domain%vars_2d(domain%var_indx(kVARS%lai)%v)%data_2d,                     &
                                domain%vars_2d(domain%var_indx(kVARS%sai)%v)%data_2d,                     &
                                domain%vars_2d(domain%var_indx(kVARS%mass_ag_grain)%v)%data_2d,           &
                                domain%vars_2d(domain%var_indx(kVARS%growing_degree_days)%v)%data_2d,     &
                                domain%vars_3d(domain%var_indx(kVARS%crop_type)%v)%data_3d,               &
                                domain%vars_2d(domain%var_indx(kVARS%crop_category)%v)%data_2di,                   &
                                domain%vars_2d(domain%var_indx(kVARS%irr_eventno_sprinkler)%v)%data_2di,           &
                                domain%vars_2d(domain%var_indx(kVARS%irr_eventno_micro)%v)%data_2di,               &
                                domain%vars_2d(domain%var_indx(kVARS%irr_eventno_flood)%v)%data_2di,               &
                                domain%vars_2d(domain%var_indx(kVARS%irr_alloc_sprinkler)%v)%data_2d,     &
                                domain%vars_2d(domain%var_indx(kVARS%irr_alloc_micro)%v)%data_2d,         &
                                domain%vars_2d(domain%var_indx(kVARS%irr_alloc_flood)%v)%data_2d,         &
                                domain%vars_2d(domain%var_indx(kVARS%irr_evap_loss_sprinkler)%v)%data_2d, &
                                domain%vars_2d(domain%var_indx(kVARS%irr_amt_sprinkler)%v)%data_2d,       &
                                domain%vars_2d(domain%var_indx(kVARS%irr_amt_micro)%v)%data_2d,           &
                                domain%vars_2d(domain%var_indx(kVARS%irr_amt_flood)%v)%data_2d,           &
                                domain%vars_2d(domain%var_indx(kVARS%evap_heat_sprinkler)%v)%data_2d,     &
                                domain%vars_2d(domain%var_indx(kVARS%temperature_2m_veg)%v)%data_2d,      &
                                domain%vars_2d(domain%var_indx(kVARS%temperature_2m_bare)%v)%data_2d,     &
                                chstarxy,                               &   !doesn't do anything -_-
                                num_soil_layers,                        &
                                options%restart%restart,             &    !restart
                                .True.,                                 &    !allowed_to_read
                                IOPT_RUN,  IOPT_CROP, IOPT_IRR, IOPT_IRRM, &
                                SF_URBAN_PHYSICS,                         &  ! urban scheme
                                ids,ide, jds,jde, kds,kde,                &
                                ims,ime, jms,jme, kms,kme,                &
                                its,ite, jts,jte, kts,kte)

  !                           TLE: GROUNDWATER OFF FOR NOW
  !                                   smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&     ! Optional groundwater
  !                                   wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &      ! Optional groundwater
  !                                   fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &      ! Optional groundwater
  !                                   rechclim,                                                             &      ! Optional groundwater
  !                                   gecros_state)                                                                ! Optional gecros crop

            domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d
            domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d=0.01
            ! where(domain%veg_type==ISWATER) domain%land_mask=kLC_WATER ! ensure VEGTYPE (land cover) and land-sea mask are consistent (BK 202208: this does not include lakes!!)
            where((domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di==ISWATER) .OR. (domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di==ISLAKE)) domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di=kLC_WATER
        endif

        if(options%physics%watersurface==kWATER_LAKE) then
        ! ____________ Lake model ______________________
        ! From WRF's /run/README.namelist:  These could at some point become namelist options in ICAR?
        ! lakedepth_default(max_dom)          = 50,      ! default lake depth (If there is no lake_depth information in the input data, then lake depth is assumed to be 50m)
        ! lake_min_elev(max_dom)              = 5,       ! minimum elevation of lakes. May be used to determine whether a water point is a lake in the absence of lake
        !                                                  category. If the landuse type includes 'lake' (i.e. Modis_lake and USGS_LAKE), this variable is of no effects.
        ! use_lakedepth (max_dom)             = 1,       ! option to use lake depth data. Lake depth data is available from 3.6 geogrid program. If one didn't process
        !                                                    the lake depth data, but this switch is set to 1, the program will stop and tell one to go back to geogrid
        !                                                     program.
        !                                     = 0, do not use lake depth data.

            if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Lake model"

            ! allocate arrays:
            if (allocated(lake_or_not)) deallocate( lake_or_not)
            if (allocated(TH2)) deallocate( TH2)
            allocate(lake_or_not(ims:ime, jms:jme))
            allocate(TH2( ims:ime, jms:jme ))

            ! ISURBAN = options%lsm%urban_category
            ISICE   = options%lsm%ice_category
            ISWATER = options%lsm%water_category
            ! MMINLU  = options%lsm%LU_Categories !"MODIFIED_IGBP_MODIS_NOAH"
            ISLAKE  = options%lsm%lake_category

            ! allocate_noah_data already sets xice_threshold, so if we are using noah (mp/lsm) leave as is.
            if(.not.(options%physics%landsurface==kLSM_NOAHMP .OR. options%physics%landsurface==kLSM_NOAH)) then
                xice_threshold = 1.0  ! allocate_noah_data sets it to 1., BUT WRF's module_physics_init.F sets xice_threshold to 0.5 .... so?
            endif

            lake_count=0
            if(ISLAKE==-1) then
                if(STD_OUT_PE .and. .not.context_change) write(*,*)  "   WARNING: no lake category in LU data: The model will try to guess lake-gridpoints. This option has not been properly tested!"
                lakeflag=0  ! If no lake cat is provided, the lake model will determine lakes based
                            ! on the criterion (ivgtyp(i,j)==iswater .and. ht(i,j)>=lake_min_elev))
            else
                lakeflag=1
                ! from WRF's module_initialize_real.F:
                DO j = jts, MIN(jde-1,jte)
                    DO i = its, MIN(ide-1,ite)
                    !    IF ( grid%lu_index(i,j) .NE. grid%islake ) THEN
                        if(domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(i,j) .NE. ISLAKE ) then
                            domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d(i,j) = 0       ! grid%lakemask(i,j) = 0
                        ELSE
                            domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d(i,j) = 1       ! grid%lakemask(i,j) = 1
                            lake_count= lake_count + 1
                        end if
                    END DO
                END DO
            endif
            ! if(options%general%debug) write(*,*)"   ",lake_count, " lake cells in image ", this_image()

            ! setlake_depth_flag and use_lakedepth flag. (They seem to be redundant, but whatever):
            if(domain%var_indx(kVARS%lake_depth)%v > 0) then
                if(STD_OUT_PE .and. .not.context_change) write(*,*) "   Using Lake depth data "
                use_lakedepth = 1
                lake_depth_flag = 1
            else
                use_lakedepth = 0
                lake_depth_flag = 0
            endif

            call lakeini( &
                IVGTYP = domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di                        &
                ,ISLTYP = domain%vars_2d(domain%var_indx(kVARS%soil_type)%v)%data_2di                      &
                ,HT=domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d                      & ! terrain height [m] if ht(i,j)>=lake_min_elev -> lake  (grid%ht in WRF)
                ,SNOW=domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d      & !i  ! SNOW in kg/m^2  (NoahLSM: SNOW liquid water-equivalent snow depth (m)
                ,lake_min_elev=5.                               & ! minimum elevation of lakes. May be used to determine whether a water point is a lake in the absence of lake category. If the landuse type includes 'lake' (i.e. Modis_lake and USGS_LAKE), this variable is of no effects.
                ,restart=options%restart%restart                & ! if restart, this (lakeini) subroutine is simply skipped.
                ,lakedepth_default=50.                          & ! default lake depth (If there is no lake_depth information in the input data, then lake depth is assumed to be 50m)
                ,lake_depth=domain%vars_2d(domain%var_indx(kVARS%lake_depth)%v)%data_2d           & !INTENT(IN)
                ,lakedepth2d=domain%vars_2d(domain%var_indx(kVARS%lakedepth2d)%v)%data_2d         & !INTENT(OUT) (will be equal to lake_depth if lake_depth data is provided in hi-res input, otherwise lakedepth_default)
                ,savedtke12d=domain%vars_2d(domain%var_indx(kVARS%savedtke12d)%v)%data_2d         & !INTENT(OUT)
                ,snowdp2d=domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d            & ! domain%snowdp2d%data_2d
                ,h2osno2d=domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d  & !domain%h2osno2d%data_2d
                ,snl2d=domain%vars_2d(domain%var_indx(kVARS%snl2d)%v)%data_2d                     & ! snowlevel 2d?
                ,t_grnd2d=domain%vars_2d(domain%var_indx(kVARS%t_grnd2d)%v)%data_2d               & ! ground temperature?
                ,t_lake3d=domain%vars_3d(domain%var_indx(kVARS%t_lake3d)%v)%data_3d               & ! lake temperature 3d
                ,lake_icefrac3d=domain%vars_3d(domain%var_indx(kVARS%lake_icefrac3d)%v)%data_3d   & ! lake ice fraction ?
                ,z_lake3d=domain%vars_3d(domain%var_indx(kVARS%z_lake3d)%v)%data_3d               & !
                ,dz_lake3d=domain%vars_3d(domain%var_indx(kVARS%dz_lake3d)%v)%data_3d             &
                ,t_soisno3d=domain%vars_3d(domain%var_indx(kVARS%t_soisno3d)%v)%data_3d           & ! temperature of both soil and snow
                ,h2osoi_ice3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_ice3d)%v)%data_3d       & !  ice lens (kg/m2)
                ,h2osoi_liq3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_liq3d)%v)%data_3d       & ! liquid water (kg/m2)
                ,h2osoi_vol3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_vol3d)%v)%data_3d       & ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
                ,z3d=domain%vars_3d(domain%var_indx(kVARS%z3d)%v)%data_3d                         & ! layer depth for snow & soil (m)
                ,dz3d=domain%vars_3d(domain%var_indx(kVARS%dz3d)%v)%data_3d                       & ! layer thickness for soil or snow (m)
                ,zi3d=domain%vars_3d(domain%var_indx(kVARS%zi3d)%v)%data_3d                                              &
                ,watsat3d=domain%vars_3d(domain%var_indx(kVARS%watsat3d)%v)%data_3d               &
                ,csol3d=domain%vars_3d(domain%var_indx(kVARS%csol3d)%v)%data_3d                   &
                ,tkmg3d=domain%vars_3d(domain%var_indx(kVARS%tkmg3d)%v)%data_3d                   &
                ,iswater=iswater,       xice=domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,           xice_threshold=xice_threshold                                              &
                ,xland=domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di                         & !-- XLAND         land mask (1 for land, 2 for water)  i/o
                ,tsk=domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d            &
                ,lakemask=domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d               & ! 2d var that says lake(1) or not lake(0)
                ,lakeflag=lakeflag                              & ! flag to read in lakemask (lakeflag=1), or to determine lakemask from ivgtyp(i,j)==iswater.and.ht(i,j)>=lake_min_elev (lakeflag=0)
                ,lake_depth_flag=lake_depth_flag,   use_lakedepth=use_lakedepth               & ! flags to use the provided lake depth data (in hi-res input domain file) or not.
                ,tkdry3d=domain%vars_3d(domain%var_indx(kVARS%tkdry3d)%v)%data_3d                  &
                ,tksatu3d=domain%vars_3d(domain%var_indx(kVARS%tksatu3d)%v)%data_3d                  &
                ,lake=lake_or_not                               & ! Logical (:,:) if gridpoint is lake or not (INTENT(OUT)) not used further?
                ,its=its, ite=ite, jts=jts, jte=jte             &
                ,ims=ims, ime=ime, jms=jms, jme=jme             &
                )
        endif

        if (options%physics%snowmodel == kSM_FSM) then
#ifdef FSM
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    SnowModel: FSM2"
            call sm_FSM_init(domain,options,context_change)
            
            !DZs already allocated in alloc_noah above
            !Hard code DZs here just for FSM, since this is also hard-coded in FSM
            DZs = [0.1,0.2,0.4,0.8]
            Zs(1) = DZs(1)/2
            do i = 2,num_soil_layers
                Zs(i) = Zs(i-1) + DZs(i)/2 + DZs(i-1)/2
            end do
            !!
            !allocate(Zs(num_soil_layers))
            !allocate(DZs(num_soil_layers))
            !!
#else
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    User asked to use FSM, but it is not compiled in this version of HICAR"
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Please de-select FSM in the namelist, or recompile HICAR with the FSM library linked"
            stop "FSM2 not compiled in this version of HICAR"
#endif
        endif

        ! defines the height of the middle of the first model level
        z_atm = domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(:,kts,:) - domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d

        update_interval=options%lsm%update_interval
        last_model_time=-999

    end subroutine lsm_init


    subroutine lsm(domain,options,dt)
        implicit none

        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real, intent(in) :: dt
        real :: lsm_dt
        integer :: i,j

        real, allocatable, dimension(:) :: solar_elevation

        if (options%physics%landsurface == 0) return

        if (last_model_time==-999) then
            if (update_interval<=dt) then
                last_model_time = domain%sim_time%seconds()-dt
            else
                last_model_time = domain%sim_time%seconds()-update_interval
            endif
        endif

        if ((domain%sim_time%seconds() - last_model_time) >= update_interval) then
            lsm_dt = domain%sim_time%seconds() - last_model_time
            last_model_time = domain%sim_time%seconds()

            if (options%physics%radiation_downScaling==1  .and. options%physics%radiation>1) then
                SW = domain%vars_2d(domain%var_indx(kVARS%shortwave_total)%v)%data_2d
            else
                SW = domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d
            endif

            ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: snow_water_equivalent max:", MAXVAL(domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d)

            ! exchange coefficients
            windspd = sqrt(domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d**2 + domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d**2)
            
            !If there is no surface layer scheme, then we should compute exchange coefficients here. Assign CHS to CHS2 and CQS2 as well...
            if (options%physics%surfacelayer ==  0) then
                if (exchange_term==1) then
                    call calc_exchange_coefficient(windspd,domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d)
                elseif (exchange_term==2) then
                    call calc_mahrt_holtslag_exchange_coefficient(windspd,domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d)
                endif
                
                domain%vars_2d(domain%var_indx(kVARS%chs2)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d
                domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d
            endif


            if((options%physics%watersurface==kWATER_SIMPLE) .or.      &
                (options%physics%watersurface==kWATER_LAKE) ) then
                    call water_simple(options,                              &
                                      domain%vars_2d(domain%var_indx(kVARS%sst)%v)%data_2d(its:ite,jts:jte),                   &
                                      domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d(its:ite,jts:jte),      &
                                      windspd(its:ite,jts:jte),                              &
                                      domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d(its:ite,kms,jts:jte),       &
                                      domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d(its:ite,kms,jts:jte),       &
                                      domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d(its:ite,jts:jte),         &
                                      domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d(its:ite,jts:jte),           &
                                      domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(its:ite,jts:jte),                     &
                                      QSFC(its:ite,jts:jte),                                 &
                                      domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d(its:ite,jts:jte),                   &
                                      domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d(its:ite,jts:jte),      &
                                      domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d(its:ite,jts:jte),   &
                                      domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di(its:ite,jts:jte),                      &
                                      its, ite, kts, kte, jts, jte)
                                !   ,domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d               & ! terrain height [m] if ht(i,j)>=lake_min_elev -> lake (in case no lake category is provided, but lake model is selected, we need to not run the simple water as well - left comment in for future reference)
            endif

            !___________________ Lake model _____________________
            ! This lake model (ported from WRF V4.4) is run for the grid cells that are defined as lake in the hi-res input file.
            ! It also is advised to supply a lake_depth parameter in the hi-res input, otherwise the default depth of 50m is used (see lakeini above)
            ! It requires the VEGPARM.TBL landuse category to be one which has a separate lake category (i.e. MODIFIED_IGBP_MODIS_NOAH, USGS-RUC or MODI-RUC).
            ! For the grid cells that are defined as water, but not as lake (i.e. oceans), the simple water model above will be run.
            if (options%physics%watersurface==kWATER_LAKE) then    ! WRF's lake model

                ! current_precipitation = (domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2dd-RAINBL)+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE  ! analogous to noah calls

                call lake( &
                    t_phy=domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d                            & !-- t_phy         temperature (K)     !Temprature at the mid points (K)
                    ,p8w=domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d(:,kms:kme,:)         & !-- p8w           pressure at full levels (Pa) ! Naming convention: 8~at => p8w reads as "p-at-w" (w=full levels)
                    ,dz8w=domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d                           & !-- dz8w          dz between full levels (m)
                    ,qvcurr=domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d                          &  !i
                    ,u_phy=domain%vars_3d(domain%var_indx(kVARS%u_mass)%v)%data_3d                                & !-- u_phy         u-velocity interpolated to theta points (m/s)
                    ,v_phy=domain%vars_3d(domain%var_indx(kVARS%v_mass)%v)%data_3d                                & !-- v_phy         v-velocity interpolated to theta points (m/s)
                    ,glw=domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d                                & !-- GLW           downward long wave flux at ground surface (W/m^2)
                    ,emiss=domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d                       & !-- EMISS         surface emissivity (between 0 and 1)
                    ,rainbl=current_precipitation                               & ! RAINBL in mm (Accumulation between PBL calls)
                    ,dtbl=lsm_dt                                                & !-- dtbl          timestep (s) or ITIMESTEP?
                    ,swdown=SW                                                  & !-- SWDOWN        downward short wave flux at ground surface (W/m^2)
                    ,albedo=ALBEDO                                              & ! albedo? fixed at 0.17?
                    ,xlat_urb2d=domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d                         & ! optional ?
                    ,z_lake3d=domain%vars_3d(domain%var_indx(kVARS%z_lake3d)%v)%data_3d                           &
                    ,dz_lake3d=domain%vars_3d(domain%var_indx(kVARS%dz_lake3d)%v)%data_3d                         &
                    ,lakedepth2d=domain%vars_2d(domain%var_indx(kVARS%lakedepth2d)%v)%data_2d                     &
                    ,watsat3d=domain%vars_3d(domain%var_indx(kVARS%watsat3d)%v)%data_3d                           &
                    ,csol3d=domain%vars_3d(domain%var_indx(kVARS%csol3d)%v)%data_3d                               &
                    ,tkmg3d=domain%vars_3d(domain%var_indx(kVARS%tkmg3d)%v)%data_3d                               &
                    ,tkdry3d=domain%vars_3d(domain%var_indx(kVARS%tkdry3d)%v)%data_3d        &
                    ,tksatu3d=domain%vars_3d(domain%var_indx(kVARS%tksatu3d)%v)%data_3d                  &
                    ,ivgtyp=domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di                                     &
                    ,HT=domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d                                  &
                    ,xland=real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di)                               & !-- XLAND         land mask (1 for land, 2 OR 0 for water)  i/o
                    ,iswater=iswater,  xice=domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,   xice_threshold=xice_threshold   &
                    ,lake_min_elev=5.                                           & ! if this value is changed, also change it in lake_ini
                    ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde       &
                    ,ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme       &
                    ,its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte       &
                    ,h2osno2d=domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d             & !domain%h2osno2d%data_2d
                    ,snowdp2d=domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d                        & ! domain%snowdp2d%data_2d
                    ,snl2d=domain%vars_2d(domain%var_indx(kVARS%snl2d)%v)%data_2d                                 &
                    ,z3d=domain%vars_3d(domain%var_indx(kVARS%z3d)%v)%data_3d                                     &
                    ,dz3d=domain%vars_3d(domain%var_indx(kVARS%dz3d)%v)%data_3d                                   &
                    ,zi3d=domain%vars_3d(domain%var_indx(kVARS%zi3d)%v)%data_3d                                   &
                    ,h2osoi_vol3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_vol3d)%v)%data_3d           &
                    ,h2osoi_liq3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_liq3d)%v)%data_3d       &
                    ,h2osoi_ice3d=domain%vars_3d(domain%var_indx(kVARS%h2osoi_ice3d)%v)%data_3d           &
                    ,t_grnd2d=domain%vars_2d(domain%var_indx(kVARS%t_grnd2d)%v)%data_2d               &
                    ,t_soisno3d=domain%vars_3d(domain%var_indx(kVARS%t_soisno3d)%v)%data_3d                                      &
                    ,t_lake3d=domain%vars_3d(domain%var_indx(kVARS%t_lake3d)%v)%data_3d                           & ! 3d lake temperature (K)
                    ,savedtke12d=domain%vars_2d(domain%var_indx(kVARS%savedtke12d)%v)%data_2d                &
                    ,lake_icefrac3d=domain%vars_3d(domain%var_indx(kVARS%lake_icefrac3d)%v)%data_3d    &
                    ,lakemask=domain%vars_2d(domain%var_indx(kVARS%lakemask)%v)%data_2d                                        &
                    ,lakeflag=lakeflag                                          &
                    ,hfx= domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d                          & !(OUT)-- HFX         upward heat flux at the surface (W/m^2)   (INTENT:OUT)
                    ,lh=domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d                              & !(OUT)-- LH          net upward latent heat flux at surface (W/m^2)
                    ,grdflx=domain%vars_2d(domain%var_indx(kVARS%ground_heat_flux)%v)%data_2d                     & !(OUT)-- GRDFLX(I,J) ground heat flux (W m-2)
                    ,tsk=domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d                        & !(OUT)-- TSK          skin temperature [K]
                    ,qfx=domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d                                     & !(OUT)-- QFX        upward moisture flux at the surface (kg/m^2/s) in
                    ,t2= domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d                          & !(OUT)-- t2         diagnostic 2-m temperature from surface layer and lsm
                    ,th2=TH2                                                    & !(OUT)-- th2        diagnostic 2-m theta from surface layer and lsm
                    ,q2=domain%vars_2d(domain%var_indx(kVARS%humidity_2m)%v)%data_2d                              & !(OUT)-- q2         diagnostic 2-m mixing ratio from surface layer and lsm
                )

            endif


            where(windspd<1) windspd=1 ! minimum wind speed to prevent the exchange coefficient from blowing up
            
            
            ! --------------------------------------------------
            ! Now handle the land surface options
            ! --------------------------------------------------
            if (options%physics%landsurface == kLSM_NOAH) then
                ! Call the Noah Land Surface Model

                ! 2m saturated mixing ratio
                do j=jms,jme
                    do i=ims,ime
                        if (domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(i,j) == kLC_LAND) then
                            QGH(i,j) = sat_mr(domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d(i,j),domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm%monthly_albedo) then
                    if (cur_vegmonth /= domain%sim_time%month) then
                        ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :)
                    endif
                endif
                if (options%lsm%monthly_vegfrac) then
                    if (cur_vegmonth /= domain%sim_time%month) then
                        VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, domain%sim_time%month, :)
                        cur_vegmonth = domain%sim_time%month
                    endif
                endif

                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: accumulated_precipitation max:", MAXVAL(domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2dd)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: RAINBL max:", MAXVAL(RAINBL)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: domain%precipitation_bucket max:", MAXVAL(domain%precipitation_bucket)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: rain_bucket max:", MAXVAL(rain_bucket)


                ! RAINBL(i,j) = [kg m-2]   RAINBL = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2dd  ! used to store last time step accumulated precip so that it can be subtracted from the current step
                current_precipitation = (domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d - rainlsm) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE

                CHS = domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d*windspd
                CHS2 = domain%vars_2d(domain%var_indx(kVARS%chs2)%v)%data_2d*windspd
                CQS2 = domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d*windspd


                call lsm_noah(domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,                &
                            domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                   &
                            domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d(:,kms:kme,:),&
                            domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,                   &
                            domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,              &
                            domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,                 &
                            domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d,                           &
                            domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d,                   &
                            domain%vars_2d(domain%var_indx(kVARS%ground_heat_flux)%v)%data_2d,              &
                            QGH,                                          &
                            GSW,                                          &
                            SW,                                           &
                            domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,                      &
                            SMSTAV,                                       &
                            domain%vars_2d(domain%var_indx(kVARS%soil_totalmoisture)%v)%data_2d,            &  ! this is not defined on lsm_init (BK 2021/03/20)
                            SFCRUNOFF,                                    &
                            UDRUNOFF,                                     &
                            domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di,                              &
                            domain%vars_2d(domain%var_indx(kVARS%soil_type)%v)%data_2di,                             &
                            ISURBAN,                                      &
                            ISICE,                                        &
                            VEGFRAC,                                      &
                            ALBEDO,                                       &
                            ALBBCK,                                       &
                            domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,                  &
                            Z0,                                           &
                            domain%vars_2d(domain%var_indx(kVARS%soil_deep_temperature)%v)%data_2d,         &
                            real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di),                       &
                            domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,                                  &
                            domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d,               &
                            EMBCK,                                        &
                            SNOWC,                                        &
                            QSFC,                                         &
                            current_precipitation,                        &  ! RAINBL
                            MMINLU,                                       &
                            num_soil_layers,                              &
                            lsm_dt,                                       &
                            DZS,                                          &
                            ITIMESTEP,                                    &
                            domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d,            &
                            domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,              &
                            domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d,         &
                            domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d,                  &
                            CHS,                                          &
                            CHS2,                                         &
                            CQS2,                                         &
                            CPM,                                          &
                            rcp,                                          &
                            SR,                                           &
                            chklowq,                                      &
                            domain%vars_2d(domain%var_indx(kVARS%lai)%v)%data_2d,                           &
                            qz0,                                          & !H
                            myj,frpcpn,                                   &
                            domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d,        &
                            domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d,                   &     !SNOWH,                                   & !H
                            SNOALB,SHDMIN,SHDMAX,                         & !I
                            SNOTIME,                                      & !?
                            ACSNOM,ACSNOW,                                & !O
                            SNOPCX,                                       & !O
                            POTEVP,                                       & !O
                            SMCREL,                                       & !O
                            XICE_THRESHOLD,                               &
                            RDLAI2D,USEMONALB,                            &
                            Ri,                                           & !I
                            NOAHRES,                                      &
                            ua_phys,flx4_2d,fvb_2d,fbur_2d,fgsn_2d,       & ! Noah UA changes
                            ids,ide, jds,jde, kds,kde,                    &
                            ims,ime, jms,jme, kms,kme,                    &
                            its,ite, jts,jte, kts,kte)

                where(domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d > options%lsm%max_swe) domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = options%lsm%max_swe
                ! now that znt (roughness_z0) has been updated, we need to recalculate terms
            else if (options%physics%landsurface == kLSM_NOAHMP) then
            ! Call the Noah-MP Land Surface Model

                ! 2m saturated mixing ratio
                do j=jms,jme
                    do i=ims,ime
                        if (domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(i,j) == kLC_LAND) then
                            QGH(i,j) = sat_mr(domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d(i,j),domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d(i,j))
                        endif
                    enddo
                enddo
                if (options%lsm%monthly_albedo) then
                    ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :)
                else
                    ALBEDO = domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, 1, :)
                endif
                if (options%lsm%monthly_vegfrac) then
                    if (cur_vegmonth /= domain%sim_time%month) then
                        VEGFRAC = domain%vars_3d(domain%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:, domain%sim_time%month, :)
                        cur_vegmonth = domain%sim_time%month
                    endif
                endif

                !more parameters
                landuse_name = options%lsm%LU_Categories            !test whether this works or if we need something separate

                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: accumulated_precipitation max:", MAXVAL(domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: RAINBL max:", MAXVAL(RAINBL)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: domain%precipitation_bucket max:", MAXVAL(domain%precipitation_bucket)
                ! if (STD_OUT_PE .and. .not.context_change) write(*,*) "    lsm start: rain_bucket max:", MAXVAL(rain_bucket)
                
                current_snow = (domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d-snowlsm)!+(domain%snowfall_bucket-snow_bucket)*kPRECIP_BUCKET_SIZE !! MJ: snowfall in kg m-2
                current_precipitation = (domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d - rainlsm) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE

!                do I = ims,ime
!                  do J = jms,jme
!                    call calc_declin(domain%sim_time%day_of_year(),real(domain%sim_time%hour),real(domain%sim_time%minute),real(domain%sim_time%second),domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d(I,J),domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d(I,J),domain%cos_zenith%data_2d(I,J))
!                  enddo
!                enddo

                allocate(solar_elevation(ims:ime))

                do j = jms,jme

                    solar_elevation  = calc_solar_elevation(date=domain%sim_time, tzone=options%rad%tzone, &
                        lon=domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d, lat=domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d, j=j, &
                        ims=ims,ime=ime,jms=jms,jme=jme,its=its,ite=ite)
                    domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d(ims:ime,j)=sin(solar_elevation(ims:ime))
                enddo

                if (options%physics%snowmodel==kSM_FSM) then
                    SR = 0.0 ! This, in combination with setting OPT_SNF to 4 in the LSM_init, will turn off snowfall partitioning in NoahMP
                    current_precipitation = current_precipitation-current_snow ! Now remove snowfall from precipitation, so we only have liquid precip going into NMP
                    
                    nmp_snowh = 0.0
                    nmp_snow = 0.0
                else
                    SR = current_snow/(current_precipitation+epsilon)
                    nmp_snowh = domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d
                    nmp_snow = domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d
                endif

                call noahmplsm(ITIMESTEP,                              &
                             domain%sim_time%year,                   &
                             domain%sim_time%day_of_year(),          &
                             domain%vars_2d(domain%var_indx(kVARS%cosine_zenith_angle)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,                  &
                             domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d,                 &
                             domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,              & ! domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d,              & !
                             lsm_dt,                                   &
                             DZS,                                      &
                             num_soil_layers,                          &
                             domain%dx,                                &
                             domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di,                          &
                             domain%vars_2d(domain%var_indx(kVARS%soil_type)%v)%data_2di,                         &
                             VEGFRAC,                                  &
                             domain%vars_2d(domain%var_indx(kVARS%vegetation_fraction_max)%v)%data_2d,   &
                             domain%vars_2d(domain%var_indx(kVARS%soil_deep_temperature)%v)%data_2d,     &
                             real(domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di),                   &
                             domain%vars_2d(domain%var_indx(kVARS%xice)%v)%data_2d,                              &
                             XICE_THRESHOLD,                           &
                             domain%vars_2d(domain%var_indx(kVARS%crop_category)%v)%data_2di,                     &  !only used if iopt_crop>0; not currently set up
                             domain%vars_2d(domain%var_indx(kVARS%date_planting)%v)%data_2d,             &  !only used if iopt_crop>0; not currently set up
                             domain%vars_2d(domain%var_indx(kVARS%date_harvest)%v)%data_2d,              &  !only used if iopt_crop>0; not currently set up
                             domain%vars_2d(domain%var_indx(kVARS%growing_season_gdd)%v)%data_2d,        &  !only used if iopt_crop>0; not currently set up
                             IDVEG, IOPT_CRS,  IOPT_BTR, IOPT_RUN,     &
                             IOPT_SFC, IOPT_FRZ, IOPT_INF, IOPT_RAD,   &
                             IOPT_ALB, IOPT_SNF, IOPT_TBOT, IOPT_STC,  &
                             IOPT_GLA, IOPT_RSF, IOPT_SOIL,IOPT_PEDO,  &
                             IOPT_CROP, IOPT_IRR, IOPT_IRRM, IZ0TLND,  &
                             SF_URBAN_PHYSICS,                         &
                             domain%vars_3d(domain%var_indx(kVARS%soil_sand_and_clay)%v)%data_3d,        &  ! only used if iopt_soil = 3
                             domain%vars_2d(domain%var_indx(kVARS%soil_texture_1)%v)%data_2d,            &  ! only used if iopt_soil = 2
                             domain%vars_2d(domain%var_indx(kVARS%soil_texture_2)%v)%data_2d,            &  ! only used if iopt_soil = 2
                             domain%vars_2d(domain%var_indx(kVARS%soil_texture_3)%v)%data_2d,            &  ! only used if iopt_soil = 2
                             domain%vars_2d(domain%var_indx(kVARS%soil_texture_4)%v)%data_2d,            &  ! only used if iopt_soil = 2
                             domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d,               &
                             domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,               &
                             domain%vars_3d(domain%var_indx(kVARS%u_mass)%v)%data_3d,                    &
                             domain%vars_3d(domain%var_indx(kVARS%v_mass)%v)%data_3d,                    &
                             SW,                                       &
                             domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d,          &  ! only used in urban modules, which are currently disabled
                             domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d,         &  ! only used in urban modules, which are currently disabled
                             domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d,                  &
                             domain%vars_3d(domain%var_indx(kVARS%pressure_interface)%v)%data_3d(:,kms:kme,:),&
                             current_precipitation,                    &
                             SR,                                       &
                             domain%vars_2d(domain%var_indx(kVARS%irr_frac_total)%v)%data_2d,            &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_frac_sprinkler)%v)%data_2d,        &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_frac_micro)%v)%data_2d,            &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_frac_flood)%v)%data_2d,            &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,          &  ! TSK
                             domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,             &  !  HFX
                             domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d,                       &
                             domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d,               &  ! LH
                             domain%vars_2d(domain%var_indx(kVARS%ground_heat_flux)%v)%data_2d,          &  ! GRDFLX
                             SMSTAV,                                   &
                             domain%vars_2d(domain%var_indx(kVARS%soil_totalmoisture)%v)%data_2d,        &
                             SFCRUNOFF, UDRUNOFF,                      &
                             ALBEDO, SNOWC,                            &
                             domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d,        &
                             domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d,    &
                             domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d,          &
                             nmp_snow,                                 &
                             nmp_snowh,                                &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_water)%v)%data_2d,              &
                             ACSNOM, ACSNOW,                           &
                             domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d, QSFC, Z0, &
                             domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d,              &
                             domain%vars_2d(domain%var_indx(kVARS%hpbl)%v)%data_2d,                      &
                             domain%vars_2d(domain%var_indx(kVARS%irr_eventno_sprinkler)%v)%data_2di,             &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_eventno_micro)%v)%data_2di,                 &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_eventno_flood)%v)%data_2di,                 &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_alloc_sprinkler)%v)%data_2d,       &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_alloc_micro)%v)%data_2d,           &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_alloc_flood)%v)%data_2d,           &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_evap_loss_sprinkler)%v)%data_2d,   &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_amt_sprinkler)%v)%data_2d,         &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_amt_micro)%v)%data_2d,             &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%irr_amt_flood)%v)%data_2d,             &  ! only used if iopt_irr > 0
                             domain%vars_2d(domain%var_indx(kVARS%evap_heat_sprinkler)%v)%data_2d,       &  ! only used if iopt_irr > 0
                             landuse_name,                             &
                             domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di,             &
                             domain%vars_2d(domain%var_indx(kVARS%veg_leaf_temperature)%v)%data_2d,      &
                             domain%vars_2d(domain%var_indx(kVARS%ground_surf_temperature)%v)%data_2d,   &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_water_ice)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_water_liquid)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_vapor_pressure)%v)%data_2d,     &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_temperature)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%coeff_momentum_drag)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d,                       &
                             domain%vars_2d(domain%var_indx(kVARS%canopy_fwet)%v)%data_2d,               &
                             domain%vars_2d(domain%var_indx(kVARS%snow_water_eq_prev)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%snow_albedo_prev)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%snowfall_ground)%v)%data_2d,           &
                             domain%vars_2d(domain%var_indx(kVARS%rainfall_ground)%v)%data_2d,           &
                             domain%vars_2d(domain%var_indx(kVARS%storage_lake)%v)%data_2d,              &
                             domain%vars_2d(domain%var_indx(kVARS%water_table_depth)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%water_aquifer)%v)%data_2d,             &
                             domain%vars_2d(domain%var_indx(kVARS%storage_gw)%v)%data_2d,                &
                             domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(:,1:3,:), &
                             domain%vars_3d(domain%var_indx(kVARS%snow_layer_depth)%v)%data_3d(:,1:7,:), &
                             domain%vars_3d(domain%var_indx(kVARS%snow_layer_ice)%v)%data_3d(:,1:3,:),   &
                             domain%vars_3d(domain%var_indx(kVARS%snow_layer_liquid_water)%v)%data_3d(:,1:3,:),&
                             domain%vars_2d(domain%var_indx(kVARS%mass_leaf)%v)%data_2d,                 &
                             domain%vars_2d(domain%var_indx(kVARS%mass_root)%v)%data_2d,                 &
                             domain%vars_2d(domain%var_indx(kVARS%mass_stem)%v)%data_2d,                 &
                             domain%vars_2d(domain%var_indx(kVARS%mass_wood)%v)%data_2d,                 &
                             domain%vars_2d(domain%var_indx(kVARS%soil_carbon_stable)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%soil_carbon_fast)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%lai)%v)%data_2d,                       &
                             domain%vars_2d(domain%var_indx(kVARS%sai)%v)%data_2d,                       &
                             domain%vars_2d(domain%var_indx(kVARS%snow_age_factor)%v)%data_2d,           &
                             domain%vars_3d(domain%var_indx(kVARS%eq_soil_moisture)%v)%data_3d,          &
                             domain%vars_2d(domain%var_indx(kVARS%smc_watertable_deep)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%recharge_deep)%v)%data_2d,             &
                             domain%vars_2d(domain%var_indx(kVARS%recharge)%v)%data_2d,                  &
                             domain%vars_2d(domain%var_indx(kVARS%mass_ag_grain)%v)%data_2d,             &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%vars_2d(domain%var_indx(kVARS%growing_degree_days)%v)%data_2d,       &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%vars_2d(domain%var_indx(kVARS%plant_growth_stage)%v)%data_2di,                &  ! currently left as zeroes; not used if iopt_crop = 0?
                             domain%vars_3d(domain%var_indx(kVARS%gecros_state)%v)%data_3d,              &  ! not set up; only used if iopt_crop = 2
                             domain%vars_2d(domain%var_indx(kVARS%temperature_2m_veg)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%temperature_2m_bare)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%mixing_ratio_2m_veg)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%mixing_ratio_2m_bare)%v)%data_2d,      &
                             domain%vars_2d(domain%var_indx(kVARS%surface_rad_temperature)%v)%data_2d,   &
              	             domain%vars_2d(domain%var_indx(kVARS%net_ecosystem_exchange)%v)%data_2d,    &
                             domain%vars_2d(domain%var_indx(kVARS%gross_primary_prod)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%net_primary_prod)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%vegetation_fraction_out)%v)%data_2d,   &
                             domain%vars_2d(domain%var_indx(kVARS%runoff_surface)%v)%data_2d,            &
                             domain%vars_2d(domain%var_indx(kVARS%runoff_subsurface)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%evap_canopy)%v)%data_2d,               &
                             domain%vars_2d(domain%var_indx(kVARS%evap_soil_surface)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%transpiration_rate)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%rad_absorbed_total)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%rad_net_longwave)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%apar)%v)%data_2d,                      &
                             domain%vars_2d(domain%var_indx(kVARS%photosynthesis_total)%v)%data_2d,      &
                             domain%vars_2d(domain%var_indx(kVARS%rad_absorbed_veg)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%rad_absorbed_bare)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%stomatal_resist_sun)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%stomatal_resist_shade)%v)%data_2d,     &
                             domain%vars_2d(domain%var_indx(kVARS%frac_between_gap)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%frac_within_gap)%v)%data_2d,           &
                             domain%vars_2d(domain%var_indx(kVARS%ground_temperature_canopy)%v)%data_2d, &
                             domain%vars_2d(domain%var_indx(kVARS%ground_temperature_bare)%v)%data_2d,   &
                             domain%vars_2d(domain%var_indx(kVARS%ch_veg)%v)%data_2d,                    &
                             domain%vars_2d(domain%var_indx(kVARS%ch_bare)%v)%data_2d,                   &
                             domain%vars_2d(domain%var_indx(kVARS%sensible_heat_veg)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%sensible_heat_canopy)%v)%data_2d,      &
                             domain%vars_2d(domain%var_indx(kVARS%sensible_heat_bare)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%evap_heat_veg)%v)%data_2d,             &
                             domain%vars_2d(domain%var_indx(kVARS%evap_heat_bare)%v)%data_2d,            &
                             domain%vars_2d(domain%var_indx(kVARS%ground_heat_veg)%v)%data_2d,           &
                             domain%vars_2d(domain%var_indx(kVARS%ground_heat_bare)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%net_longwave_veg)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%net_longwave_canopy)%v)%data_2d,       &
                             domain%vars_2d(domain%var_indx(kVARS%net_longwave_bare)%v)%data_2d,         &
                             domain%vars_2d(domain%var_indx(kVARS%transpiration_heat)%v)%data_2d,        &
                             domain%vars_2d(domain%var_indx(kVARS%evap_heat_canopy)%v)%data_2d,          &
                             domain%vars_2d(domain%var_indx(kVARS%ch_leaf)%v)%data_2d,                   &
                             domain%vars_2d(domain%var_indx(kVARS%ch_under_canopy)%v)%data_2d,           &
                             domain%vars_2d(domain%var_indx(kVARS%ch_veg_2m)%v)%data_2d,                 &
                             domain%vars_2d(domain%var_indx(kVARS%ch_bare_2m)%v)%data_2d,                &
                             domain%vars_2d(domain%var_indx(kVARS%stomatal_resist_total)%v)%data_2d,     &
                             ids,ide,  jds,jde,  kds,kde,              &
                             ims,ime,  jms,jme,  kms,kme,              &
                             its,ite,  jts,jte,  kts,kte)
                             
                             
                if (options%lsm%monthly_albedo) then
                    domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, domain%sim_time%month, :) = ALBEDO
                else
                    domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(:, 1, :) = ALBEDO
                endif

                VEGFRAC = domain%vars_2d(domain%var_indx(kVARS%vegetation_fraction_out)%v)%data_2d(:, :)*100.0

                if ( .not.(options%physics%snowmodel==kSM_FSM)) then
                    domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d = nmp_snowh
                    domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = nmp_snow
                endif

    !         TLE: OMITTING OPTIONAL PRECIP INPUTS FOR NOW
    !                         MP_RAINC, MP_RAINNC, MP_SHCV, MP_SNOW, MP_GRAUP, MP_HAIL     )
                where(domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d > options%lsm%max_swe) domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = options%lsm%max_swe
            endif

            !! MJ added: this block is for FSM as sm.
            if (options%physics%snowmodel == kSM_FSM) then
                current_precipitation = (domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d - rainlsm) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE
                current_snow = (domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d-snowlsm)!+(domain%snowfall_bucket-snow_bucket)*kPRECIP_BUCKET_SIZE !! MJ: snowfall in kg m-2
                current_rain = max(current_precipitation-current_snow,0.) !! MJ: rainfall in kg m-2
                !!
                domain%vars_2d(domain%var_indx(kVARS%windspd_10m)%v)%data_2d(its:ite,jts:jte)=windspd(its:ite,jts:jte)
                !!
#ifdef FSM
                call sm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)
#endif
                !!
                
                !do i = 1, num_soil_layers              ! soil
                !   domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d(its:ite,i,jts:jte) =  domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(its:ite,i,jts:jte)    / (1000. * DZS(i))
                !end do
                
            endif
            !!
            if (options%physics%landsurface > kLSM_BASIC) then
            
                rainlsm = domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d
                snowlsm = domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d
                
                domain%vars_2d(domain%var_indx(kVARS%longwave_up)%v)%data_2d = STBOLT * domain%vars_2d(domain%var_indx(kVARS%land_emissivity)%v)%data_2d * domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d**4
                ! accumulate soil moisture over the entire column
                domain%vars_2d(domain%var_indx(kVARS%soil_totalmoisture)%v)%data_2d = domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(:,1,:) * DZS(1) * 1000
                do i = 2,num_soil_layers
                    domain%vars_2d(domain%var_indx(kVARS%soil_totalmoisture)%v)%data_2d = domain%vars_2d(domain%var_indx(kVARS%soil_totalmoisture)%v)%data_2d + domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(:,i,:) * DZS(i) * 1000
                enddo

                if (options%physics%surfacelayer == 0) then
                    ! 2m Air T and Q are not well defined if Tskin is not coupled with the surface fluxes
                    call surface_diagnostics(domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d,          &
                                            domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d,                    &
                                            domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d,       &
                                            QSFC,                                  &
                                            domain%vars_2d(domain%var_indx(kVARS%chs2)%v)%data_2d,                   &
                                            domain%vars_2d(domain%var_indx(kVARS%cqs2)%v)%data_2d,                   &
                                            domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d,         &
                                            domain%vars_2d(domain%var_indx(kVARS%humidity_2m)%v)%data_2d,            &
                                            domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d,       &
                                            (VEGFRAC/100.0),                       &
                                            domain%vars_2d(domain%var_indx(kVARS%veg_type)%v)%data_2di,                       &
                                            domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di,                      &
                                            domain%vars_2d(domain%var_indx(kVARS%temperature_2m_veg)%v)%data_2d,     &
                                            domain%vars_2d(domain%var_indx(kVARS%temperature_2m_bare)%v)%data_2d,    &
                                            domain%vars_2d(domain%var_indx(kVARS%mixing_ratio_2m_veg)%v)%data_2d,    &
                                            domain%vars_2d(domain%var_indx(kVARS%mixing_ratio_2m_bare)%v)%data_2d)
                endif

            endif
            !!
        endif
        

    end subroutine lsm
end module land_surface
