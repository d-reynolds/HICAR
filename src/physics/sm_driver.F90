!>----------------------------------------------------------
!! This module provides a wrapper to call various snow models
!! It sets up variables specific to the physics package to be used including both
!!
!! The main entry point to the code is snow_model(domain,options,lsm_dt)
!!
!! <pre>
!! Call tree graph :
!!  sm_init->[ external initialization routines]
!!  sm->[  external sm routines]
!!  sm_finalize
!!
!! High level routine descriptions / purpose
!!   sm_init           - initializes physics package
!!   sm                - sets up and calls main physics package
!!   sm_finalize       - permits physics package cleanup (close files, deallocate memory)
!!
!! Inputs: domain, options, lsm_dt
!!      domain,options  = as defined in data_structures
!!      lsm_dt          = time step (seconds)
!! </pre>
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@epfl.ch)
!!
!!----------------------------------------------------------
module snow_model_driver
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use mod_wrf_constants,   only : KARMAN, gravity, cp, R_d, rcp, EP_1, EP_2, SVPT0, SVP1, SVP2, SVP3, EOMEG, STBOLT, p1000mb, XLV
    use ieee_arithmetic ! for debugging
    use icar_constants
#ifdef FSM
    use module_sm_FSMdrv,   only : sm_FSM_init, sm_FSM
#endif
#ifdef SNOWPACK
#ifdef SNOWPACK_FORTRAN
    use snowpack_driver, only : snowpack_init, snowpack_step
#else
    use module_sm_SNOWPACKdrv, only : sm_snowpack_init, sm_SNOWPACK
#endif
#endif
    use snow_drift,            only : snow_drift_var_request, snow_drift_init
    use module_snowslide,      only : snowslide_var_request, snowslide_init, snowslide_step
    use NoahmpIOVarType,            only : NoahmpIO_type
    use NoahmpDriverMainMod,        only : noahmp
    use ConfigVarInTransferMod,      only : ConfigVarInTransfer
    use ForcingVarInTransferMod,     only : ForcingVarInTransfer
    use EnergyVarInTransferMod,      only : EnergyVarInTransfer
    use WaterVarInTransferMod,       only : WaterVarInTransfer
    use GroundThermalPropertyMod,    only : GroundThermalProperty
    use SoilSnowTemperatureMainMod,  only : SoilSnowTemperatureMain
    use SoilSnowWaterPhaseChangeMod, only : SoilSnowWaterPhaseChange
    use SoilWaterMainMod,            only : SoilWaterMain

    implicit none

    private
    public :: sm_var_request, sm_init, snow_model   !, sm_finalize

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: restart
    integer :: update_interval
    real*8  :: last_model_time(kMAX_NESTS), next_update_time(kMAX_NESTS)
    real*8  :: last_snowslide_time(kMAX_NESTS), next_snowslide_time(kMAX_NESTS)
    real :: lsm_dt

    ! Cadence for SNOWSLIDE gravitational redistribution. Snowslide runs
    ! once per hour (at the first snow-model tick on or after each hour
    ! boundary), immediately before the snow model is invoked, so that the
    ! snow model's built-in merge/split/metamorphism passes process the
    ! fresh deposit on the same step.
    real*8, parameter :: snowslide_interval_s = 3600.0d0
    
contains

    subroutine sm_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%snowmodel == kSM_FSM) then

            call options%alloc_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg, kVARS%slope_angle, kVARS%lsm_last_snow,&
                         kVARS%lsm_last_precip, kVARS%QFX, kVARS%chs, kVARS%land_emissivity,    &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%snow_nlayers,   &
                         kVARS%shd, kVARS%meltflux_out_tstep, kVARS%meltflux_out_cumul, kVARS%Sliq_out, &
                         kVARS%windspd_10m, kVARS%bs_saltation_flux, kVARS%bs_saltation_height, kVARS%bs_saltation_concentration, &
                         kVARS%bs_suspension_flux, kVARS%bs_sublimation_flux, kVARS%bs_drift_swe_salt, kVARS%bs_drift_swe_susp, kVARS%bs_drift_swe_subl, kVARS%dSWE_slide, &
                         kVARS%shortwave_direct, kVARS%shortwave_diffuse, kVARS%ground_surf_temperature])
             
             call options%restart_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,  kVARS%snowfall, kVARS%albedo, kVARS%QFX, kVARS%land_emissivity,            &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0, kVARS%lsm_last_snow, kVARS%lsm_last_precip,      &
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%snow_nlayers, &
                         kVARS%meltflux_out_cumul ])

        else if (options%physics%snowmodel == kSM_SNOWPACK) then

            call options%alloc_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%shortwave,     &
                         kVARS%longwave, kVARS%vegetation_fraction, kVARS%canopy_water, kVARS%snow_water_equivalent,    &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%soil_deep_temperature, kVARS%roughness_z0, kVARS%ustar,        &
                         kVARS%snow_height, kVARS%lai, kVARS%temperature_2m_veg, kVARS%slope_angle, kVARS%lsm_last_snow,&
                         kVARS%lsm_last_precip, kVARS%sm_last_snow, kVARS%sm_last_precip,                              &
                         kVARS%QFX, kVARS%chs, kVARS%land_emissivity,    &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%snow_temperature, kVARS%Ds, kVARS%snow_temperature_i, kVARS%Vol_Frac_I, kVARS%Vol_Frac_W, &
                         kVARS%Vol_Frac_A, kVARS%Vol_Frac_S, kVARS%Vol_Frac_WP, kVARS%Rg, kVARS%Rb, kVARS%Dd, kVARS%Sp, kVARS%mk, &
                         kVARS%snow_nlayers, kVARS%bs_threshold_ustar, kVARS%bs_saltation_flux, kVARS%bs_saltation_height, kVARS%bs_saltation_concentration, &
                         kVARS%bs_swe_exchange, kVARS%bs_swe_erode_max, kVARS%mass_hoar, kVARS%CDot, kVARS%snow_stress, kVARS%N3, kVARS%depositionDate, kVARS%dSWE_subl, &
                         kVARS%meltflux_out_tstep, kVARS%meltflux_out_cumul, kVARS%snow_basal_heat_flux, kVARS%soil_water_content_liq, &
                         kVARS%runoff_surface, kVARS%runoff_subsurface])

             call options%restart_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,  kVARS%snowfall, kVARS%albedo, kVARS%QFX, kVARS%land_emissivity,            &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0, kVARS%lsm_last_snow, kVARS%lsm_last_precip,      &
                         kVARS%sm_last_snow, kVARS%sm_last_precip,                                                      &
                         kVARS%snow_temperature, kVARS%Ds, kVARS%snow_temperature_i, kVARS%Vol_Frac_I, kVARS%Vol_Frac_W, &
                         kVARS%Vol_Frac_A, kVARS%Vol_Frac_S, kVARS%Vol_Frac_WP, kVARS%Rg, kVARS%Rb, kVARS%Dd, kVARS%Sp, kVARS%mk, &
                         kVARS%snow_nlayers,                                                                      &
                         kVARS%mass_hoar, kVARS%CDot, kVARS%snow_stress, kVARS%N3, kVARS%depositionDate,           &
                         kVARS%meltflux_out_cumul])
        endif

        ! CRYOWRF-style blowing snow drift (works with both FSM and SNOWPACK)
        if (options%sm%suspension_layer == 1) then
            call snow_drift_var_request(options)
        endif

        ! SNOWSLIDE gravitational redistribution (works with any snow model)
        if (options%sm%snowslide == 1) then
            call snowslide_var_request(options)
        endif

    end subroutine sm_var_request


    subroutine sm_init(domain,options,context_chng)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        logical, optional,  intent(in)      :: context_chng

        logical :: context_change
        real*8 :: elapsed, eff_interval
        integer :: n_calls

        if (options%physics%snowmodel == 0) return

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .false.
        endif

        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        if (STD_OUT_PE .and. .not.context_change) write(*,*) "Initializing Snow Model"

        if (options%physics%snowmodel==kSM_FSM) then
#ifdef FSM
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    SnowModel: FSM2"
            call sm_FSM_init(domain,options,context_change)
#else
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    User asked to use FSM, but it is not compiled in this version of HICAR"
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Please de-select FSM as the snow model in the namelist, or recompile HICAR with the FSM library linked"
            stop "FSM2 not compiled in this version of HICAR"
#endif
        else if (options%physics%snowmodel==kSM_SNOWPACK) then
#ifdef SNOWPACK
#ifdef SNOWPACK_FORTRAN
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    SnowModel: Snowpack (GPU/OpenACC)"
            call snowpack_init(domain,options,context_change)
#else

            ! --- Update host: SNOWPACK-specific inputs ---
            !$acc update host( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d)

            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    SnowModel: Snowpack"
            call sm_snowpack_init(domain,options,context_change)

            ! --- Update device: SNOWPACK outputs ---
            !$acc update device( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d)

#endif
#else
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    User asked to use Snowpack, but it is not compiled in this version of HICAR"
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Please de-select Snowpack as the snow model in the namelist, or recompile HICAR with the Snowpack library linked"
            stop "Snowpack not compiled in this version of HICAR"
#endif

        endif

        ! Initialize CRYOWRF-style blowing snow drift
        if (options%sm%suspension_layer == 1) then
            call snow_drift_init(domain, options)
        endif

        ! Initialize SNOWSLIDE gravitational redistribution
        if (options%sm%snowslide == 1) then
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Initializing SNOWSLIDE gravitational redistribution"
            call snowslide_init(domain, options)
        endif

        update_interval=options%lsm%update_interval
        if (.not.(context_change)) then
            if (update_interval<=10) then
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-10
                next_update_time(domain%nest_indx) = domain%sim_time%seconds()
            else
                last_model_time(domain%nest_indx) = domain%sim_time%seconds()-update_interval
                next_update_time(domain%nest_indx) = domain%sim_time%seconds()
            endif
            if (options%restart%restart) then
                ! Determine when radiation was last called before the restart time,
                ! based on the original simulation start time and the update interval
                eff_interval = max(dble(update_interval), 10.0d0)
                ! add a small fraction of a second in the case of roundoff errors restart time
                elapsed = (options%restart%restart_time%seconds() - 0.01) - options%general%start_time%seconds()
                n_calls = int(elapsed / eff_interval)
                last_model_time(domain%nest_indx) = options%general%start_time%seconds() + n_calls * eff_interval
                next_update_time(domain%nest_indx) = last_model_time(domain%nest_indx) + eff_interval - 0.01
            endif

            ! Snowslide cadence: fire on the first snow-model tick, then every
            ! snowslide_interval_s afterwards. On restart, align to the last
            ! complete hour since start_time so timing is consistent.
            if (options%sm%snowslide == 1) then
                last_snowslide_time(domain%nest_indx) = domain%sim_time%seconds() - snowslide_interval_s
                next_snowslide_time(domain%nest_indx) = domain%sim_time%seconds()
                if (options%restart%restart) then
                    elapsed = (options%restart%restart_time%seconds() - 0.01) - options%general%start_time%seconds()
                    n_calls = int(elapsed / snowslide_interval_s)
                    last_snowslide_time(domain%nest_indx) = options%general%start_time%seconds() &
                                                            + n_calls * snowslide_interval_s
                    next_snowslide_time(domain%nest_indx) = last_snowslide_time(domain%nest_indx) &
                                                            + snowslide_interval_s - 0.01d0
                endif
            endif
        endif

    end subroutine sm_init

    subroutine snow_model(domain,options,dt,NoahmpIO_arg)
        implicit none
        type(domain_t),                  intent(inout) :: domain
        type(options_t),                 intent(in)    :: options
        real,                            intent(in)    :: dt
        type(NoahmpIO_type), optional,   intent(inout) :: NoahmpIO_arg

        real, allocatable, dimension(:,:) :: windspd
        real, allocatable, dimension(:,:) :: current_precipitation
        real, allocatable, dimension(:,:) :: current_snow
        real, allocatable, dimension(:,:) :: current_rain
        integer :: i, j

        if (options%physics%snowmodel == 0) return

        if ((domain%sim_time%seconds()) >= next_update_time(domain%nest_indx)) then
            lsm_dt = domain%sim_time%seconds() - last_model_time(domain%nest_indx)
            last_model_time(domain%nest_indx) = domain%sim_time%seconds() 
            next_update_time(domain%nest_indx) = next_update_time(domain%nest_indx) + update_interval

            allocate(current_precipitation(ims:ime,jms:jme), source=0.0)
            allocate(windspd(ims:ime,jms:jme), source=1.0)
            allocate(current_snow(ims:ime,jms:jme), source=0.0) ! MJ added 
            allocate(current_rain(ims:ime,jms:jme), source=0.0) ! MJ added 


            ! --- Update host: common snow model inputs ---
            ! For SNOWPACK we read sm_last_* (snapshots taken at the last
            ! layer deposit) instead of lsm_last_*; this lets sub-threshold
            ! snowfall accumulate across calls until enough mass exists to
            ! form a layer. For FSM the lsm_last_* path is used, since FSM
            ! has no layer-formation threshold and lsm_last_* updates each
            ! LSM tick.
            if (options%physics%snowmodel == kSM_SNOWPACK) then
                !$acc update host( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sm_last_precip)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sm_last_snow)%v)%data_2d)
            else
                !$acc update host( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%lsm_last_precip)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%lsm_last_snow)%v)%data_2d)
            endif

            associate( &
                u_10m => domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
                v_10m => domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
                precipitation => domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
                snowfall => domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d, &
                lsm_last_precip => domain%vars_2d(domain%var_indx(kVARS%lsm_last_precip)%v)%data_2d, &
                lsm_last_snow => domain%vars_2d(domain%var_indx(kVARS%lsm_last_snow)%v)%data_2d &
                )

            windspd = sqrt(u_10m**2 + v_10m**2)
            where(windspd<1) windspd=1 ! minimum wind speed to prevent the exchange coefficient from blowing up

            ! Setup the input data for the snow models. For SNOWPACK use the
            ! per-deposit snapshots so held mass below hn<0.001 m carries
            ! across calls; for FSM use the per-call lsm_last_* snapshots.
            if (options%physics%snowmodel == kSM_SNOWPACK) then
                associate(sm_last_precip => domain%vars_2d(domain%var_indx(kVARS%sm_last_precip)%v)%data_2d, &
                          sm_last_snow   => domain%vars_2d(domain%var_indx(kVARS%sm_last_snow)%v)%data_2d)
                current_precipitation = (precipitation - sm_last_precip)
                current_snow          = (snowfall      - sm_last_snow)
                end associate
            else
                current_precipitation = (precipitation - lsm_last_precip)
                current_snow          = (snowfall      - lsm_last_snow)
            endif
            current_rain = max(current_precipitation-current_snow,0.) !! MJ: rainfall in kg m-2

            ! SNOWSLIDE gravitational redistribution runs on its own hourly
            ! cadence, immediately before the snow model. This lets the snow
            ! model's built-in merge/split/metamorphism/creep passes process
            ! the fresh deposit on the same tick it was inserted.
            if (options%sm%snowslide == 1 .and. &
                domain%sim_time%seconds() >= next_snowslide_time(domain%nest_indx)) then
                last_snowslide_time(domain%nest_indx) = domain%sim_time%seconds()
                next_snowslide_time(domain%nest_indx) = next_snowslide_time(domain%nest_indx) &
                                                        + snowslide_interval_s
                call snowslide_step(domain, options)
            endif

            if (options%physics%snowmodel==kSM_FSM) then
#ifdef FSM
                ! --- Update host: FSM-specific inputs ---
                !$acc update host( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%windspd_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_flux)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_concentration)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_suspension_flux)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_sublimation_flux)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_salt)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_susp)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_subl)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%humidity_2m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%QFX)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%ground_surf_temperature)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%runoff_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%Sliq_out)%v)%data_2d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d)

                domain%vars_2d(domain%var_indx(kVARS%windspd_10m)%v)%data_2d(its:ite,jts:jte)=windspd(its:ite,jts:jte)

                call sm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)

                ! --- Update device: FSM outputs ---
                !$acc update device( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%windspd_10m)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_flux)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_saltation_concentration)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_suspension_flux)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_salt)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_susp)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%bs_drift_swe_subl)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%Sliq_out)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%QFX)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%ground_surf_temperature)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%runoff_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_cumul)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d)
#endif

            else if (options%physics%snowmodel==kSM_SNOWPACK) then
#ifdef SNOWPACK
#ifdef SNOWPACK_FORTRAN
                ! GPU SNOWPACK: data stays on device, no host/device transfers needed
                !$acc enter data copyin(current_rain, current_snow, windspd)
                call snowpack_step(domain,options,lsm_dt,current_rain,current_snow,windspd)
                !$acc exit data delete(current_rain, current_snow, windspd)

                ! Hand soil-T evolution on snow cells over to NoahMP, using the
                ! conductive heat flux at the snowpack base (G_base, computed by
                ! snowpack_step) as the upper Neumann BC.
                if (options%physics%landsurface == kLSM_NOAHMP .and. present(NoahmpIO_arg)) then
                    call noahmp_soil_only(domain, options, lsm_dt, NoahmpIO_arg)
                end if
#else

                ! --- Update host: SNOWPACK-specific inputs ---
                !$acc update host( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%pressure)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_cumul)%v)%data_2d)

                call sm_SNOWPACK(domain,options,lsm_dt,current_rain,current_snow,windspd)

                ! --- Update device: SNOWPACK outputs ---
                !$acc update device( &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%albedo)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%depositionDate)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature_i)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_I)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_W)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_A)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_S)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Vol_Frac_WP)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Rg)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Rb)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Dd)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%Sp)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%mk)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%mass_hoar)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%CDot)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_stress)%v)%data_3d, &
                !$acc   domain%vars_3d(domain%var_indx(kVARS%N3)%v)%data_3d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
                !$acc   domain%vars_2d(domain%var_indx(kVARS%meltflux_out_cumul)%v)%data_2d)
#endif
#endif
            endif
            end associate

            ! Accumulate cumulative melt runoff (works for all snow models)
            associate(meltflux_inst => domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
                      meltflux_cum  => domain%vars_2d(domain%var_indx(kVARS%meltflux_out_cumul)%v)%data_2d)
                !$acc kernels default(present)
                meltflux_cum(its:ite,jts:jte) = meltflux_cum(its:ite,jts:jte) + meltflux_inst(its:ite,jts:jte)
                meltflux_inst(its:ite,jts:jte) = 0.0 ! reset the instantaneous melt flux for the next LSM step
                !$acc end kernels
            end associate

            if (options%sm%suspension_layer == 1) then
                ! zero out bs_swe_exch since it tracked mass changes
                ! across this LSM step only.
                associate (bs_swe_exch  => domain%vars_2d(domain%var_indx(kVARS%bs_swe_exchange)%v)%data_2d, &
                           bs_swe_emax  => domain%vars_2d(domain%var_indx(kVARS%bs_swe_erode_max)%v)%data_2d)
                !$acc parallel loop collapse(2) present(bs_swe_exch, bs_swe_emax)
                do j= jts, jte
                    do i= its, ite
                        bs_swe_exch(i,j) = 0.0
                        bs_swe_emax(i,j) = 0.0
                    end do
                end do
                end associate
            endif
        endif

        end subroutine snow_model


    !>------------------------------------------------------------
    !! Run NoahMP's soil-temperature pipeline only, on snow-covered cells,
    !! using SNOWPACK's basal conductive heat flux (G_base, computed by
    !! snowpack_step) as the upper Neumann BC.
    !!
    !! Why: when SNOWPACK owns the snow column, NoahmpHICARmain runs with
    !! snow state spoofed to zero and its soil-T output is discarded for
    !! snow cells (the line-1282 gate in lsm_driver.F90). That left
    !! soil_temperature on snow-covered cells frozen in time. This helper
    !! evolves soil_temperature on those cells by calling NoahMP's
    !! GroundThermalProperty + SoilSnowTemperatureMain only — no surface
    !! energy balance, no canopy, no snow physics, no soil moisture solve.
    !!
    !! Coupling is staggered: SNOWPACK uses last-step soil_temperature(:,1,:)
    !! as Dirichlet at its base; this helper takes SNOWPACK's just-computed
    !! G_base as Neumann at the top of the soil column.
    !!
    !! For non-snow cells the soil-T pipeline is also re-run as a side
    !! effect (NoahMP routines parallelize internally over the full domain),
    !! but the result is not written back to HICAR (gated on snow_height>0)
    !! and noahmp%* state is overwritten by the next NoahmpHICARmain call's
    !! *VarInTransfer step.
    !!------------------------------------------------------------
    subroutine noahmp_soil_only(domain, options, dt, NIO)
        implicit none
        type(domain_t),       intent(inout) :: domain
        type(options_t),      intent(in)    :: options
        real,                 intent(in)    :: dt
        type(NoahmpIO_type),  intent(inout) :: NIO

        integer :: i, j, k
        integer :: nsoil, nsnow_max
        real, pointer, dimension(:,:)   :: snow_height, G_base, melt_basal
        real, pointer, dimension(:,:)   :: runoff_surface, runoff_subsurface
        real, pointer, dimension(:,:,:) :: soil_temperature, soil_moisture, soil_liq

        associate( &
            snow_height       => domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d, &
            G_base            => domain%vars_2d(domain%var_indx(kVARS%snow_basal_heat_flux)%v)%data_2d, &
            melt_basal        => domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d, &
            runoff_surface    => domain%vars_2d(domain%var_indx(kVARS%runoff_surface)%v)%data_2d, &
            runoff_subsurface => domain%vars_2d(domain%var_indx(kVARS%runoff_subsurface)%v)%data_2d, &
            soil_temperature  => domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d, &
            soil_moisture     => domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d, &
            soil_liq          => domain%vars_3d(domain%var_indx(kVARS%soil_water_content_liq)%v)%data_3d)

        nsoil     = noahmp%config%domain%NumSoilLayer
        nsnow_max = noahmp%config%domain%NumSnowLayerMax

        ! 1. Refresh NoahmpIO with HICAR's truth on snow cells AND on cells
        !    where SNOWPACK just produced meltwater (snow_height now 0 but
        !    meltflux > 0 — the meltout-step transition). Without that
        !    second condition, the meltout pulse would be dropped because
        !    the helper would skip the cell entirely.
        !    NoahmpIO currently holds NoahmpHICARmain's spoofed-no-snow
        !    output for these cells; we overwrite that with HICAR truth.
        !$acc parallel loop collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                if (snow_height(i,j) <= 0.0 .and. melt_basal(i,j) <= 0.0) cycle
                !$acc loop seq
                do k = 1, nsoil
                    NIO%TSLB (i,k,j) = soil_temperature(i,k,j)
                    NIO%SMOIS(i,k,j) = soil_moisture   (i,k,j)
                    NIO%SH2O (i,k,j) = soil_liq        (i,k,j)
                end do
                NIO%ISNOWXY(i,j) = 0
                !$acc loop seq
                do k = -nsnow_max+1, 0
                    NIO%TSNOXY (i,k,j) = 0.0
                    NIO%ZSNSOXY(i,k,j) = 0.0
                    NIO%SNICEXY(i,k,j) = 0.0
                    NIO%SNLIQXY(i,k,j) = 0.0
                end do
            end do
        end do

        ! 2. Marshal NoahmpIO -> module-level noahmp (whole-domain transfer;
        !    the routines wrap their own !$acc parallel loop internally).
        !    Config first so timestep, soil-layer thicknesses, and
        !    NumSnowLayerNeg are all sync'd to the snow-cell-corrected
        !    NoahmpIO state. Forcing is needed for TemperatureSoilBottom
        !    (the deep-ground Dirichlet temperature read from NoahmpIO%TMN)
        !    so we don't depend on whatever value the prior NoahmpHICARmain
        !    call happened to leave behind.
        call ConfigVarInTransfer (noahmp, NIO)
        call ForcingVarInTransfer(noahmp, NIO)
        call EnergyVarInTransfer (noahmp, NIO)
        call WaterVarInTransfer  (noahmp, NIO)

        ! 3. Override Neumann BC on snow cells: G_base from SNOWPACK and
        !    zero SW penetration (SNOWPACK absorbs all SW above).
        !    NumSnowLayerNeg was already set to 0 by ConfigVarInTransfer
        !    via NIO%ISNOWXY (which we set to 0 for snow cells in step 1).
        !    On a meltout step (snow_height now 0 but melt_basal > 0)
        !    G_base was zeroed by snowpack_driver's bare-ground branch, so
        !    HeatGroundTotMean = 0 here is correct: no spurious heat flux
        !    is injected for the part of the step after the snow vanished.
        !$acc parallel loop collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                if (snow_height(i,j) <= 0.0 .and. melt_basal(i,j) <= 0.0) cycle
                noahmp%energy%flux%HeatGroundTotMean(i,j)   = G_base(i,j)
                noahmp%energy%flux%RadSwPenetrateGrd(i,:,j) = 0.0
            end do
        end do

        ! 4. Soil-only physics: refresh thermal properties, solve soil-T
        !    tridiagonal system, then handle ice<->liquid phase change in
        !    the soil column. With NumSnowLayerNeg = 0, every snow loop
        !    inside SoilSnowWaterPhaseChange has bound (1,0) and skips —
        !    only soil layers see updates. All routines parallelize
        !    internally over the full (I,J) domain.
        call GroundThermalProperty   (noahmp)
        call SoilSnowTemperatureMain (noahmp)
        call SoilSnowWaterPhaseChange(noahmp)

        ! 5. Soil-water solve. Drive infiltration on snow cells from
        !    SNOWPACK's basal melt flux (kg/m^2 per LSM step ->
        !    SoilSfcInflowMean in m/s). Surface evap and transpiration are
        !    zero under snow (SNOWPACK owns surface energy; no veg).
        !
        !    SoilIce must be refreshed after phase change because
        !    SoilSnowWaterPhaseChange modifies SoilLiqWater & SoilMoisture
        !    but leaves SoilIce stale. NoahMP's normal pipeline does the
        !    same refresh at WaterMainMod.F90:114 before calling SoilWaterMain.
        !
        !    For non-snow cells the call still runs (SoilWaterMain has no
        !    per-cell gate other than IndicatorIceSfc==-1) but its result is
        !    ignored: their writeback is gated on snow_height>0, and noahmp
        !    state is reset by next step's *VarInTransfer.
        !$acc parallel loop collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                noahmp%water%flux%SoilSfcInflowMean (i,j) = 0.0
                noahmp%water%flux%EvapSoilSfcLiqMean(i,j) = 0.0
                !$acc loop seq
                do k = 1, nsoil
                    noahmp%water%state%SoilIce(i,k,j) = max(0.0, &
                        noahmp%water%state%SoilMoisture (i,k,j) - &
                        noahmp%water%state%SoilLiqWater(i,k,j))
                    noahmp%water%flux%TranspWatLossSoilMean(i,k,j) = 0.0
                end do
                ! Inject SNOWPACK basal melt as soil-surface inflow on
                ! snow cells AND on meltout-step cells (snow_height now 0
                ! but melt_basal > 0 — last basal-water pulse from the
                ! disappearing snowpack).
                if (snow_height(i,j) > 0.0 .or. melt_basal(i,j) > 0.0) then
                    ! kg/m^2 per LSM step -> m/s: divide by rho_water and dt
                    noahmp%water%flux%SoilSfcInflowMean(i,j) = melt_basal(i,j) / (1000.0 * dt)
                end if
            end do
        end do

        call SoilWaterMain(noahmp)

        ! 6. Write updated soil state back to HICAR for snow cells, reading
        !    DIRECTLY from noahmp%energy%state and noahmp%water%state.
        !
        !    We deliberately DO NOT call Energy/WaterVarOutTransfer here:
        !    every NoahMP *VarOutTransfer routine starts with
        !        if (NoahmpIO%XLAND(I,J) - 1.5 >= 0.0) cycle
        !    (e.g. EnergyVarOutTransferMod.F90:48, WaterVarOutTransferMod.F90:39),
        !    and lsm_driver.F90:1009 marks snow cells as XLAND = kLC_WATER (=2).
        !    So those out-transfers would skip the cells we care about, leaving
        !    NoahmpIO%TSLB / SMOIS for snow cells equal to whatever step 1
        !    seeded — i.e. HICAR's pre-step values — and silently discarding
        !    the soil-T and phase-change result. Reading noahmp%* directly
        !    bypasses that gate.
        !
        !    Side benefit: NoahmpIO%SMOIS / SH2O for non-snow cells stays at
        !    NoahmpHICARmain's authoritative first-call result, so there is
        !    no cumulative drift across LSM steps (soil moisture has no
        !    HICAR<->NoahmpIO round-trip equivalent to the soil-T resync at
        !    lsm_driver.F90:999).
        !
        !    The noahmp%* state for non-snow cells was perturbed by our
        !    second-call solve, but it is reset by the next LSM step's
        !    *VarInTransfer (which reads from the still-clean NoahmpIO).
        !$acc parallel loop collapse(2) default(present)
        do j = jts, jte
            do i = its, ite
                ! Same gate as steps 1, 3, 5 — snow cells AND meltout-step
                ! cells. The meltout step is the LAST step where SNOWPACK
                ! owns the cell; from the next step on, snow_height is 0
                ! AND melt_basal is 0, control returns to NoahmpHICARmain.
                if (snow_height(i,j) <= 0.0 .and. melt_basal(i,j) <= 0.0) cycle
                !$acc loop seq
                do k = 1, nsoil
                    soil_temperature(i,k,j) = noahmp%energy%state%TemperatureSoilSnow(i,k,j)
                    soil_moisture   (i,k,j) = noahmp%water%state%SoilMoisture       (i,k,j)
                    soil_liq        (i,k,j) = noahmp%water%state%SoilLiqWater       (i,k,j)
                end do
                ! NoahmpHICARmain populates runoff_surface/runoff_subsurface for
                ! non-snow cells from NoahmpIO%RUNSFXY/RUNSBXY (which the XLAND
                ! gate left stale for snow cells). Refresh from the live
                ! noahmp%water%flux values so snow cells aren't frozen at junk.
                ! Units match NoahmpHICARmain's existing wiring: mm per soil
                ! timestep (= mm per LSM step when NumSoilTimeStep = 1).
                runoff_surface   (i,j) = noahmp%water%flux%RunoffSurface   (i,j)
                runoff_subsurface(i,j) = noahmp%water%flux%RunoffSubsurface(i,j)
            end do
        end do
        end associate
        
    end subroutine noahmp_soil_only

    end module snow_model_driver