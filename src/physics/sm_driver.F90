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
    use module_sm_SNOWPACKdrv, only : sm_snowpack_init, sm_SNOWPACK
#endif

    implicit none

    private
    public :: sm_var_request, sm_init, snow_model   !, sm_finalize

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: restart

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
                         kVARS%lsm_last_precip, kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity,    &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%snow_nlayers,   &
                         kVARS%shd, kVARS%meltflux_out_tstep, kVARS%Sliq_out, &
                         kVARS%windspd_10m, kVARS%dSWE_salt, kVARS%dSWE_susp, kVARS%dSWE_blow_subl, kVARS%dSWE_subl, kVARS%dSWE_slide, &
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
                         kVARS%runoff_tstep, kVARS%snow_temperature, kVARS%Sice, kVARS%Sliq, kVARS%Ds, kVARS%fsnow, kVARS%snow_nlayers  ])

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
                         kVARS%lsm_last_precip, kVARS%QFX, kVARS%chs, kVARS%chs2, kVARS%cqs2, kVARS%land_emissivity,    &
                         kVARS%veg_type, kVARS%soil_type, kVARS%land_mask, kVARS%snowfall, kVARS%albedo,                &
                         kVARS%snow_temperature, kVARS%Ds, kVARS%snow_temperature_i, kVARS%Vol_Frac_I, kVARS%Vol_Frac_W, &
                         kVARS%Vol_Frac_A, kVARS%Vol_Frac_S, kVARS%Vol_Frac_WP, kVARS%Rg, kVARS%Rb, kVARS%Dd, kVARS%Sp, kVARS%mk, &
                         kVARS%snow_nlayers,                                                                      &
                         kVARS%mass_hoar, kVARS%CDot, kVARS%metamo, kVARS%depositionDate, kVARS%dSWE_subl])
             
             call options%restart_vars( &
                         [kVARS%sst, kVARS%water_vapor, kVARS%potential_temperature, kVARS%precipitation, kVARS%temperature, &
                         kVARS%density, kVARS%pressure_interface, kVARS%shortwave,                                      &
                         kVARS%longwave, kVARS%canopy_water, kVARS%snow_water_equivalent,                               &
                         kVARS%skin_temperature, kVARS%soil_water_content, kVARS%soil_temperature, kVARS%terrain,       &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m, kVARS%temperature_2m,        &
                         kVARS%snow_height,  kVARS%snowfall, kVARS%albedo, kVARS%QFX, kVARS%land_emissivity,            &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%longwave_up, kVARS%ground_heat_flux,          &
                         kVARS%soil_totalmoisture, kVARS%roughness_z0, kVARS%lsm_last_snow, kVARS%lsm_last_precip,      &
                         kVARS%snow_temperature, kVARS%Ds, kVARS%snow_temperature_i, kVARS%Vol_Frac_I, kVARS%Vol_Frac_W, &
                         kVARS%Vol_Frac_A, kVARS%Vol_Frac_S, kVARS%Vol_Frac_WP, kVARS%Rg, kVARS%Rb, kVARS%Dd, kVARS%Sp, kVARS%mk, &
                         kVARS%snow_nlayers,                                                                      &
                         kVARS%mass_hoar, kVARS%CDot, kVARS%metamo, kVARS%depositionDate])
        endif
    end subroutine sm_var_request


    subroutine sm_init(domain,options,context_chng)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        logical, optional,  intent(in)      :: context_chng

        logical :: context_change

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
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    SnowModel: Snowpack"
            call sm_snowpack_init(domain,options,context_change)
#else
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    User asked to use Snowpack, but it is not compiled in this version of HICAR"
            if (STD_OUT_PE .and. .not.context_change) write(*,*) "    Please de-select Snowpack as the snow model in the namelist, or recompile HICAR with the Snowpack library linked"
            stop "Snowpack not compiled in this version of HICAR"
#endif

        endif

    end subroutine sm_init

    subroutine snow_model(domain,options,lsm_dt)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options
        real,               intent(in)      :: lsm_dt

        real, allocatable, dimension(:,:) :: windspd
        real, allocatable, dimension(:,:) :: current_precipitation
        real, allocatable, dimension(:,:) :: current_snow
        real, allocatable, dimension(:,:) :: current_rain

        if (options%physics%snowmodel == 0) return


        allocate(current_precipitation(ims:ime,jms:jme), source=0.0)
        allocate(windspd(ims:ime,jms:jme), source=1.0)
        allocate(current_snow(ims:ime,jms:jme), source=0.0) ! MJ added 
        allocate(current_rain(ims:ime,jms:jme), source=0.0) ! MJ added 

        associate( &
            u_10m => domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
            v_10m => domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
            precipitation => domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
            snowfall => domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d, &
            lsm_last_precip => domain%vars_2d(domain%var_indx(kVARS%lsm_last_precip)%v)%data_2d, &
            lsm_last_snow => domain%vars_2d(domain%var_indx(kVARS%lsm_last_snow)%v)%data_2d &
            )

        ! --- Update host: common snow model inputs ---
        !$acc update host( &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%u_10m)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%v_10m)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%precipitation)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%snowfall)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%lsm_last_precip)%v)%data_2d, &
        !$acc   domain%vars_2d(domain%var_indx(kVARS%lsm_last_snow)%v)%data_2d)

        windspd = sqrt(u_10m**2 + v_10m**2)
        where(windspd<1) windspd=1 ! minimum wind speed to prevent the exchange coefficient from blowing up
        current_precipitation = (precipitation - lsm_last_precip) !+(domain%precipitation_bucket-rain_bucket)*kPRECIP_BUCKET_SIZE

        ! Setup the input data for the snow models
        current_snow = (snowfall-lsm_last_snow) !! MJ: snowfall in kg m-2
        current_rain = max(current_precipitation-current_snow,0.) !! MJ: rainfall in kg m-2

        if (options%physics%snowmodel==kSM_FSM) then
#ifdef FSM
            ! --- Update host: FSM-specific inputs ---
            !$acc update host( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
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
            !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%snow_nlayers)%v)%data_2di, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di, &
            !$acc   domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d, &
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
            !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_salt)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_susp)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%dSWE_blow_subl)%v)%data_2d, &
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
            ! EVAN'S BIG DRIVER HERE:
            ! input variables for SNOWPACK, which we need, are:

            ! - ta;       ///< Air temperature (K)
            ! - rh;       ///< Relative humidity (% or 1)
            ! - vw;       ///< Wind velocity at snow station (m s-1)
            !!! - ustar;    ///< The friction velocity (m s-1) computed in mt_MicroMet() and also used later for the MeteoHeat fluxes
            ! - z0;       ///< The roughness length computed in SnowDrift and also used later for the MeteoHeat fluxes (m)
            ! - psi_s;    ///< Stability correction for scalar heat fluxes
            ! - iswr;     ///< Incoming SHORTWAVE radiation (W m-2)
            !!! - ea;       ///< Atmospheric emissivity (1)
            ! - lw_net;   ///< Net longwave radiation (W m-2)
            ! - tss;      ///< Snow surface temperature (K)
            ! - ts0;      ///< Bottom temperatures of snow/soil pack (K)
            ! - psum;     ///< precipitation sum over the current timestep (mm)
            ! - psum_ph;  ///< precipitation phase for the current timestep (between 0 and 1, 0 is fully solid while 1 is fully liquid).
            !!! - hs;       ///< The measured height of snow (m)
            ! - geo_heat; ///< Geo heat flux (W/m^2), for the neumann lower boundary condition in the heat equation

            ! --- Update host: SNOWPACK-specific inputs ---
            !$acc update host( &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d, &
            !$acc   domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d, &
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
            !$acc   domain%vars_3d(domain%var_indx(kVARS%metamo)%v)%data_3d)

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
            !$acc   domain%vars_3d(domain%var_indx(kVARS%metamo)%v)%data_3d)
#endif
        endif

        end associate

        end subroutine snow_model

    end module snow_model_driver