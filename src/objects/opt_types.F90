module options_types

    use icar_constants,             only : kMAX_STRING_LENGTH, MAXLEVELS, MAXFILELENGTH, MAX_NUMBER_FILES, MAXVARLENGTH, kMAX_STORAGE_VARS, kMAX_NAME_LENGTH
    use time_object,                only : Time_type
    use time_delta_object,          only : time_delta_t

    implicit none


    ! ! ! !  TEST for AGL branch

    ! ------------------------------------------------
    ! type to store integer options for each physics package
    ! ------------------------------------------------
    type physics_type
        integer::microphysics
        integer::advection
        integer::boundarylayer
        integer::landsurface
        integer::surfacelayer
        integer::snowmodel
        integer::watersurface
        integer::radiation
        integer::convection
        integer::windtype
        integer::radiation_downScaling        
    end type physics_type
    
    
    ! ------------------------------------------------
    ! store wind solver and parameterization options
    ! ------------------------------------------------
    type wind_type
        logical :: Sx
        logical :: thermal
        real    :: TPI_scale
        real    :: TPI_dmax
        real    :: Sx_dmax
        real    :: Sx_scale_ang
        real    :: alpha_const
        integer :: roughness
        logical :: wind_only
        type(time_delta_t) :: update_dt  ! how often winds are updated
        integer :: wind_iterations      ! number of time to iterate for wind=3 option
        real :: smooth_wind_distance    ! distance over which to smooth the forcing wind field (m)

    end type wind_type
    
    type time_options_type
        logical :: RK3  !Logical of whether ot use RK3 time stepping for advection terms
        real :: cfl_reduction_factor    ! amount to multiple CFL by to improve stability (typically 1)
    end type time_options_type
    
    ! ------------------------------------------------
    ! store Microphysics sensitivity options
    ! ------------------------------------------------
    type mp_options_type
        real    :: Nt_c
        real    :: TNO
        real    :: am_s
        real    :: rho_g
        real    :: av_s, bv_s, fv_s, av_i
        real    :: av_g, bv_g
        real    :: Ef_si, Ef_rs, Ef_rg, Ef_ri
        real    :: C_cubes, C_sqrd
        real    :: mu_r
        real    :: t_adjust
        logical :: Ef_rw_l, EF_sw_l

        real :: update_interval  ! maximum number of seconds between updates
        integer :: top_mp_level     ! top model level to process in the microphysics
    end type mp_options_type

    ! ------------------------------------------------
    ! store Linear Theory options
    ! ------------------------------------------------
    type lt_options_type
        integer :: buffer                   ! number of grid cells to buffer around the domain MUST be >=1
        integer :: stability_window_size    ! window to average nsq over
        real    :: max_stability            ! limits on the calculated Brunt Vaisala Frequency
        real    :: min_stability            ! these may need to be a little narrower.
        logical :: variable_N               ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: smooth_nsq               ! Smooth the Calculated N^2 over vert_smooth vertical levels
        integer :: vert_smooth              ! number of model levels to smooth winds over in the vertical

        real    :: N_squared                ! static Brunt Vaisala Frequency (N^2) to use
        real    :: linear_contribution      ! fractional contribution of linear perturbation to wind field (e.g. u_hat multiplied by this)
        logical :: remove_lowres_linear     ! attempt to remove the linear mountain wave from the forcing low res model
        real    :: rm_N_squared             ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        real    :: rm_linear_contribution   ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        real    :: linear_update_fraction   ! fraction of linear perturbation to add each time step
        logical :: spatial_linear_fields    ! use a spatially varying linear wind perturbation
        logical :: linear_mask              ! use a spatial mask for the linear wind field
        logical :: nsq_calibration          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! Look up table generation parameters
        real    :: dirmax, dirmin           ! minimum and maximum directions to use in the LUT (typically 0 and 2*pi)
        real    :: spdmax, spdmin           ! minimum and maximum wind speeds to use in the LU (typically 0 and ~30)
        real    :: nsqmax, nsqmin           ! minimum and maximum brunt_vaisalla frequencies (typically ~1e-8 and 1e-3)
        integer :: n_dir_values, n_nsq_values, n_spd_values ! number of LUT bins for each parameter
        real    :: minimum_layer_size       ! Minimum vertical step to permit when computing LUT.
                                            ! If model layers are thicker, substepping will be used.

        logical :: read_LUT, write_LUT      ! options to read the LUT from disk (or write it)
        character(len=MAXFILELENGTH) :: u_LUT_Filename  ! u LUT filename to write
        character(len=MAXFILELENGTH) :: v_LUT_Filename  ! v LUT filename to write
        logical :: overwrite_lt_lut         ! if true any existing LUT file will be over written

    end type lt_options_type

    ! ------------------------------------------------
    ! store Advection options
    ! ------------------------------------------------
    type adv_options_type
        logical :: boundary_buffer          ! buffer to smooth a bit to cut down on oscillations at the border if FCT is not used
        logical :: MPDATA_FCT               ! use Flux Corrected Transport (FCT) to maintain stability and prevent any wild oscllations
        integer :: mpdata_order             ! accuracy order for MP_DATA advection scheme.
        integer :: flux_corr                ! Designates which flux-correction scheme to use
        integer :: h_order                  ! Designates which order the horizontal advection should be
        integer :: v_order                  ! Designates which order the vertical advection should be
        logical :: advect_density       ! properly incorporate density into the advection calculations.
                                        ! Doesn't play nice with linear winds

    end type adv_options_type


    ! ------------------------------------------------
    ! store Convection parameter options
    ! ------------------------------------------------
    type cu_options_type
        real :: stochastic_cu
        real :: tendency_fraction
        real :: tend_qv_fraction
        real :: tend_qc_fraction
        real :: tend_th_fraction
        real :: tend_qi_fraction
    end type cu_options_type

    ! ------------------------------------------------
    ! store PBL parameter options
    ! ------------------------------------------------
    type pbl_options_type
        integer :: ysu_topdown_pblmix
    end type pbl_options_type

    ! ------------------------------------------------
    ! store sfc parameter options
    ! ------------------------------------------------
    type sfc_options_type
        integer :: isfflx
        integer :: scm_force_flux
        integer :: iz0tlnd
        integer :: isftcflx
        real    :: sbrlim
    end type sfc_options_type

    ! ------------------------------------------------
    ! store Land Surface Model options
    ! ------------------------------------------------
    type lsm_options_type
        character (len=MAXVARLENGTH) :: LU_Categories   ! land use categories to read from VEGPARM.tbl (e.g. "USGS")
        real :: max_swe                                 ! maximum value for Snow water equivalent (excess above this is removed)
        real :: snow_den_const                          ! variable for converting snow height into SWE or visa versa when input data is incomplete 
        integer :: fsm_nsnow_max                        ! maximum number of snow layers for FSM2 to use. Set here, since it will
                                                        ! change the size of arrays elsewhere in domain_obj        
        real :: update_interval                         ! minimum time to let pass before recomputing LSM ~300s (it may be longer)  [s]
        ! the following categories will be set by default if an known LU_Category is used
        integer :: urban_category                       ! LU index value that equals "urban"
        integer :: ice_category
        integer :: water_category
        integer :: lake_category
        ! integer :: snow_category ! = ice cat
        ! use monthly vegetation fraction data, not just a single value
        logical :: monthly_vegfrac
        logical :: monthly_albedo
        integer :: sf_urban_phys
        
        integer :: nmp_dveg
        integer :: nmp_opt_crs
        integer :: nmp_opt_sfc
        integer :: nmp_opt_btr
        integer :: nmp_opt_run
        integer :: nmp_opt_infdv
        integer :: nmp_opt_frz
        integer :: nmp_opt_inf
        integer :: nmp_opt_rad
        integer :: nmp_opt_alb
        integer :: nmp_opt_snf
        integer :: nmp_opt_tbot
        integer :: nmp_opt_stc
        integer :: nmp_opt_gla
        integer :: nmp_opt_rsf
        integer :: nmp_opt_soil
        integer :: nmp_opt_pedo
        integer :: nmp_opt_crop
        integer :: nmp_opt_irr
        integer :: nmp_opt_irrm
        integer :: nmp_opt_tdrn
        integer :: noahmp_output
        real    :: nmp_soiltstep

    end type lsm_options_type

    ! ------------------------------------------------
    ! store Radiation options
    ! ------------------------------------------------
    type rad_options_type
       real    :: update_interval_rrtmg                ! how ofen to update the radiation in seconds.
                                                       ! RRTMG scheme is expensive. Default is 1800s (30 minutes)
       integer :: icloud                               ! How RRTMG interact with clouds
       integer :: cldovrlp                             ! how RRTMG considers cloud overlapping (1 = random, 2 = maximum-random, 3 = maximum, 4 = exponential, 5 = exponential-random)
       logical :: read_ghg                             ! Eihter use default green house gas mixing ratio, or read the in from file
       real    :: tzone !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe
    end type rad_options_type

    ! ------------------------------------------------
    ! store output related options
    ! ------------------------------------------------
    type output_options_type

        ! file names
        character (len=MAXFILELENGTH) :: output_file        

        real :: outputinterval          ! time between output [s]
        integer :: frames_per_outfile      ! frames (outputintervals) per out file

        ! these are the variables that need to be written and read from disk as primary output
        integer :: vars_for_output( kMAX_STORAGE_VARS ) = 0

        ! The following are NOT read from the namelist -- instead they are set/calculated from other options which ARE read from the namelist
        type(time_delta_t) :: output_dt ! store out_dt as a time delta object

    end type output_options_type

    ! ------------------------------------------------
    ! store restart related options
    ! ------------------------------------------------
    type restart_options_type

        logical :: restart              ! this is a restart run, read model conditions from a restart file
        integer :: restart_count        ! the frequency, in number of output timesteps, that we write out a restart file

        ! file names
        character (len=MAXFILELENGTH) :: restart_out_file, restart_in_file

        ! restart information
        type(Time_type) :: restart_time ! Date of the restart time step
        integer :: restart_step_in_file ! step in restart file to initialize from

    end type restart_options_type


    ! ------------------------------------------------
    ! store general simulation options
    ! ------------------------------------------------
    type general_options_type

        character (len=MAXVARLENGTH) :: version, comment, phys_suite

        logical :: debug                ! outputs a little more information at runtime (not much at present)
        logical :: interactive          ! set to true if running at the commandline to see %complete printed
        logical :: ideal                ! this is an ideal simulation, forcing will be held constant

        ! date/time parameters
        type(Time_type) :: start_time   ! Date to start running the model
        type(Time_type) :: end_time     ! End point for the model simulation
        character(len=MAXFILELENGTH) :: calendar


        ! physics parameterization options
        logical :: use_mp_options
        logical :: use_cu_options
        logical :: use_lt_options
        logical :: use_adv_options
        logical :: use_lsm_options
        logical :: use_rad_options
        logical :: use_pbl_options
        logical :: use_sfc_options
        logical :: use_wind_options

    end type general_options_type


    ! ------------------------------------------------
    ! store forcing related options
    ! ------------------------------------------------
    type forcing_options_type

        logical :: qv_is_relative_humidity! if true the input water vapor is assumed to be relative humidity instead of mixing ratio
        logical :: qv_is_spec_humidity  ! if true the input water vapor is assumed to be specific humidity instead of mixing ratio
        logical :: t_is_potential       ! if true the input temperature is interpreted as potential temperature
        logical :: z_is_geopotential    ! if true the z variable is interpreted as geopotential height
        logical :: z_is_on_interface    ! if true the z variable is interpreted as residing at model level interfaces
        logical :: time_varying_z       ! read in a new z coordinate every time step and interpolate accordingly

        real :: t_offset                ! offset to temperature because WRF outputs potential temperature-300
        logical :: limit_rh                ! impose a limit on relative humidity in the forcing data to be <=1

        real :: inputinterval           ! time between forcing steps [s]


        ! variable names from init/BC/wind/... files
        character (len=MAXVARLENGTH) :: latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,wvar, &
                                        pvar,tvar,qvvar,qcvar,qivar,qrvar,qsvar,qgvar,i2mvar,i3mvar,&
                                        qncvar,qnivar,qnrvar,qnsvar,qngvar,i2nvar,i3nvar,&
                                        i1avar,i1cvar,i2avar,i2cvar,i3avar,i3cvar,hgtvar, &
                                        pslvar, psvar, sst_var, pblhvar, &
                                        shvar,lhvar,zvar, &
                                        swdown_var, lwdown_var, &
                                        time_var

        ! The following are NOT read from the namelist -- instead they are set/calculated from other options which ARE read from the namelist
        character(len=MAXVARLENGTH) :: vars_to_read(kMAX_STORAGE_VARS)
        integer                     :: dim_list(    kMAX_STORAGE_VARS)

        character (len=MAXFILELENGTH), dimension(:), allocatable :: boundary_files
        type(time_delta_t) :: input_dt  ! store in_dt as a time delta object
        logical :: compute_z            ! flag that we need to compute z from p, this is determined from the vars specified (not read)

    end type forcing_options_type


    ! ------------------------------------------------
    ! store all model options
    ! ------------------------------------------------
    type domain_options_type

        ! file names
        character (len=MAXFILELENGTH) :: init_conditions_file                                        

        
        ! various real parameters/options
        real :: dx                      ! grid cell width [m]

        integer :: longitude_system     ! specify center for longitude system
                                        ! 0 = kPRIME_CENTERED    (-180 - 180)
                                        ! 1 = kDATELINE_CENTERED (0 - 360)        

        integer :: nz                   ! number of model vertical levels

        ! note this can't be allocatable because gfortran does not support allocatable components inside derived type coarrays...
        real, dimension(MAXLEVELS)::dz_levels ! model layer thicknesses to be read from namelist
        real    :: flat_z_height        ! height above mean ground level [m] above which z levels are flat in space
        
        logical :: sleve                ! Using a sleve space_varying_dz offers control over the decay of terrain features in the vertical grid structure. See Schär et al 2002, Leuenberger et al 2009
        integer :: terrain_smooth_windowsize
        integer :: terrain_smooth_cycles
        real    :: decay_rate_L_topo    !
        real    :: decay_rate_S_topo    !
        real    :: sleve_n              ! Additional parameter introduced by Leuenberger 2009.


        logical :: use_agl_height       ! interpolate from forcing to model layers using Z above ground level, not sea level
        real    :: agl_cap              ! height up to which AGL height is used for vertical interpolation

        ! variable names from init/BC/wind/... files
        character (len=MAXVARLENGTH) :: hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,landvar,lakedepthvar, &
                                        snowh_var, soiltype_var, soil_t_var,soil_vwc_var,swe_var,soil_deept_var, &
                                        vegtype_var,vegfrac_var, albedo_var, vegfracmax_var, lai_var, canwat_var, &
                                        linear_mask_var, nsq_calibration_var, &
                                        sinalpha_var, cosalpha_var

        character(len=MAXVARLENGTH) :: svf_var, hlm_var, slope_var, slope_angle_var, aspect_angle_var, shd_var !!MJ added


    end type domain_options_type

end module options_types
