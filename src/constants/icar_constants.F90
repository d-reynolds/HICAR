!>------------------------------------------------
!! Defines model constants (e.g. gravity, and kMAX_FILE_LENGTH)
!!
!!------------------------------------------------
module icar_constants

    implicit none

    character(len=5) :: kVERSION_STRING = "v2.1"
    
    ! Define process info and team info
    integer, parameter :: kCOMPUTE_TEAM = 1
    integer, parameter :: kIO_TEAM = 2
    logical :: STD_OUT_PE = .False.
    logical :: STD_OUT_PE_IO = .False.
    
    character(len=5), parameter :: kCHAR_NO_VAL = "-9999"
    integer,          parameter :: kINT_NO_VAL = -123456
    real,             parameter :: kREAL_NO_VAL = -123456.0
    real,             parameter :: kUNSET_REAL = -987654.0
    !Flag-value to indicate a part of a read-write buffer which was never filled
    real, parameter :: kEMPT_BUFF = -123456789.0
    
    ! string lengths
    integer, parameter :: kMAX_FILE_LENGTH =   1024
    integer, parameter :: kMAX_DIM_LENGTH  =   256
    integer, parameter :: kMAX_NAME_LENGTH =   256
    integer, parameter :: kMAX_ATTR_LENGTH =   256
    integer, parameter :: kMAX_STRING_LENGTH = 256  ! maximum length of other strings (e.g. netcdf attributes)

    ! calendar information
    character(len=20)          :: kDEFAULT_CALENDAR = "GREGORIAN"  ! default calendar type
    integer, parameter, public :: GREGORIAN=0, NOLEAP=1, THREESIXTY=2, NOCALENDAR=-1
    integer, parameter, public :: NON_VALID_YEAR = -9999


    ! maximum number of nests
    integer, parameter :: kMAX_NESTS = 10

    type var_constants_type
        SEQUENCE    ! technically SEQUENCE just requires the compiler leave them in order,
                    ! but it can also keep compilers (e.g. ifort) from padding for alignment,
                    ! as long as there is no padding we can test last_var = sizeof(kVARS)

        integer :: u
        integer :: v
        integer :: w
        integer :: w_real
        integer :: pressure
        integer :: pressure_interface
        integer :: pressure_base
        integer :: geopotential_base
        integer :: potential_temperature
        integer :: temperature
        integer :: water_vapor
        integer :: cloud_water_mass
        integer :: cloud_number
        integer :: ice_mass
        integer :: ice_number
        integer :: rain_mass
        integer :: rain_number
        integer :: snow_mass
        integer :: snow_number
        integer :: graupel_mass
        integer :: graupel_number
        integer :: ice1_a
        integer :: ice1_c
        integer :: ice2_mass
        integer :: ice2_number
        integer :: ice2_a
        integer :: ice2_c
        integer :: ice3_mass
        integer :: ice3_number
        integer :: ice3_a
        integer :: ice3_c
        integer :: precipitation
        integer :: convective_precipitation
        integer :: snowfall
        integer :: graupel
        integer :: snowfall_ground
        integer :: rainfall_ground
        integer :: lsm_last_snow
        integer :: lsm_last_precip
        integer :: exner
        integer :: nsquared
        integer :: density
        integer :: cloud_fraction
        integer :: shortwave
        integer :: shortwave_direct
        integer :: shortwave_diffuse
        integer :: shortwave_direct_above
        integer :: longwave
        integer :: albedo
        integer :: soil_albedo_dir
        integer :: soil_albedo_diff
        integer :: vegetation_fraction
        integer :: vegetation_fraction_max
        integer :: vegetation_fraction_out
        integer :: veg_type
        integer :: mass_leaf
        integer :: mass_root
        integer :: mass_stem
        integer :: mass_wood
        integer :: soil_type
        integer :: soil_texture_1
        integer :: soil_texture_2
        integer :: soil_texture_3
        integer :: soil_texture_4
        integer :: soil_sand_and_clay
        integer :: soil_carbon_stable
        integer :: soil_carbon_fast
        integer :: lai
        integer :: sai
        integer :: crop_category
        integer :: crop_type
        integer :: date_planting
        integer :: date_harvest
        integer :: growing_season_gdd
        integer :: irr_frac_total
        integer :: irr_frac_sprinkler
        integer :: irr_frac_micro
        integer :: irr_frac_flood
        integer :: irr_eventno_sprinkler
        integer :: irr_eventno_micro
        integer :: irr_eventno_flood
        integer :: irr_alloc_sprinkler
        integer :: irr_alloc_micro
        integer :: irr_alloc_flood
        integer :: irr_evap_loss_sprinkler
        integer :: irr_amt_sprinkler
        integer :: irr_amt_micro
        integer :: irr_amt_flood
        integer :: evap_heat_sprinkler
        integer :: mass_ag_grain
        integer :: growing_degree_days
        integer :: plant_growth_stage
        integer :: net_ecosystem_exchange
        integer :: gross_primary_prod
        integer :: net_primary_prod
        integer :: apar
        integer :: photosynthesis_total
        integer :: stomatal_resist_total
        integer :: stomatal_resist_sun
        integer :: stomatal_resist_shade
        integer :: gecros_state
        integer :: canopy_water
        integer :: canopy_water_ice
        integer :: canopy_water_liquid
        integer :: canopy_vapor_pressure
        integer :: canopy_temperature
        integer :: canopy_fwet
        integer :: veg_leaf_temperature
        integer :: ground_surf_temperature
        integer :: frac_between_gap
        integer :: frac_within_gap
        integer :: ground_temperature_bare
        integer :: ground_temperature_canopy
        integer :: sensible_heat
        integer :: latent_heat
        integer :: u_10m
        integer :: v_10m
        integer :: windspd_10m
        integer :: ustar
        integer :: coeff_momentum_drag
        integer :: chs
        integer :: chs2
        integer :: cqs2
        integer :: coeff_heat_exchange_3d
        integer :: coeff_momentum_exchange_3d
        integer :: QFX
        integer :: br
        integer :: mol
        integer :: psim
        integer :: psih
        integer :: fm
        integer :: fh
        integer :: surface_rad_temperature
        integer :: temperature_2m
        integer :: humidity_2m
        integer :: temperature_2m_veg
        integer :: temperature_2m_bare
        integer :: mixing_ratio_2m_veg
        integer :: mixing_ratio_2m_bare
        integer :: surface_pressure
        integer :: sea_surface_pressure
        integer :: rad_absorbed_total
        integer :: rad_absorbed_veg
        integer :: rad_absorbed_bare
        integer :: rad_net_longwave
        integer :: longwave_up
        integer :: ground_heat_flux
        integer :: evap_canopy
        integer :: evap_soil_surface
        integer :: transpiration_rate
        integer :: ch_veg
        integer :: ch_veg_2m
        integer :: ch_bare
        integer :: ch_bare_2m
        integer :: ch_under_canopy
        integer :: ch_leaf
        integer :: sensible_heat_veg
        integer :: sensible_heat_bare
        integer :: sensible_heat_canopy
        integer :: evap_heat_veg
        integer :: evap_heat_bare
        integer :: evap_heat_canopy
        integer :: transpiration_heat
        integer :: ground_heat_veg
        integer :: ground_heat_bare
        integer :: net_longwave_veg
        integer :: net_longwave_bare
        integer :: net_longwave_canopy
        integer :: runoff_surface
        integer :: runoff_subsurface
        integer :: soil_totalmoisture
        integer :: soil_deep_temperature
        integer :: water_table_depth
        integer :: water_aquifer
        integer :: storage_gw
        integer :: storage_lake
        integer :: roughness_z0
        integer :: snow_water_eq_prev
        integer :: snow_albedo_prev
        integer :: snow_layer_depth
        integer :: snow_layer_ice
        integer :: snow_layer_liquid_water
        integer :: snow_age_factor
        integer :: soil_water_content
        integer :: soil_water_content_liq
        integer :: eq_soil_moisture
        integer :: smc_watertable_deep
        integer :: recharge
        integer :: recharge_deep
        integer :: soil_temperature
        integer :: skin_temperature
        integer :: sst
        integer :: tend_qv_adv
        integer :: tend_qv_pbl
        integer :: tend_qv
        integer :: tend_th
        integer :: tend_th_pbl
        integer :: tend_qc
        integer :: tend_qc_pbl
        integer :: tend_qi
        integer :: tend_qi_pbl
        integer :: tend_qs
        integer :: tend_qr
        integer :: tend_u
        integer :: tend_v
        integer :: tend_th_lwrad
        integer :: tend_th_swrad
        integer :: u_mass
        integer :: v_mass
        integer :: re_cloud
        integer :: re_ice
        integer :: re_snow

        !ISHMAEL MP variables
        integer :: ice1_rho
        integer :: ice1_phi
        integer :: ice1_vmi
        integer :: ice2_rho
        integer :: ice2_phi
        integer :: ice2_vmi
        integer :: ice3_rho
        integer :: ice3_phi
        integer :: ice3_vmi

        integer :: wind_alpha
        integer :: froude
        integer :: blk_ri
        integer :: out_longwave_rad
        integer :: longwave_cloud_forcing
        integer :: shortwave_cloud_forcing
        integer :: land_emissivity
        integer :: temperature_interface
        integer :: runoff_tstep


        integer :: kpbl
        integer :: hpbl

        ! Lake model variables
        integer :: lake_depth
        integer :: t_lake3d
        integer :: snl2d
        integer :: t_grnd2d
        integer :: lake_icefrac3d
        integer :: z_lake3d
        integer :: dz_lake3d
        integer :: t_soisno3d
        integer :: h2osoi_ice3d
        integer :: h2osoi_liq3d
        integer :: h2osoi_vol3d
        integer :: z3d
        integer :: dz3d
        integer :: watsat3d
        integer :: csol3d
        integer :: tkmg3d
        integer :: lakemask
        integer :: xice
        integer :: zi3d
        integer :: tksatu3d
        integer :: tkdry3d
        integer :: savedtke12d
        integer :: lakedepth2d


        integer :: ivt
        integer :: iwv
        integer :: iwl
        integer :: iwi

        ! GRID VARIABLES
        integer :: z
        integer :: z_interface
        integer :: dzdx
        integer :: dzdy
        integer :: dz
        integer :: dz_interface
        integer :: advection_dz
        integer :: dzdy_v
        integer :: dzdx_u
        integer :: jacobian
        integer :: jacobian_u
        integer :: jacobian_v
        integer :: jacobian_w
        integer :: land_mask
        integer :: terrain
        integer :: latitude
        integer :: longitude
        integer :: global_terrain
        integer :: global_dz_interface
        integer :: global_z_interface
        integer :: u_latitude
        integer :: u_longitude
        integer :: v_latitude
        integer :: v_longitude
        integer :: Sx
        integer :: TPI
        integer :: neighbor_terrain
        integer :: froude_terrain
        integer :: relax_filter_2d
        integer :: relax_filter_3d
        integer :: costheta
        integer :: sintheta
        integer :: cosine_zenith_angle
        integer :: slope
        integer :: slope_angle
        integer :: aspect_angle
        integer :: svf
        integer :: shd
        integer :: hlm
        integer :: h1
        integer :: h2
        integer :: h1_u
        integer :: h1_v
        integer :: h2_u
        integer :: h2_v

        !NoahMP Variables
        integer :: wetland_sat_frac, wetland_h20_store

        !SNICAR Variables
        integer :: snicar_sn_rad, snicar_sn_fr
        integer :: snicar_bcphi, snicar_bcpho, snicar_ocphi, snicar_ocpho
        integer :: snicar_dust1, snicar_dust2, snicar_dust3, snicar_dust4, snicar_dust5
        integer :: snicar_bcphi_conc, snicar_bcpho_conc, snicar_ocphi_conc, snicar_ocpho_conc
        integer :: snicar_dust1_conc, snicar_dust2_conc, snicar_dust3_conc, snicar_dust4_conc, snicar_dust5_conc

        !General Snow model Variables
        integer :: Sice
        integer :: Sliq
        integer :: Ds
        integer :: fsnow
        integer :: snow_water_equivalent
        integer :: snow_temperature
        integer :: snow_height
        integer :: snow_nlayers

        ! FSM2trans variables
        integer :: dSWE_salt
        integer :: dSWE_susp
        integer :: dSWE_subl
        integer :: dSWE_blow_subl
        integer :: dSWE_slide
        integer :: meltflux_out_tstep
        integer :: Sliq_out

        !SNOWPACK variables
        integer :: depositionDate
        integer :: snow_temperature_i
        integer :: Vol_Frac_I
        integer :: Vol_Frac_W
        integer :: Vol_Frac_A
        integer :: Vol_Frac_S
        integer :: Vol_Frac_WP
        integer :: Rg
        integer :: Rb
        integer :: Dd
        integer :: Sp
        integer :: mk
        integer :: mass_hoar
        integer :: CDot
        integer :: metamo


        integer :: last_var
    end type var_constants_type


    type(var_constants_type) :: kVARS 



    character(len=18) :: one_d_column_dimensions(1)         = [character(len=18) :: "level"]
    character(len=18) :: two_d_dimensions(2)                = [character(len=18) :: "lon_x","lat_y"]
    character(len=18) :: two_d_t_dimensions(3)              = [character(len=18) :: "lon_x","lat_y","time"]
    character(len=18) :: two_d_u_dimensions(2)              = [character(len=18) :: "lon_u","lat_y"]
    character(len=18) :: two_d_v_dimensions(2)              = [character(len=18) :: "lon_x","lat_v"]
    character(len=18) :: two_d_global_dimensions(2)         = [character(len=18) :: "lon_x_global","lat_y_global"]
    character(len=18) :: two_d_neighbor_dimensions(2)       = [character(len=18) :: "lon_x_neighbor","lat_y_neighbor"]
    character(len=18) :: three_d_u_t_dimensions(4)          = [character(len=18) :: "lon_u","level","lat_y","time"]
    character(len=18) :: three_d_v_t_dimensions(4)          = [character(len=18) :: "lon_x","level","lat_v","time"]
    character(len=18) :: three_d_u_dimensions(3)            = [character(len=18) :: "lon_u","level","lat_y"]
    character(len=18) :: three_d_v_dimensions(3)            = [character(len=18) :: "lon_x","level","lat_v"]
    character(len=18) :: three_d_dimensions(3)              = [character(len=18) :: "lon_x","level","lat_y"]
    character(len=18) :: three_d_global_dimensions(3)       = [character(len=18) :: "lon_x_global","level","lat_y_global"]
    character(len=18) :: three_d_neighbor_dimensions(3)     = [character(len=18) :: "lon_x_neighbor","level","lat_y_neighbor"]
    character(len=18) :: three_d_global_interface_dimensions(3)       = [character(len=18) :: "lon_x_global","level_i","lat_y_global"]
    character(len=18) :: three_d_neighbor_interface_dimensions(3)     = [character(len=18) :: "lon_x_neighbor","level_i","lat_y_neighbor"]
    character(len=18) :: three_d_t_dimensions(4)            = [character(len=18) :: "lon_x","level","lat_y","time"]
    character(len=18) :: three_d_interface_dimensions(3)    = [character(len=18) :: "lon_x","level_i","lat_y"]
    character(len=18) :: three_d_t_interface_dimensions(4)  = [character(len=18) :: "lon_x","level_i","lat_y","time"]
    character(len=18) :: three_d_hlm_dimensions(3)          = [character(len=18) :: "lon_x","azimuth","lat_y"]
    character(len=18) :: three_d_t_soil_dimensions(4)       = [character(len=18) :: "lon_x","nsoil","lat_y","time"]
    character(len=18) :: three_d_t_snow_dimensions(4)       = [character(len=18) :: "lon_x","nsnow","lat_y","time"]
    character(len=18) :: three_d_t_snow_i_dimensions(4)     = [character(len=18) :: "lon_x","nsnow_i","lat_y","time"]
    character(len=18) :: three_d_t_snowsoil_dimensions(4)   = [character(len=18) :: "lon_x","nsnowsoil","lat_y","time"]
    character(len=18) :: three_d_soilcomp_dimensions(3)     = [character(len=18) :: "lon_x","nsoil_composition","lat_y"]
    character(len=18) :: three_d_crop_dimensions(3)         = [character(len=18) :: "lon_x","crop","lat_y"]
    character(len=18) :: three_d_t_gecros_dimensions(4)     = [character(len=18) :: "lon_x","gecros","lat_y","time"]
    character(len=18) :: three_d_t_month_dimensions(4)          = [character(len=18) :: "lon_x","month","lat_y","time"]
    character(len=18) :: three_d_t_lake_dimensions(4)           = [character(len=18) :: "lon_x","nlevlake","lat_y","time"]
    character(len=18) :: three_d_t_lake_soisno_dimensions(4)    = [character(len=18) :: "lon_x","nlevsoisno","lat_y","time"] !grid_lake_soisno
    character(len=18) :: three_d_t_lake_soisno_1_dimensions(4)  = [character(len=18) :: "lon_x","nlevsoisno_1","lat_y","time"]
    character(len=18) :: three_d_t_lake_soi_dimensions(4)       = [character(len=18) :: "lon_x","nlevsoi_lake","lat_y","time"] !grid_lake_soi
    character(len=18) :: four_d_azim_dimensions(4)                = [character(len=18) :: "lon_x","level","lat_y","Sx_azimuth"]


    integer, parameter :: kINTEGER_BITS     = storage_size(kINTEGER_BITS)
    integer, parameter :: kMAX_STORAGE_VARS = storage_size(kVARS) / kINTEGER_BITS

    integer, parameter :: kINTEGER          = 1
    integer, parameter :: kREAL             = 4
    integer, parameter :: kDOUBLE           = 8

    ! Initial number of output variables for which pointers are created
    integer, parameter :: kINITIAL_VAR_SIZE= 128

    ! Maximum number of dimensions
    ! Note this is defined in NetCDF, though not enforced (a file can have more than 1024 dimensions)
    integer, parameter :: kMAX_DIMENSIONS  = 1024

!>------------------------------------------------
!! Model constants (mostly string lengths)
!! ------------------------------------------------
    integer, parameter :: MAXLEVELS          =    500  ! maximum number of vertical layers (should typically be ~10-20)
    integer, parameter :: MAX_NUMBER_FILES   =  50000  ! maximum number of permitted input files (probably a bit extreme)
    character(len=54), parameter :: kOUTPUT_FMT = '("days since ",i4,"-",i2.2,"-",i2.2," ",i2.2,":00:00")'

!>------------------------------------------------
!!  Default width of coarray halos, ideally might be physics dependant (e.g. based on advection spatial order)
!! ------------------------------------------------
    integer,parameter :: kDEFAULT_HALO_SIZE = 1

!>------------------------------------------------
!! Value to accept for difference between real numbers should be as a fraction but then have to test for non-zero...
!! For some variables (cloud ice) 1e-6 is not that small, for others (pressure) it might be below precision...
!! ------------------------------------------------
    real,   parameter :: kSMALL_VALUE = 1e-6

    integer, parameter :: kMAINTAIN_LON      = 0
    integer, parameter :: kPRIME_CENTERED    = 1
    integer, parameter :: kDATELINE_CENTERED = 2
    integer, parameter :: kGUESS_LON         = 3

! ------------------------------------------------
! Physics scheme selection definitions
!
! NB: BASIC typically means "use the data from the low res model"
!     SIMPLE typically means a relatively simple formulation written for ICAR
! These could all be switched to enums too, but this makes it easy to see what number each has for the options file...
! ------------------------------------------------
    integer, parameter :: kNO_STOCHASTIC = -9999
    integer, parameter :: kCU_TIEDTKE    = 1
    integer, parameter :: kCU_NSAS       = 2
    integer, parameter :: kCU_BMJ        = 3

    integer, parameter :: kMP_THOMPSON   = 1
    integer, parameter :: kMP_SB04       = 2
    integer, parameter :: kMP_MORRISON   = 3
    integer, parameter :: kMP_WSM6       = 4
    integer, parameter :: kMP_THOMP_AER  = 5
    integer, parameter :: kMP_WSM3       = 6
    integer, parameter :: kMP_ISHMAEL    = 7
 
    integer, parameter :: kPBL_YSU         = 1

    integer, parameter :: kWATER_SIMPLE  = 1
    integer, parameter :: kWATER_LAKE    = 2

    integer, parameter :: kSFC_MM5REV    = 1

    integer, parameter :: kLSM_BASIC     = 1
    integer, parameter :: kLSM_NOAHMP    = 2
    
    integer, parameter :: kSM_FSM        = 1 !! MJ added
    integer, parameter :: kSM_SNOWPACK   = 2

    integer, parameter :: kRA_BASIC      = 1
    integer, parameter :: kRA_SIMPLE     = 2
    integer, parameter :: kRA_RRTMG      = 3
    integer, parameter :: kRA_RRTMGP     = 4
    
    integer, parameter :: kADV_STD       = 1
    integer, parameter :: kADV_MPDATA    = 2

    integer, parameter :: kFLUXCOR_MONO   = 1

    integer, parameter :: kITERATIVE_WINDS = 1

    integer, parameter :: kLC_LAND       = 1
    integer, parameter :: kLC_WATER      = 2 ! 0  ! This should maybe become an argument in the namelist if we use different hi-es files?

    ! the fixed lengths of various land-surface grids
    integer, parameter :: kSOIL_GRID_Z       = 4
    integer            :: kSNOW_GRID_Z       = 3
    integer            :: kSNOWSOIL_GRID_Z   = 7
    integer, parameter :: kCROP_GRID_Z       = 5
    integer, parameter :: kMONTH_GRID_Z      = 12
    integer, parameter :: kGECROS_GRID_Z     = 60
    integer, parameter :: kSOILCOMP_GRID_Z   = 8
    
    integer, parameter :: kLAKE_Z            = 10
    integer, parameter :: kLAKE_SOISNO_Z     = 9
    integer, parameter :: kLAKE_SOI_Z        = 4
    integer, parameter :: kLAKE_SOISNO_1_Z   = 10


    ! SNOWPACK constants
    integer, parameter :: kSNOWPACK_ATMOS_STAB_MO_HOLTSLAG = 1
    integer, parameter :: kSNOWPACK_ATMOS_STAB_MO_MICHLMAYR = 2

    integer, parameter :: kSNOWPACK_ALBEDO_PARAM_LEHNING_2 = 1
    integer, parameter :: kSNOWPACK_ALBEDO_PARAM_SCHMUCKI_OGS = 2

    integer, parameter :: kSNOWPACK_VARIANT_ANTARCTICA = 1
    integer, parameter :: kSNOWPACK_VARIANT_ALPS      = 2

    ! mm of accumulated precip before "tipping" into the bucket
    ! only performed on output operations
    integer, parameter :: kPRECIP_BUCKET_SIZE=100

! DR commented out September 2023, confusing to have these defined here AND in wrf_constants. Defer to WRF constants, since these are likely more robust
! than the values here, if/when they differ.

! ------------------------------------------------
! Physical Constants
! ------------------------------------------------
    !real, parameter :: LH_vaporization=2260000.0 ! J/kg
    !! could be calculated as 2.5E6 + (-2112.0)*temp_degC ?
    !real, parameter :: Rd  = 287   ! J/(kg K) specific gas constant for dry air
    !real, parameter :: Rw  = 461.6     ! J/(kg K) specific gas constant for moist air
    !real, parameter :: cp  = 1004.5    ! J/kg/K   specific heat capacity of moist STP air?
    !real, parameter :: gravity= 9.81   ! m/s^2    gravity
    !real, parameter :: pi  = 3.1415927 ! pi
    !real, parameter :: stefan_boltzmann = 5.67e-8 ! the Stefan-Boltzmann constant
    !real, parameter :: karman = 0.41   ! the von Karman constant
    !real, parameter :: solar_constant = 1366 ! W/m^2

    ! convenience parameters for various physics packages
    !real, parameter :: rovcp = Rd/cp
    !real, parameter :: rovg  = Rd/gravity

    ! from wrf module_model_constants
    ! parameters for calculating latent heat as a function of temperature for
    ! vaporization
    !real, parameter ::  XLV0 = 3.15E6
    !real, parameter ::  XLV1 = 2370.
    ! sublimation
    !real, parameter ::  XLS0 = 2.905E6
    !real, parameter ::  XLS1 = 259.532

    ! saturated vapor pressure parameters (?)
    !real, parameter ::  SVP1 = 0.6112
    !real, parameter ::  SVP2 = 17.67
    !real, parameter ::  SVP3 = 29.65
    !real, parameter ::  SVPT0= 273.15

    !real, parameter ::  EP1  = Rw/Rd-1.
    !real, parameter ::  EP2  = Rd/Rw

end module
