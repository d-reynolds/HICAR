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

    !Flag-value to indicate a part of a read-write buffer which was never filled
    real, parameter :: kEMPT_BUFF = -123456789.0
    
    ! string lengths
    integer, parameter :: kMAX_FILE_LENGTH =   1024
    integer, parameter :: kMAX_DIM_LENGTH  =   256
    integer, parameter :: kMAX_NAME_LENGTH =   256
    integer, parameter :: kMAX_ATTR_LENGTH =   256
    integer, parameter :: kMAX_STRING_LENGTH = 256  ! maximum length of other strings (e.g. netcdf attributes)

    ! maximum number of nests
    integer, parameter :: kMAX_NESTS = 10
    !>--------------------------------------------
    ! list of integer constants to be used when accessing various arrays that track variable allocation, usage, etc. requests
    !
    ! NOTE: IF YOU ADD TO THIS LIST BE SURE TO ADD AN INTEGER TO THE kVARS STRUCTURE CONSTRUCTOR BELOW IT!
    ! This could be transitioned to an enum... but then one can't just "use, only:kVARS"...
    ! enum, bind(C)
    !   enumerator ::  u, v, w,...
    ! end enum
    ! --------------------------------------------
    type var_constants_type
        SEQUENCE    ! technically SEQUENCE just requires the compiler leave them in order,
                    ! but it can also keep compilers (e.g. ifort) from padding for alignment,
                    ! as long as there is no padding we can test last_var = sizeof(kVARS)

        integer :: u = 1
        integer :: v = 2
        integer :: w = 3
        integer :: w_real = 4
        integer :: pressure = 5
        integer :: pressure_interface = 6
        integer :: potential_temperature = 7
        integer :: temperature = 8
        integer :: water_vapor = 9
        integer :: cloud_water_mass = 10
        integer :: cloud_number = 11
        integer :: ice_mass = 12
        integer :: ice_number = 13
        integer :: rain_mass = 14
        integer :: rain_number = 15
        integer :: snow_mass = 16
        integer :: snow_number = 17
        integer :: graupel_mass = 18
        integer :: graupel_number = 19
        integer :: ice1_a = 20
        integer :: ice1_c = 21
        integer :: ice2_mass = 22
        integer :: ice2_number = 23
        integer :: ice2_a = 24
        integer :: ice2_c = 25
        integer :: ice3_mass = 26
        integer :: ice3_number = 27
        integer :: ice3_a = 28
        integer :: ice3_c = 29
        integer :: precipitation = 30
        integer :: convective_precipitation = 31
        integer :: snowfall = 32
        integer :: graupel = 33
        integer :: snowfall_ground = 34
        integer :: rainfall_ground = 35
        integer :: exner = 36
        integer :: nsquared = 37
        integer :: density = 38
        integer :: cloud_fraction = 39
        integer :: shortwave = 40
        integer :: shortwave_direct = 41
        integer :: shortwave_diffuse = 42
        integer :: shortwave_direct_above = 43
        integer :: shortwave_total = 44
        integer :: longwave = 45
        integer :: albedo = 46
        integer :: vegetation_fraction = 47
        integer :: vegetation_fraction_max = 48
        integer :: vegetation_fraction_out = 49
        integer :: veg_type = 50
        integer :: mass_leaf = 51
        integer :: mass_root = 52
        integer :: mass_stem = 53
        integer :: mass_wood = 54
        integer :: soil_type = 55
        integer :: soil_texture_1 = 56
        integer :: soil_texture_2 = 57
        integer :: soil_texture_3 = 58
        integer :: soil_texture_4 = 59
        integer :: soil_sand_and_clay = 60
        integer :: soil_carbon_stable = 61
        integer :: soil_carbon_fast = 62
        integer :: lai = 63
        integer :: sai = 64
        integer :: crop_category = 65
        integer :: crop_type = 66
        integer :: date_planting = 67
        integer :: date_harvest = 68
        integer :: growing_season_gdd = 69
        integer :: irr_frac_total = 70
        integer :: irr_frac_sprinkler = 71
        integer :: irr_frac_micro = 72
        integer :: irr_frac_flood = 73
        integer :: irr_eventno_sprinkler = 74
        integer :: irr_eventno_micro = 75
        integer :: irr_eventno_flood = 76
        integer :: irr_alloc_sprinkler = 77
        integer :: irr_alloc_micro = 78
        integer :: irr_alloc_flood = 79
        integer :: irr_evap_loss_sprinkler = 80
        integer :: irr_amt_sprinkler = 81
        integer :: irr_amt_micro = 82
        integer :: irr_amt_flood = 83
        integer :: evap_heat_sprinkler = 84
        integer :: mass_ag_grain = 85
        integer :: growing_degree_days = 86
        integer :: plant_growth_stage = 87
        integer :: net_ecosystem_exchange = 88
        integer :: gross_primary_prod = 89
        integer :: net_primary_prod = 90
        integer :: apar = 91
        integer :: photosynthesis_total = 92
        integer :: stomatal_resist_total = 93
        integer :: stomatal_resist_sun = 94
        integer :: stomatal_resist_shade = 95
        integer :: gecros_state = 96
        integer :: canopy_water = 97
        integer :: canopy_water_ice = 98
        integer :: canopy_water_liquid = 99
        integer :: canopy_vapor_pressure = 100
        integer :: canopy_temperature = 101
        integer :: canopy_fwet = 102
        integer :: veg_leaf_temperature = 103
        integer :: ground_surf_temperature = 104
        integer :: frac_between_gap = 105
        integer :: frac_within_gap = 106
        integer :: ground_temperature_bare = 107
        integer :: ground_temperature_canopy = 108
        integer :: sensible_heat = 109
        integer :: latent_heat = 110
        integer :: u_10m = 111
        integer :: v_10m = 112
        integer :: windspd_10m = 113
        integer :: ustar = 114
        integer :: coeff_momentum_drag = 115
        integer :: chs = 116
        integer :: chs2 = 117
        integer :: cqs2 = 118
        integer :: coeff_heat_exchange_3d = 119
        integer :: coeff_momentum_exchange_3d = 120
        integer :: QFX = 121
        integer :: br = 122
        integer :: mol = 123
        integer :: psim = 124
        integer :: psih = 125
        integer :: fm = 126
        integer :: fh = 127
        integer :: surface_rad_temperature = 128
        integer :: temperature_2m = 129
        integer :: humidity_2m = 130
        integer :: temperature_2m_veg = 131
        integer :: temperature_2m_bare = 132
        integer :: mixing_ratio_2m_veg = 133
        integer :: mixing_ratio_2m_bare = 134
        integer :: surface_pressure = 135
        integer :: rad_absorbed_total = 136
        integer :: rad_absorbed_veg = 137
        integer :: rad_absorbed_bare = 138
        integer :: rad_net_longwave = 139
        integer :: longwave_up = 140
        integer :: ground_heat_flux = 141
        integer :: evap_canopy = 142
        integer :: evap_soil_surface = 143
        integer :: transpiration_rate = 144
        integer :: ch_veg = 145
        integer :: ch_veg_2m = 146
        integer :: ch_bare = 147
        integer :: ch_bare_2m = 148
        integer :: ch_under_canopy = 149
        integer :: ch_leaf = 150
        integer :: sensible_heat_veg = 151
        integer :: sensible_heat_bare = 152
        integer :: sensible_heat_canopy = 153
        integer :: evap_heat_veg = 154
        integer :: evap_heat_bare = 155
        integer :: evap_heat_canopy = 156
        integer :: transpiration_heat = 157
        integer :: ground_heat_veg = 158
        integer :: ground_heat_bare = 159
        integer :: net_longwave_veg = 160
        integer :: net_longwave_bare = 161
        integer :: net_longwave_canopy = 162
        integer :: runoff_surface = 163
        integer :: runoff_subsurface = 164
        integer :: soil_totalmoisture = 165
        integer :: soil_deep_temperature = 166
        integer :: water_table_depth = 167
        integer :: water_aquifer = 168
        integer :: storage_gw = 169
        integer :: storage_lake = 170
        integer :: roughness_z0 = 171
        integer :: snow_water_equivalent = 172
        integer :: snow_water_eq_prev = 173
        integer :: snow_albedo_prev = 174
        integer :: snow_temperature = 175
        integer :: snow_layer_depth = 176
        integer :: snow_layer_ice = 177
        integer :: snow_layer_liquid_water = 178
        integer :: snow_age_factor = 179
        integer :: snow_height = 180
        integer :: snow_nlayers = 181
        integer :: soil_water_content = 182
        integer :: soil_water_content_liq = 183
        integer :: eq_soil_moisture = 184
        integer :: smc_watertable_deep = 185
        integer :: recharge = 186
        integer :: recharge_deep = 187
        integer :: soil_temperature = 188
        integer :: skin_temperature = 189
        integer :: sst = 190
        integer :: tend_qv_adv = 191
        integer :: tend_qv_pbl = 192
        integer :: tend_qv = 193
        integer :: tend_th = 194
        integer :: tend_th_pbl = 195
        integer :: tend_qc = 196
        integer :: tend_qc_pbl = 197
        integer :: tend_qi = 198
        integer :: tend_qi_pbl = 199
        integer :: tend_qs = 200
        integer :: tend_qr = 201
        integer :: tend_u = 202
        integer :: tend_v = 203
        integer :: u_mass = 204
        integer :: v_mass = 205
        integer :: re_cloud = 206
        integer :: re_ice = 207
        integer :: re_snow = 208
        integer :: ice1_rho = 209
        integer :: ice1_phi = 210
        integer :: ice1_vmi = 211
        integer :: ice2_rho = 212
        integer :: ice2_phi = 213
        integer :: ice2_vmi = 214
        integer :: ice3_rho = 215
        integer :: ice3_phi = 216
        integer :: ice3_vmi = 217
        integer :: wind_alpha = 218
        integer :: froude = 219
        integer :: blk_ri = 220
        integer :: out_longwave_rad = 221
        integer :: longwave_cloud_forcing = 222
        integer :: shortwave_cloud_forcing = 223
        integer :: land_emissivity = 224
        integer :: temperature_interface = 225
        integer :: tend_swrad = 226
        integer :: runoff_tstep = 227
        integer :: Tsnow = 228
        integer :: Sice = 229
        integer :: Sliq = 230
        integer :: Ds = 231
        integer :: fsnow = 232
        integer :: Nsnow = 233
        integer :: dSWE_salt = 234
        integer :: dSWE_susp = 235
        integer :: dSWE_subl = 236
        integer :: dSWE_slide = 237
        integer :: meltflux_out_tstep = 238
        integer :: Sliq_out = 239
        integer :: kpbl = 240
        integer :: hpbl = 241
        integer :: lake_depth = 242
        integer :: t_lake3d = 243
        integer :: snl2d = 244
        integer :: t_grnd2d = 245
        integer :: lake_icefrac3d = 246
        integer :: z_lake3d = 247
        integer :: dz_lake3d = 248
        integer :: t_soisno3d = 249
        integer :: h2osoi_ice3d = 250
        integer :: h2osoi_liq3d = 251
        integer :: h2osoi_vol3d = 252
        integer :: z3d = 253
        integer :: dz3d = 254
        integer :: watsat3d = 255
        integer :: csol3d = 256
        integer :: tkmg3d = 257
        integer :: lakemask = 258
        integer :: xice = 259
        integer :: zi3d = 260
        integer :: tksatu3d = 261
        integer :: tkdry3d = 262
        integer :: savedtke12d = 263
        integer :: lakedepth2d = 264
        integer :: ivt = 265
        integer :: iwv = 266
        integer :: iwl = 267
        integer :: iwi = 268
        ! GRID VARIABLES
        integer :: z = 269
        integer :: z_interface = 270
        integer :: dzdx = 271
        integer :: dzdy = 272
        integer :: dz = 273
        integer :: dz_interface = 274
        integer :: advection_dz = 275
        integer :: dzdy_v = 276
        integer :: dzdx_u = 277
        integer :: jacobian = 278
        integer :: jacobian_u = 279
        integer :: jacobian_v = 280
        integer :: jacobian_w = 281
        integer :: land_mask = 282
        integer :: terrain = 283
        integer :: latitude = 284
        integer :: longitude = 285
        integer :: global_terrain = 286
        integer :: global_dz_interface = 287
        integer :: global_z_interface = 288
        integer :: u_latitude = 289
        integer :: u_longitude = 290
        integer :: v_latitude = 291
        integer :: v_longitude = 292
        integer :: Sx = 293
        integer :: TPI = 294
        integer :: neighbor_terrain = 295
        integer :: froude_terrain = 296
        integer :: relax_filter_2d = 297
        integer :: relax_filter_3d = 298
        integer :: costheta = 299
        integer :: sintheta = 300
        integer :: cosine_zenith_angle = 301
        integer :: slope = 302
        integer :: slope_angle = 303
        integer :: aspect_angle = 304
        integer :: svf = 305
        integer :: shd = 306
        integer :: hlm = 307
        integer :: h1 = 308
        integer :: h2 = 309
        integer :: h1_u = 310
        integer :: h1_v = 311
        integer :: h2_u = 312
        integer :: h2_v = 313
        integer :: last_var = 314
    end type var_constants_type


    type(var_constants_type) :: kVARS != var_constants_type(   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  &
    !                                                          11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  &
    !                                                          21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  &
    !                                                          31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  &
    !                                                          41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  &
    !                                                          51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  &
    !                                                          61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  &
    !                                                          71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  &
    !                                                          81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  &
    !                                                          91,  92,  93,  94,  95,  96,  97,  98,  99, 100,  &
    !                                                         101, 102, 103, 104, 105, 106, 107, 108, 109, 110,  &
    !                                                         111, 112, 113, 114, 115, 116, 117, 118, 119, 120,  &
    !                                                         121, 122, 123, 124, 125, 126, 127, 128, 129, 130,  &
    !                                                         131, 132, 133, 134, 135, 136, 137, 138, 139, 140,  &
    !                                                         141, 142, 143, 144, 145, 146, 147, 148, 149, 150,  &
    !                                                         151, 152, 153, 154, 155, 156, 157, 158, 159, 160,  &
    !                                                         161, 162, 163, 164, 165, 166, 167, 168, 169, 170,  &
    !                                                         171, 172, 173, 174, 175, 176, 177, 178, 179, 180,  &
    !                                                         181, 182, 183, 184, 185, 186, 187, 188, 189, 190,  &
    !                                                         191, 192, 193, 194, 195, 196, 197, 198, 199, 200,  &
    !                                                         201, 202, 203, 204, 205, 206, 207, 208, 209, 210,  &
    !                                                         211, 212, 213, 214, 215, 216, 217, 218, 219, 220,  &
    !                                                         221, 222, 223, 224, 225, 226, 227, 228, 229, 230,  &
    !                                                         231, 232, 233, 234, 235, 236, 237, 238, 239, 240,  &
    !                                                         241, 242, 243, 244, 245, 246, 247, 248, 249, 250,  &
    !                                                         251, 252, 253, 254, 255, 256, 257, 258, 259, 260,  &
    !                                                         261, 262, 263, 264, 265, 266, 267, 268, 269, 270,  &
    !                                                         271, 272, 273, 274, 275, 276, 277, 278, 279, 280,  &
    !                                                         281, 282, 283, 284, 285, 286, 287, 288, 289, 290,  &
    !                                                         291, 289, 293, 294, 295, 296, 297, 298, 299, 300,  &
    !                                                         301, 302, 303, 304, 305, 306, 307, 308, 309, 310)

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
    integer, parameter :: kLSM_NOAH      = 2
    integer, parameter :: kLSM_NOAHMP    = 3
    
    integer, parameter :: kSM_FSM        = 1 !! MJ added

    integer, parameter :: kRA_BASIC      = 1
    integer, parameter :: kRA_SIMPLE     = 2
    integer, parameter :: kRA_RRTMG      = 3

    integer, parameter :: kADV_STD       = 1
    integer, parameter :: kADV_MPDATA    = 2

    integer, parameter :: kFLUXCOR_MONO   = 1

    integer, parameter :: kWIND_LINEAR   = 1
    integer, parameter :: kOBRIEN_WINDS  = 2
    integer, parameter :: kITERATIVE_WINDS = 3
    integer, parameter :: kLINEAR_OBRIEN_WINDS = 4
    integer, parameter :: kLINEAR_ITERATIVE_WINDS = 5

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
