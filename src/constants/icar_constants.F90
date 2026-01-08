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
        integer :: pressure_base = 7
        integer :: geopotential_base = 8
        integer :: potential_temperature = 9
        integer :: temperature = 10
        integer :: water_vapor = 11
        integer :: cloud_water_mass = 12
        integer :: cloud_number = 13
        integer :: ice_mass = 14
        integer :: ice_number = 15
        integer :: rain_mass = 16
        integer :: rain_number = 17
        integer :: snow_mass = 18
        integer :: snow_number = 19
        integer :: graupel_mass = 20
        integer :: graupel_number = 21
        integer :: ice1_a = 22
        integer :: ice1_c = 23
        integer :: ice2_mass = 24
        integer :: ice2_number = 25
        integer :: ice2_a = 26
        integer :: ice2_c = 27
        integer :: ice3_mass = 28
        integer :: ice3_number = 29
        integer :: ice3_a = 30
        integer :: ice3_c = 31
        integer :: precipitation = 32
        integer :: convective_precipitation = 33
        integer :: snowfall = 34
        integer :: graupel = 35
        integer :: snowfall_ground = 36
        integer :: rainfall_ground = 37
        integer :: lsm_last_snow = 38
        integer :: lsm_last_precip = 39
        integer :: exner = 40
        integer :: nsquared = 41
        integer :: density = 42
        integer :: cloud_fraction = 43
        integer :: shortwave = 44
        integer :: shortwave_direct = 45
        integer :: shortwave_diffuse = 46
        integer :: shortwave_direct_above = 47
        integer :: longwave = 48
        integer :: albedo = 49
        integer :: vegetation_fraction = 50
        integer :: vegetation_fraction_max = 51
        integer :: vegetation_fraction_out = 52
        integer :: veg_type = 53
        integer :: mass_leaf = 54
        integer :: mass_root = 55
        integer :: mass_stem = 56
        integer :: mass_wood = 57
        integer :: soil_type = 58
        integer :: soil_texture_1 = 59
        integer :: soil_texture_2 = 60
        integer :: soil_texture_3 = 61
        integer :: soil_texture_4 = 62
        integer :: soil_sand_and_clay = 63
        integer :: soil_carbon_stable = 64
        integer :: soil_carbon_fast = 65
        integer :: lai = 66
        integer :: sai = 67
        integer :: crop_category = 68
        integer :: crop_type = 69
        integer :: date_planting = 70
        integer :: date_harvest = 71
        integer :: growing_season_gdd = 72
        integer :: irr_frac_total = 73
        integer :: irr_frac_sprinkler = 74
        integer :: irr_frac_micro = 75
        integer :: irr_frac_flood = 76
        integer :: irr_eventno_sprinkler = 77
        integer :: irr_eventno_micro = 78
        integer :: irr_eventno_flood = 79
        integer :: irr_alloc_sprinkler = 80
        integer :: irr_alloc_micro = 81
        integer :: irr_alloc_flood = 82
        integer :: irr_evap_loss_sprinkler = 83
        integer :: irr_amt_sprinkler = 84
        integer :: irr_amt_micro = 85
        integer :: irr_amt_flood = 86
        integer :: evap_heat_sprinkler = 87
        integer :: mass_ag_grain = 88
        integer :: growing_degree_days = 89
        integer :: plant_growth_stage = 90
        integer :: net_ecosystem_exchange = 91
        integer :: gross_primary_prod = 92
        integer :: net_primary_prod = 93
        integer :: apar = 94
        integer :: photosynthesis_total = 95
        integer :: stomatal_resist_total = 96
        integer :: stomatal_resist_sun = 97
        integer :: stomatal_resist_shade = 98
        integer :: gecros_state = 99
        integer :: canopy_water = 100
        integer :: canopy_water_ice = 101
        integer :: canopy_water_liquid = 102
        integer :: canopy_vapor_pressure = 103
        integer :: canopy_temperature = 104
        integer :: canopy_fwet = 105
        integer :: veg_leaf_temperature = 106
        integer :: ground_surf_temperature = 107
        integer :: frac_between_gap = 108
        integer :: frac_within_gap = 109
        integer :: ground_temperature_bare = 110
        integer :: ground_temperature_canopy = 111
        integer :: sensible_heat = 112
        integer :: latent_heat = 113
        integer :: u_10m = 114
        integer :: v_10m = 115
        integer :: windspd_10m = 116
        integer :: ustar = 117
        integer :: coeff_momentum_drag = 118
        integer :: chs = 119
        integer :: chs2 = 120
        integer :: cqs2 = 121
        integer :: coeff_heat_exchange_3d = 122
        integer :: coeff_momentum_exchange_3d = 123
        integer :: QFX = 124
        integer :: br = 125
        integer :: mol = 126
        integer :: psim = 127
        integer :: psih = 128
        integer :: fm = 129
        integer :: fh = 130
        integer :: surface_rad_temperature = 131
        integer :: temperature_2m = 132
        integer :: humidity_2m = 133
        integer :: temperature_2m_veg = 134
        integer :: temperature_2m_bare = 135
        integer :: mixing_ratio_2m_veg = 136
        integer :: mixing_ratio_2m_bare = 137
        integer :: surface_pressure = 138
        integer :: sea_surface_pressure = 139
        integer :: rad_absorbed_total = 140
        integer :: rad_absorbed_veg = 141
        integer :: rad_absorbed_bare = 142
        integer :: rad_net_longwave = 143
        integer :: longwave_up = 144
        integer :: ground_heat_flux = 145
        integer :: evap_canopy = 146
        integer :: evap_soil_surface = 147
        integer :: transpiration_rate = 148
        integer :: ch_veg = 149
        integer :: ch_veg_2m = 150
        integer :: ch_bare = 151
        integer :: ch_bare_2m = 152
        integer :: ch_under_canopy = 153
        integer :: ch_leaf = 154
        integer :: sensible_heat_veg = 155
        integer :: sensible_heat_bare = 156
        integer :: sensible_heat_canopy = 157
        integer :: evap_heat_veg = 158
        integer :: evap_heat_bare = 159
        integer :: evap_heat_canopy = 160
        integer :: transpiration_heat = 161
        integer :: ground_heat_veg = 162
        integer :: ground_heat_bare = 163
        integer :: net_longwave_veg = 164
        integer :: net_longwave_bare = 165
        integer :: net_longwave_canopy = 166
        integer :: runoff_surface = 167
        integer :: runoff_subsurface = 168
        integer :: soil_totalmoisture = 169
        integer :: soil_deep_temperature = 170
        integer :: water_table_depth = 171
        integer :: water_aquifer = 172
        integer :: storage_gw = 173
        integer :: storage_lake = 174
        integer :: roughness_z0 = 175
        integer :: snow_water_equivalent = 176
        integer :: snow_water_eq_prev = 177
        integer :: snow_albedo_prev = 178
        integer :: snow_temperature = 179
        integer :: snow_layer_depth = 180
        integer :: snow_layer_ice = 181
        integer :: snow_layer_liquid_water = 182
        integer :: snow_age_factor = 183
        integer :: snow_height = 184
        integer :: snow_nlayers = 185
        integer :: soil_water_content = 186
        integer :: soil_water_content_liq = 187
        integer :: eq_soil_moisture = 188
        integer :: smc_watertable_deep = 189
        integer :: recharge = 190
        integer :: recharge_deep = 191
        integer :: soil_temperature = 192
        integer :: skin_temperature = 193
        integer :: sst = 194
        integer :: tend_qv_adv = 195
        integer :: tend_qv_pbl = 196
        integer :: tend_qv = 197
        integer :: tend_th = 198
        integer :: tend_th_pbl = 199
        integer :: tend_qc = 200
        integer :: tend_qc_pbl = 201
        integer :: tend_qi = 202
        integer :: tend_qi_pbl = 203
        integer :: tend_qs = 204
        integer :: tend_qr = 205
        integer :: tend_u = 206
        integer :: tend_v = 207
        integer :: tend_th_lwrad = 208
        integer :: tend_th_swrad = 209
        integer :: u_mass = 210
        integer :: v_mass = 211
        integer :: re_cloud = 212
        integer :: re_ice = 213
        integer :: re_snow = 214
        integer :: ice1_rho = 215
        integer :: ice1_phi = 216
        integer :: ice1_vmi = 217
        integer :: ice2_rho = 218
        integer :: ice2_phi = 219
        integer :: ice2_vmi = 220
        integer :: ice3_rho = 221
        integer :: ice3_phi = 222
        integer :: ice3_vmi = 223
        integer :: wind_alpha = 224
        integer :: froude = 225
        integer :: blk_ri = 226
        integer :: out_longwave_rad = 227
        integer :: longwave_cloud_forcing = 228
        integer :: shortwave_cloud_forcing = 229
        integer :: land_emissivity = 230
        integer :: temperature_interface = 231
        integer :: tend_swrad = 232
        integer :: runoff_tstep = 233
        integer :: Tsnow = 234
        integer :: Sice = 235
        integer :: Sliq = 236
        integer :: Ds = 237
        integer :: fsnow = 238
        integer :: Nsnow = 239
        integer :: dSWE_salt = 240
        integer :: dSWE_susp = 241
        integer :: dSWE_subl = 242
        integer :: dSWE_slide = 243
        integer :: meltflux_out_tstep = 244
        integer :: Sliq_out = 245
        integer :: kpbl = 246
        integer :: hpbl = 247
        integer :: lake_depth = 248
        integer :: t_lake3d = 249
        integer :: snl2d = 250
        integer :: t_grnd2d = 251
        integer :: lake_icefrac3d = 252
        integer :: z_lake3d = 253
        integer :: dz_lake3d = 254
        integer :: t_soisno3d = 255
        integer :: h2osoi_ice3d = 256
        integer :: h2osoi_liq3d = 257
        integer :: h2osoi_vol3d = 258
        integer :: z3d = 259
        integer :: dz3d = 260
        integer :: watsat3d = 261
        integer :: csol3d = 262
        integer :: tkmg3d = 263
        integer :: lakemask = 264
        integer :: xice = 265
        integer :: zi3d = 266
        integer :: tksatu3d = 267
        integer :: tkdry3d = 268
        integer :: savedtke12d = 269
        integer :: lakedepth2d = 270
        integer :: ivt = 271
        integer :: iwv = 272
        integer :: iwl = 273
        integer :: iwi = 274
        ! GRID VARIABLES
        integer :: z = 275
        integer :: z_interface = 276
        integer :: dzdx = 277
        integer :: dzdy = 278
        integer :: dz = 279
        integer :: dz_interface = 280
        integer :: advection_dz = 281
        integer :: dzdy_v = 282
        integer :: dzdx_u = 283
        integer :: jacobian = 284
        integer :: jacobian_u = 285
        integer :: jacobian_v = 286
        integer :: jacobian_w = 287
        integer :: land_mask = 288
        integer :: terrain = 289
        integer :: latitude = 290
        integer :: longitude = 291
        integer :: global_terrain = 292
        integer :: global_dz_interface = 293
        integer :: global_z_interface = 294
        integer :: u_latitude = 295
        integer :: u_longitude = 296
        integer :: v_latitude = 297
        integer :: v_longitude = 298
        integer :: Sx = 299
        integer :: TPI = 300
        integer :: neighbor_terrain = 301
        integer :: froude_terrain = 302
        integer :: relax_filter_2d = 303
        integer :: relax_filter_3d = 304
        integer :: costheta = 305
        integer :: sintheta = 306
        integer :: cosine_zenith_angle = 307
        integer :: slope = 308
        integer :: slope_angle = 309
        integer :: aspect_angle = 310
        integer :: svf = 311
        integer :: shd = 312
        integer :: hlm = 313
        integer :: h1 = 314
        integer :: h2 = 315
        integer :: h1_u = 316
        integer :: h1_v = 317
        integer :: h2_u = 318
        integer :: h2_v = 319
        integer :: last_var = 320
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
    integer, parameter :: kLSM_NOAH      = 2
    integer, parameter :: kLSM_NOAHMP    = 3
    
    integer, parameter :: kSM_FSM        = 1 !! MJ added

    integer, parameter :: kRA_BASIC      = 1
    integer, parameter :: kRA_SIMPLE     = 2
    integer, parameter :: kRA_RRTMG      = 3
    integer, parameter :: kRA_RRTMGP    = 4
    
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
