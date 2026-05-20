!>------------------------------------------------------------
!! FLake lake model -- GPU-portable single-column core.
!!
!! Vendored from the official FLake LGPL distribution
!! (Mironov, DWD; file timestamps 2005-12-17).
!! Source: http://www.flake.igb-berlin.de/ (down at vendoring time; archive
!! supplied locally). The 6 constant/type modules below are verbatim except
!! for renaming `data_parameters` -> `flake_data_parameters` to avoid any
!! naming collision with host-model modules.
!!
!! GPU portability refactor (HICAR-specific changes only in MODULE flake_core):
!!  - The 46 module-level scalars that held per-column FLake state in the
!!    upstream `MODULE flake` (flake.f90:117-176) are bundled into the
!!    derived type `flake_state_t` and passed as `INTENT(INOUT) :: s` to
!!    flake_radflux and flake_driver. Inside those routines, every state-var
!!    reference `X_flk` is now `s%X_flk` (purely textual substitution; no
!!    semantic change). This makes the routines reentrant -- safe to call
!!    from a `!$acc parallel loop gang vector` over horizontal cells.
!!  - All five procedures carry `!$acc routine seq` for device dispatch.
!!  - The pure helpers (flake_buoypar, flake_snowdensity, flake_snowheatconduct)
!!    were already reentrant; only the directive is added.
!!
!! Upstream re-merge: re-run `helpers/assemble_flake_core.sh` after dropping
!! a fresh FLake tarball into ~/Downloads/src_flake_sfcflx.
!!
!! Reference: Mironov, D. V., 2008: Parameterization of lakes in numerical
!!   weather prediction. Description of a lake model. COSMO Tech. Rep. 11, DWD.
!!------------------------------------------------------------

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Data module for global parameters
!-------------------------------------------------------------------------------

MODULE flake_data_parameters

!-------------------------------------------------------------------------------
!
! Description:
!  Global parameters for the program are defined.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated intgribf, intgribc, irealgrib, iwlength and put it to data_io.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminated parameters for grid point and diagnostic calculations.
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:
!=======================================================================

IMPLICIT NONE

!=======================================================================
! Global (i.e. public) Declarations:
! Parameters for the Program:

  INTEGER, PARAMETER       ::                                         &
       ireals    = SELECTED_REAL_KIND (12,200),                       &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers

!=======================================================================

END MODULE flake_data_parameters

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_derivedtypes  

!------------------------------------------------------------------------------
!
! Description:
!
!  Derived type(s) is(are) defined.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Maximum value of the wave-length bands 
!  in the exponential decay law for the radiation flux.
!  A storage for a ten-band approximation is allocated,
!  although a smaller number of bands is actually used.
INTEGER (KIND = iintegers), PARAMETER :: & 
  nband_optic_max = 10_iintegers

!  Define TYPE "opticpar_medium"
TYPE opticpar_medium
  INTEGER (KIND = iintegers)                        ::   & 
    nband_optic                                            ! Number of wave-length bands
  REAL (KIND = ireals), DIMENSION (nband_optic_max) ::   & 
    frac_optic                                         , & ! Fractions of total radiation flux 
    extincoef_optic                                        ! Extinction coefficients 
END TYPE opticpar_medium

!==============================================================================

END MODULE flake_derivedtypes  

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_parameters

!------------------------------------------------------------------------------
!
! Description:
!
!  Values of empirical constants of the lake model FLake 
!  and of several thermodynamic parameters are set.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants 
!  in the equations for the mixed-layer depth 
!  and for the shape factor with respect to the temperature profile in the thermocline
REAL (KIND = ireals), PARAMETER ::         &
  c_cbl_1       = 0.17_ireals            , & ! Constant in the CBL entrainment equation
  c_cbl_2       = 1._ireals              , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n    = 0.5_ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_s    = 10._ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_i    = 20._ireals             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_relax_h     = 0.030_ireals           , & ! Constant in the relaxation equation for the SBL depth
  c_relax_C     = 0.0030_ireals              ! Constant in the relaxation equation for the shape factor
                                             ! with respect to the temperature profile in the thermocline

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.
REAL (KIND = ireals), PARAMETER ::         &
  C_T_min       = 0.5_ireals             , & ! Minimum value of the shape factor C_T (thermocline)
  C_T_max       = 0.8_ireals             , & ! Maximum value of the shape factor C_T (thermocline)
  Phi_T_pr0_1   = 40._ireals/3._ireals   , & ! Constant in the expression for the T shape-function derivative 
  Phi_T_pr0_2   = 20._ireals/3._ireals   , & ! Constant in the expression for the T shape-function derivative 
  C_TT_1        = 11._ireals/18._ireals  , & ! Constant in the expression for C_TT (thermocline)
  C_TT_2        = 7._ireals/45._ireals   , & ! Constant in the expression for C_TT (thermocline)
  C_B1          = 2._ireals/3._ireals    , & ! Shape factor (upper layer of bottom sediments)
  C_B2          = 3._ireals/5._ireals    , & ! Shape factor (lower layer of bottom sediments)
  Phi_B1_pr0    = 2._ireals              , & ! B1 shape-function derivative 
  C_S_lin       = 0.5_ireals             , & ! Shape factor (linear temperature profile in the snow layer)
  Phi_S_pr0_lin = 1._ireals              , & ! S shape-function derivative (linear profile) 
  C_I_lin       = 0.5_ireals             , & ! Shape factor (linear temperature profile in the ice layer)
  Phi_I_pr0_lin = 1._ireals              , & ! I shape-function derivative (linear profile) 
  Phi_I_pr1_lin = 1._ireals              , & ! I shape-function derivative (linear profile) 
  Phi_I_ast_MR  = 2._ireals              , & ! Constant in the MR2004 expression for I shape factor
  C_I_MR        = 1._ireals/12._ireals   , & ! Constant in the MR2004 expression for I shape factor
  H_Ice_max     = 3._ireals                  ! Maximum ice tickness in 
                                             ! the Mironov and Ritter (2004, MR2004) ice model [m] 

!  Security constants
REAL (KIND = ireals), PARAMETER ::         &
  h_Snow_min_flk = 1.0E-5_ireals         , & ! Minimum snow thickness [m]
  h_Ice_min_flk  = 1.0E-9_ireals         , & ! Minimum ice thickness [m]
  h_ML_min_flk   = 1.0E-2_ireals         , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk   = 1.0E+3_ireals         , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk   = 1.0E-3_ireals         , & ! Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6_ireals             ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = ireals), PARAMETER ::         &
  c_small_flk    = 1.0E-10_ireals            ! A small number

!  Thermodynamic parameters
REAL (KIND = ireals), PARAMETER ::        &
  tpl_grav          = 9.81_ireals       , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r           = 277.13_ireals     , & ! Temperature of maximum density of fresh water [K]
  tpl_T_f           = 273.15_ireals     , & ! Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05_ireals , & ! Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03_ireals    , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_I         = 9.1E+02_ireals    , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min     = 1.0E+02_ireals    , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max     = 4.0E+02_ireals    , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S   = 2.0E+02_ireals    , & ! Empirical parameter [kg m^{-4}]  
                                            ! in the expression for the snow density 
  tpl_L_f           = 3.3E+05_ireals    , & ! Latent heat of fusion [J kg^{-1}]
  tpl_c_w           = 4.2E+03_ireals    , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_c_I           = 2.1E+03_ireals    , & ! Specific heat of ice [J kg^{-1} K^{-1}]
  tpl_c_S           = 2.1E+03_ireals    , & ! Specific heat of snow [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01_ireals   , & ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I       = 2.29_ireals       , & ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min   = 0.2_ireals        , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max   = 1.5_ireals        , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S = 1.3_ireals            ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] 
                                            ! in the expression for the snow heat conductivity 

!==============================================================================

END MODULE flake_parameters

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_paramoptic_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of the optical characteristics
!  of the lake water, lake ice and snow. These reference values may be used 
!  if no information about the optical characteristics of the lake in question 
!  is available. An exponential decay law for the solar radiation flux is assumed.
!  In the simplest one-band approximation,
!  the extinction coefficient for water is set to a large value,
!  leading to the absorption of 95% of the incoming radiation 
!  within the uppermost 1 m of the lake water. 
!  The extinction coefficients for ice and snow are taken from 
!  Launiainen and Cheng (1998). The estimates for the ice correspond 
!  to the uppermost 0.1 m of the ice layer and to the clear sky conditions 
!  (see Table 2 in op. cit.).
!  Very large values of the extinction coefficients for ice and snow ("opaque")
!  can be used to prevent penetration of the solar radiation 
!  through the snow-ice cover.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters, ONLY :      &
  ireals                       , & ! KIND-type parameter for real variables 
  iintegers                        ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes, ONLY :   &
  nband_optic_max              , & ! Maximum value of the wave-length bands
  opticpar_medium                  ! Derived TYPE 

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
INTEGER (KIND = iintegers), PRIVATE :: & ! Help variable(s)
  i                                      ! DO loop index

!  Optical characteristics for water, ice and snow.
!  The simplest one-band approximation is used as a reference.
TYPE (opticpar_medium), PARAMETER ::                           & 
  opticpar_water_ref = opticpar_medium(1,                      & ! Water (reference)
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/3._ireals, (1.E+10_ireals,i=2,nband_optic_max)/))      , &
  opticpar_water_trans = opticpar_medium(2,                             & ! Transparent Water (two-band)
    (/0.10_ireals, 0.90_ireals, (0._ireals,i=3,nband_optic_max)/),      &
    (/2.0_ireals, 0.20_ireals, (1.E+10_ireals,i=3,nband_optic_max)/)) , &
!_nu  opticpar_water_trans = opticpar_medium(1,                    & ! Transparent Water (one-band)
!_nu    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
!_nu    (/0.30_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  opticpar_whiteice_ref = opticpar_medium(1,                   & ! White ice
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &   
    (/17.1_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  opticpar_blueice_ref = opticpar_medium(1,                    & ! Blue ice
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/8.4_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))     , &
  opticpar_drysnow_ref = opticpar_medium(1,                    & ! Dry snow 
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/25.0_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  opticpar_meltingsnow_ref = opticpar_medium(1,                & ! Melting snow 
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/15.0_ireals, (1.E+10_ireals,i=2,nband_optic_max)/))    , &
  opticpar_ice_opaque = opticpar_medium(1,                     & ! Opaque ice
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/1.0E+07_ireals, (1.E+10_ireals,i=2,nband_optic_max)/)) , &
  opticpar_snow_opaque = opticpar_medium(1,                    & ! Opaque snow
    (/1._ireals, (0._ireals,i=2,nband_optic_max)/),            &
    (/1.0E+07_ireals, (1.E+10_ireals,i=2,nband_optic_max)/)) 

!==============================================================================

END MODULE flake_paramoptic_ref

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_albedo_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of albedo 
!  for the lake water, lake ice and snow. 
!  As in "flake_paramoptic_ref", two ice categories, viz. white ice and blue ice,
!  and two snow categories, viz. dry snow and melting snow, are used.  
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters, ONLY :      &
  ireals                       , & ! KIND-type parameter for real variables 
  iintegers                        ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Albedo for water, ice and snow.
REAL (KIND = ireals), PARAMETER ::        &
  albedo_water_ref       = 0.07_ireals  , & ! Water
  albedo_whiteice_ref    = 0.60_ireals  , & ! White ice
  albedo_blueice_ref     = 0.10_ireals  , & ! Blue ice
  albedo_drysnow_ref     = 0.60_ireals  , & ! Dry snow 
  albedo_meltingsnow_ref = 0.10_ireals      ! Melting snow 

!  Empirical parameters.
REAL (KIND = ireals), PARAMETER :: &
  c_albice_MR = 95.6_ireals          ! Constant in the interpolation formula for 
                                     ! the ice albedo (Mironov and Ritter 2004)

!==============================================================================

END MODULE flake_albedo_ref
! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flake_configure

!------------------------------------------------------------------------------
!
! Description:
!
!  Switches and reference values of parameters 
!  that configure the lake model FLake are set.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release 
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

LOGICAL, PARAMETER :: &
  lflk_botsed_use   = .TRUE.         ! .TRUE. indicates that the bottom-sediment scheme is used
                                     ! to compute the depth penetrated by the thermal wave, 
                                     ! the temperature at this depth and the bottom heat flux.
                                     ! Otherwise, the heat flux at the water-bottom sediment interface
                                     ! is set to zero, the depth penetrated by the thermal wave 
                                     ! is set to a reference value defined below,
                                     ! and the temperature at this depth is set to 
                                     ! the temperature of maximum density of the fresh water.

REAL (KIND = ireals), PARAMETER :: &
  rflk_depth_bs_ref = 10.0_ireals    ! Reference value of the depth of the thermally active
                                     ! layer of bottom sediments [m].
                                     ! This value is used to (formally) define
                                     ! the depth penetrated by the thermal wave
                                     ! in case the bottom-sediment scheme is not used.

!==============================================================================

END MODULE flake_configure


!==============================================================================
!  HICAR GPU-portable wrapping module.
!==============================================================================

MODULE flake_core

USE flake_data_parameters, ONLY : ireals, iintegers
USE flake_derivedtypes      ! opticpar_medium, nband_optic_max
USE flake_parameters        ! physical constants, Phi_*, C_*, tpl_*, *_min_flk
USE flake_configure         ! lflk_botsed_use, rflk_depth_bs_ref

IMPLICIT NONE

!  Per-column FLake state. Bundles every scalar that the upstream
!  `MODULE flake` declared at module scope (flake.f90:117-176), so the
!  state becomes per-thread when passed as an argument.
TYPE :: flake_state_t
    ! Temperatures at the previous ("_p") and updated ("_n") time step [K]
    REAL(KIND=ireals) :: T_mnw_p_flk,  T_mnw_n_flk      ! Mean water-column temperature
    REAL(KIND=ireals) :: T_snow_p_flk, T_snow_n_flk     ! Air-snow interface temperature
    REAL(KIND=ireals) :: T_ice_p_flk,  T_ice_n_flk      ! Snow-ice or air-ice interface temperature
    REAL(KIND=ireals) :: T_wML_p_flk,  T_wML_n_flk      ! Mixed-layer temperature
    REAL(KIND=ireals) :: T_bot_p_flk,  T_bot_n_flk      ! Water-bottom sediment interface temperature
    REAL(KIND=ireals) :: T_B1_p_flk,   T_B1_n_flk       ! Bottom-of-upper-sediment-layer temperature

    ! Layer thicknesses [m]
    REAL(KIND=ireals) :: h_snow_p_flk, h_snow_n_flk
    REAL(KIND=ireals) :: h_ice_p_flk,  h_ice_n_flk
    REAL(KIND=ireals) :: h_ML_p_flk,   h_ML_n_flk
    REAL(KIND=ireals) :: H_B1_p_flk,   H_B1_n_flk

    ! Shape factors (thermocline / ice / snow)
    REAL(KIND=ireals) :: C_T_p_flk, C_T_n_flk
    REAL(KIND=ireals) :: C_TT_flk, C_Q_flk
    REAL(KIND=ireals) :: C_I_flk, C_S_flk

    ! Shape-function derivatives
    REAL(KIND=ireals) :: Phi_T_pr0_flk, Phi_I_pr0_flk, Phi_I_pr1_flk, Phi_S_pr0_flk

    ! Heat and radiation fluxes [W m^-2]
    REAL(KIND=ireals) :: Q_snow_flk, Q_ice_flk, Q_w_flk, Q_bot_flk
    REAL(KIND=ireals) :: I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk
    REAL(KIND=ireals) :: I_intm_0_h_flk, I_intm_h_D_flk
    REAL(KIND=ireals) :: Q_star_flk

    ! Velocity scales [m s^-1]
    REAL(KIND=ireals) :: u_star_w_flk, w_star_sfc_flk

    ! Rate of snow accumulation [kg m^-2 s^-1]
    REAL(KIND=ireals) :: dMsnowdt_flk
END TYPE flake_state_t

CONTAINS

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

SUBROUTINE flake_radflux ( s, depth_w, albedo_water, albedo_ice, albedo_snow, & 
                           opticpar_water, opticpar_ice, opticpar_snow )       

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean radiation flux over the mixed layer,
!  and the mean radiation flux over the thermocline.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes          ! Definitions of derived TYPEs

USE flake_parameters , ONLY : & 
  h_Snow_min_flk            , & ! Minimum snow thickness [m]
  h_Ice_min_flk             , & ! Minimum ice thickness [m]
  h_ML_min_flk                  ! Minimum mixed-layer depth [m]

!==============================================================================

IMPLICIT NONE

!$acc routine seq

TYPE(flake_state_t), INTENT(INOUT) :: s

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = ireals), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  albedo_water                      , & ! Albedo of the water surface 
  albedo_ice                        , & ! Albedo of the ice surface
  albedo_snow                           ! Albedo of the snow surface

TYPE (opticpar_medium), INTENT(IN) :: & 
  opticpar_water                    , & ! Optical characteristics of water
  opticpar_ice                      , & ! Optical characteristics of ice
  opticpar_snow                         ! Optical characteristics of snow 


!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: & ! Help variable(s)
  i                             ! DO loop index

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

  IF(s%h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
    IF(s%h_snow_p_flk.GE.h_Snow_min_flk) THEN        ! There is snow above the ice
      s%I_snow_flk = s%I_atm_flk*(1._ireals-albedo_snow) 
      s%I_bot_flk = 0._ireals
      DO i=1, opticpar_snow%nband_optic
        s%I_bot_flk = s%I_bot_flk +                    & 
        opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)*s%h_snow_p_flk) 
      END DO 
      s%I_ice_flk  = s%I_snow_flk*s%I_bot_flk
    ELSE                                           ! No snow above the ice 
      s%I_snow_flk = s%I_atm_flk  
      s%I_ice_flk  = s%I_atm_flk*(1._ireals-albedo_ice)
    END IF 
    s%I_bot_flk = 0._ireals
    DO i=1, opticpar_ice%nband_optic
      s%I_bot_flk = s%I_bot_flk +                      & 
      opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)*s%h_ice_p_flk) 
    END DO 
    s%I_w_flk      = s%I_ice_flk*s%I_bot_flk
  ELSE                                             ! No ice-snow cover
    s%I_snow_flk   = s%I_atm_flk  
    s%I_ice_flk    = s%I_atm_flk
    s%I_w_flk      = s%I_atm_flk*(1._ireals-albedo_water)
  END IF 

  IF(s%h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Radiation flux at the bottom of the mixed layer
    s%I_bot_flk = 0._ireals
    DO i=1, opticpar_water%nband_optic
      s%I_bot_flk = s%I_bot_flk +            & 
      opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*s%h_ML_p_flk) 
    END DO 
    s%I_h_flk = s%I_w_flk*s%I_bot_flk
  ELSE                                          ! Mixed-layer depth is less then a minimum value
    s%I_h_flk = s%I_w_flk
  END IF

  s%I_bot_flk = 0._ireals                         ! Radiation flux at the lake bottom
  DO i=1, opticpar_water%nband_optic
    s%I_bot_flk = s%I_bot_flk +              & 
    opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*depth_w) 
  END DO 
  s%I_bot_flk = s%I_w_flk*s%I_bot_flk

  IF(s%h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Integral-mean radiation flux over the mixed layer
    s%I_intm_0_h_flk = 0._ireals
    DO i=1, opticpar_water%nband_optic
      s%I_intm_0_h_flk = s%I_intm_0_h_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      (1._ireals - EXP(-opticpar_water%extincoef_optic(i)*s%h_ML_p_flk))
    END DO 
    s%I_intm_0_h_flk = s%I_w_flk*s%I_intm_0_h_flk/s%h_ML_p_flk
  ELSE
    s%I_intm_0_h_flk = s%I_h_flk
  END IF

  IF(s%h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN   ! Integral-mean radiation flux over the thermocline
    s%I_intm_h_D_flk = 0._ireals 
    DO i=1, opticpar_water%nband_optic
      s%I_intm_h_D_flk = s%I_intm_h_D_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      ( EXP(-opticpar_water%extincoef_optic(i)*s%h_ML_p_flk)             &
      - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
    END DO 
    s%I_intm_h_D_flk = s%I_w_flk*s%I_intm_h_D_flk/(depth_w-s%h_ML_p_flk)
  ELSE
    s%I_intm_h_D_flk = s%I_h_flk
  END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_radflux

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

SUBROUTINE flake_driver ( s, depth_w, depth_bs, T_bs, par_Coriolis,       &
                          extincoef_water_typ,                         &
                          del_time, T_sfc_p, T_sfc_n )         

!------------------------------------------------------------------------------
!
! Description:
!
!  The main driving routine of the lake model FLake 
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used 
!  for debugging purposes only.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters            ! Thermodynamic parameters and dimensionless constants of FLake

USE flake_configure             ! Switches and parameters that configure FLake

!==============================================================================

IMPLICIT NONE

!$acc routine seq

TYPE(flake_state_t), INTENT(INOUT) :: s

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = ireals), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                          , & ! Depth of the thermally active layer of bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  extincoef_water_typ               , & ! "Typical" extinction coefficient of the lake water [m^{-1}],
                                        ! used to compute the equilibrium CBL depth
  del_time                          , & ! The model time step [s]
  T_sfc_p                               ! Surface temperature at the previous time step [K]  
                                        ! (equal to either T_ice, T_snow or to T_wML)

!  Output (procedure arguments)

REAL (KIND = ireals), INTENT(OUT) ::  &
  T_sfc_n                               ! Updated surface temperature [K] 
                                        ! (equal to the updated value of either T_ice, T_snow or T_wML)


!  Local variables of type LOGICAL
LOGICAL ::          &
  l_ice_create    , & ! Switch, .TRUE. = ice does not exist but should be created
  l_snow_exists   , & ! Switch, .TRUE. = there is snow above the ice
  l_ice_meltabove     ! Switch, .TRUE. = snow/ice melting from above takes place

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = ireals) ::    &
  d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
  d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
  d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
  d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
  d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
  d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
  d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
  d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
  d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = ireals) ::    &
  N_T_mean               , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
  ZM_h_scale             , & ! The ZM96 equilibrium SBL depth scale [m] 
  conv_equil_h_scale         ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL
REAL (KIND = ireals) :: &
  h_ice_threshold     , & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
  flk_str_1           , & ! Help storage variable
  flk_str_2           , & ! Help storage variable
  R_H_icesnow         , & ! Dimensionless ratio, used to store intermediate results
  R_rho_c_icesnow     , & ! Dimensionless ratio, used to store intermediate results
  R_TI_icesnow        , & ! Dimensionless ratio, used to store intermediate results
  R_Tstar_icesnow         ! Dimensionless ratio, used to store intermediate results


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!_dm 
! Security. Set time-rate-of-change of prognostic variables to zero.
! Set prognostic variables to their values at the previous time step.
! (This is to avoid spurious changes of prognostic variables 
! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
! at the neighbouring lake points as noticed by Burkhardt Rockel.)
!_dm 

d_T_mnw_dt   = 0._ireals 
d_T_ice_dt   = 0._ireals 
d_T_bot_dt   = 0._ireals 
d_T_B1_dt    = 0._ireals 
d_h_snow_dt  = 0._ireals 
d_h_ice_dt   = 0._ireals 
d_h_ML_dt    = 0._ireals 
d_H_B1_dt    = 0._ireals 
d_C_T_dt     = 0._ireals 
s%T_snow_n_flk = s%T_snow_p_flk   
s%T_ice_n_flk  = s%T_ice_p_flk    
s%T_wML_n_flk  = s%T_wML_p_flk   
s%T_mnw_n_flk  = s%T_mnw_p_flk     
s%T_bot_n_flk  = s%T_bot_p_flk  
s%T_B1_n_flk   = s%T_B1_p_flk      
s%h_snow_n_flk = s%h_snow_p_flk 
s%h_ice_n_flk  = s%h_ice_p_flk   
s%h_ML_n_flk   = s%h_ML_p_flk    
s%H_B1_n_flk   = s%H_B1_p_flk   
s%C_T_n_flk    = s%C_T_p_flk    

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_dm
! At this point, the heat and radiation fluxes, namely,
! s%Q_snow_flk, s%Q_ice_flk, s%Q_w_flk, 
! s%I_atm_flk, s%I_snow_flk, s%I_ice_flk, s%I_w_flk, s%I_h_flk, s%I_bot_flk,     
! the mean radiation flux over the mixed layer, s%I_intm_0_h_flk, 
! and the mean radiation flux over the thermocline, s%I_intm_h_D_flk, 
! should be known.
! They are computed within "flake_interface" (or within the driving model)
! and are available to "flake_driver"
! through the above variables declared in the MODULE "flake".
! In case a lake is ice-covered, s%Q_w_flk is re-computed below.
!_dm

! Heat flux through the ice-water interface
IF(s%h_ice_p_flk.GE.h_Ice_min_flk) THEN    ! Ice exists 
  IF(s%h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Mixed-layer depth is zero, compute flux 
    s%Q_w_flk = -tpl_kappa_w*(s%T_bot_p_flk-s%T_wML_p_flk)/depth_w  ! Flux with linear T(z) 
    s%Phi_T_pr0_flk = Phi_T_pr0_1*s%C_T_p_flk-Phi_T_pr0_2         ! d\Phi(0)/d\zeta (thermocline)
    s%Q_w_flk = s%Q_w_flk*MAX(s%Phi_T_pr0_flk, 1._ireals)           ! Account for an increased d\Phi(0)/d\zeta 
  ELSE                    
    s%Q_w_flk = 0._ireals                  ! Mixed-layer depth is greater than zero, set flux to zero
  END IF   
END IF   

! A generalized heat flux scale 
s%Q_star_flk = s%Q_w_flk + s%I_w_flk + s%I_h_flk - 2._ireals*s%I_intm_0_h_flk

! Heat flux through the water-bottom sediment interface
IF(lflk_botsed_use) THEN
  s%Q_bot_flk = -tpl_kappa_w*(s%T_B1_p_flk-s%T_bot_p_flk)/MAX(s%H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
ELSE  
  s%Q_bot_flk = 0._ireals   ! The bottom-sediment scheme is not used
END IF


!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

!_dm
! Notice that a quasi-equilibrium ice-snow model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_dm

!_dev
! The dependence of snow density and of snow heat conductivity 
! on the snow thickness is accounted for parametrically.
! That is, the time derivatives of \rho_S and \kappa_S are neglected.
! The exception is the equation for the snow thickness 
! in case of snow accumulation and no melting, 
! where d\rho_S/dt is incorporated.
! Furthermore, some (presumably small) correction terms incorporating 
! the snow density and the snow heat conductivity are dropped out.
! Those terms may be included as better formulations 
! for \rho_S and \kappa_S are available.
!_dev

! Default values
l_ice_create    = .FALSE.  
l_ice_meltabove = .FALSE.  

Ice_exist: IF(s%h_ice_p_flk.LT.h_Ice_min_flk) THEN   ! Ice does not exist 

  l_ice_create = s%T_wML_p_flk.LE.(tpl_T_f+c_small_flk).AND.s%Q_w_flk.LT.0._ireals
  IF(l_ice_create) THEN                            ! Ice does not exist but should be created
    d_h_ice_dt = -s%Q_w_flk/tpl_rho_I/tpl_L_f                                  
    s%h_ice_n_flk = s%h_ice_p_flk + d_h_ice_dt*del_time                          ! Advance h_ice 
    s%T_ice_n_flk = tpl_T_f + s%h_ice_n_flk*s%Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
    d_h_snow_dt = s%dMsnowdt_flk/tpl_rho_S_min 
    s%h_snow_n_flk = s%h_snow_p_flk + d_h_snow_dt*del_time                       ! Advance h_snow
    s%Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1._ireals, s%h_ice_n_flk/H_Ice_max)       ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = s%Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(s%h_snow_n_flk) &
                * s%h_snow_n_flk/MAX(s%h_ice_n_flk, h_Ice_min_flk)
    s%T_snow_n_flk = s%T_ice_n_flk + R_H_icesnow*(s%T_ice_n_flk-tpl_T_f)           ! Snow temperature
  END IF

ELSE Ice_exist                                     ! Ice exists

  l_snow_exists = s%h_snow_p_flk.GE.h_Snow_min_flk   ! Check if there is snow above the ice

  Melting: IF(s%T_snow_p_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
                                                           ! T_snow = T_ice if snow is absent 
    IF(l_snow_exists) THEN   ! There is snow above the ice
      flk_str_1 = s%Q_snow_flk + s%I_snow_flk - s%I_ice_flk        ! Atmospheric forcing
      IF(flk_str_1.GE.0._ireals) THEN  ! Melting of snow and ice from above
        l_ice_meltabove = .TRUE.
        d_h_snow_dt = (-flk_str_1/tpl_L_f+s%dMsnowdt_flk)/flake_snowdensity(s%h_snow_p_flk)
        d_h_ice_dt  = -(s%I_ice_flk - s%I_w_flk - s%Q_w_flk)/tpl_L_f/tpl_rho_I 
      END IF 
    ELSE                     ! No snow above the ice
      flk_str_1 = s%Q_ice_flk + s%I_ice_flk - s%I_w_flk - s%Q_w_flk  ! Atmospheric forcing + heating from the water
      IF(flk_str_1.GE.0._ireals) THEN  ! Melting of ice from above, snow accumulation may occur
        l_ice_meltabove = .TRUE.
        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
        d_h_snow_dt = s%dMsnowdt_flk/tpl_rho_S_min
      END IF 
    END IF 
    IF(l_ice_meltabove) THEN  ! Melting from above takes place
      s%h_ice_n_flk  = s%h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
      s%h_snow_n_flk = s%h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
      s%T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
      s%T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
    END IF

  END IF Melting

  No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above

    d_h_snow_dt = flake_snowdensity(s%h_snow_p_flk)  
    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
     flk_str_1 = s%h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
     flk_str_1 = flk_str_1/(1._ireals-flk_str_1)
    ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
     flk_str_1 = 0._ireals
    END IF
    d_h_snow_dt = s%dMsnowdt_flk/d_h_snow_dt/(1._ireals+flk_str_1)       ! Snow accumulation
    s%h_snow_n_flk = s%h_snow_p_flk + d_h_snow_dt*del_time                         ! Advance h_snow
    
    s%Phi_I_pr0_flk = s%h_ice_p_flk/H_Ice_max                              ! h_ice relative to its maximum value
    s%C_I_flk = C_I_lin - C_I_MR*(1._ireals+Phi_I_ast_MR)*s%Phi_I_pr0_flk  ! Shape factor (ice)
    s%Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*s%Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
    s%Phi_I_pr0_flk = Phi_I_pr0_lin - s%Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)

    h_ice_threshold = MAX(1._ireals, 2._ireals*s%C_I_flk*tpl_c_I*(tpl_T_f-s%T_ice_p_flk)/tpl_L_f)
    h_ice_threshold = s%Phi_I_pr0_flk/s%C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   ! Threshold value of h_ice
    h_ice_threshold = MIN(0.9_ireals*H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
                                                                       ! h_ice(threshold) < 0.9*H_Ice_max

    IF(s%h_ice_p_flk.LT.h_ice_threshold) THEN  ! Use a quasi-equilibrium ice model

      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
        flk_str_1 = s%Q_snow_flk + s%I_snow_flk - s%I_w_flk
      ELSE                     ! Use fluxes at the air-ice interface
        flk_str_1 = s%Q_ice_flk + s%I_ice_flk - s%I_w_flk
      END IF
      d_h_ice_dt = -(flk_str_1-s%Q_w_flk)/tpl_L_f/tpl_rho_I
      s%h_ice_n_flk = s%h_ice_p_flk + d_h_ice_dt *del_time                         ! Advance h_ice
      s%T_ice_n_flk = tpl_T_f + s%h_ice_n_flk*flk_str_1/tpl_kappa_I/s%Phi_I_pr0_flk  ! Ice temperature

    ELSE                                     ! Use a complete ice model

      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-s%T_ice_p_flk)/s%h_ice_p_flk*s%Phi_I_pr0_flk
      d_h_ice_dt = (s%Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
      s%h_ice_n_flk = s%h_ice_p_flk  + d_h_ice_dt*del_time                         ! Advance h_ice

      R_TI_icesnow = tpl_c_I*(tpl_T_f-s%T_ice_p_flk)/tpl_L_f         ! Dimensionless parameter
      R_Tstar_icesnow = 1._ireals - s%C_I_flk                        ! Dimensionless parameter
      IF(l_snow_exists) THEN  ! There is snow above the ice
        R_H_icesnow = s%Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(s%h_snow_p_flk) &
                    * s%h_snow_p_flk/s%h_ice_p_flk
        R_rho_c_icesnow = flake_snowdensity(s%h_snow_p_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 
!_dev 
!_dm 
! These terms should be included as an improved understanding of the snow scheme is gained, 
! of the effect of snow density in particular. 
!_dm 
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
!_nu                        + (1._ireals+C_S_lin*s%h_snow_p_flk/s%h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow
!_dev

        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter

!_dev
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
!_nu                        + (1._ireals-R_rho_c_icesnow)*tpl_c_I*s%T_ice_p_flk/tpl_L_f
!_dev
        flk_str_2 = s%Q_snow_flk+s%I_snow_flk-s%I_w_flk                  ! Atmospheric fluxes
        flk_str_1  = s%C_I_flk*s%h_ice_p_flk + (1._ireals+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*s%h_snow_p_flk
        d_T_ice_dt = -(1._ireals-2._ireals*C_S_lin)*R_H_icesnow*(tpl_T_f-s%T_ice_p_flk)             & 
                   * tpl_c_S*s%dMsnowdt_flk                          ! Effect of snow accumulation
      ELSE                    ! No snow above the ice
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
        flk_str_2 = s%Q_ice_flk+s%I_ice_flk-s%I_w_flk                    ! Atmospheric fluxes
        flk_str_1  = s%C_I_flk*s%h_ice_p_flk
        d_T_ice_dt = 0._ireals
      END IF 
      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-s%T_ice_p_flk)/s%h_ice_p_flk*s%Phi_I_pr0_flk       &
                 * (1._ireals-R_Tstar_icesnow)                     ! Add flux due to heat conduction
      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*s%Q_w_flk            ! Add flux from water to ice
      d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
      d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
      s%T_ice_n_flk = s%T_ice_p_flk + d_T_ice_dt*del_time                          ! Advance T_ice
    END IF

    s%Phi_I_pr1_flk = MIN(1._ireals, s%h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
    s%Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*s%Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = s%Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(s%h_snow_n_flk) &
                 *s%h_snow_n_flk/MAX(s%h_ice_n_flk, h_Ice_min_flk)
    s%T_snow_n_flk = s%T_ice_n_flk + R_H_icesnow*(s%T_ice_n_flk-tpl_T_f)             ! Snow temperature

  END IF No_Melting

END IF Ice_exist   

! Security, limit h_ice by its maximum value
s%h_ice_n_flk = MIN(s%h_ice_n_flk, H_Ice_max)      

! Security, limit the ice and snow temperatures by the freezing point 
s%T_snow_n_flk = MIN(s%T_snow_n_flk, tpl_T_f)  
s%T_ice_n_flk =  MIN(s%T_ice_n_flk,  tpl_T_f)    

!_tmp
! Security, avoid too low values (these constraints are used for debugging purposes)
  s%T_snow_n_flk = MAX(s%T_snow_n_flk, 73.15_ireals)  
  s%T_ice_n_flk =  MAX(s%T_ice_n_flk,  73.15_ireals)    
!_tmp

! Remove too thin ice and/or snow
IF(s%h_ice_n_flk.LT.h_Ice_min_flk)  THEN        ! Check ice
  s%h_ice_n_flk = 0._ireals       ! Ice is too thin, remove it, and
  s%T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
  s%h_snow_n_flk = 0._ireals      ! Remove snow when there is no ice, and
  s%T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
  l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but proved to be too thin
ELSE IF(s%h_snow_n_flk.LT.h_Snow_min_flk) THEN  ! Ice exists, check snow
  s%h_snow_n_flk = 0._ireals      ! Snow is too thin, remove it, 
  s%T_snow_n_flk = s%T_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
END IF


!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF(l_ice_create) s%Q_w_flk = 0._ireals     ! Ice has just been created, set Q_w to zero
d_T_mnw_dt = (s%Q_w_flk - s%Q_bot_flk + s%I_w_flk - s%I_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w
s%T_mnw_n_flk = s%T_mnw_p_flk + d_T_mnw_dt*del_time   ! Advance T_mnw
s%T_mnw_n_flk = MAX(s%T_mnw_n_flk, tpl_T_f)           ! Limit T_mnw by the freezing point 


!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_Water: IF(s%h_ice_n_flk.GE.h_Ice_min_flk) THEN    ! Ice exists

  s%T_mnw_n_flk = MIN(s%T_mnw_n_flk, tpl_T_r) ! Limit the mean temperature under the ice by T_r 
  s%T_wML_n_flk = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 

  IF(l_ice_create) THEN                  ! Ice has just been created 
    IF(s%h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN    ! h_ML=D when ice is created 
      s%h_ML_n_flk = 0._ireals                 ! Set h_ML to zero 
      s%C_T_n_flk = C_T_min                    ! Set C_T to its minimum value 
    ELSE                                          ! h_ML<D when ice is created 
      s%h_ML_n_flk = s%h_ML_p_flk                ! h_ML remains unchanged 
      s%C_T_n_flk = s%C_T_p_flk                  ! C_T (thermocline) remains unchanged 
    END IF 
    s%T_bot_n_flk = s%T_wML_n_flk - (s%T_wML_n_flk-s%T_mnw_n_flk)/s%C_T_n_flk/(1._ireals-s%h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE IF(s%T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
    s%h_ML_n_flk = s%h_ML_p_flk                  ! h_ML remains unchanged 
    s%C_T_n_flk = s%C_T_p_flk                    ! C_T (thermocline) remains unchanged 
    s%T_bot_n_flk = s%T_wML_n_flk - (s%T_wML_n_flk-s%T_mnw_n_flk)/s%C_T_n_flk/(1._ireals-s%h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
    s%T_bot_n_flk = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
    IF(s%h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
      s%C_T_n_flk = s%C_T_p_flk                     ! C_T (thermocline) remains unchanged 
      s%h_ML_n_flk = depth_w*(1._ireals-(s%T_wML_n_flk-s%T_mnw_n_flk)/(s%T_wML_n_flk-s%T_bot_n_flk)/s%C_T_n_flk)
      s%h_ML_n_flk = MAX(s%h_ML_n_flk, 0._ireals)   ! Update the mixed-layer depth  
    ELSE                                 ! h_ML = 0 
      s%h_ML_n_flk = s%h_ML_p_flk                   ! h_ML remains unchanged 
      s%C_T_n_flk = (s%T_wML_n_flk-s%T_mnw_n_flk)/(s%T_wML_n_flk-s%T_bot_n_flk) 
      s%C_T_n_flk = MIN(C_T_max, MAX(s%C_T_n_flk, C_T_min)) ! Update the shape factor (thermocline)  
    END IF 
  END IF 

  s%T_bot_n_flk = MIN(s%T_bot_n_flk, tpl_T_r)    ! Security, limit the bottom temperature by T_r 

ELSE HTC_Water                                      ! Open water

! Generalised buoyancy flux scale and convective velocity scale
  flk_str_1 = flake_buoypar(s%T_wML_p_flk)*s%Q_star_flk/tpl_rho_w_r/tpl_c_w                    
  IF(flk_str_1.LT.0._ireals) THEN       
    s%w_star_sfc_flk = (-flk_str_1*s%h_ML_p_flk)**(1._ireals/3._ireals)  ! Convection     
  ELSE 
    s%w_star_sfc_flk = 0._ireals                                       ! Neutral or stable stratification
  END IF 

!_dm
! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
! is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_dm
  conv_equil_h_scale = -s%Q_w_flk/MAX(s%I_w_flk, c_small_flk)
  IF(conv_equil_h_scale.GT.0._ireals .AND. conv_equil_h_scale.LT.1._ireals  &
    .AND. s%T_wML_p_flk.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
    conv_equil_h_scale = SQRT(6._ireals*conv_equil_h_scale)                 &
                       + 2._ireals*conv_equil_h_scale/(1._ireals-conv_equil_h_scale)
    conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
  ELSE
    conv_equil_h_scale = 0._ireals       ! Set the equilibrium CBL depth to zero
  END IF

! Mean buoyancy frequency in the thermocline
  N_T_mean = flake_buoypar(0.5_ireals*(s%T_wML_p_flk+s%T_bot_p_flk))*(s%T_wML_p_flk-s%T_bot_p_flk)
  IF(s%h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
    N_T_mean = SQRT(N_T_mean/(depth_w-s%h_ML_p_flk))  ! Compute N                   
  ELSE 
    N_T_mean = 0._ireals                            ! h_ML=D, set N to zero
  END IF 

! The rate of change of C_T
  d_C_T_dt = MAX(s%w_star_sfc_flk, s%u_star_w_flk, u_star_min_flk)**2_iintegers
  d_C_T_dt = N_T_mean*(depth_w-s%h_ML_p_flk)**2_iintegers       &
           / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
  d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     ! Rate-of-change of C_T 

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  s%C_TT_flk = C_TT_1*s%C_T_p_flk-C_TT_2         ! C_TT, using C_T at the previous time step
  s%C_Q_flk = 2._ireals*s%C_TT_flk/s%C_T_p_flk     ! C_Q using C_T at the previous time step

  Mixing_regime: IF(flk_str_1.LT.0._ireals) THEN  ! Convective mixing 

    s%C_T_n_flk = s%C_T_p_flk + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
    s%C_T_n_flk = MIN(C_T_max, MAX(s%C_T_n_flk, C_T_min))                ! Limit C_T 
    d_C_T_dt = (s%C_T_n_flk-s%C_T_p_flk)/del_time                        ! Re-compute dC_T/dt

    IF(s%h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
      IF(s%h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Use a reduced entrainment equation (spin-up)
        d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(s%w_star_sfc_flk, c_small_flk)

!_dbg
! WRITE*, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', d_h_ML_dt*del_time
! WRITE*, '         w_*       = ', s%w_star_sfc_flk
! WRITE*, '         \beta*Q_* = ', flk_str_1
!_dbg

      ELSE                                   ! Use a complete entrainment equation 
        R_H_icesnow     = depth_w/s%h_ML_p_flk
        R_rho_c_icesnow = R_H_icesnow-1._ireals
        R_TI_icesnow    = s%C_T_p_flk/s%C_TT_flk
        R_Tstar_icesnow = (R_TI_icesnow/2._ireals-1._ireals)*R_rho_c_icesnow + 1._ireals
        d_h_ML_dt = -s%Q_star_flk*(R_Tstar_icesnow*(1._ireals+c_cbl_1)-1._ireals) - s%Q_bot_flk
        d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
        flk_str_2 = (depth_w-s%h_ML_p_flk)*(s%T_wML_p_flk-s%T_bot_p_flk)*C_TT_2/s%C_TT_flk*d_C_T_dt 
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
        flk_str_2 = s%I_bot_flk + (R_TI_icesnow-1._ireals)*s%I_h_flk - R_TI_icesnow*s%I_intm_h_D_flk
        flk_str_2 = flk_str_2 + (R_TI_icesnow-2._ireals)*R_rho_c_icesnow*(s%I_h_flk-s%I_intm_0_h_flk)
        flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
        flk_str_2 = -c_cbl_2*R_Tstar_icesnow*s%Q_star_flk/tpl_rho_w_r/tpl_c_w/MAX(s%w_star_sfc_flk, c_small_flk)
        flk_str_2 = flk_str_2 + s%C_T_p_flk*(s%T_wML_p_flk-s%T_bot_p_flk)
        d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
      END IF 
!_dm
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place.
! Physically, this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
! In reality, convective mixed-layer deepening may be completely decoupled
! from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_dm

!_dbg
! IF(d_h_ML_dt.LT.0._ireals) THEN 
!   WRITE*, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
!   WRITE*, '                d_h_ML_dt*del_time = ', MAX(d_h_ML_dt, c_small_flk)*del_time
!   WRITE*, '         u_*       = ', s%u_star_w_flk   
!   WRITE*, '         w_*       = ', s%w_star_sfc_flk
!   WRITE*, '         h_CBL_eqi = ', conv_equil_h_scale
!   WRITE*, '         ZM scale  = ', ZM_h_scale
!   WRITE*, '        s%h_ML_p_flk = ', s%h_ML_p_flk
! END IF
!   WRITE*, 'FLake: Convection, = ', d_h_ML_dt
!   WRITE*, '         Q_*       = ', s%Q_star_flk
!   WRITE*, '         \beta*Q_* = ', flk_str_1
!_dbg

      d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    
      s%h_ML_n_flk = s%h_ML_p_flk + d_h_ML_dt*del_time                       ! Update h_ML 
      s%h_ML_n_flk = MAX(h_ML_min_flk, MIN(s%h_ML_n_flk, depth_w))           ! Security, limit h_ML
    ELSE                                              ! Mixing down to the lake bottom
      s%h_ML_n_flk = depth_w
    END IF

  ELSE Mixing_regime                              ! Wind mixing

    d_h_ML_dt = MAX(s%u_star_w_flk, u_star_min_flk)                        ! The surface friction velocity
    ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2_iintegers
    ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
    ZM_h_scale = d_h_ML_dt**3_iintegers/ZM_h_scale 
    ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))        ! The ZM96 SBL depth scale 
    ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 

!_dm 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_dm 

    d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    s%h_ML_n_flk = ZM_h_scale - (ZM_h_scale-s%h_ML_p_flk)*EXP(-d_h_ML_dt)    ! Update h_ML 
    s%h_ML_n_flk = MAX(h_ML_min_flk, MIN(s%h_ML_n_flk, depth_w))             ! Limit h_ML 
    d_h_ML_dt = (s%h_ML_n_flk-s%h_ML_p_flk)/del_time                         ! Re-compute dh_ML/dt

    IF(s%h_ML_n_flk.LE.s%h_ML_p_flk)           &
      d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
    s%C_T_n_flk = s%C_T_p_flk + d_C_T_dt*del_time                            ! Update C_T
    s%C_T_n_flk = MIN(C_T_max, MAX(s%C_T_n_flk, C_T_min))                    ! Limit C_T 
    d_C_T_dt = (s%C_T_n_flk-s%C_T_p_flk)/del_time                            ! Re-compute dC_T/dt

!_dbg
! WRITE*, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
! WRITE*, '              h_CBL_eqi = ', conv_equil_h_scale
! WRITE*, '              ZM scale  = ', ZM_h_scale
! WRITE*, '              w_*       = ', s%w_star_sfc_flk
! WRITE*, '              u_*       = ', s%u_star_w_flk
! WRITE*, '             s%h_ML_p_flk = ', s%h_ML_p_flk
!_dbg

  END IF Mixing_regime

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(s%h_ML_n_flk.LE.depth_w-h_ML_min_flk) THEN       ! Mixing did not reach the bottom 

    IF(s%h_ML_n_flk.GT.s%h_ML_p_flk) THEN   ! Mixed-layer deepening 
      R_H_icesnow     = s%h_ML_p_flk/depth_w
      R_rho_c_icesnow = 1._ireals-R_H_icesnow 
      R_TI_icesnow    = 0.5_ireals*s%C_T_p_flk*R_rho_c_icesnow+s%C_TT_flk*(2._ireals*R_H_icesnow-1._ireals)
      R_Tstar_icesnow = (0.5_ireals+s%C_TT_flk-s%C_Q_flk)/R_TI_icesnow
      R_TI_icesnow    = (1._ireals-s%C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
     
      d_T_bot_dt = (s%Q_w_flk-s%Q_bot_flk+s%I_w_flk-s%I_bot_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt - s%C_T_p_flk*(s%T_wML_p_flk-s%T_bot_p_flk)*d_h_ML_dt
      d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term

      flk_str_2 = s%I_intm_h_D_flk - (1._ireals-s%C_Q_flk)*s%I_h_flk - s%C_Q_flk*s%I_bot_flk
      flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-s%h_ML_p_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term

      flk_str_2 = (1._ireals-C_TT_2*R_TI_icesnow)/s%C_T_p_flk
      flk_str_2 = flk_str_2*(s%T_wML_p_flk-s%T_bot_p_flk)*d_C_T_dt
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term
      
    ELSE                                ! Mixed-layer retreat or stationary state
      d_T_bot_dt = 0._ireals                                            ! dT_bot/dt=0
    END IF

    s%T_bot_n_flk = s%T_bot_p_flk + d_T_bot_dt*del_time                      ! Update T_bot  
    s%T_bot_n_flk = MAX(s%T_bot_n_flk, tpl_T_f)           ! Security, limit T_bot by the freezing point
    flk_str_2 = (s%T_bot_n_flk-tpl_T_r)*flake_buoypar(s%T_mnw_n_flk)
    IF(flk_str_2.LT.0._ireals) s%T_bot_n_flk = tpl_T_r  ! Security, avoid T_r crossover 
    s%T_wML_n_flk = s%C_T_n_flk*(1._ireals-s%h_ML_n_flk/depth_w)
    s%T_wML_n_flk = (s%T_mnw_n_flk-s%T_bot_n_flk*s%T_wML_n_flk)/(1._ireals-s%T_wML_n_flk)
    s%T_wML_n_flk = MAX(s%T_wML_n_flk, tpl_T_f)           ! Security, limit T_wML by the freezing point

  ELSE                                              ! Mixing down to the lake bottom 

    s%h_ML_n_flk = depth_w
    s%T_wML_n_flk = s%T_mnw_n_flk
    s%T_bot_n_flk = s%T_mnw_n_flk
    s%C_T_n_flk = C_T_min

  END IF

END IF HTC_Water


!------------------------------------------------------------------------------
!  Compute the depth of the upper layer of bottom sediments
!  and the temperature at that depth.
!------------------------------------------------------------------------------

Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  
  IF(s%H_B1_p_flk.GE.depth_bs-H_B1_min_flk) THEN   ! No T(z) maximum (no thermal wave) 
    s%H_B1_p_flk = 0._ireals                       ! Set H_B1_p to zero
    s%T_B1_p_flk = s%T_bot_p_flk                     ! Set T_B1_p to the bottom temperature
  END IF 

  flk_str_1 = 2._ireals*Phi_B1_pr0/(1._ireals-C_B1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
  h_ice_threshold = SQRT(flk_str_1)                              ! Threshold value of H_B1
  h_ice_threshold = MIN(0.9_ireals*depth_bs, h_ice_threshold)    ! Limit H_B1
  flk_str_2 = C_B2/(1._ireals-C_B2)*(T_bs-s%T_B1_p_flk)/(depth_bs-s%H_B1_p_flk)

  IF(s%H_B1_p_flk.LT.h_ice_threshold) THEN  ! Use a truncated equation for H_B1(t)
    s%H_B1_n_flk = SQRT(s%H_B1_p_flk**2_iintegers+flk_str_1)  ! Advance H_B1
    d_H_B1_dt = (s%H_B1_n_flk-s%H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
  ELSE                                    ! Use a full equation for H_B1(t)
    flk_str_1 = (s%Q_bot_flk+s%I_bot_flk)/s%H_B1_p_flk/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - (1._ireals-C_B1)*(s%T_bot_n_flk-s%T_bot_p_flk)/del_time
    d_H_B1_dt = (1._ireals-C_B1)*(s%T_bot_p_flk-s%T_B1_p_flk)/s%H_B1_p_flk + C_B1*flk_str_2
    d_H_B1_dt = flk_str_1/d_H_B1_dt
    s%H_B1_n_flk = s%H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
  END IF 
  d_T_B1_dt = flk_str_2*d_H_B1_dt
  s%T_B1_n_flk = s%T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1

!_dbg
! WRITE*, 'BS module: '
! WRITE*, '  Q_bot   = ', s%Q_bot_flk
! WRITE*, '  d_H_B1_dt = ', d_H_B1_dt
! WRITE*, '  d_T_B1_dt = ', d_T_B1_dt
! WRITE*, '  H_B1    = ', s%H_B1_n_flk
! WRITE*, '    T_bot = ', s%T_bot_n_flk
! WRITE*, '  T_B1    = ', s%T_B1_n_flk
! WRITE*, '    T_bs  = ',  T_bs
!_dbg

!_nu  
! Use a very simplistic procedure, where only the upper layer profile is used, 
! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs.
! Then, the time derivatives are zero, and the sign of the bottom heat flux depends on 
! whether T_bot is smaller or greater than T_bs.
! This is, of course, an oversimplified scheme.
!_nu  d_H_B1_dt = 0._ireals
!_nu  d_T_B1_dt = 0._ireals
!_nu  s%H_B1_n_flk = s%H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
!_nu  s%T_B1_n_flk = s%T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
!_nu  

  l_snow_exists = s%H_B1_n_flk.GE.depth_bs-H_B1_min_flk                    & ! H_B1 reached depth_bs, or
             .OR. s%H_B1_n_flk.LT.H_B1_min_flk                             & ! H_B1 decreased to zero, or
             .OR.(s%T_bot_n_flk-s%T_B1_n_flk)*(T_bs-s%T_B1_n_flk).LE.0._ireals   ! there is no T(z) maximum
  IF(l_snow_exists) THEN      
    s%H_B1_n_flk = depth_bs                     ! Set H_B1 to the depth of the thermally active layer
    s%T_B1_n_flk = T_bs                         ! Set T_B1 to the climatological temperature 
  END IF

ELSE Use_sediment                        ! The bottom-sediment scheme is not used

  s%H_B1_n_flk = rflk_depth_bs_ref              ! H_B1 is set to a reference value 
  s%T_B1_n_flk = tpl_T_r                        ! T_B1 is set to the temperature of maximum density

END IF Use_sediment


!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

! In case of unstable stratification, force mixing down to the bottom
flk_str_2 = (s%T_wML_n_flk-s%T_bot_n_flk)*flake_buoypar(s%T_mnw_n_flk)
IF(flk_str_2.LT.0._ireals) THEN 

!_dbg
! WRITE*, 'FLake: inverse (unstable) stratification !!! '
! WRITE*, '       Mixing down to the bottom is forced.'
! WRITE*, '  T_wML_p, T_wML_n ', s%T_wML_p_flk-tpl_T_f, s%T_wML_n_flk-tpl_T_f
! WRITE*, '  T_mnw_p, T_mnw_n ', s%T_mnw_p_flk-tpl_T_f, s%T_mnw_n_flk-tpl_T_f
! WRITE*, '  T_bot_p, T_bot_n ', s%T_bot_p_flk-tpl_T_f, s%T_bot_n_flk-tpl_T_f
! WRITE*, '  h_ML_p,  h_ML_n  ', s%h_ML_p_flk,          s%h_ML_n_flk
! WRITE*, '  C_T_p,   C_T_n   ', s%C_T_p_flk,           s%C_T_n_flk
!_dbg

  s%h_ML_n_flk = depth_w
  s%T_wML_n_flk = s%T_mnw_n_flk
  s%T_bot_n_flk = s%T_mnw_n_flk
  s%C_T_n_flk = C_T_min

END IF


!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF(s%h_snow_n_flk.GE.h_Snow_min_flk) THEN   
  T_sfc_n = s%T_snow_n_flk                   ! Snow exists, use the snow temperature
ELSE IF(s%h_ice_n_flk.GE.h_Ice_min_flk) THEN
  T_sfc_n = s%T_ice_n_flk                    ! Ice exists but there is no snow, use the ice temperature
ELSE 
  T_sfc_n = s%T_wML_n_flk                    ! No ice-snow cover, use the mixed-layer temperature
END IF

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_driver

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_buoypar (T_water)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the buoyancy parameter,
!  using a quadratic equation of state for the fresh-water.
!  
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_grav                  , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r                   , & ! Temperature of maximum density of fresh water [K]
  tpl_a_T                       ! Constant in the fresh-water equation of state [K^{-2}]

!==============================================================================

IMPLICIT NONE

!$acc routine seq

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) :: &
  T_water                             ! Water temperature [K]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Buoyancy parameter [m s^{-2} K^{-1}]

  flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_buoypar

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowdensity (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow density,
!  using an empirical approximation from Heise et al. (2003).
!  
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_rho_w_r               , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_S_min             , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max             , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S           , & ! Empirical parameter [kg m^{-4}] in the expression for the snow density
  c_small_flk                   ! A small number

!==============================================================================

IMPLICIT NONE

!$acc routine seq

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) :: &
  h_snow                              ! Snow thickness [m]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Snow density [kg m^{-3}]

!  Security. Ensure that the expression in () does not become negative at a very large h_snow.
  flake_snowdensity = MAX( c_small_flk, (1._ireals - h_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
  flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowdensity 

! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowheatconduct (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow heat conductivity,
!  using an empirical approximation from Heise et al. (2003).
!  
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_rho_w_r               , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_kappa_S_min           , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max           , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S             ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] in the expression for kappa_S

!==============================================================================

IMPLICIT NONE

!$acc routine seq

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = ireals), INTENT(IN) :: &
  h_snow                              ! Snow thickness [m]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

  flake_snowheatconduct = flake_snowdensity( h_snow )   ! Compute snow density
  flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
                        + h_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowheatconduct


END MODULE flake_core
