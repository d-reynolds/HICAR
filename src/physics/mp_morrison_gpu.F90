!>----------------------------------------------------------------------
!! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!!     MORRISON ET AL. (2009, MWR)
!!
!! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1
!!
!! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
!! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!!     OF VERLINDE AND COTTON (1993)
!! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!!     IN LOW REFLECTIVITY REGIONS
!! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
!! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)
!! 6) INCLUDE WRF CONSTANTS PER SUGGESTION OF JIMY
!!
!! bug fix, 5/12/10
!! 7) bug fix for saturation vapor pressure in low pressure, to avoid division by zero
!! 8) include 'EP2' WRF constant for saturation mixing ratio calculation, instead of hardwire constant
!!
!! CHANGES FOR V3.3
!! 1) MODIFICATION FOR COUPLING WITH WRF-CHEM (PREDICTED DROPLET NUMBER CONCENTRATION) AS AN OPTION
!! 2) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
!! 3) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
!! 4) CLEAN UP OF COMMENTS IN THE CODE
!!
!! additional minor bug fixes and small changes, 5/30/2011
!! minor revisions by A. Ackerman April 2011:
!! 1) replaced kinematic with dynamic viscosity
!! 2) replaced scaling by air density for cloud droplet sedimentation
!!    with viscosity-dependent Stokes expression
!! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
!! 4) corrected typo in 2nd digit of ventilation constant F2R
!!
!! additional fixes:
!! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
!! 6) NPRACS IS NOT SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS
!! 7) fix for switch for running w/o graupel/hail (cloud ice and snow only)
!!
!! hm bug fix 3/16/12
!!
!! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted
!!
!! WRFV3.5
!! hm/A. Ackerman bug fix 11/08/12
!!
!! 1) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass
!!    collected by rain
!! 2) minor changes to some comments
!! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!!    cm-3. This was done to address the problem of excessive and persistent
!!    anvil cirrus produced by the scheme.
!!
!! CHANGES FOR WRFV3.5.1
!! 1) added output for snow+cloud ice and graupel time step and accumulated
!!    surface precipitation
!! 2) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!!    Hailong Wang (PNNL)
!! 3) very minor fix to immersion freezing rate formulation (negligible impact)
!! 4) clarifications to code comments
!! 5) minor change to shedding of rain, remove limit so that the number of
!!    collected drops can smaller than number of shed drops
!! 6) change of specific heat of liquid water from 4218 to 4187 J/kg/K
!!
!! CHANGES FOR WRFV3.6.1
!! 1) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!!    from the calculation of PSMLT and PGMLT
!! 2) redundant initialization of PSMLT (non answer-changing)
!!
!! ----------------------------------------------------------------------
!!
!! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
!! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
!! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.
!!----------------------------------------------------------------------
!
!WRF:MODEL_LAYER:PHYSICS
!



MODULE MODULE_MP_MORR_TWO_MOMENT_gpu
   ! USE     module_wrf_error
   ! USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
   ! USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT
   ! USE module_mp_radar

! USE WRF PHYSICS CONSTANTS taken from ICAR data_structures module
  use mod_wrf_constants, ONLY: CP=>cp, G=>gravity, R => r_d, RV => r_v, EP_2
!  USE module_state_description

   IMPLICIT NONE
   private
   PUBLIC  ::  MP_MORR_TWO_MOMENT_gpu, MORR_TWO_MOMENT_INIT_gpu

   REAL, PARAMETER :: PI = 3.1415926535897932384626434
   REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297

   PUBLIC  ::  POLYSVP

   PRIVATE :: GAMMA, DERF1
   PRIVATE :: PI, SQRTPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! IACT = 3, ACTIVATION CALCULATED IN MODULE_MIXACTIVATE

     INTEGER, PRIVATE ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION
! !!!NOTE: PREDICTED DROPLET CONCENTRATION NOT AVAILABLE IN THIS VERSION
! CONTACT HUGH MORRISON (morrison@ucar.edu) FOR FURTHER INFORMATION

     INTEGER, PRIVATE ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     REAL, PRIVATE ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     INTEGER, PRIVATE ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS

     INTEGER, PRIVATE ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
!             NON-EQULIBRIUM SUPERSATURATION,
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

     INTEGER, PRIVATE ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

     INTEGER, PRIVATE ::  ISUB

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     INTEGER, PRIVATE ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     INTEGER, PRIVATE ::  IHAIL

! CLOUD MICROPHYSICS CONSTANTS

     REAL, PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
!     REAL, PRIVATE ::      R           ! GAS CONSTANT FOR AIR
!     REAL, PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
!     REAL, PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     REAL, PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     REAL, PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
     REAL, PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     REAL, PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
     REAL, PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
     REAL, PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     REAL, PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     REAL, PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     REAL, PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
     REAL, PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
!     REAL, PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
     REAL, PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     REAL, PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     REAL, PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     REAL, PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     REAL, PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! hm, add for V3.2
     REAL, PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     REAL, PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
     REAL, PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     REAL, PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     REAL, PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
     REAL, PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
     REAL, PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     REAL, PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     REAL, PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     REAL, PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     REAL, PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
     REAL, PRIVATE ::      BACT    ! ACTIVATION PARAMETER
     REAL, PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     REAL, PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     REAL, PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     REAL, PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     REAL, PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     REAL, PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     REAL, PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     REAL, PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING

! CONSTANTS TO IMPROVE EFFICIENCY

     REAL, PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     REAL, PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     REAL, PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     REAL, PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     REAL, PRIVATE :: CONS41
      !$acc declare create(PI, SQRTPI, QSMALL,RHOW, dcs, cons1, cons12, &
      !$acc   IACT, INUM, NDCNST, ILIQ, INUC, IBASE, ISUB, IGRAUP, IHAIL, &
      !$acc   AI, AC, AS, AR, AG, BI, BC, BS, BR, BG, RHOSU, RHOI, RHOSN, RHOG, &
      !$acc   AIMM, BIMM, ECR, MI0, MG0, F1S, F2S, F1R, F2R, EII, ECI, RIN, CPW, &
      !$acc   CI, DI, CS, DS, CG, DG, C1, K1, MW, OSM, VI, EPSM, RHOA, MAP, MA, RR, &
      !$acc   BACT, RM1, RM2, NANEW1, NANEW2, SIG1, SIG2, F11, F12, F21, F22, MMULT, &
      !$acc   LAMMAXI, LAMMINI, LAMMAXR, LAMMINR, LAMMAXS, LAMMINS, LAMMAXG, LAMMING, &
      !$acc   CONS2, CONS3, CONS4, CONS5, CONS6, CONS7, CONS8, CONS9, CONS10, &
      !$acc   CONS11, CONS13, CONS14, CONS15, CONS16, CONS17, CONS18, CONS19, CONS20, &
      !$acc   CONS21, CONS22, CONS23, CONS24, CONS25, CONS26, CONS27, CONS28, CONS29, &
      !$acc   CONS30, CONS31, CONS32, CONS33, CONS34, CONS35, CONS36, CONS37, CONS38, &
      !$acc   CONS39, CONS40, CONS41)


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MORR_TWO_MOMENT_INIT_gpu(hail_opt) ! RAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      INTEGER, INTENT(IN):: hail_opt ! RAS

      integer n,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM IS AUTOMATICALLY SET TO 0 FOR WRF-CHEM BELOW,
! ALLOWING PREDICTION OF DROPLET CONCENTRATION
! THUS, THIS PARAMETER SHOULD NOT BE CHANGED HERE
! AND SHOULD BE LEFT TO 1

      INUM = 1

! SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)
! IF NO COUPLING WITH WRF-CHEM

      NDCNST = 250.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE, THE FOLLOWING OPTIONS RELATED TO DROPLET ACTIVATION
! (IACT, IBASE, ISUB) ARE NOT AVAILABLE IN CURRENT VERSION
! FOR WRF-CHEM, DROPLET ACTIVATION IS PERFORMED
! IN 'MIX_ACTIVATE', NOT IN MICROPHYSICS SCHEME


! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

      IACT = 2

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER,
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      IBASE = 2

! INCLUDE SUB-GRID VERTICAL VELOCITY (standard deviation of w) IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! currently, sub-grid w is constant of 0.5 m/s (not coupled with PBL/turbulence scheme)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      ISUB = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      ILIQ = 0

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      INUC = 0

! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
! IGRAUP = 1, NO GRAUPEL/HAIL

      IGRAUP = 0

! HM ADDED 11/7/07
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

      !IHAIL = 0 !changed to namelist option (hail_opt) by RAS
      ! Check if namelist option is feasible, otherwise default to graupel - RAS
      IF (hail_opt .eq. 1) THEN
         IHAIL = 1
      ELSE
         IHAIL = 0
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
         IF (IHAIL.EQ.0) THEN
	 AG = 19.3
	 BG = 0.37
         ELSE ! (MATSUN AND HUGGINS 1980)
         AG = 114.5
         BG = 0.5
         END IF

! CONSTANTS AND PARAMETERS
!         R = 287.15
!         RV = 461.5
!         CP = 1005.
         RHOSU = 85000./(287.15*273.15)
         RHOW = 997.
         RHOI = 500.
         RHOSN = 100.
         IF (IHAIL.EQ.0) THEN
	 RHOG = 400.
         ELSE
         RHOG = 900.
         END IF
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
	 MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
!         F2R = 0.32
! fix 053011
         F2R = 0.308
!         G = 9.806
         QSMALL = 1.E-14
         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
! hm, 7/23/13
!         CPW = 4218.
         CPW = 4187.

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA

         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
         LAMMINR = 1./2800.E-6

         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = 0.4
              C1 = 120.

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000.

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
         RR = 8.3187
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE
! (see morrison et al. 2007, JGR)
! MODE 1

         RM1 = 0.052E-6
         SIG1 = 2.04
         NANEW1 = 72.2E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = 1.3E-6
         SIG2 = 2.5
         NANEW2 = 1.8E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)

! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

!+---+-----------------------------------------------------------------+
!..Set these variables needed for computing radar reflectivity.  These
!.. get used within radar_init to create other variables used in the
!.. radar module.

        !  xam_r = PI*RHOW/6.
        !  xbm_r = 3.
        !  xmu_r = 0.
        !  xam_s = CS
        !  xbm_s = DS
        !  xmu_s = 0.
        !  xam_g = CG
        !  xbm_g = DG
        !  xmu_g = 0.
         !
        !  call radar_init
!+---+-----------------------------------------------------------------+

      !$acc update device(QSMALL,RHOW, dcs, cons1, cons12, &
      !$acc   IACT, INUM, NDCNST, ILIQ, INUC, IBASE, ISUB, IGRAUP, IHAIL, &
      !$acc   AI, AC, AS, AR, AG, BI, BC, BS, BR, BG, RHOSU, RHOI, RHOSN, RHOG, &
      !$acc   AIMM, BIMM, ECR, MI0, MG0, F1S, F2S, F1R, F2R, EII, ECI, RIN, CPW, &
      !$acc   CI, DI, CS, DS, CG, DG, C1, K1, MW, OSM, VI, EPSM, RHOA, MAP, MA, RR, &
      !$acc   BACT, RM1, RM2, NANEW1, NANEW2, SIG1, SIG2, F11, F12, F21, F22, MMULT, &
      !$acc   LAMMAXI, LAMMINI, LAMMAXR, LAMMINR, LAMMAXS, LAMMINS, LAMMAXG, LAMMING, &
      !$acc   CONS2, CONS3, CONS4, CONS5, CONS6, CONS7, CONS8, CONS9, CONS10, &
      !$acc   CONS11, CONS13, CONS14, CONS15, CONS16, CONS17, CONS18, CONS19, CONS20, &
      !$acc   CONS21, CONS22, CONS23, CONS24, CONS25, CONS26, CONS27, CONS28, CONS29, &
      !$acc   CONS30, CONS31, CONS32, CONS33, CONS34, CONS35, CONS36, CONS37, CONS38, &
      !$acc   CONS39, CONS40, CONS41)

END SUBROUTINE MORR_TWO_MOMENT_INIT_gpu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO)
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE.
! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MP_MORR_TWO_MOMENT_gpu(ITIMESTEP,                       &
                TH, QV, QC, QR, QI, QS, QG, NI, NS, NR, NG, &
                RHO_IN, PII, P, DT_IN, DZ, W,          &
                RAINNC, RAINNCV, SR,                    &
                SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,    & ! hm added 7/13/13
                EFFC, EFFI, EFFS,                       & ! particle radiuses for radiation
                refl_10cm, diagflag, do_radar_ref,      & ! GT added for reflectivity calcs
                qrcuten, qscuten, qicuten & ! mu           & ! hm added
!               ,F_QNDROP, qndrop                        & ! hm added, wrf-chem
               ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
!jdf		   ,C2PREC3D,CSED3D,
		,ISED3D,SSED3D                                 &
!               ,GSED3D,RSED3D & ! HM ADD, WRF-CHEM
!               ,rainprod, evapprod                      &
!		   ,QLSINK,PRECR,PRECI,PRECS,PRECG &        ! HM ADD, WRF-CHEM
                                            )

! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SNOWNC - accumulated grid-scale snow plus cloud ice (mm)
! SNOWNCV - one time step grid scale snow plus cloud ice (mm/time step)
! GRAUPELNC - accumulated grid-scale graupel (mm)
! GRAUPELNCV - one time step grid scale graupel (mm/time step)
! SR - one time step mass ratio of snow to total precip
! qrcuten, rain tendency from parameterized cumulus convection
! qscuten, snow tendency from parameterized cumulus convection
! qicuten, cloud ice tendency from parameterized cumulus convection

! variables below currently not in use, not coupled to PBL or radiation codes
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! HM, ADDED FOR WRF-CHEM COUPLING
! QLSINK - TENDENCY OF CLOUD WATER TO RAIN, SNOW, GRAUPEL (KG/KG/S)
! CSED,ISED,SSED,GSED,RSED - SEDIMENTATION FLUXES (KG/M^2/S) FOR CLOUD WATER, ICE, SNOW, GRAUPEL, RAIN
! PRECI,PRECS,PRECG,PRECR - SEDIMENTATION FLUXES (KG/M^2/S) FOR ICE, SNOW, GRAUPEL, RAIN

! rainprod - total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)
! evapprod - tendency of evaporation of rain (kg kg-1 s-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reflectivity currently not included!!!!
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte
! Temporary changed from INOUT to IN

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qv, qc, qr, qi, qs, qg, ni, ns, nr, TH, NG, EFFC, EFFI, EFFS, SSED3D, ISED3D
!jdf                      qndrop ! hm added, wrf-chem
   ! REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: qndrop
!jdf  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT):: CSED3D, &
   ! REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: QLSINK, &
   !                        rainprod, evapprod, &
   !                        PRECI,PRECS,PRECG,PRECR ! HM, WRF-CHEM
!, effcs, effis

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz, RHO_IN, w !, tke, nctend, nitend,kzh
   REAL, INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: ITIMESTEP

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR, &
! hm added 7/13/13
                          SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
                          refl_10cm

   ! REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

   ! LOCAL VARIABLES

   REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                      EFFR, EFFG

! add cumulus tendencies

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
      qrcuten, qscuten, qicuten
   ! REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN):: &
   !    mu

  ! LOGICAL, INTENT(IN), OPTIONAL ::                F_QNDROP  ! wrf-chem
  LOGICAL :: flag_qndrop  ! wrf-chem
  integer :: iinum ! wrf-chem

! wrf-chem
   REAL, DIMENSION(kts:kte) :: nc1d, nc_tend1d
   !REAL, DIMENSION(kts:kte) :: rainprod1d, evapprod1d
! HM add reflectivity
   !REAL, DIMENSION(kts:kte) :: dBZ

   REAL, dimension(its:ite,jts:jte) :: PRECPRT1D, SNOWRT1D, SNOWPRT1D, GRPLPRT1D ! hm added 7/13/13

   INTEGER I,K,J

   REAL DT

   LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
   INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!START OF DECLARATIONS FROM INLINED FUNCTION!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)

      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  T3D                ! TEMPERATURE (I,K,J)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
! below for wrf-chem
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  nc3d
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  nc3dten

! HM ADDED GRAUPEL VARIABLES
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)

! hm add cumulus tendencies for precip
        REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   qrcu1d
        REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   qscu1d
        REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   qicu1d

! OUTPUT VARIABLES

      !   REAL, DIMENSION(KTS:KTE) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
      !   REAL, DIMENSION(KTS:KTE) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
      !   REAL, DIMENSION(KTS:KTE) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
      !   REAL, DIMENSION(KTS:KTE) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
      !   REAL, DIMENSION(KTS:KTE) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)
	  REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   C2PREC,CSED,ISED,SSED,GSED,RSED

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

      !   REAL DT         ! MODEL TIME STEP (SEC)

!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DRITS:ITE,KTS:KTE,JTS:JTE
! MICROPHYSICAL PRITS:ITE,KTS:KTE,JTS:JTE
     REAL ::  NSUBC     ! LOSS OF NC DURING EVAP
     REAL ::  NSUBI     ! LOSS OF NI DURING SUB.
     REAL ::  NSUBS     ! LOSS OF NS DURING SUB.
     REAL ::  NSUBR     ! LOSS OF NR DURING EVAP
     REAL ::  PRD       ! DEP CLOUD ICE
     REAL ::  PRE       ! EVAP OF RAIN
     REAL ::  PRDS      ! DEP SNOW
     REAL ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     REAL ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     REAL ::  PRA       ! ACCRETION DROPLETS BY RAIN
     REAL ::  PRC       ! AUTOCONVERSION DROPLETS
     REAL ::  PCC       ! COND/EVAP DROPLETS
     REAL ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     REAL ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
     REAL ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     REAL ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
     REAL ::  NSAGG     ! SELF-COLLECTION OF SNOW
     REAL ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     REAL ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
     REAL ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE BY SNOW
     REAL ::  PRCI      ! CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
     REAL ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     REAL ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     REAL ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     REAL ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     REAL ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     REAL ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     REAL ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     REAL ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
     REAL ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     REAL ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
     REAL ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     REAL ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     REAL ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     REAL ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
     REAL ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     REAL ::  NSMLTS    ! CHANGE N MELTING SNOW
     REAL ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
! HM ADDED 12/13/06
     REAL ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     REAL ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     REAL ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     REAL ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     REAL ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     REAL ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     REAL ::  EPRD      ! SUBLIMATION CLOUD ICE
     REAL ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED 12/13/06
     REAL ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     REAL ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     REAL ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL ::  PRDG    ! DEP OF GRAUPEL
     REAL ::  EPRDG    ! SUB OF GRAUPEL
     REAL ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     REAL ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     REAL ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     REAL ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     REAL ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     REAL ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
     REAL ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     REAL ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     REAL ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     REAL ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     REAL ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     REAL ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
! TIME-VARYING ATMOSPHERIC PARAMETERS

     REAL ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     REAL ::   EVS   ! SATURATION VAPOR PRESSURE
     REAL ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     REAL ::   QVS   ! SATURATION MIXING RATIO
     REAL ::   QVI   ! ICE SATURATION MIXING RATIO
     REAL ::   QVQVS ! SAUTRATION RATIO
     REAL ::   QVQVSI! ICE SATURAION RATIO
     REAL ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   MU    ! VISCOCITY OF AIR
     REAL ::   SC    ! SCHMIDT NUMBER
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   XLF   ! LATENT HEAT OF FREEZING
     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE)  ::   RHO   ! AIR DENSITY
     REAL ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
     REAL ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     REAL    NACNT                    ! NUMBER OF CONTACT IN
     REAL    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
     REAL    COFFI                    ! ICE AUTOCONVERSION PARAMETER

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      REAL UNI, UMI,UMR
      REAL RGVM
      REAL FALTNDR,FALTNDI,FALTNDNI,RHO2
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   DUMQS,DUMFNS
      REAL UMS,UNS
      REAL FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::    DUMC,DUMFNC
      REAL UNC,UMC,UNG,UMG
      REAL FALTNDC,FALTNDNC
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   DUMFNR,FALOUTNR
      REAL FALTNDNR

      ! the only arrays which need to be stored
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::    FR, FI, FC, FNI,FG,FNG, FS,FNS, FNR, FNC
      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::   FALOUTS,FALOUTNS,FALOUTG,FALOUTNG, FALOUTR,FALOUTI,FALOUTNI, FALOUTC,FALOUTNC

! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      REAL, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) ::    AIN,ARN,ASN,ACN,AGN

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL GAMMA,      ! EULER GAMMA FUNCTION
!      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     REAL DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

! PROGNOSTIC SUPERSATURATION

     REAL DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     REAL DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
     REAL EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     REAL EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     REAL EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     REAL EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     REAL TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     REAL TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     REAL TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     REAL TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     REAL TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     REAL DUMACT,DUM3

! COUNTING/INDEX VARIABLES

     INTEGER N, MAXN ! ,I
     INTEGER, DIMENSION(ITS:ITE,KMS:KME,JTS:JTE) ::   NSTEP ! NUMBER OF CLOUD DROPLETS IN CELL
     INTEGER, DIMENSION(ITS:ITE,JTS:JTE) ::   NSTEP_FLAT ! NUMBER OF CLOUD DROPLETS IN CELL

! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
! LTRUE, SWITCH = 0, NO HYDROMETEORS IN CELL,
!               = 1, HYDROMETEORS IN CELL
     INTEGER, DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE)  ::   LTRUE, LTRUE_COL   ! AIR DENSITY

! DROPLET ACTIVATION/FREEZING AEROSOL


     REAL    CT      ! DROPLET ACTIVATION PARAMETER
     REAL    TEMP1   ! DUMMY TEMPERATURE
     REAL    SAT1    ! DUMMY SATURATION
     REAL    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     REAL    KEL     ! KELVIN PARAMETER
     REAL    KC2     ! TOTAL ICE NUCLEATION RATE

       REAL CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     REAL DUMQI,DUMNI,DC0,DS0,DG0
     REAL DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     REAL WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      REAL ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

        REAL AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

! DUMMY SIZE DISTRIBUTION PARAMETERS

        REAL DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

        INTEGER IDROP

! FOR WRF-CHEM
! #if (WRF_CHEM == 1)
!     REAL, DIMENSION(KTS:KTE), INTENT(INOUT) :: rainprod, evapprod
! #endif
    REAL, DIMENSION(KTS:KTE)                :: tqimelt ! melting of cloud ice (tendency)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS FROM INLINED FUNCTION!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! below for wrf-chem
   flag_qndrop = .false.
   ! IF ( PRESENT ( f_qndrop ) ) flag_qndrop = f_qndrop
!!!!!!!!!!!!!!!!!!!!!!

   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN

   iinum=1

   !$acc data present(TH, QV, QC, QR, QI, QS, QG, NI, NS, NR, NG, &
   !$acc             RHO_IN, PII, P, DZ, W,              &
   !$acc             RAINNC, RAINNCV, SR,                    &
   !$acc             SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,    & 
   !$acc             EFFC, EFFI, EFFS,                       &
   !$acc             qrcuten, qscuten, qicuten,              & 
   !$acc 		   ISED3D,SSED3D)                          &
   !$acc      create(EFFR, EFFG, &
   !$acc             QC3DTEN, QI3DTEN, QNI3DTEN, QR3DTEN, &
   !$acc             NI3DTEN, NS3DTEN, NR3DTEN, &
   !$acc             QC3D, QI3D, QR3D, NI3D, NS3D, NR3D, QNI3D, &
   !$acc             T3DTEN, QV3DTEN, T3D, QV3D, PRES, &
   !$acc             DZQ, RHO, &
   !$acc             QG3DTEN, NG3DTEN, QG3D, NG3D, &
   !$acc             qrcu1d, qscu1d, qicu1d, &
   !$acc             QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN, LTRUE_COL, NSTEP, NSTEP_FLAT, &
   !$acc             nc1d, nc_tend1d, C2PREC,CSED,ISED,SSED,GSED,RSED, &
   !$acc             lamg,acn,abi,arn,ain,agn,ltrue,n0s,qvqvs,qvi,pgam, &
   !$acc             qvqvsi,dv,eis,cdist1,xlf,xxlv,xxls,nc3d,lams,sc,asn, &
   !$acc             n0g,cpm,lamr,n0rr,lami,n0i,kap,nc3dten,mu,lamc,ab,qvs,dap,evs, &
   !$acc             precprt1d,snowrt1d,snowprt1d,grplprt1d, &
   !$acc             fnc,fc,fi,dumr,fni,fnr,fs,fr,fng,fg,fns,dumg, &
   !$acc             dumi,dumqs,dumfng,dumfnr,dumfni,dumc,dumfns,dumfnc, &
   !$acc             faloutns,faloutg,faloutng,falouts,faloutr,faloutnr,faloutni,falouti,faloutc,faloutnc)

   ! default firstprivate to handle all of the module level variables that are defined in init (I hope this works as intented...)
   ! !$omp parallel default(shared) &
   ! !$omp private(i,j,k) &
   ! !$omp private(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,NI3DTEN,NS3DTEN,NR3DTEN,T3DTEN,QV3DTEN,nc_tend1d) &
   ! !$omp private(QC3D,QI3D,QNI3D,QR3D,NI3D,NS3D,NR3D,QG3D,NG3D,QG3DTEN,NG3DTEN) &
   ! !$omp private(T3D,QV3D,PRES,DZQ,qrcu1d,qscu1d,qicu1d,nc1d,iinum) &
   ! !$omp private(PRECPRT1D,SNOWPRT1D,GRPLPRT1D) &
   ! !$omp private(QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN) &
   ! !$omp shared(QC,QI,QS,QR,NI,NS,NR,QG,NG,T,QV,P,DZ,W,qrcuten,qscuten,qicuten) &
   ! !$omp shared(TH,EFFC,EFFI,EFFS,EFFR,EFFG,RAINNC,RAINNCV,SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,SR) &
   ! !$omp firstprivate(IHAIL,IGRAUP,ISUB,IBASE,INUC,ILIQ,NDCNST,INUM,IACT,DT) &
   ! !$omp firstprivate(its, ite, jts, jte, kts,kte)
   ! !$omp do schedule(dynamic)

   !$acc parallel loop gang collapse(2)
   do j=jts,jte      ! j loop (north-south)
   do i=its,ite      ! i loop (east-west)
          ! INITIALIZE PRECIP AND SNOW RATES
          PRECPRT1D(i,j) = 0.
          SNOWRT1D(i,j) = 0.
          ! hm added 7/13/13
          SNOWPRT1D(i,j) = 0.
          GRPLPRT1D(i,j) = 0.
    enddo
    enddo

   !$acc parallel loop gang vector collapse(3)
   do j=jts,jte      ! j loop (north-south)
   do k=kts,kte   ! k loop (vertical)
   do i=its,ite      ! i loop (east-west)

          QC3DTEN(I,K,J)  = 0.
          QI3DTEN(I,K,J)  = 0.
          QNI3DTEN(I,K,J) = 0.
          QR3DTEN(I,K,J)  = 0.
          NI3DTEN(I,K,J)  = 0.
          NS3DTEN(I,K,J)  = 0.
          NR3DTEN(I,K,J)  = 0.
          T3DTEN(I,K,J)   = 0.
          QV3DTEN(I,K,J)  = 0.
          nc_tend1d = 0. ! wrf-chem

          QC3D(I,K,J)       = QC(i,k,j)
          QI3D(I,K,J)       = QI(i,k,j)
          QNI3D(I,K,J)       = QS(i,k,j)
          QR3D(I,K,J)       = QR(i,k,j)

          NI3D(I,K,J)       = NI(i,k,j)
          NS3D(I,K,J)       = NS(i,k,j)
          NR3D(I,K,J)       = NR(i,k,j)
! HM ADD GRAUPEL
          QG3D(I,K,J)       = QG(I,K,j)
          NG3D(I,K,J)       = NG(I,K,j)
          QG3DTEN(I,K,J)  = 0.
          NG3DTEN(I,K,J)  = 0.

      !     EFFC1D(I,K,J)     = EFFC(i,k,j)
      !     EFFS1D(I,K,J)     = EFFS(i,k,j)
      !     EFFI1D(I,K,J)     = EFFI(i,k,j)
          
          T3D(I,K,J)        = TH(i,k,j)*PII(i,k,j)
          QV3D(I,K,J)       = QV(i,k,j)
          PRES(I,K,J)        = P(i,k,j)
          DZQ(I,K,J)       = DZ(i,k,j)
! add cumulus tendencies, decouple from mu
          qrcu1d(I,K,J)     = qrcuten(i,k,j) ! /mu(i,j) ! not coupled with mu in ICAR
          qscu1d(I,K,J)     = qscuten(i,k,j) ! /mu(i,j) ! not coupled with mu in ICAR
          qicu1d(I,K,J)     = qicuten(i,k,j) ! /mu(i,j) ! not coupled with mu in ICAR
          nc1d=0. ! temporary placeholder, set to constant in microphysics subroutine
      ENDDO
      ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!START OF INLINED FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ! SUBROUTINE MORR_TWO_MOMENT_MICRO(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,         &
      !  NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NI3D,NS3D,NR3D,              &
      !  T3DTEN,QV3DTEN,T3D,QV3D,PRES,DZQ,W3D,WVAR,PRECRT,SNOWRT,            &
      !  SNOWPRT,GRPLPRT,                & ! hm added 7/13/13
      !  EFFC,EFFI,EFFS,EFFR,DT,                                                   &
      !                                       IMS,IME, JMS,JME, KMS,KME,           &
      !                                       ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
      !                   QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,qrcu1d,qscu1d, qicu1d,    &
      !                   QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
      !                   nc3d,nc3dten,iinum, & ! wrf-chem
	! 			c2prec,CSED,ISED,SSED,GSED,RSED  &  ! hm added, wrf-chem
! #if (WRF_CHEM == 1)
!         ,rainprod, evapprod &
! #endif
                        ! )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS; MORRISON AND PINTO 2005 JAS.
! ADDITIONAL CHANGES ARE DESCRIBED IN DETAIL BY MORRISON, THOMPSON, TATARSKII (MWR, SUBMITTED)

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'MORR_TWO_MOMENT'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS


! comment lines for wrf-chem since these are intent(in) in that case
!       REAL, DIMENSION(KTS:KTE) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
!       REAL, DIMENSION(KTS:KTE) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)



! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
   !$acc parallel loop gang vector collapse(3) async(1)
   do j=jts,jte      ! j loop (north-south)
   DO K = KTS,KTE
   do i=its,ite      ! i loop (east-west)

! SET LTRUE INITIALLY TO 0

               LTRUE(I,K,J) = 0
               NSTEP(I,K,J) = 1

! NC3DTEN LOCAL ARRAY INITIALIZED
               NC3DTEN(I,K,J) = 0.
! INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO

		C2PREC(I,K,J)=0.
		CSED(I,K,J)=0.
		ISED(I,K,J)=0.
		SSED(I,K,J)=0.
		GSED(I,K,J)=0.
		RSED(I,K,J)=0.

! #if (WRF_CHEM == 1)
!          rainprod(I,K,J) = 0.
!          evapprod(I,K,J) = 0.
!          tqimelt(I,K,J)  = 0.
!          PRC      = 0.
!          PRA      = 0.
! #endif

! LATENT HEAT OF VAPORATION

            XXLV(I,K,J) = 3.1484E6-2370.*T3D(I,K,J)

! LATENT HEAT OF SUBLIMATION

            XXLS(I,K,J) = 3.15E6-2370.*T3D(I,K,J)+0.3337E6

            CPM(I,K,J) = CP*(1.+0.887*QV3D(I,K,J))

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS = min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J),0))   ! PA
            EIS = min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS.GT.EVS) EIS = EVS

            QVS = EP_2*EVS/(PRES(I,K,J)-EVS)
            QVI = EP_2*EIS/(PRES(I,K,J)-EIS)

            QVQVS = QV3D(I,K,J)/QVS
            QVQVSI = QV3D(I,K,J)/QVI

! AIR DENSITY

            RHO(I,K,J) = PRES(I,K,J)/(R*T3D(I,K,J))

! ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
! FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON

            IF (QRCU1D(I,K,J).GE.1.E-10) THEN
            DUM=1.8e5*(QRCU1D(I,K,J)*DT/(PI*RHOW*RHO(I,K,J)**3))**0.25
            NR3D(I,K,J)=NR3D(I,K,J)+DUM
            END IF
            IF (QSCU1D(I,K,J).GE.1.E-10) THEN
            DUM=3.e5*(QSCU1D(I,K,J)*DT/(CONS1*RHO(I,K,J)**3))**(1./(DS+1.))
            NS3D(I,K,J)=NS3D(I,K,J)+DUM
            END IF
            IF (QICU1D(I,K,J).GE.1.E-10) THEN
            DUM=QICU1D(I,K,J)*DT/(CI*(80.E-6)**DI)
            NI3D(I,K,J)=NI3D(I,K,J)+DUM
            END IF

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm modify 7/0/09 change limit to 1.e-8

             IF (QVQVS.LT.0.9) THEN
               IF (QR3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QR3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QR3D(I,K,J)*XXLV(I,K,J)/CPM(I,K,J)
                  QR3D(I,K,J)=0.
               END IF
               IF (QC3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QC3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QC3D(I,K,J)*XXLV(I,K,J)/CPM(I,K,J)
                  QC3D(I,K,J)=0.
               END IF
             END IF

             IF (QVQVSI.LT.0.9) THEN
               IF (QI3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QI3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QI3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QI3D(I,K,J)=0.
               END IF
               IF (QNI3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QNI3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QNI3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QNI3D(I,K,J)=0.
               END IF
               IF (QG3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QG3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QG3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QG3D(I,K,J)=0.
               END IF
             END IF

! HEAT OF FUSION

            XLF(I,K,J) = XXLS(I,K,J)-XXLV(I,K,J)

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(I,K,J).LT.QSMALL) THEN
         QC3D(I,K,J) = 0.
         NC3D(I,K,J) = 0.
         EFFC(I,K,J) = 0.
       END IF
       IF (QR3D(I,K,J).LT.QSMALL) THEN
         QR3D(I,K,J) = 0.
         NR3D(I,K,J) = 0.
         EFFR(I,K,J) = 0.
       END IF
       IF (QI3D(I,K,J).LT.QSMALL) THEN
         QI3D(I,K,J) = 0.
         NI3D(I,K,J) = 0.
         EFFI(I,K,J) = 0.
       END IF
       IF (QNI3D(I,K,J).LT.QSMALL) THEN
         QNI3D(I,K,J) = 0.
         NS3D(I,K,J) = 0.
         EFFS(I,K,J) = 0.
       END IF
       IF (QG3D(I,K,J).LT.QSMALL) THEN
         QG3D(I,K,J) = 0.
         NG3D(I,K,J) = 0.
         EFFG(I,K,J) = 0.
       END IF

! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

      QRSTEN(I,K,J) = 0.
      QISTEN(I,K,J) = 0.
      QNISTEN(I,K,J) = 0.
      QCSTEN(I,K,J) = 0.
      QGSTEN(I,K,J) = 0.

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

! fix 053011
            MU(I,K,J) = 1.496E-6*T3D(I,K,J)**1.5/(T3D(I,K,J)+120.)

! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(I,K,J))**0.54

! fix 053011
!            AIN(I,K,J) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction
            AIN(I,K,J) = (RHOSU/RHO(I,K,J))**0.35*AI
            ARN(I,K,J) = DUM*AR
            ASN(I,K,J) = DUM*AS
!            ACN(I,K,J) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(I,K,J) = G*RHOW/(18.*MU(I,K,J))
! HM ADD GRAUPEL 8/28/06
            AGN(I,K,J) = DUM*AG

!hm 4/7/09 bug fix, initialize lami to prevent later division by zero
            LAMI(I,K,J)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

      IF (.not. (QC3D(I,K,J).LT.QSMALL.AND.QI3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.QSMALL &
            .AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.QSMALL.AND.( &
            (T3D(I,K,J).LT.273.15.AND.QVQVSI.LT.0.999).OR. &
            (T3D(I,K,J).GE.273.15.AND.QVQVS.LT.0.999)))) THEN

! THERMAL CONDUCTIVITY FOR AIR

! fix 053011
            KAP = 1.414E3*MU(I,K,J)

! DIFFUSIVITY OF WATER VAPOR

            DV = 8.794E-5*T3D(I,K,J)**1.81/PRES(I,K,J)

! SCHMIT NUMBER

! fix 053011
            SC = MU(I,K,J)/(RHO(I,K,J)*DV)

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(I,K,J)**2)

            DQSDT = XXLV(I,K,J)*QVS/DUM
            DQSIDT =  XXLS(I,K,J)*QVI/DUM

            ABI = 1.+DQSIDT*XXLS(I,K,J)/CPM(I,K,J)
            AB = 1.+DQSDT*XXLV(I,K,J)/CPM(I,K,J)

!
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

            IF (T3D(I,K,J).GE.273.15) THEN

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(I,K,J)=NDCNST*1.E6/RHO(I,K,J)
         END IF

! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
       IF (QNI3D(I,K,J).LT.1.E-6) THEN
          QR3D(I,K,J)=QR3D(I,K,J)+QNI3D(I,K,J)
          NR3D(I,K,J)=NR3D(I,K,J)+NS3D(I,K,J)
          T3D(I,K,J)=T3D(I,K,J)-QNI3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
          QNI3D(I,K,J) = 0.
          NS3D(I,K,J) = 0.
       END IF
       IF (QG3D(I,K,J).LT.1.E-6) THEN
          QR3D(I,K,J)=QR3D(I,K,J)+QG3D(I,K,J)
          NR3D(I,K,J)=NR3D(I,K,J)+NG3D(I,K,J)
          T3D(I,K,J)=T3D(I,K,J)-QG3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
          QG3D(I,K,J) = 0.
          NG3D(I,K,J) = 0.
       END IF

      IF (.not. (QC3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.1.E-8.AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.1.E-8)) THEN

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NS3D(I,K,J) = MAX(0.,NS3D(I,K,J))
      NC3D(I,K,J) = MAX(0.,NC3D(I,K,J))
      NR3D(I,K,J) = MAX(0.,NR3D(I,K,J))
      NG3D(I,K,J) = MAX(0.,NG3D(I,K,J))

!......................................................................
! RAIN

      IF (QR3D(I,K,J).GE.QSMALL) THEN
      LAMR(I,K,J) = (PI*RHOW*NR3D(I,K,J)/QR3D(I,K,J))**(1./3.)
      N0RR(I,K,J) = NR3D(I,K,J)*LAMR(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(I,K,J).LT.LAMMINR) THEN

      LAMR(I,K,J) = LAMMINR

      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      ELSE IF (LAMR(I,K,J).GT.LAMMAXR) THEN
      LAMR(I,K,J) = LAMMAXR
      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(I,K,J).GE.QSMALL) THEN

         DUM = PRES(I,K,J)/(287.15*T3D(I,K,J))
         PGAM(I,K,J)=0.0005714*(NC3D(I,K,J)/1.E6*DUM)+0.2714
         PGAM(I,K,J)=1./(PGAM(I,K,J)**2)-1.
         PGAM(I,K,J)=MAX(PGAM(I,K,J),2.)
         PGAM(I,K,J)=MIN(PGAM(I,K,J),10.)

! CALCULATE LAMC

      LAMC(I,K,J) = (CONS26*NC3D(I,K,J)*GAMMA(PGAM(I,K,J)+4.)/   &
                 (QC3D(I,K,J)*GAMMA(PGAM(I,K,J)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(I,K,J)+1.)/60.E-6
      LAMMAX = (PGAM(I,K,J)+1.)/1.E-6

      IF (LAMC(I,K,J).LT.LAMMIN) THEN
      LAMC(I,K,J) = LAMMIN

      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26
      ELSE IF (LAMC(I,K,J).GT.LAMMAX) THEN
      LAMC(I,K,J) = LAMMAX

      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(I,K,J).GE.QSMALL) THEN
      LAMS(I,K,J) = (CONS1*NS3D(I,K,J)/QNI3D(I,K,J))**(1./DS)
      N0S(I,K,J) = NS3D(I,K,J)*LAMS(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(I,K,J).LT.LAMMINS) THEN
      LAMS(I,K,J) = LAMMINS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1

      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)

      ELSE IF (LAMS(I,K,J).GT.LAMMAXS) THEN

      LAMS(I,K,J) = LAMMAXS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1

      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(I,K,J).GE.QSMALL) THEN
      LAMG(I,K,J) = (CONS2*NG3D(I,K,J)/QG3D(I,K,J))**(1./DG)
      N0G(I,K,J) = NG3D(I,K,J)*LAMG(I,K,J)

! ADJUST VARS

      IF (LAMG(I,K,J).LT.LAMMING) THEN
      LAMG(I,K,J) = LAMMING
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)

      ELSE IF (LAMG(I,K,J).GT.LAMMAXG) THEN

      LAMG(I,K,J) = LAMMAXG
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)
      END IF
      END IF

!.....................................................................
! ZERO OUT PROCESS RATES

            PRC = 0.
            NPRC = 0.
            NPRC1 = 0.
            PRA = 0.
            NPRA = 0.
            NRAGG = 0.
            NSMLTS = 0.
            NSMLTR = 0.
            EVPMS = 0.
            PCC = 0.
            PRE = 0.
            NSUBC = 0.
            NSUBR = 0.
            PRACG = 0.
            NPRACG = 0.
            PSMLT = 0.
            PGMLT = 0.
            EVPMG = 0.
            PRACS = 0.
            NPRACS = 0.
            NGMLTG = 0.
            NGMLTR = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(I,K,J).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC=1350.*QC3D(I,K,J)**2.47*  &
           (NC3D(I,K,J)/1.e6*RHO(I,K,J))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1 = PRC/CONS29
        NPRC = PRC/(QC3D(I,K,J)/NC3D(I,K,J))

! hm bug fix 3/20/12
                NPRC = MIN(NPRC,NC3D(I,K,J)/DT)
                NPRC1 = MIN(NPRC1,NPRC)

         END IF

!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
! FORMULA FROM IKAWA AND SAITO (1991)

         IF (QR3D(I,K,J).GE.1.E-8.AND.QNI3D(I,K,J).GE.1.E-8) THEN

            UMS = ASN(I,K,J)*CONS3/(LAMS(I,K,J)**BS)
            UMR = ARN(I,K,J)*CONS4/(LAMR(I,K,J)**BR)
            UNS = ASN(I,K,J)*CONS5/LAMS(I,K,J)**BS
            UNR = ARN(I,K,J)*CONS6/LAMR(I,K,J)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/RHO(I,K,J))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! hm fix, 2/12/13
! for above freezing conditions to get accelerated melting of snow,
! we need collection of rain by snow (following Lin et al. 1983)
!            PRACS = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
!                  0.08*UMS*UMR)**0.5*RHO(I,K,J)*                     &
!                 N0RR(I,K,J)*N0S(I,K,J)/LAMS(I,K,J)**3*                    &
!                  (5./(LAMS(I,K,J)**3*LAMR(I,K,J))+                    &
!                  2./(LAMS(I,K,J)**2*LAMR(I,K,J)**2)+                  &
!                  0.5/(LAMS(I,K,J)*LAMR(I,K,J)**3)))

            PRACS = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(I,K,J)*                      &
                  N0RR(I,K,J)*N0S(I,K,J)/LAMR(I,K,J)**3*                              &
                  (5./(LAMR(I,K,J)**3*LAMS(I,K,J))+                    &
                  2./(LAMR(I,K,J)**2*LAMS(I,K,J)**2)+                  &
                  0.5/(LAMR(I,K,J)*LAMS(I,K,J)**3)))

! fix 053011, npracs no longer subtracted from snow
!            NPRACS = CONS32*RHO(I,K,J)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(I,K,J)*N0S(I,K,J)*              &
!                (1./(LAMR(I,K,J)**3*LAMS(I,K,J))+                      &
!                 1./(LAMR(I,K,J)**2*LAMS(I,K,J)**2)+                   &
!                 1./(LAMR(I,K,J)*LAMS(I,K,J)**3))

         END IF

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

         IF (QR3D(I,K,J).GE.1.E-8.AND.QG3D(I,K,J).GE.1.E-8) THEN

            UMG = AGN(I,K,J)*CONS7/(LAMG(I,K,J)**BG)
            UMR = ARN(I,K,J)*CONS4/(LAMR(I,K,J)**BR)
            UNG = AGN(I,K,J)*CONS8/LAMG(I,K,J)**BG
            UNR = ARN(I,K,J)*CONS6/LAMR(I,K,J)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/RHO(I,K,J))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(I,K,J)*                      &
                  N0RR(I,K,J)*N0G(I,K,J)/LAMR(I,K,J)**3*                              &
                  (5./(LAMR(I,K,J)**3*LAMG(I,K,J))+                    &
                  2./(LAMR(I,K,J)**2*LAMG(I,K,J)**2)+				   &
				  0.5/(LAMR(I,K,J)*LAMG(I,K,J)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC

            DUM = PRACG/5.2E-7

            NPRACG = CONS32*RHO(I,K,J)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(I,K,J)*N0G(I,K,J)*              &
                (1./(LAMR(I,K,J)**3*LAMG(I,K,J))+                      &
                 1./(LAMR(I,K,J)**2*LAMG(I,K,J)**2)+                   &
                 1./(LAMR(I,K,J)*LAMG(I,K,J)**3))

! hm 7/15/13, remove limit so that the number of collected drops can smaller than
! number of shed drops
!            NPRACG=MAX(NPRACG-DUM,0.)
            NPRACG=NPRACG-DUM

	    END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(I,K,J).GE.1.E-8 .AND. QC3D(I,K,J).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(I,K,J)*QR3D(I,K,J))
           PRA = 67.*(DUM)**1.15
           NPRA = PRA/(QC3D(I,K,J)/NC3D(I,K,J))

         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(I,K,J).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(I,K,J).lt.dum1) then
            dum=1.
            else if (1./lamr(I,K,J).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(I,K,J)-dum1))
            end if
!            NRAGG = -8.*NR3D(I,K,J)*QR3D(I,K,J)*RHO(I,K,J)
            NRAGG = -5.78*dum*NR3D(I,K,J)*QR3D(I,K,J)*RHO(I,K,J)
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

      IF (QR3D(I,K,J).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(I,K,J)*RHO(I,K,J)*DV*                           &
                   (F1R/(LAMR(I,K,J)*LAMR(I,K,J))+                       &
                    F2R*(ARN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS9/                   &
                (LAMR(I,K,J)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

           IF (QV3D(I,K,J).LT.QVS) THEN
              PRE = EPSR*(QV3D(I,K,J)-QVS)/AB
              PRE = MIN(PRE,0.)
           ELSE
              PRE = 0.
           END IF

!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          IF (QNI3D(I,K,J).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(I,K,J)*T3D(I,K,J)*PRACS
             DUM = -CPW/XLF(I,K,J)*(T3D(I,K,J)-273.15)*PRACS

! hm fix 1/20/15
!             PSMLT=2.*PI*N0S(I,K,J)*KAP*(273.15-T3D(I,K,J))/       &
!                    XLF(I,K,J)*RHO(I,K,J)*(F1S/(LAMS(I,K,J)*LAMS(I,K,J))+        &
!                    F2S*(ASN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
!                    SC**(1./3.)*CONS10/                   &
!                   (LAMS(I,K,J)**CONS35))+DUM
             PSMLT=2.*PI*N0S(I,K,J)*KAP*(273.15-T3D(I,K,J))/       &
                    XLF(I,K,J)*(F1S/(LAMS(I,K,J)*LAMS(I,K,J))+        &
                    F2S*(ASN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS10/                   &
                   (LAMS(I,K,J)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

      IF (QVQVS.LT.1.) THEN
        EPSS = 2.*PI*N0S(I,K,J)*RHO(I,K,J)*DV*                            &
                   (F1S/(LAMS(I,K,J)*LAMS(I,K,J))+                       &
                    F2S*(ASN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS10/                   &
               (LAMS(I,K,J)**CONS35))
! hm fix 8/4/08
        EVPMS = (QV3D(I,K,J)-QVS)*EPSS/AB
        EVPMS = MAX(EVPMS,PSMLT)
        PSMLT = PSMLT-EVPMS
      END IF
      END IF

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          IF (QG3D(I,K,J).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(I,K,J)*T3D(I,K,J)*PRACG
             DUM = -CPW/XLF(I,K,J)*(T3D(I,K,J)-273.15)*PRACG

! hm fix 1/20/15
!             PGMLT=2.*PI*N0G(I,K,J)*KAP*(273.15-T3D(I,K,J))/ 		 &
!                    XLF(I,K,J)*RHO(I,K,J)*(F1S/(LAMG(I,K,J)*LAMG(I,K,J))+                &
!                    F2S*(AGN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
!                    SC**(1./3.)*CONS11/                   &
!                   (LAMG(I,K,J)**CONS36))+DUM
             PGMLT=2.*PI*N0G(I,K,J)*KAP*(273.15-T3D(I,K,J))/ 		 &
                    XLF(I,K,J)*(F1S/(LAMG(I,K,J)*LAMG(I,K,J))+                &
                    F2S*(AGN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS11/                   &
                   (LAMG(I,K,J)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

      IF (QVQVS.LT.1.) THEN
        EPSG = 2.*PI*N0G(I,K,J)*RHO(I,K,J)*DV*                                &
                   (F1S/(LAMG(I,K,J)*LAMG(I,K,J))+                               &
                    F2S*(AGN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS11/                   &
               (LAMG(I,K,J)**CONS36))
! hm fix 8/4/08
        EVPMG = (QV3D(I,K,J)-QVS)*EPSG/AB
        EVPMG = MAX(EVPMG,PGMLT)
        PGMLT = PGMLT-EVPMG
      END IF
      END IF

! HM, V3.2
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

      PRACG = 0.
      PRACS = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

      DUM = (PRC+PRA)*DT

      IF (DUM.GT.QC3D(I,K,J).AND.QC3D(I,K,J).GE.QSMALL) THEN

        RATIO = QC3D(I,K,J)/DUM

        PRC = PRC*RATIO
        PRA = PRA*RATIO

        END IF

! CONSERVATION OF SNOW

        DUM = (-PSMLT-EVPMS+PRACS)*DT

        IF (DUM.GT.QNI3D(I,K,J).AND.QNI3D(I,K,J).GE.QSMALL) THEN

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
        RATIO = QNI3D(I,K,J)/DUM

        PSMLT = PSMLT*RATIO
        EVPMS = EVPMS*RATIO
        PRACS = PRACS*RATIO

        END IF

! CONSERVATION OF GRAUPEL

        DUM = (-PGMLT-EVPMG+PRACG)*DT

        IF (DUM.GT.QG3D(I,K,J).AND.QG3D(I,K,J).GE.QSMALL) THEN

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
        RATIO = QG3D(I,K,J)/DUM

        PGMLT = PGMLT*RATIO
        EVPMG = EVPMG*RATIO
        PRACG = PRACG*RATIO

        END IF

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

        DUM = (-PRACS-PRACG-PRE-PRA-PRC+PSMLT+PGMLT)*DT

        IF (DUM.GT.QR3D(I,K,J).AND.QR3D(I,K,J).GE.QSMALL) THEN

        RATIO = (QR3D(I,K,J)/DT+PRACS+PRACG+PRA+PRC-PSMLT-PGMLT)/ &
                        (-PRE)
        PRE = PRE*RATIO

        END IF

!....................................

      QV3DTEN(I,K,J) = QV3DTEN(I,K,J)+(-PRE-EVPMS-EVPMG)

      T3DTEN(I,K,J) = T3DTEN(I,K,J)+(PRE*XXLV(I,K,J)+(EVPMS+EVPMG)*XXLS(I,K,J)+&
                    (PSMLT+PGMLT-PRACS-PRACG)*XLF(I,K,J))/CPM(I,K,J)

      QC3DTEN(I,K,J) = QC3DTEN(I,K,J)+(-PRA-PRC)
      QR3DTEN(I,K,J) = QR3DTEN(I,K,J)+(PRE+PRA+PRC-PSMLT-PGMLT+PRACS+PRACG)
      QNI3DTEN(I,K,J) = QNI3DTEN(I,K,J)+(PSMLT+EVPMS-PRACS)
      QG3DTEN(I,K,J) = QG3DTEN(I,K,J)+(PGMLT+EVPMG-PRACG)
! fix 053011
!      NS3DTEN(I,K,J) = NS3DTEN(I,K,J)-NPRACS
! HM, bug fix 5/12/08, npracg is subtracted from nr not ng
!      NG3DTEN(I,K,J) = NG3DTEN(I,K,J)
      NC3DTEN(I,K,J) = NC3DTEN(I,K,J)+ (-NPRA-NPRC)
      NR3DTEN(I,K,J) = NR3DTEN(I,K,J)+ (NPRC1+NRAGG-NPRACG)

! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

	C2PREC(I,K,J) = PRA+PRC
      IF (PRE.LT.0.) THEN
         DUM = PRE*DT/QR3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NSUBR = DUM*NR3D(I,K,J)/DT
      END IF

        IF (EVPMS+PSMLT.LT.0.) THEN
         DUM = (EVPMS+PSMLT)*DT/QNI3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NSMLTS = DUM*NS3D(I,K,J)/DT
        END IF
        IF (PSMLT.LT.0.) THEN
          DUM = PSMLT*DT/QNI3D(I,K,J)
          DUM = MAX(-1.0,DUM)
          NSMLTR = DUM*NS3D(I,K,J)/DT
        END IF
        IF (EVPMG+PGMLT.LT.0.) THEN
         DUM = (EVPMG+PGMLT)*DT/QG3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NGMLTG = DUM*NG3D(I,K,J)/DT
        END IF
        IF (PGMLT.LT.0.) THEN
          DUM = PGMLT*DT/QG3D(I,K,J)
          DUM = MAX(-1.0,DUM)
          NGMLTR = DUM*NG3D(I,K,J)/DT
        END IF

         NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+(NSMLTS)
         NG3DTEN(I,K,J) = NG3DTEN(I,K,J)+(NGMLTG)
         NR3DTEN(I,K,J) = NR3DTEN(I,K,J)+(NSUBR-NSMLTR-NGMLTR)

      ENDIF !(QC3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.1.E-8.AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.1.E-8) old CONTINUE 300  

!  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(I,K,J)+DT*T3DTEN(I,K,J)
      DUMQV = QV3D(I,K,J)+DT*QV3DTEN(I,K,J)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J)+DT*T3DTEN(I,K,J),0))
      DUMQSS = EP_2*dum/(PRES(I,K,J)-dum)
      DUMQC = QC3D(I,K,J)+DT*QC3DTEN(I,K,J)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC = DUMS/(1.+XXLV(I,K,J)**2*DUMQSS/(CPM(I,K,J)*RV*DUMT**2))/DT
      IF (PCC*DT+DUMQC.LT.0.) THEN
           PCC = -DUMQC/DT
      END IF

      QV3DTEN(I,K,J) = QV3DTEN(I,K,J)-PCC
      T3DTEN(I,K,J) = T3DTEN(I,K,J)+PCC*XXLV(I,K,J)/CPM(I,K,J)
      QC3DTEN(I,K,J) = QC3DTEN(I,K,J)+PCC

! #if (WRF_CHEM == 1)
!          evapprod(I,K,J) = - PRE - EVPMS - EVPMG
!          rainprod(I,K,J) = PRA + PRC + tqimelt(I,K,J)
! #endif

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC.LT.0.) THEN
!        DUM = PCC*DT/QC3D(I,K,J)
!           DUM = MAX(-1.,DUM)
!        NSUBC = DUM*NC3D(I,K,J)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(I,K,J) = NC3DTEN(I,K,J)+NSUBC

!.....................................................................
!.....................................................................
         ELSE  ! TEMPERATURE < 273.15

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(I,K,J)=NDCNST*1.E6/RHO(I,K,J)
         END IF

! CALCULATE SIZE DISTRIBUTION PARAMETERS
! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(I,K,J) = MAX(0.,NI3D(I,K,J))
      NS3D(I,K,J) = MAX(0.,NS3D(I,K,J))
      NC3D(I,K,J) = MAX(0.,NC3D(I,K,J))
      NR3D(I,K,J) = MAX(0.,NR3D(I,K,J))
      NG3D(I,K,J) = MAX(0.,NG3D(I,K,J))

!......................................................................
! CLOUD ICE

      IF (QI3D(I,K,J).GE.QSMALL) THEN
         LAMI(I,K,J) = (CONS12*                 &
              NI3D(I,K,J)/QI3D(I,K,J))**(1./DI)
         N0I(I,K,J) = NI3D(I,K,J)*LAMI(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(I,K,J).LT.LAMMINI) THEN

      LAMI(I,K,J) = LAMMINI

      N0I(I,K,J) = LAMI(I,K,J)**4*QI3D(I,K,J)/CONS12

      NI3D(I,K,J) = N0I(I,K,J)/LAMI(I,K,J)
      ELSE IF (LAMI(I,K,J).GT.LAMMAXI) THEN
      LAMI(I,K,J) = LAMMAXI
      N0I(I,K,J) = LAMI(I,K,J)**4*QI3D(I,K,J)/CONS12

      NI3D(I,K,J) = N0I(I,K,J)/LAMI(I,K,J)
      END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(I,K,J).GE.QSMALL) THEN
      LAMR(I,K,J) = (PI*RHOW*NR3D(I,K,J)/QR3D(I,K,J))**(1./3.)
      N0RR(I,K,J) = NR3D(I,K,J)*LAMR(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(I,K,J).LT.LAMMINR) THEN

      LAMR(I,K,J) = LAMMINR

      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      ELSE IF (LAMR(I,K,J).GT.LAMMAXR) THEN
      LAMR(I,K,J) = LAMMAXR
      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(I,K,J).GE.QSMALL) THEN

         DUM = PRES(I,K,J)/(287.15*T3D(I,K,J))
         PGAM(I,K,J)=0.0005714*(NC3D(I,K,J)/1.E6*DUM)+0.2714
         PGAM(I,K,J)=1./(PGAM(I,K,J)**2)-1.
         PGAM(I,K,J)=MAX(PGAM(I,K,J),2.)
         PGAM(I,K,J)=MIN(PGAM(I,K,J),10.)

! CALCULATE LAMC

      LAMC(I,K,J) = (CONS26*NC3D(I,K,J)*GAMMA(PGAM(I,K,J)+4.)/   &
                 (QC3D(I,K,J)*GAMMA(PGAM(I,K,J)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(I,K,J)+1.)/60.E-6
      LAMMAX = (PGAM(I,K,J)+1.)/1.E-6

      IF (LAMC(I,K,J).LT.LAMMIN) THEN
      LAMC(I,K,J) = LAMMIN

      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26
      ELSE IF (LAMC(I,K,J).GT.LAMMAX) THEN
      LAMC(I,K,J) = LAMMAX
      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26

      END IF

! TO CALCULATE DROPLET FREEZING

        CDIST1(I,K,J) = NC3D(I,K,J)/GAMMA(PGAM(I,K,J)+1.)

      END IF

!......................................................................
! SNOW

      IF (QNI3D(I,K,J).GE.QSMALL) THEN
      LAMS(I,K,J) = (CONS1*NS3D(I,K,J)/QNI3D(I,K,J))**(1./DS)
      N0S(I,K,J) = NS3D(I,K,J)*LAMS(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(I,K,J).LT.LAMMINS) THEN
      LAMS(I,K,J) = LAMMINS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1

      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)

      ELSE IF (LAMS(I,K,J).GT.LAMMAXS) THEN

      LAMS(I,K,J) = LAMMAXS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1

      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(I,K,J).GE.QSMALL) THEN
      LAMG(I,K,J) = (CONS2*NG3D(I,K,J)/QG3D(I,K,J))**(1./DG)
      N0G(I,K,J) = NG3D(I,K,J)*LAMG(I,K,J)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(I,K,J).LT.LAMMING) THEN
      LAMG(I,K,J) = LAMMING
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)

      ELSE IF (LAMG(I,K,J).GT.LAMMAXG) THEN

      LAMG(I,K,J) = LAMMAXG
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)
      END IF
      END IF

!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC = 0.
            NNUCCC = 0.
            PRC = 0.
            NPRC = 0.
            NPRC1 = 0.
            NSAGG = 0.
            PSACWS = 0.
            NPSACWS = 0.
            PSACWI = 0.
            NPSACWI = 0.
            PRACS = 0.
            NPRACS = 0.
            NMULTS = 0.
            QMULTS = 0.
            NMULTR = 0.
            QMULTR = 0.
            NMULTG = 0.
            QMULTG = 0.
            NMULTRG = 0.
            QMULTRG = 0.
            MNUCCR = 0.
            NNUCCR = 0.
            PRA = 0.
            NPRA = 0.
            NRAGG = 0.
            PRCI = 0.
            NPRCI = 0.
            PRAI = 0.
            NPRAI = 0.
            NNUCCD = 0.
            MNUCCD = 0.
            PCC = 0.
            PRE = 0.
            PRD = 0.
            PRDS = 0.
            EPRD = 0.
            EPRDS = 0.
            NSUBC = 0.
            NSUBI = 0.
            NSUBS = 0.
            NSUBR = 0.
            PIACR = 0.
            NIACR = 0.
            PRACI = 0.
            PIACRS = 0.
            NIACRS = 0.
            PRACIS = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG = 0.
            PSACR = 0.
	    PSACWG = 0.
	    PGSACW = 0.
            PGRACS = 0.
	    PRDG = 0.
	    EPRDG = 0.
	    NPRACG = 0.
	    NPSACWG = 0.
	    NSCNG = 0.
 	    NGRACS = 0.
	    NSUBG = 0.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
        IF (QC3D(I,K,J).GE.QSMALL .AND. T3D(I,K,J).LT.269.15) THEN

! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

           NACNT = EXP(-2.80+0.262*(273.15-T3D(I,K,J)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(273.15-T3D(I,K,J)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(273.15-T3D(I,K,J)))

! CONTACT FREEZING

! MEAN FREE PATH

            DUM = 7.37*T3D(I,K,J)/(288.*10.*PRES(I,K,J))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

            DAP(I,K,J) = CONS37*T3D(I,K,J)*(1.+DUM/RIN)/MU(I,K,J)

           MNUCCC = CONS38*DAP(I,K,J)*NACNT*EXP(LOG(CDIST1(I,K,J))+   &
                   LOG(GAMMA(PGAM(I,K,J)+5.))-4.*LOG(LAMC(I,K,J)))
           NNUCCC = 2.*PI*DAP(I,K,J)*NACNT*CDIST1(I,K,J)*           &
                    GAMMA(PGAM(I,K,J)+2.)/                         &
                    LAMC(I,K,J)

! IMMERSION FREEZING (BIGG 1953)

!           MNUCCC = MNUCCC+CONS39*                   &
!                  EXP(LOG(CDIST1(I,K,J))+LOG(GAMMA(7.+PGAM(I,K,J)))-6.*LOG(LAMC(I,K,J)))*             &
!                   EXP(AIMM*(273.15-T3D(I,K,J)))

!           NNUCCC = NNUCCC+                                  &
!            CONS40*EXP(LOG(CDIST1(I,K,J))+LOG(GAMMA(PGAM(I,K,J)+4.))-3.*LOG(LAMC(I,K,J)))              &
!                *EXP(AIMM*(273.15-T3D(I,K,J)))

! hm 7/15/13 fix for consistency w/ original formula
           MNUCCC = MNUCCC+CONS39*                   &
                  EXP(LOG(CDIST1(I,K,J))+LOG(GAMMA(7.+PGAM(I,K,J)))-6.*LOG(LAMC(I,K,J)))*             &
                   (EXP(AIMM*(273.15-T3D(I,K,J)))-1.)

           NNUCCC = NNUCCC+                                  &
            CONS40*EXP(LOG(CDIST1(I,K,J))+LOG(GAMMA(PGAM(I,K,J)+4.))-3.*LOG(LAMC(I,K,J)))              &
                *(EXP(AIMM*(273.15-T3D(I,K,J)))-1.)

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

           NNUCCC = MIN(NNUCCC,NC3D(I,K,J)/DT)

        END IF

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(I,K,J).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC=1350.*QC3D(I,K,J)**2.47*  &
           (NC3D(I,K,J)/1.e6*RHO(I,K,J))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1 = PRC/CONS29
        NPRC = PRC/(QC3D(I,K,J)/NC3D(I,K,J))

! hm bug fix 3/20/12
                NPRC = MIN(NPRC,NC3D(I,K,J)/DT)
                NPRC1 = MIN(NPRC1,NPRC)

         END IF

!.......................................................................
! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

         IF (QNI3D(I,K,J).GE.1.E-8) THEN
             NSAGG = CONS15*ASN(I,K,J)*RHO(I,K,J)**            &
            ((2.+BS)/3.)*QNI3D(I,K,J)**((2.+BS)/3.)*                  &
            (NS3D(I,K,J)*RHO(I,K,J))**((4.-BS)/3.)/                       &
            (RHO(I,K,J))
         END IF

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

         IF (QNI3D(I,K,J).GE.1.E-8 .AND. QC3D(I,K,J).GE.QSMALL) THEN

           PSACWS = CONS13*ASN(I,K,J)*QC3D(I,K,J)*RHO(I,K,J)*               &
                  N0S(I,K,J)/                        &
                  LAMS(I,K,J)**(BS+3.)
           NPSACWS = CONS13*ASN(I,K,J)*NC3D(I,K,J)*RHO(I,K,J)*              &
                  N0S(I,K,J)/                        &
                  LAMS(I,K,J)**(BS+3.)

         END IF

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

         IF (QG3D(I,K,J).GE.1.E-8 .AND. QC3D(I,K,J).GE.QSMALL) THEN

           PSACWG = CONS14*AGN(I,K,J)*QC3D(I,K,J)*RHO(I,K,J)*               &
                  N0G(I,K,J)/                        &
                  LAMG(I,K,J)**(BG+3.)
           NPSACWG = CONS14*AGN(I,K,J)*NC3D(I,K,J)*RHO(I,K,J)*              &
                  N0G(I,K,J)/                        &
                  LAMG(I,K,J)**(BG+3.)
	    END IF

!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

         IF (QI3D(I,K,J).GE.1.E-8 .AND. QC3D(I,K,J).GE.QSMALL) THEN

! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

            IF (1./LAMI(I,K,J).GE.100.E-6) THEN

           PSACWI = CONS16*AIN(I,K,J)*QC3D(I,K,J)*RHO(I,K,J)*               &
                  N0I(I,K,J)/                        &
                  LAMI(I,K,J)**(BI+3.)
           NPSACWI = CONS16*AIN(I,K,J)*NC3D(I,K,J)*RHO(I,K,J)*              &
                  N0I(I,K,J)/                        &
                  LAMI(I,K,J)**(BI+3.)
           END IF
         END IF

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         IF (QR3D(I,K,J).GE.1.E-8.AND.QNI3D(I,K,J).GE.1.E-8) THEN

            UMS = ASN(I,K,J)*CONS3/(LAMS(I,K,J)**BS)
            UMR = ARN(I,K,J)*CONS4/(LAMR(I,K,J)**BR)
            UNS = ASN(I,K,J)*CONS5/LAMS(I,K,J)**BS
            UNR = ARN(I,K,J)*CONS6/LAMR(I,K,J)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/RHO(I,K,J))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(I,K,J)*                      &
                  N0RR(I,K,J)*N0S(I,K,J)/LAMR(I,K,J)**3*                              &
                  (5./(LAMR(I,K,J)**3*LAMS(I,K,J))+                    &
                  2./(LAMR(I,K,J)**2*LAMS(I,K,J)**2)+                  &
                  0.5/(LAMR(I,K,J)*LAMS(I,K,J)**3)))

            NPRACS = CONS32*RHO(I,K,J)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(I,K,J)*N0S(I,K,J)*              &
                (1./(LAMR(I,K,J)**3*LAMS(I,K,J))+                      &
                 1./(LAMR(I,K,J)**2*LAMS(I,K,J)**2)+                   &
                 1./(LAMR(I,K,J)*LAMS(I,K,J)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS = MIN(PRACS,QR3D(I,K,J)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! HM MODIFY FOR WRFV3.1
!            IF (IHAIL.EQ.0) THEN
            IF (QNI3D(I,K,J).GE.0.1E-3.AND.QR3D(I,K,J).GE.0.1E-3) THEN
            PSACR = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(I,K,J)*                     &
                 N0RR(I,K,J)*N0S(I,K,J)/LAMS(I,K,J)**3*                               &
                  (5./(LAMS(I,K,J)**3*LAMR(I,K,J))+                    &
                  2./(LAMS(I,K,J)**2*LAMR(I,K,J)**2)+                  &
                  0.5/(LAMS(I,K,J)*LAMR(I,K,J)**3)))
            END IF
!            END IF

         END IF

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990,
! USED BY REISNER ET AL 1998
         IF (QR3D(I,K,J).GE.1.E-8.AND.QG3D(I,K,J).GE.1.E-8) THEN

            UMG = AGN(I,K,J)*CONS7/(LAMG(I,K,J)**BG)
            UMR = ARN(I,K,J)*CONS4/(LAMR(I,K,J)**BR)
            UNG = AGN(I,K,J)*CONS8/LAMG(I,K,J)**BG
            UNR = ARN(I,K,J)*CONS6/LAMR(I,K,J)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/RHO(I,K,J))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(I,K,J)*                      &
                  N0RR(I,K,J)*N0G(I,K,J)/LAMR(I,K,J)**3*                              &
                  (5./(LAMR(I,K,J)**3*LAMG(I,K,J))+                    &
                  2./(LAMR(I,K,J)**2*LAMG(I,K,J)**2)+				   &
				  0.5/(LAMR(I,K,J)*LAMG(I,K,J)**3)))

            NPRACG = CONS32*RHO(I,K,J)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(I,K,J)*N0G(I,K,J)*              &
                (1./(LAMR(I,K,J)**3*LAMG(I,K,J))+                      &
                 1./(LAMR(I,K,J)**2*LAMG(I,K,J)**2)+                   &
                 1./(LAMR(I,K,J)*LAMG(I,K,J)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG = MIN(PRACG,QR3D(I,K,J)/DT)

	    END IF

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         IF (QNI3D(I,K,J).GE.0.1E-3) THEN
         IF (QC3D(I,K,J).GE.0.5E-3.OR.QR3D(I,K,J).GE.0.1E-3) THEN
         IF (PSACWS.GT.0..OR.PRACS.GT.0.) THEN
            IF (T3D(I,K,J).LT.270.16 .AND. T3D(I,K,J).GT.265.16) THEN

               IF (T3D(I,K,J).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(I,K,J).LE.270.16.AND.T3D(I,K,J).GT.268.16)  THEN
                  FMULT = (270.16-T3D(I,K,J))/2.
               ELSE IF (T3D(I,K,J).GE.265.16.AND.T3D(I,K,J).LE.268.16)   THEN
                  FMULT = (T3D(I,K,J)-265.16)/3.
               ELSE IF (T3D(I,K,J).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               IF (PSACWS.GT.0.) THEN
                  NMULTS = 35.E4*PSACWS*FMULT*1000.
                  QMULTS = NMULTS*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS = MIN(QMULTS,PSACWS)
                  PSACWS = PSACWS-QMULTS

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACS.GT.0.) THEN
                   NMULTR = 35.E4*PRACS*FMULT*1000.
                   QMULTR = NMULTR*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR = MIN(QMULTR,PRACS)

                   PRACS = PRACS-QMULTR

               END IF

            END IF
         END IF
         END IF
         END IF

!.......................................................................
! RIME-SPLINTERING - GRAUPEL
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

!         IF (IHAIL.EQ.0) THEN
! v1.4
         IF (QG3D(I,K,J).GE.0.1E-3) THEN
         IF (QC3D(I,K,J).GE.0.5E-3.OR.QR3D(I,K,J).GE.0.1E-3) THEN
         IF (PSACWG.GT.0..OR.PRACG.GT.0.) THEN
            IF (T3D(I,K,J).LT.270.16 .AND. T3D(I,K,J).GT.265.16) THEN

               IF (T3D(I,K,J).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(I,K,J).LE.270.16.AND.T3D(I,K,J).GT.268.16)  THEN
                  FMULT = (270.16-T3D(I,K,J))/2.
               ELSE IF (T3D(I,K,J).GE.265.16.AND.T3D(I,K,J).LE.268.16)   THEN
                  FMULT = (T3D(I,K,J)-265.16)/3.
               ELSE IF (T3D(I,K,J).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               IF (PSACWG.GT.0.) THEN
                  NMULTG = 35.E4*PSACWG*FMULT*1000.
                  QMULTG = NMULTG*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG = MIN(QMULTG,PSACWG)
                  PSACWG = PSACWG-QMULTG

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACG.GT.0.) THEN
                   NMULTRG = 35.E4*PRACG*FMULT*1000.
                   QMULTRG = NMULTRG*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG = MIN(QMULTRG,PRACG)
                   PRACG = PRACG-QMULTRG

               END IF
               END IF
               END IF
            END IF
            END IF
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

!           IF (IHAIL.EQ.0) THEN
	   IF (PSACWS.GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(I,K,J).GE.0.1E-3.AND.QC3D(I,K,J).GE.0.5E-3) THEN

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
	     PGSACW = MIN(PSACWS,CONS17*DT*N0S(I,K,J)*QC3D(I,K,J)*QC3D(I,K,J)* &
                          ASN(I,K,J)*ASN(I,K,J)/ &
                           (RHO(I,K,J)*LAMS(I,K,J)**(2.*BS+2.)))

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
	     DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW,0.)

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
	     NSCNG = DUM/MG0*RHO(I,K,J)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG = MIN(NSCNG,NS3D(I,K,J)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS = PSACWS - PGSACW
             END IF
	   END IF

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

	   IF (PRACS.GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(I,K,J).GE.0.1E-3.AND.QR3D(I,K,J).GE.0.1E-3) THEN
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
	      DUM = CONS18*(4./LAMS(I,K,J))**3*(4./LAMS(I,K,J))**3 &
                   /(CONS18*(4./LAMS(I,K,J))**3*(4./LAMS(I,K,J))**3+ &
                   CONS19*(4./LAMR(I,K,J))**3*(4./LAMR(I,K,J))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
	      PGRACS = (1.-DUM)*PRACS
            NGRACS = (1.-DUM)*NPRACS
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS = MIN(NGRACS,NR3D(I,K,J)/DT)
            NGRACS = MIN(NGRACS,NS3D(I,K,J)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS = PRACS - PGRACS
            NPRACS = NPRACS - NGRACS
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR=PSACR*(1.-DUM)
            END IF
	   END IF
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         IF (T3D(I,K,J).LT.269.15.AND.QR3D(I,K,J).GE.QSMALL) THEN

! IMMERSION FREEZING (BIGG 1953)
!            MNUCCR = CONS20*NR3D(I,K,J)*EXP(AIMM*(273.15-T3D(I,K,J)))/LAMR(I,K,J)**3 &
!                 /LAMR(I,K,J)**3

!            NNUCCR = PI*NR3D(I,K,J)*BIMM*EXP(AIMM*(273.15-T3D(I,K,J)))/LAMR(I,K,J)**3

! hm fix 7/15/13 for consistency w/ original formula
            MNUCCR = CONS20*NR3D(I,K,J)*(EXP(AIMM*(273.15-T3D(I,K,J)))-1.)/LAMR(I,K,J)**3 &
                 /LAMR(I,K,J)**3

            NNUCCR = PI*NR3D(I,K,J)*BIMM*(EXP(AIMM*(273.15-T3D(I,K,J)))-1.)/LAMR(I,K,J)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR = MIN(NNUCCR,NR3D(I,K,J)/DT)

         END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(I,K,J).GE.1.E-8 .AND. QC3D(I,K,J).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(I,K,J)*QR3D(I,K,J))
           PRA = 67.*(DUM)**1.15
           NPRA = PRA/(QC3D(I,K,J)/NC3D(I,K,J))

         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(I,K,J).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(I,K,J).lt.dum1) then
            dum=1.
            else if (1./lamr(I,K,J).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(I,K,J)-dum1))
            end if
!            NRAGG = -8.*NR3D(I,K,J)*QR3D(I,K,J)*RHO(I,K,J)
            NRAGG = -5.78*dum*NR3D(I,K,J)*QR3D(I,K,J)*RHO(I,K,J)
         END IF

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

         IF (QI3D(I,K,J).GE.1.E-8 .AND.QVQVSI.GE.1.) THEN

!           COFFI = 2./LAMI(I,K,J)
!           IF (COFFI.GE.DCS) THEN
              NPRCI = CONS21*(QV3D(I,K,J)-QVI)*RHO(I,K,J)                         &
                *N0I(I,K,J)*EXP(-LAMI(I,K,J)*DCS)*DV/ABI
              PRCI = CONS22*NPRCI
              NPRCI = MIN(NPRCI,NI3D(I,K,J)/DT)

!           END IF
         END IF

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         IF (QNI3D(I,K,J).GE.1.E-8 .AND. QI3D(I,K,J).GE.QSMALL) THEN
            PRAI = CONS23*ASN(I,K,J)*QI3D(I,K,J)*RHO(I,K,J)*N0S(I,K,J)/     &
                     LAMS(I,K,J)**(BS+3.)
            NPRAI = CONS23*ASN(I,K,J)*NI3D(I,K,J)*                                       &
                  RHO(I,K,J)*N0S(I,K,J)/                                 &
                  LAMS(I,K,J)**(BS+3.)
            NPRAI=MIN(NPRAI,NI3D(I,K,J)/DT)
         END IF

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

         IF (QR3D(I,K,J).GE.1.E-8.AND.QI3D(I,K,J).GE.1.E-8.AND.T3D(I,K,J).LE.273.15) THEN

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            IF (QR3D(I,K,J).GE.0.1E-3) THEN
            NIACR=CONS24*NI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J) &
                /LAMR(I,K,J)**(BR+3.)*RHO(I,K,J)
            PIACR=CONS25*NI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J) &
                /LAMR(I,K,J)**(BR+3.)/LAMR(I,K,J)**3*RHO(I,K,J)
            PRACI=CONS24*QI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J)/ &
                LAMR(I,K,J)**(BR+3.)*RHO(I,K,J)
            NIACR=MIN(NIACR,NR3D(I,K,J)/DT)
            NIACR=MIN(NIACR,NI3D(I,K,J)/DT)
            ELSE
            NIACRS=CONS24*NI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J) &
                /LAMR(I,K,J)**(BR+3.)*RHO(I,K,J)
            PIACRS=CONS25*NI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J) &
                /LAMR(I,K,J)**(BR+3.)/LAMR(I,K,J)**3*RHO(I,K,J)
            PRACIS=CONS24*QI3D(I,K,J)*N0RR(I,K,J)*ARN(I,K,J)/ &
                LAMR(I,K,J)**(BR+3.)*RHO(I,K,J)
            NIACRS=MIN(NIACRS,NR3D(I,K,J)/DT)
            NIACRS=MIN(NIACRS,NI3D(I,K,J)/DT)
            END IF
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

         IF (INUC.EQ.0) THEN

! add threshold according to Greg Thomspon

         if ((QVQVS.GE.0.999.and.T3D(I,K,J).le.265.15).or. &
              QVQVSI.ge.1.08) then

! hm, modify dec. 5, 2006, replace with cooper curve
      kc2 = 0.005*exp(0.304*(273.15-T3D(I,K,J)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
      kc2 = min(kc2,500.e3)
      kc2=MAX(kc2/RHO(I,K,J),0.)  ! convert to kg-1

          IF (KC2.GT.NI3D(I,K,J)+NS3D(I,K,J)+NG3D(I,K,J)) THEN
             NNUCCD = (KC2-NI3D(I,K,J)-NS3D(I,K,J)-NG3D(I,K,J))/DT
             MNUCCD = NNUCCD*MI0
          END IF

          END IF

          ELSE IF (INUC.EQ.1) THEN

          IF (T3D(I,K,J).LT.273.15.AND.QVQVSI.GT.1.) THEN

             KC2 = 0.16*1000./RHO(I,K,J)  ! CONVERT FROM L-1 TO KG-1
          IF (KC2.GT.NI3D(I,K,J)+NS3D(I,K,J)+NG3D(I,K,J)) THEN
             NNUCCD = (KC2-NI3D(I,K,J)-NS3D(I,K,J)-NG3D(I,K,J))/DT
             MNUCCD = NNUCCD*MI0
          END IF
          END IF

         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!  101      CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

! NO VENTILATION FOR CLOUD ICE

        IF (QI3D(I,K,J).GE.QSMALL) THEN

         EPSI = 2.*PI*N0I(I,K,J)*RHO(I,K,J)*DV/(LAMI(I,K,J)*LAMI(I,K,J))

      ELSE
         EPSI = 0.
      END IF

      IF (QNI3D(I,K,J).GE.QSMALL) THEN
        EPSS = 2.*PI*N0S(I,K,J)*RHO(I,K,J)*DV*                            &
                   (F1S/(LAMS(I,K,J)*LAMS(I,K,J))+                       &
                    F2S*(ASN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS10/                   &
               (LAMS(I,K,J)**CONS35))
      ELSE
      EPSS = 0.
      END IF

      IF (QG3D(I,K,J).GE.QSMALL) THEN
        EPSG = 2.*PI*N0G(I,K,J)*RHO(I,K,J)*DV*                                &
                   (F1S/(LAMG(I,K,J)*LAMG(I,K,J))+                               &
                    F2S*(AGN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS11/                   &
               (LAMG(I,K,J)**CONS36))


      ELSE
      EPSG = 0.
      END IF

      IF (QR3D(I,K,J).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(I,K,J)*RHO(I,K,J)*DV*                           &
                   (F1R/(LAMR(I,K,J)*LAMR(I,K,J))+                       &
                    F2R*(ARN(I,K,J)*RHO(I,K,J)/MU(I,K,J))**0.5*                      &
                    SC**(1./3.)*CONS9/                   &
                (LAMR(I,K,J)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
              IF (QI3D(I,K,J).GE.QSMALL) THEN
              DUM=(1.-EXP(-LAMI(I,K,J)*DCS)*(1.+LAMI(I,K,J)*DCS))
              PRD = EPSI*(QV3D(I,K,J)-QVI)/ABI*DUM
              ELSE
              DUM=0.
              END IF
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              IF (QNI3D(I,K,J).GE.QSMALL) THEN
              PRDS = EPSS*(QV3D(I,K,J)-QVI)/ABI+ &
                EPSI*(QV3D(I,K,J)-QVI)/ABI*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
              ELSE
              PRD = PRD+EPSI*(QV3D(I,K,J)-QVI)/ABI*(1.-DUM)
              END IF
! VAPOR DPEOSITION ON GRAUPEL
              PRDG = EPSG*(QV3D(I,K,J)-QVI)/ABI

! NO CONDENSATION ONTO RAIN, ONLY EVAP

           IF (QV3D(I,K,J).LT.QVS) THEN
              PRE = EPSR*(QV3D(I,K,J)-QVS)/AB
              PRE = MIN(PRE,0.)
           ELSE
              PRE = 0.
           END IF

! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME

           DUM = (QV3D(I,K,J)-QVI)/DT

           FUDGEF = 0.9999
           SUM_DEP = PRD+PRDS+MNUCCD+PRDG

           IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
               MNUCCD = FUDGEF*MNUCCD*DUM/SUM_DEP
               PRD = FUDGEF*PRD*DUM/SUM_DEP
               PRDS = FUDGEF*PRDS*DUM/SUM_DEP
	       PRDG = FUDGEF*PRDG*DUM/SUM_DEP
           ENDIF

! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

           IF (PRD.LT.0.) THEN
              EPRD=PRD
              PRD=0.
           END IF
           IF (PRDS.LT.0.) THEN
              EPRDS=PRDS
              PRDS=0.
           END IF
           IF (PRDG.LT.0.) THEN
              EPRDG=PRDG
              PRDG=0.
           END IF

!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

      IF (ILIQ.EQ.1) THEN
      MNUCCC=0.
      NNUCCC=0.
      MNUCCR=0.
      NNUCCR=0.
      MNUCCD=0.
      NNUCCD=0.
      END IF

! ****SENSITIVITY - NO GRAUPEL
      IF (IGRAUP.EQ.1) THEN
            PRACG = 0.
            PSACR = 0.
	    PSACWG = 0.
	    PRDG = 0.
	    EPRDG = 0.
            EVPMG = 0.
            PGMLT = 0.
	    NPRACG = 0.
	    NPSACWG = 0.
	    NSCNG = 0.
 	    NGRACS = 0.
	    NSUBG = 0.
	    NGMLTG = 0.
            NGMLTR = 0.
! fix 053011
            PIACRS=PIACRS+PIACR
            PIACR = 0.
! fix 070713
	    PRACIS=PRACIS+PRACI
	    PRACI = 0.
	    PSACWS=PSACWS+PGSACW
	    PGSACW = 0.
	    PRACS=PRACS+PGRACS
	    PGRACS = 0.
       END IF

! CONSERVATION OF QC

      DUM = (PRC+PRA+MNUCCC+PSACWS+PSACWI+QMULTS+PSACWG+PGSACW+QMULTG)*DT

      IF (DUM.GT.QC3D(I,K,J).AND.QC3D(I,K,J).GE.QSMALL) THEN
        RATIO = QC3D(I,K,J)/DUM

        PRC = PRC*RATIO
        PRA = PRA*RATIO
        MNUCCC = MNUCCC*RATIO
        PSACWS = PSACWS*RATIO
        PSACWI = PSACWI*RATIO
        QMULTS = QMULTS*RATIO
        QMULTG = QMULTG*RATIO
        PSACWG = PSACWG*RATIO
	PGSACW = PGSACW*RATIO
        END IF

! CONSERVATION OF QI

      DUM = (-PRD-MNUCCC+PRCI+PRAI-QMULTS-QMULTG-QMULTR-QMULTRG &
                -MNUCCD+PRACI+PRACIS-EPRD-PSACWI)*DT

      IF (DUM.GT.QI3D(I,K,J).AND.QI3D(I,K,J).GE.QSMALL) THEN

        RATIO = (QI3D(I,K,J)/DT+PRD+MNUCCC+QMULTS+QMULTG+QMULTR+QMULTRG+ &
                     MNUCCD+PSACWI)/ &
                      (PRCI+PRAI+PRACI+PRACIS-EPRD)

        PRCI = PRCI*RATIO
        PRAI = PRAI*RATIO
        PRACI = PRACI*RATIO
        PRACIS = PRACIS*RATIO
        EPRD = EPRD*RATIO

        END IF

! CONSERVATION OF QR

      DUM=((PRACS-PRE)+(QMULTR+QMULTRG-PRC)+(MNUCCR-PRA)+ &
             PIACR+PIACRS+PGRACS+PRACG)*DT

      IF (DUM.GT.QR3D(I,K,J).AND.QR3D(I,K,J).GE.QSMALL) THEN

        RATIO = (QR3D(I,K,J)/DT+PRC+PRA)/ &
             (-PRE+QMULTR+QMULTRG+PRACS+MNUCCR+PIACR+PIACRS+PGRACS+PRACG)

        PRE = PRE*RATIO
        PRACS = PRACS*RATIO
        QMULTR = QMULTR*RATIO
        QMULTRG = QMULTRG*RATIO
        MNUCCR = MNUCCR*RATIO
        PIACR = PIACR*RATIO
        PIACRS = PIACRS*RATIO
        PGRACS = PGRACS*RATIO
        PRACG = PRACG*RATIO

        END IF

! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

        IF (IGRAUP.EQ.0) THEN

      DUM = (-PRDS-PSACWS-PRAI-PRCI-PRACS-EPRDS+PSACR-PIACRS-PRACIS)*DT

      IF (DUM.GT.QNI3D(I,K,J).AND.QNI3D(I,K,J).GE.QSMALL) THEN

        RATIO = (QNI3D(I,K,J)/DT+PRDS+PSACWS+PRAI+PRCI+PRACS+PIACRS+PRACIS)/(-EPRDS+PSACR)

       EPRDS = EPRDS*RATIO
       PSACR = PSACR*RATIO

       END IF

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN

      DUM = (-PRDS-PSACWS-PRAI-PRCI-PRACS-EPRDS+PSACR-PIACRS-PRACIS-MNUCCR)*DT

      IF (DUM.GT.QNI3D(I,K,J).AND.QNI3D(I,K,J).GE.QSMALL) THEN

       RATIO = (QNI3D(I,K,J)/DT+PRDS+PSACWS+PRAI+PRCI+PRACS+PIACRS+PRACIS+MNUCCR)/(-EPRDS+PSACR)

       EPRDS = EPRDS*RATIO
       PSACR = PSACR*RATIO

       END IF

       END IF

! CONSERVATION OF QG

      DUM = (-PSACWG-PRACG-PGSACW-PGRACS-PRDG-MNUCCR-EPRDG-PIACR-PRACI-PSACR)*DT

      IF (DUM.GT.QG3D(I,K,J).AND.QG3D(I,K,J).GE.QSMALL) THEN

        RATIO = (QG3D(I,K,J)/DT+PSACWG+PRACG+PGSACW+PGRACS+PRDG+MNUCCR+PSACR+&
                  PIACR+PRACI)/(-EPRDG)

       EPRDG = EPRDG*RATIO

      END IF

! TENDENCIES

      QV3DTEN(I,K,J) = QV3DTEN(I,K,J)+(-PRE-PRD-PRDS-MNUCCD-EPRD-EPRDS-PRDG-EPRDG)

! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
      T3DTEN(I,K,J) = T3DTEN(I,K,J)+(PRE                                 &
               *XXLV(I,K,J)+(PRD+PRDS+                            &
                MNUCCD+EPRD+EPRDS+PRDG+EPRDG)*XXLS(I,K,J)+         &
               (PSACWS+PSACWI+MNUCCC+MNUCCR+                      &
                QMULTS+QMULTG+QMULTR+QMULTRG+PRACS &
                +PSACWG+PRACG+PGSACW+PGRACS+PIACR+PIACRS)*XLF(I,K,J))/CPM(I,K,J)

      QC3DTEN(I,K,J) = QC3DTEN(I,K,J)+                                      &
                 (-PRA-PRC-MNUCCC+PCC-                  &
                  PSACWS-PSACWI-QMULTS-QMULTG-PSACWG-PGSACW)
      QI3DTEN(I,K,J) = QI3DTEN(I,K,J)+                                      &
         (PRD+EPRD+PSACWI+MNUCCC-PRCI-                                 &
                  PRAI+QMULTS+QMULTG+QMULTR+QMULTRG+MNUCCD-PRACI-PRACIS)
      QR3DTEN(I,K,J) = QR3DTEN(I,K,J)+                                      &
                 (PRE+PRA+PRC-PRACS-MNUCCR-QMULTR-QMULTRG &
             -PIACR-PIACRS-PRACG-PGRACS)

      IF (IGRAUP.EQ.0) THEN

      QNI3DTEN(I,K,J) = QNI3DTEN(I,K,J)+                                    &
           (PRAI+PSACWS+PRDS+PRACS+PRCI+EPRDS-PSACR+PIACRS+PRACIS)
      NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+(NSAGG+NPRCI-NSCNG-NGRACS+NIACRS)
      QG3DTEN(I,K,J) = QG3DTEN(I,K,J)+(PRACG+PSACWG+PGSACW+PGRACS+ &
                    PRDG+EPRDG+MNUCCR+PIACR+PRACI+PSACR)
      NG3DTEN(I,K,J) = NG3DTEN(I,K,J)+(NSCNG+NGRACS+NNUCCR+NIACR)

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      ELSE IF (IGRAUP.EQ.1) THEN

      QNI3DTEN(I,K,J) = QNI3DTEN(I,K,J)+                                    &
           (PRAI+PSACWS+PRDS+PRACS+PRCI+EPRDS-PSACR+PIACRS+PRACIS+MNUCCR)
      NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+(NSAGG+NPRCI-NSCNG-NGRACS+NIACRS+NNUCCR)

      END IF

      NC3DTEN(I,K,J) = NC3DTEN(I,K,J)+(-NNUCCC-NPSACWS                &
            -NPRA-NPRC-NPSACWI-NPSACWG)

      NI3DTEN(I,K,J) = NI3DTEN(I,K,J)+                                      &
       (NNUCCC-NPRCI-NPRAI+NMULTS+NMULTG+NMULTR+NMULTRG+ &
               NNUCCD-NIACR-NIACRS)

      NR3DTEN(I,K,J) = NR3DTEN(I,K,J)+(NPRC1-NPRACS-NNUCCR      &
                   +NRAGG-NIACR-NIACRS-NPRACG-NGRACS)

! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

	C2PREC(I,K,J) = PRA+PRC+PSACWS+QMULTS+QMULTG+PSACWG+ &
       PGSACW+MNUCCC+PSACWI
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(I,K,J)+DT*T3DTEN(I,K,J)
      DUMQV = QV3D(I,K,J)+DT*QV3DTEN(I,K,J)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J)+DT*T3DTEN(I,K,J),0))
      DUMQSS = EP_2*dum/(PRES(I,K,J)-dum)
      DUMQC = QC3D(I,K,J)+DT*QC3DTEN(I,K,J)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC = DUMS/(1.+XXLV(I,K,J)**2*DUMQSS/(CPM(I,K,J)*RV*DUMT**2))/DT
      IF (PCC*DT+DUMQC.LT.0.) THEN
           PCC = -DUMQC/DT
      END IF

      QV3DTEN(I,K,J) = QV3DTEN(I,K,J)-PCC
      T3DTEN(I,K,J) = T3DTEN(I,K,J)+PCC*XXLV(I,K,J)/CPM(I,K,J)
      QC3DTEN(I,K,J) = QC3DTEN(I,K,J)+PCC

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC.LT.0.) THEN
!        DUM = PCC*DT/QC3D(I,K,J)
!           DUM = MAX(-1.,DUM)
!        NSUBC = DUM*NC3D(I,K,J)/DT
!     END IF

      IF (EPRD.LT.0.) THEN
         DUM = EPRD*DT/QI3D(I,K,J)
            DUM = MAX(-1.,DUM)
         NSUBI = DUM*NI3D(I,K,J)/DT
      END IF
      IF (EPRDS.LT.0.) THEN
         DUM = EPRDS*DT/QNI3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NSUBS = DUM*NS3D(I,K,J)/DT
      END IF
      IF (PRE.LT.0.) THEN
         DUM = PRE*DT/QR3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NSUBR = DUM*NR3D(I,K,J)/DT
      END IF
      IF (EPRDG.LT.0.) THEN
         DUM = EPRDG*DT/QG3D(I,K,J)
           DUM = MAX(-1.,DUM)
         NSUBG = DUM*NG3D(I,K,J)/DT
      END IF

!        nsubr=0.
!        nsubs=0.
!        nsubg=0.

! UPDATE TENDENCIES

!        NC3DTEN(I,K,J) = NC3DTEN(I,K,J)+NSUBC
         NI3DTEN(I,K,J) = NI3DTEN(I,K,J)+NSUBI
         NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+NSUBS
         NG3DTEN(I,K,J) = NG3DTEN(I,K,J)+NSUBG
         NR3DTEN(I,K,J) = NR3DTEN(I,K,J)+NSUBR

! #if (WRF_CHEM == 1)
!          evapprod(I,K,J) = - PRE - EPRDS - EPRDG
!          rainprod(I,K,J) = PRA + PRC + PSACWS + PSACWG + PGSACW &
!                        + PRAI + PRCI + PRACI + PRACIS + &
!                        + PRDS + PRDG
! #endif

         END IF !!!!!! TEMPERATURE

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE(I,K,J) = 1

      ENDIF !(QC3D(I,K,J).LT.QSMALL.AND.QI3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.QSMALL &
            !      .AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.QSMALL.AND.( &
            !      (T3D(I,K,J).LT.273.15.AND.QVQVSI.LT.0.999).OR. &
            !      (T3D(I,K,J).GE.273.15.AND.QVQVS.LT.0.999))) -- old CONTINUE 200
        END DO
      END DO
      END DO

      !calculate LTRUE_COL(I,K,J), where LTRUE_COL=1 if LTRUE=1 anywhere in the column
      !$acc parallel loop gang vector collapse(3) wait(1) async(2)
      do j=jts,jte
         do k = kts, kte
            do i=its,ite
               LTRUE_COL(i,k,j)=maxval(LTRUE(i,:,j))
            end do
         end do
      end do
! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................
   !$acc parallel loop gang vector collapse(3) async(3) wait(2)
   do j=jts,jte      ! j loop (north-south)
      do k = kts, kte
         do i=its,ite      ! i loop (east-west)

       IF (LTRUE_COL(I,K,J).EQ.0) CYCLE
        DUMI(I,K,J) = QI3D(I,K,J)+QI3DTEN(I,K,J)*DT
        DUMQS(I,K,J) = QNI3D(I,K,J)+QNI3DTEN(I,K,J)*DT
        DUMR(I,K,J) = QR3D(I,K,J)+QR3DTEN(I,K,J)*DT
        DUMFNI(I,K,J) = NI3D(I,K,J)+NI3DTEN(I,K,J)*DT
        DUMFNS(I,K,J) = NS3D(I,K,J)+NS3DTEN(I,K,J)*DT
        DUMFNR(I,K,J) = NR3D(I,K,J)+NR3DTEN(I,K,J)*DT
        DUMC(I,K,J) = QC3D(I,K,J)+QC3DTEN(I,K,J)*DT
        DUMFNC(I,K,J) = NC3D(I,K,J)+NC3DTEN(I,K,J)*DT
	DUMG(I,K,J) = QG3D(I,K,J)+QG3DTEN(I,K,J)*DT
	DUMFNG(I,K,J) = NG3D(I,K,J)+NG3DTEN(I,K,J)*DT

! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (iinum.EQ.1) THEN
        DUMFNC(I,K,J) = NC3D(I,K,J)
        END IF

! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
      DUMFNI(I,K,J) = MAX(0.,DUMFNI(I,K,J))
      DUMFNS(I,K,J) = MAX(0.,DUMFNS(I,K,J))
      DUMFNC(I,K,J) = MAX(0.,DUMFNC(I,K,J))
      DUMFNR(I,K,J) = MAX(0.,DUMFNR(I,K,J))
      DUMFNG(I,K,J) = MAX(0.,DUMFNG(I,K,J))

!......................................................................
! CLOUD ICE

      IF (DUMI(I,K,J).GE.QSMALL) THEN
        DLAMI = (CONS12*DUMFNI(I,K,J)/DUMI(I,K,J))**(1./DI)
        DLAMI=MAX(DLAMI,LAMMINI)
        DLAMI=MIN(DLAMI,LAMMAXI)
      END IF
!......................................................................
! RAIN

      IF (DUMR(I,K,J).GE.QSMALL) THEN
        DLAMR = (PI*RHOW*DUMFNR(I,K,J)/DUMR(I,K,J))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      END IF
!......................................................................
! CLOUD DROPLETS

      IF (DUMC(I,K,J).GE.QSMALL) THEN
         DUM = PRES(I,K,J)/(287.15*T3D(I,K,J))
         PGAM(I,K,J)=0.0005714*(NC3D(I,K,J)/1.E6*DUM)+0.2714
         PGAM(I,K,J)=1./(PGAM(I,K,J)**2)-1.
         PGAM(I,K,J)=MAX(PGAM(I,K,J),2.)
         PGAM(I,K,J)=MIN(PGAM(I,K,J),10.)

        DLAMC = (CONS26*DUMFNC(I,K,J)*GAMMA(PGAM(I,K,J)+4.)/(DUMC(I,K,J)*GAMMA(PGAM(I,K,J)+1.)))**(1./3.)
        LAMMIN = (PGAM(I,K,J)+1.)/60.E-6
        LAMMAX = (PGAM(I,K,J)+1.)/1.E-6
        DLAMC=MAX(DLAMC,LAMMIN)
        DLAMC=MIN(DLAMC,LAMMAX)
      END IF
!......................................................................
! SNOW

      IF (DUMQS(I,K,J).GE.QSMALL) THEN
        DLAMS = (CONS1*DUMFNS(I,K,J)/ DUMQS(I,K,J))**(1./DS)
        DLAMS=MAX(DLAMS,LAMMINS)
        DLAMS=MIN(DLAMS,LAMMAXS)
      END IF
!......................................................................
! GRAUPEL

      IF (DUMG(I,K,J).GE.QSMALL) THEN
        DLAMG = (CONS2*DUMFNG(I,K,J)/ DUMG(I,K,J))**(1./DG)
        DLAMG=MAX(DLAMG,LAMMING)
        DLAMG=MIN(DLAMG,LAMMAXG)
      END IF

!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      IF (DUMC(I,K,J).GE.QSMALL) THEN
      UNC =  ACN(I,K,J)*GAMMA(1.+BC+PGAM(I,K,J))/ (DLAMC**BC*GAMMA(PGAM(I,K,J)+1.))
      UMC = ACN(I,K,J)*GAMMA(4.+BC+PGAM(I,K,J))/  (DLAMC**BC*GAMMA(PGAM(I,K,J)+4.))
      ELSE
      UMC = 0.
      UNC = 0.
      END IF

      IF (DUMI(I,K,J).GE.QSMALL) THEN
      UNI =  AIN(I,K,J)*CONS27/DLAMI**BI
      UMI = AIN(I,K,J)*CONS28/(DLAMI**BI)
      ELSE
      UMI = 0.
      UNI = 0.
      END IF

      IF (DUMR(I,K,J).GE.QSMALL) THEN
      UNR = ARN(I,K,J)*CONS6/DLAMR**BR
      UMR = ARN(I,K,J)*CONS4/(DLAMR**BR)
      ELSE
      UMR = 0.
      UNR = 0.
      END IF

      IF (DUMQS(I,K,J).GE.QSMALL) THEN
      UMS = ASN(I,K,J)*CONS3/(DLAMS**BS)
      UNS = ASN(I,K,J)*CONS5/DLAMS**BS
      ELSE
      UMS = 0.
      UNS = 0.
      END IF

      IF (DUMG(I,K,J).GE.QSMALL) THEN
      UMG = AGN(I,K,J)*CONS7/(DLAMG**BG)
      UNG = AGN(I,K,J)*CONS8/DLAMG**BG
      ELSE
      UMG = 0.
      UNG = 0.
      END IF

! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
        dum=(rhosu/RHO(I,K,J))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
! fix 053011
! fix for correction by AA 4/6/11
        UMI=MIN(UMI,1.2*(rhosu/RHO(I,K,J))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/RHO(I,K,J))**0.35)
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)

      FR(I,K,J) = UMR
      FI(I,K,J) = UMI
      FNI(I,K,J) = UNI
      FS(I,K,J) = UMS
      FNS(I,K,J) = UNS
      FNR(I,K,J) = UNR
      FC(I,K,J) = UMC
      FNC(I,K,J) = UNC
      FG(I,K,J) = UMG
      FNG(I,K,J) = UNG

      ENDDO
      ENDDO
      ENDDO
! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

   !$acc parallel loop gang collapse(2) async(4) wait(3)
   do j=jts,jte      ! j loop (north-south)
   do i=its,ite      ! i loop (east-west)

       IF (LTRUE_COL(I,KTS,J).EQ.0) CYCLE
      !$acc loop seq
      DO K = KTE-1,KTS,-1

	! IF (K.LE.KTE-1) THEN
        IF (FR(I,K,J).LT.1.E-10) THEN
	FR(I,K,J)=FR(I,K+1,J)
	END IF
        IF (FI(I,K,J).LT.1.E-10) THEN
	FI(I,K,J)=FI(I,K+1,J)
	END IF
        IF (FNI(I,K,J).LT.1.E-10) THEN
	FNI(I,K,J)=FNI(I,K+1,J)
	END IF
        IF (FS(I,K,J).LT.1.E-10) THEN
	FS(I,K,J)=FS(I,K+1,J)
	END IF
        IF (FNS(I,K,J).LT.1.E-10) THEN
	FNS(I,K,J)=FNS(I,K+1,J)
	END IF
        IF (FNR(I,K,J).LT.1.E-10) THEN
	FNR(I,K,J)=FNR(I,K+1,J)
	END IF
        IF (FC(I,K,J).LT.1.E-10) THEN
	FC(I,K,J)=FC(I,K+1,J)
	END IF
        IF (FNC(I,K,J).LT.1.E-10) THEN
	FNC(I,K,J)=FNC(I,K+1,J)
	END IF
        IF (FG(I,K,J).LT.1.E-10) THEN
	FG(I,K,J)=FG(I,K+1,J)
	END IF
        IF (FNG(I,K,J).LT.1.E-10) THEN
	FNG(I,K,J)=FNG(I,K+1,J)
	END IF
	! END IF ! K LE KTE-1

      END DO

      END DO
      END DO

! CALCULATE NUMBER OF SPLIT TIME STEPS
   !$acc parallel loop gang vector collapse(3) async(5) wait(4)
   do j=jts,jte      ! j loop (north-south)
   DO K = KTS,KTE
      do i=its,ite      ! i loop (east-west)

       IF (LTRUE_COL(I,K,J).EQ.0) CYCLE

      RGVM = MAX(FR(I,K,J),FI(I,K,J),FS(I,K,J),FC(I,K,J),FNI(I,K,J),FNR(I,K,J),FNS(I,K,J),FNC(I,K,J),FG(I,K,J),FNG(I,K,J))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEP(I,K,J) = MAX(INT(RGVM*DT/DZQ(I,K,J)+1.),1)

! MULTIPLY VARIABLES BY RHO(I,K,J)
      DUM1 = RHO(I,K,J)
      DUMR(I,K,J) = DUMR(I,K,J)*DUM1
      DUMI(I,K,J) = DUMI(I,K,J)*DUM1
      DUMFNI(I,K,J) = DUMFNI(I,K,J)*DUM1
      DUMQS(I,K,J) = DUMQS(I,K,J)*DUM1
      DUMFNS(I,K,J) = DUMFNS(I,K,J)*DUM1
      DUMFNR(I,K,J) = DUMFNR(I,K,J)*DUM1
      DUMC(I,K,J) = DUMC(I,K,J)*DUM1
      DUMFNC(I,K,J) = DUMFNC(I,K,J)*DUM1
      DUMG(I,K,J) = DUMG(I,K,J)*DUM1
      DUMFNG(I,K,J) = DUMFNG(I,K,J)*DUM1

      END DO
      END DO
      END DO

      !$acc parallel loop gang vector collapse(2) wait(5)
      do j = jts,jte
      do i = its,ite
            NSTEP_FLAT(i,j) = MAXVAL(NSTEP(i,:,j))
      enddo
      enddo

      MAXN = 0
      !$acc parallel loop gang vector collapse(3) copy(MAXN)
      do j = jts,jte
      do k = kts,kte
      do i = its,ite
            NSTEP(i,k,j) = NSTEP_FLAT(i,j)
            MAXN = max(MAXN, NSTEP(i,k,j))
      enddo
      enddo
      enddo


  DO N = 1,MAXN

   !$acc parallel loop gang vector collapse(3)
   do j=jts,jte      ! j loop (north-south)
   DO K = KTS,KTE
   do i=its,ite      ! i loop (east-west)

      IF (LTRUE_COL(I,K,J).EQ.0 .or. N > NSTEP(I,K,J)) cycle

      ! if (LTRUE(I,K,J).EQ.0) CYCLE  !NO HYDROMETEORS CALCULATED FOR THIS CELL
      FALOUTR(I,K,J) = FR(I,K,J)*DUMR(I,K,J)
      FALOUTI(I,K,J) = FI(I,K,J)*DUMI(I,K,J)
      FALOUTNI(I,K,J) = FNI(I,K,J)*DUMFNI(I,K,J)
      FALOUTS(I,K,J) = FS(I,K,J)*DUMQS(I,K,J)
      FALOUTNS(I,K,J) = FNS(I,K,J)*DUMFNS(I,K,J)
      FALOUTNR(I,K,J) = FNR(I,K,J)*DUMFNR(I,K,J)
      FALOUTC(I,K,J) = FC(I,K,J)*DUMC(I,K,J)
      FALOUTNC(I,K,J) = FNC(I,K,J)*DUMFNC(I,K,J)
      FALOUTG(I,K,J) = FG(I,K,J)*DUMG(I,K,J)
      FALOUTNG(I,K,J) = FNG(I,K,J)*DUMFNG(I,K,J)
      END DO
      END DO
      END DO

      !$acc parallel loop gang vector collapse(2)
      do j=jts,jte      ! j loop (north-south)
      do i=its,ite      ! i loop (east-west)

            IF (LTRUE_COL(I,KTE,J).EQ.0 .or. N > NSTEP(I,KTE,J)) cycle

            DUM1 = 1/DZQ(I,KTE,J)
            FALTNDR = FALOUTR(I,KTE,J)*DUM1
            FALTNDI = FALOUTI(I,KTE,J)*DUM1
            FALTNDNI = FALOUTNI(I,KTE,J)*DUM1
            FALTNDS = FALOUTS(I,KTE,J)*DUM1
            FALTNDNS = FALOUTNS(I,KTE,J)*DUM1
            FALTNDNR = FALOUTNR(I,KTE,J)*DUM1
            FALTNDC = FALOUTC(I,KTE,J)*DUM1
            FALTNDNC = FALOUTNC(I,KTE,J)*DUM1
            FALTNDG = FALOUTG(I,KTE,J)*DUM1
            FALTNDNG = FALOUTNG(I,KTE,J)*DUM1
      ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

            DUM2 = 1/NSTEP(I,KTE,J)/RHO(I,KTE,J)
            QRSTEN(I,KTE,J) = QRSTEN(I,KTE,J)-FALTNDR*DUM2
            QISTEN(I,KTE,J) = QISTEN(I,KTE,J)-FALTNDI*DUM2
            NI3DTEN(I,KTE,J) = NI3DTEN(I,KTE,J)-FALTNDNI*DUM2
            QNISTEN(I,KTE,J) = QNISTEN(I,KTE,J)-FALTNDS*DUM2
            NS3DTEN(I,KTE,J) = NS3DTEN(I,KTE,J)-FALTNDNS*DUM2
            NR3DTEN(I,KTE,J) = NR3DTEN(I,KTE,J)-FALTNDNR*DUM2
            QCSTEN(I,KTE,J) = QCSTEN(I,KTE,J)-FALTNDC*DUM2
            NC3DTEN(I,KTE,J) = NC3DTEN(I,KTE,J)-FALTNDNC*DUM2
            QGSTEN(I,KTE,J) = QGSTEN(I,KTE,J)-FALTNDG*DUM2
            NG3DTEN(I,KTE,J) = NG3DTEN(I,KTE,J)-FALTNDNG*DUM2

            DUMT = DT/NSTEP(I,KTE,J)
            DUMR(I,KTE,J) = DUMR(I,KTE,J)+FALTNDR*DUMT
            DUMI(I,KTE,J) = DUMI(I,KTE,J)+FALTNDI*DUMT
            DUMFNI(I,KTE,J) = DUMFNI(I,KTE,J)+FALTNDNI*DUMT
            DUMQS(I,KTE,J) = DUMQS(I,KTE,J)+FALTNDS*DUMT
            DUMFNS(I,KTE,J) = DUMFNS(I,KTE,J)+FALTNDNS*DUMT
            DUMFNR(I,KTE,J) = DUMFNR(I,KTE,J)+FALTNDNR*DUMT
            DUMC(I,KTE,J) = DUMC(I,KTE,J)+FALTNDC*DUMT
            DUMFNC(I,KTE,J) = DUMFNC(I,KTE,J)+FALTNDNC*DUMT
            DUMG(I,KTE,J) = DUMG(I,KTE,J)+FALTNDG*DUMT
            DUMFNG(I,KTE,J) = DUMFNG(I,KTE,J)+FALTNDNG*DUMT

      ENDDO
      enddo

      !$acc parallel loop gang vector collapse(3)
      do j=jts,jte      ! j loop (north-south)
      DO K = KTE-1,KTS,-1
      do i=its,ite      ! i loop (east-west)

            IF (LTRUE_COL(I,K,J).EQ.0 .or. N > NSTEP(I,K,J)) cycle

            DUM1 = 1/DZQ(I,K,J)

            FALTNDR = (FALOUTR(I,K+1,J)-FALOUTR(I,K,J))*DUM1
            FALTNDI = (FALOUTI(I,K+1,J)-FALOUTI(I,K,J))*DUM1
            FALTNDNI = (FALOUTNI(I,K+1,J)-FALOUTNI(I,K,J))*DUM1
            FALTNDS = (FALOUTS(I,K+1,J)-FALOUTS(I,K,J))*DUM1
            FALTNDNS = (FALOUTNS(I,K+1,J)-FALOUTNS(I,K,J))*DUM1
            FALTNDNR = (FALOUTNR(I,K+1,J)-FALOUTNR(I,K,J))*DUM1
            FALTNDC = (FALOUTC(I,K+1,J)-FALOUTC(I,K,J))*DUM1
            FALTNDNC = (FALOUTNC(I,K+1,J)-FALOUTNC(I,K,J))*DUM1
            FALTNDG = (FALOUTG(I,K+1,J)-FALOUTG(I,K,J))*DUM1
            FALTNDNG = (FALOUTNG(I,K+1,J)-FALOUTNG(I,K,J))*DUM1

      ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES
            DUM2 = 1/NSTEP(I,K,J)/RHO(I,K,J)

            QRSTEN(I,K,J) = QRSTEN(I,K,J)+FALTNDR*DUM2
            QISTEN(I,K,J) = QISTEN(I,K,J)+FALTNDI*DUM2
            NI3DTEN(I,K,J) = NI3DTEN(I,K,J)+FALTNDNI*DUM2
            QNISTEN(I,K,J) = QNISTEN(I,K,J)+FALTNDS*DUM2
            NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+FALTNDNS*DUM2
            NR3DTEN(I,K,J) = NR3DTEN(I,K,J)+FALTNDNR*DUM2
            QCSTEN(I,K,J) = QCSTEN(I,K,J)+FALTNDC*DUM2
            NC3DTEN(I,K,J) = NC3DTEN(I,K,J)+FALTNDNC*DUM2
            QGSTEN(I,K,J) = QGSTEN(I,K,J)+FALTNDG*DUM2
            NG3DTEN(I,K,J) = NG3DTEN(I,K,J)+FALTNDNG*DUM2

            DUMT = DT/NSTEP(I,K,J)
            DUMR(I,K,J) = DUMR(I,K,J)+FALTNDR*DUMT
            DUMI(I,K,J) = DUMI(I,K,J)+FALTNDI*DUMT
            DUMFNI(I,K,J) = DUMFNI(I,K,J)+FALTNDNI*DUMT
            DUMQS(I,K,J) = DUMQS(I,K,J)+FALTNDS*DUMT
            DUMFNS(I,K,J) = DUMFNS(I,K,J)+FALTNDNS*DUMT
            DUMFNR(I,K,J) = DUMFNR(I,K,J)+FALTNDNR*DUMT
            DUMC(I,K,J) = DUMC(I,K,J)+FALTNDC*DUMT
            DUMFNC(I,K,J) = DUMFNC(I,K,J)+FALTNDNC*DUMT
            DUMG(I,K,J) = DUMG(I,K,J)+FALTNDG*DUMT
            DUMFNG(I,K,J) = DUMFNG(I,K,J)+FALTNDNG*DUMT

      ! FOR WRF-CHEM, NEED PRECIP RATES (UNITS OF KG/M^2/S)
            CSED(I,K,J)=CSED(I,K,J)+FALOUTC(I,K,J)/NSTEP(I,K,J)
            ISED(I,K,J)=ISED(I,K,J)+FALOUTI(I,K,J)/NSTEP(I,K,J)
            SSED(I,K,J)=SSED(I,K,J)+FALOUTS(I,K,J)/NSTEP(I,K,J)
            GSED(I,K,J)=GSED(I,K,J)+FALOUTG(I,K,J)/NSTEP(I,K,J)
            RSED(I,K,J)=RSED(I,K,J)+FALOUTR(I,K,J)/NSTEP(I,K,J)
      END DO
      END DO
      END DO

! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
! OF LIQUID WATER CANCELS THIS FACTOR OF 1000
      
      !$acc parallel loop gang vector collapse(2)
      do j=jts,jte      ! j loop (north-south)
      do i=its,ite      ! i loop (east-west)
        IF (LTRUE_COL(I,KTS,J).EQ.0 .or. N > NSTEP(I,KTS,J)) cycle

        PRECPRT1D(I,J) = PRECPRT1D(I,J)+(FALOUTR(I,KTS,J)+FALOUTC(I,KTS,J)+FALOUTS(I,KTS,J)+FALOUTI(I,KTS,J)+FALOUTG(I,KTS,J))  &
                     *DT/NSTEP_FLAT(I,J)
        SNOWRT1D(I,J) = SNOWRT1D(I,J)+(FALOUTS(I,KTS,J)+FALOUTI(I,KTS,J)+FALOUTG(I,KTS,J))*DT/NSTEP_FLAT(I,J)
! hm added 7/13/13
        SNOWPRT1D(I,J) = SNOWPRT1D(I,J)+(FALOUTI(I,KTS,J)+FALOUTS(I,KTS,J))*DT/NSTEP_FLAT(I,J)
        GRPLPRT1D(I,J) = GRPLPRT1D(I,J)+(FALOUTG(I,KTS,J))*DT/NSTEP_FLAT(I,J)
      END DO
      END DO

      END DO

   !$acc parallel loop gang vector collapse(3) async(6) wait(5)
   do j=jts,jte      ! j loop (north-south)
   do k=kts,kte
   do i=its,ite      ! i loop (east-west)
        if (LTRUE(I,K,J).EQ.0) CYCLE !NO HYDROMETEORS CALCULATED FOR THIS CELL

! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

        QR3DTEN(I,K,J)=QR3DTEN(I,K,J)+QRSTEN(I,K,J)
        QI3DTEN(I,K,J)=QI3DTEN(I,K,J)+QISTEN(I,K,J)
        QC3DTEN(I,K,J)=QC3DTEN(I,K,J)+QCSTEN(I,K,J)
        QG3DTEN(I,K,J)=QG3DTEN(I,K,J)+QGSTEN(I,K,J)
        QNI3DTEN(I,K,J)=QNI3DTEN(I,K,J)+QNISTEN(I,K,J)

! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

!hm 4/7/09 bug fix
!        IF (QI3D(I,K,J).GE.QSMALL.AND.T3D(I,K,J).LT.273.15) THEN
        IF (QI3D(I,K,J).GE.QSMALL.AND.T3D(I,K,J).LT.273.15.AND.LAMI(I,K,J).GE.1.E-10) THEN
        IF (1./LAMI(I,K,J).GE.2.*DCS) THEN
           QNI3DTEN(I,K,J) = QNI3DTEN(I,K,J)+QI3D(I,K,J)/DT+ QI3DTEN(I,K,J)
           NS3DTEN(I,K,J) = NS3DTEN(I,K,J)+NI3D(I,K,J)/DT+   NI3DTEN(I,K,J)
           QI3DTEN(I,K,J) = -QI3D(I,K,J)/DT
           NI3DTEN(I,K,J) = -NI3D(I,K,J)/DT
        END IF
        END IF

! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration

          QC3D(I,K,J)        = QC3D(I,K,J)+QC3DTEN(I,K,J)*DT
          QI3D(I,K,J)        = QI3D(I,K,J)+QI3DTEN(I,K,J)*DT
          QNI3D(I,K,J)        = QNI3D(I,K,J)+QNI3DTEN(I,K,J)*DT
          QR3D(I,K,J)        = QR3D(I,K,J)+QR3DTEN(I,K,J)*DT
          NC3D(I,K,J)        = NC3D(I,K,J)+NC3DTEN(I,K,J)*DT
          NI3D(I,K,J)        = NI3D(I,K,J)+NI3DTEN(I,K,J)*DT
          NS3D(I,K,J)        = NS3D(I,K,J)+NS3DTEN(I,K,J)*DT
          NR3D(I,K,J)        = NR3D(I,K,J)+NR3DTEN(I,K,J)*DT

          IF (IGRAUP.EQ.0) THEN
          QG3D(I,K,J)        = QG3D(I,K,J)+QG3DTEN(I,K,J)*DT
          NG3D(I,K,J)        = NG3D(I,K,J)+NG3DTEN(I,K,J)*DT
          END IF

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
          T3D(I,K,J)         = T3D(I,K,J)+T3DTEN(I,K,J)*DT
          QV3D(I,K,J)        = QV3D(I,K,J)+QV3DTEN(I,K,J)*DT

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS = min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J),0))   ! PA
            EIS = min(0.99*pres(I,K,J),POLYSVP(T3D(I,K,J),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS.GT.EVS) EIS = EVS

            QVS = EP_2*EVS/(PRES(I,K,J)-EVS)
            QVI = EP_2*EIS/(PRES(I,K,J)-EIS)

            QVQVS = QV3D(I,K,J)/QVS
            QVQVSI = QV3D(I,K,J)/QVI

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm 7/9/09 change limit to 1.e-8

             IF (QVQVS.LT.0.9) THEN
               IF (QR3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QR3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QR3D(I,K,J)*XXLV(I,K,J)/CPM(I,K,J)
                  QR3D(I,K,J)=0.
               END IF
               IF (QC3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QC3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QC3D(I,K,J)*XXLV(I,K,J)/CPM(I,K,J)
                  QC3D(I,K,J)=0.
               END IF
             END IF

             IF (QVQVSI.LT.0.9) THEN
               IF (QI3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QI3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QI3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QI3D(I,K,J)=0.
               END IF
               IF (QNI3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QNI3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QNI3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QNI3D(I,K,J)=0.
               END IF
               IF (QG3D(I,K,J).LT.1.E-8) THEN
                  QV3D(I,K,J)=QV3D(I,K,J)+QG3D(I,K,J)
                  T3D(I,K,J)=T3D(I,K,J)-QG3D(I,K,J)*XXLS(I,K,J)/CPM(I,K,J)
                  QG3D(I,K,J)=0.
               END IF
             END IF

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(I,K,J).LT.QSMALL) THEN
         QC3D(I,K,J) = 0.
         NC3D(I,K,J) = 0.
         EFFC(I,K,J) = 0.
       END IF
       IF (QR3D(I,K,J).LT.QSMALL) THEN
         QR3D(I,K,J) = 0.
         NR3D(I,K,J) = 0.
         EFFR(I,K,J) = 0.
       END IF
       IF (QI3D(I,K,J).LT.QSMALL) THEN
         QI3D(I,K,J) = 0.
         NI3D(I,K,J) = 0.
         EFFI(I,K,J) = 0.
       END IF
       IF (QNI3D(I,K,J).LT.QSMALL) THEN
         QNI3D(I,K,J) = 0.
         NS3D(I,K,J) = 0.
         EFFS(I,K,J) = 0.
       END IF
       IF (QG3D(I,K,J).LT.QSMALL) THEN
         QG3D(I,K,J) = 0.
         NG3D(I,K,J) = 0.
         EFFG(I,K,J) = 0.
       END IF

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

      IF (.not.(QC3D(I,K,J).LT.QSMALL.AND.QI3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.QSMALL &
           .AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.QSMALL)) THEN

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN

        IF (QI3D(I,K,J).GE.QSMALL.AND.T3D(I,K,J).GE.273.15) THEN
           QR3D(I,K,J) = QR3D(I,K,J)+QI3D(I,K,J)
           T3D(I,K,J) = T3D(I,K,J)-QI3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
! #if (WRF_CHEM == 1)
!            tqimelt(I,K,J)=QI3D(I,K,J)/DT
! #endif
           QI3D(I,K,J) = 0.
           NR3D(I,K,J) = NR3D(I,K,J)+NI3D(I,K,J)
           NI3D(I,K,J) = 0.
        END IF

! ****SENSITIVITY - NO ICE
      IF (ILIQ.EQ.0) THEN

! HOMOGENEOUS FREEZING OF CLOUD WATER

        IF (T3D(I,K,J).LE.233.15.AND.QC3D(I,K,J).GE.QSMALL) THEN
           QI3D(I,K,J)=QI3D(I,K,J)+QC3D(I,K,J)
           T3D(I,K,J)=T3D(I,K,J)+QC3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
           QC3D(I,K,J)=0.
           NI3D(I,K,J)=NI3D(I,K,J)+NC3D(I,K,J)
           NC3D(I,K,J)=0.
        END IF

! HOMOGENEOUS FREEZING OF RAIN

        IF (IGRAUP.EQ.0) THEN

        IF (T3D(I,K,J).LE.233.15.AND.QR3D(I,K,J).GE.QSMALL) THEN
           QG3D(I,K,J) = QG3D(I,K,J)+QR3D(I,K,J)
           T3D(I,K,J) = T3D(I,K,J)+QR3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
           QR3D(I,K,J) = 0.
           NG3D(I,K,J) = NG3D(I,K,J)+ NR3D(I,K,J)
           NR3D(I,K,J) = 0.
        END IF

        ELSE IF (IGRAUP.EQ.1) THEN

        IF (T3D(I,K,J).LE.233.15.AND.QR3D(I,K,J).GE.QSMALL) THEN
           QNI3D(I,K,J) = QNI3D(I,K,J)+QR3D(I,K,J)
           T3D(I,K,J) = T3D(I,K,J)+QR3D(I,K,J)*XLF(I,K,J)/CPM(I,K,J)
           QR3D(I,K,J) = 0.
           NS3D(I,K,J) = NS3D(I,K,J)+NR3D(I,K,J)
           NR3D(I,K,J) = 0.
        END IF

        END IF

      ENDIF !(ILIQ.EQ.1) old CONTINUE 778

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(I,K,J) = MAX(0.,NI3D(I,K,J))
      NS3D(I,K,J) = MAX(0.,NS3D(I,K,J))
      NC3D(I,K,J) = MAX(0.,NC3D(I,K,J))
      NR3D(I,K,J) = MAX(0.,NR3D(I,K,J))
      NG3D(I,K,J) = MAX(0.,NG3D(I,K,J))

!......................................................................
! CLOUD ICE

      IF (QI3D(I,K,J).GE.QSMALL) THEN
         LAMI(I,K,J) = (CONS12*                 &
              NI3D(I,K,J)/QI3D(I,K,J))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(I,K,J).LT.LAMMINI) THEN

      LAMI(I,K,J) = LAMMINI

      N0I(I,K,J) = LAMI(I,K,J)**4*QI3D(I,K,J)/CONS12

      NI3D(I,K,J) = N0I(I,K,J)/LAMI(I,K,J)
      ELSE IF (LAMI(I,K,J).GT.LAMMAXI) THEN
      LAMI(I,K,J) = LAMMAXI
      N0I(I,K,J) = LAMI(I,K,J)**4*QI3D(I,K,J)/CONS12

      NI3D(I,K,J) = N0I(I,K,J)/LAMI(I,K,J)
      END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(I,K,J).GE.QSMALL) THEN
      LAMR(I,K,J) = (PI*RHOW*NR3D(I,K,J)/QR3D(I,K,J))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(I,K,J).LT.LAMMINR) THEN

      LAMR(I,K,J) = LAMMINR

      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      ELSE IF (LAMR(I,K,J).GT.LAMMAXR) THEN
      LAMR(I,K,J) = LAMMAXR
      N0RR(I,K,J) = LAMR(I,K,J)**4*QR3D(I,K,J)/(PI*RHOW)

      NR3D(I,K,J) = N0RR(I,K,J)/LAMR(I,K,J)
      END IF

      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(I,K,J).GE.QSMALL) THEN

         DUM = PRES(I,K,J)/(287.15*T3D(I,K,J))
         PGAM(I,K,J)=0.0005714*(NC3D(I,K,J)/1.E6*DUM)+0.2714
         PGAM(I,K,J)=1./(PGAM(I,K,J)**2)-1.
         PGAM(I,K,J)=MAX(PGAM(I,K,J),2.)
         PGAM(I,K,J)=MIN(PGAM(I,K,J),10.)

! CALCULATE LAMC

      LAMC(I,K,J) = (CONS26*NC3D(I,K,J)*GAMMA(PGAM(I,K,J)+4.)/   &
                 (QC3D(I,K,J)*GAMMA(PGAM(I,K,J)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(I,K,J)+1.)/60.E-6
      LAMMAX = (PGAM(I,K,J)+1.)/1.E-6

      IF (LAMC(I,K,J).LT.LAMMIN) THEN
      LAMC(I,K,J) = LAMMIN
      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26

      ELSE IF (LAMC(I,K,J).GT.LAMMAX) THEN
      LAMC(I,K,J) = LAMMAX
      NC3D(I,K,J) = EXP(3.*LOG(LAMC(I,K,J))+LOG(QC3D(I,K,J))+              &
                LOG(GAMMA(PGAM(I,K,J)+1.))-LOG(GAMMA(PGAM(I,K,J)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(I,K,J).GE.QSMALL) THEN
      LAMS(I,K,J) = (CONS1*NS3D(I,K,J)/QNI3D(I,K,J))**(1./DS)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(I,K,J).LT.LAMMINS) THEN
      LAMS(I,K,J) = LAMMINS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1

      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)

      ELSE IF (LAMS(I,K,J).GT.LAMMAXS) THEN

      LAMS(I,K,J) = LAMMAXS
      N0S(I,K,J) = LAMS(I,K,J)**4*QNI3D(I,K,J)/CONS1
      NS3D(I,K,J) = N0S(I,K,J)/LAMS(I,K,J)
      END IF

      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(I,K,J).GE.QSMALL) THEN
      LAMG(I,K,J) = (CONS2*NG3D(I,K,J)/QG3D(I,K,J))**(1./DG)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(I,K,J).LT.LAMMING) THEN
      LAMG(I,K,J) = LAMMING
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)

      ELSE IF (LAMG(I,K,J).GT.LAMMAXG) THEN

      LAMG(I,K,J) = LAMMAXG
      N0G(I,K,J) = LAMG(I,K,J)**4*QG3D(I,K,J)/CONS2

      NG3D(I,K,J) = N0G(I,K,J)/LAMG(I,K,J)
      END IF

      END IF

      ENDIF !(QC3D(I,K,J).LT.QSMALL.AND.QI3D(I,K,J).LT.QSMALL.AND.QNI3D(I,K,J).LT.QSMALL &
            !.AND.QR3D(I,K,J).LT.QSMALL.AND.QG3D(I,K,J).LT.QSMALL) old CONTINUE 500
! CALCULATE EFFECTIVE RADIUS

      IF (QI3D(I,K,J).GE.QSMALL) THEN
         EFFI(I,K,J) = 3./LAMI(I,K,J)/2.*1.E6
      ELSE
         EFFI(I,K,J) = 25.
      END IF

      IF (QNI3D(I,K,J).GE.QSMALL) THEN
         EFFS(I,K,J) = 3./LAMS(I,K,J)/2.*1.E6
      ELSE
         EFFS(I,K,J) = 25.
      END IF

      IF (QR3D(I,K,J).GE.QSMALL) THEN
         EFFR(I,K,J) = 3./LAMR(I,K,J)/2.*1.E6
      ELSE
         EFFR(I,K,J) = 25.
      END IF

      IF (QC3D(I,K,J).GE.QSMALL) THEN
            EFFC(I,K,J) = GAMMA(PGAM(I,K,J)+4.)/                        &
             GAMMA(PGAM(I,K,J)+3.)/LAMC(I,K,J)/2.*1.E6
      ELSE
      EFFC(I,K,J) = 25.
      END IF

      IF (QG3D(I,K,J).GE.QSMALL) THEN
         EFFG(I,K,J) = 3./LAMG(I,K,J)/2.*1.E6
      ELSE
         EFFG(I,K,J) = 25.
      END IF

! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
!          NI3D(I,K,J) = MIN(NI3D(I,K,J),10.E6/RHO(I,K,J))
! HM, 12/28/12, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
! OF EXCESSIVE AND PERSISTENT ANVIL
! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
          NI3D(I,K,J) = MIN(NI3D(I,K,J),0.3E6/RHO(I,K,J))

! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
          IF (iinum.EQ.0.AND.IACT.EQ.2) THEN
          NC3D(I,K,J) = MIN(NC3D(I,K,J),(NANEW1+NANEW2)/RHO(I,K,J))
          END IF
! SWITCH FOR CONSTANT DROPLET NUMBER
          IF (iinum.EQ.1) THEN
! CHANGE NDCNST FROM CM-3 TO KG-1
             NC3D(I,K,J) = NDCNST*1.E6/RHO(I,K,J)
          END IF

      END DO !!! K LOOP
      END DO
      END DO

! ALL DONE !!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF INLINED FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !$acc parallel loop gang vector collapse(3) async(7) wait(6)
   do j=jts,jte      ! j loop (north-south)
   do k=kts,kte
   do i=its,ite      ! i loop (east-west)

! hm, add tendencies to update global variables
! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,k,j)        = QC3D(I,K,J)
          QI(i,k,j)        = QI3D(I,K,J)
          QS(i,k,j)        = QNI3D(I,K,J)
          QR(i,k,j)        = QR3D(I,K,J)
          NI(i,k,j)        = NI3D(I,K,J)
          NS(i,k,j)        = NS3D(I,K,J)
          NR(i,k,j)        = NR3D(I,K,J)
	      QG(I,K,j)        = QG3D(I,K,J)
          NG(I,K,j)        = NG3D(I,K,J)

          TH(I,K,J)        = T3D(I,K,J)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          QV(i,k,j)        = QV3D(I,K,J)

          ISED3D(i,k,j) = ISED(I,K,J)
          SSED3D(i,k,j) = SSED(I,K,J)
! wrf-chem
        !   IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
        !      qndrop(i,k,j) = nc1d(I,K,J)
!jdf         CSED3D(I,K,J) = CSED(I,K,J)
        !   END IF
        !   IF ( PRESENT( QLSINK ) ) THEN
        !      if(qc(i,k,j)>1.e-10) then
        !         QLSINK(I,K,J)  = C2PREC(I,K,J)/QC(I,K,J)
        !      else
        !         QLSINK(I,K,J)  = 0.0
        !      endif
        !   END IF
        !   IF ( PRESENT( PRECR ) ) PRECR(I,K,J) = RSED(I,K,J)
        !   IF ( PRESENT( PRECI ) ) PRECI(I,K,J) = ISED(I,K,J)
        !   IF ( PRESENT( PRECS ) ) PRECS(I,K,J) = SSED(I,K,J)
        !   IF ( PRESENT( PRECG ) ) PRECG(I,K,J) = GSED(I,K,J)
! EFFECTIVE RADIUS FOR RADIATION CODE (currently not coupled)
! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
          EFFC(I,K,J)     = MIN(EFFC(I,K,J),50.)*1.E-6
          EFFC(I,K,J)     = MAX(EFFC(I,K,J),2.5)*1.E-6
          EFFI(I,K,J)     = MIN(EFFI(I,K,J),125.)*1.E-6
          EFFI(I,K,J)     = MAX(EFFI(I,K,J),5.)*1.E-6
          EFFS(I,K,J)     = MIN(EFFS(I,K,J),1000.)*1.E-6
          EFFS(I,K,J)     = MAX(EFFS(I,K,J),10.)*1.E-6
!          EFFCS(I,K,J)     = MIN(EFFC(I,K,J),50.)
!          EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),1.)
!          EFFIS(I,K,J)     = MIN(EFFI(I,K,J),130.)
!          EFFIS(I,K,J)     = MAX(EFFIS(I,K,J),13.)

! #if ( WRF_CHEM == 1)
!            IF ( PRESENT( rainprod ) ) rainprod(i,k,j) = rainprod1d(I,K,J)
!            IF ( PRESENT( evapprod ) ) evapprod(i,k,j) = evapprod1d(I,K,J)
! #endif

      end do

! hm modified so that m2005 precip variables correctly match wrf precip variables
!+---+-----------------------------------------------------------------+
        !  IF ( PRESENT (diagflag) ) THEN
        !  if (diagflag .and. do_radar_ref == 1) then
        !   call refl10cm_hm (QV3D, QR3D, NR3D, QNI3D, NS3D, QG3D, NG3D,   &
        !               t1d, PRES, dBZ, kts, kte, i, j)
        !   do k = kts, kte
        !      refl_10cm(i,k,j) = MAX(-35., dBZ(I,K,J))
        !   enddo
        !  endif
        !  ENDIF
!+---+-----------------------------------------------------------------+

   end do
   end do

   !$acc parallel loop gang vector collapse(2) wait(7)
   do j=jts,jte      ! j loop (north-south)
   do i=its,ite      ! i loop (east-west)
      RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D(i,j)
      RAINNCV(i,j) = PRECPRT1D(i,j)
! hm, added 7/13/13
      SNOWNC(i,j) = SNOWNC(I,J)+SNOWPRT1D(i,j)
      SNOWNCV(i,j) = SNOWPRT1D(i,j)
      GRAUPELNC(i,j) = GRAUPELNC(I,J)+GRPLPRT1D(i,j)
      GRAUPELNCV(i,j) = GRPLPRT1D(i,j)
      SR(i,j) = SNOWRT1D(i,j)/(PRECPRT1D(i,j)+1.E-12)
   enddo
   enddo
   ! !$omp end do
   ! !$omp end parallel
   !$acc end data

END SUBROUTINE MP_MORR_TWO_MOMENT_gpu

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FUNCTION POLYSVP (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
      !$acc routine seq
      IMPLICIT NONE

      REAL DUM
      REAL, INTENT(IN) :: T
      INTEGER, INTENT(IN) :: TYPE
! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1) THEN

!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!          LOG10(6.1071))*100.


      dt = max(-80.,t-273.16)
      polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
      polysvp = polysvp*100.

      END IF

! LIQUID

      IF (TYPE.EQ.0) THEN

       dt = max(-80.,t-273.16)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.

!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-									&
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!             LOG10(1013.246))*100.

         END IF


      END FUNCTION POLYSVP

!------------------------------------------------------------------------------

      REAL FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL, INTENT(IN) :: X
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


      REAL FUNCTION DERF1(X)
      IMPLICIT NONE
      REAL, INTENT(IN) :: X
      REAL, DIMENSION(0 : 64) :: A, B
      REAL W,T,Y
      INTEGER K,I
      DATA A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      DATA (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      DATA (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      DATA (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      DATA (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      DATA (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF1 = Y
      END FUNCTION DERF1

!+---+-----------------------------------------------------------------+

!       subroutine refl10cm_hm (QV3D, QR3D, NR3D, QNI3D, NS3D, QG3D, NG3D, &
!                       t1d, PRES, dBZ, kts, kte, ii, jj)
!
!       IMPLICIT NONE
!
! !..Sub arguments
!       INTEGER, INTENT(IN):: kts, kte, ii, jj
!       REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
!                       QV3D, QR3D, NR3D, QNI3D, NS3D, QG3D, NG3D, t1d, PRES
!       REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
!
! !..Local variables
!       REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
!       REAL, DIMENSION(kts:kte):: rr, nr, rs, ns, rg, ng
!
!       DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, ilams
!       DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_g, N0_s
!       DOUBLE PRECISION:: lamr, lamg, lams
!       LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg
!
!       REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
!       DOUBLE PRECISION:: fmelt_s, fmelt_g
!       DOUBLE PRECISION:: cback, x, eta, f_d
!
!       INTEGER:: i, k, k_0, kbot, n
!       LOGICAL:: melti
!
! !+---+
!
!       do k = kts, kte
!          dBZ(I,K,J) = -35.0
!       enddo
!
! !+---+-----------------------------------------------------------------+
! !..Put column of data into local arrays.
! !+---+-----------------------------------------------------------------+
!       do k = kts, kte
!          temp(I,K,J) = T3D(I,K,J)
!          qv(I,K,J) = MAX(1.E-10, QV3D(I,K,J))
!          pres(I,K,J) = PRES(I,K,J)
!          RHO = 0.622*pres(I,K,J)/(R*temp(I,K,J)*(qv(I,K,J)+0.622))
!
!          if (QR3D(I,K,J) .gt. 1.E-9) then
!             rr(I,K,J) = QR3D(I,K,J)*RHO
!             nr(I,K,J) = NR3D(I,K,J)*RHO
!             lamr = (xam_r*xcrg(3)*xorg2*nr(I,K,J)/rr(I,K,J))**xobmr
!             ilamr(I,K,J) = 1./lamr
!             N0_r(I,K,J) = nr(I,K,J)*xorg2*lamr**xcre(2)
!             L_qr(I,K,J) = .true.
!          else
!             rr(I,K,J) = 1.E-12
!             nr(I,K,J) = 1.E-12
!             L_qr(I,K,J) = .false.
!          endif
!
!          if (QNI3D(I,K,J) .gt. 1.E-9) then
!             rs(I,K,J) = QNI3D(I,K,J)*RHO
!             ns(I,K,J) = NS3D(I,K,J)*RHO
!             lams = (xam_s*xcsg(3)*xosg2*ns(I,K,J)/rs(I,K,J))**xobms
!             ilams(I,K,J) = 1./lams
!             N0_s(I,K,J) = ns(I,K,J)*xosg2*lams**xcse(2)
!             L_qs(I,K,J) = .true.
!          else
!             rs(I,K,J) = 1.E-12
!             ns(I,K,J) = 1.E-12
!             L_qs(I,K,J) = .false.
!          endif
!
!          if (QG3D(I,K,J) .gt. 1.E-9) then
!             rg(I,K,J) = QG3D(I,K,J)*RHO
!             ng(I,K,J) = NG3D(I,K,J)*RHO
!             lamg = (xam_g*xcgg(3)*xogg2*ng(I,K,J)/rg(I,K,J))**xobmg
!             ilamg(I,K,J) = 1./lamg
!             N0_g(I,K,J) = ng(I,K,J)*xogg2*lamg**xcge(2)
!             L_qg(I,K,J) = .true.
!          else
!             rg(I,K,J) = 1.E-12
!             ng(I,K,J) = 1.E-12
!             L_qg(I,K,J) = .false.
!          endif
!       enddo
!
! !+---+-----------------------------------------------------------------+
! !..Locate K-level of start of melting (k_0 is level above).
! !+---+-----------------------------------------------------------------+
!       melti = .false.
!       k_0 = kts
!       do k = kte-1, kts, -1
!          if ( (temp(I,K,J).gt.273.15) .and. L_qr(I,K,J)                         &
!                                   .and. (L_qs(I,K+1,J).or.L_qg(I,K+1,J)) ) then
!             k_0 = MAX(k+1, k_0)
!             melti=.true.
!             goto 195
!          endif
!       enddo
!  195  continue
!
! !+---+-----------------------------------------------------------------+
! !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
! !.. and non-water-coated snow and graupel when below freezing are
! !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
! !+---+-----------------------------------------------------------------+
!
!       do k = kts, kte
!          ze_rain(I,K,J) = 1.e-22
!          ze_snow(I,K,J) = 1.e-22
!          ze_graupel(I,K,J) = 1.e-22
!          if (L_qr(I,K,J)) ze_rain(I,K,J) = N0_r(I,K,J)*xcrg(4)*ilamr(I,K,J)**xcre(4)
!          if (L_qs(I,K,J)) ze_snow(I,K,J) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
!                                  * (xam_s/900.0)*(xam_s/900.0)          &
!                                  * N0_s(I,K,J)*xcsg(4)*ilams(I,K,J)**xcse(4)
!          if (L_qg(I,K,J)) ze_graupel(I,K,J) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
!                                     * (xam_g/900.0)*(xam_g/900.0)       &
!                                     * N0_g(I,K,J)*xcgg(4)*ilamg(I,K,J)**xcge(4)
!       enddo
!
! !+---+-----------------------------------------------------------------+
! !..Special case of melting ice (snow/graupel) particles.  Assume the
! !.. ice is surrounded by the liquid water.  Fraction of meltwater is
! !.. extremely simple based on amount found above the melting level.
! !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
! !.. routines).
! !+---+-----------------------------------------------------------------+
!
!       if (melti .and. k_0.ge.kts+1) then
!        do k = k_0-1, kts, -1
!
! !..Reflectivity contributed by melting snow
!           if (L_qs(I,K,J) .and. L_qs(k_0) ) then
!            fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(I,K,J)/rs(k_0), 0.99d0))
!            eta = 0.d0
!            lams = 1./ilams(I,K,J)
!            do n = 1, nrbins
!               x = xam_s * xxDs(n)**xbm_s
!               call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
!                     fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
!                     CBACK, mixingrulestring_s, matrixstring_s,          &
!                     inclusionstring_s, hoststring_s,                    &
!                     hostmatrixstring_s, hostinclusionstring_s)
!               f_d = N0_s(I,K,J)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
!               eta = eta + f_d * CBACK * simpson(n) * xdts(n)
!            enddo
!            ze_snow(I,K,J) = SNGL(lamda4 / (pi5 * K_w) * eta)
!           endif
!
!
! !..Reflectivity contributed by melting graupel
!
!           if (L_qg(I,K,J) .and. L_qg(k_0) ) then
!            fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(I,K,J)/rg(k_0), 0.99d0))
!            eta = 0.d0
!            lamg = 1./ilamg(I,K,J)
!            do n = 1, nrbins
!               x = xam_g * xxDg(n)**xbm_g
!               call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
!                     fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
!                     CBACK, mixingrulestring_g, matrixstring_g,          &
!                     inclusionstring_g, hoststring_g,                    &
!                     hostmatrixstring_g, hostinclusionstring_g)
!               f_d = N0_g(I,K,J)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
!               eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
!            enddo
!            ze_graupel(I,K,J) = SNGL(lamda4 / (pi5 * K_w) * eta)
!           endif
!
!        enddo
!       endif
!
!       do k = kte, kts, -1
!          dBZ(I,K,J) = 10.*log10((ze_rain(I,K,J)+ze_snow(I,K,J)+ze_graupel(I,K,J))*1.d18)
!       enddo
!
!
!       end subroutine refl10cm_hm
!
!+---+-----------------------------------------------------------------+

END MODULE module_mp_morr_two_moment_gpu
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
