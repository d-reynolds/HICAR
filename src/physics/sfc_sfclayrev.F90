!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_sfclayrev
  use icar_constants,    only : STD_OUT_PE

 REAL    , PARAMETER ::  VCONVC=1.
 REAL    , PARAMETER ::  CZO=0.0185
 REAL    , PARAMETER ::  OZO=1.59E-5
 REAL    , SAVE      ::  SBRLIM=250.
 
 REAL,   DIMENSION(0:1000 ),SAVE :: psim_stab,psim_unstab,psih_stab,psih_unstab
!$acc declare create(psim_stab,psim_unstab,psih_stab,psih_unstab)

CONTAINS

!-------------------------------------------------------------------
   SUBROUTINE SFCLAYREV(U3D,V3D,T3D,QV3D,P3D,dz8w,                    &
                     CP,G,ROVCP,R,XLV,PSFCPA,CHS,CHS2,CQS2,CPM,      &
                     ZNT,UST,PBLH,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                     FM,FH,                                        &
                     XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                     U10,V10,TH2,T2,Q2,                            &
                     GZ1OZ0,WSPD,BR,ISFFLX,DX,                     &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,EOMEG,STBOLT,                          &
                     P1000mb,                                      &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte,                    &
                     ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,          &
                     shalwater_z0,                                 & 
                     scm_force_flux                                )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!   Changes in V3.7 over water surfaces:
!          1. for ZNT/Cd, replacing constant OZO with 0.11*1.5E-5/UST(I,J)
!             the COARE 3.5 (Edson et al. 2013) formulation is also available
!          2. for VCONV, reducing magnitude by half
!          3. for Ck, replacing Carlson-Boland with COARE 3
!-------------------------------------------------------------------
!-- U3D         3D u-velocity interpolated to theta points (m/s)
!-- V3D         3D v-velocity interpolated to theta points (m/s)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- dz8w        dz between full levels (m)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G           acceleration due to gravity (m/s^2)
!-- ROVCP       R/CP
!-- R           gas constant for dry air (J/kg/K)
!-- XLV         latent heat of vaporization for water (J/kg)
!-- PSFCPA      surface pressure (Pa)
!-- ZNT         roughness length (m)
!-- UST         u* in similarity theory (m/s)
!-- USTM        u* in similarity theory (m/s) without vconv correction
!               used to couple with TKE scheme
!-- PBLH        PBL height from previous time (m)
!-- MAVAIL      surface moisture availability (between 0 and 1)
!-- ZOL         z/L height over Monin-Obukhov length
!-- MOL         T* (similarity theory) (K)
!-- REGIME      flag indicating PBL regime (stable, unstable, etc.)
!-- PSIM        similarity stability function for momentum
!-- PSIH        similarity stability function for heat
!-- FM          integrated stability function for momentum
!-- FH          integrated stability function for heat
!-- XLAND       land mask (1 for land, 2 for water)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          net upward latent heat flux at surface (W/m^2)
!-- TSK         surface temperature (K)
!-- FLHC        exchange coefficient for heat (W/m^2/K)
!-- FLQC        exchange coefficient for moisture (kg/m^2/s)
!-- CHS         heat/moisture exchange coefficient for LSM (m/s)
!-- QGH         lowest-level saturated mixing ratio
!-- QSFC        ground saturated mixing ratio
!-- U10         diagnostic 10m u wind
!-- V10         diagnostic 10m v wind
!-- TH2         diagnostic 2m theta (K)
!-- T2          diagnostic 2m temperature (K)
!-- Q2          diagnostic 2m mixing ratio (kg/kg)
!-- GZ1OZ0      log(z/z0) where z0 is roughness length
!-- WSPD        wind speed at lowest model level (m/s)
!-- BR          bulk Richardson number in surface layer
!-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
!-- DX          horizontal grid size (m)
!-- SVP1        constant for saturation vapor pressure (kPa)
!-- SVP2        constant for saturation vapor pressure (dimensionless)
!-- SVP3        constant for saturation vapor pressure (K)
!-- SVPT0       constant for saturation vapor pressure (K)
!-- EP1         constant for virtual temperature (R_v/R_d - 1) (dimensionless)
!-- EP2         constant for specific humidity calculation 
!               (R_d/R_v) (dimensionless)
!-- KARMAN      Von Karman constant
!-- EOMEG       angular velocity of earth's rotation (rad/s)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- isftcflx    =0, (Charnock and Carlson-Boland); =1, AHW Ck, Cd, =2 Garratt
!-- iz0tlnd     =0 Carlson-Boland, =1 Czil_new
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
      INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
!                                                               
      INTEGER,  INTENT(IN )   ::        ISFFLX
      REAL,     INTENT(IN )   ::        SVP1,SVP2,SVP3,SVPT0
      REAL,     INTENT(IN )   ::        EP1,EP2,KARMAN,EOMEG,STBOLT
      REAL,     INTENT(IN )   ::        P1000mb
!
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           dz8w
                                        
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           QV3D, &
                                                              P3D, &
                                                              T3D

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::             MAVAIL, &
                                                             PBLH, &
                                                              TSK

      INTEGER,     DIMENSION( ims:ime, jms:jme )                 , &
                INTENT(IN   )               ::             XLAND

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )               ::                U10, &
                                                              V10, &
                                                              TH2, &
                                                               T2, &
                                                               Q2

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)               ::             REGIME, &
                                                              HFX, &
                                                              QFX, &
                                                               LH, &
                                                             QSFC, &
                                                          MOL,RMOL
!m the following 5 are change to memory size
!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                 GZ1OZ0,WSPD,BR, &
                                                  PSIM,PSIH,FM,FH

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                            U3D, &
                                                              V3D
                                        
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::               PSFCPA

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                            ZNT, &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                             CHS2, &
                                                             CQS2, &
                                                              CHS

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                      FLHC,FLQC

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)   ::                                 &
                                                              QGH
                                    
      REAL,     INTENT(IN   )               ::   CP,G,ROVCP,R,XLV,DX
 
      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(OUT)     ::                  ck,cka,cd,cda

      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(INOUT)   ::                           USTM

      INTEGER,  OPTIONAL,  INTENT(IN )   ::     ISFTCFLX, IZ0TLND
      INTEGER,  OPTIONAL,  INTENT(IN )   ::     SCM_FORCE_FLUX

      INTEGER,  INTENT(IN )   ::     shalwater_z0 
!      REAL,     INTENT(IN )   ::     shalwater_depth ! not used in this routine...
!      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
!                INTENT(INOUT)   ::                               water_depth

! LOCAL VARS

      REAL,     DIMENSION( its:ite, jts:jte ) ::                       water_depth

      ! REAL,     DIMENSION( its:ite ) ::                    dz8w1d

! LOCAL VARS

      REAL,     DIMENSION( its:ite, jts:jte )        ::    RHOX, &
                                                           THX,QX, &
                                                           DENOMQ, &
                                                          DENOMQ2, &
                                                          DENOMT2

    REAL :: ZA, THVX, ZQKL, ZQKLP1, PSIH2, PSIM2, PSIH10, PSIM10, &
             WSPDI, GZ2OZ0, GZ10OZ0, GOVRTH, TGDSA, SCR3, SCR4, THGB, PSFC, &
             pq, pq2, pq10                                                       
!                                                      

      INTEGER                               ::                 KL

      INTEGER ::  N,I,J,K,KK,L,NZOL,NK,NZOL2,NZOL10

      REAL    ::  PL,THCON,TVCON,E1
      REAL    ::  ZL,TSKV,DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
      REAL    ::  DTG,PSIX,DTTHX,PSIX10,PSIT,PSIT2,PSIQ,PSIQ2,PSIQ10
      REAL    ::  FLUXC,VSGD,Z0Q,VISC,RESTAR,CZIL,GZ0OZQ,GZ0OZT
      REAL    ::  ZW, ZN1, ZN2
!
! .... paj ...
!
      REAL    :: zolzz,zol0
!     REAL    :: zolri,zolri2
!     REAL    :: psih_stable,psim_stable,psih_unstable,psim_unstable
!     REAL    :: psih_stable_full,psim_stable_full,psih_unstable_full,psim_unstable_full
      REAL    :: zl2,zl10,z0t

      REAL,     PARAMETER     ::        XKA=2.4E-5
      REAL,     PARAMETER     ::        PRT=1.

      ! DO J=jts,jte
      !   DO i=its,ite
      !     dz8w1d(I,J) = dz8w(i,1,j)
      !   ENDDO
   
      ! do J=jts,jte
      !   DO i=its,ite
      !      U1D(I,J) =U3D(i,1,j)
      !      V1D(I,J) =V3D(i,1,j)
      !      QV1D(I,J)=QV3D(i,1,j)
      !      P1D(I,J) =P3D(i,1,j)
      !      T1D(I,J) =T3D(i,1,j)
          !  water_depth(i,j) = 0.0
      !   ENDDO
      ! ENDDO

        !  Sending array starting locations of optional variables may cause
        !  troubles, so we explicitly change the call.

        ! CALL SFCLAYREV1D(J,U1D,V1D,T1D,QV1D,P1D,dz8w1d,               &
        !         CP,G,ROVCP,R,XLV,PSFC(ims,j),CHS(ims,j),CHS2(ims,j),&
        !         CQS2(ims,j),CPM(ims,j),PBLH(ims,j), RMOL(ims,j),   &
        !         ZNT(ims,j),UST(ims,j),MAVAIL(ims,j),ZOL(ims,j),    &
        !         MOL(ims,j),REGIME(ims,j),PSIM(ims,j),PSIH(ims,j),  &
        !         FM(ims,j),FH(ims,j),                               &
        !         XLAND(ims,j),HFX(ims,j),QFX(ims,j),TSK(ims,j),     &
        !         U10(ims,j),V10(ims,j),TH2(ims,j),T2(ims,j),        &
        !         Q2(ims,j),FLHC(ims,j),FLQC(ims,j),QGH(ims,j),      &
        !         QSFC(ims,j),LH(ims,j),                             &
        !         GZ1OZ0(ims,j),WSPD(ims,j),BR(ims,j),ISFFLX,DX,     &
        !         SVP1,SVP2,SVP3,SVPT0,EP1,EP2,KARMAN,EOMEG,STBOLT,  &
        !         P1000mb,                                           &
        !         shalwater_z0,water_depth,                          & 
        !         ids,ide, jds,jde, kds,kde,                         &
        !         ims,ime, jms,jme, kms,kme,                         &
        !         its,ite, jts,jte, kts,kte                          &
        !         ,isftcflx,iz0tlnd,scm_force_flux                   &
!#if ( EM_CORE == 1 )
!                USTM(ims,j),CK(ims,j),CKA(ims,j),                  &
!                CD(ims,j),CDA(ims,j)                               &
!#endif
                                                                  !  )
      ! ENDDO


  !  END SUBROUTINE SFCLAYREV


!-------------------------------------------------------------------
!    SUBROUTINE SFCLAYREV1D(J,UX,VX,T1D,QV1D,P1D,dz8w1d,                &
!                      CP,G,ROVCP,R,XLV,PSFCPA,CHS,CHS2,CQS2,CPM,PBLH,RMOL, &
!                      ZNT,UST,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH,FM,FH,&
!                      XLAND,HFX,QFX,TSK,                            &
!                      U10,V10,TH2,T2,Q2,FLHC,FLQC,QGH,              &
!                      QSFC,LH,GZ1OZ0,WSPD,BR,ISFFLX,DX,             &
!                      SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
!                      KARMAN,EOMEG,STBOLT,                          &
!                      P1000mb,                                      &
!                      shalwater_z0,water_depth,                     & 
!                      ids,ide, jds,jde, kds,kde,                    &
!                      ims,ime, jms,jme, kms,kme,                    &
!                      its,ite, jts,jte, kts,kte,                    &
!                      isftcflx, iz0tlnd,scm_force_flux,                            &
!                      ustm,ck,cka,cd,cda                            )
! !-------------------------------------------------------------------
!       IMPLICIT NONE
!-------------------------------------------------------------------

!       INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
!                                         ims,ime, jms,jme, kms,kme, &
!                                         its,ite, jts,jte, kts,kte, &
!                                         J
! !                                                               
!       INTEGER,  INTENT(IN )   ::        ISFFLX
!       REAL,     INTENT(IN )   ::        SVP1,SVP2,SVP3,SVPT0
!       REAL,     INTENT(IN )   ::        EP1,EP2,KARMAN,EOMEG,STBOLT
!       REAL,     INTENT(IN )   ::        P1000mb

! !
!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(IN   )               ::             MAVAIL, &
!                                                              PBLH, &
!                                                             XLAND, &
!                                                               TSK
! !
!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(IN   )               ::             PSFCPA

!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(INOUT)               ::             REGIME, &
!                                                               HFX, &
!                                                               QFX, &
!                                                          MOL,RMOL
! !m the following 5 are changed to memory size---
! !
!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(INOUT)   ::                 GZ1OZ0,WSPD,BR, &
!                                                   PSIM,PSIH,FM,FH

!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(INOUT)   ::                            ZNT, &
!                                                               ZOL, &
!                                                               UST, &
!                                                               CPM, &
!                                                              CHS2, &
!                                                              CQS2, &
!                                                               CHS

!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(INOUT)   ::                      FLHC,FLQC

!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(INOUT)   ::                                 &
!                                                          QSFC,QGH

!       REAL,     DIMENSION( ims:ime )                             , &
!                 INTENT(OUT)     ::                        U10,V10, &
!                                                      TH2,T2,Q2,LH

                                    
!       REAL,     INTENT(IN   )               ::   CP,G,ROVCP,R,XLV,DX

!       INTEGER,  INTENT(IN )   ::     shalwater_z0 
!       REAL,     DIMENSION( ims:ime ), INTENT(IN)  :: water_depth 
! ! MODULE-LOCAL VARIABLES, DEFINED IN SUBROUTINE SFCLAY
!       REAL,     DIMENSION( its:ite ),  INTENT(IN   )   ::  dz8w1d

!       REAL,     DIMENSION( its:ite ),  INTENT(IN   )   ::      UX, &
!                                                                VX, &
!                                                              QV1D, &
!                                                               P1D, &
!                                                               T1D


!-------------------------------------------------------------------
      KL=kte

    !$acc data present(  &
    !$acc     U3D,V3D,T3D,QV3D,P3D,dz8w,                        &
    !$acc     PSFCPA,CHS,CHS2,CQS2,CPM,      &
    !$acc     ZNT,UST,PBLH,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
    !$acc     FM,FH,                                        &
    !$acc     HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
    !$acc     U10,V10,TH2,T2,Q2,                            &
    !$acc     GZ1OZ0,WSPD,BR,                     &
    !$acc     ustm,ck,cka,cd,cda, xland)          &
    !$acc create(water_depth,THX,QX, &
    !$acc        DENOMQ,DENOMQ2,DENOMT2, &
    !$acc        RHOX)


    !$acc parallel loop gang vector collapse(2)
      do j = jts,jte
      DO i=its,ite
! PSFCPA cb
         PSFC=PSFCPA(I,J)/1000.
      ! ENDDO
!                                                      
!----CONVERT GROUND TEMPERATURE TO POTENTIAL TEMPERATURE:  
!                                                            
      ! DO I=its,ite                                   
        TGDSA=TSK(I,J)                                    
! PSFC cb
!        THGB=TSK(I,J)*(100./PSFC)**ROVCP                
        THGB=TSK(I,J)*(P1000mb/PSFCPA(I,J))**ROVCP   
      ! ENDDO                                               
!                                                            
!-----DECOUPLE FLUX-FORM VARIABLES TO GIVE U,V,T,THETA,THETA-VIR.,
!     T-VIR., QV, AND QC AT CROSS POINTS AND AT KTAU-1.  
!                                                                 
!     *** NOTE ***                                           
!         THE BOUNDARY WINDS MAY NOT BE ADEQUATELY AFFECTED BY FRICTION,         
!         SO USE ONLY INTERIOR VALUES OF UX AND VX TO CALCULATE 
!         TENDENCIES.                             
!                                                           

!     DO 24 I=its,ite
!        UX(I,J)=U1D(I,J)
!        VX(I,J)=V1D(I,J)
                                                   
!.....SCR3(I,K) STORE TEMPERATURE,                           
!     SCR4(I,K) STORE VIRTUAL TEMPERATURE.                                       
                                                                                 
      ! DO I=its,ite
! PL cb
         PL=P3D(I,1,J)/1000.
         SCR3=T3D(I,1,J)                                                   
!         THCON=(100./PL)**ROVCP                                                 
         THCON=(P1000mb*0.001/PL)**ROVCP
         THX(I,J)=SCR3*THCON                                               
         SCR4=SCR3                                                    
         THVX=THX(I,J)                                                     
         QX(I,J)=0.                                                             
  !  ENDDO                                                                 
!                                                                                
      ! DO I=its,ite
         water_depth(i,j) = 0.0

         QGH(I,J)=0.                                                                
         FLHC(I,J)=0.                                                               
         FLQC(I,J)=0.                                                               
         CPM(I,J)=CP                                                                
      ! ENDDO
!                                                                                
!     IF(IDRY.EQ.1)GOTO 80                                                   
      ! DO I=its,ite
         QX(I,J)=QV3D(I,1,J)                                                    
         TVCON=(1.+EP1*QX(I,J))                                      
         THVX=THX(I,J)*TVCON                                               
         SCR4=SCR3*TVCON                                              
      ! ENDDO                                                                 
!                                                                                
      ! DO I=its,ite
        E1=SVP1*EXP(SVP2*(TGDSA-SVPT0)/(TGDSA-SVP3))                       
!  for land points QSFC can come from previous time step
        !if(xland(I,J).gt.1.5.or.qsfc(I,J).le.0.0)
        QSFC(I,J)=EP2*E1/(PSFC-E1)                                                 
! QGH CHANGED TO USE LOWEST-LEVEL AIR TEMP CONSISTENT WITH MYJSFC CHANGE
! Q2SAT = QGH IN LSM
        E1=SVP1*EXP(SVP2*(T3D(I,1,J)-SVPT0)/(T3D(I,1,J)-SVP3))                       
        PL=P3D(I,1,J)/1000.
        QGH(I,J)=EP2*E1/(PL-E1)                                                 
        CPM(I,J)=CP*(1.+0.8*QX(I,J))                                   
      ! ENDDO                                                                   
  !  80 CONTINUE
                                                                                 
!-----COMPUTE THE HEIGHT OF FULL- AND HALF-SIGMA LEVELS ABOVE GROUND             
!     LEVEL, AND THE LAYER THICKNESSES.                                          
                                                                                 
      ! DO I=its,ite
        ZQKLP1=0.
        RHOX(I,J)=PSFC*1000./(R*SCR4)                                       
      ! ENDDO                                                                   
!                                                                                
      ! DO I=its,ite                                                   
           ZQKL=dz8w(I,1,J)+ZQKLP1
      ! ENDDO                                                                 
!                                                                                
      ! DO I=its,ite
         ZA=0.5*(ZQKL+ZQKLP1)                                        
      ! ENDDO                                                               
!                                                                                
      ! DO I=its,ite
        GOVRTH=G/THX(I,J)                                                    
      ! ENDDO                                                                   
                                                                                 
!-----CALCULATE BULK RICHARDSON NO. OF SURFACE LAYER, ACCORDING TO               
!     AKB(1976), EQ(12).                                                         
                   
      ! DO I=its,ite
        GZ1OZ0(I,J)=ALOG((ZA+ZNT(I,J))/ZNT(I,J))   ! log((z+z0)/z0)                                     
        GZ2OZ0=ALOG((2.+ZNT(I,J))/ZNT(I,J))      ! log((2+z0)/z0)                           
        GZ10OZ0=ALOG((10.+ZNT(I,J))/ZNT(I,J))    ! log((10+z0)z0)                    
        IF((XLAND(I,J)-1.5).GE.0)THEN                                            
          ZL=ZNT(I,J)                                                            
        ELSE                                                                     
          ZL=0.01                                                                
        ENDIF                                                                    
        WSPD(I,J)=SQRT(U3D(I,1,J)*U3D(I,1,J)+V3D(I,1,J)*V3D(I,1,J))                        

        TSKV=THGB*(1.+EP1*QSFC(I,J))                     
        DTHVDZ=(THVX-TSKV)                                                 
!  Convective velocity scale Vc and subgrid-scale velocity Vsg
!  following Beljaars (1994, QJRMS) and Mahrt and Sun (1995, MWR)
!                                ... HONG Aug. 2001
!
!       VCONV = 0.25*sqrt(g/tskv*pblh(I,J)*dthvm)
!      Use Beljaars over land, old MM5 (Wyngaard) formula over water
        if (xland(I,J).lt.1.5) then
        fluxc = max(hfx(I,J)/rhox(I,J)/cp                    &
              + ep1*tskv*qfx(I,J)/rhox(I,J),0.)
        VCONV = vconvc*(g/TGDSA*pblh(I,J)*fluxc)**.33
        else
        IF(-DTHVDZ.GE.0)THEN
          DTHVM=-DTHVDZ
        ELSE
          DTHVM=0.
        ENDIF
!       VCONV = 2.*SQRT(DTHVM)
! V3.7: reducing contribution in calm conditions
        VCONV = SQRT(DTHVM)
        endif
! Mahrt and Sun low-res correction
        VSGD = 0.32 * (max(dx/5000.-1.,0.))**.33
        WSPD(I,J)=SQRT(WSPD(I,J)*WSPD(I,J)+VCONV*VCONV+vsgd*vsgd)
        WSPD(I,J)=AMAX1(WSPD(I,J),0.1)
        BR(I,J)=GOVRTH*ZA*DTHVDZ/(WSPD(I,J)*WSPD(I,J))                        
!  IF PREVIOUSLY UNSTABLE, DO NOT LET INTO REGIMES 1 AND 2
        IF(MOL(I,J).LT.0.)BR(I,J)=AMIN1(BR(I,J),0.0)
!jdf
        RMOL(I,J)=-GOVRTH*DTHVDZ*ZA*KARMAN
!jdf

      ! ENDDO                                                                   

!                                                                                
!-----DIAGNOSE BASIC PARAMETERS FOR THE APPROPRIATED STABILITY CLASS:            
!                                                                                
!                                                                                
!     THE STABILITY CLASSES ARE DETERMINED BY BR (BULK RICHARDSON NO.)           
!     AND HOL (HEIGHT OF PBL/MONIN-OBUKHOV LENGTH).                              
!                                                                                
!     CRITERIA FOR THE CLASSES ARE AS FOLLOWS:                                   
!                                                                                
!        1. BR .GE. 0.0;                                                         
!               REPRESENTS NIGHTTIME STABLE CONDITIONS (REGIME=1),               
!                                                                                
!        3. BR .EQ. 0.0                                                          
!               REPRESENTS FORCED CONVECTION CONDITIONS (REGIME=3),              
!                                                                                
!        4. BR .LT. 0.0                                                          
!               REPRESENTS FREE CONVECTION CONDITIONS (REGIME=4).                
!                                                                                
!CCCCC                                                                           

      ! DO I=its,ite
!                                                                           
      if (br(I,J).gt.0) then
        if (br(I,J).gt.SBRLIM) then
        zol(I,J)=zolri(SBRLIM,ZA,ZNT(I,J))
        else
        zol(I,J)=zolri(br(I,J),ZA,ZNT(I,J))
        endif
      endif
!
      if (br(I,J).lt.0) then
       IF(UST(I,J).LT.0.001)THEN
          ZOL(I,J)=BR(I,J)*GZ1OZ0(I,J)
        ELSE
        if (br(I,J).lt.-SBRLIM) then
        zol(I,J)=zolri(-SBRLIM,ZA,ZNT(I,J))
        else
        zol(I,J)=zolri(br(I,J),ZA,ZNT(I,J))
        endif
       ENDIF
      endif!
! ... paj: compute integrated similarity functions.
!
        zolzz=zol(I,J)*(ZA+znt(I,J))/ZA ! (z+z0/L
        zol10=zol(I,J)*(10.+znt(I,J))/ZA   ! (10+z0)/L
        zol2=zol(I,J)*(2.+znt(I,J))/ZA     ! (2+z0)/L
        zol0=zol(I,J)*znt(I,J)/ZA          ! z0/L
        ZL2=(2.)/ZA*ZOL(I,J)             ! 2/L      
        ZL10=(10.)/ZA*ZOL(I,J)           ! 10/L

        IF((XLAND(I,J)-1.5).LT.0.)THEN
        ZL=(0.01)/ZA*ZOL(I,J)   ! (0.01)/L     
        ELSE
        ZL=ZOL0                     ! z0/L
        ENDIF
!                                                                                
!-----CLASS 1; STABLE (NIGHTTIME) CONDITIONS:                                    
!
       IF(BR(I,J).GT.0.) THEN
          REGIME(I,J)=1.
  !
  ! ... paj: psim and psih. Follows Cheng and Brutsaert 2005 (CB05).
  !
          psim(I,J)=psim_stable(zolzz)-psim_stable(zol0)
          psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
  !
          PSIM10=psim_stable(zol10)-psim_stable(zol0)
          PSIH10=psih_stable(zol10)-psih_stable(zol0)
  !
          PSIM2=psim_stable(zol2)-psim_stable(zol0)
          PSIH2=psih_stable(zol2)-psih_stable(zol0)
  !
  ! ... paj: preparations to compute PSIQ. Follows CB05+Carlson Boland JAM 1978.
  !
          pq=psih_stable(zol(I,J))-psih_stable(zl)
          pq2=psih_stable(zl2)-psih_stable(zl)
          pq10=psih_stable(zl10)-psih_stable(zl)
  !
  !       1.0 over Monin-Obukhov length
          RMOL(I,J)=ZOL(I,J)/ZA 
        ELSEIF(BR(I,J).EQ.0) THEN
  !                                                                                
  !                                                                                
  !-----CLASS 3; FORCED CONVECTION:                                                
  !                                                                                
          REGIME(I,J)=3.                                                           
          PSIM(I,J)=0.0                                                              
          PSIH(I,J)=PSIM(I,J)                                                          
          PSIM10=0.                                                   
          PSIH10=PSIM10                                           
          PSIM2=0.                                                  
          PSIH2=PSIM2                                           
  !
  ! paj: preparations to compute PSIQ.
  !
          pq=PSIH(I,J)
          pq2=PSIH2
          pq10=0.
  !
          ZOL(I,J)=0.                                             
          RMOL(I,J) = ZOL(I,J)/ZA  

        ELSEIF(BR(I,J).LT.0) THEN
  !                                                                                
  !-----CLASS 4; FREE CONVECTION:                                                  
  !                                                                                
          REGIME(I,J)=4.                                                           
  !
  ! ... paj: PSIM and PSIH ...
  !
          psim(I,J)=psim_unstable(zolzz)-psim_unstable(zol0)
          psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
  !
          PSIM10=psim_unstable(zol10)-psim_unstable(zol0)
          PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
  !
          PSIM2=psim_unstable(zol2)-psim_unstable(zol0)
          PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
  !
  ! ... paj: preparations to compute PSIQ 
  !
          pq=psih_unstable(zol(I,J))-psih_unstable(zl)
          pq2=psih_unstable(zl2)-psih_unstable(zl)
          pq10=psih_unstable(zl10)-psih_unstable(zl)
  !
  !---LIMIOT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND HIGH ROUGHNESS            
  !---  THIS PREVENTS DENOMINATOR IN FLUXES FROM GETTING TOO SMALL                 
          PSIH(I,J)=AMIN1(PSIH(I,J),0.9*GZ1OZ0(I,J))
          PSIM(I,J)=AMIN1(PSIM(I,J),0.9*GZ1OZ0(I,J))
          PSIH2=AMIN1(PSIH2,0.9*GZ2OZ0)
          PSIM10=AMIN1(PSIM10,0.9*GZ10OZ0)
  !
  ! AHW: mods to compute ck, cd
          PSIH10=AMIN1(PSIH10,0.9*GZ10OZ0)

          RMOL(I,J) = ZOL(I,J)/ZA  

      ENDIF     
      ! ENDDO                                                              
!                                                                                
!-----COMPUTE THE FRICTIONAL VELOCITY:                                           
!     ZA(1982) EQS(2.60),(2.61).                                                 
!                                                                                
      ! DO I=its,ite
        DTG=THX(I,J)-THGB                                                   
        PSIX=GZ1OZ0(I,J)-PSIM(I,J)                                                   
        PSIX10=GZ10OZ0-PSIM10

!     LOWER LIMIT ADDED TO PREVENT LARGE FLHC IN SOIL MODEL
!     ACTIVATES IN UNSTABLE CONDITIONS WITH THIN LAYERS OR HIGH Z0
!       PSIT=AMAX1(GZ1OZ0(I,J)-PSIH(I,J),2.)
       PSIT=GZ1OZ0(I,J)-PSIH(I,J)
       PSIT2=GZ2OZ0-PSIH2
!
        IF((XLAND(I,J)-1.5).GE.0)THEN                                            
          ZL=ZNT(I,J)                                                            
        ELSE                                                                     
          ZL=0.01                                                                
        ENDIF                                                                    
!
        PSIQ=ALOG(KARMAN*UST(I,J)*ZA/XKA+ZA/ZL)-pq
        PSIQ2=ALOG(KARMAN*UST(I,J)*2./XKA+2./ZL)-pq2

! AHW: mods to compute ck, cd
        PSIQ10=ALOG(KARMAN*UST(I,J)*10./XKA+10./ZL)-pq10

! V3.7: using Fairall 2003 to compute z0q and z0t over water:
!       adapted from module_sf_mynn.F
        IF ( (XLAND(I,J)-1.5).GE.0. ) THEN
              VISC=(1.32+0.009*(SCR3-273.15))*1.E-5
              RESTAR=UST(I,J)*ZNT(I,J)/VISC
              Z0T = (5.5e-5)*(RESTAR**(-0.60))
              Z0T = MIN(Z0T,1.0e-4)
              Z0T = MAX(Z0T,2.0e-9)
              Z0Q = Z0T

! following paj:
           zolzz=zol(I,J)*(ZA+z0t)/ZA    ! (z+z0t)/L
           zol10=zol(I,J)*(10.+z0t)/ZA   ! (10+z0t)/L
           zol2=zol(I,J)*(2.+z0t)/ZA     ! (2+z0t)/L
           zol0=zol(I,J)*z0t/ZA          ! z0t/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
              PSIT=ALOG((ZA+z0t)/Z0t)-PSIH(I,J)
              PSIT2=ALOG((2.+z0t)/Z0t)-PSIH2

           zolzz=zol(I,J)*(ZA+z0q)/ZA    ! (z+z0q)/L
           zol10=zol(I,J)*(10.+z0q)/ZA   ! (10+z0q)/L
           zol2=zol(I,J)*(2.+z0q)/ZA     ! (2+z0q)/L
           zol0=zol(I,J)*z0q/ZA          ! z0q/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
!
              PSIQ=ALOG((ZA+z0q)/Z0q)-PSIH(I,J)
              PSIQ2=ALOG((2.+z0q)/Z0q)-PSIH2
              PSIQ10=ALOG((10.+z0q)/Z0q)-PSIH10
        ENDIF

        IF ( PRESENT(ISFTCFLX) ) THEN
           IF ( ISFTCFLX.EQ.1 .AND. (XLAND(I,J)-1.5).GE.0. ) THEN
! v3.1
!             Z0Q = 1.e-4 + 1.e-3*(MAX(0.,UST(I,J)-1.))**2
! hfip1
!             Z0Q = 0.62*2.0E-5/UST(I,J) + 1.E-3*(MAX(0.,UST(I,J)-1.5))**2
! v3.2
              Z0Q = 1.e-4
!
! ... paj: recompute psih for z0q
!
           zolzz=zol(I,J)*(ZA+z0q)/ZA    ! (z+z0q)/L
           zol10=zol(I,J)*(10.+z0q)/ZA   ! (10+z0q)/L
           zol2=zol(I,J)*(2.+z0q)/ZA     ! (2+z0q)/L
           zol0=zol(I,J)*z0q/ZA          ! z0q/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
!
              PSIQ=ALOG((ZA+z0q)/Z0Q)-PSIH(I,J)
              PSIT=PSIQ
              PSIQ2=ALOG((2.+z0q)/Z0Q)-PSIH2
              PSIQ10=ALOG((10.+z0q)/Z0Q)-PSIH10
              PSIT2=PSIQ2
           ENDIF
          IF ( ISFTCFLX.EQ.2 .AND. (XLAND(I,J)-1.5).GE.0. ) THEN
! AHW: Garratt formula: Calculate roughness Reynolds number
!        Kinematic viscosity of air (linear approc to
!                 temp dependence at sea level)
! GZ0OZT and GZ0OZQ are based off formulas from Brutsaert (1975), which
! Garratt (1992) used with values of k = 0.40, Pr = 0.71, and Sc = 0.60
              VISC=(1.32+0.009*(SCR3-273.15))*1.E-5
!!            VISC=1.5E-5
              RESTAR=UST(I,J)*ZNT(I,J)/VISC
              GZ0OZT=0.40*(7.3*SQRT(SQRT(RESTAR))*SQRT(0.71)-5.)
!
! ... paj: compute psih for z0t for temperature ...
!
              z0t=znt(I,J)/exp(GZ0OZT)
!
           zolzz=zol(I,J)*(ZA+z0t)/ZA    ! (z+z0t)/L
           zol10=zol(I,J)*(10.+z0t)/ZA   ! (10+z0t)/L
           zol2=zol(I,J)*(2.+z0t)/ZA     ! (2+z0t)/L
           zol0=zol(I,J)*z0t/ZA          ! z0t/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
!
!              PSIT=GZ1OZ0(I,J)-PSIH(I,J)+RESTAR2
!              PSIT2=GZ2OZ0-PSIH2+RESTAR2
              PSIT=ALOG((ZA+z0t)/Z0t)-PSIH(I,J)
              PSIT2=ALOG((2.+z0t)/Z0t)-PSIH2
!
              GZ0OZQ=0.40*(7.3*SQRT(SQRT(RESTAR))*SQRT(0.60)-5.)
              z0q=znt(I,J)/exp(GZ0OZQ)
!
           zolzz=zol(I,J)*(ZA+z0q)/ZA    ! (z+z0q)/L
           zol10=zol(I,J)*(10.+z0q)/ZA   ! (10+z0q)/L
           zol2=zol(I,J)*(2.+z0q)/ZA     ! (2+z0q)/L
           zol0=zol(I,J)*z0q/ZA          ! z0q/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
!
              PSIQ=ALOG((ZA+z0q)/Z0q)-PSIH(I,J)
              PSIQ2=ALOG((2.+z0q)/Z0q)-PSIH2
              PSIQ10=ALOG((10.+z0q)/Z0q)-PSIH10
!              PSIQ=GZ1OZ0(I,J)-PSIH(I,J)+2.28*SQRT(SQRT(RESTAR))-2.
!              PSIQ2=GZ2OZ0-PSIH2+2.28*SQRT(SQRT(RESTAR))-2.
!              PSIQ10=GZ10OZ0-PSIH(I,J)+2.28*SQRT(SQRT(RESTAR))-2.
           ENDIF
        ENDIF
        IF(PRESENT(ck) .and. PRESENT(cd) .and. PRESENT(cka) .and. PRESENT(cda)) THEN
           Ck(I,J)=(karman/psix10)*(karman/psiq10)
           Cd(I,J)=(karman/psix10)*(karman/psix10)
           Cka(I,J)=(karman/psix)*(karman/psiq)
           Cda(I,J)=(karman/psix)*(karman/psix)
        ENDIF
        IF ( PRESENT(IZ0TLND) ) THEN
           IF ( IZ0TLND.GE.1 .AND. (XLAND(I,J)-1.5).LE.0. ) THEN
              ZL=ZNT(I,J)
!             CZIL RELATED CHANGES FOR LAND
              VISC=(1.32+0.009*(SCR3-273.15))*1.E-5
              RESTAR=UST(I,J)*ZL/VISC
!             Modify CZIL according to Chen & Zhang, 2009 if iz0tlnd = 1
!             If iz0tlnd = 2, use traditional value

              IF ( IZ0TLND.EQ.1 ) THEN
                 CZIL = 10.0 ** ( -0.40 * ( ZL / 0.07 ) )
              ELSE IF ( IZ0TLND.EQ.2 ) THEN
                 CZIL = 0.1 
              END IF
!
! ... paj: compute phish for z0t over land
!
              z0t=znt(I,J)/exp(CZIL*KARMAN*SQRT(RESTAR))
!
           zolzz=zol(I,J)*(ZA+z0t)/ZA    ! (z+z0t)/L
           zol10=zol(I,J)*(10.+z0t)/ZA   ! (10+z0t)/L
           zol2=zol(I,J)*(2.+z0t)/ZA     ! (2+z0t)/L
           zol0=zol(I,J)*z0t/ZA          ! z0t/L
!
              if (zol(I,J).gt.0.) then
              psih(I,J)=psih_stable(zolzz)-psih_stable(zol0)
              PSIH10=psih_stable(zol10)-psih_stable(zol0)
              PSIH2=psih_stable(zol2)-psih_stable(zol0)
              else
                if (zol(I,J).eq.0) then
                psih(I,J)=0.
                PSIH10=0.
                PSIH2=0.
                else
                psih(I,J)=psih_unstable(zolzz)-psih_unstable(zol0)
                PSIH10=psih_unstable(zol10)-psih_unstable(zol0)
                PSIH2=psih_unstable(zol2)-psih_unstable(zol0)
                endif
              endif
!
              PSIQ=ALOG((ZA+z0t)/Z0t)-PSIH(I,J)
              PSIQ2=ALOG((2.+z0t)/Z0t)-PSIH2
              PSIT=PSIQ
              PSIT2=PSIQ2
!
!              PSIT=GZ1OZ0(I,J)-PSIH(I,J)+CZIL*KARMAN*SQRT(RESTAR)
!              PSIQ=GZ1OZ0(I,J)-PSIH(I,J)+CZIL*KARMAN*SQRT(RESTAR)
!              PSIT2=GZ2OZ0-PSIH2+CZIL*KARMAN*SQRT(RESTAR)
!              PSIQ2=GZ2OZ0-PSIH2+CZIL*KARMAN*SQRT(RESTAR)

           ENDIF
        ENDIF
! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE 
        UST(I,J)=0.5*UST(I,J)+0.5*KARMAN*WSPD(I,J)/PSIX                                             
! TKE coupling: compute ust without vconv for use in tke scheme
        WSPDI=SQRT(U3D(I,1,J)*U3D(I,1,J)+V3D(I,1,J)*V3D(I,1,J))
        IF ( PRESENT(USTM) ) THEN
        USTM(I,J)=0.5*USTM(I,J)+0.5*KARMAN*WSPDI/PSIX
        ENDIF

        U10(I,J)=U3D(I,1,J)*PSIX10/PSIX                                    
        V10(I,J)=V3D(I,1,J)*PSIX10/PSIX                                   
        TH2(I,J)=THGB+DTG*PSIT2/PSIT                                
        Q2(I,J)=QSFC(I,J)+(QX(I,J)-QSFC(I,J))*PSIQ2/PSIQ                   
        T2(I,J) = TH2(I,J)*(PSFCPA(I,J)/P1000mb)**ROVCP                     
!                                                                                
        IF((XLAND(I,J)-1.5).LT.0.)THEN                                            
          UST(I,J)=AMAX1(UST(I,J),0.001)
        ENDIF                                                                    
        MOL(I,J)=KARMAN*DTG/PSIT/PRT                              
        DENOMQ(I,J)=PSIQ
        DENOMQ2(I,J)=PSIQ2
        DENOMT2(I,J)=PSIT2
        FM(I,J)=PSIX
        FH(I,J)=PSIT
  ! ENDDO                                                                   
!                                                                                
  ! 335 CONTINUE                                                                   
                                                                                  
!-----COMPUTE THE SURFACE SENSIBLE AND LATENT HEAT FLUXES:                       
      IF ( PRESENT(SCM_FORCE_FLUX) ) THEN
         IF (.NOT.(SCM_FORCE_FLUX.EQ.1)) THEN
          ! DO i=its,ite
            QFX(I,J)=0.                                                              
            HFX(I,J)=0.                                                              
          ! ENDDO
        ENDIF
      ENDIF
    ENDDO
    ENDDO

    IF (.not.(ISFFLX.EQ.0)) THEN                                           
                                                                                 
!-----OVER WATER, ALTER ROUGHNESS LENGTH (ZNT) ACCORDING TO WIND (UST).          
      !$acc parallel loop gang vector collapse(2)
      DO J=jts,jte                                              
      DO I=its,ite
        IF((XLAND(I,J)-1.5).GE.0)THEN                                            
!         ZNT(I,J)=CZO*UST(I,J)*UST(I,J)/G+OZO                                   
          ! PSH - formulation for depth-dependent roughness from
          ! ... Jimenez and Dudhia, 2018
          IF ( shalwater_z0 .eq. 1 ) THEN
             ZNT(I,J) = depth_dependent_z0(water_depth(I,J),ZNT(I,J),UST(I,J))
          ELSE
             ! Since V3.7 (ref: EC Physics document for Cy36r1)
             ZNT(I,J)=CZO*UST(I,J)*UST(I,J)/G+0.11*1.5E-5/UST(I,J)
             ! V3.9: Add limit as in isftcflx = 1,2
             ZNT(I,J)=MIN(ZNT(I,J),2.85e-3)
          ENDIF
! COARE 3.5 (Edson et al. 2013)
!         CZC = 0.0017*WSPD(I,J)-0.005
!         CZC = min(CZC,0.028)
!         ZNT(I,J)=CZC*UST(I,J)*UST(I,J)/G+0.11*1.5E-5/UST(I,J)
! AHW: change roughness length, and hence the drag coefficients Ck and Cd
          IF ( PRESENT(ISFTCFLX) ) THEN
             IF ( ISFTCFLX.NE.0 ) THEN
!               ZNT(I,J)=10.*exp(-9.*UST(I,J)**(-.3333))
!               ZNT(I,J)=10.*exp(-9.5*UST(I,J)**(-.3333))
!               ZNT(I,J)=ZNT(I,J) + 0.11*1.5E-5/AMAX1(UST(I,J),0.01)
!               ZNT(I,J)=0.011*UST(I,J)*UST(I,J)/G+OZO
!               ZNT(I,J)=MAX(ZNT(I,J),3.50e-5)
! AHW 2012:
                ZW  = MIN((UST(I,J)/1.06)**(0.3),1.0)
                ZN1 = 0.011*UST(I,J)*UST(I,J)/G + OZO
                ZN2 = 10.*exp(-9.5*UST(I,J)**(-.3333)) + &
                       0.11*1.5E-5/AMAX1(UST(I,J),0.01)
                ZNT(I,J)=(1.0-ZW) * ZN1 + ZW * ZN2
                ZNT(I,J)=MIN(ZNT(I,J),2.85e-3)
                ZNT(I,J)=MAX(ZNT(I,J),1.27e-7)
             ENDIF
          ENDIF
          ZL = ZNT(I,J)
        ELSE
          ZL = 0.01
        ENDIF                                                                    
        FLQC(I,J)=RHOX(I,J)*MAVAIL(I,J)*UST(I,J)*KARMAN/DENOMQ(I,J)
!       FLQC(I,J)=RHOX(I,J)*MAVAIL(I,J)*UST(I,J)*KARMAN/(   &
!               ALOG(KARMAN*UST(I,J)*ZA/XKA+ZA/ZL)-PSIH(I,J))
        DTTHX=ABS(THX(I,J)-THGB)                                            
        IF(DTTHX.GT.1.E-5)THEN                                                   
          FLHC(I,J)=CPM(I,J)*RHOX(I,J)*UST(I,J)*MOL(I,J)/(THX(I,J)-THGB)          
!         write(*,1001)FLHC(I,J),CPM(I,J),RHOX(I,J),UST(I,J),MOL(I,J),THX(I,J),THGB,I
 1001   format(f8.5,2x,f12.7,2x,f12.10,2x,f12.10,2x,f13.10,2x,f12.8,f12.8,2x,i3)
        ELSE                                                                     
          FLHC(I,J)=0.                                                             
        ENDIF                                                                    !                                                                                
!-----COMPUTE SURFACE MOIST FLUX:                                                
!                                                                                
!                                                                                
     IF ( PRESENT(SCM_FORCE_FLUX) ) THEN
        IF (.not.(SCM_FORCE_FLUX.EQ.1)) THEN
            QFX(I,J)=FLQC(I,J)*(QSFC(I,J)-QX(I,J))                                     
            QFX(I,J)=AMAX1(QFX(I,J),0.)                                            
            LH(I,J)=XLV*QFX(I,J)
          ! ENDDO                                                                  
    !-----COMPUTE SURFACE HEAT FLUX:                                                 
    !                                                                                
          ! DO I=its,ite
            IF(XLAND(I,J)-1.5.GT.0.)THEN                                           
              HFX(I,J)=FLHC(I,J)*(THGB-THX(I,J)) 
    !         IF ( PRESENT(ISFTCFLX) ) THEN
    !            IF ( ISFTCFLX.NE.0 ) THEN
    ! AHW: add dissipative heating term (commented out in 3.6.1)
    !               HFX(I,J)=HFX(I,J)+RHOX(I,J)*USTM(I,J)*USTM(I,J)*WSPDI
    !            ENDIF
    !         ENDIF 
            ELSEIF(XLAND(I,J)-1.5.LT.0.)THEN                                       
              HFX(I,J)=FLHC(I,J)*(THGB-THX(I,J))                                
              HFX(I,J)=AMAX1(HFX(I,J),-250.)                                       
            ENDIF         
      ENDIF
    ENDIF

         IF((XLAND(I,J)-1.5).GE.0)THEN
           ZL=ZNT(I,J)
         ELSE
           ZL=0.01
         ENDIF
!v3.1.1
!         CHS(I,J)=UST(I,J)*KARMAN/(ALOG(KARMAN*UST(I,J)*ZA &
!                /XKA+ZA/ZL)-PSIH(I,J))
         CHS(I,J)=UST(I,J)*KARMAN/DENOMQ(I,J)
!        GZ2OZ0=ALOG(2./ZNT(I,J))
!        PSIM2=-10.*GZ2OZ0
!        PSIM2=AMAX1(PSIM2,-10.)
!        PSIH2=PSIM2
! v3.1.1
!         CQS2(I,J)=UST(I,J)*KARMAN/(ALOG(KARMAN*UST(I,J)*2.0  &
!               /XKA+2.0/ZL)-PSIH2)
!         CHS2(I,J)=UST(I,J)*KARMAN/(GZ2OZ0-PSIH2)
         CQS2(I,J)=UST(I,J)*KARMAN/DENOMQ2(I,J)
         CHS2(I,J)=UST(I,J)*KARMAN/DENOMT2(I,J)
      ENDDO
      ENDDO                                                     
  ! 410 CONTINUE                                                                   
!jdf
!     DO I=its,ite
!       IF(UST(I,J).GE.0.1) THEN
!         RMOL(I,J)=RMOL(I,J)*(-FLHC(I,J))/(UST(I,J)*UST(I,J)*UST(I,J))
!       ELSE
!         RMOL(I,J)=RMOL(I,J)*(-FLHC(I,J))/(0.1*0.1*0.1)
!       ENDIF
!     ENDDO
!jdf

!                                 
    ENDIF

  !$acc end data

   END SUBROUTINE SFCLAYREV

!====================================================================
   SUBROUTINE sfclayrevinit(ims,ime,jms,jme,                    &
                            its,ite,jts,jte, sbrlim_in)!,        &
!                            bathymetry_flag, shalwater_z0,      &
!                            shalwater_depth, water_depth,       &
!                            xland,LakeModel,lake_depth,lakemask )

    INTEGER                   ::      N
    REAL                      ::      zolf

   INTEGER, INTENT(IN)      ::   ims,ime,jms,jme,its,ite,jts,jte
   REAL, OPTIONAL, INTENT(IN)      ::   sbrlim_in

!   INTEGER, INTENT(IN)      ::   shalwater_z0 
!   REAL,    INTENT(IN)      ::   shalwater_depth 
!   INTEGER, INTENT(IN)      ::   bathymetry_flag 
!   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)    ::  water_depth
!   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   )    ::  xland 
!   INTEGER         ::     LakeModel
!   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   )    ::  lake_depth
!   REAL,    DIMENSION( ims:ime, jms:jme )                   ::  lakemask  

    if (present(sbrlim_in)) SBRLIM = sbrlim_in

    !$acc parallel loop
    DO N=0,1000
! stable function tables
       zolf = float(n)*0.01
       psim_stab(n)=psim_stable_full(zolf)
       psih_stab(n)=psih_stable_full(zolf)
 
! unstable function tables
       zolf = -float(n)*0.01
       psim_unstab(n)=psim_unstable_full(zolf)
       psih_unstab(n)=psih_unstable_full(zolf)

    ENDDO

!    IF ( shalwater_z0 .EQ. 1 ) THEN
!       CALL shalwater_init(ims,ime,jms,jme,            &
!                 its,ite,jts,jte,                      &
!                 bathymetry_flag, shalwater_z0,        &
!                 shalwater_depth, water_depth,         &
!                 xland,LakeModel,lake_depth,lakemask   )
!    END IF

   END SUBROUTINE sfclayrevinit

   SUBROUTINE shalwater_init(ims,ime,jms,jme,                    &
                             its,ite,jts,jte,                    &
                             bathymetry_flag, shalwater_z0,      &
                             shalwater_depth, water_depth,       &
                             xland,LakeModel,lake_depth,lakemask )

   INTEGER, INTENT(IN)      ::   ims,ime,jms,jme,its,ite,jts,jte
   INTEGER, INTENT(IN)      ::   shalwater_z0 
   REAL,    INTENT(IN)      ::   shalwater_depth 
   INTEGER, INTENT(IN)      ::   bathymetry_flag 
   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)    ::  water_depth
   INTEGER,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   )    ::  xland 
   INTEGER         ::     LakeModel
   REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   )    ::  lake_depth
   REAL,    DIMENSION( ims:ime, jms:jme )                   ::  lakemask  

   ! Local 
   LOGICAL :: overwrite_water_depth
   INTEGER :: i, j
   overwrite_water_depth = .False.

   IF ( bathymetry_flag .eq. 1 ) THEN
      IF ( shalwater_depth .LE. 0.0 ) THEN
         IF ( LakeModel .ge. 1 ) THEN
            DO j = jts,jte
               DO i = its,ite
                  IF ( lakemask(i,j) .EQ. 1 ) THEN
                     water_depth(i,j) = lake_depth(i,j)
                  END IF
               END DO
            END DO
         END IF
      ELSE
         overwrite_water_depth = .True.
      END IF
   ELSE
      IF ( shalwater_depth .GT. 0.0 ) THEN
         overwrite_water_depth = .True.
      ELSE
         if (STD_OUT_PE) write(*,*) "No bathymetry data detected and shalwater_depth not greater than 0.0."
         stop
      END IF
   END IF

   IF (overwrite_water_depth) THEN 
      DO j = jts,jte
         DO i = its,ite
            IF((XLAND(i,j)-1.5).GE.0)THEN
               water_depth(i,j) = shalwater_depth
            ELSE
               water_depth(i,j) = -2.0 
            END IF
         END DO
      END DO
   END IF

   END SUBROUTINE shalwater_init

      function zolri(ri,z,z0)
      !$acc routine seq

      real, intent(in) :: ri
      real, intent(in)    :: z, z0
      real    :: x1, x2, fx1, fx2, zolri
      integer :: iter

!
      if (ri.lt.0.)then
        x1=-5.
        x2=0.
      else
        x1=0.
        x2=5.
      endif
!
      fx1=zolri2(x1,ri,z,z0)
      fx2=zolri2(x2,ri,z,z0)
      iter = 0
      Do While (abs(x1 - x2) > 0.01)
      if (iter .eq. 10) return
! check added for potential divide by zero (2019/11)
      if(fx1.eq.fx2)return
      if(abs(fx2).lt.abs(fx1))then
        x1=x1-fx1/(fx2-fx1)*(x2-x1)
        fx1=zolri2(x1,ri,z,z0)
        zolri=x1
      else
        x2=x2-fx2/(fx2-fx1)*(x2-x1)
        fx2=zolri2(x2,ri,z,z0)
        zolri=x2
      endif
!
      iter = iter + 1
      enddo
!

      return
      end function

!
! -----------------------------------------------------------------------
!
      function zolri2(zol2,ri2,z,z0)
      !$acc routine seq
      real, intent(inout) :: zol2
      real, intent(in)    :: ri2, z, z0
      real  :: zol3, zol20, psix2, psih2, zolri2
!
      if(zol2*ri2 .lt. 0.)zol2=0.  ! limit zol2 - must be same sign as ri2
!
      zol20=zol2*z0/z ! z0/L
      zol3=zol2+zol20 ! (z+z0)/L
!
      if (ri2.lt.0) then
      psix2=log((z+z0)/z0)-(psim_unstable(zol3)-psim_unstable(zol20))
      psih2=log((z+z0)/z0)-(psih_unstable(zol3)-psih_unstable(zol20))
      else
      psix2=log((z+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
      psih2=log((z+z0)/z0)-(psih_stable(zol3)-psih_stable(zol20))
      endif
!
      zolri2=zol2*psih2/psix2**2-ri2
!
      return
      end function
!
! ... integrated similarity functions ...
!
      function psim_stable_full(zolf)
      !$acc routine seq
      real   :: psim_stable_full
      real, intent(in)    :: zolf

        psim_stable_full=-6.1*log(zolf+(1+zolf**2.5)**(1./2.5))
      return
      end function

      function psih_stable_full(zolf)
      real   :: psih_stable_full
      real, intent(in)    :: zolf

        psih_stable_full=-5.3*log(zolf+(1+zolf**1.1)**(1./1.1))
      return
      end function
      
      function psim_unstable_full(zolf)
      !$acc routine seq
      real   :: x, psimk, ym, psimc, psim_unstable_full
      real, intent(in)    :: zolf

        x=(1.-16.*zolf)**.25
        psimk=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*ATAN(1.)
!
        ym=(1.-10.*zolf)**0.33
        psimc=(3./2.)*log((ym**2.+ym+1.)/3.)-sqrt(3.)*ATAN((2.*ym+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
!
        psim_unstable_full=(psimk+zolf**2*(psimc))/(1+zolf**2.)

      return
      end function

      function psih_unstable_full(zolf)
      !$acc routine seq
      real   :: y, psihk, yh, psihc, psih_unstable_full
      real, intent(in)    :: zolf

        y=(1.-16.*zolf)**.5
        psihk=2.*log((1+y)/2.)
!
        yh=(1.-34.*zolf)**0.33
        psihc=(3./2.)*log((yh**2.+yh+1.)/3.)-sqrt(3.)*ATAN((2.*yh+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
!
        psih_unstable_full=(psihk+zolf**2*(psihc))/(1+zolf**2.)

      return
      end function

! look-up table functions
      function psim_stable(zolf)
      !$acc routine seq
      integer :: nzol
      real    :: rzol, psim_stable
      real, intent(in)    :: zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_stable = psim_stab(nzol) + rzol*(psim_stab(nzol+1)-psim_stab(nzol))
        else
           psim_stable = psim_stable_full(zolf)
        endif
      return
      end function

      function psih_stable(zolf)
      !$acc routine seq
      integer :: nzol
      real    :: rzol, psih_stable
      real, intent(in)    :: zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_stable = psih_stab(nzol) + rzol*(psih_stab(nzol+1)-psih_stab(nzol))
        else
           psih_stable = psih_stable_full(zolf)
        endif
      return
      end function
      
      function psim_unstable(zolf)
      !$acc routine seq
      integer :: nzol
      real    :: rzol, psim_unstable
      real, intent(in)    :: zolf
      
        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_unstable = psim_unstab(nzol) + rzol*(psim_unstab(nzol+1)-psim_unstab(nzol))
        else
           psim_unstable = psim_unstable_full(zolf)
        endif
      return
      end function

      function psih_unstable(zolf)
      !$acc routine seq
      integer :: nzol
      real    :: rzol, psih_unstable
      real, intent(in)    :: zolf
      
        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_unstable = psih_unstab(nzol) + rzol*(psih_unstab(nzol+1)-psih_unstab(nzol))
        else
           psih_unstable = psih_unstable_full(zolf)
        endif
      return
      end function

      function depth_dependent_z0(water_depth,z0,UST)
      real :: depth_b, effective_depth, depth_dependent_z0
      real, intent(in):: water_depth, z0, UST
      
         IF ( water_depth .lt. 10.0 ) THEN
           effective_depth = 10.0
         ELSEIF ( water_depth .gt. 100.0 ) THEN
           effective_depth = 100.0
         ELSE
           effective_depth = water_depth
         ENDIF
 
         depth_b = 1 / 30.0 * log (1260.0 / effective_depth)
         depth_dependent_z0 = exp((2.7 * ust - 1.8 / depth_b) / (ust + 0.17 / depth_b) )
         depth_dependent_z0 = MIN(depth_dependent_z0,0.1)
     return
     end function
!-------------------------------------------------------------------          

END MODULE module_sf_sfclayrev

!
! ----------------------------------------------------------
!