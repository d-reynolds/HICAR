!-----------------------------------------------------------------------
! Factorial Snow Model interface 
!
!-----------------------------------------------------------------------
module FSM_interface

	use LANDUSE

    use MODCONF, only: SNTRAN,SNSLID

!-----------------------------------------------------------------------
! Meteorological driving variables
!-----------------------------------------------------------------------
    use DRIVING, only: &
      year,              &! Year
      month,             &! Month of year
      day,               &! Day of month
      hour,              &! Hour of day
      dt,                &! Timestep (s)
      LW,                &! Incoming longwave radiation (W/m2)
      Ps,                &! Surface pressure (Pa)
      Qa,                &! Specific humidity (kg/kg)
      Rf,                &! Rainfall rate (kg/m2/s)
      Sdif,              &! Diffuse shortwave radiation (W/m^2)
      Sdir,              &! Direct-beam shortwave radiation (W/m^2)
      Sf,                &! Snowfall rate (kg/m2/s)
      Sf24h,             &! Snowfall 24hr (kg/m2)
      Ta,                &! Air temperature (K)
      Ua,                &! Wind speed (m/s)
      Udir                ! Wind direction

    use PARAMETERS, only : &
      z0sn              ! Snow roughness length (m)

!-----------------------------------------------------------------------
! Model state variables  
!-----------------------------------------------------------------------
    use STATE_VARIABLES, only: &
      Tsrf,              &! Surface skin temperature (K)
      Tsnow,             &! Snow layer temperatures (K)
      Sice,              &! Ice content of snow layers (kg/m^2)
      Sliq,              &! Liquid content of snow layers (kg/m^2)
      Ds,                &! Snow layer thicknesses (m)
      fsnow,             &! Surface skin temperature (K)
      Nsnow,             &! Number of snow layers
      Tsoil,             &! Soil layer temperatures (K)
      albs,              &! Snow albedo
      theta               ! Volumetric moisture content of soil layers
!-----------------------------------------------------------------------
! output to HICAR variables which re decaled in PHYSICS subroutine but not in MODULES 
!-----------------------------------------------------------------------
    use MODULES_interface
           
    use SNOWTRAN3D_interface, only : SNOWTRAN3D_setup, SNOWTRAN3D_salt_start, SNOWTRAN3D_salt, SNOWTRAN3D_salt_end, SNOWTRAN3D_susp_start, SNOWTRAN3D_susp, SNOWTRAN3D_susp_end, SNOWTRAN3D_accum

    implicit none
    
    private
    !!
    public :: FSM_SETUP,FSM_DRIVE,FSM_PHYSICS,FSM_SNOWSLIDE, FSM_SNOWSLIDE_END, FSM_CUMULATE_SD, FSM_SNOWTRAN_SETUP, FSM_SNOWTRAN_SALT_START, FSM_SNOWTRAN_SALT, FSM_SNOWTRAN_SALT_END, FSM_SNOWTRAN_SUSP_START, FSM_SNOWTRAN_SUSP, FSM_SNOWTRAN_SUSP_END, FSM_SNOWTRAN_ACCUM
    
    public :: Nx_HICAR, Ny_HICAR, zH_HICAR, NNsmax_HICAR, NNsoil_HICAR, lat_HICAR,lon_HICAR,terrain_HICAR,dx_HICAR,slope_HICAR,&
              shd_HICAR, SNTRAN, SNSLID, LHN_ON, LFOR_HN, &
              NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,NOSHDTN,NALRADT,&
              NSNTRAN,NSNSLID,NSNOLAY,NHISWET,NCHECKS, DDs_min,DDs_surflay,                           &
              CTILE,rtthresh, LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT

    !!
    public ::            &
      year,              &
      month,             &
      day,               &
      hour,              &
      dt,                &
      LW,                &
      Ps,                &
      Qa,                &
      Rf,                &
      Sdif,              &
      Sdir,              &
      Sf,                &
      Sf24h,             &
      Ta,                &
      Ua,                &
      Udir
    public ::  	&
      Tsrf,              &
      Tsnow,             &
      Sice,              &
      Sliq,              &
      Ds,                &
      fsnow,             &
      Nsnow,             &
      Tsoil,             &
      albs,              &
      theta,             &
      z0sn
    public ::  	&
      firstit,     &
      Esrf_,       &
      Gsoil_,      &
      H_,          &
      LE_,         &
      Melt_,       &
      Rnet_,       &
      Roff_,       &
      snowdepth_,  &
      SWE_,        &  
      KH_,         &  
      meltflux_out_, &
      Sliq_out_,     &
      dSWE_salt_,      &
      dSWE_susp_,      &
      dSWE_subl_,      &
      dSWE_slide_!,     &
!      Qsalt_u,       &
!      Qsalt_v
      
    
    integer :: Nx_HICAR, Ny_HICAR, NNsmax_HICAR, NNsoil_HICAR

    integer :: &
        NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,NOSHDTN,NALRADT,&
        NSNTRAN,NSNSLID,NSNOLAY,NHISWET,NCHECKS
        
    logical :: LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT
        
    logical :: LHN_ON, LFOR_HN

    real, allocatable:: &
      zH_HICAR(:,:),    & ! height of first model level in HICAR
      lat_HICAR(:,:),   & ! lat from HICAR based WGS84
      lon_HICAR(:,:),   & ! lon from HICAR based WGS84
      slope_HICAR(:,:), & ! terrain slope
      shd_HICAR(:,:),   & ! snow holding depth
      terrain_HICAR(:,:)
 
    real :: & !!!!!!!!---->
      DDs_min,DDs_surflay, & ! minimum snow layer thickness and thickness of surface layer
      rtthresh,            & ! threshold for tile
      dx_HICAR            !   The grid spacing of the high-resolution data in HICAR
      
    character(len=20) :: CTILE
 contains
    
    !!!!!!!!!
    subroutine FSM_SETUP()
        implicit none
        call SETUP_interface()
    end subroutine FSM_SETUP
    !!
    subroutine FSM_DRIVE()
        implicit none
        call DRIVE_interface()
    end subroutine FSM_DRIVE    
    !!
     subroutine FSM_PHYSICS()
        implicit none  
        call PHYSICS_interface()
    end subroutine FSM_PHYSICS  
    
    subroutine FSM_CUMULATE_SD()
        implicit none
        
        call CUMULATE_SD_interface()
        
    end subroutine FSM_CUMULATE_SD


    subroutine FSM_SNOWSLIDE(snowdepth0,Sice0,snowdepth0_buff,Sice0_buff,aval,first_it,dSWE_slide)
        implicit none
        
        real, intent(inout) :: &
          snowdepth0(Nx_HICAR,Ny_HICAR), &
          Sice0(Nx_HICAR,Ny_HICAR),      &
          dSWE_slide(Nx_HICAR,Ny_HICAR)
        real, intent(out) :: &
          snowdepth0_buff(Nx_HICAR,Ny_HICAR), &
          Sice0_buff(Nx_HICAR,Ny_HICAR)

        logical, intent(inout) :: &
          first_it,               &
          aval(Nx_HICAR,Ny_HICAR)

        if (SNSLID==1) then
        
            !Run over all cells, not removing snow from image border
            call SNOWSLIDE_interface(snowdepth0,Sice0,dSWE_slide,first_it,.False.,aval)

            !Save Sice0 and snowdepth0 to output buffers to hand to HICAR for exchange
            snowdepth0_buff = snowdepth0
            Sice0_buff = Sice0
            first_it = .False.
            !Run over image border cells to adjust snowdepth0 and sice0
            call SNOWSLIDE_interface(snowdepth0,Sice0,dSWE_slide,first_it,.True.,aval)
        endif
        
    end subroutine FSM_SNOWSLIDE
    
    subroutine FSM_SNOWSLIDE_END(snowdepth0,Sice0)
        implicit none
        
        real, intent(in) :: &
          snowdepth0(Nx_HICAR,Ny_HICAR), &
          Sice0(Nx_HICAR,Ny_HICAR)
          
        ! Accumulation of new snow, calculation of snow cover fraction and relayering
        call SNOW_LAYERING(snowdepth0,Sice0)

        !call CUMULATE_SD_interface()

    end subroutine FSM_SNOWSLIDE_END


    subroutine FSM_SNOWTRAN_SETUP()
        implicit none

        if (SNTRAN==1) call SNOWTRAN3D_setup()
        
    end subroutine FSM_SNOWTRAN_SETUP

    subroutine FSM_SNOWTRAN_SALT_START(Qs_u, Qs_v)
        implicit none
        real, intent(inout) :: &
          Qs_u(Nx_HICAR,Ny_HICAR), &
          Qs_v(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            call SNOWTRAN3D_salt_start(Qs_u, Qs_v)    
        endif 
        
    end subroutine FSM_SNOWTRAN_SALT_START


    subroutine FSM_SNOWTRAN_SALT(Qs_u, Qs_v)
        implicit none
        real, intent(inout) :: &
          Qs_u(Nx_HICAR,Ny_HICAR), &
          Qs_v(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            call SNOWTRAN3D_salt(Qs_u, Qs_v)    
        endif 
        
    end subroutine FSM_SNOWTRAN_SALT
    
    subroutine FSM_SNOWTRAN_SALT_END(Qs_u, Qs_v)
        implicit none
        real, intent(inout) :: &
          Qs_u(Nx_HICAR,Ny_HICAR), &
          Qs_v(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            call SNOWTRAN3D_salt_end(Qs_u, Qs_v)    
        endif 
        
    end subroutine FSM_SNOWTRAN_SALT_END

    subroutine FSM_SNOWTRAN_SUSP_START(Qs_u, Qs_v)
        implicit none
        real, intent(inout) :: &
          Qs_u(Nx_HICAR,Ny_HICAR), &
          Qs_v(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            call SNOWTRAN3D_susp_start(Qs_u, Qs_v)    
        endif 
        
    end subroutine FSM_SNOWTRAN_SUSP_START

    subroutine FSM_SNOWTRAN_SUSP(Qs_u, Qs_v)
        implicit none
        real, intent(inout) :: &
          Qs_u(Nx_HICAR,Ny_HICAR), &
          Qs_v(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            call SNOWTRAN3D_susp(Qs_u, Qs_v)    
        endif 
        
    end subroutine FSM_SNOWTRAN_SUSP

    subroutine FSM_SNOWTRAN_SUSP_END()
        implicit none

        if (SNTRAN==1) then
            call SNOWTRAN3D_susp_end()    
        endif 
        
    end subroutine FSM_SNOWTRAN_SUSP_END

    
    subroutine FSM_SNOWTRAN_ACCUM()
        implicit none
        real :: &
          dSWE_salt(Nx_HICAR,Ny_HICAR), &
          dSWE_susp(Nx_HICAR,Ny_HICAR), &
          dSWE_subl(Nx_HICAR,Ny_HICAR)

        if (SNTRAN==1) then
            dSWE_salt(:,:) = 0.
            dSWE_susp(:,:) = 0.
            dSWE_subl(:,:) = 0.

            call SNOWTRAN3D_accum(dSWE_salt,dSWE_susp,dSWE_subl)
        
            call CUMULATE_SNOWTRAN3D_interface(dSWE_salt,dSWE_susp,dSWE_subl)
        endif
        
    end subroutine FSM_SNOWTRAN_ACCUM

end module FSM_interface
