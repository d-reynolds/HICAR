
!!
!!----------------------------------------------------------
module module_sf_FSMdrv
    use module_sf_noahdrv,   only : lsm_noah, lsm_noah_init
    use module_water_simple, only : water_simple
    use mod_atm_utilities,   only : sat_mr
    use time_object,         only : Time_type
    use mod_wrf_constants,   only : piconst, XLS
    use icar_constants
    use iso_fortran_env,       only : output_unit
    use options_interface,   only : options_t
    use variable_interface,  only : variable_t
    use domain_interface,    only : domain_t
    use io_routines,         only : io_write, io_read, io_add_attribute
    use FSM_interface , only:  FSM_SETUP,FSM_DRIVE,FSM_PHYSICS, FSM_SNOWSLIDE, FSM_SNOWSLIDE_END, FSM_CUMULATE_SD, FSM_SNOWTRAN_SETUP, FSM_SNOWTRAN_SALT_START, FSM_SNOWTRAN_SALT, FSM_SNOWTRAN_SALT_END, FSM_SNOWTRAN_SUSP_START, FSM_SNOWTRAN_SUSP, FSM_SNOWTRAN_SUSP_END, FSM_SNOWTRAN_ACCUM
    use FSM_interface , only:  Nx_HICAR, Ny_HICAR, NNsmax_HICAR, NNsoil_HICAR, lat_HICAR, &
                               lon_HICAR,terrain_HICAR,dx_HICAR,slope_HICAR,shd_HICAR, LHN_ON, LFOR_HN, &
                               NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,NOSHDTN,NALRADT,&
                               NSNTRAN,NSNSLID,NSNOLAY,NCHECKS, LZ0PERT, LWCPERT, LFSPERT, LALPERT, LSLPERT, CTILE, rtthresh, DDs_min, DDs_surflay
    use FSM_interface, only: &
      year,          &
      month,         &
      day,           &
      hour,          &
      dt,            &
      LW,            &
      Ps,            &
      Qa,            &
      Rf,            &
      Sdif,          &
      Sdir,          &
      Sf,            &
      Sf24h,         &
      Ta,            &
      Ua,            &
      Udir
    use FSM_interface, only: &
      Esrf_,         &
      Gsoil_,        &
      H_,            &
      LE_,           &
      Melt_,         &
      Rnet_,         &
      Roff_,         &
      snowdepth_,    &
      SWE_,          &  
      KH_,           &  
      meltflux_out_, &
      Sliq_out_,     &
      dSWE_salt_,      &
      dSWE_susp_,      &
      dSWE_subl_,      &
      dSWE_slide_!,     &
      !Qs_u,       &
      !Qs_v

    use FSM_interface, only: &
      firstit,       &
      Tsrf,          &
      Tsnow,         &
      Sice,          &
      Sliq,          &
      Ds,            &
      fsnow,         &
      histowet,      &
      Nsnow,         &
      Tsoil,         &
      albs,          &
      theta,         &
      z0sn
    use FSM_interface, only: SNTRAN, SNSLID
    
    implicit none

    private
    public :: sm_FSM_init,sm_FSM
   
    type(variable_t) :: Qs_u_var, Qs_v_var, SD_0_var, Sice_0_var  
    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions
    integer :: hour_old

    type(Time_type) :: last_output
    real            :: last_snowslide
   
    real, allocatable :: &
        z0_bare(:,:),    &       !bare ground z0 before getting covered in snow
        Sf24_tracker(:,:,:), &    !24h snow fall tracker
        ! snowfall_sum(:,:),       & !aggregated per output interval
        ! rainfall_sum(:,:),       & !aggregated per output interval
        Roff_sum(:,:),           & !aggregated per output interval
        meltflux_out_sum(:,:)      !aggregated per output interval

contains
 
    subroutine sm_FSM_init(domain,options, context_chng)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        logical, optional, intent(in) :: context_chng
        integer :: i,j,k, hj, hi, i_s, i_e, j_s, j_e
        logical :: context_change
        !!

        if (present(context_chng)) then
            context_change = context_chng
        else
            context_change = .False.
        endif
        
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        
        !if (SNTRAN+SNSLID > 0) then
        its = domain%grid%its-1
        ite = domain%grid%ite+1
        jts = domain%grid%jts-1
        jte = domain%grid%jte+1
        !else
        !    its = domain%grid%its
        !    ite = domain%grid%ite
        !    jts = domain%grid%jts
        !    jte = domain%grid%jte
        !endif
        !!       
        
        Ny_HICAR=ite-its+1
        Nx_HICAR=jte-jts+1
        
        if (options%restart%restart) then
            last_output = options%restart%restart_time
        else if (context_change) then
            last_output = domain%next_output - domain%output_dt
        else
            last_output = options%general%start_time
        endif
        
        last_snowslide = 4000

        !!
        if (allocated(lat_HICAR)) deallocate(lat_HICAR)
        if (allocated(lon_HICAR)) deallocate(lon_HICAR)
        if (allocated(terrain_HICAR)) deallocate(terrain_HICAR)
        if (allocated(slope_HICAR)) deallocate(slope_HICAR)
        if (allocated(shd_HICAR)) deallocate(shd_HICAR)
        if (allocated(z0_bare)) deallocate(z0_bare)

        allocate(lat_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(lon_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(terrain_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(slope_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(shd_HICAR(Nx_HICAR,Ny_HICAR))
        allocate(z0_bare(domain%grid%its:domain%grid%ite,domain%grid%jts:domain%grid%jte))

        lat_HICAR=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d(its:ite,jts:jte))
        lon_HICAR=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%longitude)%v)%data_2d(its:ite,jts:jte))
        terrain_HICAR=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d(its:ite,jts:jte))
        if (domain%var_indx(kVARS%slope_angle)%v > 0) then
            slope_HICAR=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%slope_angle)%v)%data_2d(its:ite,jts:jte))*180.0/piconst !Convert from radians to degrees
        else
            slope_HICAR = 0.
        endif
        if (domain%var_indx(kVARS%shd)%v > 0) then
            shd_HICAR=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%shd)%v)%data_2d(its:ite,jts:jte))
        else
            shd_HICAR = 10000.
        endif
        
        z0_bare = 0.01
        dx_HICAR=domain%dx
        NNsoil_HICAR=options%lsm%num_soil_layers
        NNsmax_HICAR=options%sm%fsm_nsnow_max
        DDs_min=options%sm%fsm_ds_min
        DDs_surflay=options%sm%fsm_ds_surflay

        !copy over namelist options to FSM2trans
        NALBEDO=options%sm%fsm_albedo
        NCANMOD=options%sm%fsm_canmod
        NCONDCT=options%sm%fsm_condct
        NDENSTY=options%sm%fsm_densty
        NEXCHNG=options%sm%fsm_exchng
        NHYDROL=options%sm%fsm_hydrol
        NSNFRAC=options%sm%fsm_snfrac
        NRADSBG=options%sm%fsm_radsbg
        NZOFFST=options%sm%fsm_zoffst
        NOSHDTN=options%sm%fsm_oshdtn
        NALRADT=options%sm%fsm_alradt
        NSNTRAN=options%sm%fsm_sntran
        NSNSLID=options%sm%fsm_snslid
        NSNOLAY=options%sm%fsm_snolay
        NCHECKS=options%sm%fsm_checks

        LZ0PERT=options%sm%fsm_z0pert
        LWCPERT=options%sm%fsm_wcpert
        LFSPERT=options%sm%fsm_fspert
        LALPERT=options%sm%fsm_alpert
        LSLPERT=options%sm%fsm_slpert

        ! Hard coding these tile values for now
        ! not sure how they are supposed to be used?
        CTILE='open'
        rtthresh=0.1

        LHN_ON=options%sm%fsm_hn_on
        LFOR_HN=options%sm%fsm_for_hn

        !!
        if (allocated(Esrf_)) deallocate(Esrf_)
        if (allocated(Gsoil_)) deallocate(Gsoil_)
        if (allocated(H_)) deallocate(H_)
        if (allocated(LE_)) deallocate(LE_)
        if (allocated(Melt_)) deallocate(Melt_)
        if (allocated(Rnet_)) deallocate(Rnet_)
        if (allocated(Roff_)) deallocate(Roff_)
        if (allocated(snowdepth_)) deallocate(snowdepth_)
        if (allocated(SWE_)) deallocate(SWE_)
        if (allocated(KH_)) deallocate(KH_)
        if (allocated(meltflux_out_)) deallocate(meltflux_out_)
        if (allocated(Sliq_out_)) deallocate(Sliq_out_)
        if (allocated(dSWE_salt_)) deallocate(dSWE_salt_)
        if (allocated(dSWE_susp_)) deallocate(dSWE_susp_)
        if (allocated(dSWE_subl_)) deallocate(dSWE_subl_)
        if (allocated(dSWE_slide_)) deallocate(dSWE_slide_)
        if (allocated(Sf24_tracker)) deallocate(Sf24_tracker)
        if (allocated(Roff_sum)) deallocate(Roff_sum)
        if (allocated(meltflux_out_sum)) deallocate(meltflux_out_sum)

        allocate(Esrf_(Nx_HICAR,Ny_HICAR)); Esrf_=0.
        allocate(Gsoil_(Nx_HICAR,Ny_HICAR));Gsoil_=0.
        allocate(H_(Nx_HICAR,Ny_HICAR)); H_=0.
        allocate(LE_(Nx_HICAR,Ny_HICAR)); LE_=0.
        allocate(Melt_(Nx_HICAR,Ny_HICAR)); Melt_=0.
        allocate(Rnet_(Nx_HICAR,Ny_HICAR)); Rnet_=0.
        allocate(Roff_(Nx_HICAR,Ny_HICAR)); Roff_=0.
        allocate(snowdepth_(Nx_HICAR,Ny_HICAR)); snowdepth_=0.
        allocate(SWE_(Nx_HICAR,Ny_HICAR)); SWE_=0.0
        allocate(KH_(Nx_HICAR,Ny_HICAR)); KH_=0.0
        allocate(meltflux_out_(Nx_HICAR,Ny_HICAR)); meltflux_out_=0.
        allocate(Sliq_out_(Nx_HICAR,Ny_HICAR)); Sliq_out_=0.
        allocate(dSWE_salt_(Nx_HICAR,Ny_HICAR)); dSWE_salt_=0.
        allocate(dSWE_susp_(Nx_HICAR,Ny_HICAR)); dSWE_susp_=0.
        allocate(dSWE_subl_(Nx_HICAR,Ny_HICAR)); dSWE_subl_=0.
        allocate(dSWE_slide_(Nx_HICAR,Ny_HICAR)); dSWE_slide_=0.
        allocate(Sf24_tracker(24,Nx_HICAR,Ny_HICAR));Sf24_tracker=0.
        allocate(Roff_sum(Nx_HICAR,Ny_HICAR)); Roff_sum=0.
        allocate(meltflux_out_sum(Nx_HICAR,Ny_HICAR)); meltflux_out_sum=0.   
        !!

        call FSM_SETUP()

        if (SNTRAN+SNSLID > 0) then
            call Qs_u_var%initialize(domain%grid2d)
            call Qs_v_var%initialize(domain%grid2d)
            call SD_0_var%initialize(domain%grid2d)
            call Sice_0_var%initialize(domain%grid2d)
        endif

        !!        
        !! MJ added this block to read in while we use restart file:
        if (options%restart%restart .or. context_change) then

            !! giving feedback to HICAR
            Tsrf = TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d(its:ite,jts:jte))
            if (options%lsm%monthly_albedo) then
                albs = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(its:ite, domain%sim_time%month, jts:jte))
            else
                albs = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(its:ite, 1, jts:jte))
            endif
            
            fsnow = TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(its:ite,jts:jte))
            Nsnow = TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(its:ite,jts:jte))                        
            !!
            do i=1,NNsmax_HICAR
                Tsnow(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(its:ite,i,jts:jte))
                Sice(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d(its:ite,i,jts:jte))
                Sliq(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d(its:ite,i,jts:jte))
                Ds(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(its:ite,i,jts:jte))
            enddo
            do i=1,kSOIL_GRID_Z
                Tsoil(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d(its:ite,i,jts:jte))
                theta(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(its:ite,i,jts:jte))
            enddo

            call exch_FSM_state_vars(domain,corners_in=.True.)
        !If this is not a restart or a context change (i.e. pure intialization), see if we were 
        ! initilaized with any snow already present...
        else

            j_s = 2
            i_s = 2
            j_e = Nx_HICAR-1
            i_e = Ny_HICAR-1
    
            do j=j_s,j_e
                do i=i_s,i_e
                    hj = j-j_s+domain%jts
                    hi = i-i_s+domain%its

                    ! The state variables in FSM to set are Ds (array of snow depth by layer)
                    ! Sice (array of ice content by layer), Sliq (array of liquid water content by layer),
                    ! Nsnow (number of snow layers), fsnow (fraction of snow cover)
                    if ( (domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(hi,hj) > 0) .and. .not.(options%physics%watersurface==kWATER_SIMPLE .and. domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(hi,hj)==kLC_WATER)) then

                        !find out how many even layers of snow we have
                        do k = 1,(NNsmax_HICAR-1)
                            if (domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(hi,hj)/(k+1) < DDs_min) exit
                        enddo

                        Nsnow(j,i)     = 1 
                        fsnow(j,i)     = 0.95
                        domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(hi,hj) = Nsnow(j,i)
                        domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(hi,hj) = fsnow(j,i)

                        Ds(1:k,j,i)      = domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(hi,hj)/k  
                        Sice(1:k,j,i)    = domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(hi,hj)/k 
                        Sliq(1:k,j,i)    = 0.0
                        Tsnow(1:k,j,i)   = 273.15
                        Tsrf(j,i)      = 273.15
                        histowet(1:k,j,i) = 0.0

                        domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(hi,1:k,hj) = Tsnow(1:k,j,i)
                        domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d(hi,1:k,hj) = Sliq(1:k,j,i) 
                        domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d(hi,1:k,hj) = Sice(1:k,j,i) 
                        domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(hi,1:k,hj) = Ds(1:k,j,i)
            
                    endif
                enddo
            enddo
            
            call exch_FSM_state_vars(domain,corners_in=.True.)

        endif

        !Test if restart was not succesful, or if we were not passed FSM restart...
        !do j = 1, Ny_HICAR
        !    do i = 1, Nx_HICAR
        !        if (.not.(fsnow(i,j)*SUM(Ds(:,i,j))==domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(j-1+domain%its,i-1+domain%jts))) then
        !            albs(i,j) = 0.75
        !            fsnow(i,j) = 1
        !            Nsnow (i,j)= min(6,options%lsm%fsm_nsnow_max)
        !            Tsnow(:,i,j) = 268
        !            Tsoil(:,i,j) = 273
        !            theta(:,i,j) = 0.2
        !            Sliq(:,i,j) = 0.0
        !            Ds(:,i,j) = domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(j-1+domain%its,i-1+domain%jts)/Nsnow(i,j)
        !            Sice(:,i,j) = domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(j-1+domain%its,i-1+domain%jts)/Nsnow(i,j)
        !        endif
        !    end do
        !end do
        !SYNC ALL
    end subroutine sm_FSM_init

    subroutine sm_FSM(domain,options,lsm_dt,current_rain,current_snow,windspd)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        !!       
        real, intent(in) :: &
            lsm_dt          ! Timestep (s) from HICAR for lsm
        real, dimension(ims:ime,jms:jme), intent(in) :: &   !! Note that: when the input is an array with whatever indexing either using (:,:) or (N,M), it always starts from index 1 in all dims
            current_rain,                 &! rainfall in kg m-2
            current_snow,                 &! snowfall in kg m-2
            windspd                        ! Wind speed (m/s)

        integer :: i,j,k, hj, hi, i_s, i_e, j_s, j_e
        real :: Delta_t
        real, dimension(its:ite,jts:jte) :: SWE_pre
        real, dimension(Nx_HICAR,Ny_HICAR) :: SD_0, Sice_0, SD_0_buff, Sice_0_buff, Qs_u, Qs_v
        logical, dimension(Nx_HICAR,Ny_HICAR) :: aval
        logical :: first_SLIDE, do_snowslide
        character(len=1024) :: filename
        
        !if (SNTRAN+SNSLID > 0) then
        do_snowslide = .False.
        if (SNSLID > 0) then
            do_snowslide = (last_snowslide > 3600.) !Hard code to just be called every hour...
        endif

        j_s = 2
        i_s = 2
        j_e = Nx_HICAR-1
        i_e = Ny_HICAR-1
        !else
        !    j_s = 1
        !    i_s = 1
        !    j_e = Ny_HICAR
        !    i_e = Nx_HICAR
        !endif
        
        !!
        !! giving the required input from HICAR to FSM
        year=real(domain%sim_time%year)
        month=real(domain%sim_time%month)
        day=real(domain%sim_time%day)
        hour=real(domain%sim_time%hour)
        dt=lsm_dt
        
        call exch_FSM_state_vars(domain,corners_in=.True.)

        LW=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%longwave)%v)%data_2d(its:ite,jts:jte))
        Ps=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%surface_pressure)%v)%data_2d(its:ite,jts:jte))
        Rf=TRANSPOSE(current_rain(its:ite,jts:jte))
        Sdir=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%shortwave_direct)%v)%data_2d(its:ite,jts:jte))  !Sdir=domain%vars_2d(domain%var_indx(kVARS%shortwave)%v)%data_2d(its:ite,jts:jte)
        Sdif=TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%shortwave_diffuse)%v)%data_2d(its:ite,jts:jte)) !Sdif=0.0
        Sf=TRANSPOSE(current_snow(its:ite,jts:jte))

        ! if (hour /= hour_old) then
        !     Sf24_tracker(INT(hour),:,:) = Sf

        !     hour_old = hour
        ! else
        !     Sf24_tracker(INT(hour),:,:) = Sf24_tracker(INT(hour),:,:) + Sf
        ! endif
        Sf24h=0.0 !sum(Sf24_tracker,dim=1)
        
        ! This should always be 2m temperature and humidity
        Ta= TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%temperature_2m)%v)%data_2d(its:ite,jts:jte))
        Qa= TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%humidity_2m)%v)%data_2d(its:ite,jts:jte))

        ! This should always be 10m wind speed
        Ua=TRANSPOSE(windspd(its:ite,jts:jte))
        
        !CAUTION -- Udir hardcoded to use domain%u/v instead of v_10m/u_10m. This is done since we need to access 
        ! these variables, one index out from the tile indices (For FSM). However, the LSM only updates the roughness length
        ! within the tile indices, which is needed to compute 10m wind speeds. As a cheep work around, and since we are 
        ! only doing simulations with first level thickness =20m, they are hard coded to lowest level of 3d arrays here
        
        Udir = ATAN2(TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(its:ite,domain%grid%kms,jts:jte)), &
                     TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(its:ite,domain%grid%kms,jts:jte))) 
        Udir = 90 - (Udir * 180/piconst + 180)
        where(Udir<0) Udir=Udir+360
        
        if ((domain%sim_time%seconds()   >=    last_output%seconds())) then
            !If we are the first call since the last output, reset the per-output counters
            domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d = 0.
            domain%vars_2d(domain%var_indx(kVARS%dSWE_salt)%v)%data_2d  = 0.
            domain%vars_2d(domain%var_indx(kVARS%dSWE_susp)%v)%data_2d  = 0.
            domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d  = 0.
            last_output = domain%next_output
        endif
        
        SWE_pre = domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(its:ite,jts:jte)
                !!  
        !! FSM processing      
        call FSM_DRIVE()
        call FSM_PHYSICS()
                        
        !Call Snowtran here -- must be done here and not in FSM since we need control over the parallelization of the routine
        if (SNTRAN > 0) then        
            ! call exch_FSM_state_vars(domain)
            
            call FSM_SNOWTRAN_SETUP()
            
            !--------------------------SALTATION--------------------------
            call FSM_SNOWTRAN_SALT_START(Qs_u,Qs_v)
            
            !First guess for fluxes
            call FSM_SNOWTRAN_SALT(Qs_u,Qs_v)
            
            !Exchange fluxes between processes
            call exch_SNTRAN_Qs(domain,Qs_u,Qs_v)
            !Recalculate fluxes with intermediate values from neighbors
            call FSM_SNOWTRAN_SALT(Qs_u,Qs_v)
            
            !Exchange fluxes between processes
            call exch_SNTRAN_Qs(domain,Qs_u,Qs_v)
            !Recalculate fluxes with intermediate values from neighbors
            call FSM_SNOWTRAN_SALT(Qs_u,Qs_v)
            
            call FSM_SNOWTRAN_SALT_END(Qs_u,Qs_v)
            !------------------------END SALTATION------------------------

            !--------------------------SUSPENSION--------------------------
            call FSM_SNOWTRAN_SUSP_START(Qs_u,Qs_v)
            
            !First guess for fluxes
            call FSM_SNOWTRAN_SUSP(Qs_u,Qs_v)
            
            !Exchange fluxes between processes
            call exch_SNTRAN_Qs(domain,Qs_u,Qs_v)
            !Recalculate fluxes with intermediate values from neighbors
            call FSM_SNOWTRAN_SUSP(Qs_u,Qs_v)
            
            !Exchange fluxes between processes
            call exch_SNTRAN_Qs(domain,Qs_u,Qs_v)
            !Recalculate fluxes with intermediate values from neighbors
            call FSM_SNOWTRAN_SUSP(Qs_u,Qs_v)

            call FSM_SNOWTRAN_SUSP_END()
            !------------------------END SUSPENSION------------------------

        endif

        
        !Call Snowslide here -- must be done here and not in FSM since we need control over the parallelization of the routine
        if (do_snowslide .and. SNSLID > 0) then

            last_snowslide = last_snowslide + lsm_dt

            if (do_snowslide) then
                SD_0 = 0.0
                Sice_0 = 0.0
                SD_0_buff = 0.0
                Sice_0_buff = 0.0

                dSWE_slide_ = 0.0
                first_SLIDE = .True.
                aval = .False.
                !Snowslide needs the corner snow depth information from corner neighbor processes
                ! call exch_FSM_state_vars(domain,corners_in=.True.)
                
                do i=1,10
                    call FSM_SNOWSLIDE(SD_0,Sice_0,SD_0_buff,Sice_0_buff,aval,first_SLIDE,dSWE_slide_)
                    
                    ! Copy interior buffer for exchange
                    call exch_SLIDE_buffers(domain,SD_0_buff,Sice_0_buff)
                    
                    !Now, where the buffer is positive for halo cells, record that an avalanche was passed here
                    aval = .False.
                    where(SD_0_buff > SD_0) aval=.True.
                    aval(j_s:j_e,i_s:i_e) = .False.
                    
                    !Must accumulate slide changes here, since we will loop over calls to snowslide
                    domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = &
                            domain%vars_2d(domain%var_indx(kVARS%dSWE_slide)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte) + TRANSPOSE(dSWE_slide_(j_s:j_e,i_s:i_e))

                    ! trade over buffer to halo cells, in case avalanching snow was passed on previous iteration.
                    SD_0(1,:) = SD_0_buff(1,:)
                    SD_0(Nx_HICAR,:) = SD_0_buff(Nx_HICAR,:)
                    SD_0(:,1) = SD_0_buff(:,1)
                    SD_0(:,Ny_HICAR) = SD_0_buff(:,Ny_HICAR)

                    Sice_0(1,:) = Sice_0_buff(1,:)
                    Sice_0(Nx_HICAR,:) = Sice_0_buff(Nx_HICAR,:)
                    Sice_0(:,1) = Sice_0_buff(:,1)
                    Sice_0(:,Ny_HICAR) = Sice_0_buff(:,Ny_HICAR)


                enddo
                call FSM_SNOWSLIDE_END(SD_0,Sice_0)
                last_snowslide = 0
            endif
        endif
        
        !Accumulate all the SNTRAN fluxes calculated above
        ! Do this here so that the grid and halo are still consistent for SNSLID if it is called
        if (SNTRAN > 0) call FSM_SNOWTRAN_ACCUM()

        !The End.
        call FSM_CUMULATE_SD()
        !! giving feedback to HICAR -- should only be done for snow-covered cells, or cells which were just snowed on
        do j=j_s,j_e
            do i=i_s,i_e
                hj = j-j_s+domain%jts
                hi = i-i_s+domain%its
                if ( SWE_(j,i)+SWE_pre(hi,hj) > 0 .and. .not.(options%physics%watersurface==kWATER_SIMPLE .and. domain%vars_2d(domain%var_indx(kVARS%land_mask)%v)%data_2di(hi,hj)==kLC_WATER)) then
                    !If we are covering pixel for the first time, save current (bare) roughness length for later
                    if (domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(hi,hj)==0) then
                        z0_bare(hi,hj) = domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d(hi,hj)
                    endif    
                    domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d(hi,hj)=H_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d(hi,hj)=LE_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(hi,hj)=SWE_(j,i)
                    if (options%lsm%monthly_albedo) then
                        domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(hi, domain%sim_time%month, hj) = albs(j,i)
                    else
                        domain%vars_3d(domain%var_indx(kVARS%albedo)%v)%data_3d(hi, 1, hj) = albs(j,i)
                    endif
                    !
                    domain%vars_2d(domain%var_indx(kVARS%skin_temperature)%v)%data_2d(hi,hj)=Tsrf(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%snow_height)%v)%data_2d(hi,hj)=snowdepth_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%Sliq_out)%v)%data_2d(hi,hj)=Sliq_out_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(hi,hj)=fsnow(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(hi,hj)=Nsnow(j,i)
                    !
                    domain%vars_2d(domain%var_indx(kVARS%qfx)%v)%data_2d(hi,hj)=Esrf_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%chs)%v)%data_2d(hi,hj)=KH_(j,i)
                    domain%vars_2d(domain%var_indx(kVARS%roughness_z0)%v)%data_2d(hi,hj) = z0sn*fsnow(j,i)+(1-fsnow(j,i))*z0_bare(hi,hj)
                    
                    do k=1,kSNOW_GRID_Z
                        domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(hi,k,hj) = Tsnow(k,j,i)
                        domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d(hi,k,hj) = Sice(k,j,i)
                        domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d(hi,k,hj) = Sliq(k,j,i)
                        domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(hi,k,hj) = Ds(k,j,i)
                    enddo
                    do k=1,kSOIL_GRID_Z
                        domain%vars_3d(domain%var_indx(kVARS%soil_temperature)%v)%data_3d(hi,k,hj) = Tsoil(k,j,i)
                        domain%vars_3d(domain%var_indx(kVARS%soil_water_content)%v)%data_3d(hi,k,hj)=theta(k,j,i)
                    enddo
                
                    if (SNTRAN>0) then
                        ! Convert to rate 1/s
                        domain%vars_2d(domain%var_indx(kVARS%dSWE_salt)%v)%data_2d(hi,hj)=domain%vars_2d(domain%var_indx(kVARS%dSWE_salt)%v)%data_2d(hi,hj) + dSWE_salt_(j,i)
                        domain%vars_2d(domain%var_indx(kVARS%dSWE_susp)%v)%data_2d(hi,hj)= domain%vars_2d(domain%var_indx(kVARS%dSWE_susp)%v)%data_2d(hi,hj) + dSWE_susp_(j,i)
                        domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d(hi,hj)= domain%vars_2d(domain%var_indx(kVARS%dSWE_subl)%v)%data_2d(hi,hj) + dSWE_subl_(j,i)
                        !Add sublimated snow to latent heat flux. 
                        !Sometimes FSM returns NaN values for blowing snow sublimation, so mask those out here
                        if (abs(dSWE_subl_(j,i))>1) dSWE_subl_(j,i) = 0.0
                        domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d(hi,hj)   = domain%vars_2d(domain%var_indx(kVARS%latent_heat)%v)%data_2d(hi,hj)   + (-dSWE_subl_(j,i))*XLS/dt 
                        domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d(hi,hj) = domain%vars_2d(domain%var_indx(kVARS%sensible_heat)%v)%data_2d(hi,hj) + ( dSWE_subl_(j,i))*XLS/dt 
                    endif
                endif
            enddo
        enddo
        
        !If we just melted out the snow in this step, set this variable. This is to ensure good hand-off from snow-covered to not with NoahMP
        where(SWE_pre(domain%its:domain%ite,domain%jts:domain%jte) > 0 .and. &
                domain%vars_2d(domain%var_indx(kVARS%snow_water_equivalent)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte)==0) 
            domain%vars_2d(domain%var_indx(kVARS%ground_surf_temperature)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte)=273.15
        end where
        
        ! Let FSM know that it has done the first iteration 
        if (firstit==1) firstit=0
        
        !These are FSM2 Diagnostic/output vars, so we can update them everywhere
        !domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(its:ite,jts:jte)=Nsnow
        !!
        !snowfall_sum=snowfall_sum+Sf
        !rainfall_sum=rainfall_sum+Rf
        Roff_sum=Roff_sum+Roff_
        meltflux_out_sum=meltflux_out_sum+meltflux_out_
        

        Delta_t=mod(domain%sim_time%seconds(),options%output%output_dt%seconds())
        if ( abs(options%output%output_dt%seconds()-(Delta_t+dt)) <= 1.e-3 ) then
        !    if (STD_OUT_PE) write(*,*) "resetting/aggregating vars e.g. runoff during t-1->t"!, Delta_t,Delta_t+dt
        !    !!
        !    domain%vars_2d(domain%var_indx(kVARS%rainfall_tstep)%v)%data_2d(its:ite,jts:jte)=rainfall_sum
        !    domain%vars_2d(domain%var_indx(kVARS%snowfall_tstep)%v)%data_2d(its:ite,jts:jte)=snowfall_sum
           domain%vars_2d(domain%var_indx(kVARS%runoff_tstep)%v)%data_2d(its:ite,jts:jte)=TRANSPOSE(Roff_sum)
           domain%vars_2d(domain%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d(its:ite,jts:jte)=TRANSPOSE(meltflux_out_sum)
        !    !! reseting the container to zero for next output interval
        !    rainfall_sum = 0.
        !    snowfall_sum = 0.                
           Roff_sum = 0.
           meltflux_out_sum = 0.
        endif
        
        !!
        !SYNC ALL  
    end subroutine sm_FSM
    
    
    subroutine exch_FSM_state_vars(domain,corners_in)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        logical, optional, intent(in) :: corners_in
        
        integer :: i
        logical :: corners
        
        corners=.False.
        if (present(corners_in)) corners=corners_in
        
        domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(fsnow(2:Nx_HICAR-1,2:Ny_HICAR-1))
        domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(Nsnow(2:Nx_HICAR-1,2:Ny_HICAR-1))                        
        !!
        do i=1,NNsmax_HICAR
            domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Tsnow(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Sice(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Sliq(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
            domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(domain%its:domain%ite,i,domain%jts:domain%jte) = TRANSPOSE(Ds(i,2:Nx_HICAR-1,2:Ny_HICAR-1))
        enddo

        call domain%halo%exch_var(domain%vars_2d(domain%var_indx(kVARS%fsnow)%v))
        call domain%halo%exch_var(domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v))
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v))
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Sice)%v))
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Sliq)%v))
        call domain%halo%exch_var(domain%vars_3d(domain%var_indx(kVARS%Ds)%v))

        !call domain%halo%batch_exch(domain%exch_vars, domain%adv_vars, two_d=.True.)
        !call domain%halo%batch_exch(domain%exch_vars, domain%adv_vars, two_d=.False.,exch_var_only=.True.)      

        if (corners) call domain%halo%exch_var(domain%vars_2d(domain%var_indx(kVARS%Ds)%v),corners=corners)
        if (corners) call domain%halo%exch_var(domain%vars_2d(domain%var_indx(kVARS%fsnow)%v),corners=corners)

        fsnow = TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%fsnow)%v)%data_2d(its:ite,jts:jte))
        Nsnow = TRANSPOSE(domain%vars_2d(domain%var_indx(kVARS%Nsnow)%v)%data_2d(its:ite,jts:jte))                        
        !!
        do i=1,NNsmax_HICAR
            Tsnow(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%snow_temperature)%v)%data_3d(its:ite,i,jts:jte))
            Sice(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Sice)%v)%data_3d(its:ite,i,jts:jte))
            Sliq(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Sliq)%v)%data_3d(its:ite,i,jts:jte))
            Ds(i,:,:) = TRANSPOSE(domain%vars_3d(domain%var_indx(kVARS%Ds)%v)%data_3d(its:ite,i,jts:jte))
        enddo
        
    end subroutine exch_FSM_state_vars
    
    
    subroutine exch_SNTRAN_Qs(domain, Qs_u, Qs_v)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        real, dimension(Nx_HICAR,Ny_HICAR), intent(inout) :: Qs_u, Qs_v

        Qs_u_var%data_2d = 0.0
        Qs_v_var%data_2d = 0.0

        Qs_u_var%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(Qs_u(2:Nx_HICAR-1,2:Ny_HICAR-1))
        Qs_v_var%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(Qs_v(2:Nx_HICAR-1,2:Ny_HICAR-1))
        
        call domain%halo%exch_var(Qs_u_var)
        call domain%halo%exch_var(Qs_v_var)

        Qs_u = TRANSPOSE(Qs_u_var%data_2d(its:ite,jts:jte))
        Qs_v = TRANSPOSE(Qs_v_var%data_2d(its:ite,jts:jte))

    end subroutine exch_SNTRAN_Qs

    subroutine exch_SLIDE_buffers(domain,SD_0,Sice_0)
        implicit none
        
        type(domain_t), intent(inout) :: domain
        
        real, dimension(Nx_HICAR,Ny_HICAR), intent(inout) :: SD_0, Sice_0

        SD_0_var%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(SD_0(2:Nx_HICAR-1,2:Ny_HICAR-1))
        Sice_0_var%data_2d(domain%its:domain%ite,domain%jts:domain%jte) = TRANSPOSE(Sice_0(2:Nx_HICAR-1,2:Ny_HICAR-1))
        
        call domain%halo%exch_var(SD_0_var)
        call domain%halo%exch_var(Sice_0_var)
        call domain%halo%exch_var(SD_0_var, corners=.True.)
        call domain%halo%exch_var(Sice_0_var, corners=.True.)

        SD_0 = TRANSPOSE(SD_0_var%data_2d(its:ite,jts:jte))
        Sice_0 = TRANSPOSE(Sice_0_var%data_2d(its:ite,jts:jte))


    end subroutine exch_SLIDE_buffers


!!
end module module_sf_FSMdrv
