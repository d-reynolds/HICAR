submodule(options_interface) options_implementation

    use icar_constants
    use mod_wrf_constants,          only : piconst
    use io_routines,                only : io_newunit, check_variable_present, check_file_exists
    use time_delta_object,          only : time_delta_t
    use time_object,                only : Time_type
    use string,                     only : str
    use model_tracking,             only : print_model_diffs
    use output_metadata,            only : get_varname, get_varindx
    use namelist_utils,             only : set_nml_var, set_nml_var_default, set_namelist, write_nml_file_end
    implicit none

contains


    !> ----------------------------------------------------------------------------
    !!  Read all namelists from the options file specified on the command line
    !!
    !!  Reads the commandline (or uses default icar_options.nml filename)
    !!  Checks that the version of the options file matches the version of the code
    !!  Reads each namelist successively, all options are stored in supplied options object
    !!
    !! ----------------------------------------------------------------------------

    !> -------------------------------
    !! Initialize an options object
    !!
    !!
    !! -------------------------------
    module subroutine init(this, namelist_file, n_indx, info_only, gen_nml)
        implicit none
        class(options_t),   intent(inout)  :: this
        character(len=*),   intent(in)     :: namelist_file
        integer,            intent(in)     :: n_indx
        logical,            intent(in)     :: info_only, gen_nml

!       reads a series of options from a namelist file and stores them in the
!       options data structure
        integer :: i

        this%nest_indx = n_indx
        if (STD_OUT_PE .and. n_indx==1 .and. .not.(info_only .or. gen_nml)) write(*,*) "  Using options file = ", trim(namelist_file)
        if (gen_nml) call set_namelist(namelist_file)

        call general_namelist(         namelist_file,   this%general, n_indx, info_only=info_only, gen_nml=gen_nml)
        call restart_namelist(         namelist_file,   this, n_indx, info_only=info_only, gen_nml=gen_nml)
        call domain_namelist(          namelist_file,   this%domain, n_indx, info_only=info_only, gen_nml=gen_nml)
        call forcing_namelist(         namelist_file,   this, n_indx, info_only=info_only, gen_nml=gen_nml)
        call physics_namelist(         namelist_file,   this%physics, n_indx, info_only=info_only, gen_nml=gen_nml)
        call time_parameters_namelist( namelist_file,   this%time, n_indx, info_only=info_only, gen_nml=gen_nml)
        call lt_parameters_namelist(   namelist_file,   this%lt, this%general%use_lt_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call mp_parameters_namelist(   namelist_file,   this%mp, this%general%use_mp_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call adv_parameters_namelist(  namelist_file,   this%adv, this%general%use_adv_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call sm_parameters_namelist(   namelist_file,   this%sm, this%general%use_sm_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call lsm_parameters_namelist(  namelist_file,   this%lsm, this%general%use_lsm_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call cu_parameters_namelist(   namelist_file,   this%cu, this%general%use_cu_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call rad_parameters_namelist(  namelist_file,   this%rad, this%general%use_rad_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call pbl_parameters_namelist(  namelist_file,   this%pbl, this%general%use_pbl_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call sfc_parameters_namelist(  namelist_file,   this%sfc, this%general%use_sfc_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call wind_namelist(            namelist_file,   this%wind, this%general%use_wind_options, n_indx, info_only=info_only, gen_nml=gen_nml)
        call output_namelist(          namelist_file,   this%output, n_indx, info_only=info_only, gen_nml=gen_nml)

        ! If this run was just done to output the namelist options, stop now
        if (info_only .or. gen_nml) then
            if (gen_nml) then
                if (STD_OUT_PE) write(*,*) 'Default namelist written to file: ', trim(namelist_file)
                call write_nml_file_end()
            endif
            stop
        endif

        if (this%general%phys_suite /= '') call set_phys_suite(this)
        
        if (this%restart%restart) this%general%start_time = this%restart%restart_time

        call default_var_requests(this)

        if (n_indx == 1) call version_check(this%general)

    end subroutine init

    !> -------------------------------
    !! Check options for consistency among different nests
    !!
    !! -------------------------------
    subroutine inter_nest_options_check(options)
        implicit none
        type(options_t), intent(inout) :: options(:)

        integer :: i, n, child_nest_indx
        integer :: new_restartinterval
        type(Time_type) :: parent_restart_time, child_restart_time, new_child_restart_time
        type(time_delta_t) :: nest_start_offset, helper_delta, new_restart_delta, tmp_delta

        do i = 1, size(options)
            !See if we are the head of a nest chain, and if we have children
            if (options(i)%general%parent_nest == 0 .and. size(options(i)%general%child_nests) > 0) then
                call helper_delta%set(options(i)%restart%restart_count*options(i)%output%outputinterval)
                parent_restart_time = options(i)%general%start_time + helper_delta

                !Loop through all child nests in the chain
                do n = i+1, size(options)

                    child_nest_indx = n
                    ! If we have hit another root nest, then we are done with this chain
                    if (options(child_nest_indx)%general%parent_nest == 0) exit

                    if (child_nest_indx <= i) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: Child nest index ",child_nest_indx," is less than the parent nest index ",i
                        stop
                    endif

                    !! All forcing%inputinterval should be the same
                    if (.not.(options(i)%forcing%inputinterval == options(child_nest_indx)%forcing%inputinterval)) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: Forcing dt for nest ", child_nest_indx, " does not match the parent nest"
                        stop
                    endif

                    !! Start and end times of a nest should be within the span of the parent nest run time
                    if (options(i)%general%start_time > options(child_nest_indx)%general%start_time) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: Start time of nest ", child_nest_indx, " is before the start time of the parent nest"
                        stop
                    endif
                    if (options(i)%general%end_time < options(child_nest_indx)%general%end_time) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: End time of nest ", child_nest_indx, " is after the end time of the parent nest"
                        stop
                    endif

                    !! Ensure that the child nest dx is actually smaller than the parent nest
                    if (options(i)%domain%dx <= options(child_nest_indx)%domain%dx) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: dx of nest ", child_nest_indx, " is not smaller than the parent nest."
                        if (STD_OUT_PE) write(*,*) "  ERROR: This is sort of useless, so it was probably done in error"
                        stop
                    endif

                    !! Check that, for a given chain of nests, the time at which restart files are output is the same for all nests
                    !Compute when the first restart time is for the child nest
                    call helper_delta%set(options(child_nest_indx)%restart%restart_count*options(child_nest_indx)%output%outputinterval)
                    child_restart_time = options(child_nest_indx)%general%start_time + helper_delta
                    nest_start_offset = options(child_nest_indx)%general%start_time - options(i)%general%start_time

                    ! If the two times are different, warn the user
                    if ((parent_restart_time+nest_start_offset) /= child_restart_time) then
                        !Try to set the restart interval for the child nest to match the parent nest restart output time
                        new_restart_delta = parent_restart_time+nest_start_offset-options(child_nest_indx)%general%start_time
                        new_restartinterval = (new_restart_delta%seconds())/options(child_nest_indx)%output%outputinterval
                        
                        call helper_delta%set(new_restartinterval*options(child_nest_indx)%output%outputinterval)
                        child_restart_time = options(child_nest_indx)%general%start_time + helper_delta
                        if ((parent_restart_time+nest_start_offset) == child_restart_time) then
                            call set_nml_var(options(child_nest_indx)%restart%restart_count, new_restartinterval, 'restartinterval')
                            if (STD_OUT_PE) write(*,*) "  ATTENTION: Restart interval for nest ", child_nest_indx, " has been set to match the"
                            if (STD_OUT_PE) write(*,*) "  ATTENTION: first restart date of the parent nest."
                            if (STD_OUT_PE) write(*,*) "  ATTENTION: Restart interval for nest ", child_nest_indx, " is now: ", new_restartinterval
                            if (STD_OUT_PE) write(*,*) "  ATTENTION: First restart time for nest ", child_nest_indx, " is now: ", trim(child_restart_time%as_string())
                        else
                            !If this does not work, error out and report to the user
                            if (STD_OUT_PE) write(*,*) "  ERROR: Restart interval for nest ", child_nest_indx, " does not match the parent nest"
                            if (STD_OUT_PE) write(*,*) "  ERROR: Restart interval for nest ", child_nest_indx, " is: ", options(child_nest_indx)%restart%restart_count
                            if (STD_OUT_PE) write(*,*) "  ERROR: First restart time for nest ", child_nest_indx, " is: ", trim(child_restart_time%as_string())
                            if (STD_OUT_PE) write(*,*) "  ERROR: This time is calculated as start_time + restartinterval*outputinterval"
                            child_restart_time = parent_restart_time+nest_start_offset
                            if (STD_OUT_PE) write(*,*) "  ERROR: First restart time for nest ", child_nest_indx, " should be: ", trim(child_restart_time%as_string())
                            if (STD_OUT_PE) write(*,*) "  ERROR: Please adjust the restartinterval and outputinterval for this nest to match the parent nest"
                            stop
                        end if
                    end if

                    !! Check that, for a given chain of nests, the start times for all nests are separated by an integer
                    ! multiple of the input interval

                    ! Should take the difference between child and parent start times, and divide by the input interval
                    ! If this is not an integer, then the start times are not separated by an integer multiple of the input interval
                    tmp_delta = options(child_nest_indx)%general%start_time - options(i)%general%start_time
                    if (mod(tmp_delta%seconds(), options(child_nest_indx)%forcing%inputinterval) /= 0) then
                        if (STD_OUT_PE) write(*,*) "  ERROR: Start time for nest ", child_nest_indx, " is not an integer multiple"
                        if (STD_OUT_PE) write(*,*) "  ERROR: of the inputinterval from the start time of the parent nest"
                        stop
                    endif
                enddo
            endif
        enddo

    end subroutine inter_nest_options_check

    !> -------------------------------
    !! Checks options in the options data structure for consistency
    !!
    !! Stops or prints a large warning depending on warning level requested and error found
    !!
    !! -------------------------------
    subroutine check(this)
        ! Minimal error checking on option settings
        implicit none
        class(options_t), intent(inout)::this

        integer :: i

        ! Check that boundary conditions file exists
        if (trim(this%domain%init_conditions_file) /= '') then
            call check_file_exists(trim(this%domain%init_conditions_file), message='A boundary conditions file does not exist.')
        endif

        ! Check that the output and restart folders exists
        if (trim(this%output%output_folder) /= '') then
            call check_file_exists(trim(this%output%output_folder), message='Output folder does not exist.')
        endif
        if (trim(this%restart%restart_folder) /= '') then
            call check_file_exists(trim(this%restart%restart_folder), message='Restart folder does not exist.')
        endif

        !clean output var list
        do i=1, size(this%output%vars_for_output)
            if ((this%output%vars_for_output(i)+this%vars_for_restart(i) > 0) .and. (this%vars_to_allocate(i) <= 0)) then
                !if (STD_OUT_PE) write(*,*) 'variable ',trim(get_varname(this%vars_to_allocate(i))),' requested for output/restart, but was not allocated by one of the modules'
                if (this%output%vars_for_output(i) > 0) this%output%vars_for_output(i) = 0
                if (this%vars_for_restart(i) > 0) this%vars_for_restart(i) = 0
            endif
        enddo

        if (this%mp%top_mp_level < 0) this%mp%top_mp_level = this%domain%nz + this%mp%top_mp_level

        ! In read wind options, we set update_dt to be the FREQUENCY, not the actual dt. Compute the actual dt here
        call this%wind%update_dt%set(seconds=this%forcing%input_dt%seconds()/this%wind%update_dt%seconds())

        !Perform checks
        if (this%wind%smooth_wind_distance.eq.(-9999)) then
            this%wind%smooth_wind_distance=this%domain%dx*2
            if (STD_OUT_PE) write(*,*) "  Default smoothing distance = dx*2 = ", this%wind%smooth_wind_distance
        elseif (this%wind%smooth_wind_distance<0) then
            write(*,*) "  Wind smoothing must be a positive number"
            write(*,*) "  this%wind%smooth_wind_distance = ",this%wind%smooth_wind_distance
            this%wind%smooth_wind_distance = this%domain%dx*2
        endif

        !If user does not define this option, then let's Do The Right Thing
        if (this%lsm%nmp_opt_sfc == -1) then
            !If user has not turned on the surface layer scheme, then we set this to the recommended default value of 1
            if (this%physics%surfacelayer == 0) then
                this%lsm%nmp_opt_sfc = 1
            !If user has turned on the surface layer scheme, then let's opt to use the surface layer scheme's surface exchange coefficients
            else if (this%physics%surfacelayer == 1) then
                this%lsm%nmp_opt_sfc = 3
            endif
        endif

        if (.not.(this%physics%radiation == kRA_RRTMG)) this%pbl%ysu_topdown_pblmix = 0

        ! Change z length of snow arrays here, since we need to change their size for the output arrays, which are set in
        ! output_options_namelist
        if (this%physics%snowmodel==kSM_FSM) then
            kSNOW_GRID_Z = this%sm%fsm_nsnow_max
            kSNOWSOIL_GRID_Z = kSNOW_GRID_Z+kSOIL_GRID_Z
    endif

        ! if using a real LSM, feedback will probably keep hot-air from getting even hotter, so not likely a problem
        if ((this%physics%landsurface>0).and.(this%physics%boundarylayer==0)) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "  WARNING, Using surfaces fluxes (lsm>0) without a PBL scheme may overheat the surface and CRASH the model."
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
        endif

        ! if using a real LSM, feedback will probably keep hot-air from getting even hotter, so not likely a problem
        if ((this%physics%surfacelayer==0).and.(this%physics%boundarylayer>0)) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  ERROR, a surface layer scheme is required when using a PBL scheme,"
            if (STD_OUT_PE) write(*,*) "  ERROR, set sfc > 0 in the namelist."
        endif

        ! prior to v 0.9.3 this was assumed, so throw a warning now just in case.
        if ((this%forcing%z_is_geopotential .eqv. .False.).and. &
            (this%forcing%zvar=="PH")) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "  WARNING z variable is not assumed to be geopotential height when it is 'PH'."
            if (STD_OUT_PE) write(*,*) "  WARNING If z is geopotential, set z_is_geopotential=True in the namelist."
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
        endif
        if ((this%forcing%z_is_geopotential .eqv. .True.).and. &
            (this%forcing%z_is_on_interface .eqv. .False.)) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "  WARNING geopotential height is no longer assumed to be on model interface levels."
            if (STD_OUT_PE) write(*,*) "  WARNING To interpolate geopotential, set z_is_on_interface=True in the namelist. "
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
        endif
        
        !! MJ added
        if ((this%physics%radiation_downScaling==1).and.(this%physics%radiation==0)) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  STOP STOP STOP"
            if (STD_OUT_PE) write(*,*) "  STOP, Running radiation_downScaling=1 cannot not be used with rad=0"
            if (STD_OUT_PE) write(*,*) "  STOP STOP STOP"
            stop
        endif
        
        if (this%time%RK3) then
            if (max(this%adv%h_order,this%adv%v_order)==5 .and. this%time%cfl_reduction_factor > 1.4) then
                if (STD_OUT_PE) write(*,*) "  CFL reduction factor should be less than 1.4 when horder or vorder = 5, limiting to 1.4"
                this%time%cfl_reduction_factor = min(1.4,this%time%cfl_reduction_factor)
            elseif (max(this%adv%h_order,this%adv%v_order)==3 .and. this%time%cfl_reduction_factor > 1.6) then
                if (STD_OUT_PE) write(*,*) "  CFL reduction factor should be less than 1.6 when horder or vorder = 3, limiting to 1.6"
                this%time%cfl_reduction_factor = min(1.6,this%time%cfl_reduction_factor)
            endif
        else
            if (this%time%cfl_reduction_factor > 1.0) then   
                if (STD_OUT_PE) write(*,*) "  CFL reduction factor should be less than 1.0 when RK3=.False., limiting to 1.0"
                this%time%cfl_reduction_factor = min(1.0,this%time%cfl_reduction_factor)
            endif
        endif
        
        if (this%wind%alpha_const > 0) then
            if (this%wind%alpha_const > 1.0) then
                if (STD_OUT_PE) write(*,*) "  Alpha currently limited to values between 0.01 and 1.0, setting to 1.0"
                this%wind%alpha_const = 1.0
            else if (this%wind%alpha_const < 0.01) then
                if (STD_OUT_PE) write(*,*) "  Alpha currently limited to values between 0.01 and 1.0, setting to 0.01"
                this%wind%alpha_const = 0.01
            endif
        endif
        
        ! should warn user if lsm is run without radiation
        if ((this%physics%landsurface>kLSM_BASIC .or. this%physics%snowmodel>0).and.(this%physics%radiation==0)) then
            if (STD_OUT_PE) write(*,*) "  "
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "  WARNING, Using land surface model without radiation input does not make sense."
            if (STD_OUT_PE) write(*,*) "  WARNING WARNING WARNING"
        endif

        if (this%physics%landsurface>0) then
            this%sfc%isfflx = 1
            this%sfc%scm_force_flux = 1
        endif
        
        !Allow for microphysics precipitation partitioning with NoahMP if using a snow model
        if (this%physics%snowmodel>0) then
            this%lsm%nmp_opt_snf = 4
        endif

        ! check if the last entry in dz_levels is zero, which would indicate that nz is larger than the number
        ! of entries in dz_levels, or that the user passed bad data
        if (this%domain%nz > 1) then
            if (this%domain%dz_levels(this%domain%nz) == 0) then
                if (STD_OUT_PE) write(*,*) "  nz is larger than the number of entries in dz_levels, or the last entry in dz_levels is zero."
                stop
            endif
        endif

        ! check if start time is before end time
        if (this%general%start_time >= this%general%end_time) then
            if (STD_OUT_PE) write(*,*) "  Start time must be before end time"
            stop
        endif

        ! check if restart_time is between start and end time
        if (this%restart%restart) then
            if (this%restart%restart_time < this%general%start_time .or. &
                this%restart%restart_time > this%general%end_time) then
                if (STD_OUT_PE) write(*,*) "  Restart time must be between start and end time"
                stop
            endif
        endif

        ! Check if supporting files exist, if they are needed by physics modules
        if (this%physics%landsurface==kLSM_NOAHMP) then
            if (STD_OUT_PE) write(*,*) '  NoahMP LSM turned on, checking for supporting files...'
            call check_file_exists('GENPARM.TBL', message='GENPARM.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('MPTABLE.TBL', message='MPTABLE.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('SOILPARM.TBL', message='SOILPARM.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('VEGPARM.TBL', message='VEGPARM.TBL file does not exist. This should be in the same directory as the namelist.')
        endif
        if (this%physics%radiation==kRA_RRTMG) then
            if (STD_OUT_PE) write(*,*) '  RRTMG radiation turned on, checking for supporting files...'
            call check_file_exists('rrtmg_support/forrefo_1.nc', message='At least one of the RRTMG supporting files does not exist. These files should be in a folder "rrtmg_support" placed in the same directory as the namelist.')
        endif
        if (this%physics%microphysics==kMP_ISHMAEL) then
            if (STD_OUT_PE) write(*,*) '  ISHMAEL microphysics turned on, checking for supporting files...'
            call check_file_exists('mp_support/ishmael_gamma_tab.nc', message='At least one of the ISHMAEL supporting files does not exist. These files should be in a folder "mp_support" placed in the same directory as the namelist.')
        endif

        ! check that input variables are specified when relevant physics options are chosen
        if (this%physics%landsurface > kLSM_BASIC) then
            call require_var(this%forcing%latvar, 'latvar', 'This variable is required when running with this LSM option')
            !call require_var(this%forcing%sst_var, 'sst_var')
            !call require_var(this%forcing%time_var, 'time_var')
        endif

        if (this%physics%snowmodel == kSM_FSM) then
            !call require_var(this%forcing%qcvar, 'qcvar', 'This variable is required when running with the FSM snowmodel option')
            !call require_var(this%forcing%qngvar, 'qngvar')
            !call require_var(this%forcing%qnsvar, 'qnsvar')
        endif

        if (trim(this%lsm%LU_Categories)=="USGS") then
            if((this%physics%watersurface==kWATER_LAKE) .AND. (STD_OUT_PE)) then
                write(*,*) "  WARNING: Lake model selected (water=2), but USGS LU-categories has no lake category"
            endif
        elseif (trim(this%lsm%LU_Categories)=="NLCD40") then
            if(this%physics%watersurface==kWATER_LAKE) write(*,*) "  WARNING: Lake model selected (water=2), but NLCD40 LU-categories has no lake category"
        endif

        ! There needs to be a unique domain file for each nest. Additionally, dx needs to be set for each nest. Check for these here.
        if (trim(this%domain%init_conditions_file)=="") then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error: 'init_conditions_file' must be set in the domain namelist for each nest"
            if (STD_OUT_PE) write(*,*) "  Error: missing for nest: ", this%nest_indx
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        if (this%domain%dx<=0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error: 'dx' must be set in the domain namelist for each nest"
            if (STD_OUT_PE) write(*,*) "  Error: missing for nest: ", this%nest_indx
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif


        ! -------------------------------------
        ! Restart Checks
        ! -------------------------------------
        ! Check that the restart interval is a multiple of the input interval

    end subroutine check

    !> -------------------------------
    !! Read physics options to use from a namelist file
    !!
    !! -------------------------------
    subroutine physics_namelist(filename, phys_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),intent(in)     :: filename
        type(physics_type), intent(inout) :: phys_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml
        integer :: name_unit, rc
        !variables to be used in the namelist
        integer, dimension(kMAX_NESTS) :: pbl, lsm, mp, sfc, sm, water, rad, conv, adv, wind, radiation_downscaling
        logical :: print_info, gennml
        !define the namelist
        namelist /physics/ pbl, lsm, sfc, sm, water, mp, rad, conv, adv, wind, radiation_downscaling

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(pbl, 'pbl', print_info, gennml)
        call set_nml_var_default(lsm, 'lsm', print_info, gennml)
        call set_nml_var_default(sfc, 'sfc', print_info, gennml)
        call set_nml_var_default(sm, 'sm', print_info, gennml)
        call set_nml_var_default(water, 'water', print_info, gennml)
        call set_nml_var_default(mp, 'mp', print_info, gennml)
        call set_nml_var_default(rad, 'rad', print_info, gennml)
        call set_nml_var_default(conv, 'conv', print_info, gennml)
        call set_nml_var_default(adv, 'adv', print_info, gennml)
        call set_nml_var_default(wind, 'wind', print_info, gennml)
        call set_nml_var_default(radiation_downscaling, 'radiation_downscaling', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        !read the namelist
        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=physics)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'physics' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        !store options
        call set_nml_var(phys_options%boundarylayer, pbl(n_indx), 'pbl', pbl(1))
        call set_nml_var(phys_options%landsurface, lsm(n_indx), 'lsm', lsm(1))
        call set_nml_var(phys_options%surfacelayer, sfc(n_indx), 'sfc', sfc(1))
        call set_nml_var(phys_options%snowmodel, sm(n_indx), 'sm', sm(1))
        call set_nml_var(phys_options%watersurface, water(n_indx), 'water', water(1))
        call set_nml_var(phys_options%microphysics, mp(n_indx), 'mp',mp(1))
        call set_nml_var(phys_options%radiation, rad(n_indx), 'rad', rad(1))
        call set_nml_var(phys_options%convection, conv(n_indx), 'conv', conv(1))
        call set_nml_var(phys_options%advection, adv(n_indx), 'adv', adv(1))
        call set_nml_var(phys_options%windtype, wind(n_indx), 'wind', wind(1))
        call set_nml_var(phys_options%radiation_downScaling, radiation_downscaling(n_indx), 'radiation_downscaling', radiation_downscaling(1))

    end subroutine physics_namelist


    !> -------------------------------
    !! Check that a required input variable is present
    !!
    !! If not present, halt the program
    !!
    !! -------------------------------
    subroutine require_var(inputvar, var_name, message)
        implicit none
        character(len=*), intent(in) :: inputvar
        character(len=*), intent(in) :: var_name
        character(len=*), optional, intent(in) :: message

        if (trim(inputvar)=="") then
            if (STD_OUT_PE) write(*,*) "  Variable: ",trim(var_name), " is required."
            if (STD_OUT_PE .and. present(message)) write(*,*) "  Variable: ",trim(var_name), " ",trim(message)
            stop
        endif

    end subroutine require_var

    !> -------------------------------
    !! Initialize the variable names to be written to standard output
    !!
    !! Reads the output_list namelist
    !!
    !! -------------------------------
    subroutine output_namelist(filename, output_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(output_options_type), intent(inout) :: output_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, i, j, status, var_indx
        integer :: frames_per_outfile(kMAX_NESTS)
        real    :: outputinterval(kMAX_NESTS)
        logical :: file_exists, print_info, gennml

        ! Local variables
        character(len=kMAX_FILE_LENGTH) :: output_folder(kMAX_NESTS)
        character(len=kMAX_NAME_LENGTH) :: output_vars(kMAX_STORAGE_VARS, kMAX_NESTS)
        
        namelist /output/ output_vars, outputinterval, frames_per_outfile, output_folder

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(outputinterval, 'outputinterval', print_info, gennml)
        call set_nml_var_default(frames_per_outfile, 'frames_per_outfile', print_info, gennml)
        call set_nml_var_default(output_folder, 'output_folder', print_info, gennml)
        call set_nml_var_default(output_vars, 'output_vars', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc, nml=output)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'output' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        do j=1, kMAX_STORAGE_VARS
            if (trim(output_vars(j, n_indx)) /= "" .and. trim(output_vars(j, n_indx)) /= kCHAR_NO_VAL ) then

                !get the var index for this output variable name
                var_indx = get_varindx(trim(output_vars(j, n_indx)))
                if (var_indx <= kMAX_STORAGE_VARS) call add_to_varlist(output_options%vars_for_output, [var_indx])
            endif
        enddo        

        call set_nml_var(output_options%output_folder, output_folder(n_indx), 'output_folder', output_folder(1))
        call set_nml_var(output_options%outputinterval, outputinterval(n_indx), 'outputinterval', outputinterval(1))
        call output_options%output_dt%set(seconds=outputinterval(n_indx))
        call set_nml_var(output_options%frames_per_outfile, frames_per_outfile(n_indx), 'frames_per_outfile', frames_per_outfile(1))

    end subroutine output_namelist

    !> -------------------------------
    !! Initialize the variable names to be written to standard output
    !!
    !! Reads the restart_list namelist
    !!
    !! -------------------------------
    subroutine restart_namelist(filename, options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),intent(in)    :: filename
        type(options_t), intent(inout) :: options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        integer    :: restartinterval(kMAX_NESTS)

        logical :: file_exists, print_info, gennml, restart_run(kMAX_NESTS)

        ! Local variables
        character(len=kMAX_FILE_LENGTH) :: restart_folder(kMAX_NESTS)
        character(len=kMAX_FILE_LENGTH) :: restart_date(kMAX_NESTS)    ! date to restart

        namelist /restart/  restartinterval, restart_folder, restart_date, restart_run

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(restartinterval, 'restartinterval', print_info, gennml)
        call set_nml_var_default(restart_folder, 'restart_folder', print_info, gennml)
        call set_nml_var_default(restart_date, 'restart_date', print_info, gennml)
        call set_nml_var_default(restart_run, 'restart_run', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc, nml=restart)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'restart' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        ! Restart interval for this particular nest must be set to be the number of output intervals for this nest
        ! which results in the same restart timestep for the first nest
        call set_nml_var(options%restart%restart_count, restartinterval(n_indx), 'restartinterval', restartinterval(1))
        call set_nml_var(options%restart%restart_folder, restart_folder(n_indx), 'restart_folder', restart_folder(1))
        call set_nml_var(options%restart%restart, restart_run(n_indx), 'restart_run', restart_run(1))

        !If the user did not ask for a restart run, leave the function now
        if (.not.(options%restart%restart)) return
        
        ! calculate the modified julian day for th restart date given
        call options%restart%restart_time%init(options%general%calendar)
        if (restart_date(1)=="") then
            if (STD_OUT_PE) write(*,*) "  ERROR: restart_date must be specified in the namelist"
            stop
        else
            call options%restart%restart_time%set(restart_date(1))
        endif
        
        
    end subroutine restart_namelist

    !> -------------------------------
    !! Initialize the variable names to be read
    !!
    !! Reads the var_list namelist
    !!
    !! -------------------------------
    subroutine forcing_namelist(filename, options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),     intent(in)    :: filename
        type(options_t),      intent(inout) :: options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, i, j, nfiles
        logical :: compute_p, print_info, gennml, no_check
        logical, dimension(kMAX_NESTS) :: limit_rh, z_is_geopotential, z_is_on_interface,&
                   time_varying_z, t_is_potential, qv_is_spec_humidity, &
                   qv_is_relative_humidity
        real, dimension(kMAX_NESTS)    :: t_offset, inputinterval

        character(len=kMAX_FILE_LENGTH) :: forcing_file_list
        character(len=kMAX_FILE_LENGTH), allocatable :: boundary_files(:)

        character(len=kMAX_NAME_LENGTH) :: latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,wvar,zvar,  &
                                        pvar,tvar,qvvar,qcvar,qivar,qrvar,qgvar,qsvar,            &
                                        qncvar,qnivar,qnrvar,qngvar,qnsvar,hgtvar,shvar,lhvar,pblhvar,  &
                                        i2mvar, i3mvar, i2nvar, i3nvar, i1avar, i2avar, i3avar, i1cvar, i2cvar, i3cvar, &
                                        psvar, pslvar, swdown_var, lwdown_var, sst_var, time_var

        namelist /forcing/ forcing_file_list, inputinterval, t_offset, limit_rh, z_is_geopotential, z_is_on_interface, time_varying_z, &
                            t_is_potential, qv_is_relative_humidity, qv_is_spec_humidity, &
                            pvar,tvar,qvvar,qcvar,qivar,qrvar,qgvar,qsvar,qncvar,qnivar,qnrvar,qngvar,qnsvar,&
                            i2mvar, i3mvar, i2nvar, i3nvar, i1avar, i2avar, i3avar, i1cvar, i2cvar, i3cvar, &
                            hgtvar,shvar,lhvar,pblhvar,   &
                            latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,wvar,zvar, &
                            psvar, pslvar, swdown_var, lwdown_var, sst_var, time_var

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.    
        if (present(gen_nml)) gennml = gen_nml

        ! Default values for forcing options
        allocate(boundary_files(MAX_NUMBER_FILES))

        call set_nml_var_default(forcing_file_list, 'forcing_file_list', print_info, gennml)
        call set_nml_var_default(t_offset, 't_offset', print_info, gennml)
        call set_nml_var_default(inputinterval, 'inputinterval', print_info, gennml)
        call set_nml_var_default(limit_rh, 'limit_rh', print_info, gennml)
        call set_nml_var_default(z_is_geopotential, 'z_is_geopotential', print_info, gennml)
        call set_nml_var_default(z_is_on_interface, 'z_is_on_interface', print_info, gennml)
        call set_nml_var_default(time_varying_z, 'time_varying_z', print_info, gennml)
        call set_nml_var_default(t_is_potential, 't_is_potential', print_info, gennml)
        call set_nml_var_default(qv_is_relative_humidity, 'qv_is_relative_humidity', print_info, gennml)
        call set_nml_var_default(qv_is_spec_humidity, 'qv_is_spec_humidity', print_info, gennml)
        call set_nml_var_default(latvar, 'latvar', print_info, gennml)
        call set_nml_var_default(lonvar, 'lonvar', print_info, gennml)
        call set_nml_var_default(hgtvar, 'hgtvar', print_info, gennml)
        call set_nml_var_default(zvar, 'zvar', print_info, gennml)
        call set_nml_var_default(uvar, 'uvar', print_info, gennml)
        call set_nml_var_default(ulat, 'ulat', print_info, gennml)
        call set_nml_var_default(ulon, 'ulon', print_info, gennml)
        call set_nml_var_default(vvar, 'vvar', print_info, gennml)
        call set_nml_var_default(vlat, 'vlat', print_info, gennml)
        call set_nml_var_default(vlon, 'vlon', print_info, gennml)
        call set_nml_var_default(wvar, 'wvar', print_info, gennml)
        call set_nml_var_default(pslvar, 'pslvar', print_info, gennml)
        call set_nml_var_default(psvar, 'psvar', print_info, gennml)
        call set_nml_var_default(pvar, 'pvar', print_info, gennml)
        call set_nml_var_default(tvar, 'tvar', print_info, gennml)
        call set_nml_var_default(qvvar, 'qvvar', print_info, gennml)
        call set_nml_var_default(qcvar, 'qcvar', print_info, gennml)
        call set_nml_var_default(qivar, 'qivar', print_info, gennml)
        call set_nml_var_default(qrvar, 'qrvar', print_info, gennml)
        call set_nml_var_default(qsvar, 'qsvar', print_info, gennml)
        call set_nml_var_default(qgvar, 'qgvar', print_info, gennml)
        call set_nml_var_default(i2mvar, 'i2mvar', print_info, gennml)
        call set_nml_var_default(i3mvar, 'i3mvar', print_info, gennml)
        call set_nml_var_default(qncvar, 'qncvar', print_info, gennml)
        call set_nml_var_default(qnivar, 'qnivar', print_info, gennml)
        call set_nml_var_default(qnrvar, 'qnrvar', print_info, gennml)
        call set_nml_var_default(qnsvar, 'qnsvar', print_info, gennml)
        call set_nml_var_default(qngvar, 'qngvar', print_info, gennml)
        call set_nml_var_default(i2nvar, 'i2nvar', print_info, gennml)
        call set_nml_var_default(i3nvar, 'i3nvar', print_info, gennml)
        call set_nml_var_default(i1avar, 'i1avar', print_info, gennml)
        call set_nml_var_default(i2avar, 'i2avar', print_info, gennml)
        call set_nml_var_default(i3avar, 'i3avar', print_info, gennml)
        call set_nml_var_default(i1cvar, 'i1cvar', print_info, gennml)
        call set_nml_var_default(i2cvar, 'i2cvar', print_info, gennml)
        call set_nml_var_default(i3cvar, 'i3cvar', print_info, gennml)
        call set_nml_var_default(shvar, 'shvar', print_info, gennml)
        call set_nml_var_default(lhvar, 'lhvar', print_info, gennml)
        call set_nml_var_default(swdown_var, 'swdown_var', print_info, gennml)
        call set_nml_var_default(lwdown_var, 'lwdown_var', print_info, gennml)
        call set_nml_var_default(sst_var, 'sst_var', print_info, gennml)
        call set_nml_var_default(pblhvar, 'pblhvar', print_info, gennml)
        call set_nml_var_default(time_var, 'time_var', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=forcing)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'forcing' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        compute_p = .False.
        if ((pvar=="") .and. ((pslvar/="") .or. (psvar/=""))) compute_p = .True.
        if (compute_p) then
            if ((pslvar == "").and.(hgtvar == "")) then
                write(*,*) "  ERROR: if surface pressure is used to compute air pressure, then surface height must be specified"
                error stop
            endif
        endif

        options%forcing%compute_z = .False.
        if ((zvar=="") .and. ((pslvar/="") .or. (psvar/=""))) options%forcing%compute_z = .True.
        if (options%forcing%compute_z) then
            if (pvar=="") then
                if (STD_OUT_PE) write(*,*) "  ERROR: either pressure (pvar) or atmospheric level height (zvar) must be specified"
                error stop
            endif
        endif

        call set_nml_var(options%forcing%t_offset, t_offset(n_indx), 't_offset', t_offset(1))
        call set_nml_var(options%forcing%limit_rh, limit_rh(n_indx), 'limit_rh', limit_rh(1))
        call set_nml_var(options%forcing%z_is_geopotential, z_is_geopotential(n_indx), 'z_is_geopotential', z_is_geopotential(1))
        call set_nml_var(options%forcing%z_is_on_interface, z_is_on_interface(n_indx), 'z_is_on_interface', z_is_on_interface(1))
        call set_nml_var(options%forcing%time_varying_z, time_varying_z(n_indx), 'time_varying_z', time_varying_z(1))
        call set_nml_var(options%forcing%t_is_potential, t_is_potential(n_indx), 't_is_potential', t_is_potential(1))
        call set_nml_var(options%forcing%qv_is_relative_humidity, qv_is_relative_humidity(n_indx), 'qv_is_relative_humidity', qv_is_relative_humidity(1))
        call set_nml_var(options%forcing%qv_is_spec_humidity, qv_is_spec_humidity(n_indx), 'qv_is_spec_humidity', qv_is_spec_humidity(1))
        call set_nml_var(options%forcing%inputinterval, inputinterval(n_indx), 'inputinterval', inputinterval(1))

        call options%forcing%input_dt%set(seconds=options%forcing%inputinterval)

        call require_var(lonvar, "Longitude")
        call require_var(latvar, "Latitude")
        call require_var(hgtvar, "Terrain Height")
        call require_var(zvar, "Verticle Level Height")
        call require_var(uvar, "U winds")
        call require_var(vvar, "V winds")
        call require_var(tvar, "Temperature")
        call require_var(qvvar, "Water Vapor Mixing Ratio")
        call require_var(time_var, "Time")

        call check_file_exists(forcing_file_list, message="Forcing file list does not exist.")

        ! See if we should check for these variables in the boundary conditions files
        no_check = (options%general%parent_nest > 0)

        if (.not.(no_check)) then
            nfiles = read_forcing_file_names(forcing_file_list, boundary_files)

            if (nfiles==0) then
                stop "No boundary conditions files specified."
            endif

            allocate(options%forcing%boundary_files(nfiles))
            options%forcing%boundary_files(1:nfiles) = boundary_files(1:nfiles)
            deallocate(boundary_files)
        endif

        ! NOTE: temperature must be the first of the forcing variables read
        options%forcing%vars_to_read(:) = ""
        options%forcing%dim_list(:) = 0
        i = 1
        call set_nml_var(options%forcing%tvar, tvar, 'tvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%latvar, latvar, 'latvar', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%lonvar, lonvar, 'lonvar', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%hgtvar, hgtvar, 'hgtvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%uvar, uvar, 'uvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%ulat, ulat, 'ulat', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%ulon, ulon, 'ulon', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%vvar, vvar, 'vvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%vlat, vlat, 'vlat', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%vlon, vlon, 'vlon', options%forcing, no_check=no_check)
        call set_nml_var(options%forcing%wvar, wvar, 'wvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%pslvar, pslvar, 'pslvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%psvar, psvar, 'psvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qvvar, qvvar, 'qvvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qcvar, qcvar, 'qcvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qivar, qivar, 'qivar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qrvar, qrvar, 'qrvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qgvar, qgvar, 'qgvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qsvar, qsvar, 'qsvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qncvar, qncvar, 'qncvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qnivar, qnivar, 'qnivar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qnrvar, qnrvar, 'qnrvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qngvar, qngvar, 'qngvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%qnsvar, qnsvar, 'qnsvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i2mvar, i2mvar, 'i2mvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i3mvar, i3mvar, 'i3mvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i2nvar, i2nvar, 'i2nvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i3nvar, i3nvar, 'i3nvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i1avar, i1avar, 'i1avar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i2avar, i2avar, 'i2avar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i3avar, i3avar, 'i3avar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i1cvar, i1cvar, 'i1cvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i2cvar, i2cvar, 'i2cvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%i3cvar, i3cvar, 'i3cvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%shvar, shvar, 'shvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%lhvar, lhvar, 'lhvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%swdown_var, swdown_var, 'swdown_var', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%lwdown_var, lwdown_var, 'lwdown_var', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%sst_var, sst_var, 'sst_var', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%pblhvar, pblhvar, 'pblhvar', options%forcing, i, no_check=no_check)
        call set_nml_var(options%forcing%time_var, time_var, 'time_var', options%forcing, no_check=no_check)

        ! vertical coordinate
        ! if (options%forcing%time_varying_z) then
        if (options%forcing%compute_z) then
            zvar = "height_computed"
            options%forcing%zvar        = zvar; options%forcing%vars_to_read(i) = zvar;      options%forcing%dim_list(i) = -3;    i = i + 1
        else
            call set_nml_var(options%forcing%zvar, zvar, 'zvar', options%forcing, i, no_check=no_check)
        endif


        if (compute_p) then
            pvar = "air_pressure_computed"
            options%forcing%pvar        = pvar  ; options%forcing%vars_to_read(i) = pvar;       options%forcing%dim_list(i) = -3;   i = i + 1
        else
            call set_nml_var(options%forcing%pvar, pvar, 'pvar', options%forcing, i, no_check=no_check)
        endif


    end subroutine forcing_namelist


    !> -------------------------------
    !! Initialize the main parameter options
    !!
    !! Reads parameters for the ICAR simulation
    !! These include setting flags that request other namelists be read
    !!
    !! -------------------------------
    subroutine general_namelist(filename, gen_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(general_options_type), intent(inout) :: gen_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, count, i
        integer, dimension(kMAX_NESTS) :: parent_nest, nests
        logical :: print_info, gennml

        ! parameters to read
        logical, dimension(kMAX_NESTS) :: interactive, debug, &
                                            use_mp_options, use_lt_options, use_adv_options, use_lsm_options, &
                                            use_cu_options, use_rad_options, use_pbl_options, use_sfc_options, &
                                            use_wind_options, use_sm_options

        character(len=kMAX_FILE_LENGTH), dimension(kMAX_NESTS) :: start_date, end_date, calendar, version, comment, phys_suite
        character(len=kMAX_NAME_LENGTH)                        :: start_date_checked, end_date_checked
        namelist /general/    debug, interactive, calendar,          &
                              version, comment, phys_suite,                 &
                              start_date, end_date, &
                              nests, parent_nest, &
                              use_mp_options,     &
                              use_lt_options,     &
                              use_lsm_options,    &
                              use_sm_options,     &
                              use_adv_options,    &
                              use_cu_options,     &
                              use_rad_options,    &
                              use_sfc_options,    &
                              use_wind_options,   &
                              use_pbl_options

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(version, 'version', print_info, gennml)
        call set_nml_var_default(comment, 'comment', print_info, gennml)
        call set_nml_var_default(debug, 'debug', print_info, gennml)
        call set_nml_var_default(interactive, 'interactive', print_info, gennml)
        call set_nml_var_default(calendar, 'calendar', print_info, gennml)
        call set_nml_var_default(start_date, 'start_date', print_info, gennml)
        call set_nml_var_default(end_date, 'end_date', print_info, gennml)
        call set_nml_var_default(phys_suite, 'phys_suite', print_info, gennml)
        call set_nml_var_default(nests, 'nests', print_info, gennml)
        call set_nml_var_default(parent_nest, 'parent_nest', print_info, gennml)
        call set_nml_var_default(use_mp_options, 'use_mp_options', print_info, gennml)
        call set_nml_var_default(use_lt_options, 'use_lt_options', print_info, gennml)
        call set_nml_var_default(use_adv_options, 'use_adv_options', print_info, gennml)
        call set_nml_var_default(use_cu_options, 'use_cu_options', print_info, gennml)
        call set_nml_var_default(use_lsm_options, 'use_lsm_options', print_info, gennml)
        call set_nml_var_default(use_sm_options, 'use_sm_options', print_info, gennml)
        call set_nml_var_default(use_rad_options, 'use_rad_options', print_info, gennml)
        call set_nml_var_default(use_pbl_options, 'use_pbl_options', print_info, gennml)
        call set_nml_var_default(use_sfc_options, 'use_sfc_options', print_info, gennml)
        call set_nml_var_default(use_wind_options, 'use_wind_options', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=general)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'general' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        call set_nml_var(gen_options%calendar, calendar(n_indx), 'calendar', calendar(1))
        call set_nml_var(gen_options%version, version(n_indx), 'version', version(1))
        call set_nml_var(gen_options%comment, comment(n_indx), 'comment', comment(1))
        call set_nml_var(gen_options%phys_suite, phys_suite(n_indx), 'phys_suite',  phys_suite(1))
        call set_nml_var(gen_options%debug, debug(n_indx), 'debug', debug(1))
        call set_nml_var(gen_options%interactive, interactive(n_indx), 'interactive', interactive(1)) 
        call set_nml_var(gen_options%nests, nests(n_indx), 'nests',  nests(1))
        call set_nml_var(gen_options%use_mp_options, use_mp_options(n_indx), 'use_mp_options', use_mp_options(1))
        call set_nml_var(gen_options%use_lt_options, use_lt_options(n_indx), 'use_lt_options', use_lt_options(1))
        call set_nml_var(gen_options%use_adv_options, use_adv_options(n_indx), 'use_adv_options', use_adv_options(1))
        call set_nml_var(gen_options%use_lsm_options, use_lsm_options(n_indx), 'use_lsm_options', use_lsm_options(1))
        call set_nml_var(gen_options%use_sm_options, use_sm_options(n_indx), 'use_sm_options', use_sm_options(1))
        call set_nml_var(gen_options%use_cu_options, use_cu_options(n_indx), 'use_cu_options', use_cu_options(1))
        call set_nml_var(gen_options%use_rad_options, use_rad_options(n_indx), 'use_rad_options', use_rad_options(1))
        call set_nml_var(gen_options%use_pbl_options, use_pbl_options(n_indx), 'use_pbl_options', use_pbl_options(1))
        call set_nml_var(gen_options%use_sfc_options, use_sfc_options(n_indx), 'use_sfc_options', use_sfc_options(1))
        call set_nml_var(gen_options%use_wind_options, use_wind_options(n_indx), 'use_wind_options', use_wind_options(1))
        call set_nml_var(start_date_checked, start_date(n_indx), 'start_date', start_date(1))
        call set_nml_var(end_date_checked, end_date(n_indx), 'end_date', end_date(1))

        if (trim(start_date_checked)/="") then
            call gen_options%start_time%init(gen_options%calendar)
            call gen_options%start_time%set(start_date_checked)
        else
            stop 'start date must be supplied in namelist'
        endif
        if (trim(end_date_checked)/="") then
            call gen_options%end_time%init(gen_options%calendar)
            call gen_options%end_time%set(end_date_checked)
        else
            stop 'end date must be supplied in namelist'
        endif

        if (parent_nest(n_indx) >= n_indx) then
            if (STD_OUT_PE) write(*,*) "  ERROR for nest ", n_indx, ": parent nest must be less than or equal to the current nest"
            if (STD_OUT_PE) write(*,*) "  ERROR for nest ", n_indx, ": parent nest is ", parent_nest(n_indx)
            stop
        else
            call set_nml_var(gen_options%parent_nest, parent_nest(n_indx), 'parent_nest')
        endif

        count = 0
        do i = 1, gen_options%nests
            if (parent_nest(i) == n_indx) count = count + 1
        end do

        allocate(gen_options%child_nests(count))
        if (count == 0) return

        count = 0
        do i = 1, gen_options%nests
            if (parent_nest(i) == n_indx) then
                count = count + 1
                gen_options%child_nests(count) = i
            endif
        end do


        ! options are updated when complete
    end subroutine general_namelist


    !> -------------------------------
    !! Set up model levels, either read from a namelist, or from a default set of values
    !!
    !! Reads the z_info namelist or sets default values
    !!
    !! -------------------------------
    subroutine domain_namelist(filename, domain_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(domain_options_type), intent(inout) :: domain_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, this_level
        logical :: print_info, gennml

        real, dimension(MAXLEVELS, kMAX_NESTS) :: dz_levels
        logical, dimension(kMAX_NESTS) :: sleve, use_agl_height

        real, dimension(kMAX_NESTS) :: dx, flat_z_height, decay_rate_L_topo, decay_rate_S_topo, sleve_n, agl_cap, max_agl_height
        integer, dimension(kMAX_NESTS) :: nz, longitude_system, terrain_smooth_windowsize, terrain_smooth_cycles

        character(len=kMAX_FILE_LENGTH) :: init_conditions_file(kMAX_NESTS)

        character(len=kMAX_NAME_LENGTH), dimension(kMAX_NESTS) :: landvar,lakedepthvar,hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                                        snowh_var, soiltype_var, soil_t_var,soil_vwc_var,swe_var, soil_deept_var,           &
                                        vegtype_var,vegfrac_var, vegfracmax_var, albedo_var, lai_var, canwat_var, linear_mask_var, nsq_calibration_var,  &
                                        sinalpha_var, cosalpha_var, svf_var, hlm_var, slope_var, slope_angle_var, &
                                        aspect_angle_var, shd_var  !!MJ added

        namelist /domain/ dx, nz, longitude_system, init_conditions_file, &
                            landvar,lakedepthvar, snowh_var, agl_cap, use_agl_height, &
                            hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                            soiltype_var, soil_t_var,soil_vwc_var,swe_var,soil_deept_var,           &
                            vegtype_var,vegfrac_var, vegfracmax_var, albedo_var, lai_var, canwat_var, linear_mask_var, nsq_calibration_var,  &
                            sinalpha_var, cosalpha_var, svf_var, hlm_var, slope_var, slope_angle_var, aspect_angle_var, shd_var, & !! MJ added
                            dz_levels, flat_z_height, sleve, terrain_smooth_windowsize, terrain_smooth_cycles, decay_rate_L_topo, decay_rate_S_topo, sleve_n


        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(init_conditions_file, 'init_conditions_file', print_info, gennml)
        call set_nml_var_default(dx, 'dx', print_info, gennml)
        call set_nml_var_default(longitude_system, 'longitude_system', print_info, gennml)
        call set_nml_var_default(nz, 'nz', print_info, gennml)
        call set_nml_var_default(flat_z_height, 'flat_z_height', print_info, gennml)
        call set_nml_var_default(sleve, 'sleve', print_info, gennml)
        call set_nml_var_default(terrain_smooth_windowsize, 'terrain_smooth_windowsize', print_info, gennml)
        call set_nml_var_default(terrain_smooth_cycles, 'terrain_smooth_cycles', print_info, gennml)
        call set_nml_var_default(decay_rate_L_topo, 'decay_rate_L_topo', print_info, gennml)
        call set_nml_var_default(decay_rate_S_topo, 'decay_rate_S_topo', print_info, gennml)
        call set_nml_var_default(sleve_n, 'sleve_n', print_info, gennml)
        call set_nml_var_default(use_agl_height, 'use_agl_height', print_info, gennml)
        call set_nml_var_default(agl_cap, 'agl_cap', print_info, gennml)

        call set_nml_var_default(hgt_hi, 'hgt_hi', print_info, gennml)
        call set_nml_var_default(landvar, 'landvar', print_info, gennml)
        call set_nml_var_default(lakedepthvar, 'lakedepthvar', print_info, gennml)
        call set_nml_var_default(lat_hi, 'lat_hi', print_info, gennml)
        call set_nml_var_default(lon_hi, 'lon_hi', print_info, gennml)
        call set_nml_var_default(ulat_hi, 'ulat_hi', print_info, gennml)
        call set_nml_var_default(ulon_hi, 'ulon_hi', print_info, gennml)
        call set_nml_var_default(vlat_hi, 'vlat_hi', print_info, gennml)
        call set_nml_var_default(vlon_hi, 'vlon_hi', print_info, gennml)
        call set_nml_var_default(soiltype_var, 'soiltype_var', print_info, gennml)
        call set_nml_var_default(soil_t_var, 'soil_t_var', print_info, gennml)
        call set_nml_var_default(soil_vwc_var, 'soil_vwc_var', print_info, gennml)
        call set_nml_var_default(swe_var, 'swe_var', print_info, gennml)
        call set_nml_var_default(snowh_var, 'snowh_var', print_info, gennml)
        call set_nml_var_default(soil_deept_var, 'soil_deept_var', print_info, gennml)
        call set_nml_var_default(vegtype_var, 'vegtype_var', print_info, gennml)
        call set_nml_var_default(vegfrac_var, 'vegfrac_var', print_info, gennml)
        call set_nml_var_default(vegfracmax_var, 'vegfracmax_var', print_info, gennml)
        call set_nml_var_default(albedo_var, 'albedo_var', print_info, gennml)
        call set_nml_var_default(lai_var, 'lai_var', print_info, gennml)
        call set_nml_var_default(canwat_var, 'canwat_var', print_info, gennml)
        call set_nml_var_default(linear_mask_var, 'linear_mask_var', print_info, gennml)
        call set_nml_var_default(nsq_calibration_var, 'nsq_calibration_var', print_info, gennml)
        call set_nml_var_default(sinalpha_var, 'sinalpha_var', print_info, gennml)
        call set_nml_var_default(cosalpha_var, 'cosalpha_var', print_info, gennml)

        call set_nml_var_default(svf_var, 'svf_var', print_info, gennml)
        call set_nml_var_default(hlm_var, 'hlm_var', print_info, gennml)
        call set_nml_var_default(slope_var, 'slope_var', print_info, gennml)
        call set_nml_var_default(slope_angle_var, 'slope_angle_var', print_info, gennml)
        call set_nml_var_default(aspect_angle_var, 'aspect_angle_var', print_info, gennml)
        call set_nml_var_default(shd_var, 'shd_var', print_info, gennml)

        call set_nml_var_default(dz_levels, 'dz_levels', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=domain)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'domain' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        sleve(n_indx) = sleve(1)
        use_agl_height(n_indx) = use_agl_height(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=domain)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'domain' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif


        call set_nml_var(domain_options%dx, dx(n_indx), 'dx', dx(1))
        call set_nml_var(domain_options%nz, nz(n_indx), 'nz', nz(1))
        call set_nml_var(domain_options%longitude_system, longitude_system(n_indx), 'longitude_system', longitude_system(1))
        call set_nml_var(domain_options%init_conditions_file, init_conditions_file(n_indx), 'init_conditions_file', init_conditions_file(1))
        call set_nml_var(domain_options%flat_z_height, flat_z_height(n_indx), 'flat_z_height', flat_z_height(1))
        call set_nml_var(domain_options%sleve, sleve(n_indx), 'sleve', sleve(1))
        call set_nml_var(domain_options%terrain_smooth_windowsize, terrain_smooth_windowsize(n_indx), 'terrain_smooth_windowsize', terrain_smooth_windowsize(1))
        call set_nml_var(domain_options%terrain_smooth_cycles, terrain_smooth_cycles(n_indx), 'terrain_smooth_cycles', terrain_smooth_cycles(1))
        call set_nml_var(domain_options%decay_rate_L_topo, decay_rate_L_topo(n_indx), 'decay_rate_L_topo', decay_rate_L_topo(1))
        call set_nml_var(domain_options%decay_rate_S_topo, decay_rate_S_topo(n_indx), 'decay_rate_S_topo', decay_rate_S_topo(1))
        call set_nml_var(domain_options%sleve_n, sleve_n(n_indx), 'sleve_n', sleve_n(1))
        call set_nml_var(domain_options%use_agl_height, use_agl_height(n_indx), 'use_agl_height', use_agl_height(1))
        call set_nml_var(domain_options%agl_cap, agl_cap(n_indx), 'agl_cap', agl_cap(1))

        ! NOTE: hgt_hi has to be the first of the variables read
        call set_nml_var(domain_options%hgt_hi, hgt_hi(n_indx), 'hgt_hi',domain_options, hgt_hi(1))
        call set_nml_var(domain_options%landvar, landvar(n_indx), 'landvar',domain_options, landvar(1))
        call set_nml_var(domain_options%lakedepthvar, lakedepthvar(n_indx), 'lakedepthvar',domain_options, lakedepthvar(1))
        call set_nml_var(domain_options%lat_hi, lat_hi(n_indx), 'lat_hi',domain_options, lat_hi(1))
        call set_nml_var(domain_options%lon_hi, lon_hi(n_indx), 'lon_hi',domain_options, lon_hi(1))
        call set_nml_var(domain_options%ulat_hi, ulat_hi(n_indx), 'ulat_hi',domain_options, ulat_hi(1))
        call set_nml_var(domain_options%ulon_hi, ulon_hi(n_indx), 'ulon_hi',domain_options, ulon_hi(1))
        call set_nml_var(domain_options%vlat_hi, vlat_hi(n_indx), 'vlat_hi',domain_options, vlat_hi(1))
        call set_nml_var(domain_options%vlon_hi, vlon_hi(n_indx), 'vlon_hi',domain_options, vlon_hi(1))
        call set_nml_var(domain_options%soiltype_var, soiltype_var(n_indx), 'soiltype_var',domain_options, soiltype_var(1))
        call set_nml_var(domain_options%soil_t_var, soil_t_var(n_indx), 'soil_t_var',domain_options, soil_t_var(1))
        call set_nml_var(domain_options%soil_vwc_var, soil_vwc_var(n_indx), 'soil_vwc_var',domain_options, soil_vwc_var(1))
        call set_nml_var(domain_options%swe_var, swe_var(n_indx), 'swe_var',domain_options, swe_var(1))
        call set_nml_var(domain_options%snowh_var, snowh_var(n_indx), 'snowh_var',domain_options, snowh_var(1))
        call set_nml_var(domain_options%soil_deept_var, soil_deept_var(n_indx), 'soil_deept_var',domain_options, soil_deept_var(1))
        call set_nml_var(domain_options%vegtype_var, vegtype_var(n_indx), 'vegtype_var',domain_options, vegtype_var(1))
        call set_nml_var(domain_options%vegfrac_var, vegfrac_var(n_indx), 'vegfrac_var',domain_options, vegfrac_var(1))
        call set_nml_var(domain_options%vegfracmax_var, vegfracmax_var(n_indx), 'vegfracmax_var',domain_options, vegfracmax_var(1))
        call set_nml_var(domain_options%albedo_var, albedo_var(n_indx), 'albedo_var',domain_options, albedo_var(1))
        call set_nml_var(domain_options%lai_var, lai_var(n_indx), 'lai_var',domain_options, lai_var(1))
        call set_nml_var(domain_options%canwat_var, canwat_var(n_indx), 'canwat_var',domain_options, canwat_var(1))
        call set_nml_var(domain_options%linear_mask_var, linear_mask_var(n_indx), 'linear_mask_var',domain_options, linear_mask_var(1))
        call set_nml_var(domain_options%nsq_calibration_var, nsq_calibration_var(n_indx), 'nsq_calibration_var',domain_options, nsq_calibration_var(1))
        call set_nml_var(domain_options%sinalpha_var, sinalpha_var(n_indx), 'sinalpha_var',domain_options, sinalpha_var(1))
        call set_nml_var(domain_options%cosalpha_var, cosalpha_var(n_indx), 'cosalpha_var',domain_options, cosalpha_var(1))

        call set_nml_var(domain_options%svf_var, svf_var(n_indx), 'svf_var',domain_options, svf_var(1))
        call set_nml_var(domain_options%hlm_var, hlm_var(n_indx), 'hlm_var',domain_options, hlm_var(1))
        call set_nml_var(domain_options%slope_var, slope_var(n_indx), 'slope_var',domain_options, slope_var(1))
        call set_nml_var(domain_options%slope_angle_var, slope_angle_var(n_indx), 'slope_angle_var',domain_options, slope_angle_var(1))
        call set_nml_var(domain_options%aspect_angle_var, aspect_angle_var(n_indx), 'aspect_angle_var',domain_options, aspect_angle_var(1))
        call set_nml_var(domain_options%shd_var, shd_var(n_indx), 'shd_var',domain_options, shd_var(1))

        ! dx and init_conditions_file are required, check that they are set
        call require_var(domain_options%init_conditions_file, "Initial Conditions file")
        if (domain_options%dx <= 0) then
            if (STD_OUT_PE) write(*,*) "  ERROR: dx must be specified in namelist"
            if (STD_OUT_PE) write(*,*) "  ERROR: dx not specified for domain: ", n_indx
            stop
        endif

        allocate(domain_options%dz_levels(domain_options%nz))
        
        ! if nz wasn't specified in the namelist, we assume a HUGE number of levels
        ! so now we have to figure out what the actual number of levels read was
        if (ALL(dz_levels(:,n_indx)==kREAL_NO_VAL)) then
            if (STD_OUT_PE) write(*,*) "  WARNING: dz_levels not specified in namelist for domain: ", n_indx
            if (n_indx == 1) then
                stop 'dz_levels must be specified in namelist for at least the first domain'
            else
                if (STD_OUT_PE) write(*,*) "  WARNING: Copying over the values from the first domain"
            endif
            dz_levels(:,n_indx) = dz_levels(:,1)
        endif
        if (domain_options%nz==MAXLEVELS) then
            do this_level=1,MAXLEVELS-1
                if (dz_levels(this_level+1,n_indx)<=0) then
                    domain_options%nz=this_level
                    exit
                endif
            end do
            domain_options%nz=this_level
        endif

        call set_nml_var(domain_options%dz_levels(1:domain_options%nz), dz_levels(1:domain_options%nz,n_indx), 'dz_levels')

        call require_var(domain_options%lat_hi, "High-res Lat")
        call require_var(domain_options%lon_hi, "High-res Lon")
        call require_var(domain_options%hgt_hi, "High-res HGT")

    end subroutine domain_namelist


    !> -------------------------------
    !! Initialize the microphysics options
    !!
    !! Reads the mp_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine mp_parameters_namelist(mp_filename, mp_options, use_mp_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)    :: mp_filename
        type(mp_options_type), intent(inout) :: mp_options
        logical, intent(in) :: use_mp_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        logical :: print_info, gennml
        integer :: name_unit, rc

        real, dimension(kMAX_NESTS)    :: Nt_c, TNO, am_s, rho_g, av_s, bv_s, fv_s, av_g, bv_g, av_i, &
                                          Ef_si, Ef_rs, Ef_rg, Ef_ri, C_cubes, C_sqrd, mu_r, t_adjust
        logical, dimension(kMAX_NESTS) :: Ef_rw_l, EF_sw_l
        integer :: top_mp_level(kMAX_NESTS)
        real    :: update_interval_mp(kMAX_NESTS)

        namelist /mp_parameters/ Nt_c, TNO, am_s, rho_g, av_s, bv_s, fv_s, av_g, bv_g, av_i,    &   ! thompson microphysics parameters
                                Ef_si, Ef_rs, Ef_rg, Ef_ri,                                     &   ! thompson microphysics parameters
                                C_cubes, C_sqrd, mu_r, Ef_rw_l, Ef_sw_l, t_adjust,              &   ! thompson microphysics parameters
                                top_mp_level, update_interval_mp

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(Nt_c, 'Nt_c', print_info, gennml)
        call set_nml_var_default(TNO, 'TNO', print_info, gennml)
        call set_nml_var_default(am_s, 'am_s', print_info, gennml)
        call set_nml_var_default(rho_g, 'rho_g', print_info, gennml)
        call set_nml_var_default(av_s, 'av_s', print_info, gennml)
        call set_nml_var_default(bv_s, 'bv_s', print_info, gennml)
        call set_nml_var_default(fv_s, 'fv_s', print_info, gennml)
        call set_nml_var_default(av_g, 'av_g', print_info, gennml)
        call set_nml_var_default(bv_g, 'bv_g', print_info, gennml)
        call set_nml_var_default(av_i, 'av_i', print_info, gennml)
        call set_nml_var_default(Ef_si, 'Ef_si', print_info, gennml)
        call set_nml_var_default(Ef_rs, 'Ef_rs', print_info, gennml)
        call set_nml_var_default(Ef_rg, 'Ef_rg', print_info, gennml)
        call set_nml_var_default(Ef_ri, 'Ef_ri', print_info, gennml)
        call set_nml_var_default(C_cubes, 'C_cubes', print_info, gennml)
        call set_nml_var_default(C_sqrd, 'C_sqrd', print_info, gennml)
        call set_nml_var_default(mu_r, 'mu_r', print_info, gennml)
        call set_nml_var_default(t_adjust, 't_adjust', print_info, gennml)
        call set_nml_var_default(Ef_rw_l, 'Ef_rw_l', print_info, gennml)
        call set_nml_var_default(Ef_sw_l, 'Ef_sw_l', print_info, gennml)
        call set_nml_var_default(top_mp_level, 'top_mp_level', print_info, gennml)
        call set_nml_var_default(update_interval_mp, 'update_interval_mp', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read in the namelist
        if (use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit,iostat=rc, nml=mp_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'mp_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        Ef_rw_l(n_indx) = Ef_rw_l(1)
        Ef_sw_l(n_indx) = Ef_sw_l(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        ! read in the namelist
        if (use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit,iostat=rc, nml=mp_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'mp_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(mp_options%Nt_c, Nt_c(n_indx), 'Nt_c', Nt_c(1))
        call set_nml_var(mp_options%TNO, TNO(n_indx), 'TNO', TNO(1))
        call set_nml_var(mp_options%am_s, am_s(n_indx), 'am_s', am_s(1))
        call set_nml_var(mp_options%rho_g, rho_g(n_indx), 'rho_g', rho_g(1))
        call set_nml_var(mp_options%av_s, av_s(n_indx), 'av_s', av_s(1))
        call set_nml_var(mp_options%bv_s, bv_s(n_indx), 'bv_s', bv_s(1))
        call set_nml_var(mp_options%fv_s, fv_s(n_indx), 'fv_s', fv_s(1))

        call set_nml_var(mp_options%av_g, av_g(n_indx), 'av_g', av_g(1))
        call set_nml_var(mp_options%bv_g, bv_g(n_indx), 'bv_g', bv_g(1))
        call set_nml_var(mp_options%av_i, av_i(n_indx), 'av_i', av_i(1))
        call set_nml_var(mp_options%Ef_si, Ef_si(n_indx), 'Ef_si', Ef_si(1))
        call set_nml_var(mp_options%Ef_rs, Ef_rs(n_indx), 'Ef_rs', Ef_rs(1))
        call set_nml_var(mp_options%Ef_rg, Ef_rg(n_indx), 'Ef_rg', Ef_rg(1))
        call set_nml_var(mp_options%Ef_ri, Ef_ri(n_indx), 'Ef_ri', Ef_ri(1))
        call set_nml_var(mp_options%mu_r, mu_r(n_indx), 'mu_r', mu_r(1))
        call set_nml_var(mp_options%t_adjust, t_adjust(n_indx), 't_adjust', t_adjust(1))
        call set_nml_var(mp_options%C_cubes, C_cubes(n_indx), 'C_cubes', C_cubes(1))
        call set_nml_var(mp_options%C_sqrd, C_sqrd(n_indx), 'C_sqrd', C_sqrd(1))
        call set_nml_var(mp_options%Ef_rw_l, Ef_rw_l(n_indx), 'Ef_rw_l', Ef_rw_l(1))
        call set_nml_var(mp_options%Ef_sw_l, Ef_sw_l(n_indx), 'Ef_sw_l', Ef_sw_l(1))

        call set_nml_var(mp_options%update_interval, update_interval_mp(n_indx), 'update_interval_mp', update_interval_mp(1))
        call set_nml_var(mp_options%top_mp_level, top_mp_level(n_indx), 'top_mp_level', top_mp_level(1))

    end subroutine mp_parameters_namelist


    !> -------------------------------
    !! Initialize the Linear Theory options
    !!
    !! Reads the lt_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lt_parameters_namelist(filename, lt_options, use_lt_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(lt_options_type), intent(inout) :: lt_options
        logical, intent(in) :: use_lt_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        integer :: vert_smooth(kMAX_NESTS)
        logical :: variable_N(kMAX_NESTS)           ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: smooth_nsq(kMAX_NESTS)               ! Smooth the Calculated N^2 over vert_smooth vertical levels
        integer :: buffer(kMAX_NESTS)                   ! number of grid cells to buffer around the domain MUST be >=1
        integer :: stability_window_size(kMAX_NESTS)    ! window to average nsq over
        real    :: max_stability(kMAX_NESTS)            ! limits on the calculated Brunt Vaisala Frequency
        real    :: min_stability(kMAX_NESTS)            ! these may need to be a little narrower.
        real    :: linear_contribution(kMAX_NESTS)      ! multiplier on uhat,vhat before adding to u,v
        real    :: linear_update_fraction(kMAX_NESTS)   ! controls the rate at which the linearfield updates (should be calculated as f(in_dt))

        real    :: N_squared(kMAX_NESTS)                ! static Brunt Vaisala Frequency (N^2) to use
        logical :: remove_lowres_linear(kMAX_NESTS)     ! attempt to remove the linear mountain wave from the forcing low res model
        real    :: rm_N_squared(kMAX_NESTS)             ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        real    :: rm_linear_contribution(kMAX_NESTS)   ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        logical :: spatial_linear_fields(kMAX_NESTS)    ! use a spatially varying linear wind perturbation
        logical :: linear_mask(kMAX_NESTS)              ! use a spatial mask for the linear wind field
        logical :: nsq_calibration(kMAX_NESTS)          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! Look up table generation parameters
        real, dimension(kMAX_NESTS)    :: dirmax, dirmin
        real, dimension(kMAX_NESTS)    :: spdmax, spdmin
        real, dimension(kMAX_NESTS)    :: nsqmax, nsqmin
        integer, dimension(kMAX_NESTS) :: n_dir_values, n_nsq_values, n_spd_values
        real, dimension(kMAX_NESTS)    :: minimum_layer_size       ! Minimum vertical step to permit when computing LUT.
                                            ! If model layers are thicker, substepping will be used.

        ! parameters to control reading from or writing an LUT file
        logical, dimension(kMAX_NESTS) :: read_LUT, write_LUT
        character(len=kMAX_FILE_LENGTH), dimension(kMAX_NESTS) :: u_LUT_Filename, v_LUT_Filename, LUT_Filename
        logical :: overwrite_lt_lut(kMAX_NESTS)

        ! define the namelist
        namelist /lt_parameters/ variable_N, smooth_nsq, buffer, stability_window_size, max_stability, min_stability, &
                                 linear_contribution, linear_update_fraction, N_squared, vert_smooth, &
                                 remove_lowres_linear, rm_N_squared, rm_linear_contribution, &
                                 spatial_linear_fields, linear_mask, nsq_calibration, minimum_layer_size, &
                                 dirmax, dirmin, spdmax, spdmin, nsqmax, nsqmin, n_dir_values, n_nsq_values, n_spd_values, &
                                 read_LUT, write_LUT, u_LUT_Filename, v_LUT_Filename, overwrite_lt_lut, LUT_Filename


        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(variable_N, 'variable_N', print_info, gennml)
        call set_nml_var_default(smooth_nsq, 'smooth_nsq', print_info, gennml)
        call set_nml_var_default(buffer, 'buffer', print_info, gennml)
        call set_nml_var_default(stability_window_size, 'stability_window_size', print_info, gennml)
        call set_nml_var_default(max_stability, 'max_stability', print_info, gennml)
        call set_nml_var_default(min_stability, 'min_stability', print_info, gennml)
        call set_nml_var_default(vert_smooth, 'vert_smooth', print_info, gennml)
        call set_nml_var_default(N_squared, 'N_squared', print_info, gennml)
        call set_nml_var_default(linear_contribution, 'linear_contribution', print_info, gennml)
        call set_nml_var_default(remove_lowres_linear, 'remove_lowres_linear', print_info, gennml)
        call set_nml_var_default(rm_N_squared, 'rm_N_squared', print_info, gennml)
        call set_nml_var_default(rm_linear_contribution, 'rm_linear_contribution', print_info, gennml)
        call set_nml_var_default(linear_update_fraction, 'linear_update_fraction', print_info, gennml)
        call set_nml_var_default(spatial_linear_fields, 'spatial_linear_fields', print_info, gennml)
        call set_nml_var_default(linear_mask, 'linear_mask', print_info, gennml)
        call set_nml_var_default(nsq_calibration, 'nsq_calibration', print_info, gennml)
        call set_nml_var_default(dirmax, 'dirmax', print_info, gennml)
        call set_nml_var_default(dirmin, 'dirmin', print_info, gennml)
        call set_nml_var_default(spdmax, 'spdmax', print_info, gennml)
        call set_nml_var_default(spdmin, 'spdmin', print_info, gennml)
        call set_nml_var_default(nsqmax, 'nsqmax', print_info, gennml)
        call set_nml_var_default(nsqmin, 'nsqmin', print_info, gennml)
        call set_nml_var_default(n_dir_values, 'n_dir_values', print_info, gennml)
        call set_nml_var_default(n_nsq_values, 'n_nsq_values', print_info, gennml)
        call set_nml_var_default(n_spd_values, 'n_spd_values', print_info, gennml)
        call set_nml_var_default(minimum_layer_size, 'minimum_layer_size', print_info, gennml)
        call set_nml_var_default(read_LUT, 'read_LUT', print_info, gennml)
        call set_nml_var_default(write_LUT, 'write_LUT', print_info, gennml)
        call set_nml_var_default(u_LUT_Filename, 'u_LUT_Filename', print_info, gennml)
        call set_nml_var_default(v_LUT_Filename, 'v_LUT_Filename', print_info, gennml)
        call set_nml_var_default(LUT_Filename, 'LUT_Filename', print_info, gennml)
        call set_nml_var_default(overwrite_lt_lut, 'overwrite_lt_lut', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_lt_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lt_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'lt_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        variable_N(n_indx) = variable_N(1)
        smooth_nsq(n_indx) = smooth_nsq(1)
        remove_lowres_linear(n_indx) = remove_lowres_linear(1)
        spatial_linear_fields(n_indx) = spatial_linear_fields(1)
        linear_mask(n_indx) = linear_mask(1)
        nsq_calibration(n_indx) = nsq_calibration(1)
        read_LUT(n_indx) = read_LUT(1)
        write_LUT(n_indx) = write_LUT(1)
        overwrite_lt_lut(n_indx) = overwrite_lt_lut(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        ! read the namelist options
        if (use_lt_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lt_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'lt_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(lt_options%variable_N, variable_N(n_indx), 'variable_N', variable_N(1))
        call set_nml_var(lt_options%smooth_nsq, smooth_nsq(n_indx), 'smooth_nsq', smooth_nsq(1))
        call set_nml_var(lt_options%buffer, buffer(n_indx), 'buffer', buffer(1))
        call set_nml_var(lt_options%stability_window_size, stability_window_size(n_indx), 'stability_window_size', stability_window_size(1))
        call set_nml_var(lt_options%max_stability, max_stability(n_indx), 'max_stability', max_stability(1))
        call set_nml_var(lt_options%min_stability, min_stability(n_indx), 'min_stability', min_stability(1))
        call set_nml_var(lt_options%vert_smooth, vert_smooth(n_indx), 'vert_smooth', vert_smooth(1))
        call set_nml_var(lt_options%N_squared, N_squared(n_indx), 'N_squared', N_squared(1))
        call set_nml_var(lt_options%linear_contribution, linear_contribution(n_indx), 'linear_contribution', linear_contribution(1))
        call set_nml_var(lt_options%remove_lowres_linear, remove_lowres_linear(n_indx), 'remove_lowres_linear')
        call set_nml_var(lt_options%rm_N_squared, rm_N_squared(n_indx), 'rm_N_squared', rm_N_squared(1))
        call set_nml_var(lt_options%rm_linear_contribution, rm_linear_contribution(n_indx), 'rm_linear_contribution', rm_linear_contribution(1))
        call set_nml_var(lt_options%linear_update_fraction, linear_update_fraction(n_indx), 'linear_update_fraction', linear_update_fraction(1))
        call set_nml_var(lt_options%spatial_linear_fields, spatial_linear_fields(n_indx), 'spatial_linear_fields', spatial_linear_fields(1))
        call set_nml_var(lt_options%linear_mask, linear_mask(n_indx), 'linear_mask', linear_mask(1))
        call set_nml_var(lt_options%nsq_calibration, nsq_calibration(n_indx), 'nsq_calibration', nsq_calibration(1))
        call set_nml_var(lt_options%dirmax, dirmax(n_indx), 'dirmax', dirmax(1))
        call set_nml_var(lt_options%dirmin, dirmin(n_indx), 'dirmin', dirmin(1))
        call set_nml_var(lt_options%spdmax, spdmax(n_indx), 'spdmax', spdmax(1))
        call set_nml_var(lt_options%spdmin, spdmin(n_indx), 'spdmin', spdmin(1))
        call set_nml_var(lt_options%nsqmax, nsqmax(n_indx), 'nsqmax', nsqmax(1))
        call set_nml_var(lt_options%nsqmin, nsqmin(n_indx), 'nsqmin', nsqmin(1))
        call set_nml_var(lt_options%n_dir_values, n_dir_values(n_indx), 'n_dir_values', n_dir_values(1))
        call set_nml_var(lt_options%n_nsq_values, n_nsq_values(n_indx), 'n_nsq_values', n_nsq_values(1))
        call set_nml_var(lt_options%n_spd_values, n_spd_values(n_indx), 'n_spd_values', n_spd_values(1))
        call set_nml_var(lt_options%minimum_layer_size, minimum_layer_size(n_indx), 'minimum_layer_size', minimum_layer_size(1))
        call set_nml_var(lt_options%read_LUT, read_LUT(n_indx), 'read_LUT', read_LUT(1))
        call set_nml_var(lt_options%write_LUT, write_LUT(n_indx), 'write_LUT', write_LUT(1))
        call set_nml_var(lt_options%u_LUT_Filename, u_LUT_Filename(n_indx), 'u_LUT_Filename', u_LUT_Filename(1))
        call set_nml_var(lt_options%v_LUT_Filename, v_LUT_Filename(n_indx), 'v_LUT_Filename', v_LUT_Filename(1))
        call set_nml_var(lt_options%overwrite_lt_lut, overwrite_lt_lut(n_indx), 'overwrite_lt_lut', overwrite_lt_lut(1))


    end subroutine lt_parameters_namelist


    !> -------------------------------
    !! Initialize the advection options
    !!
    !! Reads the adv_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine adv_parameters_namelist(filename, adv_options, use_adv_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(adv_options_type), intent(inout) :: adv_options
        logical, intent(in) :: use_adv_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        logical :: boundary_buffer(kMAX_NESTS)   ! apply some smoothing to the x and y boundaries in MPDATA
        logical :: MPDATA_FCT(kMAX_NESTS) ! use the flux corrected transport option in MPDATA
        logical :: advect_density(kMAX_NESTS)
        ! MPDATA order of correction (e.g. 1st=upwind, 2nd=classic, 3rd=better)
        integer, dimension(kMAX_NESTS) :: mpdata_order, flux_corr, h_order, v_order
        
        ! define the namelist
        namelist /adv_parameters/ boundary_buffer, MPDATA_FCT, mpdata_order, flux_corr, h_order, v_order, advect_density

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(boundary_buffer, 'boundary_buffer', print_info, gennml)
        call set_nml_var_default(advect_density, 'advect_density', print_info, gennml)
        call set_nml_var_default(MPDATA_FCT, 'MPDATA_FCT', print_info, gennml)
        call set_nml_var_default(mpdata_order, 'mpdata_order', print_info, gennml)
        call set_nml_var_default(flux_corr, 'flux_corr', print_info, gennml)
        call set_nml_var_default(h_order, 'h_order', print_info, gennml)
        call set_nml_var_default(v_order, 'v_order', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_adv_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=adv_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'adv_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        boundary_buffer(n_indx) = boundary_buffer(1)
        MPDATA_FCT(n_indx) = MPDATA_FCT(1)
        advect_density(n_indx) = advect_density(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        ! read the namelist options
        if (use_adv_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=adv_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'adv_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(adv_options%boundary_buffer, boundary_buffer(n_indx), 'boundary_buffer', boundary_buffer(1))
        call set_nml_var(adv_options%MPDATA_FCT, MPDATA_FCT(n_indx), 'MPDATA_FCT', MPDATA_FCT(1))
        call set_nml_var(adv_options%mpdata_order, mpdata_order(n_indx), 'mpdata_order', mpdata_order(1))
        call set_nml_var(adv_options%flux_corr, flux_corr(n_indx), 'flux_corr', flux_corr(1))
        call set_nml_var(adv_options%h_order, h_order(n_indx), 'h_order', h_order(1))
        call set_nml_var(adv_options%v_order, v_order(n_indx), 'v_order', v_order(1))
        call set_nml_var(adv_options%advect_density, advect_density(n_indx), 'advect_density', advect_density(1))

    end subroutine adv_parameters_namelist


    !> -------------------------------
    !! Initialize the PBL options
    !!
    !! Reads the pbl_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine pbl_parameters_namelist(filename, pbl_options, use_pbl_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(pbl_options_type), intent(inout) :: pbl_options
        logical, intent(in) :: use_pbl_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        integer :: ysu_topdown_pblmix(kMAX_NESTS) ! controls if radiative, top-down mixing in YSU scheme is turned on
        
        ! define the namelist
        namelist /pbl_parameters/ ysu_topdown_pblmix

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(ysu_topdown_pblmix, 'ysu_topdown_pblmix', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_pbl_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=pbl_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'pbl_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(pbl_options%ysu_topdown_pblmix, ysu_topdown_pblmix(n_indx), 'ysu_topdown_pblmix', ysu_topdown_pblmix(1))
    end subroutine pbl_parameters_namelist

    !> -------------------------------
    !! Initialize the surface layer options
    !!
    !! Reads the sfc_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine sfc_parameters_namelist(filename, sfc_options, use_sfc_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(sfc_options_type), intent(inout) :: sfc_options
        logical, intent(in) :: use_sfc_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        integer, dimension(kMAX_NESTS) :: isfflx, scm_force_flux, iz0tlnd, isftcflx
        real    :: sbrlim(kMAX_NESTS)

        logical :: print_info, gennml
        ! define the namelist
        namelist /sfc_parameters/ isfflx, scm_force_flux, iz0tlnd, sbrlim, isftcflx
         
        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(isfflx, 'isfflx', print_info, gennml)
        call set_nml_var_default(scm_force_flux, 'scm_force_flux', print_info, gennml)
        call set_nml_var_default(iz0tlnd, 'iz0tlnd', print_info, gennml)
        call set_nml_var_default(isftcflx, 'isftcflx', print_info, gennml)
        call set_nml_var_default(sbrlim, 'sbrlim', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_sfc_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=sfc_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'sfc_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        
        call set_nml_var(sfc_options%isfflx, isfflx(n_indx), 'isfflx', isfflx(1))
        call set_nml_var(sfc_options%scm_force_flux, scm_force_flux(n_indx), 'scm_force_flux', scm_force_flux(1))
        call set_nml_var(sfc_options%iz0tlnd, iz0tlnd(n_indx), 'iz0tlnd', iz0tlnd(1))
        call set_nml_var(sfc_options%isftcflx, isftcflx(n_indx), 'isftcflx', isftcflx(1))
        call set_nml_var(sfc_options%sbrlim, sbrlim(n_indx), 'sbrlim', sbrlim(1))

    end subroutine sfc_parameters_namelist


    !> -------------------------------
    !! Initialize the convection scheme options
    !!
    !! Reads the cu_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine cu_parameters_namelist(filename, cu_options, use_cu_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(cu_options_type), intent(inout) :: cu_options
        logical, intent(in) :: use_cu_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        real, dimension(kMAX_NESTS) :: tendency_fraction, tend_qv_fraction, tend_qc_fraction, tend_th_fraction, tend_qi_fraction, &
                                       stochastic_cu


        ! define the namelist
        namelist /cu_parameters/ tendency_fraction, tend_qv_fraction, tend_qc_fraction, tend_th_fraction, tend_qi_fraction, &
                                 stochastic_cu

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(stochastic_cu, 'stochastic_cu', print_info, gennml)
        call set_nml_var_default(tendency_fraction, 'tendency_fraction', print_info, gennml)
        call set_nml_var_default(tend_qv_fraction, 'tend_qv_fraction', print_info, gennml)
        call set_nml_var_default(tend_qc_fraction, 'tend_qc_fraction', print_info, gennml)
        call set_nml_var_default(tend_th_fraction, 'tend_th_fraction', print_info, gennml)
        call set_nml_var_default(tend_qi_fraction, 'tend_qi_fraction', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_cu_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=cu_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'cu_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        ! if not set separately, default to the global tendency setting
        if (tend_qv_fraction(n_indx) < 0) tend_qv_fraction(n_indx) = tendency_fraction(n_indx)
        if (tend_qc_fraction(n_indx) < 0) tend_qc_fraction(n_indx) = tendency_fraction(n_indx)
        if (tend_th_fraction(n_indx) < 0) tend_th_fraction(n_indx) = tendency_fraction(n_indx)
        if (tend_qi_fraction(n_indx) < 0) tend_qi_fraction(n_indx) = tendency_fraction(n_indx)

        ! store everything in the cu_options structure
        call set_nml_var(cu_options%tendency_fraction, tendency_fraction(n_indx), 'tendency_fraction', tendency_fraction(1))
        call set_nml_var(cu_options%tend_qv_fraction, tend_qv_fraction(n_indx), 'tend_qv_fraction', tend_qv_fraction(1))
        call set_nml_var(cu_options%tend_qc_fraction, tend_qc_fraction(n_indx), 'tend_qc_fraction', tend_qc_fraction(1))
        call set_nml_var(cu_options%tend_th_fraction, tend_th_fraction(n_indx), 'tend_th_fraction', tend_th_fraction(1))
        call set_nml_var(cu_options%tend_qi_fraction, tend_qi_fraction(n_indx), 'tend_qi_fraction', tend_qi_fraction(1))
        call set_nml_var(cu_options%stochastic_cu, stochastic_cu(n_indx), 'stochastic_cu', stochastic_cu(1))

    end subroutine cu_parameters_namelist


    !> -------------------------------
    !! Initialize the land surface model options
    !!
    !! Reads the lsm_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lsm_parameters_namelist(filename, lsm_options, use_lsm_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(lsm_options_type), intent(inout) :: lsm_options
        logical, intent(in) :: use_lsm_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        character(len=kMAX_NAME_LENGTH) :: LU_Categories(kMAX_NESTS) ! Category definitions (e.g. USGS, MODIFIED_IGBP_MODIS_NOAH)
        real    :: max_swe(kMAX_NESTS)
        real    :: snow_den_const(kMAX_NESTS)                    ! variable for converting snow height into SWE or visa versa when input data is incomplete 

        logical :: monthly_vegfrac(kMAX_NESTS)                   ! read in 12 months of vegfrac data
        logical :: monthly_albedo(kMAX_NESTS)                    ! same for albedo (requires vegfrac be monthly)
        real :: update_interval_lsm(kMAX_NESTS)                  ! minimum number of seconds between LSM updates
        integer :: urban_category(kMAX_NESTS)                    ! index that defines the urban category in LU_Categories
        integer :: ice_category(kMAX_NESTS)                      ! index that defines the ice category in LU_Categories
        integer :: water_category(kMAX_NESTS)                    ! index that defines the water category in LU_Categories
        integer :: sf_urban_phys(kMAX_NESTS)
        integer :: num_soil_layers(kMAX_NESTS)
        real    :: nmp_soiltstep(kMAX_NESTS)
        integer, dimension(kMAX_NESTS) :: nmp_dveg, nmp_opt_crs, nmp_opt_sfc, nmp_opt_btr, nmp_opt_run, nmp_opt_frz, nmp_opt_inf, nmp_opt_rad, nmp_opt_alb, nmp_opt_snf, nmp_opt_tbot, nmp_opt_stc, nmp_opt_gla, nmp_opt_rsf, nmp_opt_soil, nmp_opt_pedo, nmp_opt_crop, nmp_opt_irr, nmp_opt_irrm, nmp_opt_tdrn, noahmp_output

        integer :: lake_category(kMAX_NESTS)                    ! index that defines the lake category in (some) LU_Categories

        ! define the namelist
        namelist /lsm_parameters/ LU_Categories, update_interval_lsm, &
                                  urban_category, ice_category, water_category, lake_category, snow_den_const,&
                                  monthly_vegfrac, monthly_albedo, max_swe,  nmp_dveg,   &
                                  nmp_opt_crs, nmp_opt_sfc, nmp_opt_btr, nmp_opt_run, nmp_opt_frz, &
                                  nmp_opt_inf, nmp_opt_rad, nmp_opt_alb, nmp_opt_snf, nmp_opt_tbot,           &
                                  nmp_opt_stc, nmp_opt_gla, nmp_opt_rsf, nmp_opt_soil, nmp_opt_pedo,          &
                                  nmp_opt_crop, nmp_opt_irr, nmp_opt_irrm, nmp_opt_tdrn, nmp_soiltstep,       &
                                  sf_urban_phys, noahmp_output, num_soil_layers !! MJ added

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(LU_Categories, 'LU_Categories', print_info, gennml)
        call set_nml_var_default(update_interval_lsm, 'update_interval_lsm', print_info, gennml)
        call set_nml_var_default(monthly_vegfrac, 'monthly_vegfrac', print_info, gennml)
        call set_nml_var_default(num_soil_layers, 'num_soil_layers', print_info, gennml)

        call set_nml_var_default(monthly_albedo, 'monthly_albedo', print_info, gennml)
        call set_nml_var_default(urban_category, 'urban_category', print_info, gennml)
        call set_nml_var_default(ice_category, 'ice_category', print_info, gennml)
        call set_nml_var_default(water_category, 'water_category', print_info, gennml)
        call set_nml_var_default(lake_category, 'lake_category', print_info, gennml)
        call set_nml_var_default(snow_den_const, 'snow_den_const', print_info, gennml)
        call set_nml_var_default(max_swe, 'max_swe', print_info, gennml)
        call set_nml_var_default(sf_urban_phys, 'sf_urban_phys', print_info, gennml)
        call set_nml_var_default(nmp_dveg, 'nmp_dveg', print_info, gennml)
        call set_nml_var_default(nmp_opt_crs, 'nmp_opt_crs', print_info, gennml)
        call set_nml_var_default(nmp_opt_sfc, 'nmp_opt_sfc', print_info, gennml)
        call set_nml_var_default(nmp_opt_btr, 'nmp_opt_btr', print_info, gennml)
        call set_nml_var_default(nmp_opt_run, 'nmp_opt_run', print_info, gennml)
        call set_nml_var_default(nmp_opt_frz, 'nmp_opt_frz', print_info, gennml)
        call set_nml_var_default(nmp_opt_inf, 'nmp_opt_inf', print_info, gennml)
        call set_nml_var_default(nmp_opt_rad, 'nmp_opt_rad', print_info, gennml)
        call set_nml_var_default(nmp_opt_alb, 'nmp_opt_alb', print_info, gennml)
        call set_nml_var_default(nmp_opt_snf, 'nmp_opt_snf', print_info, gennml)
        call set_nml_var_default(nmp_opt_tbot, 'nmp_opt_tbot', print_info, gennml)
        call set_nml_var_default(nmp_opt_stc, 'nmp_opt_stc', print_info, gennml)
        call set_nml_var_default(nmp_opt_gla, 'nmp_opt_gla', print_info, gennml)
        call set_nml_var_default(nmp_opt_rsf, 'nmp_opt_rsf', print_info, gennml)
        call set_nml_var_default(nmp_opt_soil, 'nmp_opt_soil', print_info, gennml)

        call set_nml_var_default(nmp_opt_pedo, 'nmp_opt_pedo', print_info, gennml)
        call set_nml_var_default(nmp_opt_crop, 'nmp_opt_crop', print_info, gennml)
        call set_nml_var_default(nmp_opt_irr, 'nmp_opt_irr', print_info, gennml)
        call set_nml_var_default(nmp_opt_irrm, 'nmp_opt_irrm', print_info, gennml)
        call set_nml_var_default(nmp_opt_tdrn, 'nmp_opt_tdrn', print_info, gennml)
        call set_nml_var_default(nmp_soiltstep, 'nmp_soiltstep', print_info, gennml)
        call set_nml_var_default(noahmp_output, 'noahmp_output', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_lsm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lsm_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'lsm_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        monthly_albedo(n_indx) = monthly_albedo(1)
        monthly_vegfrac(n_indx) = monthly_vegfrac(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        ! read the namelist options
        if (use_lsm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lsm_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'lsm_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_default_LU_categories(LU_Categories(n_indx), urban_category(n_indx), ice_category(n_indx), water_category(n_indx), lake_category(n_indx))

        call set_nml_var(lsm_options%LU_Categories, LU_Categories(n_indx), 'LU_Categories', LU_Categories(1))
        call set_nml_var(lsm_options%update_interval, update_interval_lsm(n_indx), 'update_interval_lsm', update_interval_lsm(1))
        call set_nml_var(lsm_options%monthly_vegfrac, monthly_vegfrac(n_indx), 'monthly_vegfrac', monthly_vegfrac(1))
        call set_nml_var(lsm_options%num_soil_layers, num_soil_layers(n_indx), 'num_soil_layers', num_soil_layers(1))

        call set_nml_var(lsm_options%monthly_albedo, monthly_albedo(n_indx), 'monthly_albedo', monthly_albedo(1))
        call set_nml_var(lsm_options%urban_category, urban_category(n_indx), 'urban_category', urban_category(1))
        call set_nml_var(lsm_options%ice_category, ice_category(n_indx), 'ice_category', ice_category(1))
        call set_nml_var(lsm_options%water_category, water_category(n_indx), 'water_category', water_category(1))
        call set_nml_var(lsm_options%lake_category, lake_category(n_indx), 'lake_category', lake_category(1))
        call set_nml_var(lsm_options%snow_den_const, snow_den_const(n_indx), 'snow_den_const', snow_den_const(1))
        call set_nml_var(lsm_options%max_swe, max_swe(n_indx), 'max_swe', max_swe(1))
        call set_nml_var(lsm_options%sf_urban_phys, sf_urban_phys(n_indx), 'sf_urban_phys', sf_urban_phys(1))
        call set_nml_var(lsm_options%nmp_dveg, nmp_dveg(n_indx), 'nmp_dveg', nmp_dveg(1))

        call set_nml_var(lsm_options%nmp_opt_crs, nmp_opt_crs(n_indx), 'nmp_opt_crs', nmp_opt_crs(1))
        call set_nml_var(lsm_options%nmp_opt_sfc, nmp_opt_sfc(n_indx), 'nmp_opt_sfc', nmp_opt_sfc(1))
        call set_nml_var(lsm_options%nmp_opt_btr, nmp_opt_btr(n_indx), 'nmp_opt_btr', nmp_opt_btr(1))
        call set_nml_var(lsm_options%nmp_opt_run, nmp_opt_run(n_indx), 'nmp_opt_run', nmp_opt_run(1))
        call set_nml_var(lsm_options%nmp_opt_frz, nmp_opt_frz(n_indx), 'nmp_opt_frz', nmp_opt_frz(1))
        call set_nml_var(lsm_options%nmp_opt_inf, nmp_opt_inf(n_indx), 'nmp_opt_inf', nmp_opt_inf(1))
        call set_nml_var(lsm_options%nmp_opt_rad, nmp_opt_rad(n_indx), 'nmp_opt_rad', nmp_opt_rad(1))
        call set_nml_var(lsm_options%nmp_opt_alb, nmp_opt_alb(n_indx), 'nmp_opt_alb', nmp_opt_alb(1))
        call set_nml_var(lsm_options%nmp_opt_snf, nmp_opt_snf(n_indx), 'nmp_opt_snf', nmp_opt_snf(1))
        call set_nml_var(lsm_options%nmp_opt_tbot, nmp_opt_tbot(n_indx), 'nmp_opt_tbot', nmp_opt_tbot(1))
        call set_nml_var(lsm_options%nmp_opt_stc, nmp_opt_stc(n_indx), 'nmp_opt_stc', nmp_opt_stc(1))
        call set_nml_var(lsm_options%nmp_opt_gla, nmp_opt_gla(n_indx), 'nmp_opt_gla', nmp_opt_gla(1))
        call set_nml_var(lsm_options%nmp_opt_rsf, nmp_opt_rsf(n_indx), 'nmp_opt_rsf', nmp_opt_rsf(1))
        call set_nml_var(lsm_options%nmp_opt_soil, nmp_opt_soil(n_indx), 'nmp_opt_soil', nmp_opt_soil(1))
        call set_nml_var(lsm_options%nmp_opt_pedo, nmp_opt_pedo(n_indx), 'nmp_opt_pedo', nmp_opt_pedo(1))
        call set_nml_var(lsm_options%nmp_opt_crop, nmp_opt_crop(n_indx), 'nmp_opt_crop', nmp_opt_crop(1))
        call set_nml_var(lsm_options%nmp_opt_irr, nmp_opt_irr(n_indx), 'nmp_opt_irr', nmp_opt_irr(1))
        call set_nml_var(lsm_options%nmp_opt_irrm, nmp_opt_irrm(n_indx), 'nmp_opt_irrm', nmp_opt_irrm(1))
        call set_nml_var(lsm_options%nmp_opt_tdrn, nmp_opt_tdrn(n_indx), 'nmp_opt_tdrn', nmp_opt_tdrn(1))
        call set_nml_var(lsm_options%nmp_soiltstep, nmp_soiltstep(n_indx), 'nmp_soiltstep', nmp_soiltstep(1))
        call set_nml_var(lsm_options%noahmp_output, noahmp_output(n_indx), 'noahmp_output', noahmp_output(1))

        
    end subroutine lsm_parameters_namelist

    subroutine sm_parameters_namelist(filename, sm_options, use_sm_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(sm_options_type), intent(inout) :: sm_options
        logical, intent(in) :: use_sm_options
        integer, intent(in) :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        integer :: fsm_nsnow_max(kMAX_NESTS)      ! maximum number of snow layers in the FSM2trans snow model
        integer, dimension(kMAX_NESTS) :: fsm_albedo, fsm_canmod, fsm_checks, fsm_condct, fsm_densty, fsm_exchng, &
                   fsm_hydrol, fsm_radsbg, fsm_snfrac, fsm_snolay, fsm_snslid, fsm_sntran, fsm_zoffst
        real, dimension(kMAX_NESTS)    :: fsm_ds_min, fsm_ds_surflay
        logical, dimension(kMAX_NESTS) :: fsm_hn_on, fsm_for_hn

        ! define the namelist
        namelist /sm_parameters/ fsm_nsnow_max, fsm_albedo, fsm_canmod, fsm_checks, fsm_condct, fsm_densty, fsm_exchng, &
                                 fsm_hydrol, fsm_radsbg, fsm_snfrac, fsm_snolay, fsm_snslid, fsm_sntran, fsm_zoffst, &
                                 fsm_ds_min, fsm_ds_surflay, fsm_hn_on, fsm_for_hn


        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(fsm_nsnow_max, 'fsm_nsnow_max', print_info, gennml)
        call set_nml_var_default(fsm_albedo, 'fsm_albedo', print_info, gennml)
        call set_nml_var_default(fsm_canmod, 'fsm_canmod', print_info, gennml)
        call set_nml_var_default(fsm_checks, 'fsm_checks', print_info, gennml)
        call set_nml_var_default(fsm_condct, 'fsm_condct', print_info, gennml)
        call set_nml_var_default(fsm_densty, 'fsm_densty', print_info, gennml)
        call set_nml_var_default(fsm_exchng, 'fsm_exchng', print_info, gennml)
        call set_nml_var_default(fsm_hydrol, 'fsm_hydrol', print_info, gennml)
        call set_nml_var_default(fsm_radsbg, 'fsm_radsbg', print_info, gennml)
        call set_nml_var_default(fsm_snfrac, 'fsm_snfrac', print_info, gennml)
        call set_nml_var_default(fsm_snolay, 'fsm_snolay', print_info, gennml)
        call set_nml_var_default(fsm_snslid, 'fsm_snslid', print_info, gennml)
        call set_nml_var_default(fsm_sntran, 'fsm_sntran', print_info, gennml)
        call set_nml_var_default(fsm_zoffst, 'fsm_zoffst', print_info, gennml)
        call set_nml_var_default(fsm_ds_min, 'fsm_ds_min', print_info, gennml)
        call set_nml_var_default(fsm_ds_surflay, 'fsm_ds_surflay', print_info, gennml)
        
        call set_nml_var_default(fsm_hn_on, 'fsm_hn_on', print_info, gennml)
        call set_nml_var_default(fsm_for_hn, 'fsm_for_hn', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_sm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=sm_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'sm_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        fsm_hn_on(n_indx) = fsm_hn_on(1)
        fsm_for_hn(n_indx) = fsm_for_hn(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        ! read the namelist options
        if (use_sm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=sm_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'sm_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(sm_options%fsm_nsnow_max, fsm_nsnow_max(n_indx), 'fsm_nsnow_max', fsm_nsnow_max(1))
        call set_nml_var(sm_options%fsm_albedo, fsm_albedo(n_indx), 'fsm_albedo', fsm_albedo(1))
        call set_nml_var(sm_options%fsm_canmod, fsm_canmod(n_indx), 'fsm_canmod', fsm_canmod(1))
        call set_nml_var(sm_options%fsm_checks, fsm_checks(n_indx), 'fsm_checks', fsm_checks(1))
        call set_nml_var(sm_options%fsm_condct, fsm_condct(n_indx), 'fsm_condct', fsm_condct(1))
        call set_nml_var(sm_options%fsm_densty, fsm_densty(n_indx), 'fsm_densty', fsm_densty(1))
        call set_nml_var(sm_options%fsm_exchng, fsm_exchng(n_indx), 'fsm_exchng', fsm_exchng(1))
        call set_nml_var(sm_options%fsm_hydrol, fsm_hydrol(n_indx), 'fsm_hydrol', fsm_hydrol(1))
        call set_nml_var(sm_options%fsm_radsbg, fsm_radsbg(n_indx), 'fsm_radsbg', fsm_radsbg(1))
        call set_nml_var(sm_options%fsm_snfrac, fsm_snfrac(n_indx), 'fsm_snfrac', fsm_snfrac(1))
        call set_nml_var(sm_options%fsm_snolay, fsm_snolay(n_indx), 'fsm_snolay', fsm_snolay(1))
        call set_nml_var(sm_options%fsm_snslid, fsm_snslid(n_indx), 'fsm_snslid', fsm_snslid(1))
        call set_nml_var(sm_options%fsm_sntran, fsm_sntran(n_indx), 'fsm_sntran', fsm_sntran(1))
        call set_nml_var(sm_options%fsm_zoffst, fsm_zoffst(n_indx), 'fsm_zoffst', fsm_zoffst(1))
        call set_nml_var(sm_options%fsm_ds_min, fsm_ds_min(n_indx), 'fsm_ds_min', fsm_ds_min(1))
        call set_nml_var(sm_options%fsm_ds_surflay, fsm_ds_surflay(n_indx), 'fsm_ds_surflay', fsm_ds_surflay(1))

        call set_nml_var(sm_options%fsm_hn_on, fsm_hn_on(n_indx), 'fsm_hn_on')
        call set_nml_var(sm_options%fsm_for_hn, fsm_for_hn(n_indx), 'fsm_for_hn')
        
    end subroutine sm_parameters_namelist
    !> -------------------------------
    !! Initialize the radiation model options
    !!
    !! Reads the rad_parameters namelist or sets default values
    !! -------------------------------
    subroutine rad_parameters_namelist(filename, rad_options, use_rad_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(rad_options_type), intent(inout) :: rad_options
        logical, intent(in)  :: use_rad_options
        integer, intent(in)  :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        real    :: update_interval_rrtmg(kMAX_NESTS)             ! minimum number of seconds between RRTMG updates
        integer :: icloud(kMAX_NESTS)                            ! how RRTMG interacts with clouds
        integer :: cldovrlp(kMAX_NESTS)                          ! how RRTMG considers cloud overlapping
        logical :: read_ghg(kMAX_NESTS)
        real    :: tzone(kMAX_NESTS) !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe
        ! define the namelist
        namelist /rad_parameters/ update_interval_rrtmg, icloud, read_ghg, cldovrlp, tzone !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(update_interval_rrtmg, 'update_interval_rrtmg', print_info, gennml)
        call set_nml_var_default(icloud, 'icloud', print_info, gennml)
        call set_nml_var_default(cldovrlp, 'cldovrlp', print_info, gennml)
        call set_nml_var_default(read_ghg, 'read_ghg', print_info, gennml)
        call set_nml_var_default(tzone, 'tzone', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (use_rad_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit, iostat=rc, nml=rad_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'rad_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        read_ghg(n_indx) = read_ghg(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        if (use_rad_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit, iostat=rc, nml=rad_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'rad_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(rad_options%update_interval_rrtmg, update_interval_rrtmg(n_indx), 'update_interval_rrtmg', update_interval_rrtmg(1))
        call set_nml_var(rad_options%icloud, icloud(n_indx), 'icloud', icloud(1))
        call set_nml_var(rad_options%cldovrlp, cldovrlp(n_indx), 'cldovrlp', cldovrlp(1))
        call set_nml_var(rad_options%read_ghg, read_ghg(n_indx), 'read_ghg', read_ghg(1))
        call set_nml_var(rad_options%tzone, tzone(n_indx), 'tzone', tzone(1))
        
    end subroutine rad_parameters_namelist
    
    
    subroutine wind_namelist(filename, wind_options, use_wind_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(wind_type), intent(inout) :: wind_options
        logical, intent(in)           :: use_wind_options
        integer, intent(in)           :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml
        
        integer :: name_unit, rc, update_frequency_checked     ! logical unit number for namelist
        logical :: print_info, gennml
        !Define parameters
        integer, dimension(kMAX_NESTS) :: wind_iterations, update_frequency
        logical, dimension(kMAX_NESTS) :: Sx, thermal, wind_only
        real, dimension(kMAX_NESTS)    :: Sx_dmax, Sx_scale_ang, TPI_scale, TPI_dmax, alpha_const, smooth_wind_distance
        
        !Make name-list
        namelist /wind/ Sx, thermal, wind_only, Sx_dmax, Sx_scale_ang, TPI_scale, TPI_dmax, alpha_const, &
                        update_frequency, smooth_wind_distance, wind_iterations
        
        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(Sx, 'Sx', print_info, gennml)
        call set_nml_var_default(thermal, 'thermal', print_info, gennml)
        call set_nml_var_default(wind_only, 'wind_only', print_info, gennml)
        call set_nml_var_default(Sx_dmax, 'Sx_dmax', print_info, gennml)
        call set_nml_var_default(Sx_scale_ang, 'Sx_scale_ang', print_info, gennml)
        call set_nml_var_default(TPI_scale, 'TPI_scale', print_info, gennml)
        call set_nml_var_default(TPI_dmax, 'TPI_dmax', print_info, gennml)
        call set_nml_var_default(alpha_const, 'alpha_const', print_info, gennml)
        call set_nml_var_default(smooth_wind_distance, 'smooth_wind_distance', print_info, gennml)
        call set_nml_var_default(wind_iterations, 'wind_iterations', print_info, gennml)
        call set_nml_var_default(update_frequency, 'update_frequency', print_info, gennml)
        
        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        !Read namelist file
        if (use_wind_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=wind)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'wind' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        Sx(n_indx) = Sx(1)
        thermal(n_indx) = thermal(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        !Read namelist file
        if (use_wind_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=wind)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
                if (STD_OUT_PE) write(*,*) "  Error reading 'wind' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "  --------------------------------"
            endif
        endif

        call set_nml_var(wind_options%Sx, Sx(n_indx), 'Sx', Sx(1))
        call set_nml_var(wind_options%thermal, thermal(n_indx), 'thermal', thermal(1))
        call set_nml_var(wind_options%wind_only, wind_only(n_indx), 'wind_only', wind_only(1))
        call set_nml_var(wind_options%Sx_dmax, Sx_dmax(n_indx), 'Sx_dmax', Sx_dmax(1))
        call set_nml_var(wind_options%TPI_dmax, TPI_dmax(n_indx), 'TPI_dmax', TPI_dmax(1))
        call set_nml_var(wind_options%TPI_scale, TPI_scale(n_indx), 'TPI_scale', TPI_scale(1))
        call set_nml_var(wind_options%Sx_scale_ang, Sx_scale_ang(n_indx), 'Sx_scale_ang', Sx_scale_ang(1))
        call set_nml_var(wind_options%alpha_const, alpha_const(n_indx), 'alpha_const', alpha_const(1))
        call set_nml_var(wind_options%wind_iterations, wind_iterations(n_indx), 'wind_iterations', wind_iterations(1))
        call set_nml_var(wind_options%smooth_wind_distance, smooth_wind_distance(n_indx), 'smooth_wind_distance', smooth_wind_distance(1))
        call set_nml_var(update_frequency_checked, update_frequency(n_indx), 'update_frequency', update_frequency(1))

        call wind_options%update_dt%set(seconds=update_frequency_checked)

    end subroutine wind_namelist


    subroutine time_parameters_namelist(filename, time_options, n_indx, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(time_options_type), intent(inout) :: time_options
        integer, intent(in)           :: n_indx
        logical, intent(in), optional  :: info_only, gen_nml
        
        integer :: name_unit, rc                            ! logical unit number for namelist
        !Define parameters
        real :: cfl_reduction_factor(kMAX_NESTS)    
        logical :: RK3(kMAX_NESTS)

        logical :: print_info, gennml
        
        !Make name-list
        namelist /time_parameters/ cfl_reduction_factor, RK3
        
        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(cfl_reduction_factor, 'cfl_reduction_factor', print_info, gennml)
        call set_nml_var_default(RK3, 'RK3', print_info, gennml)
        
        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        !Read namelist file
        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=time_parameters)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'time_parameters' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif
        ! Copy the first value of logical variables -- this way we can have a user_default value if the value for this nest was not explicitly set
        RK3(n_indx) = RK3(1)
        ! Now read namelist again, -- if the value of the logical option is set in the namelist, it will be set to the user set value again
        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=time_parameters)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            if (STD_OUT_PE) write(*,*) "  Error reading 'time_parameters' namelist"
            if (STD_OUT_PE) write(*,*) "  --------------------------------"
            stop
        endif

        call set_nml_var(time_options%cfl_reduction_factor, cfl_reduction_factor(n_indx), 'cfl_reduction_factor', cfl_reduction_factor(1))
        call set_nml_var(time_options%RK3, RK3(n_indx), 'RK3')

    end subroutine time_parameters_namelist
    
    
    subroutine set_phys_suite(options)
        implicit none
        type(options_t),    intent(inout) :: options
    
        select case (options%general%phys_suite)
            case('HICAR')
                if (options%domain%dx >= 1000 .and. STD_OUT_PE) then
                    write(*,*) '------------------------------------------------'
                    write(*,*) 'WARNING: Setting HICAR namelist options'
                    write(*,*) 'When user has selected dx => 1000.'
                    write(*,*) 'High-resolution wind solver will be used,'
                    write(*,*) 'which may not be appropriate for this resolution'
                    write(*,*) '------------------------------------------------'
                endif
                !Add base HICAR options here
                
                !options%physics%boundarylayer = at some point, add a scale aware / LES turbulence scheme
                
                options%physics%windtype = 4
                options%physics%convection = 0
                options%wind%Sx = .True.
                options%time%RK3 = .True.
                options%domain%use_agl_height = .True.
                options%domain%agl_cap = 800
                options%domain%sleve = .True.
                options%adv%advect_density = .True.

        end select
    
    end subroutine
    
    !> ----------------------------------------------------------------------------
    !!  Read in the name of the boundary condition files from a text file
    !!
    !!  @param      filename        The name of the text file to read
    !!  @param[out] forcing_files   An array to store the filenames in
    !!  @retval     nfiles          The number of files read.
    !!
    !! ----------------------------------------------------------------------------
    function read_forcing_file_names(filename, forcing_files) result(nfiles)
        implicit none
        character(len=*) :: filename
        character(len=kMAX_FILE_LENGTH), dimension(MAX_NUMBER_FILES) :: forcing_files
        integer :: nfiles
        integer :: file_unit
        integer :: i, error
        logical :: first_file_exists, last_file_exists
        character(len=kMAX_FILE_LENGTH) :: temporary_file

        open(unit=io_newunit(file_unit), file=filename)
        i=0
        error=0
        do while (error==0)
            read(file_unit, *, iostat=error) temporary_file
            if (error==0) then
                i=i+1
                forcing_files(i) = temporary_file
            endif
        enddo
        close(file_unit)
        nfiles = i
        ! print out a summary
        if (STD_OUT_PE) write(*,*) "    Boundary conditions files to be used:"
        if (nfiles>10) then
            if (STD_OUT_PE) write(*,*) "      nfiles=", trim(str(nfiles)), ", too many to print."
            if (STD_OUT_PE) write(*,*) "      First file:", trim(forcing_files(1))
            if (STD_OUT_PE) write(*,*) "      Last file: ", trim(forcing_files(nfiles))
        else
            do i=1,nfiles
                if (STD_OUT_PE) write(*,*) "        ",trim(forcing_files(i))
            enddo
        endif

        ! Check that the options file actually exists
        INQUIRE(file=trim(forcing_files(1)), exist=first_file_exists)
        INQUIRE(file=trim(forcing_files(nfiles)), exist=last_file_exists)

        ! if options file does not exist, print an error and quit
        if (.not.first_file_exists .or. .not.last_file_exists) then
            if (.not.first_file_exists .and. STD_OUT_PE) write(*,*) "  The first forcing file does not exist = ", trim(forcing_files(1))
            if (.not.last_file_exists .and. STD_OUT_PE) write(*,*) "  The last forcing file does not exist = ", trim(forcing_files(nfiles))

            ! stop "At least the first or last forcing file does not exist. Only the first and last file are checked, please check the rest."
        endif


    end function read_forcing_file_names


    !> -------------------------------
    !! Sets the default value for each of three land use categories depending on the LU_Categories input
    !!
    !! -------------------------------
    subroutine set_default_LU_categories(LU_Categories, urban_category, ice_category, water_category, lake_category)
        ! if various LU categories were not defined in the namelist (i.e. they == -1) then attempt
        ! to define default values for them based on the LU_Categories variable supplied.
        implicit none
        integer, intent(inout) :: urban_category, ice_category, water_category, lake_category
        character(len=kMAX_NAME_LENGTH), intent(in) :: LU_Categories

        if (trim(LU_Categories)=="MODIFIED_IGBP_MODIS_NOAH") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            if (lake_category==-1) lake_category = 21

        elseif (trim(LU_Categories)=="USGS") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = -1
            if (water_category==-1) water_category = 16
            ! if (lake_category==-1) lake_category = 16  ! No separate lake category!

        elseif (trim(LU_Categories)=="USGS-RUC") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = 24
            if (water_category==-1) water_category = 16
            if (lake_category==-1) lake_category = 28
            ! also note, lakes_category = 28
            ! write(*,*) "  WARNING: not handling lake category (28)"

        elseif (trim(LU_Categories)=="MODI-RUC") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            if (lake_category==-1) lake_category = 21
            ! also note, lakes_category = 21
            ! write(*,*) "  WARNING: not handling lake category (21)"

        elseif (trim(LU_Categories)=="NLCD40") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15 ! and 22?
            ! if (water_category==-1) water_category = 17 ! and 21 'Open Water'
            write(*,*) "  WARNING: not handling all varients of categories (e.g. permanent_snow=15 is, but permanent_snow_ice=22 is not)"
        endif

    end subroutine set_default_LU_categories

    !> -------------------------------
    !! Add variables needed by all domains to the list of requested variables
    !!
    !! -------------------------------
    subroutine default_var_requests(options)
        type(options_t) :: options
        
        ! List the variables that are required to be allocated for any domain
        call options%alloc_vars(                                                    &
                     [kVARS%z,                      kVARS%z_interface,              &
                      kVARS%dz,                     kVARS%dz_interface,             &
                      kVARS%dzdx,                   kVARS%dzdy,                     &
                      kVARS%u,                      kVARS%v,                        &
                      kVARS%w,                      kVARS%w_real,                   &
                      kVARS%surface_pressure,       kVARS%roughness_z0,             &
                      kVARS%terrain,                kVARS%pressure,                 &
                      kVARS%temperature,            kVARS%pressure_interface,       &
                      kVARS%exner,                  kVARS%potential_temperature,    &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude,              &
                      kVARS%temperature_interface,  kVars%density])

        ! List the variables that are required for any restart
        call options%restart_vars(                                                  &
                     [kVARS%z,                                                      &
                      kVARS%terrain,                kVARS%potential_temperature,    &
                      kVARS%latitude,               kVARS%longitude,                &
                      kVARS%u,                      kVARS%v,                        &
                      kVARS%w,                      kVARS%w_real,                   &
                      kVARS%u_latitude,             kVARS%u_longitude,              &
                      kVARS%v_latitude,             kVARS%v_longitude               ])

    end subroutine default_var_requests

    !> -------------------------------
    !! Add list of new variables to a list of variables
    !!
    !! Adds one to the associated index of the list and returns an error
    !! Sets Error/=0 if any of the variables suggested are outside the bounds of the list
    !!
    !! -------------------------------
    subroutine add_to_varlist(varlist, varids, error)
        implicit none
        integer, intent(inout)  :: varlist(:)
        integer, intent(in)     :: varids(:)
        integer, intent(out), optional  :: error

        integer :: i, ierr

        ierr=0
        do i=1,size(varids)
            if (varids(i) <= size(varlist)) then
                varlist( varids(i) ) = varlist( varids(i) ) + 1
            else
                if (STD_OUT_PE) write(*,*) "  WARNING: trying to add var outside of permitted list:",varids(i), size(varlist)
                ierr=1
            endif
        enddo

        if (present(error)) error=ierr

    end subroutine add_to_varlist


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be allocated
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine alloc_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_to_allocate,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_to_allocate,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine alloc_vars

    !> -------------------------------
    !! Overwrites options%forcing struct to expect "forcing" data from the parent nest
    !!
    !! Mostly just changes forcing var names
    !!
    !! -------------------------------
    module subroutine setup_synthetic_forcing(this)
        implicit none
        class(options_t),  intent(inout):: this
        integer :: ierr, i

        this%forcing%qv_is_relative_humidity = .false.
        this%forcing%qv_is_spec_humidity = .false.  
        this%forcing%t_is_potential = .True.
        this%forcing%z_is_geopotential = .False.
        this%forcing%z_is_on_interface = .False.
        this%forcing%time_varying_z = .False.
        this%forcing%t_offset = 0.0
        this%forcing%compute_z = .False.

        ! Now set the forcing variable names -- these are the same as the domain variable names
        ! NOTE: temperature must be the first of the forcing variables read
        this%forcing%vars_to_read(:) = ""
        this%forcing%dim_list(:) = 0
        i = 1
        call set_nml_var(this%forcing%zvar, get_varname( kVARS%z ), 'zvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%uvar, get_varname( kVARS%u ), 'uvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%vvar, get_varname( kVARS%v ), 'vvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%wvar, get_varname( kVARS%w_real ), 'wvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%qvvar, get_varname( kVARS%water_vapor ), 'qvvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%tvar, get_varname( kVARS%potential_temperature ), 'tvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%pvar, get_varname( kVARS%pressure ), 'pvar', this%forcing, i, no_check=.True.)
        call set_nml_var(this%forcing%latvar, get_varname( kVARS%latitude ), 'latvar', this%forcing, no_check=.True.)
        call set_nml_var(this%forcing%lonvar, get_varname( kVARS%longitude ), 'lonvar', this%forcing, no_check=.True.)
        call set_nml_var(this%forcing%hgtvar, get_varname( kVARS%terrain ), 'hgtvar', this%forcing, no_check=.True.)

        if (0<this%vars_to_allocate( kVARS%cloud_water) ) call set_nml_var(this%forcing%qcvar, get_varname( kVARS%cloud_water ), 'qcvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%cloud_ice) ) call set_nml_var(this%forcing%qivar, get_varname( kVARS%cloud_ice ), 'qivar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%rain_in_air) ) call set_nml_var(this%forcing%qrvar, get_varname( kVARS%rain_in_air ), 'qrvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%graupel_in_air) ) call set_nml_var(this%forcing%qgvar, get_varname( kVARS%graupel_in_air ), 'qgvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%snow_in_air) ) call set_nml_var(this%forcing%qsvar, get_varname( kVARS%snow_in_air ), 'qsvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%cloud_number_concentration) ) call set_nml_var(this%forcing%qncvar, get_varname( kVARS%cloud_number_concentration ), 'qncvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice_number_concentration) ) call set_nml_var(this%forcing%qnivar, get_varname( kVARS%ice_number_concentration ), 'qnivar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%rain_number_concentration) ) call set_nml_var(this%forcing%qnrvar, get_varname( kVARS%rain_number_concentration ), 'qnrvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%graupel_number_concentration) ) call set_nml_var(this%forcing%qngvar, get_varname( kVARS%graupel_number_concentration ), 'qngvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%snow_number_concentration) ) call set_nml_var(this%forcing%qnsvar, get_varname( kVARS%snow_number_concentration ), 'qnsvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice2_mass) ) call set_nml_var(this%forcing%i2mvar, get_varname( kVARS%ice2_mass ), 'i2mvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice3_mass) ) call set_nml_var(this%forcing%i3mvar, get_varname( kVARS%ice3_mass ), 'i3mvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice2_number) ) call set_nml_var(this%forcing%i2nvar, get_varname( kVARS%ice2_number ), 'i2nvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice3_number) ) call set_nml_var(this%forcing%i3nvar, get_varname( kVARS%ice3_number ), 'i3nvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice1_a) ) call set_nml_var(this%forcing%i1avar, get_varname( kVARS%ice1_a ), 'i1avar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice2_a) ) call set_nml_var(this%forcing%i2avar, get_varname( kVARS%ice2_a ), 'i2avar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice3_a) ) call set_nml_var(this%forcing%i3avar, get_varname( kVARS%ice3_a ), 'i3avar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice1_c) ) call set_nml_var(this%forcing%i1cvar, get_varname( kVARS%ice1_c ), 'i1cvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice2_c) ) call set_nml_var(this%forcing%i2cvar, get_varname( kVARS%ice2_c ), 'i2cvar', this%forcing, i, no_check=.True.)
        if (0<this%vars_to_allocate( kVARS%ice3_c) ) call set_nml_var(this%forcing%i3cvar, get_varname( kVARS%ice3_c ), 'i3cvar', this%forcing, i, no_check=.True.)
        
    end subroutine setup_synthetic_forcing


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be output in a restart file
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine restart_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_for_restart,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_for_restart,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine


    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be advected
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine advect_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_to_advect,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_to_advect,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine

    !> -------------------------------
    !! Add a set of variable(s) to the internal list of variables to be exchanged-only
    !!
    !! Sets error /= 0 if an error occurs in add_to_varlist
    !!
    !! -------------------------------
    module subroutine exch_vars(this, input_vars, var_idx, error)
        class(options_t),  intent(inout):: this
        integer, optional, intent(in)  :: input_vars(:)
        integer, optional, intent(in)  :: var_idx
        integer, optional, intent(out) :: error

        integer :: ierr

        ierr=0
        if (present(var_idx)) then
            call add_to_varlist(this%vars_to_exch,[var_idx], ierr)
        endif

        if (present(input_vars)) then
            call add_to_varlist(this%vars_to_exch,input_vars, ierr)
        endif

        if (present(error)) error=ierr

    end subroutine


    !> -------------------------------
    !! Check the version number in the namelist file and compare to the current model version
    !!
    !! If the namelist version doesn't match, print the differences between that version and this
    !! and STOP execution
    !!
    !! -------------------------------
    subroutine version_check(options)
        type(general_options_type),intent(inout)  :: options


        if (options%version.ne.kVERSION_STRING) then
            if (STD_OUT_PE) write(*,*) "  Model version does not match namelist version"
            if (STD_OUT_PE) write(*,*) "    Model version: ",kVERSION_STRING
            if (STD_OUT_PE) write(*,*) "    Namelist version: ",trim(options%version)
            call print_model_diffs(options%version)
            stop
        endif

        if (STD_OUT_PE) write(*,*) "    Model version: ",trim(options%version)

    end subroutine version_check


end submodule
