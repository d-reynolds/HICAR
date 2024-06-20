submodule(options_interface) options_implementation

    use icar_constants
    use mod_wrf_constants,          only : piconst
    use io_routines,                only : io_newunit, check_variable_present, check_file_exists
    use time_io,                    only : find_timestep_in_file
    use time_delta_object,          only : time_delta_t
    use time_object,                only : Time_type
    use string,                     only : str
    use model_tracking,             only : print_model_diffs

    use convection,                 only : cu_var_request
    use land_surface,               only : lsm_var_request
    use surface_layer,              only : sfc_var_request
    use radiation,                  only : ra_var_request
    use microphysics,               only : mp_var_request
    use advection,                  only : adv_var_request
    use wind,                       only : wind_var_request
    use planetary_boundary_layer,   only : pbl_var_request

    use output_metadata,            only : get_varname
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
    module subroutine init(this, namelist_file, info_only, gen_nml, only_namelist_check)
        implicit none
        class(options_t),   intent(inout)  :: this
        character(len=*),   intent(in)     :: namelist_file
        logical,            intent(in)     :: info_only, gen_nml, only_namelist_check

!       reads a series of options from a namelist file and stores them in the
!       options data structure
        integer :: i

        if (STD_OUT_PE .and. .not.(info_only .or. gen_nml)) write(*,*) "Using options file = ", trim(namelist_file)
        if (gen_nml) call set_namelist(namelist_file)

        call general_namelist(      namelist_file,   this%general, info_only=info_only, gen_nml=gen_nml)
        call domain_namelist(       namelist_file,   this%domain, info_only=info_only, gen_nml=gen_nml)
        call forcing_namelist(      namelist_file,   this%forcing, info_only=info_only, gen_nml=gen_nml)

        call physics_namelist(      namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call time_parameters_namelist(         namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call lt_parameters_namelist(    namelist_file,    this, info_only=info_only, gen_nml=gen_nml)
        call mp_parameters_namelist(    namelist_file,    this, info_only=info_only, gen_nml=gen_nml)
        call adv_parameters_namelist(   namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call lsm_parameters_namelist(   namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call cu_parameters_namelist(    namelist_file,    this, info_only=info_only, gen_nml=gen_nml)
        call rad_parameters_namelist(   namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call pbl_parameters_namelist(   namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call sfc_parameters_namelist(   namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        
        call wind_namelist(         namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call output_namelist(       namelist_file,   this, info_only=info_only, gen_nml=gen_nml)
        call restart_namelist(      namelist_file,   this, info_only=info_only, gen_nml=gen_nml)

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

        call collect_var_requests(this)

        ! check for any inconsistencies in the options requested
        call options_check(this)

        ! If this run was just done to check the namelist options, stop now
        if (only_namelist_check) then
            if (STD_OUT_PE) write(*,*) 'Namelist options check complete.'
            stop
        endif
    end subroutine init

    subroutine collect_var_requests(options)
        type(options_t) :: options

        call collect_physics_requests(options)
        call default_var_requests(options)

    end subroutine collect_var_requests

    !> -------------------------------
    !! Call all physics driver var_request routines
    !!
    !! var_request routines allow physics modules to requested
    !! which variables they need to have allocated, advected, and written in restart files
    !!
    !! -------------------------------
    subroutine collect_physics_requests(options)
        type(options_t) :: options

        call ra_var_request(options)
        call lsm_var_request(options)
        call sfc_var_request(options)
        call pbl_var_request(options)
        call cu_var_request(options)
        call mp_var_request(options)
        call adv_var_request(options)
        call wind_var_request(options)

    end subroutine collect_physics_requests


    !> -------------------------------
    !! Checks options in the options data structure for consistency
    !!
    !! Stops or prints a large warning depending on warning level requested and error found
    !!
    !! -------------------------------
    subroutine options_check(options)
        ! Minimal error checking on option settings
        implicit none
        type(options_t), intent(inout)::options
        integer :: i

        !clean output var list
        do i=1, size(options%output%vars_for_output)
            if ((options%output%vars_for_output(i)+options%vars_for_restart(i) > 0) .and. (options%vars_to_allocate(i) <= 0)) then
                !if (STD_OUT_PE) write(*,*) 'variable ',trim(get_varname(options%vars_to_allocate(i))),' requested for output/restart, but was not allocated by one of the modules'
                if (options%output%vars_for_output(i) > 0) options%output%vars_for_output(i) = 0
                if (options%vars_for_restart(i) > 0) options%vars_for_restart(i) = 0
            endif
        enddo

        ! if using a real LSM, feedback will probably keep hot-air from getting even hotter, so not likely a problem
        if ((options%physics%landsurface>0).and.(options%physics%boundarylayer==0)) then
            if (STD_OUT_PE) write(*,*) ""
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "WARNING, Using surfaces fluxes (lsm>0) without a PBL scheme may overheat the surface and CRASH the model."
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
        endif

        ! prior to v 0.9.3 this was assumed, so throw a warning now just in case.
        if ((options%forcing%z_is_geopotential .eqv. .False.).and. &
            (options%forcing%zvar=="PH")) then
            if (STD_OUT_PE) write(*,*) ""
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "WARNING z variable is not assumed to be geopotential height when it is 'PH'."
            if (STD_OUT_PE) write(*,*) "WARNING If z is geopotential, set z_is_geopotential=True in the namelist."
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
        endif
        if ((options%forcing%z_is_geopotential .eqv. .True.).and. &
            (options%forcing%z_is_on_interface .eqv. .False.)) then
            if (STD_OUT_PE) write(*,*) ""
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "WARNING geopotential height is no longer assumed to be on model interface levels."
            if (STD_OUT_PE) write(*,*) "WARNING To interpolate geopotential, set z_is_on_interface=True in the namelist. "
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
        endif
        
        !! MJ added
        if ((options%physics%radiation_downScaling==1).and.(options%physics%radiation==0)) then
            if (STD_OUT_PE) write(*,*) ""
            if (STD_OUT_PE) write(*,*) "STOP STOP STOP"
            if (STD_OUT_PE) write(*,*) "STOP, Running radiation_downScaling=1 cannot not be used with rad=0"
            if (STD_OUT_PE) write(*,*) "STOP STOP STOP"
            stop
        endif
        
        if (options%time%RK3) then
            if (max(options%adv%h_order,options%adv%v_order)==5 .and. options%time%cfl_reduction_factor > 1.4) then
                if (STD_OUT_PE) write(*,*) "CFL reduction factor should be less than 1.4 when horder or vorder = 5, limiting to 1.4"
                options%time%cfl_reduction_factor = min(1.4,options%time%cfl_reduction_factor)
            elseif (max(options%adv%h_order,options%adv%v_order)==3 .and. options%time%cfl_reduction_factor > 1.6) then
                if (STD_OUT_PE) write(*,*) "CFL reduction factor should be less than 1.6 when horder or vorder = 3, limiting to 1.6"
                options%time%cfl_reduction_factor = min(1.6,options%time%cfl_reduction_factor)
            endif
        else
            if (options%time%cfl_reduction_factor > 1.0) then   
                if (STD_OUT_PE) write(*,*) "CFL reduction factor should be less than 1.0 when RK3=.False., limiting to 1.0"
                options%time%cfl_reduction_factor = min(1.0,options%time%cfl_reduction_factor)
            endif
        endif
        
        if (options%wind%alpha_const > 0) then
            if (options%wind%alpha_const > 1.0) then
                if (STD_OUT_PE) write(*,*) "Alpha currently limited to values between 0.01 and 1.0, setting to 1.0"
                options%wind%alpha_const = 1.0
            else if (options%wind%alpha_const < 0.01) then
                if (STD_OUT_PE) write(*,*) "Alpha currently limited to values between 0.01 and 1.0, setting to 0.01"
                options%wind%alpha_const = 0.01
            endif
        endif
        
        ! should warn user if lsm is run without radiation
        if ((options%physics%landsurface>kLSM_BASIC .or. options%physics%snowmodel>0).and.(options%physics%radiation==0)) then
            if (STD_OUT_PE) write(*,*) ""
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
            if (STD_OUT_PE) write(*,*) "WARNING, Using land surface model without radiation input does not make sense."
            if (STD_OUT_PE) write(*,*) "WARNING WARNING WARNING"
        endif

        if (options%physics%landsurface>0) then
            options%sfc%isfflx = 1
            options%sfc%scm_force_flux = 1
        endif
        
        !Allow for microphysics precipitation partitioning with NoahMP if using a snow model
        if (options%physics%snowmodel>0) then
            options%lsm%nmp_opt_snf = 4
        endif

        ! check if the last entry in dz_levels is zero, which would indicate that nz is larger than the number
        ! of entries in dz_levels, or that the user passed bad data
        if (options%domain%nz > 1) then
            if (options%domain%dz_levels(options%domain%nz) == 0) then
                if (STD_OUT_PE) write(*,*) "nz is larger than the number of entries in dz_levels, or the last entry in dz_levels is zero."
                stop
            endif
        endif

        ! check if start time is before end time
        if (options%general%start_time >= options%general%end_time) then
            if (STD_OUT_PE) write(*,*) "Start time must be before end time"
            stop
        endif

        ! check if restart_time is between start and end time
        if (options%restart%restart) then
            if (options%restart%restart_time < options%general%start_time .or. &
                options%restart%restart_time > options%general%end_time) then
                if (STD_OUT_PE) write(*,*) "Restart time must be between start and end time"
                stop
            endif
        endif



        ! Check if supporting files exist, if they are needed by physics modules
        if (options%physics%landsurface==kLSM_NOAHMP) then
            if (STD_OUT_PE) write(*,*) '  NoahMP LSM turned on, checking for supporting files...'
            call check_file_exists('GENPARM.TBL', message='GENPARM.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('MPTABLE.TBL', message='MPTABLE.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('SOILPARM.TBL', message='SOILPARM.TBL file does not exist. This should be in the same directory as the namelist.')
            call check_file_exists('VEGPARM.TBL', message='VEGPARM.TBL file does not exist. This should be in the same directory as the namelist.')
        endif
        if (options%physics%radiation==kRA_RRTMG) then
            if (STD_OUT_PE) write(*,*) '  RRTMG radiation turned on, checking for supporting files...'
            call check_file_exists('rrtmg_support/forrefo_1.nc', message='At least one of the RRTMG supporting files does not exist. These files should be in a folder "rrtmg_support" placed in the same directory as the namelist.')
        endif
        if (options%physics%microphysics==kMP_ISHMAEL) then
            if (STD_OUT_PE) write(*,*) '  ISHMAEL microphysics turned on, checking for supporting files...'
            call check_file_exists('mp_support/ishmael_gamma_tab.nc', message='At least one of the ISHMAEL supporting files does not exist. These files should be in a folder "mp_support" placed in the same directory as the namelist.')
        endif

    end subroutine options_check



    !> -------------------------------
    !! Read physics options to use from a namelist file
    !!
    !! -------------------------------
    subroutine physics_namelist(filename,options, info_only, gen_nml)
        implicit none
        character(len=*),intent(in)     :: filename
        type(options_t), intent(inout)  :: options
        logical, intent(in), optional  :: info_only, gen_nml
        integer :: name_unit, rc
        !variables to be used in the namelist
        integer :: pbl, lsm, sfc, sm, water, mp, rad, conv, adv, wind, radiation_downscaling
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
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'physics' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

        !store options
        call set_nml_var(options%physics%boundarylayer, pbl, 'pbl')
        call set_nml_var(options%physics%landsurface, lsm, 'lsm')
        call set_nml_var(options%physics%surfacelayer, sfc, 'sfc')
        call set_nml_var(options%physics%snowmodel, sm, 'sm')
        call set_nml_var(options%physics%watersurface, water, 'water')
        call set_nml_var(options%physics%microphysics, mp, 'mp')
        call set_nml_var(options%physics%radiation, rad, 'rad')
        call set_nml_var(options%physics%convection, conv, 'conv')
        call set_nml_var(options%physics%advection, adv, 'adv')
        call set_nml_var(options%physics%windtype, wind, 'wind')
        call set_nml_var(options%physics%radiation_downScaling, radiation_downscaling, 'radiation_downscaling')


    end subroutine physics_namelist


    !> -------------------------------
    !! Check that a required input variable is present
    !!
    !! If not present, halt the program
    !!
    !! -------------------------------
    subroutine require_var(inputvar, var_name)
        implicit none
        character(len=*), intent(in) :: inputvar
        character(len=*), intent(in) :: var_name

        if (trim(inputvar)=="") then
            if (STD_OUT_PE) write(*,*) "Variable: ",trim(var_name), " is required."
            stop
        endif

    end subroutine require_var

    !> -------------------------------
    !! Initialize the variable names to be written to standard output
    !!
    !! Reads the output_list namelist
    !!
    !! -------------------------------
    subroutine output_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(options_t),              intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, i, j, status, frames_per_outfile
        real    :: outputinterval
        logical :: file_exists, print_info, gennml

        ! Local variables
        character(len=kMAX_FILE_LENGTH) :: output_file
        character(len=kMAX_NAME_LENGTH), allocatable :: output_vars(:)
        
        namelist /output/ output_vars, outputinterval, frames_per_outfile, output_file

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        allocate(output_vars(kMAX_STORAGE_VARS))

        call set_nml_var_default(outputinterval, 'outputinterval', print_info, gennml)
        call set_nml_var_default(frames_per_outfile, 'frames_per_outfile', print_info, gennml)
        call set_nml_var_default(output_file, 'output_file', print_info, gennml)
        call set_nml_var_default(output_vars, 'output_vars', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc, nml=output)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'output' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

        do j=1, kMAX_STORAGE_VARS
            if (trim(output_vars(j)) /= "") then
                do i=1, kMAX_STORAGE_VARS
                    if (trim(get_varname(i)) == trim(output_vars(j))) then
                        call add_to_varlist(options%output%vars_for_output, [i])
                    endif
                enddo
            endif
        enddo        

        call set_nml_var(options%output%output_file, output_file, 'output_file')
        call set_nml_var(options%output%outputinterval, outputinterval, 'outputinterval')
        call options%output%output_dt%set(seconds=outputinterval)
        call set_nml_var(options%output%frames_per_outfile, frames_per_outfile, 'frames_per_outfile')

    end subroutine output_namelist

    !> -------------------------------
    !! Initialize the variable names to be written to standard output
    !!
    !! Reads the restart_list namelist
    !!
    !! -------------------------------
    subroutine restart_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(options_t),              intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        integer    :: restartinterval
        logical :: file_exists, restart_run, print_info, gennml

        ! Local variables
        character(len=kMAX_FILE_LENGTH) :: restart_out_file, restart_in_file
        
        integer :: restart_step                         ! time step relative to the start of the restart file
        character(len=MAXFILELENGTH) :: restart_date    ! date to restart
        type(Time_type) :: restart_time, time_at_step   ! restart date as a modified julian day        

        namelist /restart/  restartinterval, restart_out_file, restart_in_file, restart_date, restart_run

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(restartinterval, 'restartinterval', print_info, gennml)
        call set_nml_var_default(restart_out_file, 'restart_out_file', print_info, gennml)
        call set_nml_var_default(restart_in_file, 'restart_in_file', print_info, gennml)
        call set_nml_var_default(restart_date, 'restart_date', print_info, gennml)
        call set_nml_var_default(restart_run, 'restart_run', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc, nml=restart)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'restart' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

        call set_nml_var(options%restart%restart_count, restartinterval, 'restartinterval')
        call set_nml_var(options%restart%restart_out_file, restart_out_file, 'restart_out_file')
        call set_nml_var(options%restart%restart_in_file, restart_in_file, 'restart_in_file')
        call set_nml_var(options%restart%restart, restart_run, 'restart_run')

        !If the user did not ask for a restart run, leave the function now
        if (.not.(options%restart%restart)) return
        
        ! calculate the modified julian day for th restart date given
        call restart_time%init(options%general%calendar)
        if (restart_date=="") then
            if (STD_OUT_PE) write(*,*) "ERROR: restart_date must be specified in the namelist"
            stop
        else
            call restart_time%set(restart_date)
        endif
                              
        write(options%restart%restart_in_file, '(A,A,".nc")')    &
                trim(restart_in_file),   &
                trim(restart_time%as_string('(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'))

        call check_file_exists(options%restart%restart_in_file, message='Restart file does not exist.')

        ! find the time step that most closely matches the requested restart time (<=)
        restart_step = find_timestep_in_file(options%restart%restart_in_file, 'time', restart_time, time_at_step)

        ! check to see if we actually udpated the restart date and print if in a more verbose mode
        if (options%general%debug) then
            if (restart_time /= time_at_step) then
                if (STD_OUT_PE) write(*,*) " updated restart date: ", trim(time_at_step%as_string())
            endif
        endif

        restart_time = time_at_step

        ! save the parameters in the master options structure
        options%restart%restart_step_in_file = restart_step
        options%restart%restart_time         = restart_time


        if (options%general%debug) then
            if (STD_OUT_PE) write(*,*) " ------------------ "
            if (STD_OUT_PE) write(*,*) "RESTART INFORMATION"
            if (STD_OUT_PE) write(*,*) "mjd",         options%restart%restart_time%mjd()
            if (STD_OUT_PE) write(*,*) "date:",       trim(options%restart%restart_time%as_string())
            if (STD_OUT_PE) write(*,*) "file:",   trim(options%restart%restart_in_file)
            if (STD_OUT_PE) write(*,*) "forcing step",options%restart%restart_step_in_file
            if (STD_OUT_PE) write(*,*) " ------------------ "
        endif

    end subroutine restart_namelist

    !> -------------------------------
    !! Initialize the variable names to be read
    !!
    !! Reads the var_list namelist
    !!
    !! -------------------------------
    subroutine forcing_namelist(filename,options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(forcing_options_type), intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, i, j, nfiles
        logical :: compute_p, print_info, gennml
        logical :: limit_rh, z_is_geopotential, z_is_on_interface,&
                   time_varying_z, t_is_potential, qv_is_spec_humidity, &
                   qv_is_relative_humidity
        real    :: t_offset, inputinterval

        character(len=MAXFILELENGTH) :: forcing_file_list
        character(len=MAXFILELENGTH), allocatable :: boundary_files(:)

        character(len=MAXVARLENGTH) :: latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon,wvar,zvar,  &
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
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'forcing' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

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

        nfiles = read_forcing_file_names(forcing_file_list, boundary_files)

        if (nfiles==0) then
            stop "No boundary conditions files specified."
        endif

        allocate(options%boundary_files(nfiles))
        options%boundary_files(1:nfiles) = boundary_files(1:nfiles)
        deallocate(boundary_files)

        compute_p = .False.
        if ((pvar=="") .and. ((pslvar/="") .or. (psvar/=""))) compute_p = .True.
        if (compute_p) then
            if ((pslvar == "").and.(hgtvar == "")) then
                write(*,*) "ERROR: if surface pressure is used to compute air pressure, then surface height must be specified"
                error stop
            endif
        endif

        options%compute_z = .False.
        if ((zvar=="") .and. ((pslvar/="") .or. (psvar/=""))) options%compute_z = .True.
        if (options%compute_z) then
            if (pvar=="") then
                if (STD_OUT_PE) write(*,*) "ERROR: either pressure (pvar) or atmospheric level height (zvar) must be specified"
                error stop
            endif
        endif

        call set_nml_var(options%t_offset, t_offset, 't_offset')
        call set_nml_var(options%limit_rh, limit_rh, 'limit_rh')
        call set_nml_var(options%z_is_geopotential, z_is_geopotential, 'z_is_geopotential')
        call set_nml_var(options%z_is_on_interface, z_is_on_interface, 'z_is_on_interface')
        call set_nml_var(options%time_varying_z, time_varying_z, 'time_varying_z')
        call set_nml_var(options%t_is_potential, t_is_potential, 't_is_potential')
        call set_nml_var(options%qv_is_relative_humidity, qv_is_relative_humidity, 'qv_is_relative_humidity')
        call set_nml_var(options%qv_is_spec_humidity, qv_is_spec_humidity, 'qv_is_spec_humidity')
        call set_nml_var(options%inputinterval, inputinterval, 'inputinterval')

        call options%input_dt%set(seconds=inputinterval)

        ! NOTE: temperature must be the first of the forcing variables read
        options%vars_to_read(:) = ""
        i = 1
        call set_nml_var(options%tvar, tvar, 'tvar', options, i)
        call set_nml_var(options%latvar, latvar, 'latvar', options)
        call set_nml_var(options%lonvar, lonvar, 'lonvar', options)
        call set_nml_var(options%hgtvar, hgtvar, 'hgtvar', options, i)
        call set_nml_var(options%uvar, uvar, 'uvar', options, i)
        call set_nml_var(options%ulat, ulat, 'ulat', options)
        call set_nml_var(options%ulon, ulon, 'ulon', options)
        call set_nml_var(options%vvar, vvar, 'vvar', options, i)
        call set_nml_var(options%vlat, vlat, 'vlat', options)
        call set_nml_var(options%vlon, vlon, 'vlon', options)
        call set_nml_var(options%wvar, wvar, 'wvar', options, i)
        call set_nml_var(options%pslvar, pslvar, 'pslvar', options, i)
        call set_nml_var(options%psvar, psvar, 'psvar', options, i)
        call set_nml_var(options%qvvar, qvvar, 'qvvar', options, i)
        call set_nml_var(options%qcvar, qcvar, 'qcvar', options, i)
        call set_nml_var(options%qivar, qivar, 'qivar', options, i)
        call set_nml_var(options%qrvar, qrvar, 'qrvar', options, i)
        call set_nml_var(options%qgvar, qgvar, 'qgvar', options, i)
        call set_nml_var(options%qsvar, qsvar, 'qsvar', options, i)
        call set_nml_var(options%qncvar, qncvar, 'qncvar', options, i)
        call set_nml_var(options%qnivar, qnivar, 'qnivar', options, i)
        call set_nml_var(options%qnrvar, qnrvar, 'qnrvar', options, i)
        call set_nml_var(options%qngvar, qngvar, 'qngvar', options, i)
        call set_nml_var(options%qnsvar, qnsvar, 'qnsvar', options, i)
        call set_nml_var(options%i2mvar, i2mvar, 'i2mvar', options, i)
        call set_nml_var(options%i3mvar, i3mvar, 'i3mvar', options, i)
        call set_nml_var(options%i2nvar, i2nvar, 'i2nvar', options, i)
        call set_nml_var(options%i3nvar, i3nvar, 'i3nvar', options, i)
        call set_nml_var(options%i1avar, i1avar, 'i1avar', options, i)
        call set_nml_var(options%i2avar, i2avar, 'i2avar', options, i)
        call set_nml_var(options%i3avar, i3avar, 'i3avar', options, i)
        call set_nml_var(options%i1cvar, i1cvar, 'i1cvar', options, i)
        call set_nml_var(options%i2cvar, i2cvar, 'i2cvar', options, i)
        call set_nml_var(options%i3cvar, i3cvar, 'i3cvar', options, i)
        call set_nml_var(options%shvar, shvar, 'shvar', options, i)
        call set_nml_var(options%lhvar, lhvar, 'lhvar', options, i)
        call set_nml_var(options%swdown_var, swdown_var, 'swdown_var', options, i)
        call set_nml_var(options%lwdown_var, lwdown_var, 'lwdown_var', options, i)
        call set_nml_var(options%sst_var, sst_var, 'sst_var', options, i)
        call set_nml_var(options%pblhvar, pblhvar, 'pblhvar', options, i)
        call set_nml_var(options%time_var, time_var, 'time_var', options)

        ! vertical coordinate
        ! if (options%time_varying_z) then
        if (options%compute_z) then
            zvar = "height_computed"
            options%zvar        = zvar; options%vars_to_read(i) = zvar;      options%dim_list(i) = -3;    i = i + 1
        else
            call set_nml_var(options%zvar, zvar, 'zvar', options, i)
        endif


        if (compute_p) then
            pvar = "air_pressure_computed"
            options%pvar        = pvar  ; options%vars_to_read(i) = pvar;       options%dim_list(i) = -3;   i = i + 1
        else
            call set_nml_var(options%pvar, pvar, 'pvar', options, i)
        endif


    end subroutine forcing_namelist


    !> -------------------------------
    !! Initialize the main parameter options
    !!
    !! Reads parameters for the ICAR simulation
    !! These include setting flags that request other namelists be read
    !!
    !! -------------------------------
    subroutine general_namelist(filename,options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(general_options_type), intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        ! parameters to read
        logical :: interactive, debug, &
                   use_mp_options, use_lt_options, use_adv_options, use_lsm_options, &
                   use_cu_options, use_rad_options, use_pbl_options, use_sfc_options, &
                   use_wind_options

        character(len=MAXFILELENGTH) :: calendar, start_date, end_date
        character(len=MAXVARLENGTH)  :: version, comment, phys_suite

        namelist /general/    debug, interactive, calendar,          &
                              version, comment, phys_suite,                 &
                              start_date, end_date,         &
                              use_mp_options,     &
                              use_lt_options,     &
                              use_lsm_options,    &
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
        call set_nml_var_default(use_mp_options, 'use_mp_options', print_info, gennml)
        call set_nml_var_default(use_lt_options, 'use_lt_options', print_info, gennml)
        call set_nml_var_default(use_adv_options, 'use_adv_options', print_info, gennml)
        call set_nml_var_default(use_cu_options, 'use_cu_options', print_info, gennml)
        call set_nml_var_default(use_lsm_options, 'use_lsm_options', print_info, gennml)
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
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'general' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

        call set_nml_var(options%calendar, calendar, 'calendar')
        call set_nml_var(options%version, version, 'version')
        call set_nml_var(options%comment, comment, 'comment')
        call set_nml_var(options%phys_suite, phys_suite, 'phys_suite')
        call set_nml_var(options%debug, debug, 'debug')
        call set_nml_var(options%interactive, interactive, 'interactive')
        call set_nml_var(options%use_mp_options, use_mp_options, 'use_mp_options')
        call set_nml_var(options%use_lt_options, use_lt_options, 'use_lt_options')
        call set_nml_var(options%use_adv_options, use_adv_options, 'use_adv_options')
        call set_nml_var(options%use_lsm_options, use_lsm_options, 'use_lsm_options')
        call set_nml_var(options%use_cu_options, use_cu_options, 'use_cu_options')
        call set_nml_var(options%use_rad_options, use_rad_options, 'use_rad_options')
        call set_nml_var(options%use_pbl_options, use_pbl_options, 'use_pbl_options')
        call set_nml_var(options%use_sfc_options, use_sfc_options, 'use_sfc_options')
        call set_nml_var(options%use_wind_options, use_wind_options, 'use_wind_options')

        call version_check(options)

        if (trim(start_date)/="") then
            call options%start_time%init(calendar)
            call options%start_time%set(start_date)
        else
            stop 'start date must be supplied in namelist'
        endif
        if (trim(end_date)/="") then
            call options%end_time%init(calendar)
            call options%end_time%set(end_date)
        else
            stop 'end date must be supplied in namelist'
        endif


        ! options are updated when complete
    end subroutine general_namelist


    !> -------------------------------
    !! Set up model levels, either read from a namelist, or from a default set of values
    !!
    !! Reads the z_info namelist or sets default values
    !!
    !! -------------------------------
    subroutine domain_namelist(filename,options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(domain_options_type), intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc, this_level
        logical :: print_info, gennml

        real, allocatable, dimension(:) :: dz_levels
        logical :: sleve, use_agl_height

        real :: dx, flat_z_height, decay_rate_L_topo, decay_rate_S_topo, sleve_n, agl_cap, max_agl_height
        integer :: nz, longitude_system, terrain_smooth_windowsize, terrain_smooth_cycles

        character(len=MAXFILELENGTH) :: init_conditions_file

        character(len=MAXVARLENGTH) :: landvar,lakedepthvar,hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi,           &
                                        snowh_var, soiltype_var, soil_t_var,soil_vwc_var,swe_var, soil_deept_var,           &
                                        vegtype_var,vegfrac_var, vegfracmax_var, albedo_var, lai_var, canwat_var, linear_mask_var, nsq_calibration_var,  &
                                        sinalpha_var, cosalpha_var


        character(len=MAXVARLENGTH) :: svf_var, hlm_var, slope_var, slope_angle_var, aspect_angle_var, shd_var  !!MJ added

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

        allocate(dz_levels(MAXLEVELS))
        call set_nml_var_default(dz_levels, 'dz_levels', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        open(io_newunit(name_unit), file=filename)
        read(name_unit,iostat=rc,nml=domain)
        close(name_unit)
        if (rc /= 0) then
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'domain' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif

        call require_var(lat_hi, "High-res Lat")
        call require_var(lon_hi, "High-res Lon")
        call require_var(hgt_hi, "High-res HGT")

        call set_nml_var(options%dx, dx, 'dx')
        call set_nml_var(options%nz, nz, 'nz')
        call set_nml_var(options%longitude_system, longitude_system, 'longitude_system')
        call set_nml_var(options%init_conditions_file, init_conditions_file, 'init_conditions_file')
        call set_nml_var(options%flat_z_height, flat_z_height, 'flat_z_height')
        call set_nml_var(options%sleve, sleve, 'sleve')
        call set_nml_var(options%terrain_smooth_windowsize, terrain_smooth_windowsize, 'terrain_smooth_windowsize')
        call set_nml_var(options%terrain_smooth_cycles, terrain_smooth_cycles, 'terrain_smooth_cycles')
        call set_nml_var(options%decay_rate_L_topo, decay_rate_L_topo, 'decay_rate_L_topo')
        call set_nml_var(options%decay_rate_S_topo, decay_rate_S_topo, 'decay_rate_S_topo')
        call set_nml_var(options%sleve_n, sleve_n, 'sleve_n')
        call set_nml_var(options%use_agl_height, use_agl_height, 'use_agl_height')
        call set_nml_var(options%agl_cap, agl_cap, 'agl_cap')

        ! NOTE: hgt_hi has to be the first of the variables read
        call set_nml_var(options%hgt_hi, hgt_hi, 'hgt_hi',options)
        call set_nml_var(options%landvar, landvar, 'landvar',options)
        call set_nml_var(options%lakedepthvar, lakedepthvar, 'lakedepthvar',options)
        call set_nml_var(options%lat_hi, lat_hi, 'lat_hi',options)
        call set_nml_var(options%lon_hi, lon_hi, 'lon_hi',options)
        call set_nml_var(options%ulat_hi, ulat_hi, 'ulat_hi',options)
        call set_nml_var(options%ulon_hi, ulon_hi, 'ulon_hi',options)
        call set_nml_var(options%vlat_hi, vlat_hi, 'vlat_hi',options)
        call set_nml_var(options%vlon_hi, vlon_hi, 'vlon_hi',options)
        call set_nml_var(options%soiltype_var, soiltype_var, 'soiltype_var',options)
        call set_nml_var(options%soil_t_var, soil_t_var, 'soil_t_var',options)
        call set_nml_var(options%soil_vwc_var, soil_vwc_var, 'soil_vwc_var',options)
        call set_nml_var(options%swe_var, swe_var, 'swe_var',options)
        call set_nml_var(options%snowh_var, snowh_var, 'snowh_var',options)
        call set_nml_var(options%soil_deept_var, soil_deept_var, 'soil_deept_var',options)
        call set_nml_var(options%vegtype_var, vegtype_var, 'vegtype_var',options)
        call set_nml_var(options%vegfrac_var, vegfrac_var, 'vegfrac_var',options)
        call set_nml_var(options%vegfracmax_var, vegfracmax_var, 'vegfracmax_var',options)
        call set_nml_var(options%albedo_var, albedo_var, 'albedo_var',options)
        call set_nml_var(options%lai_var, lai_var, 'lai_var',options)
        call set_nml_var(options%canwat_var, canwat_var, 'canwat_var',options)
        call set_nml_var(options%linear_mask_var, linear_mask_var, 'linear_mask_var',options)
        call set_nml_var(options%nsq_calibration_var, nsq_calibration_var, 'nsq_calibration_var',options)
        call set_nml_var(options%sinalpha_var, sinalpha_var, 'sinalpha_var',options)
        call set_nml_var(options%cosalpha_var, cosalpha_var, 'cosalpha_var',options)

        call set_nml_var(options%svf_var, svf_var, 'svf_var',options)
        call set_nml_var(options%hlm_var, hlm_var, 'hlm_var',options)
        call set_nml_var(options%slope_var, slope_var, 'slope_var',options)
        call set_nml_var(options%slope_angle_var, slope_angle_var, 'slope_angle_var',options)
        call set_nml_var(options%aspect_angle_var, aspect_angle_var, 'aspect_angle_var',options)
        call set_nml_var(options%shd_var, shd_var, 'shd_var',options)


        options%nz = nz
        call set_nml_var(options%nz, nz, 'nz')
        ! if nz wasn't specified in the namelist, we assume a HUGE number of levels
        ! so now we have to figure out what the actual number of levels read was
        if (options%nz==MAXLEVELS) then
            do this_level=1,MAXLEVELS-1
                if (dz_levels(this_level+1)<=0) then
                    options%nz=this_level
                    exit
                endif
            end do
            options%nz=this_level
        endif

        call set_nml_var(options%dz_levels(1:options%nz), dz_levels(1:options%nz), 'dz_levels')
        deallocate(dz_levels)

    end subroutine domain_namelist


    !> -------------------------------
    !! Initialize the microphysics options
    !!
    !! Reads the mp_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine mp_parameters_namelist(mp_filename,options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)    :: mp_filename
        type(options_t),    intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml

        logical :: print_info, gennml
        integer :: name_unit, rc

        real    :: Nt_c, TNO, am_s, rho_g, av_s, bv_s, fv_s, av_g, bv_g, av_i
        real    :: Ef_si, Ef_rs, Ef_rg, Ef_ri
        real    :: C_cubes, C_sqrd, mu_r, t_adjust
        logical :: Ef_rw_l, EF_sw_l
        integer :: top_mp_level
        real    :: update_interval_mp

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
        if (options%general%use_mp_options) then
            open(io_newunit(name_unit), file=mp_filename)
            read(name_unit,iostat=rc, nml=mp_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'mp_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        if (top_mp_level < 0) top_mp_level = options%domain%nz + top_mp_level

        call set_nml_var(options%mp%Nt_c, Nt_c, 'Nt_c')
        call set_nml_var(options%mp%TNO, TNO, 'TNO')
        call set_nml_var(options%mp%am_s, am_s, 'am_s')
        call set_nml_var(options%mp%rho_g, rho_g, 'rho_g')
        call set_nml_var(options%mp%av_s, av_s, 'av_s')
        call set_nml_var(options%mp%bv_s, bv_s, 'bv_s')
        call set_nml_var(options%mp%fv_s, fv_s, 'fv_s')

        call set_nml_var(options%mp%av_g, av_g, 'av_g')
        call set_nml_var(options%mp%bv_g, bv_g, 'bv_g')
        call set_nml_var(options%mp%av_i, av_i, 'av_i')
        call set_nml_var(options%mp%Ef_si, Ef_si, 'Ef_si')
        call set_nml_var(options%mp%Ef_rs, Ef_rs, 'Ef_rs')
        call set_nml_var(options%mp%Ef_rg, Ef_rg, 'Ef_rg')
        call set_nml_var(options%mp%Ef_ri, Ef_ri, 'Ef_ri')
        call set_nml_var(options%mp%mu_r, mu_r, 'mu_r')
        call set_nml_var(options%mp%t_adjust, t_adjust, 't_adjust')
        call set_nml_var(options%mp%C_cubes, C_cubes, 'C_cubes')
        call set_nml_var(options%mp%C_sqrd, C_sqrd, 'C_sqrd')
        call set_nml_var(options%mp%Ef_rw_l, Ef_rw_l, 'Ef_rw_l')
        call set_nml_var(options%mp%Ef_sw_l, Ef_sw_l, 'Ef_sw_l')

        call set_nml_var(options%mp%update_interval, update_interval_mp, 'update_interval_mp')
        call set_nml_var(options%mp%top_mp_level, top_mp_level, 'top_mp_level')

    end subroutine mp_parameters_namelist


    !> -------------------------------
    !! Initialize the Linear Theory options
    !!
    !! Reads the lt_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lt_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options
        logical, intent(in), optional  :: info_only, gen_nml

        integer :: name_unit, rc
        logical :: print_info, gennml

        integer :: vert_smooth
        logical :: variable_N           ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: smooth_nsq               ! Smooth the Calculated N^2 over vert_smooth vertical levels
        integer :: buffer                   ! number of grid cells to buffer around the domain MUST be >=1
        integer :: stability_window_size    ! window to average nsq over
        real    :: max_stability            ! limits on the calculated Brunt Vaisala Frequency
        real    :: min_stability            ! these may need to be a little narrower.
        real    :: linear_contribution      ! multiplier on uhat,vhat before adding to u,v
        real    :: linear_update_fraction   ! controls the rate at which the linearfield updates (should be calculated as f(in_dt))

        real    :: N_squared                ! static Brunt Vaisala Frequency (N^2) to use
        logical :: remove_lowres_linear     ! attempt to remove the linear mountain wave from the forcing low res model
        real    :: rm_N_squared             ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        real    :: rm_linear_contribution   ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        logical :: spatial_linear_fields    ! use a spatially varying linear wind perturbation
        logical :: linear_mask              ! use a spatial mask for the linear wind field
        logical :: nsq_calibration          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! Look up table generation parameters
        real    :: dirmax, dirmin
        real    :: spdmax, spdmin
        real    :: nsqmax, nsqmin
        integer :: n_dir_values, n_nsq_values, n_spd_values
        real    :: minimum_layer_size       ! Minimum vertical step to permit when computing LUT.
                                            ! If model layers are thicker, substepping will be used.

        ! parameters to control reading from or writing an LUT file
        logical :: read_LUT, write_LUT
        character(len=MAXFILELENGTH) :: u_LUT_Filename, v_LUT_Filename, LUT_Filename
        logical :: overwrite_lt_lut

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
        if (options%general%use_lt_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lt_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'lt_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        call set_nml_var(options%lt%variable_N, variable_N, 'variable_N')
        call set_nml_var(options%lt%smooth_nsq, smooth_nsq, 'smooth_nsq')
        call set_nml_var(options%lt%buffer, buffer, 'buffer')
        call set_nml_var(options%lt%stability_window_size, stability_window_size, 'stability_window_size')
        call set_nml_var(options%lt%max_stability, max_stability, 'max_stability')
        call set_nml_var(options%lt%min_stability, min_stability, 'min_stability')
        call set_nml_var(options%lt%vert_smooth, vert_smooth, 'vert_smooth')
        call set_nml_var(options%lt%N_squared, N_squared, 'N_squared')
        call set_nml_var(options%lt%linear_contribution, linear_contribution, 'linear_contribution')
        call set_nml_var(options%lt%remove_lowres_linear, remove_lowres_linear, 'remove_lowres_linear')
        call set_nml_var(options%lt%rm_N_squared, rm_N_squared, 'rm_N_squared')
        call set_nml_var(options%lt%rm_linear_contribution, rm_linear_contribution, 'rm_linear_contribution')
        call set_nml_var(options%lt%linear_update_fraction, linear_update_fraction, 'linear_update_fraction')
        call set_nml_var(options%lt%spatial_linear_fields, spatial_linear_fields, 'spatial_linear_fields')
        call set_nml_var(options%lt%linear_mask, linear_mask, 'linear_mask')
        call set_nml_var(options%lt%nsq_calibration, nsq_calibration, 'nsq_calibration')
        call set_nml_var(options%lt%dirmax, dirmax, 'dirmax')
        call set_nml_var(options%lt%dirmin, dirmin, 'dirmin')
        call set_nml_var(options%lt%spdmax, spdmax, 'spdmax')
        call set_nml_var(options%lt%spdmin, spdmin, 'spdmin')
        call set_nml_var(options%lt%nsqmax, nsqmax, 'nsqmax')
        call set_nml_var(options%lt%nsqmin, nsqmin, 'nsqmin')
        call set_nml_var(options%lt%n_dir_values, n_dir_values, 'n_dir_values')
        call set_nml_var(options%lt%n_nsq_values, n_nsq_values, 'n_nsq_values')
        call set_nml_var(options%lt%n_spd_values, n_spd_values, 'n_spd_values')
        call set_nml_var(options%lt%minimum_layer_size, minimum_layer_size, 'minimum_layer_size')
        call set_nml_var(options%lt%read_LUT, read_LUT, 'read_LUT')
        call set_nml_var(options%lt%write_LUT, write_LUT, 'write_LUT')
        call set_nml_var(options%lt%u_LUT_Filename, u_LUT_Filename, 'u_LUT_Filename')
        call set_nml_var(options%lt%v_LUT_Filename, v_LUT_Filename, 'v_LUT_Filename')
        call set_nml_var(options%lt%overwrite_lt_lut, overwrite_lt_lut, 'overwrite_lt_lut')


    end subroutine lt_parameters_namelist


    !> -------------------------------
    !! Initialize the advection options
    !!
    !! Reads the adv_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine adv_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options
        logical, intent(in), optional  :: info_only, gen_nml

        type(adv_options_type)::adv_options
        integer :: name_unit, rc
        logical :: print_info, gennml

        logical :: boundary_buffer   ! apply some smoothing to the x and y boundaries in MPDATA
        logical :: MPDATA_FCT ! use the flux corrected transport option in MPDATA
        logical :: advect_density
        ! MPDATA order of correction (e.g. 1st=upwind, 2nd=classic, 3rd=better)
        integer :: mpdata_order, flux_corr, h_order, v_order
        
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
        if (options%general%use_adv_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=adv_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'adv_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        call set_nml_var(adv_options%boundary_buffer, boundary_buffer, 'boundary_buffer')
        call set_nml_var(adv_options%MPDATA_FCT, MPDATA_FCT, 'MPDATA_FCT')
        call set_nml_var(adv_options%mpdata_order, mpdata_order, 'mpdata_order')
        call set_nml_var(adv_options%flux_corr, flux_corr, 'flux_corr')
        call set_nml_var(adv_options%h_order, h_order, 'h_order')
        call set_nml_var(adv_options%v_order, v_order, 'v_order')
        call set_nml_var(adv_options%advect_density, advect_density, 'advect_density')

        ! copy the data back into the global options data structure
        options%adv = adv_options
    end subroutine adv_parameters_namelist


    !> -------------------------------
    !! Initialize the PBL options
    !!
    !! Reads the pbl_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine pbl_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options
        logical, intent(in), optional  :: info_only, gen_nml

        type(pbl_options_type)::pbl_options
        integer :: name_unit, rc
        logical :: print_info, gennml

        integer :: ysu_topdown_pblmix ! controls if radiative, top-down mixing in YSU scheme is turned on
        
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
        if (options%general%use_pbl_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=pbl_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'pbl_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        if (.not.(options%physics%radiation == kRA_RRTMG)) ysu_topdown_pblmix = 0

        call set_nml_var(pbl_options%ysu_topdown_pblmix, ysu_topdown_pblmix, 'ysu_topdown_pblmix')

        ! copy the data back into the global options data structure
        options%pbl = pbl_options
    end subroutine pbl_parameters_namelist

    !> -------------------------------
    !! Initialize the surface layer options
    !!
    !! Reads the sfc_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine sfc_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options
        logical, intent(in), optional  :: info_only, gen_nml

        type(sfc_options_type)::sfc_options
        integer :: name_unit, rc, isfflx, scm_force_flux, iz0tlnd, isftcflx
        real    :: sbrlim
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
        if (options%general%use_sfc_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=sfc_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'sfc_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif
        
        call set_nml_var(sfc_options%isfflx, isfflx, 'isfflx')
        call set_nml_var(sfc_options%scm_force_flux, scm_force_flux, 'scm_force_flux')
        call set_nml_var(sfc_options%iz0tlnd, iz0tlnd, 'iz0tlnd')
        call set_nml_var(sfc_options%isftcflx, isftcflx, 'isftcflx')
        call set_nml_var(sfc_options%sbrlim, sbrlim, 'sbrlim')

        ! copy the data back into the global options data structure
        options%sfc = sfc_options
    end subroutine sfc_parameters_namelist


    !> -------------------------------
    !! Initialize the convection scheme options
    !!
    !! Reads the cu_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine cu_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout):: options
        logical, intent(in), optional  :: info_only, gen_nml

        type(cu_options_type) :: cu_options
        integer :: name_unit, rc
        logical :: print_info, gennml

        real :: tendency_fraction, tend_qv_fraction, tend_qc_fraction, tend_th_fraction, tend_qi_fraction
        real :: stochastic_cu


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
        if (options%general%use_cu_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=cu_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'cu_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        ! if not set separately, default to the global tendency setting
        if (tend_qv_fraction < 0) tend_qv_fraction = tendency_fraction
        if (tend_qc_fraction < 0) tend_qc_fraction = tendency_fraction
        if (tend_th_fraction < 0) tend_th_fraction = tendency_fraction
        if (tend_qi_fraction < 0) tend_qi_fraction = tendency_fraction

        ! store everything in the cu_options structure
        call set_nml_var(cu_options%tendency_fraction, tendency_fraction, 'tendency_fraction')
        call set_nml_var(cu_options%tend_qv_fraction, tend_qv_fraction, 'tend_qv_fraction')
        call set_nml_var(cu_options%tend_qc_fraction, tend_qc_fraction, 'tend_qc_fraction')
        call set_nml_var(cu_options%tend_th_fraction, tend_th_fraction, 'tend_th_fraction')
        call set_nml_var(cu_options%tend_qi_fraction, tend_qi_fraction, 'tend_qi_fraction')
        call set_nml_var(cu_options%stochastic_cu, stochastic_cu, 'stochastic_cu')

        ! copy the data back into the global options data structure
        options%cu = cu_options
    end subroutine cu_parameters_namelist


    !> -------------------------------
    !! Initialize the land surface model options
    !!
    !! Reads the lsm_parameters namelist or sets default values
    !!
    !! -------------------------------
    subroutine lsm_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout)::options
        logical, intent(in), optional  :: info_only, gen_nml

        type(lsm_options_type) :: lsm_options
        integer :: name_unit, rc
        logical :: print_info, gennml

        character(len=MAXVARLENGTH) :: LU_Categories ! Category definitions (e.g. USGS, MODIFIED_IGBP_MODIS_NOAH)
        real    :: max_swe
        real    :: snow_den_const                    ! variable for converting snow height into SWE or visa versa when input data is incomplete 

        logical :: monthly_vegfrac                   ! read in 12 months of vegfrac data
        logical :: monthly_albedo                    ! same for albedo (requires vegfrac be monthly)
        real :: update_interval_lsm                  ! minimum number of seconds between LSM updates
        integer :: urban_category                    ! index that defines the urban category in LU_Categories
        integer :: ice_category                      ! index that defines the ice category in LU_Categories
        integer :: water_category                    ! index that defines the water category in LU_Categories
        integer :: sf_urban_phys
        real    :: nmp_soiltstep
        integer :: fsm_nsnow_max, nmp_dveg, nmp_opt_crs, nmp_opt_sfc, nmp_opt_btr, nmp_opt_run, nmp_opt_frz, nmp_opt_inf, nmp_opt_rad, nmp_opt_alb, nmp_opt_snf, nmp_opt_tbot, nmp_opt_stc, nmp_opt_gla, nmp_opt_rsf, nmp_opt_soil, nmp_opt_pedo, nmp_opt_crop, nmp_opt_irr, nmp_opt_irrm, nmp_opt_tdrn, noahmp_output

        integer :: lake_category                    ! index that defines the lake category in (some) LU_Categories

        ! define the namelist
        namelist /lsm_parameters/ LU_Categories, update_interval_lsm, &
                                  urban_category, ice_category, water_category, lake_category, snow_den_const,&
                                  monthly_vegfrac, monthly_albedo, max_swe,  nmp_dveg,   &
                                  nmp_opt_crs, nmp_opt_sfc, nmp_opt_btr, nmp_opt_run, nmp_opt_frz, &
                                  nmp_opt_inf, nmp_opt_rad, nmp_opt_alb, nmp_opt_snf, nmp_opt_tbot,           &
                                  nmp_opt_stc, nmp_opt_gla, nmp_opt_rsf, nmp_opt_soil, nmp_opt_pedo,          &
                                  nmp_opt_crop, nmp_opt_irr, nmp_opt_irrm, nmp_opt_tdrn, nmp_soiltstep,       &
                                  sf_urban_phys, noahmp_output, fsm_nsnow_max !! MJ added

        print_info = .False.
        if (present(info_only)) print_info = info_only

        gennml = .False.
        if (present(gen_nml)) gennml = gen_nml

        call set_nml_var_default(LU_Categories, 'LU_Categories', print_info, gennml)
        call set_nml_var_default(update_interval_lsm, 'update_interval_lsm', print_info, gennml)
        call set_nml_var_default(monthly_vegfrac, 'monthly_vegfrac', print_info, gennml)

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
        call set_nml_var_default(fsm_nsnow_max, 'fsm_nsnow_max', print_info, gennml)

        ! If this is just a verbose print run, exit here so we don't need a namelist
        if (print_info .or. gennml) return

        ! read the namelist options
        if (options%general%use_lsm_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=lsm_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'lsm_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        !If user does not define this option, then let's Do The Right Thing
        if (nmp_opt_sfc == -1) then
            !If user has not turned on the surface layer scheme, then we set this to the recommended default value of 1
            if (options%physics%surfacelayer == 0) then
                nmp_opt_sfc = 1
            !If user has turned on the surface layer scheme, then let's opt to use the surface layer scheme's surface exchange coefficients
            else if (options%physics%surfacelayer == 1) then
                nmp_opt_sfc = 3
            endif
        endif

        call set_default_LU_categories(options, urban_category, ice_category, water_category, LU_Categories, lake_category)

        call set_nml_var(lsm_options%LU_Categories, LU_Categories, 'LU_Categories')
        call set_nml_var(lsm_options%update_interval, update_interval_lsm, 'update_interval_lsm')
        call set_nml_var(lsm_options%monthly_vegfrac, monthly_vegfrac, 'monthly_vegfrac')
        call set_nml_var(lsm_options%monthly_albedo, monthly_albedo, 'monthly_albedo')
        call set_nml_var(lsm_options%urban_category, urban_category, 'urban_category')
        call set_nml_var(lsm_options%ice_category, ice_category, 'ice_category')
        call set_nml_var(lsm_options%water_category, water_category, 'water_category')
        call set_nml_var(lsm_options%lake_category, lake_category, 'lake_category')
        call set_nml_var(lsm_options%snow_den_const, snow_den_const, 'snow_den_const')
        call set_nml_var(lsm_options%max_swe, max_swe, 'max_swe')
        call set_nml_var(lsm_options%sf_urban_phys, sf_urban_phys, 'sf_urban_phys')
        call set_nml_var(lsm_options%nmp_dveg, nmp_dveg, 'nmp_dveg')

        call set_nml_var(lsm_options%nmp_opt_crs, nmp_opt_crs, 'nmp_opt_crs')
        call set_nml_var(lsm_options%nmp_opt_sfc, nmp_opt_sfc, 'nmp_opt_sfc')
        call set_nml_var(lsm_options%nmp_opt_btr, nmp_opt_btr, 'nmp_opt_btr')
        call set_nml_var(lsm_options%nmp_opt_run, nmp_opt_run, 'nmp_opt_run')
        call set_nml_var(lsm_options%nmp_opt_frz, nmp_opt_frz, 'nmp_opt_frz')
        call set_nml_var(lsm_options%nmp_opt_inf, nmp_opt_inf, 'nmp_opt_inf')
        call set_nml_var(lsm_options%nmp_opt_rad, nmp_opt_rad, 'nmp_opt_rad')
        call set_nml_var(lsm_options%nmp_opt_alb, nmp_opt_alb, 'nmp_opt_alb')
        call set_nml_var(lsm_options%nmp_opt_snf, nmp_opt_snf, 'nmp_opt_snf')
        call set_nml_var(lsm_options%nmp_opt_tbot, nmp_opt_tbot, 'nmp_opt_tbot')
        call set_nml_var(lsm_options%nmp_opt_stc, nmp_opt_stc, 'nmp_opt_stc')
        call set_nml_var(lsm_options%nmp_opt_gla, nmp_opt_gla, 'nmp_opt_gla')
        call set_nml_var(lsm_options%nmp_opt_rsf, nmp_opt_rsf, 'nmp_opt_rsf')
        call set_nml_var(lsm_options%nmp_opt_soil, nmp_opt_soil, 'nmp_opt_soil')
        call set_nml_var(lsm_options%nmp_opt_pedo, nmp_opt_pedo, 'nmp_opt_pedo')
        call set_nml_var(lsm_options%nmp_opt_crop, nmp_opt_crop, 'nmp_opt_crop')
        call set_nml_var(lsm_options%nmp_opt_irr, nmp_opt_irr, 'nmp_opt_irr')
        call set_nml_var(lsm_options%nmp_opt_irrm, nmp_opt_irrm, 'nmp_opt_irrm')
        call set_nml_var(lsm_options%nmp_opt_tdrn, nmp_opt_tdrn, 'nmp_opt_tdrn')
        call set_nml_var(lsm_options%nmp_soiltstep, nmp_soiltstep, 'nmp_soiltstep')
        call set_nml_var(lsm_options%noahmp_output, noahmp_output, 'noahmp_output')
        call set_nml_var(lsm_options%fsm_nsnow_max, fsm_nsnow_max, 'fsm_nsnow_max')

        
        ! copy the data back into the global options data structure
        options%lsm = lsm_options

        ! Change z length of snow arrays here, since we need to change their size for the output arrays, which are set in
        ! output_options_namelist
        if (options%physics%snowmodel==kSM_FSM) then
                kSNOW_GRID_Z = options%lsm%fsm_nsnow_max
                kSNOWSOIL_GRID_Z = kSNOW_GRID_Z+kSOIL_GRID_Z
        endif

    end subroutine lsm_parameters_namelist


    !> -------------------------------
    !! Initialize the radiation model options
    !!
    !! Reads the rad_parameters namelist or sets default values
    !! -------------------------------
    subroutine rad_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),   intent(in)   :: filename
        type(options_t),    intent(inout)::options
        logical, intent(in), optional  :: info_only, gen_nml

        type(rad_options_type) :: rad_options
        integer :: name_unit, rc
        logical :: print_info, gennml

        real    :: update_interval_rrtmg             ! minimum number of seconds between RRTMG updates
        integer :: icloud                            ! how RRTMG interacts with clouds
        integer :: cldovrlp                          ! how RRTMG considers cloud overlapping
        logical :: read_ghg
        real    :: tzone !! MJ adedd,tzone is UTC Offset and 1 here for centeral Erupe
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
        if (options%general%use_rad_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit, iostat=rc, nml=rad_parameters)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'rad_parameters' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        call set_nml_var(rad_options%update_interval_rrtmg, update_interval_rrtmg, 'update_interval_rrtmg')
        call set_nml_var(rad_options%icloud, icloud, 'icloud')
        call set_nml_var(rad_options%cldovrlp, cldovrlp, 'cldovrlp')
        call set_nml_var(rad_options%read_ghg, read_ghg, 'read_ghg')
        call set_nml_var(rad_options%tzone, tzone, 'tzone')

        ! copy the data back into the global options data structure
        options%rad = rad_options
        
    end subroutine rad_parameters_namelist
    
    
    subroutine wind_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(options_t),    intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml
        
        integer :: name_unit, update_frequency             ! logical unit number for namelist
        logical :: print_info, gennml
        !Define parameters
        integer :: wind_iterations, rc
        logical :: Sx, wind_only, thermal
        real    :: Sx_dmax, Sx_scale_ang, TPI_scale, TPI_dmax, alpha_const, smooth_wind_distance
        
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
        if (options%general%use_wind_options) then
            open(io_newunit(name_unit), file=filename)
            read(name_unit,iostat=rc,nml=wind)
            close(name_unit)
            if (rc /= 0) then
                if (STD_OUT_PE) write(*,*) "--------------------------------"
                if (STD_OUT_PE) write(*,*) "Error reading 'wind' namelist, continuing with defaults"
                if (STD_OUT_PE) write(*,*) "--------------------------------"
            endif
        endif

        !Perform checks
        if (smooth_wind_distance.eq.(-9999)) then
            smooth_wind_distance=options%domain%dx*2
            if (STD_OUT_PE) write(*,*) " Default smoothing distance = dx*2 = ", smooth_wind_distance
        elseif (smooth_wind_distance<0) then
            write(*,*) " Wind smoothing must be a positive number"
            write(*,*) " smooth_wind_distance = ",smooth_wind_distance
            smooth_wind_distance = options%domain%dx*2
        endif

        call set_nml_var(options%wind%smooth_wind_distance, smooth_wind_distance, 'smooth_wind_distance')
        call set_nml_var(options%wind%Sx, Sx, 'Sx')
        call set_nml_var(options%wind%thermal, thermal, 'thermal')
        call set_nml_var(options%wind%wind_only, wind_only, 'wind_only')
        call set_nml_var(options%wind%Sx_dmax, Sx_dmax, 'Sx_dmax')
        call set_nml_var(options%wind%TPI_dmax, TPI_dmax, 'TPI_dmax')
        call set_nml_var(options%wind%TPI_scale, TPI_scale, 'TPI_scale')
        call set_nml_var(options%wind%Sx_scale_ang, Sx_scale_ang, 'Sx_scale_ang')
        call set_nml_var(options%wind%alpha_const, alpha_const, 'alpha_const')
        call set_nml_var(options%wind%wind_iterations, wind_iterations, 'wind_iterations')

        call options%wind%update_dt%set(seconds=options%forcing%input_dt%seconds()/update_frequency)

    end subroutine wind_namelist


    subroutine time_parameters_namelist(filename, options, info_only, gen_nml)
        implicit none
        character(len=*),             intent(in)    :: filename
        type(options_t),              intent(inout) :: options
        logical, intent(in), optional  :: info_only, gen_nml
        
        integer :: name_unit, rc                            ! logical unit number for namelist
        !Define parameters
        real :: cfl_reduction_factor    
        logical :: RK3
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
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            if (STD_OUT_PE) write(*,*) "Error reading 'time_parameters' namelist"
            if (STD_OUT_PE) write(*,*) "--------------------------------"
            stop
        endif
                
        call set_nml_var(options%time%cfl_reduction_factor, cfl_reduction_factor, 'cfl_reduction_factor')
        call set_nml_var(options%time%RK3, RK3, 'RK3')

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
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_FILES) :: forcing_files
        integer :: nfiles
        integer :: file_unit
        integer :: i, error
        logical :: first_file_exists, last_file_exists
        character(len=MAXFILELENGTH) :: temporary_file

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
        if (STD_OUT_PE) write(*,*) "  Boundary conditions files to be used:"
        if (nfiles>10) then
            if (STD_OUT_PE) write(*,*) "    nfiles=", trim(str(nfiles)), ", too many to print."
            if (STD_OUT_PE) write(*,*) "    First file:", trim(forcing_files(1))
            if (STD_OUT_PE) write(*,*) "    Last file: ", trim(forcing_files(nfiles))
        else
            do i=1,nfiles
                if (STD_OUT_PE) write(*,*) "      ",trim(forcing_files(i))
            enddo
        endif

        ! Check that the options file actually exists
        INQUIRE(file=trim(forcing_files(1)), exist=first_file_exists)
        INQUIRE(file=trim(forcing_files(nfiles)), exist=last_file_exists)

        ! if options file does not exist, print an error and quit
        if (.not.first_file_exists .or. .not.last_file_exists) then
            if (.not.first_file_exists .and. STD_OUT_PE) write(*,*) "The first forcing file does not exist = ", trim(forcing_files(1))
            if (.not.last_file_exists .and. STD_OUT_PE) write(*,*) "The last forcing file does not exist = ", trim(forcing_files(nfiles))

            stop "At least the first or last forcing file does not exist. Only the first and last file are checked, please check the rest."
        endif


    end function read_forcing_file_names


    !> -------------------------------
    !! Sets the default value for each of three land use categories depending on the LU_Categories input
    !!
    !! -------------------------------
    subroutine set_default_LU_categories(options, urban_category, ice_category, water_category, LU_Categories, lake_category)
        ! if various LU categories were not defined in the namelist (i.e. they == -1) then attempt
        ! to define default values for them based on the LU_Categories variable supplied.
        implicit none
        type(options_t),    intent(inout)::options
        integer, intent(inout) :: urban_category, ice_category, water_category, lake_category
        character(len=MAXVARLENGTH), intent(in) :: LU_Categories

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
            if((options%physics%watersurface==kWATER_LAKE) .AND. (STD_OUT_PE)) then
                write(*,*) "WARNING: Lake model selected (water=2), but USGS LU-categories has no lake category"
            endif

        elseif (trim(LU_Categories)=="USGS-RUC") then
            if (urban_category==-1) urban_category = 1
            if (ice_category==-1)   ice_category = 24
            if (water_category==-1) water_category = 16
            if (lake_category==-1) lake_category = 28
            ! also note, lakes_category = 28
            ! write(*,*) "WARNING: not handling lake category (28)"

        elseif (trim(LU_Categories)=="MODI-RUC") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15
            if (water_category==-1) water_category = 17
            if (lake_category==-1) lake_category = 21
            ! also note, lakes_category = 21
            ! write(*,*) "WARNING: not handling lake category (21)"

        elseif (trim(LU_Categories)=="NLCD40") then
            if (urban_category==-1) urban_category = 13
            if (ice_category==-1)   ice_category = 15 ! and 22?
            ! if (water_category==-1) water_category = 17 ! and 21 'Open Water'
            if(options%physics%watersurface==kWATER_LAKE) write(*,*) "WARNING: Lake model selected (water=2), but NLCD40 LU-categories has no lake category"
            write(*,*) "WARNING: not handling all varients of categories (e.g. permanent_snow=15 is, but permanent_snow_ice=22 is not)"
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
                if (STD_OUT_PE) write(*,*) "WARNING: trying to add var outside of permitted list:",varids(i), size(varlist)
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
            if (STD_OUT_PE) write(*,*) "Model version does not match namelist version"
            if (STD_OUT_PE) write(*,*) "  Model version: ",kVERSION_STRING
            if (STD_OUT_PE) write(*,*) "  Namelist version: ",trim(options%version)
            call print_model_diffs(options%version)
            stop
        endif

        if (STD_OUT_PE) write(*,*) "  Model version: ",trim(options%version)

    end subroutine version_check


end submodule
