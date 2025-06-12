!>------------------------------------------------------------
!!  Implementation of domain object
!!
!!  implements all domain type bound procedures
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(domain_interface) domain_implementation
    use assertions_mod,       only : assert, assertions
    use mod_atm_utilities,    only : exner_function, update_pressure, compute_ivt, compute_iq
    use icar_constants
    use string,               only : str
    use io_routines,          only : io_write, io_read
    use geo,                  only : geo_lut, geo_interp, geo_interp2d, standardize_coordinates
    use array_utilities,      only : array_offset_x, array_offset_y, smooth_array, smooth_array_2d, make_2d_x, make_2d_y
    use vertical_interpolation,only : vinterp, vLUT
    use output_metadata,            only : get_varname, get_varmeta, get_varindx
    use mod_wrf_constants,    only : gravity, R_d, KARMAN
    use iso_fortran_env

    implicit none

    interface setup
        module procedure setup_var
    end interface
    
    real, parameter::deg2rad=0.017453293 !2*pi/360
    integer :: FILTER_WIDTH = 7

    ! primary public routines : init, get_initial_conditions, halo_send, halo_retrieve, or halo_exchange
contains


    !> -------------------------------
    !! Initialize the size of the domain
    !!
    !! -------------------------------
    module subroutine init(this, options, nest_indx)
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options
        integer, intent(in) :: nest_indx
        
        ! Call the parent type's init procedure
        call this%init_flow_obj(options, nest_indx)

        this%dx = options%domain%dx

        call read_domain_shape(this, options)
        
        call create_variables(this, options)
        
        call initialize_core_variables(this, options)  ! split into several subroutines?

        call read_land_variables(this, options)

        call set_var_lists(this, options)

        call init_relax_filters(this,options)

        call init_batch_exch(this)

    end subroutine

    module subroutine init_batch_exch(this)
        implicit none
        class(domain_t), intent(inout) :: this

        type(variable_t) :: tmp_var
        integer :: i


        ! pack adv_vars with the variables to advect
        do i = 1,size(this%adv_vars)
            tmp_var = get_varmeta(get_varindx(this%adv_vars(i)%n))
            if (tmp_var%three_d) then
                this%n_adv_3d = this%n_adv_3d + 1
            end if
        end do
        this%n_adv_2d = size(this%adv_vars) - this%n_adv_3d

        ! pack adv_vars with the variables to advect
        do i = 1,size(this%exch_vars)
            tmp_var = get_varmeta(get_varindx(this%exch_vars(i)%n))
            if (tmp_var%three_d) then
                this%n_exch_3d = this%n_exch_3d + 1
            end if
        end do
        this%n_exch_2d = size(this%exch_vars) - this%n_exch_3d

        call this%halo%init(this%exch_vars, this%adv_vars, this%grid, this%compute_comms)

    end subroutine init_batch_exch

    module subroutine batch_exch(this, two_d, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: two_d, exch_only


        if (present(two_d)) then
            if (two_d) then
                call this%halo_2d_send_batch()
                call this%halo_2d_retrieve_batch()    
            endif
        else
            if (present(exch_only)) then
                call this%halo_3d_send_batch(exch_only=exch_only)
                call this%halo_3d_retrieve_batch(exch_only=exch_only)
            else
                call this%halo_3d_send_batch()
                call this%halo_3d_retrieve_batch()    
            end if
        end if

    end subroutine batch_exch

    module subroutine halo_3d_send_batch(this, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: exch_only


        if (present(exch_only)) then
            call this%halo%halo_3d_send_batch(this%exch_vars, this%adv_vars, this%vars_3d, exch_var_only=exch_only)
        else
            call this%halo%halo_3d_send_batch(this%exch_vars, this%adv_vars, this%vars_3d)
        end if

    end subroutine halo_3d_send_batch

    module subroutine halo_3d_retrieve_batch(this, exch_only)
        implicit none
        class(domain_t), intent(inout) :: this
        logical, optional,   intent(in) :: exch_only

        if (present(exch_only)) then
            call this%halo%halo_3d_retrieve_batch(this%exch_vars, this%adv_vars, this%vars_3d, exch_var_only=exch_only)
        else
            call this%halo%halo_3d_retrieve_batch(this%exch_vars, this%adv_vars, this%vars_3d)
        end if

    end subroutine halo_3d_retrieve_batch

    module subroutine halo_2d_send_batch(this)
        implicit none
        class(domain_t), intent(inout) :: this

        if (this%n_exch_2d + this%n_adv_2d == 0) return

        call this%halo%halo_2d_send_batch(this%exch_vars, this%adv_vars, this%vars_2d)

    end subroutine halo_2d_send_batch

    module subroutine halo_2d_retrieve_batch(this)
        implicit none
        class(domain_t), intent(inout) :: this

        if (this%n_exch_2d + this%n_adv_2d == 0) return

        call this%halo%halo_2d_retrieve_batch(this%exch_vars, this%adv_vars, this%vars_2d)

    end subroutine halo_2d_retrieve_batch

    !> -------------------------------
    !! Release memory associated with the domain
    !!
    !! -------------------------------
    module subroutine release(this)
        class(domain_t), intent(inout) :: this

        call this%halo%finalize()

    end subroutine
    
    subroutine set_var_lists(this, options)
        class(domain_t), intent(inout) :: this
        type(options_t), intent(in)    :: options

        type(variable_t) :: tmp_var
        integer :: var_list(kMAX_STORAGE_VARS), i, n_vars, var_indx

        !Advection variables -- these are exchanged AND advected
        n_vars = 0
        do i = 1, size(kADV_VARS)
            var_indx = get_varindx(trim(kADV_VARS(i)))
            if (this%var_indx(var_indx)%v > 0) n_vars = n_vars + 1
        enddo
        allocate(this%adv_vars(n_vars))
        n_vars = 0

        do i = 1, size(kADV_VARS)
            var_indx = get_varindx(trim(kADV_VARS(i)))
            if (this%var_indx(var_indx)%v > 0) then
                n_vars = n_vars + 1
                this%adv_vars(n_vars)%n = this%var_indx(var_indx)%n
                this%adv_vars(n_vars)%v = this%var_indx(var_indx)%v
            endif
        enddo

        n_vars = 0
        do i = 1, size(kEXCH_VARS)
            var_indx = get_varindx(trim(kEXCH_VARS(i)))
            if (this%var_indx(var_indx)%v > 0) n_vars = n_vars + 1
        enddo
        allocate(this%exch_vars(n_vars))
        n_vars = 0

        do i = 1, size(kEXCH_VARS)
            var_indx = get_varindx(trim(kEXCH_VARS(i)))
            if (this%var_indx(var_indx)%v > 0) then
                n_vars = n_vars + 1
                this%exch_vars(n_vars)%n = this%var_indx(var_indx)%n
                this%exch_vars(n_vars)%v = this%var_indx(var_indx)%v
            endif
        enddo

        var_list = options%output%vars_for_output + options%vars_for_restart
        ! call sort_by_kVARS(var_list)
        n_vars = 0
        do i = 1, kMAX_STORAGE_VARS
            if (0<var_list(i)) then
                ! get the variable meta data defined in var_defs.F90
                tmp_var = get_varmeta(i)
                
                if (tmp_var%name == "") cycle
                n_vars = n_vars + 1
                this%vars_to_out(i)%n = this%var_indx(i)%n
                this%vars_to_out(i)%v = this%var_indx(i)%v
            endif
        enddo

        ! Sort the output vars according to their ordering in kVARS. This lets the
        ! above lines be in any arbitrary order.
        ! call this%vars_to_out%sort_by_kVARS()
    end subroutine set_var_lists

    !> -------------------------------
    !! Set up the initial conditions for the domain
    !!
    !! This includes setting up all of the geographic interpolation now that we have the forcing grid
    !! and interpolating the first time step of forcing data on to the high res domain grids
    !!
    !! -------------------------------
    module subroutine get_initial_conditions(this, forcing, options)
      implicit none
      class(domain_t),  intent(inout) :: this
      type(boundary_t), intent(inout) :: forcing
      type(options_t),  intent(in)    :: options

      integer :: i

      ! create geographic lookup table for domain
      call setup_geo_interpolation(this, forcing, options)

        ! for all variables with a forcing_var /= "", get forcing, interpolate to local domain
      call this%interpolate_forcing(forcing)

      call this%enforce_limits()

    end subroutine


    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    module subroutine diagnostic_update(this, forcing_update)
        implicit none
        class(domain_t),  intent(inout)   :: this
        logical, intent(in), optional    :: forcing_update
        integer :: i, j, k
        real :: qsum
        logical :: forcing_update_only
        real, dimension(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme) :: mod_temp_3d
        real, dimension(this%ims:this%ime, this%jms:this%jme) :: surf_temp_1, surf_temp_2, surf_temp_3

        forcing_update_only = .False.
        if (present(forcing_update)) forcing_update_only = forcing_update

        associate(ims => this%ims, ime => this%ime,                             &
                  jms => this%jms, jme => this%jme,                             &
                  ids => this%ids, ide => this%ide,                             &
                  jds => this%jds, jde => this%jde,                             &
                  kms => this%kms, kme => this%kme,                             &
                  its => this%its, ite => this%ite,                             &
                  jts => this%jts, jte => this%jte,                             &
                  exner                 => this%vars_3d(this%var_indx(kVARS%exner)%v)%data_3d,                  &
                  pressure              => this%vars_3d(this%var_indx(kVARS%pressure)%v)%data_3d,               &
                  pressure_i            => this%vars_3d(this%var_indx(kVARS%pressure_interface)%v)%data_3d,     &
                  dz_i                  => this%vars_3d(this%var_indx(kVARS%dz_interface)%v)%data_3d,           &
                  dz_mass               => this%vars_3d(this%var_indx(kVARS%dz)%v)%data_3d,                &
                  psfc                  => this%vars_2d(this%var_indx(kVARS%surface_pressure)%v)%data_2d,       &
                  density               => this%vars_3d(this%var_indx(kVARS%density)%v)%data_3d,                &
                  temperature           => this%vars_3d(this%var_indx(kVARS%temperature)%v)%data_3d,            &
                  qv                    => this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%data_3d,            &  
                  temperature_i         => this%vars_3d(this%var_indx(kVARS%temperature_interface)%v)%data_3d,  &
                  u                     => this%vars_3d(this%var_indx(kVARS%u)%v)%data_3d,                      &
                  v                     => this%vars_3d(this%var_indx(kVARS%v)%v)%data_3d,                      &
                  u_mass                => this%vars_3d(this%var_indx(kVARS%u_mass)%v)%data_3d,                 &
                  v_mass                => this%vars_3d(this%var_indx(kVARS%v_mass)%v)%data_3d,                 &
                  potential_temperature => this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%data_3d )

                !Calculation of density
                if (forcing_update_only) then
                    exner = exner_function(this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d)
                    temperature = potential_temperature * exner
                    density =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d / (R_d * temperature*(1+qv)) ! kg/m^3
                else
                    exner = exner_function(pressure)
                    temperature = potential_temperature * exner
                    density =  pressure / (R_d * temperature*(1+qv)) ! kg/m^3
                endif
        

        ! differences between forcing data at the boudnary and the internal model state can lead to strong discontinuities in temperature and qv
        ! these then affect density, leading to discontinuities in density, and thus winds. So, here we set the density for points on the "frame"
        ! to be the same as the density in the first cell within the physics region of the domain
                if (ims==ids) then
                    do j = jms,jme
                        do k = kms,kme
                            do i = ims,its-1
                                temperature(i,k,j) = potential_temperature(its,k,j) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,j))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,j))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif
                if (ime==ide) then
                    do j = jms,jme
                        do k = kms,kme
                            do i = ite+1,ime
                                temperature(i,k,j) = potential_temperature(ite,k,j) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,j))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,j))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif
                if (jms==jds) then
                    do j = jms,jts-1
                        do k = kms,kme
                            do i = ims,ime
                                temperature(i,k,j) = potential_temperature(i,k,jts) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(i,k,jts))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(i,k,jts))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif
                if (jme==jde) then
                    do j = jte+1,jme
                        do k = kms,kme
                            do i = ims,ime
                                temperature(i,k,j) = potential_temperature(i,k,jte) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(i,k,jte))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(i,k,jte))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif

                if (ims==ids .and. jms==jds) then
                    do j = jms,jts-1
                        do k = kms,kme
                            do i = ims,its-1
                                temperature(i,k,j) = potential_temperature(its,k,jts) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,jts))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,jts))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif

                if (ims==ids .and. jme==jde) then
                    do j = jte+1,jme
                        do k = kms,kme
                            do i = ims,its-1
                                temperature(i,k,j) = potential_temperature(its,k,jte) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,jte))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(its,k,jte))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif

                if (ime==ide .and. jms==jds) then
                    do j = jms,jts-1
                        do k = kms,kme
                            do i = ite+1,ime
                                temperature(i,k,j) = potential_temperature(ite,k,jts) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,jts))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,jts))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif

                if (ime==ide .and. jme==jde) then
                    do j = jte+1,jme
                        do k = kms,kme
                            do i = ite+1,ime
                                temperature(i,k,j) = potential_temperature(ite,k,jte) * exner(i,k,j)
                                if (forcing_update_only) then
                                    density(i,k,j) =  this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,jte))) ! kg/m^3
                                else
                                    density(i,k,j) =  pressure(i,k,j) / (R_d * temperature(i,k,j)*(1+qv(ite,k,jte))) ! kg/m^3
                                endif
                            enddo
                        enddo
                    enddo
                endif
        
                if (forcing_update_only) return

        
        temperature_i(ims:ime,kms,jms:jme) = temperature(ims:ime,kms,jms:jme) + (temperature(ims:ime,kms,jms:jme) - temperature(ims:ime,kms+1,jms:jme)) * 0.5
        pressure_i(ims:ime,kms,jms:jme) = pressure(ims:ime,kms,jms:jme) + (pressure(ims:ime,kms,jms:jme) - pressure(ims:ime,kms+1,jms:jme)) * 0.5

        do j = jms,jme
            do k = kms+1,kme
                do i = ims,ime
                    pressure_i(i,k,j) = (dz_i(i,k,j)*pressure(i,k-1,j)+dz_i(i,k-1,j)*pressure(i,k,j))/((dz_i(i,k-1,j)+dz_i(i,k,j)))
                    temperature_i(i,k,j) = (dz_i(i,k,j)*temperature(i,k-1,j)+dz_i(i,k-1,j)*temperature(i,k,j))/((dz_i(i,k-1,j)+dz_i(i,k,j)))
                enddo
            enddo
        enddo
        temperature_i(ims:ime,kme+1,jms:jme) = temperature(ims:ime,kme,jms:jme) + (temperature(ims:ime,kme,jms:jme) - temperature(ims:ime,kme-1,jms:jme)) * 0.5
        pressure_i(ims:ime,kme+1,jms:jme) = pressure(ims:ime,kme,jms:jme) + (pressure(ims:ime,kme,jms:jme) - pressure(ims:ime,kme-1,jms:jme)) * 0.5

        if (this%var_indx(kVARS%u_mass)%v > 0) then
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        u_mass(i,k,j) = (u(i+1,k,j) + u(i,k,j)) * 0.5
                        v_mass(i,k,j) = (v(i,k,j+1) + v(i,k,j)) * 0.5
                    enddo
                enddo
            enddo
        endif
                
        if (this%var_indx(kVARS%surface_pressure)%v > 0) then
            psfc = pressure_i(ims:ime, kms, jms:jme)
        endif
        if (this%var_indx(kVARS%ivt)%v > 0) then
            call compute_ivt(this%vars_2d(this%var_indx(kVARS%ivt)%v)%data_2d, qv, u_mass, v_mass, pressure_i(:,kms:kme,:))
        endif
        if (this%var_indx(kVARS%iwv)%v > 0) then
            call compute_iq(this%vars_2d(this%var_indx(kVARS%iwv)%v)%data_2d, qv, pressure_i(:,kms:kme,:))
        endif
        if (this%var_indx(kVARS%iwl)%v > 0) then
            mod_temp_3d = 0
            if (this%var_indx(kVARS%cloud_water_mass)%v > 0) mod_temp_3d = mod_temp_3d + this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%data_3d(i,k,j)
            if (this%var_indx(kVARS%rain_mass)%v > 0) mod_temp_3d = mod_temp_3d + this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%data_3d(i,k,j)
            call compute_iq(this%vars_2d(this%var_indx(kVARS%iwl)%v)%data_2d, mod_temp_3d, pressure_i(:,kms:kme,:))
        endif
        if (this%var_indx(kVARS%iwi)%v > 0) then
            mod_temp_3d = 0
            if (this%var_indx(kVARS%ice_mass)%v > 0) mod_temp_3d = mod_temp_3d + this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%data_3d(i,k,j)
            if (this%var_indx(kVARS%snow_mass)%v > 0) mod_temp_3d = mod_temp_3d + this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%data_3d(i,k,j)
            if (this%var_indx(kVARS%graupel_mass)%v > 0) mod_temp_3d = mod_temp_3d + this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%data_3d(i,k,j)
            call compute_iq(this%vars_2d(this%var_indx(kVARS%iwi)%v)%data_2d, mod_temp_3d, pressure_i(:,kms:kme,:))
        endif
        
        ! temporary constant
        if (this%var_indx(kVARS%roughness_z0)%v > 0) then
            ! use log-law of the wall to convert from first model level to surface
            surf_temp_1 = karman / log((this%vars_3d(this%var_indx(kVARS%z)%v)%data_3d(ims:ime,kms,jms:jme) - this%vars_2d(this%var_indx(kVARS%terrain)%v)%data_2d(ims:ime,jms:jme)) / this%vars_2d(this%var_indx(kVARS%roughness_z0)%v)%data_2d(ims:ime,jms:jme))
            ! use log-law of the wall to convert from surface to 10m height
            surf_temp_2 = log(10.0 / this%vars_2d(this%var_indx(kVARS%roughness_z0)%v)%data_2d(ims:ime,jms:jme)) / karman
        endif

        if (this%var_indx(kVARS%u_10m)%v > 0) then
            surf_temp_3                         = u_mass      (ims:ime,kms,jms:jme) * surf_temp_1
            this%vars_2d(this%var_indx(kVARS%u_10m)%v)%data_2d(ims:ime,jms:jme) = surf_temp_3     * surf_temp_2
            surf_temp_3                         = v_mass      (ims:ime,kms,jms:jme) * surf_temp_1
            this%vars_2d(this%var_indx(kVARS%v_10m)%v)%data_2d(ims:ime,jms:jme) = surf_temp_3     * surf_temp_2
        endif

        if (this%var_indx(kVARS%ustar)%v > 0) then
            ! now calculate master ustar based on U and V combined in quadrature
            this%vars_2d(this%var_indx(kVARS%ustar)%v)%data_2d(its:ite,jts:jte) = sqrt(u_mass(its:ite,kms,jts:jte)**2 + v_mass(its:ite,kms,jts:jte)**2) * surf_temp_1(its:ite,jts:jte)
        endif
        
        end associate

    end subroutine diagnostic_update

   
    !> -------------------------------
    !! Allocate and or initialize all domain variables if they have been requested
    !!
    !! -------------------------------
    subroutine create_variables(this, opt)
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: opt
        integer :: i,j

        integer :: ims, ime, jms, jme, kms, kme
        type(variable_t) :: tmp_var, forcing_var_hi
        type(grid_t)     :: grid
        integer :: n_dims, n_one_d, n_two_d, n_three_d, n_four_d, n_forcing_var

        if (STD_OUT_PE) print *,"  Initializing variables"
        if (STD_OUT_PE) flush(output_unit)
        n_one_d = 0
        n_two_d = 0
        n_three_d = 0
        n_four_d = 0
        n_forcing_var = 0
        ! get counts of 1d, 2d, 3d, and 4d variables
        do i = 1, kMAX_STORAGE_VARS
            if (0<opt%vars_to_allocate(i)) then
                ! get the variable meta data defined in var_defs.F90
                tmp_var = get_varmeta(i,opt)

                if (tmp_var%one_d) then
                    n_one_d = n_one_d + 1
                else if (tmp_var%two_d) then
                    n_two_d = n_two_d + 1
                else if (tmp_var%three_d) then
                    n_three_d = n_three_d + 1
                else if (tmp_var%four_d) then
                    n_four_d = n_four_d + 1
                endif
                                        
                if (tmp_var%forcing_var /= "") then
                    n_forcing_var = n_forcing_var + 1
                endif

            endif
        enddo

        allocate(this%vars_1d(n_one_d))
        allocate(this%vars_2d(n_two_d))
        allocate(this%vars_3d(n_three_d))
        allocate(this%vars_4d(n_four_d))
        allocate(this%forcing_hi(n_forcing_var))

        !reset to be used as index counters
        n_one_d = 0
        n_two_d = 0
        n_three_d = 0
        n_four_d = 0
        n_forcing_var = 0
        do i = 1, kMAX_STORAGE_VARS
            if (0<opt%vars_to_allocate(i)) then
                ! get the variable meta data defined in var_defs.F90
                tmp_var = get_varmeta(i,opt)
                
                ! test if variable even has an entry
                if (tmp_var%name == "") cycle 

                n_dims = size(tmp_var%dimensions)

                if (tmp_var%one_d) then
                    grid = this%column_grid
                else if (tmp_var%two_d) then
                    if (tmp_var%dimensions(1) == "lon_x_global") then
                        grid = this%global_grid_2d
                    elseif (tmp_var%dimensions(1) == "lon_x_neighbor") then
                        grid = this%neighbor_grid_2d
                    else
                        grid = this%grid2d
                        if (tmp_var%xstag > 0) grid = this%u_grid2d
                        if (tmp_var%ystag > 0) grid = this%v_grid2d
                    endif
                else if (tmp_var%three_d) then
                    if (tmp_var%dimensions(1) == "lon_x_global") then
                        if (tmp_var%dimensions(2) == "level_i") then
                            grid = this%global_grid8w
                        else if (tmp_var%dimensions(2) == "level") then
                            grid = this%global_grid
                        endif
                    elseif (tmp_var%dimensions(1) == "lon_x_neighbor") then
                        if (tmp_var%dimensions(2) == "level_i") then
                            grid = this%neighbor_grid8w
                        else if (tmp_var%dimensions(2) == "level") then
                            grid = this%neighbor_grid
                        endif
                    else
                        select case (tmp_var%dimensions(2))
                            case ("level")
                                grid = this%grid
                                if (tmp_var%xstag > 0) grid = this%u_grid
                                if (tmp_var%ystag > 0) grid = this%v_grid
                            case ("level_i")
                                grid = this%grid8w
                            case ("nsoil")
                                grid = this%grid_soil
                            case ("nsnow")
                                grid = this%grid_snow
                            case ("nsnowsoil")
                                grid = this%grid_snowsoil
                            case ("nsoil_composition")
                                grid = this%grid_soilcomp
                            case ("crop")
                                grid = this%grid_croptype
                            case ("gecros")
                                grid = this%grid_gecros
                            case ("month")
                                grid = this%grid_monthly
                            case ("nlevlake")
                                grid = this%grid_lake
                            case ("nlevsoisno")
                                grid = this%grid_lake_soisno
                            case ("nlevsoisno_1")
                                grid = this%grid_lake_soisno_1
                            case ("nlevsoi_lake")
                                grid = this%grid_lake_soi
                            case ("azimuth")
                                grid = this%grid_hlm
                        end select
                    endif
                else if (tmp_var%four_d) then
                    grid = this%grid_Sx
                endif            

                ! if we are using the linear wind solver, we need dz and z information for the global grid. 
                ! but to save memory, global_z/dz is setup by default for the neighbor grid. Change here
                if ( (opt%physics%windtype == kWIND_LINEAR) .or. (opt%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
                     (opt%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                    if (i==kVARS%h1 .or. i==kVARS%h2) then
                        grid = this%global_grid_2d
                    else if (i==kVARS%global_z_interface) then
                        grid = this%global_grid8w
                    else if (i==kVARS%global_dz_interface) then
                        grid = this%global_grid
                    endif
                endif

                ! test if forcing var is empty
                if (trim(tmp_var%forcing_var) == '') then
                    call setup(tmp_var, grid, dtype=tmp_var%dtype)
                else
                    call setup(tmp_var, grid, forcing_var=trim(tmp_var%forcing_var), force_boundaries=tmp_var%force_boundaries, dtype=tmp_var%dtype)

                    n_forcing_var = n_forcing_var + 1
                    call forcing_var_hi%initialize( grid, forcing_var=trim(tmp_var%forcing_var))
                    forcing_var_hi%name = trim(tmp_var%name)

                    this%forcing_var_indx(i)%n = trim(tmp_var%name)
                    this%forcing_var_indx(i)%v = n_forcing_var

                    this%forcing_hi(this%forcing_var_indx(i)%v) = forcing_var_hi
                endif

                this%var_indx(i)%n = trim(tmp_var%name)
                if (tmp_var%one_d) then
                    n_one_d = n_one_d + 1
                    this%var_indx(i)%v = n_one_d
                    this%vars_1d(this%var_indx(i)%v) = tmp_var
                else if (tmp_var%two_d) then
                    n_two_d = n_two_d + 1
                    this%var_indx(i)%v = n_two_d
                    this%vars_2d(this%var_indx(i)%v) = tmp_var
                else if (tmp_var%three_d) then
                    n_three_d = n_three_d + 1
                    this%var_indx(i)%v = n_three_d
                    this%vars_3d(this%var_indx(i)%v) = tmp_var
                else if (tmp_var%four_d) then
                    n_four_d = n_four_d + 1
                    this%var_indx(i)%v = n_four_d
                    this%vars_4d(this%var_indx(i)%v) = tmp_var
                endif
            endif
        enddo


        if (0<opt%vars_to_allocate( kVARS%tend_qv_adv) )                allocate(this%tend%qv_adv(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),   source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qv_pbl) )                allocate(this%tend%qv_pbl(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),   source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qv) )                    allocate(this%tend%qv(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_th) )                    allocate(this%tend%th(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_th_pbl) )                allocate(this%tend%th_pbl(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qc_pbl) )                allocate(this%tend%qc_pbl(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qi_pbl) )                allocate(this%tend%qi_pbl(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qc) )                    allocate(this%tend%qc(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qi) )                    allocate(this%tend%qi(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qs) )                    allocate(this%tend%qs(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_qr) )                    allocate(this%tend%qr(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),       source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_u) )                     allocate(this%tend%u(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),        source=0.0)
        if (0<opt%vars_to_allocate( kVARS%tend_v) )                     allocate(this%tend%v(this%ims:this%ime, this%kms:this%kme, this%jms:this%jme),        source=0.0)

    end subroutine

    !> -------------------------------
    !! Setup a regular variable.
    !!
    !! Initializes the variable
    !! including the forcing_variable if it was set
    !! and adds that variable to the list of variables that has forcing data if the list is supplied
    !! and the forcing_var is both present and not blank ("")
    !!
    !! -------------------------------
    subroutine setup_var(var, grid, forcing_var, force_boundaries, dtype)
        implicit none
        type(variable_t),   intent(inout) :: var
        type(grid_t),       intent(in)    :: grid
        character(len=*),   intent(in),   optional :: forcing_var
        logical,            intent(in),   optional :: force_boundaries
        integer,            intent(in),   optional :: dtype

        if (present(forcing_var)) then
            if (present(dtype)) then
                call var%initialize(grid, forcing_var=forcing_var, dtype=dtype)
            else
                call var%initialize(grid, forcing_var=forcing_var)
            endif
            var%forcing_var = forcing_var
            if (present(force_boundaries)) var%force_boundaries = force_boundaries
        else

            if (present(dtype)) then
                call var%initialize(grid, dtype=dtype)
            else
                call var%initialize(grid)
            endif
        endif

    end subroutine

    !> ---------------------------------
    !! Read the core model variables from disk
    !!
    !! Reads Terrain, lat, lon and u/v lat/lon on the high-res domain grid
    !! Passing data between images and disk is handled by io_read
    !!
    !! ---------------------------------
    subroutine read_core_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options
        real, allocatable :: temporary_data(:,:), temp_offset(:,:)

        ! Read the terrain data
        call io_read(options%domain%init_conditions_file,   &
                       options%domain%hgt_hi,                 &
                       temporary_data)
        this%vars_2d(this%var_indx(kVARS%terrain)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        
        this%vars_2d(this%var_indx(kVARS%neighbor_terrain)%v)%data_2d = temporary_data(this%ihs:this%ihe, this%jhs:this%jhe)

        allocate(temp_offset(1:this%grid%ide+1,1:this%grid%jde+1))

        if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
             (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                this%vars_2d(this%var_indx(kVARS%global_terrain)%v)%data_2d = temporary_data
        end if

        !while we have global terrain loaded, pass to split_topography
        if (options%domain%sleve) call split_topography(this, temporary_data, options)  ! here h1 and h2 are calculated

        ! Read the latitude data
        call io_read(options%domain%init_conditions_file,   &
                       options%domain%lat_hi,                 &
                       temporary_data)

        call make_2d_y(temporary_data, this%grid%ims, this%grid%ime)
        this%vars_2d(this%var_indx(kVARS%latitude)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        ! allocate(this%grid_vars(this%var_indx(kVARS%latitude_global)%v), source=temporary_data)

        ! Read the longitude data
        call io_read(options%domain%init_conditions_file,   &
                       options%domain%lon_hi,                 &
                       temporary_data)
        call make_2d_x(temporary_data, this%grid%jms, this%grid%jme)
        this%vars_2d(this%var_indx(kVARS%longitude)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
        ! allocate(this%grid_vars(this%var_indx(kVARS%longitude_global)%v), source=temporary_data)

        !-----------------------------------------
        !
        ! Handle staggered lat/lon grids, straightfoward if ulat/ulon are supplied
        ! If not, then read in mass grid lat/lon and stagger them
        !
        !-----------------------------------------
        ! Read the u-grid longitude data if specified, other wise interpolate from mass grid
        if (options%domain%ulon_hi /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%ulon_hi,                &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde)
            this%vars_2d(this%var_indx(kVARS%u_longitude)%v)%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lon_hi,                 &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde)
            call array_offset_x(temporary_data, temp_offset)
            this%vars_2d(this%var_indx(kVARS%u_longitude)%v)%data_2d = temp_offset(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        endif

        ! Read the u-grid latitude data if specified, other wise interpolate from mass grid
        if (options%domain%ulat_hi /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%ulat_hi,                &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide+1)
            this%vars_2d(this%var_indx(kVARS%u_latitude)%v)%data_2d = temporary_data(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lat_hi,                 &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide+1)
            call array_offset_x(temporary_data, temp_offset)
            this%vars_2d(this%var_indx(kVARS%u_latitude)%v)%data_2d = temp_offset(this%u_grid%ims:this%u_grid%ime,this%u_grid%jms:this%u_grid%jme)
        endif

        ! Read the v-grid longitude data if specified, other wise interpolate from mass grid
        if (options%domain%vlon_hi /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%vlon_hi,                &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde+1)
            this%vars_2d(this%var_indx(kVARS%v_longitude)%v)%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lon_hi,                 &
                           temporary_data)

            call make_2d_y(temporary_data, 1, this%jde+1)
            call array_offset_y(temporary_data, temp_offset)
            this%vars_2d(this%var_indx(kVARS%v_longitude)%v)%data_2d = temp_offset(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        endif

        ! Read the v-grid latitude data if specified, other wise interpolate from mass grid
        if (options%domain%vlat_hi /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%vlat_hi,                &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide)
            this%vars_2d(this%var_indx(kVARS%v_latitude)%v)%data_2d = temporary_data(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        else
            ! load the mass grid data again to get the full grid
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lat_hi,                 &
                           temporary_data)

            call make_2d_x(temporary_data, 1, this%ide)
            call array_offset_y(temporary_data, temp_offset)
            this%vars_2d(this%var_indx(kVARS%v_latitude)%v)%data_2d = temp_offset(this%v_grid%ims:this%v_grid%ime,this%v_grid%jms:this%v_grid%jme)
        endif

        if (STD_OUT_PE) write(*,*) "  Finished reading core domain variables"
        if (STD_OUT_PE) flush(output_unit)

    end subroutine


    !> -------------------------------
    !! Setup a single Geographic structure given a latitude, longitude, and z array
    !!
    !! -------------------------------
    subroutine setup_geo(geo, latitude, longitude, longitude_system, z)
        implicit none
        type(interpolable_type),  intent(inout) :: geo
        real,                     intent(in)    :: latitude(:,:)
        real,                     intent(in)    :: longitude(:,:)
        integer,                  intent(in)    :: longitude_system
        real, optional,           intent(in)    :: z(:,:,:)
        if (allocated(geo%lat)) deallocate(geo%lat)
        allocate( geo%lat, source=latitude)

        if (allocated(geo%lon)) deallocate(geo%lon)
        allocate( geo%lon, source=longitude)

        if (present(z)) then
            if (allocated(geo%z)) deallocate(geo%z)
            allocate( geo%z, source=z)
        endif
        ! This makes 2D variables out of lat/lon if they come in as 1D variables
        ! This also puts the longitudes onto a 0-360 if they are -180-180 (important for Alaska)
        ! Though if working in Europe the -180-180 grid is better ideally the optimal value should be checked.
        ! and good luck if you want to work over the poles...
        call standardize_coordinates(geo, longitude_system)

    end subroutine


    function find_flat_model_level(options, nz, dz) result(max_level)
        implicit none
        type(options_t), intent(in) :: options
        integer,         intent(in) :: nz
        real,            intent(in) :: dz(:)
        integer :: max_level

        integer :: j
        real :: height

        if (options%domain%flat_z_height > nz) then
            if (STD_OUT_PE) write(*,*) "    Treating flat_z_height as specified in meters above mean terrain height: ", options%domain%flat_z_height," meters"
            height = 0
            do j = 1, nz
                if (height <= options%domain%flat_z_height) then
                    height = height + dz(j)
                    max_level = j
                endif
            enddo

        elseif (options%domain%flat_z_height <= 0) then
            if (STD_OUT_PE) write(*,*) "    Treating flat_z_height as counting levels down from the model top: ", options%domain%flat_z_height," levels"
            max_level = nz + options%domain%flat_z_height

        else
            if (STD_OUT_PE) write(*,*) "    Treating flat_z_height as counting levels up from the ground: ", options%domain%flat_z_height," levels"
            max_level = options%domain%flat_z_height
        endif
        if (STD_OUT_PE) flush(output_unit)

    end function find_flat_model_level


    !> -------------------------------
    !! Setup the SLEVE vertical grid structure.
    !!   This basically entails 2 transformations: First a linear one so that sum(dz) ranges from 0 to smooth_height H.
    !!   (boundary cnd (3) in Schär et al 2002)  Next, the nonlinear SLEVE transformation
    !!    eqn (2) from Leuenberger et al 2009 z_sleve = Z + terrain * sinh((H/s)**n - (Z/s)**n) / SINH((H/s)**n) (for both smallscale and largescale terrain)
    !!   Here H is the model top or (flat_z_height in m), s controls how fast the terrain decays
    !!   and n controls the compression throughout the column (this last factor was added by Leuenberger et al 2009)
    !!   References: Leuenberger et al 2009 "A Generalization of the SLEVE Vertical Coordinate"
    !!               Schär et al 2002 "A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models"
    !!
    !! N.B. flat dz height != 0 makes little sense here? But works (?)
    !! -------------------------------
    subroutine setup_sleve(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:), gamma_n(:), neighbor_jacobian(:,:,:), neighbor_z(:,:,:)
        integer :: i, max_level
        real :: s, n, s1, s2, gamma, gamma_min
        real :: b1_i, b1_mass, db1_i, db1_mass, b2_i, b2_mass, db2_i, db2_mass
        
        associate(ims => this%ims,      ime => this%ime,                        &
            jms => this%jms,      jme => this%jme,                        &
            kms => this%kms,      kme => this%kme,                        &
            z                     => this%vars_3d(this%var_indx(kVARS%z)%v)%data_3d,                      &
            z_u                   => this%geo_u%z,                        &
            z_v                   => this%geo_v%z,                        &
            z_interface           => this%vars_3d(this%var_indx(kVARS%z_interface)%v)%data_3d,            &
            nz                    => options%domain%nz,               &
            dz                    => options%domain%dz_levels,        &
            dz_mass               => this%vars_3d(this%var_indx(kVARS%dz)%v)%data_3d,                &
            dz_interface          => this%vars_3d(this%var_indx(kVARS%dz_interface)%v)%data_3d,           &
            terrain               => this%vars_2d(this%var_indx(kVARS%terrain)%v)%data_2d,                &
            h1                    => this%vars_2d(this%var_indx(kVARS%h1)%v)%data_2d,                &
            h2                    => this%vars_2d(this%var_indx(kVARS%h2)%v)%data_2d,                &
            h1_u                  => this%vars_2d(this%var_indx(kVARS%h1_u)%v)%data_2d,                &
            h2_u                  => this%vars_2d(this%var_indx(kVARS%h2_u)%v)%data_2d,                &
            h1_v                  => this%vars_2d(this%var_indx(kVARS%h1_v)%v)%data_2d,                &
            h2_v                  => this%vars_2d(this%var_indx(kVARS%h2_v)%v)%data_2d,                &
            global_z_interface    => this%vars_3d(this%var_indx(kVARS%global_z_interface)%v)%data_3d,             &
            global_dz_interface   => this%vars_3d(this%var_indx(kVARS%global_dz_interface)%v)%data_3d,            &
            jacobian_u            => this%vars_3d(this%var_indx(kVARS%jacobian_u)%v)%data_3d,                     &
            jacobian_v            => this%vars_3d(this%var_indx(kVARS%jacobian_v)%v)%data_3d,                     &
            jacobian_w            => this%vars_3d(this%var_indx(kVARS%jacobian_w)%v)%data_3d,                     &
            dzdx                  => this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d,                   &
            dzdy                  => this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d,                   &
            dzdx_u                => this%vars_3d(this%var_indx(kVARS%dzdx_u)%v)%data_3d,                         &
            dzdy_v                => this%vars_3d(this%var_indx(kVARS%dzdy_v)%v)%data_3d,                         &
            jacobian              => this%vars_3d(this%var_indx(kVARS%jacobian)%v)%data_3d,                       &
            smooth_height         => this%smooth_height)

            ! Still not 100% convinced this works well in cases other than flat_z_height = 0 (w sleve). So for now best to keep at 0 when using sleve?
            max_level = nz !find_flat_model_level(options, nz, dz)

            smooth_height = sum(dz(1:max_level))!+dz(max_level)*0.5

            ! Terminology from Schär et al 2002, Leuenberger 2009: (can be simpliied later on, but for clarity)
            s1 = smooth_height / options%domain%decay_rate_L_topo
            s2 = smooth_height / options%domain%decay_rate_S_topo
            n  =  options%domain%sleve_n 

            ! Scale dz with smooth_height/sum(dz(1:max_level)) before calculating sleve levels.
            ! dz_scl(:)   =   dz(1:nz) !*  smooth_height / sum(dz(1:max_level))  ! this leads to a jump in dz thickness at max_level+1. Not sure if this is a problem.

            ! - - -   calculate invertibility parameter gamma (Schär et al 2002 eqn 20):  - - - - - -
            gamma  =  1  -  MAXVAL(h1)/s1 * COSH(smooth_height/s1)/SINH(smooth_height/s1) &
                          - MAXVAL(h2)/s2 * COSH(smooth_height/s2)/SINH(smooth_height/s2)

            ! with the new (leuenberger et al 2010) Sleve formulation, the inveribiltiy criterion is as follows:
            ! ( Although an argument could be made to calculate this on the offset (u/v) grid b/c that is most
            !   relevant for advection? In reality this is probably a sufficient approximation, as long as we
            !   aren't pushing the gamma factor too close to zero )
            allocate(gamma_n(this%kds : this%kde+1))
            i=kms
            gamma_n(i) =  1                                                     &
                - MAXVAL(h1) * n/(s1**n)                                        &
                * COSH((smooth_height/s1)**n) / SINH((smooth_height/s1)**n)     &
                - MAXVAL(h2) * n/(s2**n)                                        &
                * COSH((smooth_height/s2)**n) / SINH((smooth_height/s2)**n)

            do i = this%grid%kds, this%grid%kde
                gamma_n(i+1)  =  1                                    &    ! # for i != kds !!
                - MAXVAL(h1) * n/(s1**n) * sum(dz(1:i))**(n-1)                                             &
                * COSH((smooth_height/s1)**n -(sum(dz(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)    &
                - MAXVAL(h2) * n/(s2**n) *  sum(dz(1:i))**(n-1)                                            &
                * COSH((smooth_height/s2)**n -(sum(dz(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
            enddo

            if (n==1) then
                gamma_min = gamma
            else
                gamma_min = MINVAL(gamma_n)
            endif


            ! For reference: COSMO1 operational setting (but model top is at ~22000 masl):
            !    Decay Rate for Large-Scale Topography: svc1 = 10000.0000
            !    Decay Rate for Small-Scale Topography: svc2 =  3300.0000
            if ((STD_OUT_PE)) then
                write(*,*) "    Using a SLEVE coordinate with a Decay height for Large-Scale Topography: (s1) of ", s1, " m."
                write(*,*) "    Using a SLEVE coordinate with a Decay height for Small-Scale Topography: (s2) of ", s2, " m."
                write(*,*) "    Using a sleve_n of ", options%domain%sleve_n
                write(*,*) "    Smooth height is ", smooth_height, "m.a.s.l     (model top ", sum(dz(1:nz)), "m.a.s.l.)"
                write(*,*) "    invertibility parameter gamma is: ", gamma_min
                if(gamma_min <= 0) write(*,*) " CAUTION: coordinate transformation is not invertible (gamma <= 0 ) !!! reduce decay rate(s), and/or increase flat_z_height!"
                ! if(options%general%debug)  write(*,*) "   (for (debugging) reference: 'gamma(n=1)'= ", gamma,")"
                write(*,*) ""
                flush(output_unit)
            endif

            ! use temp to store global z-interface so that global-jacobian can be calculated

            if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
                 (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                 
                allocate(temp(this%ids:this%ide, this%kds:this%kde+1, this%jds:this%jde))
                temp(:,kms,:)   = this%vars_2d(this%var_indx(kVARS%global_terrain)%v)%data_2d
            else
                allocate(temp(this%ihs:this%ihe, this%khs:this%khe+1, this%jhs:this%jhe))
                temp(:,kms,:)   = this%vars_2d(this%var_indx(kVARS%neighbor_terrain)%v)%data_2d
            endif
            
            allocate(neighbor_jacobian(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
            allocate(neighbor_z(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
                                    
            ! - - - - -  k levels  - - - - -
            do i = this%grid%kms, this%grid%kme

                if (i==kms) then
                    b1_i = SINH( (smooth_height/s1)**n - (dz(i)/s1)**n ) / SINH((smooth_height/s1)**n)
                    b2_i = SINH( (smooth_height/s2)**n - (dz(i)/s2)**n ) / SINH((smooth_height/s2)**n)
                    b1_mass = SINH( (smooth_height/s1)**n -  ( (dz(i)/2) /s1)**n ) / SINH((smooth_height/s1)**n)
                    b2_mass = SINH( (smooth_height/s2)**n -  ( (dz(i)/2) /s2)**n ) / SINH((smooth_height/s2)**n)

                    db1_i = -n/(s1**n) * dz(i)**(n-1) * COSH((smooth_height/s1)**n - & 
                            (dz(i)/s1)**n ) / SINH((smooth_height/s1)**n)
                    db2_i = -n/(s2**n) * dz(i)**(n-1) * COSH((smooth_height/s2)**n - & 
                            (dz(i)/s2)**n ) / SINH((smooth_height/s2)**n)

                    db1_mass = -n/(s1**n) * (dz(i)/2)**(n-1) * COSH((smooth_height/s1)**n - &
                            ((dz(i)/2)/s1)**n ) / SINH((smooth_height/s1)**n)
                    db2_mass = -n/(s2**n) * (dz(i)/2)**(n-1) * COSH((smooth_height/s2)**n - &
                            ((dz(i)/2)/s2)**n ) / SINH((smooth_height/s2)**n)

                    temp(:,i+1,:)  = dz(i) + h1*b1_i + h2*b2_i

                    global_dz_interface(:,i,:)  =  temp(:,i+1,:) - temp(:,i,:)  ! same for higher k
                    global_z_interface(:,i,:)  = temp(:,i,:)

                    dz_mass(:,i,:)       = global_dz_interface(ims:ime,i,jms:jme) / 2           ! Diff for k=1            

                    ! ! - - - - -   u/v grid calculations for lowest level (i=kms)  - - - - -
                    ! ! for the u and v grids, z(1) was already initialized with terrain.
                    ! ! but the first level needs to be offset, and the rest of the levels need to be created
                    ! ! BK: So if z_u is already offset in the u dir, but not in the z dir, we can say that
                    ! !     z_u(:,1,:) is the terrain on the u grid, and it needs to be offset in the z-dir
                    ! !     to reach mass levels (so by dz[i]/2)

                    neighbor_z(:,i,:)  = (dz(i)/2)  + h1(this%ihs:this%ihe,this%jhs:this%jhe)*b1_mass + &
                                                          h2(this%ihs:this%ihe,this%jhs:this%jhe)*b2_mass
                    z_u(:,i,:)   = (dz(i)/2) + h1_u*b1_mass + h2_u*b2_mass
                    z_v(:,i,:)   = (dz(i)/2) + h1_v*b1_mass + h2_v*b2_mass

                else if(i>kms) then
                    if(i<=max_level) then

                        b1_i = SINH( (smooth_height/s1)**n - (sum(dz(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)
                        b2_i = SINH( (smooth_height/s2)**n - (sum(dz(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
                        b1_mass = SINH( (smooth_height/s1)**n -  ( (sum(dz(1:(i-1)))+dz(i)/2) /s1)**n ) / SINH((smooth_height/s1)**n)
                        b2_mass = SINH( (smooth_height/s2)**n -  ( (sum(dz(1:(i-1)))+dz(i)/2) /s2)**n ) / SINH((smooth_height/s2)**n)


                        db1_i = -n/(s1**n) * sum(dz(1:i))**(n-1) * COSH((smooth_height/s1)**n - & 
                                (sum(dz(1:i))/s1)**n ) / SINH((smooth_height/s1)**n)
                        db2_i = -n/(s2**n) * sum(dz(1:i))**(n-1) * COSH((smooth_height/s2)**n - & 
                                (sum(dz(1:i))/s2)**n ) / SINH((smooth_height/s2)**n)
                                
                        db1_mass = -n/(s1**n) * (sum(dz(1:(i-1)))+dz(i)/2)**(n-1) * COSH((smooth_height/s1)**n - &
                                ((sum(dz(1:(i-1)))+dz(i)/2)/s1)**n ) / SINH((smooth_height/s1)**n)
                        db2_mass = -n/(s2**n) * (sum(dz(1:(i-1)))+dz(i)/2)**(n-1) * COSH((smooth_height/s2)**n - &
                                ((sum(dz(1:(i-1)))+dz(i)/2)/s2)**n ) / SINH((smooth_height/s2)**n)

                        temp(:,i+1,:)  = sum(dz(1:i)) + h1*b1_i + h2*b2_i 
                        global_dz_interface(:,i,:)  =  temp(:,i+1,:) - temp(:,i,:)

                        global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)

                        neighbor_z(:,i,:)  = (sum(dz(1:(i-1))) + dz(i)/2)  + h1(this%ihs:this%ihe,this%jhs:this%jhe)*b1_mass + &
                                                                                     h2(this%ihs:this%ihe,this%jhs:this%jhe)*b2_mass  
                        z_u(:,i,:)   = (sum(dz(1:(i-1))) + dz(i)/2) + h1_u*b1_mass + h2_u*b2_mass  
                        z_v(:,i,:)   = (sum(dz(1:(i-1))) + dz(i)/2) + h1_v*b1_mass + h2_v*b2_mass  

                        if ( ANY(global_z_interface(:,i,:)<0) ) then   ! Eror catching. Probably good to engage.
                            if (STD_OUT_PE) then
                                write(*,*) "Error: global_z_interface below zero (for level  ",i,")"
                                write(*,*)  "min max global_z_interface: ",MINVAL(global_z_interface(:,i,:)),MAXVAL(global_z_interface(:,i,:))
                                error stop
                            endif
                        else if ( ANY(global_z_interface(:,i,:)<=0.01) ) then
                            write(*,*) "WARNING: global_z_interface very low (at level ",i,")"
                        endif

                    else ! above the flat_z_height
                        b1_i = 0
                        b2_i = 0
                        b1_mass = 0
                        b2_mass = 0

                        db1_i = 0
                        db2_i = 0
                        db1_mass = 0
                        db2_mass = 0

                        global_dz_interface(:,i,:) =  dz(i)
                        global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
                        !  if (i/=this%grid%kme)   z_interface(:,i+1,:) = z_interface(:,i,:) + dz(i) ! (dz(i) + dz( i) )/2 !test in icar_s5T
                        z_u(:,i,:)  = z_u(:,i-1,:)  + (dz(i) + dz(i-1))*0.5 ! zr_u only relevant for first i above max level, aferwards both zr_u(i) AND zr_u(i-1) ar
                        z_v(:,i,:)  = z_v(:,i-1,:)  + (dz(i) + dz(i-1))*0.5
                        neighbor_z(:,i,:) =  neighbor_z(:,i-1,:) + (dz(i-1) + dz(i))*0.5

                    endif
                    dz_mass(:,i,:)   =  global_dz_interface(ims:ime,i-1,jms:jme) / 2  +  global_dz_interface(ims:ime,i,jms:jme) / 2
                endif ! if (i>kms)
                
                neighbor_jacobian(:,i,:) = 1 + h1(this%ihs:this%ihe,this%jhs:this%jhe)*db1_mass + h2(this%ihs:this%ihe,this%jhs:this%jhe)*db2_mass
                jacobian_u(:,i,:) = 1 + h1_u*db1_mass + h2_u*db2_mass
                jacobian_v(:,i,:) = 1 + h1_v*db1_mass + h2_v*db2_mass

                jacobian_w(:,i,:) = 1 + h1(ims:ime,jms:jme)*db1_i + h2(ims:ime,jms:jme)*db2_i
                
                dzdx(ims+1:ime-1,i,:) = (b1_mass*(h1(ims+2:ime,jms:jme)-h1(ims:ime-2,jms:jme)) + b2_mass*(h2(ims+2:ime,jms:jme)-h2(ims:ime-2,jms:jme)))/(2*this%dx)
                dzdx(ims,i,:)   = (-neighbor_z(ims+2,i,jms:jme) + 4*neighbor_z(ims+1,i,jms:jme) - 3*neighbor_z(ims,i,jms:jme) )/(2*this%dx)
                dzdx(ime,i,:)   = (neighbor_z(ime-2,i,jms:jme) - 4*neighbor_z(ime-1,i,jms:jme) + 3*neighbor_z(ime,i,jms:jme) )/(2*this%dx)

                dzdy(:,i,jms+1:jme-1) = (b1_mass*(h1(ims:ime,jms+2:jme)-h1(ims:ime,jms:jme-2)) + b2_mass*(h2(ims:ime,jms+2:jme)-h2(ims:ime,jms:jme-2)))/(2*this%dx)
                dzdy(:,i,jms)   = (-neighbor_z(ims:ime,i,jms+2) + 4*neighbor_z(ims:ime,i,jms+1) - 3*neighbor_z(ims:ime,i,jms) )/(2*this%dx)
                dzdy(:,i,jme)   = (neighbor_z(ims:ime,i,jme-2) - 4*neighbor_z(ims:ime,i,jme-1) + 3*neighbor_z(ims:ime,i,jme) )/(2*this%dx)
                
                dzdx_u(ims+1:ime,i,:) = (b1_mass*(h1(ims+1:ime,jms:jme)-h1(ims:ime-1,jms:jme))   + b2_mass*(h2(ims+1:ime,jms:jme)-h2(ims:ime-1,jms:jme)))/(this%dx)
                dzdy_v(:,i,jms+1:jme) = (b1_mass*(h1(ims:ime,jms+1:jme)-h1(ims:ime,jms:jme-1))   + b2_mass*(h2(ims:ime,jms+1:jme)-h2(ims:ime,jms:jme-1)))/(this%dx)

                dzdx_u(ims+1:ime,i,:) = (dzdx(ims:ime-1,i,:)+dzdx(ims+1:ime,i,:))*0.5
                dzdx_u(ims,i,:)   = dzdx(ims,i,:)*1.5 - dzdx(ims+1,i,:)*0.5
                dzdx_u(ime+1,i,:)   = dzdx(ime,i,:)*1.5 - dzdx(ime-1,i,:)*0.5

                dzdy_v(:,i,jms+1:jme) = (dzdy(:,i,jms:jme-1)+dzdy(:,i,jms+1:jme))*0.5
                dzdy_v(:,i,jms)   = dzdy(:,i,jms)*1.5 - dzdy(:,i,jms+1)*0.5
                dzdy_v(:,i,jme+1)   = dzdy(:,i,jme)*1.5 - dzdy(:,i,jme-1)*0.5
                
            enddo
            
            !Finishing touch
            i=kme+1
            global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
            
            jacobian_w(:,this%kme,:) = 1.0


            ! this is on the subset grid:
            dz_interface = global_dz_interface(ims:ime,:,jms:jme)
            z_interface  = global_z_interface(ims:ime,:,jms:jme)
            z            = neighbor_z(ims:ime,:,jms:jme)
            jacobian     = neighbor_jacobian(ims:ime,:,jms:jme)

        end associate
        ! call array_offset_x(neighbor_jacobian(this%ims:this%ime,:,this%jms:this%jme), temp)
        ! this%vars_3d(this%var_indx(kVARS%jacobian_u)%v)%data_3d(this%ims:this%ime+1,:,this%jms:this%jme) = temp
        ! call array_offset_y(neighbor_jacobian(this%ims:this%ime,:,this%jms:this%jme), temp)
        ! this%vars_3d(this%var_indx(kVARS%jacobian_v)%v)%data_3d(this%ims:this%ime,:,this%jms:this%jme+1) = temp

    end subroutine setup_sleve



    !> -------------------------------
    !! Setup the vertical grid structure, in case SLEVE coordinates are not used.
    !!    This means either constant vertical height, or a simple terrain following coordinate (Gal-Chen)
    !!
    !! --------------------------------
    subroutine setup_simple_z(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:), temp_offset(:,:), global_jacobian(:,:,:)
        integer :: i, max_level

        associate(  ims => this%ims,      ime => this%ime,                        &
                    jms => this%jms,      jme => this%jme,                        &
                    kms => this%kms,      kme => this%kme,                        &
                    z                     => this%vars_3d(this%var_indx(kVARS%z)%v)%data_3d,                      &
                    z_u                   => this%geo_u%z,                        &
                    z_v                   => this%geo_v%z,                        &
                    z_interface           => this%vars_3d(this%var_indx(kVARS%z_interface)%v)%data_3d,            &
                    nz                    => options%domain%nz,               &
                    dz                    => options%domain%dz_levels,        &
                    dz_mass               => this%vars_3d(this%var_indx(kVARS%dz)%v)%data_3d,                &
                    dz_interface          => this%vars_3d(this%var_indx(kVARS%dz_interface)%v)%data_3d,           &
                    terrain               => this%vars_2d(this%var_indx(kVARS%terrain)%v)%data_2d,                &
                    global_z_interface    => this%vars_3d(this%var_indx(kVARS%global_z_interface)%v)%data_3d,             &
                    global_dz_interface   => this%vars_3d(this%var_indx(kVARS%global_dz_interface)%v)%data_3d,            &
                    neighbor_terrain      => this%vars_2d(this%var_indx(kVARS%neighbor_terrain)%v)%data_2d,               &
                    global_terrain        => this%vars_2d(this%var_indx(kVARS%global_terrain)%v)%data_2d,               &
                    jacobian_u            => this%vars_3d(this%var_indx(kVARS%jacobian_u)%v)%data_3d,                     &
                    jacobian_v            => this%vars_3d(this%var_indx(kVARS%jacobian_v)%v)%data_3d,                     &
                    jacobian_w            => this%vars_3d(this%var_indx(kVARS%jacobian_w)%v)%data_3d,                     &
                    dzdy                  => this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d,                           &
                    jacobian              => this%vars_3d(this%var_indx(kVARS%jacobian)%v)%data_3d,                       &
                    smooth_height         => this%smooth_height)

            ! Start with a separate calculation for the lowest model level z=1
            i = this%grid%kms
            
            
            if ( (options%physics%windtype == kWIND_LINEAR) .or. (options%physics%windtype == kLINEAR_OBRIEN_WINDS) .or. &
                 (options%physics%windtype == kLINEAR_ITERATIVE_WINDS) ) then
                global_z_interface(:,i,:)   = global_terrain
                allocate(global_jacobian(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde))
            else
                global_z_interface(:,i,:)   = neighbor_terrain
                allocate(global_jacobian(this%ihs:this%ihe, this%khs:this%khe, this%jhs:this%jhe))
            endif


            max_level = nz !find_flat_model_level(options, nz, dz)

            smooth_height = sum(dz(1:max_level))

            jacobian(:,i,:) = (smooth_height - terrain) / smooth_height ! sum(dz(1:max_level))
            global_jacobian(:,i,:) = (smooth_height - global_z_interface(:,i,:) ) /smooth_height !sum(dz(1:max_level))

            dz_mass(:,i,:)      = dz(i) / 2 * jacobian(:,i,:)
            dz_interface(:,i,:) = dz(i) * jacobian(:,i,:)
            z(:,i,:)            = terrain + dz_mass(:,i,:)
            z_interface(:,i,:)  = terrain

            global_dz_interface(:,i,:) = dz(i) * global_jacobian(:,i,:)
            
            ! Now the higher (k!=1) levels can be calculated:
            do i = this%grid%kms+1, this%grid%kme
                if (i<=max_level) then
                    jacobian(:,i,:) = jacobian(:,i-1,:)
                    global_jacobian(:,i,:) = global_jacobian(:,i-1,:)
                else
                    jacobian(:,i,:) = 1
                    global_jacobian(:,i,:) = 1
                endif

                dz_mass(:,i,:)     = (dz(i)/2 * jacobian(:,i,:) + dz(i-1)/2 * jacobian(:,i-1,:))
                dz_interface(:,i,:)= dz(i) * jacobian(:,i,:)
                z(:,i,:)           = z(:,i-1,:)           + dz_mass(:,i,:)
                z_interface(:,i,:) = z_interface(:,i-1,:) + dz_interface(:,i-1,:)

                global_dz_interface(:,i,:) = dz(i) * global_jacobian(:,i,:)
                global_z_interface(:,i,:)  = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)


            enddo

            i = this%grid%kme + 1
            global_z_interface(:,i,:) = global_z_interface(:,i-1,:) + global_dz_interface(:,i-1,:)
            
            if (allocated(temp)) deallocate(temp)
            allocate(temp(this%ihs:this%ihe+1, this%khs:this%khe, this%jhs:this%jhe+1))

            temp(this%ihs,:,this%jhs:this%jhe) = global_jacobian(this%ihs,:,this%jhs:this%jhe)
            temp(this%ihe+1,:,this%jhs:this%jhe) = global_jacobian(this%ihe,:,this%jhs:this%jhe)
            temp(this%ihs+1:this%ihe,:,this%jhs:this%jhe) = (global_jacobian(this%ihs+1:this%ihe,:,this%jhs:this%jhe) + &
                                                                global_jacobian(this%ihs:this%ihe-1,:,this%jhs:this%jhe))/2
            jacobian_u = temp(ims:ime+1,:,jms:jme)

            temp(this%ihs:this%ihe,:,this%jhs) = global_jacobian(this%ihs:this%ihe,:,this%jhs)
            temp(this%ihs:this%ihe,:,this%jhe+1) = global_jacobian(this%ihs:this%ihe,:,this%jhe)
            temp(this%ihs:this%ihe,:,this%jhs+1:this%jhe) = (global_jacobian(this%ihs:this%ihe,:,this%jhs+1:this%jhe) + &
                                                global_jacobian(this%ihs:this%ihe,:,this%jhs:this%jhe-1))/2
            jacobian_v = temp(ims:ime,:,jms:jme+1)

            jacobian_w(:,this%kme,:) = 1.0 !jacobian(:,this%kme,:)
            jacobian_w(:,this%kms:this%kme-1,:) = (dz_interface(:,this%kms+1:this%kme,:)* jacobian(:,this%kms:this%kme-1,:) + &
                                                   dz_interface(:,this%kms:this%kme-1,:)* jacobian(:,this%kms+1:this%kme,:))/ &
                                                                                (dz_interface(:,this%kms:this%kme-1,:)+dz_interface(:,this%kms+1:this%kme,:))
                                                                                

            call array_offset_x(terrain, temp_offset)
            z_u(:,1,:) = temp_offset
            call array_offset_y(terrain, temp_offset)
            z_v(:,1,:) = temp_offset

            z_u(:,1,:)          = z_u(:,1,:) + dz(1) / 2 * jacobian_u(:,1,:)
            z_v(:,1,:)          = z_v(:,1,:) + dz(1) / 2 * jacobian_v(:,1,:)

            do i = this%grid%kms+1, this%grid%kme
                z_u(:,i,:) = z_u(:,i-1,:)  + ((dz(i)/2 * jacobian_u(:,i,:) + dz(i-1)/2 * jacobian_u(:,i-1,:)))
                z_v(:,i,:) = z_v(:,i-1,:)  + ((dz(i)/2 * jacobian_v(:,i,:) + dz(i-1)/2 * jacobian_v(:,i-1,:)))  
            enddo
                                                                                
            call setup_dzdxy(this, options, global_jacobian)

        end associate

    end subroutine setup_simple_z



    !> -------------------------------
    !! Initialize various domain variables, mostly z, dz, etc.
    !!
    !! -------------------------------
    subroutine initialize_core_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:,:)
        integer :: i, j

        call read_core_variables(this, options)

        !setup geo_u/v here, because their z arrays will be calculated in the setup methods below
        call setup_geo(this%geo_u,   this%vars_2d(this%var_indx(kVARS%u_latitude)%v)%data_2d,   this%vars_2d(this%var_indx(kVARS%u_longitude)%v)%data_2d, options%domain%longitude_system)
        call setup_geo(this%geo_v,   this%vars_2d(this%var_indx(kVARS%v_latitude)%v)%data_2d,   this%vars_2d(this%var_indx(kVARS%v_longitude)%v)%data_2d, options%domain%longitude_system)
        allocate( this%geo_u%z(this%u_grid%ims:this%u_grid%ime, this%u_grid%nz, this%u_grid%jms:this%u_grid%jme))
        allocate( this%geo_v%z(this%v_grid%ims:this%v_grid%ime, this%v_grid%nz, this%v_grid%jms:this%v_grid%jme))

        do i=this%grid%kms, this%grid%kme
            this%vars_3d(this%var_indx(kVARS%advection_dz)%v)%data_3d(:,i,:) = options%domain%dz_levels(i)
        enddo

        ! Setup the vertical grid structure, either as a SLEVE coordinate, or a more 'simple' vertical structure:
        if (options%domain%sleve) then
            call setup_sleve(this, options)
        else
            ! This will set up either a Gal-Chen terrainfollowing coordinate, or no terrain following.
            call setup_simple_z(this, options)
        endif
        
        call setup_geo(this%geo,   this%vars_2d(this%var_indx(kVARS%latitude)%v)%data_2d,   this%vars_2d(this%var_indx(kVARS%longitude)%v)%data_2d, options%domain%longitude_system,  this%vars_3d(this%var_indx(kVARS%z)%v)%data_3d)
        call setup_geo(this%geo_agl,   this%vars_2d(this%var_indx(kVARS%latitude)%v)%data_2d,   this%vars_2d(this%var_indx(kVARS%longitude)%v)%data_2d, options%domain%longitude_system,  this%vars_3d(this%var_indx(kVARS%z)%v)%data_3d)

        call setup_grid_rotations(this, options)
        

    end subroutine initialize_core_variables
        
    subroutine setup_grid_rotations(this,options)
        type(domain_t),  intent(inout) :: this
        type(options_t), intent(in)    :: options

        integer :: i, j, i_s, i_e, j_s, j_e, smooth_loops
        integer :: starti, endi, smooth_window_size
        double precision :: dist, dlat, dlon

        real, allocatable :: lat(:,:), lon(:,:), costheta(:,:), sintheta(:,:)

        
        if (options%domain%sinalpha_var /= "") then
            i_s = this%ims
            i_e = this%ime
            j_s = this%jms
            j_e = this%jme

            if (STD_OUT_PE) print*, "Reading Sinalpha/cosalpha"
            if (STD_OUT_PE) flush(output_unit)

            call io_read(options%domain%init_conditions_file, options%domain%sinalpha_var, lon)
            this%vars_2d(this%var_indx(kVARS%sintheta)%v)%data_2d = lon(i_s:i_e, j_s:j_e)

            call io_read(options%domain%init_conditions_file, options%domain%cosalpha_var, lon)
            this%vars_2d(this%var_indx(kVARS%costheta)%v)%data_2d = lon(i_s:i_e, j_s:j_e)

        else

            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lat_hi,                 &
                           lat)
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lon_hi,                 &
                           lon)

            smooth_window_size = 50

            i_s = this%ids!max(this%ims-smooth_window_size,this%ids)
            i_e = this%ide!min(this%ime+smooth_window_size,this%ide)
            j_s = this%jds!max(this%jms-smooth_window_size,this%jms)
            j_e = this%jde!min(this%jme+smooth_window_size,this%jde)

            allocate(sintheta(i_s:i_e,j_s:j_e))
            allocate(costheta(i_s:i_e,j_s:j_e))

            do j = j_s, j_e
                do i = i_s, i_e
                    ! in case we are in the first or last grid, reset boundaries
                    starti = max(this%ids, i-2)
                    endi   = min(this%ide, i+2)

                    ! change in latitude
                    dlat = DBLE(lat(endi,j) - lat(starti,j))
                    ! change in longitude
                    dlon = DBLE(lon(endi,j) - lon(starti,j)) * cos(deg2rad*DBLE(lat(i,j)))
                    !if (abs(dlat) > 1) write(*,*) 'dlat:  ', dlat, '  ', ims, '  ', ime, '  ', jms, '  ', jme
                    !if (abs(dlon) > 1) write(*,*) 'dlon:  ', dlon, '  ', ims, '  ', ime, '  ', jms, '  ', jme
                    
                    ! distance between two points
                    dist = sqrt(DBLE(dlat)**2 + DBLE(dlon)**2) 

                    ! sin/cos of angles for use in rotating fields later
                    costheta(i, j) = abs(dlon / dist)
                    sintheta(i, j) =  (-1) * dlat / dist

                enddo
            enddo
            
            !Smooth cos/sin in case there are jumps from the lat/lon grids (more likely at low resolutions)
            smooth_loops = int(1000/this%dx)
            
            do i=1,smooth_loops
             call smooth_array_2d( costheta , windowsize  =  4)!int((ime-ims)/5))
             call smooth_array_2d( sintheta , windowsize  =  4)!int((ime-ims)/5))
            enddo
            this%vars_2d(this%var_indx(kVARS%costheta)%v)%data_2d(this%ims:this%ime,this%jms:this%jme) = costheta(this%ims:this%ime,this%jms:this%jme)
            this%vars_2d(this%var_indx(kVARS%sintheta)%v)%data_2d(this%ims:this%ime,this%jms:this%jme) = sintheta(this%ims:this%ime,this%jms:this%jme)
             
        endif
        if (options%general%debug .and.(STD_OUT_PE)) then
            print*, ""
            print*, "Domain Geometry"
            print*, "MAX / MIN SIN(theta) (ideally 0)"
            print*, "   ", maxval(this%vars_2d(this%var_indx(kVARS%sintheta)%v)%data_2d), minval(this%vars_2d(this%var_indx(kVARS%sintheta)%v)%data_2d)
            print*, "MAX / MIN COS(theta) (ideally 1)"
            print*, "   ", maxval(this%vars_2d(this%var_indx(kVARS%costheta)%v)%data_2d), minval(this%vars_2d(this%var_indx(kVARS%costheta)%v)%data_2d)
            print*, ""
            flush(output_unit)
        endif


    end subroutine setup_grid_rotations

    subroutine setup_dzdxy(this, options, neighbor_jacobian)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options
        real, allocatable, intent(in)   :: neighbor_jacobian(:,:,:)
        
        real, allocatable :: neighbor_z(:,:,:)
        real, allocatable :: neighbor_dzdx(:,:,:)
        real, allocatable :: neighbor_dzdy(:,:,:)
        integer :: i

        allocate(neighbor_z( this% ihs : this% ihe, this% khs : this% khe, this% jhs : this% jhe) )
        allocate(neighbor_dzdx( this% ihs : this% ihe+1, this% khs : this% khe, this% jhs : this% jhe) )
        allocate(neighbor_dzdy( this% ihs : this% ihe, this% khs : this% khe, this% jhs : this% jhe+1) )

        if (this%var_indx(kVARS%neighbor_terrain)%v > 0) then
                neighbor_z(:,1,:) = this%vars_2d(this%var_indx(kVARS%neighbor_terrain)%v)%data_2d + (options%domain%dz_levels(1)/2)*neighbor_jacobian(:,1,:)
        else
                neighbor_z(:,1,:) = this%vars_2d(this%var_indx(kVARS%global_terrain)%v)%data_2d(this%ihs:this%ihe, this%jhs:this%jhe) + (options%domain%dz_levels(1)/2)*neighbor_jacobian(:,1,:)
        endif
        do i=2,this%khe
            neighbor_z(:,i,:) = neighbor_z(:,i-1,:) + (((options%domain%dz_levels(i)) / 2)*neighbor_jacobian(:,i,:)) + &
                                                  (((options%domain%dz_levels(i-1)) / 2)*neighbor_jacobian(:,i-1,:))
        enddo

        neighbor_dzdx = 0
        neighbor_dzdy = 0

        !For dzdx
        neighbor_dzdx(this%ihs+1:this%ihe-1,:,:) = (neighbor_z(this%ihs+2:this%ihe,:,:) - &
                                                           neighbor_z(this%ihs:this%ihe-2,:,:))/(2*this%dx)
                                                                                                          
        neighbor_dzdx(this%ihs,:,:) = (-3*neighbor_z(this%ihs,:,:) + &
                                          4*neighbor_z(this%ihs+1,:,:) - neighbor_z(this%ihs+2,:,:)) / (2*this%dx)
                                          
        neighbor_dzdx(this%ihe,:,:) = (3*neighbor_z(this%ihe,:,:) - &
                                         4*neighbor_z(this%ihe-1,:,:) + neighbor_z(this%ihe-2,:,:)) / (2*this%dx)
        
        this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d(:,:,:) = neighbor_dzdx(this%ims:this%ime,:,this%jms:this%jme)
        
        

        neighbor_dzdx(this%ihs+1:this%ihe,:,:) = (neighbor_z(this%ihs+1:this%ihe,:,:) - neighbor_z(this%ihs:this%ihe-1,:,:))/this%dx
        neighbor_dzdx(this%ihs,:,:) = neighbor_dzdx(this%ihs+1,:,:) 
        neighbor_dzdx(this%ihe+1,:,:) = neighbor_dzdx(this%ihe,:,:)
        this%vars_3d(this%var_indx(kVARS%dzdx_u)%v)%data_3d(this%ims+1:this%ime,:,:) = neighbor_dzdx(this%ims+1:this%ime,:,this%jms:this%jme)
        
        this%vars_3d(this%var_indx(kVARS%dzdx_u)%v)%data_3d(this%ims,:,:)   = this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d(this%ims,:,:)*1.5 - this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d(this%ims+1,:,:)*0.5
        this%vars_3d(this%var_indx(kVARS%dzdx_u)%v)%data_3d(this%ime+1,:,:) = this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d(this%ime,:,:)*1.5 - this%vars_3d(this%var_indx(kVARS%dzdx)%v)%data_3d(this%ime-1,:,:)*0.5
        
        
        !For dzdy
        neighbor_dzdy(:,:,this%jhs+1:this%jhe-1) = (neighbor_z(:,:,this%jhs+2:this%jhe) - &
                                                           neighbor_z(:,:,this%jhs:this%jhe-2))/(2*this%dx)
        neighbor_dzdy(:,:,this%jhs) = (-3*neighbor_z(:,:,this%jms) + &
                                          4*neighbor_z(:,:,this%jms+1) - neighbor_z(:,:,this%jms+2)) / (2*this%dx)
                                          
        neighbor_dzdy(:,:,this%jhe) = (3*neighbor_z(:,:,this%jhe) - &
                                         4*neighbor_z(:,:,this%jhe-1) + neighbor_z(:,:,this%jhe-2)) / (2*this%dx)
        this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d(:,:,:) = neighbor_dzdy(this%ims:this%ime,:,this%jms:this%jme)


        neighbor_dzdy(:,:,this%jhs+1:this%jhe) = (neighbor_z(:,:,this%jhs+1:this%jhe) - neighbor_z(:,:,this%jhs:this%jhe-1))/this%dx
        neighbor_dzdy(:,:,this%jhs) = neighbor_dzdy(:,:,this%jhs+1) 
        neighbor_dzdy(:,:,this%jhe+1) = neighbor_dzdy(:,:,this%jhe)
                
        this%vars_3d(this%var_indx(kVARS%dzdy_v)%v)%data_3d(:,:,this%jms+1:this%jme) = neighbor_dzdy(this%ims:this%ime,:,this%jms+1:this%jme)
        this%vars_3d(this%var_indx(kVARS%dzdy_v)%v)%data_3d(:,:,this%jms)   = this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d(:,:,this%jms)*1.5 - this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d(:,:,this%jms+1)*0.5
        this%vars_3d(this%var_indx(kVARS%dzdy_v)%v)%data_3d(:,:,this%jme+1) = this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d(:,:,this%jme)*1.5 - this%vars_3d(this%var_indx(kVARS%dzdy)%v)%data_3d(:,:,this%jme-1)*0.5

    end subroutine setup_dzdxy


    !> -------------------------------
    !!  Separate the terrain into large scale and small scale terrain for SLEVE coordinate calculation
    !!  h(x,y) = h_1(x,y) + h_2(x,y) ;
    !!  where the subscripts 1 and 2 refer to large-scale and small-scale contributions, respectively.
    !!  The large-scale contribution h1 can be obtained from the full topography by an appropriate smoothing operation.
    !!
    !!  The smoothing is done over the entire (non-parallelized terrain, i.e. ids-ide). Afterwards the relevant variables
    !!  are subset to the respective paralellized grids. This is not the most efficient, but it makes the smoothing easier.
    !!
    !> -------------------------------

    subroutine split_topography(this, global_terr, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        real, dimension(this%ids:this%ide,this%jds:this%jde), intent(in)   :: global_terr
        type(options_t), intent(in)     :: options

        real, allocatable :: temp(:,:)
        integer :: i

                            
        associate(h1                    => this%vars_2d(this%var_indx(kVARS%h1)%v)%data_2d,                             &
                  h2                    => this%vars_2d(this%var_indx(kVARS%h2)%v)%data_2d,                             &
                  h1_u                  => this%vars_2d(this%var_indx(kVARS%h1_u)%v)%data_2d,                           &
                  h2_u                  => this%vars_2d(this%var_indx(kVARS%h2_u)%v)%data_2d,                           &
                  h1_v                  => this%vars_2d(this%var_indx(kVARS%h1_v)%v)%data_2d,                           &
                  h2_v                  => this%vars_2d(this%var_indx(kVARS%h2_v)%v)%data_2d,                           &
                  terrain               => this%vars_2d(this%var_indx(kVARS%terrain)%v)%data_2d)


        if ((STD_OUT_PE)) then
          write(*,*) "  Setting up the SLEVE vertical coordinate:"
          write(*,*) "    Smoothing large-scale terrain (h1) with a windowsize of ", &
                  options%domain%terrain_smooth_windowsize, " for ",        &
                  options%domain%terrain_smooth_cycles, " smoothing cylces."
        endif


        ! create a separate variable that will be smoothed later on:
        allocate(temp(this%ids:this%ide,this%jds:this%jde))
        temp =  global_terr

        ! Smooth the terrain to attain the large-scale contribution h1 (_u/v):
        do i =1,options%domain%terrain_smooth_cycles
          call smooth_array( temp, windowsize  =  options%domain%terrain_smooth_windowsize)
        enddo
        
        if (this%var_indx(kVARS%global_terrain)%v > 0) then
                h1   =  temp
                h2   =  global_terr - h1
        else
                h1   =  temp(this%ihs:this%ihe,this%jhs:this%jhe)
                h2   =  global_terr(this%ihs:this%ihe,this%jhs:this%jhe) - h1
        endif
        ! offset the global terrain for the h_(u/v) calculations:
        deallocate(temp)
        allocate(temp(this%ids:this%ide+1,this%jds:this%jde))
        call array_offset_x(global_terr, temp)
        !temp(this%ids,this%jds:this%jde) = temp(this%ids+1,this%jds:this%jde)
        !temp(this%ide+1,this%jds:this%jde) = temp(this%ide,this%jds:this%jde)
        
        h2_u = temp(this%u_grid2d%ims:this%u_grid2d%ime, this%u_grid2d%jms:this%u_grid2d%jme)
        do i =1,options%domain%terrain_smooth_cycles
          call smooth_array( temp, windowsize = options%domain%terrain_smooth_windowsize)
        enddo
        
        h1_u = temp(this%u_grid2d%ims:this%u_grid2d%ime, this%u_grid2d%jms:this%u_grid2d%jme)
        h2_u =  h2_u  - h1_u

        
        ! offset the global terrain for the h_(u/v) calculations:
        deallocate(temp)
        allocate(temp(this%ids:this%ide,this%jds:this%jde+1))
        call array_offset_y(global_terr, temp)
        !temp(this%ids:this%ide,this%jds) = temp(this%ids:this%ide,this%jds+1)
        !temp(this%ids:this%ide,this%jde+1) = temp(this%ids:this%ide,this%jde)
        h2_v = temp(this%v_grid2d%ims:this%v_grid2d%ime, this%v_grid2d%jms:this%v_grid2d%jme)
        
        do i =1,options%domain%terrain_smooth_cycles
          call smooth_array( temp, windowsize = options%domain%terrain_smooth_windowsize)
        enddo
        
        h1_v = temp(this%v_grid2d%ims:this%v_grid2d%ime, this%v_grid2d%jms:this%v_grid2d%jme)
        h2_v =  h2_v  - h1_v
        
        !if (STD_OUT_PE) then
        !   write(*,*) "       Max of full topography", MAXVAL(neighbor_terrain )
        !   write(*,*) "       Max of large-scale topography (h1)  ", MAXVAL(h1)
        !   write(*,*) "       Max of small-scale topography (h2)  ", MAXVAL(h2)
        !end if

        end associate

    end subroutine split_topography


    



    subroutine read_land_variables(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        integer :: i, nsoil
        real, allocatable :: temporary_data(:,:), temporary_data_3d(:,:,:)
        real :: soil_thickness(20)
        real :: init_surf_temp

        soil_thickness = 1.0
        soil_thickness(1:4) = [0.1, 0.2, 0.5, 1.0]
        init_surf_temp = 280

        if (STD_OUT_PE) write (*,*) "Reading Land Variables"
        if (STD_OUT_PE) flush(output_unit)

        if (this%var_indx(kVARS%soil_water_content)%v > 0) then
            nsoil = size(this%vars_3d(this%var_indx(kVARS%soil_water_content)%v)%data_3d, 2)
        elseif (this%var_indx(kVARS%soil_temperature)%v > 0) then
            nsoil = size(this%vars_3d(this%var_indx(kVARS%soil_temperature)%v)%data_3d, 2)
        endif

        if (options%domain%landvar /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%landvar,         &
                           temporary_data)
            if (this%var_indx(kVARS%land_mask)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%land_mask)%v)%data_2di = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                where(this%vars_2d(this%var_indx(kVARS%land_mask)%v)%data_2di==0) this%vars_2d(this%var_indx(kVARS%land_mask)%v)%data_2di = kLC_WATER  ! To ensure conisitency. land_mask can be 0 or 2 for water, enforce a single value.
            else 
                if (this%var_indx(kVARS%land_mask)%v > 0) then
                    this%vars_2d(this%var_indx(kVARS%land_mask)%v)%data_2di = kLC_LAND
                endif
            endif
        endif

        if ((options%physics%watersurface==kWATER_LAKE) .AND.(options%domain%lakedepthvar /= "")) then
            if (STD_OUT_PE) write(*,*) "   reading lake depth data from hi-res file"

            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lakedepthvar,         &
                           temporary_data)
            if (this%var_indx(kVARS%lake_depth)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%lake_depth)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif

        endif

        if (options%domain%soiltype_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%soiltype_var,         &
                           temporary_data)
            if (this%var_indx(kVARS%soil_type)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%soil_type)%v)%data_2di = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%soil_type)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%soil_type)%v)%data_2di = 3
            endif
        endif

        if (options%domain%cropcategory_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%cropcategory_var,         &
                           temporary_data)
            if (this%var_indx(kVARS%crop_category)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%crop_category)%v)%data_2di = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%crop_category)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%crop_category)%v)%data_2di = 0
            endif
        endif


        if (options%domain%soil_deept_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%soil_deept_var,       &
                           temporary_data)
            if (this%var_indx(kVARS%soil_deep_temperature)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)

                if (minval(temporary_data)< 200) then
                    if (STD_OUT_PE) write(*,*) "WARNING, VERY COLD SOIL TEMPERATURES SPECIFIED:", minval(temporary_data)
                    if (STD_OUT_PE) write(*,*) trim(options%domain%init_conditions_file),"  ",trim(options%domain%soil_deept_var)
                endif
                if (minval(this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d)< 200) then
                    where(this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d<200) this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d=init_surf_temp ! <200 is just broken, set to mean annual air temperature at mid-latidudes
                endif
            endif
        else
            if (this%var_indx(kVARS%soil_deep_temperature)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d = init_surf_temp
            endif
        endif

        if (options%domain%soil_t_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%soil_t_var,           &
                           temporary_data_3d)
            if (this%var_indx(kVARS%soil_temperature)%v > 0) then
                do i=1,nsoil
                    this%vars_3d(this%var_indx(kVARS%soil_temperature)%v)%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                enddo
                if (options%domain%soil_deept_var == "") then
                    if (this%var_indx(kVARS%soil_deep_temperature)%v > 0) then
                        this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d = this%vars_3d(this%var_indx(kVARS%soil_temperature)%v)%data_3d(:,nsoil,:)
                    endif
                endif
            endif

        else
            if (this%var_indx(kVARS%soil_temperature)%v > 0) then
                if (this%var_indx(kVARS%soil_deep_temperature)%v > 0) then
                    do i=1,nsoil
                        this%vars_3d(this%var_indx(kVARS%soil_temperature)%v)%data_3d(:,i,:) = this%vars_2d(this%var_indx(kVARS%soil_deep_temperature)%v)%data_2d
                    enddo
                endif
            endif
        endif


        if (options%domain%swe_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%swe_var,         &
                           temporary_data)
            if (this%var_indx(kVARS%snow_water_equivalent)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%snow_water_equivalent)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = 0
            endif
        endif

        if (options%domain%snowh_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%snowh_var,         &
                           temporary_data)
            if (this%var_indx(kVARS%snow_height)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%snow_height)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%snow_height)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%snow_height)%v)%data_2d = 0
            endif
        endif
        
        if ( (this%var_indx(kVARS%snow_height)%v > 0) .and. (this%var_indx(kVARS%snow_water_equivalent)%v > 0)) then
            !Do check if we read in SWE but not snow height -- convert with user supplied constant density
            if (options%domain%swe_var /= "" .and. options%domain%snowh_var == "") then
                this%vars_2d(this%var_indx(kVARS%snow_height)%v)%data_2d = this%vars_2d(this%var_indx(kVARS%snow_water_equivalent)%v)%data_2d/options%lsm%snow_den_const
            endif
            
            !Do check if we read in snow height but not SWE -- convert with user supplied constant density
            if (options%domain%snowh_var /= "" .and. options%domain%swe_var == "") then
                this%vars_2d(this%var_indx(kVARS%snow_water_equivalent)%v)%data_2d = this%vars_2d(this%var_indx(kVARS%snow_height)%v)%data_2d*options%lsm%snow_den_const
            endif
        endif



        if (options%domain%soil_vwc_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%soil_vwc_var,         &
                           temporary_data_3d)
            if (this%var_indx(kVARS%soil_water_content)%v > 0) then
                do i=1,nsoil
                    this%vars_3d(this%var_indx(kVARS%soil_water_content)%v)%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                enddo
            endif

        else
            if (this%var_indx(kVARS%soil_water_content)%v > 0) then
                this%vars_3d(this%var_indx(kVARS%soil_water_content)%v)%data_3d = 0.4
            endif
        endif

        if (options%domain%vegtype_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%vegtype_var,          &
                           temporary_data)
            if (this%var_indx(kVARS%veg_type)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%veg_type)%v)%data_2di = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%veg_type)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%veg_type)%v)%data_2di = 7
            endif
        endif

        if (options%domain%albedo_var /= "") then
            if (options%lsm%monthly_albedo) then
                call io_read(options%domain%init_conditions_file,   &
                            options%domain%albedo_var,          &
                            temporary_data_3d)

                if (this%var_indx(kVARS%albedo)%v > 0) then
                    do i=1,size(this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d, 2)
                        this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme,i)
                    enddo

                    if (maxval(temporary_data_3d) > 1) then
                        if (STD_OUT_PE) write(*,*) "Changing input ALBEDO % to fraction"
                        this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d = this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d / 100
                    endif
                endif
            else
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%albedo_var,          &
                               temporary_data)
                if (this%var_indx(kVARS%albedo)%v > 0) then
                    do i=1,size(this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d, 2)
                        this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d(:,i,:) = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                    enddo

                    if (maxval(temporary_data) > 1) then
                        if (STD_OUT_PE) write(*,*) "Changing input ALBEDO % to fraction"
                        this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d = this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d / 100
                    endif
                endif
            endif

        else
            if (this%var_indx(kVARS%albedo)%v > 0) then
                this%vars_3d(this%var_indx(kVARS%albedo)%v)%data_3d = 0.17
            endif
        endif


        if (options%domain%vegfrac_var /= "") then
            if (options%lsm%monthly_albedo) then
                call io_read(options%domain%init_conditions_file,   &
                            options%domain%vegfrac_var,          &
                            temporary_data_3d)

                if (this%var_indx(kVARS%vegetation_fraction)%v > 0) then
                    do i=1,size(this%vars_3d(this%var_indx(kVARS%vegetation_fraction)%v)%data_3d, 2)
                        this%vars_3d(this%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme,i)
                    enddo
                endif
            else
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%vegfrac_var,          &
                               temporary_data)
                if (this%var_indx(kVARS%vegetation_fraction)%v > 0) then
                    do i=1,size(this%vars_3d(this%var_indx(kVARS%vegetation_fraction)%v)%data_3d, 2)
                        this%vars_3d(this%var_indx(kVARS%vegetation_fraction)%v)%data_3d(:,i,:) = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                    enddo
                endif
            endif

        else
            if (this%var_indx(kVARS%vegetation_fraction)%v > 0) then
                this%vars_3d(this%var_indx(kVARS%vegetation_fraction)%v)%data_3d = 60.
            endif
        endif

        if (this%var_indx(kVARS%soil_totalmoisture)%v > 0) then
            this%vars_2d(this%var_indx(kVARS%soil_totalmoisture)%v)%data_2d = 0
            if (this%var_indx(kVARS%soil_water_content)%v > 0) then
                do i=1, nsoil
                    this%vars_2d(this%var_indx(kVARS%soil_totalmoisture)%v)%data_2d = this%vars_2d(this%var_indx(kVARS%soil_totalmoisture)%v)%data_2d + this%vars_3d(this%var_indx(kVARS%soil_water_content)%v)%data_3d(:,i,:) * soil_thickness(i) * 1000 !! MJ added
                enddo
            endif
        endif

        if (options%domain%vegfracmax_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%vegfracmax_var,       &
                           temporary_data)
            if (this%var_indx(kVARS%vegetation_fraction_max)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%vegetation_fraction_max)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%vegetation_fraction_max)%v > 0) then
                if (STD_OUT_PE) write(*,*) "    VEGMAX not specified; using default value of 0.8"
                this%vars_2d(this%var_indx(kVARS%vegetation_fraction_max)%v)%data_2d = 80.
            endif
        endif

        if (options%domain%lai_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%lai_var,              &
                           temporary_data)
            if (this%var_indx(kVARS%lai)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%lai)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%lai)%v > 0) then
                if (STD_OUT_PE) write(*,*) "    LAI not specified; using default value of 1"
                this%vars_2d(this%var_indx(kVARS%lai)%v)%data_2d = 1
            endif
        endif

        if (options%domain%canwat_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%canwat_var,              &
                           temporary_data)
            if (this%var_indx(kVARS%canopy_water)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%canopy_water)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        else
            if (this%var_indx(kVARS%canopy_water)%v > 0) then
                if (STD_OUT_PE) write(*,*) "    CANWAT not specified; using default value of 0"
                this%vars_2d(this%var_indx(kVARS%canopy_water)%v)%data_2d = 0
            endif
        endif
        if (options%domain%shd_var /= "") then
            call io_read(options%domain%init_conditions_file,   &
                           options%domain%shd_var,       &
                           temporary_data)
            if (this%var_indx(kVARS%shd)%v > 0) then
                this%vars_2d(this%var_indx(kVARS%shd)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
            endif
        endif

        if (options%physics%radiation_downScaling==1) then
            !!
            if (options%domain%hlm_var /= "") then
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%hlm_var,           &
                               temporary_data_3d)
                if (this%var_indx(kVARS%hlm)%v > 0) then
                    do i=1,90
                        this%vars_3d(this%var_indx(kVARS%hlm)%v)%data_3d(:,i,:) = temporary_data_3d(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme, i)
                        !if (STD_OUT_PE) write(*,*),"hlm ", i, this%vars_3d(this%var_indx(kVARS%hlm)%v)%data_3d(this%grid%its,i,this%grid%jts)
                    enddo
                endif
            else  
                stop "hlm_var not specified in domain file, but required for radiation downscaling"
            endif
            !!
            if (options%domain%svf_var /= "") then
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%svf_var,         &
                               temporary_data)
                if (this%var_indx(kVARS%svf)%v > 0) then
                    this%vars_2d(this%var_indx(kVARS%svf)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            else
                stop "svf_var not specified in domain file, but required for radiation downscaling"
            endif            
            !!
            if (options%domain%slope_var /= "") then
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%slope_var,         &
                               temporary_data)
                if (this%var_indx(kVARS%slope)%v > 0) then
                    this%vars_2d(this%var_indx(kVARS%slope)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            else
                stop "slope_var not specified in domain file, but required for radiation downscaling"
            endif            
            !!
            if (options%domain%slope_angle_var /= "") then
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%slope_angle_var,         &
                               temporary_data)
                if (this%var_indx(kVARS%slope_angle)%v > 0) then
                    this%vars_2d(this%var_indx(kVARS%slope_angle)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            else
                stop "slope_angle_var not specified in domain file, but required for radiation downscaling"
            endif
            !!
            if (options%domain%aspect_angle_var /= "") then
                call io_read(options%domain%init_conditions_file,   &
                               options%domain%aspect_angle_var,         &
                               temporary_data)
                if (this%var_indx(kVARS%aspect_angle)%v > 0) then
                    this%vars_2d(this%var_indx(kVARS%aspect_angle)%v)%data_2d = temporary_data(this%grid%ims:this%grid%ime, this%grid%jms:this%grid%jme)
                endif
            else
                stop "aspect_angle_var not specified in domain file, but required for radiation downscaling"
            endif
        endif

        ! these will all be udpated by either forcing data or the land model, but initialize to sensible values to avoid breaking other initialization routines
        if (this%var_indx(kVARS%skin_temperature)%v > 0) this%vars_2d(this%var_indx(kVARS%skin_temperature)%v)%data_2d = init_surf_temp
        if (this%var_indx(kVARS%sst)%v > 0) this%vars_2d(this%var_indx(kVARS%sst)%v)%data_2d = init_surf_temp
        if (this%var_indx(kVARS%roughness_z0)%v > 0) this%vars_2d(this%var_indx(kVARS%roughness_z0)%v)%data_2d = 0.001
        if (this%var_indx(kVARS%sensible_heat)%v > 0) this%vars_2d(this%var_indx(kVARS%sensible_heat)%v)%data_2d=0
        if (this%var_indx(kVARS%latent_heat)%v > 0) this%vars_2d(this%var_indx(kVARS%latent_heat)%v)%data_2d=0
        if (this%var_indx(kVARS%u_10m)%v > 0) this%vars_2d(this%var_indx(kVARS%u_10m)%v)%data_2d=0
        if (this%var_indx(kVARS%v_10m)%v > 0) this%vars_2d(this%var_indx(kVARS%v_10m)%v)%data_2d=0

        if (this%var_indx(kVARS%windspd_10m)%v > 0) this%vars_2d(this%var_indx(kVARS%windspd_10m)%v)%data_2d=0
        if (this%var_indx(kVARS%temperature_2m)%v > 0) this%vars_2d(this%var_indx(kVARS%temperature_2m)%v)%data_2d=init_surf_temp
        if (this%var_indx(kVARS%humidity_2m)%v > 0) this%vars_2d(this%var_indx(kVARS%humidity_2m)%v)%data_2d=0.001
        if (this%var_indx(kVARS%surface_pressure)%v > 0) this%vars_2d(this%var_indx(kVARS%surface_pressure)%v)%data_2d=102000
        if (this%var_indx(kVARS%longwave_up)%v > 0) this%vars_2d(this%var_indx(kVARS%longwave_up)%v)%data_2d=0
        if (this%var_indx(kVARS%ground_heat_flux)%v > 0) this%vars_2d(this%var_indx(kVARS%ground_heat_flux)%v)%data_2d=0
        if (this%var_indx(kVARS%veg_leaf_temperature)%v > 0) this%vars_2d(this%var_indx(kVARS%veg_leaf_temperature)%v)%data_2d=init_surf_temp
        if (this%var_indx(kVARS%ground_surf_temperature)%v > 0) this%vars_2d(this%var_indx(kVARS%ground_surf_temperature)%v)%data_2d=init_surf_temp
        if (this%var_indx(kVARS%canopy_vapor_pressure)%v > 0) this%vars_2d(this%var_indx(kVARS%canopy_vapor_pressure)%v)%data_2d=2000
        if (this%var_indx(kVARS%canopy_temperature)%v > 0) this%vars_2d(this%var_indx(kVARS%canopy_temperature)%v)%data_2d=init_surf_temp
        if (this%var_indx(kVARS%coeff_momentum_drag)%v > 0) this%vars_2d(this%var_indx(kVARS%coeff_momentum_drag)%v)%data_2d=0.01
        if (this%var_indx(kVARS%chs)%v > 0) this%vars_2d(this%var_indx(kVARS%chs)%v)%data_2d=0.01
        if (this%var_indx(kVARS%chs2)%v > 0) this%vars_2d(this%var_indx(kVARS%chs2)%v)%data_2d=0.01
        if (this%var_indx(kVARS%cqs2)%v > 0) this%vars_2d(this%var_indx(kVARS%cqs2)%v)%data_2d=0.01

        if (this%var_indx(kVARS%QFX)%v > 0) this%vars_2d(this%var_indx(kVARS%QFX)%v)%data_2d=0.0
        if (this%var_indx(kVARS%br)%v > 0) this%vars_2d(this%var_indx(kVARS%br)%v)%data_2d=0.0
        if (this%var_indx(kVARS%mol)%v > 0) this%vars_2d(this%var_indx(kVARS%mol)%v)%data_2d=0.0
        if (this%var_indx(kVARS%psim)%v > 0) this%vars_2d(this%var_indx(kVARS%psim)%v)%data_2d=0.0
        if (this%var_indx(kVARS%psih)%v > 0) this%vars_2d(this%var_indx(kVARS%psih)%v)%data_2d=0.0
        if (this%var_indx(kVARS%fm)%v > 0) this%vars_2d(this%var_indx(kVARS%fm)%v)%data_2d=0.0
        if (this%var_indx(kVARS%fh)%v > 0) this%vars_2d(this%var_indx(kVARS%fh)%v)%data_2d=0.0

        if (this%var_indx(kVARS%hpbl)%v > 0) this%vars_2d(this%var_indx(kVARS%hpbl)%v)%data_2d=100.0
        if (this%var_indx(kVARS%coeff_heat_exchange_3d)%v > 0) this%vars_3d(this%var_indx(kVARS%coeff_heat_exchange_3d)%v)%data_3d=0.01
        if (this%var_indx(kVARS%coeff_momentum_exchange_3d)%v > 0) this%vars_3d(this%var_indx(kVARS%coeff_momentum_exchange_3d)%v)%data_3d=0.01
        if (this%var_indx(kVARS%canopy_fwet)%v > 0) this%vars_2d(this%var_indx(kVARS%canopy_fwet)%v)%data_2d=0
        if (this%var_indx(kVARS%snow_water_eq_prev)%v > 0) this%vars_2d(this%var_indx(kVARS%snow_water_eq_prev)%v)%data_2d=0
        if (this%var_indx(kVARS%snow_albedo_prev)%v > 0) this%vars_2d(this%var_indx(kVARS%snow_albedo_prev)%v)%data_2d=0.65
        if (this%var_indx(kVARS%storage_lake)%v > 0) this%vars_2d(this%var_indx(kVARS%storage_lake)%v)%data_2d=0

        if (this%var_indx(kVARS%ustar)%v > 0)                  this%vars_2d(this%var_indx(kVARS%ustar)%v)%data_2d=0.1
        if (this%var_indx(kVARS%irr_eventno_sprinkler)%v > 0)  this%vars_2d(this%var_indx(kVARS%irr_eventno_sprinkler)%v)%data_2di=0
        if (this%var_indx(kVARS%irr_eventno_micro)%v > 0)      this%vars_2d(this%var_indx(kVARS%irr_eventno_micro)%v)%data_2di=0
        if (this%var_indx(kVARS%irr_eventno_flood)%v > 0)      this%vars_2d(this%var_indx(kVARS%irr_eventno_flood)%v)%data_2di=0
        if (this%var_indx(kVARS%plant_growth_stage)%v > 0)     this%vars_2d(this%var_indx(kVARS%plant_growth_stage)%v)%data_2di=0
        if (this%var_indx(kVARS%kpbl)%v > 0)                   this%vars_2d(this%var_indx(kVARS%kpbl)%v)%data_2di=0

        if (this%var_indx(kVARS%runoff_tstep)%v > 0)        this%vars_2d(this%var_indx(kVARS%runoff_tstep)%v)%data_2d=0.
        if (this%var_indx(kVARS%snow_temperature)%v > 0)    this%vars_3d(this%var_indx(kVARS%snow_temperature)%v)%data_3d=273.15
        if (this%var_indx(kVARS%Sice)%v > 0)                this%vars_3d(this%var_indx(kVARS%Sice)%v)%data_3d=0.
        if (this%var_indx(kVARS%Sliq)%v > 0)                this%vars_3d(this%var_indx(kVARS%Sliq)%v)%data_3d=0.
        if (this%var_indx(kVARS%Ds)%v > 0)                  this%vars_3d(this%var_indx(kVARS%Ds)%v)%data_3d=0.
        if (this%var_indx(kVARS%fsnow)%v > 0)               this%vars_2d(this%var_indx(kVARS%fsnow)%v)%data_2d=0.
        if (this%var_indx(kVARS%Nsnow)%v > 0)               this%vars_2d(this%var_indx(kVARS%Nsnow)%v)%data_2d=0.
        if (this%var_indx(kVARS%dSWE_salt)%v > 0)             this%vars_2d(this%var_indx(kVARS%dSWE_salt)%v)%data_2d=0.
        if (this%var_indx(kVARS%dSWE_susp)%v > 0)             this%vars_2d(this%var_indx(kVARS%dSWE_susp)%v)%data_2d=0.
        if (this%var_indx(kVARS%dSWE_subl)%v > 0)             this%vars_2d(this%var_indx(kVARS%dSWE_subl)%v)%data_2d=0.
        if (this%var_indx(kVARS%dSWE_slide)%v > 0)            this%vars_2d(this%var_indx(kVARS%dSWE_slide)%v)%data_2d=0.

        !!
        if (this%var_indx(kVARS%meltflux_out_tstep)%v > 0)  this%vars_2d(this%var_indx(kVARS%meltflux_out_tstep)%v)%data_2d=0.
        if (this%var_indx(kVARS%shortwave_direct)%v > 0)  this%vars_2d(this%var_indx(kVARS%shortwave_direct)%v)%data_2d=0.
        if (this%var_indx(kVARS%shortwave_diffuse)%v > 0)  this%vars_2d(this%var_indx(kVARS%shortwave_diffuse)%v)%data_2d=0.
        if (this%var_indx(kVARS%shortwave_direct_above)%v > 0)  this%vars_2d(this%var_indx(kVARS%shortwave_direct_above)%v)%data_2d=0.
        if (this%var_indx(kVARS%Sliq_out)%v > 0)  this%vars_2d(this%var_indx(kVARS%Sliq_out)%v)%data_2d=0.

    end subroutine read_land_variables

    !> -------------------------------
    !! Read in the shape of the domain required and setup the grid objects
    !!
    !! -------------------------------
    subroutine read_domain_shape(this, options)
        implicit none
        class(domain_t), intent(inout)  :: this
        type(options_t), intent(in)     :: options

        real, allocatable :: temporary_data(:,:)
        integer :: nx_global, ny_global, nz_global, nsmooth, adv_order, my_index

        nsmooth = max(1, int(options%wind%smooth_wind_distance / options%domain%dx))
        if (options%wind%smooth_wind_distance == 0.0) nsmooth = 0
        this%nsmooth = nsmooth
        if ((STD_OUT_PE).and.(options%general%debug)) write(*,*) "number of gridcells to smooth = ",nsmooth
        ! This doesn't need to read in this variable, it could just request the dimensions
        ! but this is not a performance sensitive part of the code (for now)
        call io_read(options%domain%init_conditions_file,   &
                     options%domain%hgt_hi,                 &
                     temporary_data)

        nx_global = size(temporary_data,1)
        ny_global = size(temporary_data,2)
        nz_global = options%domain%nz
        
        adv_order = options%adv%h_order
        
        !If we are using the monotonic flux limiter, it is necesarry to calculate the fluxes one location deep into the
        !halo. Thus, we need one extra cell in each halo direction to support the finite difference stencil
        !This is achieved here by artificially inflating the adv_order which is passed to the grid setup
        if (options%adv%flux_corr==kFLUXCOR_MONO) adv_order = adv_order+2
        
        !If using MPDATA, we need a halo of size 2 to support the difference stencil
        if (options%physics%advection==kADV_MPDATA) adv_order = 4
        
        if (this%compute_comms == MPI_COMM_NULL) then
            my_index = 1
        else
            call MPI_Comm_rank(this%compute_comms, my_index)
            ! MPI returns rank, which is 0-indexed
            my_index = my_index + 1
        endif

        call this%grid%set_grid_dimensions(     nx_global, ny_global, nz_global, image=my_index, comms=this%compute_comms, adv_order=adv_order)
        call this%grid8w%set_grid_dimensions(   nx_global, ny_global, nz_global+1, image=my_index, comms=this%compute_comms, adv_order=adv_order)

        call this%u_grid%set_grid_dimensions( nx_global, ny_global, nz_global, image=my_index, comms=this%compute_comms, adv_order=adv_order, nx_extra = 1)
        call this%v_grid%set_grid_dimensions( nx_global, ny_global, nz_global, image=my_index, comms=this%compute_comms, adv_order=adv_order, ny_extra = 1)

        ! for 2D mass variables
        call this%grid2d%set_grid_dimensions( nx_global, ny_global, 0, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)

        ! setup a 2D lat/lon grid extended by nsmooth grid cells so that smoothing can take place "across" images
        ! This just sets up the fields to interpolate u and v to so that the input data are handled on an extended
        ! grid.  They are then subset to the u_grid and v_grids above before actual use.
        call this%u_grid2d%set_grid_dimensions(     nx_global, ny_global, 0, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order, nx_extra = 1)

        ! handle the v-grid too
        call this%v_grid2d%set_grid_dimensions(     nx_global, ny_global, 0, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order, ny_extra = 1)
        
        call this%column_grid%set_grid_dimensions(               0,         0, nz_global, image=my_index, comms=this%compute_comms, adv_order=adv_order) !! MJ added
        call this%grid_soil%set_grid_dimensions(         nx_global, ny_global, kSOIL_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_snow%set_grid_dimensions(         nx_global, ny_global, kSNOW_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_snowsoil%set_grid_dimensions(     nx_global, ny_global, kSNOWSOIL_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_soilcomp%set_grid_dimensions(     nx_global, ny_global, kSOILCOMP_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_gecros%set_grid_dimensions(       nx_global, ny_global, kGECROS_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_croptype%set_grid_dimensions(     nx_global, ny_global, kCROP_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_monthly%set_grid_dimensions(      nx_global, ny_global, kMONTH_GRID_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_lake%set_grid_dimensions(         nx_global, ny_global, kLAKE_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_lake_soisno%set_grid_dimensions(  nx_global, ny_global, kLAKE_SOISNO_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_lake_soi%set_grid_dimensions(     nx_global, ny_global, kLAKE_SOI_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_lake_soisno_1%set_grid_dimensions(nx_global, ny_global, kLAKE_SOISNO_1_Z, image=my_index, comms=this%compute_comms, global_nz=nz_global, adv_order=adv_order)
        call this%grid_hlm%set_grid_dimensions(     nx_global, ny_global, 90, image=my_index, comms=this%compute_comms, adv_order=adv_order) !! MJ added
        call this%grid_Sx%set_grid_dimensions(     nx_global, ny_global, nz_global, 72, image=my_index, comms=this%compute_comms, adv_order=adv_order) !! MJ added

        call this%global_grid_2d%set_grid_dimensions(   nx_global, ny_global, 0)
        call this%global_grid%set_grid_dimensions(   nx_global, ny_global, nz_global)
        call this%global_grid8w%set_grid_dimensions(   nx_global, ny_global, nz_global+1)

        ! We need to manually set the neighbor grid bounds below, this is just to initialize the grid
        call this%neighbor_grid_2d%set_grid_dimensions(   nx_global, ny_global, 0)
        call this%neighbor_grid%set_grid_dimensions(   nx_global, ny_global, nz_global)
        call this%neighbor_grid8w%set_grid_dimensions(   nx_global, ny_global, nz_global+1)

        this%ximg = this%grid%ximg
        this%ximages = this%grid%ximages
        this%yimg = this%grid%yimg
        this%yimages = this%grid%yimages

        this%north_boundary = (this%grid%yimg == this%grid%yimages)
        this%south_boundary = (this%grid%yimg == 1)
        this%east_boundary  = (this%grid%ximg == this%grid%ximages)
        this%west_boundary  = (this%grid%ximg == 1)

        this%ims = this%grid%ims; this%its = this%grid%its; this%ids = this%grid%ids
        this%ime = this%grid%ime; this%ite = this%grid%ite; this%ide = this%grid%ide
        this%kms = this%grid%kms; this%kts = this%grid%kts; this%kds = this%grid%kds
        this%kme = this%grid%kme; this%kte = this%grid%kte; this%kde = this%grid%kde
        this%jms = this%grid%jms; this%jts = this%grid%jts; this%jds = this%grid%jds
        this%jme = this%grid%jme; this%jte = this%grid%jte; this%jde = this%grid%jde
        
        !Calculate neighborhood indexes. These are used to store terrain in the local neighborhood for non-local wind calculations
        this%neighborhood_max = max(nsmooth,8)
        
        !Considering blocking terrain...
        if (options%physics%windtype == kLINEAR_ITERATIVE_WINDS .or. options%physics%windtype ==  kITERATIVE_WINDS) then
            this%neighborhood_max = int(max(4000.0/this%dx,1.0))
        endif
        
        !Considering TPI...
        if (options%wind%Sx) then
            this%neighborhood_max = max(this%neighborhood_max,floor(max(1.0,(options%wind%TPI_dmax+options%wind%Sx_dmax)/this%dx)))
        endif
        
        this%ihs=max(this%grid%ims-this%neighborhood_max,this%grid%ids); this%ihe=min(this%grid%ime+this%neighborhood_max,this%grid%ide)
        this%jhs=max(this%grid%jms-this%neighborhood_max,this%grid%jds); this%jhe=min(this%grid%jme+this%neighborhood_max,this%grid%jde)
        this%khs=this%grid%kms;                                          this%khe=this%grid%kme

        this%neighbor_grid_2d%ims=this%ihs; this%neighbor_grid_2d%ime=this%ihe
        this%neighbor_grid_2d%jms=this%jhs; this%neighbor_grid_2d%jme=this%jhe
        this%neighbor_grid%ims=this%ihs; this%neighbor_grid%ime=this%ihe
        this%neighbor_grid%jms=this%jhs; this%neighbor_grid%jme=this%jhe
        this%neighbor_grid8w%ims=this%ihs; this%neighbor_grid8w%ime=this%ihe
        this%neighbor_grid8w%jms=this%jhs; this%neighbor_grid8w%jme=this%jhe


    end subroutine

    !> -------------------------------
    !! Check that a set of variables is within realistic bounds (i.e. >0)
    !!
    !! Need to add more variables to the list
    !!
    !! -------------------------------
    module subroutine enforce_limits(this,update_in)
      class(domain_t), intent(inout) :: this
      logical, optional,  intent(in) :: update_in

      logical update
      update = .False.
      if (present(update_in)) update = update_in

      if (update) then
        if (this%var_indx(kVARS%water_vapor)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%dqdt_3d)      ) then 
                where(this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%dqdt_3d < 0)           this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%dqdt_3d = 0
            endif
        endif

        if (this%var_indx(kVARS%potential_temperature)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%dqdt_3d) ) then
                where(this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%dqdt_3d < 0) this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%dqdt_3d = 0
            endif
        endif

        if (this%var_indx(kVARS%cloud_water_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%dqdt_3d) ) then
                where(this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%dqdt_3d < 0)      this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%dqdt_3d = 0
            endif
        endif

        if (this%var_indx(kVARS%cloud_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%cloud_number)%v)%dqdt_3d)    ) then
                where(this%vars_3d(this%var_indx(kVARS%cloud_number)%v)%dqdt_3d < 0)          this%vars_3d(this%var_indx(kVARS%cloud_number)%v)%dqdt_3d = 0
            endif 
        endif
        if (this%var_indx(kVARS%ice_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%dqdt_3d)  ) then
                where(this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%dqdt_3d < 0)        this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%dqdt_3d = 0
            endif
        endif 
        if (this%var_indx(kVARS%ice_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice_number)%v)%dqdt_3d)) then
                where(this%vars_3d(this%var_indx(kVARS%ice_number)%v)%dqdt_3d < 0)      this%vars_3d(this%var_indx(kVARS%ice_number)%v)%dqdt_3d = 0
            endif
        endif 
        if (this%var_indx(kVARS%rain_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%dqdt_3d)       ) then
                where(this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%dqdt_3d < 0)             this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%dqdt_3d = 0
            endif
        endif 
        if (this%var_indx(kVARS%rain_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%rain_number)%v)%dqdt_3d)     ) then
                where(this%vars_3d(this%var_indx(kVARS%rain_number)%v)%dqdt_3d < 0)           this%vars_3d(this%var_indx(kVARS%rain_number)%v)%dqdt_3d = 0
            endif 
        endif
        if (this%var_indx(kVARS%snow_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%dqdt_3d)       ) then
                where(this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%dqdt_3d < 0)             this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%dqdt_3d = 0
            endif 
        endif
        if (this%var_indx(kVARS%snow_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%snow_number)%v)%dqdt_3d)     ) then
                where(this%vars_3d(this%var_indx(kVARS%snow_number)%v)%dqdt_3d < 0)           this%vars_3d(this%var_indx(kVARS%snow_number)%v)%dqdt_3d = 0
            endif 
        endif
        if (this%var_indx(kVARS%graupel_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%dqdt_3d)    ) then
                where(this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%dqdt_3d < 0)          this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%graupel_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%graupel_number)%v)%dqdt_3d)  ) then
                where(this%vars_3d(this%var_indx(kVARS%graupel_number)%v)%dqdt_3d < 0)        this%vars_3d(this%var_indx(kVARS%graupel_number)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice1_a)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice1_a)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice1_a)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice1_a)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice1_c)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice1_c)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice1_c)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice1_c)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice2_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice2_mass)%v)%dqdt_3d)       ) then
                where(this%vars_3d(this%var_indx(kVARS%ice2_mass)%v)%dqdt_3d < 0)             this%vars_3d(this%var_indx(kVARS%ice2_mass)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice2_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice2_number)%v)%dqdt_3d)     ) then
                where(this%vars_3d(this%var_indx(kVARS%ice2_number)%v)%dqdt_3d < 0)           this%vars_3d(this%var_indx(kVARS%ice2_number)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice2_a)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice2_a)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice2_a)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice2_a)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice2_c)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice2_c)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice2_c)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice2_c)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice3_mass)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice3_mass)%v)%dqdt_3d)       ) then
                where(this%vars_3d(this%var_indx(kVARS%ice3_mass)%v)%dqdt_3d < 0)             this%vars_3d(this%var_indx(kVARS%ice3_mass)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice3_number)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice3_number)%v)%dqdt_3d)     ) then
                where(this%vars_3d(this%var_indx(kVARS%ice3_number)%v)%dqdt_3d < 0)           this%vars_3d(this%var_indx(kVARS%ice3_number)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice3_a)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice3_a)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice3_a)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice3_a)%v)%dqdt_3d = 0
            endif
        endif
        if (this%var_indx(kVARS%ice3_c)%v > 0) then
            if (allocated(this%vars_3d(this%var_indx(kVARS%ice3_c)%v)%dqdt_3d)          ) then
                where(this%vars_3d(this%var_indx(kVARS%ice3_c)%v)%dqdt_3d < 0)                this%vars_3d(this%var_indx(kVARS%ice3_c)%v)%dqdt_3d = 0
            endif
        endif
    else
        if (this%var_indx(kVARS%water_vapor)%v > 0           ) where(this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%data_3d < 0)             this%vars_3d(this%var_indx(kVARS%water_vapor)%v)%data_3d = 0
        if (this%var_indx(kVARS%potential_temperature)%v > 0 ) where(this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%data_3d < 0)   this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%data_3d = 0
        if (this%var_indx(kVARS%cloud_water_mass)%v > 0      ) where(this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%data_3d < 0)        this%vars_3d(this%var_indx(kVARS%cloud_water_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%cloud_number)%v > 0          ) where(this%vars_3d(this%var_indx(kVARS%cloud_number)%v)%data_3d < 0)            this%vars_3d(this%var_indx(kVARS%cloud_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice_mass)%v > 0        ) where(this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%data_3d < 0)          this%vars_3d(this%var_indx(kVARS%ice_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice_number)%v > 0      ) where(this%vars_3d(this%var_indx(kVARS%ice_number)%v)%data_3d < 0)        this%vars_3d(this%var_indx(kVARS%ice_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%rain_mass)%v > 0             ) where(this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%data_3d < 0)               this%vars_3d(this%var_indx(kVARS%rain_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%rain_number)%v > 0           ) where(this%vars_3d(this%var_indx(kVARS%rain_number)%v)%data_3d < 0)             this%vars_3d(this%var_indx(kVARS%rain_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%snow_mass)%v > 0             ) where(this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%data_3d < 0)               this%vars_3d(this%var_indx(kVARS%snow_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%snow_number)%v > 0           ) where(this%vars_3d(this%var_indx(kVARS%snow_number)%v)%data_3d < 0)             this%vars_3d(this%var_indx(kVARS%snow_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%graupel_mass)%v > 0          ) where(this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%data_3d < 0)            this%vars_3d(this%var_indx(kVARS%graupel_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%graupel_number)%v > 0        ) where(this%vars_3d(this%var_indx(kVARS%graupel_number)%v)%data_3d < 0)          this%vars_3d(this%var_indx(kVARS%graupel_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice1_a)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice1_a)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice1_a)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice1_c)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice1_c)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice1_c)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice2_mass)%v > 0             ) where(this%vars_3d(this%var_indx(kVARS%ice2_mass)%v)%data_3d < 0)               this%vars_3d(this%var_indx(kVARS%ice2_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice2_number)%v > 0           ) where(this%vars_3d(this%var_indx(kVARS%ice2_number)%v)%data_3d < 0)             this%vars_3d(this%var_indx(kVARS%ice2_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice2_a)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice2_a)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice2_a)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice2_c)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice2_c)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice2_c)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice3_mass)%v > 0             ) where(this%vars_3d(this%var_indx(kVARS%ice3_mass)%v)%data_3d < 0)               this%vars_3d(this%var_indx(kVARS%ice3_mass)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice3_number)%v > 0           ) where(this%vars_3d(this%var_indx(kVARS%ice3_number)%v)%data_3d < 0)             this%vars_3d(this%var_indx(kVARS%ice3_number)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice3_a)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice3_a)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice3_a)%v)%data_3d = 0
        if (this%var_indx(kVARS%ice3_c)%v > 0                ) where(this%vars_3d(this%var_indx(kVARS%ice3_c)%v)%data_3d < 0)                  this%vars_3d(this%var_indx(kVARS%ice3_c)%v)%data_3d = 0

      endif
    end subroutine


    !> -------------------------------
    !! Setup the Geographic look up tables for interpolating a given forcing data set to each of the grids
    !!
    !! -------------------------------
    subroutine setup_geo_interpolation(this, forcing, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(options_t), intent(in)     :: options

        type(interpolable_type) :: forc_u_from_mass, forc_v_from_mass
        type(variable_t) :: temporary_data
        integer :: nx, ny, nz, i, j, k, ims, ime, jms, jme
        real, allocatable, dimension(:,:) :: AGL_cap, AGL_u_cap, AGL_v_cap, AGL_n, AGL_u_n, AGL_v_n

        ! this%geo and forcing%geo have to be of class interpolable
        ! which means they must contain lat, lon, z, geolut, and vLUT components

        call geo_LUT(this%geo_agl,forcing%geo_agl)
        call geo_LUT(this%geo,forcing%geo)
        call geo_LUT(this%geo_u,  forcing%geo_u)
        call geo_LUT(this%geo_v,  forcing%geo_v)

        if (allocated(forcing%z)) then  ! In case of external 2D forcing data, skip the VLUTs.

            ! See if we did not set forcing z at initialization of boundary object. This can happen
            ! if we are a nest, and the parent domain's global z field was not yet set at the time of
            ! initialization.
            if (.not.forcing%z_is_set) then
                temporary_data = forcing%variables%get_var(options%forcing%zvar)
                forcing%z = temporary_data%data_3d
            endif

            forc_u_from_mass%lat = forcing%geo%lat
            forc_u_from_mass%lon = forcing%geo%lon
            forc_v_from_mass%lat = forcing%geo%lat
            forc_v_from_mass%lon = forcing%geo%lon

            call geo_LUT(this%geo_u, forc_u_from_mass)
            
            call geo_LUT(this%geo_v, forc_v_from_mass)
            
            nz = ubound(forcing%z,  2)
            ims = lbound(this%geo_u%z,1)
            ime = ubound(this%geo_u%z,1)
            jms = lbound(this%geo_u%z,3)
            jme = ubound(this%geo_u%z,3)
            allocate(forcing%geo_u%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            

            ims = lbound(this%geo_v%z,1)
            ime = ubound(this%geo_v%z,1)
            jms = lbound(this%geo_v%z,3)
            jme = ubound(this%geo_v%z,3)            
            allocate(forcing%geo_v%z(ims:ime, forcing%kts:forcing%kte, jms:jme))
            

            ims = lbound(this%geo%z,1)
            ime = ubound(this%geo%z,1)
            jms = lbound(this%geo%z,3)
            jme = ubound(this%geo%z,3)
            allocate(forcing%geo%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            
            allocate(forcing%geo_agl%z(ims:ime, forcing%kts:forcing%kte, jms:jme))            

            call geo_interp(forcing%geo%z, forcing%z, forcing%geo%geolut)
            call vLUT(this%geo,   forcing%geo)

            call geo_interp(forcing%geo_agl%z, forcing%z, forcing%geo%geolut)
            call geo_interp(forcing%geo_u%z, forcing%z, forc_u_from_mass%geolut)
            call geo_interp(forcing%geo_v%z, forcing%z, forc_v_from_mass%geolut)
            

            if (options%domain%use_agl_height) then
                
                nx = size(this%geo_agl%z, 1)
                ny = size(this%geo_agl%z, 3)
                allocate(AGL_n(nx,ny))
                allocate(AGL_cap(nx,ny))

                nx = size(this%geo_u%z, 1)
                ny = size(this%geo_u%z, 3)
                allocate(AGL_u_n(nx,ny))
                allocate(AGL_u_cap(nx,ny))
                
                nx = size(this%geo_v%z, 1)
                ny = size(this%geo_v%z, 3)
                allocate(AGL_v_n(nx,ny))
                allocate(AGL_v_cap(nx,ny))


                AGL_cap = forcing%geo_agl%z(:,1,:)+real(options%domain%agl_cap)
                where (AGL_cap <= (this%geo_agl%z(:,1,:)+200)) AGL_cap = this%geo_agl%z(:,1,:)+200
                
                AGL_u_cap = forcing%geo_u%z(:,1,:)+real(options%domain%agl_cap)
                where (AGL_u_cap <= (this%geo_u%z(:,1,:)+200)) AGL_u_cap = this%geo_u%z(:,1,:)+200
                
                AGL_v_cap = forcing%geo_v%z(:,1,:)+real(options%domain%agl_cap)
                where (AGL_v_cap <= (this%geo_v%z(:,1,:)+200)) AGL_v_cap = this%geo_v%z(:,1,:)+200

                !Do AGL interpolation for forcing geo z's
                do k=size(forcing%geo_agl%z, 2),1,-1
                    AGL_u_n = (AGL_u_cap-forcing%geo_u%z(:,k,:))/max(abs(AGL_u_cap-forcing%geo_u%z(:,1,:)),0.00001)
                    AGL_v_n = (AGL_v_cap-forcing%geo_v%z(:,k,:))/max(abs(AGL_v_cap-forcing%geo_v%z(:,1,:)),0.00001)
                    AGL_n = (AGL_cap-forcing%geo_agl%z(:,k,:))/max(abs(AGL_cap-forcing%geo_agl%z(:,1,:)),0.00001)

                    where (AGL_n < 0.0) AGL_n = 0.0
                    where (AGL_u_n < 0.0) AGL_u_n = 0.0
                    where (AGL_v_n < 0.0) AGL_v_n = 0.0
                    
                    forcing%geo_u%z(:,k,:) = forcing%geo_u%z(:,k,:)-forcing%geo_u%z(:,1,:)*AGL_u_n
                    forcing%geo_v%z(:,k,:) = forcing%geo_v%z(:,k,:)-forcing%geo_v%z(:,1,:)*AGL_v_n
                    forcing%geo_agl%z(:,k,:) = forcing%geo_agl%z(:,k,:)-forcing%geo_agl%z(:,1,:)*AGL_n
                enddo
                ! Step in reverse so that the bottom level is preserved until it is no longer needed
                ! Do AGL interpolation for domain grid
                
                do k=size(this%geo_agl%z,   2),1,-1
                    ! Multiply subtraction of base-topography by a factor that scales from 1 at surface to 0 at AGL_cap height
                    AGL_u_n = (AGL_u_cap-this%geo_u%z(:,k,:))/max(abs(AGL_u_cap-this%geo_u%z(:,1,:)),0.00001)
                    AGL_v_n = (AGL_v_cap-this%geo_v%z(:,k,:))/max(abs(AGL_v_cap-this%geo_v%z(:,1,:)),0.00001)
                    AGL_n = (AGL_cap-this%geo_agl%z(:,k,:))/max(abs(AGL_cap-this%geo_agl%z(:,1,:)),0.00001)

                    where (AGL_n < 0.0) AGL_n = 0.0
                    where (AGL_u_n < 0.0) AGL_u_n = 0.0
                    where (AGL_v_n < 0.0) AGL_v_n = 0.0
                    
                    this%geo_u%z(:,k,:) = this%geo_u%z(:,k,:)-this%geo_u%z(:,1,:)*AGL_u_n
                    this%geo_v%z(:,k,:) = this%geo_v%z(:,k,:)-this%geo_v%z(:,1,:)*AGL_v_n
                    this%geo_agl%z(:,k,:) = this%geo_agl%z(:,k,:)-this%geo_agl%z(:,1,:)*AGL_n
                enddo
            endif

            call vLUT(this%geo_agl,   forcing%geo_agl)
            call vLUT(this%geo_u, forcing%geo_u)
            call vLUT(this%geo_v, forcing%geo_v)
                        
            ! check that the forcing z is higher than the domain z
            ! if ( maxval(forcing%z) < maxval(this%geo%z) ) then
            !     write(*,*) "ERROR: Forcing or parent-nest z is lower than domain z."
            !     write(*,*) "ERROR: Check earlier output during domain initialization"
            !     write(*,*) "ERROR: to ensure that nested domains fit within their parent domain."
            !     write(*,*) "ERROR: Otherwise, ensure that the vertical extent of all domains"
            !     write(*,*) "ERROR: fits within the forcing domain."
            !     stop
            ! endif

        end if
        

    end subroutine

    subroutine init_relax_filters(this,options)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(options_t),    intent(in)     :: options
        integer :: hs, k, i
        real, dimension(FILTER_WIDTH) :: rs, rs_r
        logical :: corner
        !Setup relaxation filters, start with 2D then expand for 3D version
        
        
        associate( relax_filter => this%vars_2d(this%var_indx(kVARS%relax_filter_2d)%v)%data_2d, relax_filter_3d => this%vars_3d(this%var_indx(kVARS%relax_filter_3d)%v)%data_3d)

        corner = ((this%west_boundary .or. this%east_boundary) .and. (this%north_boundary .or. this%south_boundary))

        hs = this%grid%halo_size

        !relaxation boundary -- set to be 7 for default
        FILTER_WIDTH = min(FILTER_WIDTH,(this%ime-this%ims-hs),(this%jme-this%jms-hs))
        
        if (options%forcing%relax_filters) then
            rs = (/0.9, 0.75, 0.6, 0.5, 0.4, 0.25, 0.1 /)
            rs_r = (/0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9/)
        else
            rs = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
            rs_r = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
        endif

        relax_filter = 0.0
        
        if (this%west_boundary) then
            relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
            do k=this%jms,this%jme
                relax_filter(this%ims+hs:this%ims+hs+FILTER_WIDTH-1,k) = rs(1:FILTER_WIDTH)
            enddo
        endif
        if (this%east_boundary) then
            relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0
            do k=this%jms,this%jme
                relax_filter(this%ime-hs-FILTER_WIDTH+1:this%ime-hs,k) = rs_r(1:FILTER_WIDTH)        
            enddo
        endif
        if (this%north_boundary) then
            relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0
            do k=this%ims,this%ime
                relax_filter(k,this%jme-hs-FILTER_WIDTH+1:this%jme-hs) = rs_r(1:FILTER_WIDTH)
            enddo
        endif
        if (this%south_boundary) then
            relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0
            do k=this%ims,this%ime
                relax_filter(k,this%jms+hs:this%jms+hs+FILTER_WIDTH-1) = rs(1:FILTER_WIDTH)
            enddo
        endif
        if (this%north_boundary .and. this%west_boundary) then
            relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
            relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0

            do i = 1, FILTER_WIDTH
                do k = 1, FILTER_WIDTH
                    relax_filter(this%ims+hs+i-1,this%jme-hs-k+1) = rs(min(i,k))
                enddo
            enddo
        endif
        if (this%north_boundary .and. this%east_boundary) then
            relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0
            relax_filter(this%ims:this%ime,this%jme-hs+1:this%jme) = 1.0

            do i = 1, FILTER_WIDTH
                do k = 1, FILTER_WIDTH
                    relax_filter(this%ime-hs-i+1,this%jme-hs-k+1) = rs(min(i,k))
                enddo
            enddo
        endif
        if (this%south_boundary .and. this%west_boundary) then

            relax_filter(this%ims:this%ims+hs-1,this%jms:this%jme) = 1.0
            relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0

            do i = 1, FILTER_WIDTH
                do k = 1, FILTER_WIDTH
                    relax_filter(this%ims+hs+i-1,this%jms+hs+k-1) = rs(min(i,k))
                enddo
            enddo
        endif
        if (this%south_boundary .and. this%east_boundary) then

            relax_filter(this%ime-hs+1:this%ime,this%jms:this%jme) = 1.0
            relax_filter(this%ims:this%ime,this%jms:this%jms+hs-1) = 1.0

            do i = 1, FILTER_WIDTH
                do k = 1, FILTER_WIDTH
                    relax_filter(this%ime-hs-i+1,this%jms+hs+k-1) = rs(min(i,k))
                enddo
            enddo
        endif

        do k=this%kms,this%kme
            relax_filter_3d(this%ims:this%ime,k,this%jms:this%jme) = relax_filter
        enddo
        
        end associate

    end subroutine init_relax_filters
    
    
    !> -------------------------------
    !! Update the dQdt fields for all forced variables which force the whole domain
    !! Forced variables which force just the boundary are handeled by a similar function called on the boundary object
    !! 
    !! For domain-forced variables, this routine is the partner of apply_forcing below.
    !! update_delta_fields normalizes the difference by the time step of that difference field
    !! apply_forcing multiplies that /second value and multiplies it by the current time step before adding it
    !!
    !! -------------------------------
    module subroutine update_delta_fields(this)
        implicit none
        class(domain_t),    intent(inout) :: this

        type(time_delta_t)  :: dt
        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_update
        integer :: i, var_indx

        dt = this%next_input - this%sim_time

        ! check if the difference between the simulation time and next_input is less than an input_dt
        ! if so, this signals that we are in between two input times, so advance the state of the variables%data_3d
        ! to the simulation time

        !include "-1" to accomodate rounding errors
        if (dt%seconds() < (this%input_dt%seconds()-1)) call dt%set(seconds= (this%input_dt%seconds() - dt%seconds()) )

        ! Now iterate through the dictionary as long as there are more elements present
        do i = 1,size(this%forcing_hi)
            !Update delta fields on the high-resolution forcing varaibles...

            if (this%forcing_hi(i)%two_d) then
                this%forcing_hi(i)%dqdt_2d = (this%forcing_hi(i)%dqdt_2d - this%forcing_hi(i)%data_2d) / dt%seconds()
            else if (this%forcing_hi(i)%three_d) then
                this%forcing_hi(i)%dqdt_3d = (this%forcing_hi(i)%dqdt_3d - this%forcing_hi(i)%data_3d) / dt%seconds()
            endif

            ! now update delta fields for domain variables
            var_indx = get_varindx(trim(this%forcing_hi(i)%name))
            var_to_update = get_varmeta(var_indx)

            if (var_to_update%force_boundaries) cycle

            if (var_to_update%two_d) then
                this%vars_2d(this%var_indx(var_indx)%v)%dqdt_2d = (this%vars_2d(this%var_indx(var_indx)%v)%dqdt_2d - this%vars_2d(this%var_indx(var_indx)%v)%data_2d) / dt%seconds()
            else if (var_to_update%three_d) then
                this%vars_3d(this%var_indx(var_indx)%v)%dqdt_3d = (this%vars_3d(this%var_indx(var_indx)%v)%dqdt_3d - this%vars_3d(this%var_indx(var_indx)%v)%data_3d) / dt%seconds()
            endif
            
        enddo


        ! w has to be handled separately because it is the only variable that can be updated using the delta fields but is not
        ! actually read from disk. Note that if we move to balancing winds every timestep, then it doesn't matter.
        ! this%vars_3d(this%var_indx(kVARS%w)%v)%dqdt_3d = (this%vars_3d(this%var_indx(kVARS%w)%v)%dqdt_3d - this%vars_3d(this%var_indx(kVARS%w)%v)%data_3d) / dt%seconds()

    end subroutine


    !> -------------------------------
    !! Add the forcing update to boundaries and internal diagnosed fields
    !!
    !! This routine is the partner of update_delta_fields above.
    !! update_delta_fields normalizes the difference by the time step of that difference field
    !! apply forcing multiplies that /second value and multiplies it by the current time step before adding it
    !!
    !! -------------------------------
    module subroutine apply_forcing(this, options, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(options_t), intent(in)       :: options
        real, intent(in)                  :: dt
        integer :: ims, ime, jms, jme
        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_update
        integer :: i, k, j, var_indx, n
        real    :: dt_h
        logical :: do_boundary, is_wind, is_w_real
        
        !calculate dt in units of hours
        dt_h = dt/3600.0

        do n = 1,size(this%forcing_hi)

            var_indx = get_varindx(trim(this%forcing_hi(n)%name))
            var_to_update = get_varmeta(var_indx)            
            is_w_real = (this%var_indx(var_indx)%v == this%var_indx(kVARS%w_real)%v)
            is_wind = (this%var_indx(var_indx)%v == this%var_indx(kVARS%u)%v) .or. (this%var_indx(var_indx)%v == this%var_indx(kVARS%v)%v) .or. is_w_real

            if (var_to_update%two_d) then
                ims = this%vars_2d(this%var_indx(var_indx)%v)%grid%ims
                ime = this%vars_2d(this%var_indx(var_indx)%v)%grid%ime
                jms = this%vars_2d(this%var_indx(var_indx)%v)%grid%jms
                jme = this%vars_2d(this%var_indx(var_indx)%v)%grid%jme
    
                do_boundary = (ims < this%ids+this%grid%halo_size+FILTER_WIDTH) .or. (ime > this%ide-this%grid%halo_size-FILTER_WIDTH) .or. &
                    (jms < this%jds+this%grid%halo_size+FILTER_WIDTH) .or. (jme > this%jde-this%grid%halo_size-FILTER_WIDTH)

                ! apply forcing throughout the domain for 2D diagnosed variables (e.g. SST, SW)
                if (.not.(var_to_update%force_boundaries)) then
                    do concurrent (j = jms:jme, i = ims:ime)
                        this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) = this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) + (this%vars_2d(this%var_indx(var_indx)%v)%dqdt_2d(i,j) * dt)
                        this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) = max(this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j),0.0)
                    enddo
                else if (do_boundary) then
                    if (any(this%vars_2d(this%var_indx(kVARS%relax_filter_2d)%v)%data_2d > 0.0)) then
                        !Update forcing data to current time step
                        do concurrent (j = jms:jme, i = ims:ime)
                            if (this%vars_2d(this%var_indx(kVARS%relax_filter_2d)%v)%data_2d(i,j) > 0.0) then
                                this%forcing_hi(n)%data_2d(i,j) = this%forcing_hi(n)%data_2d(i,j) + (this%forcing_hi(n)%dqdt_2d(i,j) * dt)
                                this%forcing_hi(n)%data_2d(i,j) = max(this%forcing_hi(n)%data_2d(i,j),0.0)

                                if (this%vars_2d(this%var_indx(kVARS%relax_filter_2d)%v)%data_2d(i,j) == 1.0) then
                                    this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) = this%forcing_hi(n)%data_2d(i,j)
                                else
                                    this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) = this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) + &
                                                    (this%vars_2d(this%var_indx(kVARS%relax_filter_2d)%v)%data_2d(i,j) * dt_h) * &
                                                    (this%forcing_hi(n)%data_2d(i,j) - this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j))

                                    this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j) = max(this%vars_2d(this%var_indx(var_indx)%v)%data_2d(i,j),0.0)
                                endif
                            endif
                        enddo
                    endif
                endif 

            else if (var_to_update%three_d) then
                ims = this%vars_3d(this%var_indx(var_indx)%v)%grid%ims
                ime = this%vars_3d(this%var_indx(var_indx)%v)%grid%ime
                jms = this%vars_3d(this%var_indx(var_indx)%v)%grid%jms
                jme = this%vars_3d(this%var_indx(var_indx)%v)%grid%jme

                !see if we are on the boundary of the domain
                do_boundary = (ims < this%ids+this%grid%halo_size+FILTER_WIDTH) .or. (ime > this%ide-this%grid%halo_size-FILTER_WIDTH) .or. &
                               (jms < this%jds+this%grid%halo_size+FILTER_WIDTH) .or. (jme > this%jde-this%grid%halo_size-FILTER_WIDTH)

                ! only apply forcing data on the boundaries for advected scalars (e.g. temperature, humidity)
                ! applying forcing to the edges has already been handeled when updating dqdt using the relaxation filter
                if (.not.(var_to_update%force_boundaries)) then
                    do concurrent (j = jms:jme, k = this%kms:this%kme, i = ims:ime)
                        this%forcing_hi(n)%data_3d(i,k,j)    = this%forcing_hi(n)%data_3d(i,k,j) + (this%forcing_hi(n)%dqdt_3d(i,k,j) * dt)
                        if (.not.(is_wind)) this%forcing_hi(n)%data_3d(i,k,j) = max(this%forcing_hi(n)%data_3d(i,k,j),0.0)

                        if (.not.(is_w_real)) this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) = this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) + &
                                                        (this%vars_3d(this%var_indx(var_indx)%v)%dqdt_3d(i,k,j) * dt)
                                                        
                        if (.not.(is_wind)) this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) = max(this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j),0.0)
                    enddo
                else if (do_boundary) then
                    if (any(this%vars_3d(this%var_indx(kVARS%relax_filter_3d)%v)%data_3d > 0.0)) then
                        !Update forcing data to current time step
                        do concurrent (j = jms:jme, k = this%kms:this%kme, i = ims:ime)
                            if (this%vars_3d(this%var_indx(kVARS%relax_filter_3d)%v)%data_3d(i,k,j) > 0.0) then
                                this%forcing_hi(n)%data_3d(i,k,j) = this%forcing_hi(n)%data_3d(i,k,j) + (this%forcing_hi(n)%dqdt_3d(i,k,j) * dt)
                                this%forcing_hi(n)%data_3d(i,k,j) = max(this%forcing_hi(n)%data_3d(i,k,j),0.0)

                                if (this%vars_3d(this%var_indx(kVARS%relax_filter_3d)%v)%data_3d(i,k,j) == 1.0) then
                                    this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) = this%forcing_hi(n)%data_3d(i,k,j)
                                else
                                    this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) = this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) + &
                                                    (this%vars_3d(this%var_indx(kVARS%relax_filter_3d)%v)%data_3d(i,k,j) * dt_h) * &
                                                    (this%forcing_hi(n)%data_3d(i,k,j) - this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j))

                                    this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j) = max(this%vars_3d(this%var_indx(var_indx)%v)%data_3d(i,k,j),0.0)
                                endif
                            endif
                        enddo
                    endif
                endif
            endif
        enddo

        ! w has to be handled separately because it is the only variable that can be updated using the delta fields but is not
        ! actually read from disk. Note that if we move to balancing winds every timestep, then it doesn't matter.
        ! if (.not.(options%adv%advect_density)) then
        !     do concurrent (j = jms:jme, k = this%kms:this%kme, i = ims:ime)
        !         this%vars_3d(this%var_indx(kVARS%w)%v)%data_3d(i,k,j) = this%vars_3d(this%var_indx(kVARS%w)%v)%data_3d(i,k,j) + (this%vars_3d(this%var_indx(kVARS%w)%v)%dqdt_3d(i,k,j) * dt)
        !     enddo
        ! endif


    end subroutine


    !> -------------------------------
    !! Loop through all variables for which forcing data have been supplied and interpolate the forcing data to the domain
    !!
    !! -------------------------------
    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        logical,          intent(in),   optional :: update

        ! internal field always present for value of optional "update"
        logical :: update_only
        ! temporary to hold the variable to be interpolated to
        type(variable_t) :: var_to_interpolate
        ! temporary to hold the forcing variable to be interpolated from
        type(variable_t) :: input_data

        ! number of layers has to be used when subsetting for update_pressure (for now)
        integer :: nz, i, var_indx, pressure_indx, pot_temp_indx
        logical :: var_is_u, var_is_v, var_is_pressure, var_is_potential_temp, agl_interp

        update_only = .False.
        if (present(update)) update_only = update


        ! Now iterate through the dictionary as long as there are more elements present
        do i = 1,size(this%forcing_hi)

            var_indx = get_varindx(trim(this%forcing_hi(i)%name))
            var_to_interpolate = get_varmeta(var_indx)

            if (var_to_interpolate%two_d) then
                var_to_interpolate = this%vars_2d(this%var_indx(var_indx)%v)
            else if (var_to_interpolate%three_d) then
                var_to_interpolate = this%vars_3d(this%var_indx(var_indx)%v)
            endif

            ! get the associated forcing data
            input_data = forcing%variables%get_var(var_to_interpolate%forcing_var)
            ! interpolate
            if (var_to_interpolate%two_d) then
                if (update_only) then
                    call geo_interp2d(this%forcing_hi(i)%dqdt_2d, input_data%data_2d, forcing%geo%geolut)
                    !If this variable is forcing the whole domain, we can copy the next forcing step directly over to domain
                    if (.not.(var_to_interpolate%force_boundaries)) this%vars_2d(this%var_indx(var_indx)%v)%dqdt_2d = this%forcing_hi(i)%dqdt_2d
                else
                    call geo_interp2d(this%forcing_hi(i)%data_2d, input_data%data_2d, forcing%geo%geolut)
                    !If this is an initialization step, copy high res directly over to domain
                    this%vars_2d(this%var_indx(var_indx)%v)%data_2d = this%forcing_hi(i)%data_2d
                endif
            else

                var_is_pressure = (trim(var_to_interpolate%forcing_var) == trim(this%vars_3d(this%var_indx(kVARS%pressure)%v)%forcing_var))
                var_is_potential_temp = (trim(var_to_interpolate%forcing_var) == trim(this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%forcing_var))
                var_is_u = (trim(var_to_interpolate%forcing_var) == trim(this%vars_3d(this%var_indx(kVARS%u)%v)%forcing_var))
                var_is_v = (trim(var_to_interpolate%forcing_var) == trim(this%vars_3d(this%var_indx(kVARS%v)%v)%forcing_var))
                !If we are dealing with anything but pressure and temperature (basically mass/number species), consider height above ground
                !for interpolation. If the user has not selected AGL interpolation in the namelist, this will result in standard z-interpolation
                agl_interp = .not.(var_is_pressure .or. var_is_potential_temp)

                ! if just updating, use the dqdt variable otherwise use the 3D variable
                if (update_only) then
                    call interpolate_variable(this%forcing_hi(i)%dqdt_3d, input_data, forcing, this, &
                                    interpolate_agl_in=agl_interp, var_is_u=var_is_u, var_is_v=var_is_v, nsmooth=this%nsmooth)
                    !If this variable is forcing the whole domain, we can copy the next forcing step directly over to domain
                    if (.not.(var_to_interpolate%force_boundaries).and..not.var_is_u.and..not.var_is_v) this%vars_3d(this%var_indx(var_indx)%v)%dqdt_3d = this%forcing_hi(i)%dqdt_3d
                else
                    call interpolate_variable(this%forcing_hi(i)%data_3d, input_data, forcing, this, &
                                    interpolate_agl_in=agl_interp, var_is_u=var_is_u, var_is_v=var_is_v, nsmooth=this%nsmooth)
                    !If this is an initialization step, copy high res directly over to domain
                    this%vars_3d(this%var_indx(var_indx)%v)%data_3d = this%forcing_hi(i)%data_3d
                endif
                if (var_is_pressure) pressure_indx = i
                if (var_is_potential_temp) pot_temp_indx = i
            endif
            call forcing%variables%add_var(var_to_interpolate%forcing_var, input_data)
        enddo

        !Adjust potential temperature (first) and pressure (second) to account for points below forcing grid
        !Only domain-wide-forced variables are updated with the domain dqdt_3d
        if (update_only) then
            call adjust_pressure_temp(this%forcing_hi(pressure_indx)%dqdt_3d,this%forcing_hi(pot_temp_indx)%dqdt_3d, forcing%geo%z, this%geo%z)
            this%vars_3d(this%var_indx(kVARS%pressure)%v)%dqdt_3d = this%forcing_hi(pressure_indx)%dqdt_3d
            this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%dqdt_3d = this%forcing_hi(pot_temp_indx)%dqdt_3d
        else
            call adjust_pressure_temp(this%forcing_hi(pressure_indx)%data_3d,this%forcing_hi(pot_temp_indx)%data_3d, forcing%geo%z, this%geo%z)
            this%vars_3d(this%var_indx(kVARS%pressure)%v)%data_3d = this%forcing_hi(pressure_indx)%data_3d
            this%vars_3d(this%var_indx(kVARS%potential_temperature)%v)%data_3d = this%forcing_hi(pot_temp_indx)%data_3d
        endif

        !Ensure that input data for hydrometeors after interpolation have been forced to 0-minimum
        call this%enforce_limits(update_in=update_only)

        !Perform a diagnostic_update to ensure that all diagnostic variables are set for the new forcing data
        !This will be overwriten as soon as we enter the physics loop, but it is necesery to compute density
        !For the future step so that the wind solver uses both future winds, and future density.
        call this%diagnostic_update(forcing_update=update_only)

    end subroutine

    subroutine adjust_pressure_temp(pressure, potential_temp, input_z, output_z)
        implicit none
        real, intent(inout), dimension(:,:,:) :: pressure, potential_temp
        real, intent(in), dimension(:,:,:) :: input_z, output_z !> z on the forcing and ICAR model levels [m]
        integer :: i,j,k, nz, nx, ny
        real    :: t, p_guess, dz, H
        
        !For all output_z less than input_z, extrapolate downwards based on lapse rate of -6.5C/km
        
        nx = size(potential_temp, 1)
        nz = size(potential_temp, 2)
        ny = size(potential_temp, 3)

        do j = 1, ny
            do i = 1, nx
                do k = 1, nz
                    if (input_z(i,1,j) > output_z(i,k,j)) then
                        
                        !From vertical interpolation, potential_temperature and pressure will be kept constant when below the grid
                        !So the current values at these below-indices reflect the temp/pressure of the closest forcing grid cell
                    
                        dz = input_z(i,1,j)-output_z(i,k,j)
                        
                        !Assume lapse rate of -6.5ºC/1km
                        potential_temp(i,k,j) = potential_temp(i,k,j) + 6.5*dz/1000.0
                        
                        !estimate pressure difference 1100 Pa for each 100m difference for exner function
                        H = 29.3 * (potential_temp(i,k,j) * exner_function(pressure(i,k,j)))
                        pressure(i,k,j) = pressure(i,k,j) * exp(dz/H)
                        !pressure(i,k,j) = pressure(i,k,j) * exp( ((gravity/R_d) * dz) / t )
                    else
                        exit
                    endif
                end do
            enddo
        enddo


    end subroutine adjust_pressure_temp

    !> -------------------------------
    !! Adjust a 3d pressure field from the forcing data to the ICAR model grid
    !!
    !! Because the GCM grid can be very different from the ICAR grid, we first roughly match up
    !! the GCM level that is closest to the ICAR level. This has to be done grid cell by gridcell.
    !! This still is not ideal, in that it has already subset the GCM levels to the same number as are in ICAR
    !! If the GCM has a LOT of fine layers ICAR will not be getting layers higher up in the atmosphere.
    !! It would be nice to first use vinterp to get as close as we can, then update pressure only for grid cells below.
    !! Uses update_pressure to make a final adjustment (including below the lowest model level).
    !!
    !! -------------------------------
    subroutine adjust_pressure(pressure, input_z, output_z, potential_temperature)
        implicit none
        real, intent(inout), dimension(:,:,:) :: pressure !> Pressure on the forcing model levels [Pa]
        real, intent(in), dimension(:,:,:) :: input_z, output_z !> z on the forcing and ICAR model levels [m]
        real, intent(in), dimension(:,:,:) :: potential_temperature !> potential temperature of the forcing data [K]

        ! store a temporary copy of P and Z from the forcing data after selecting the closest GCM level to the ICAR data
        real, allocatable, dimension(:,:,:) :: temp_z, temp_p, temp_t
        ! loop counter variables
        integer :: k, nz, in_z_idx
        integer :: i,j, nx, ny

        allocate(temp_z, temp_p, temp_t, mold=pressure)

        nx = size(pressure, 1)
        nz = size(pressure, 2)
        ny = size(pressure, 3)

        do j = 1, ny
            do i = 1, nx
                ! keep track of the nearest z level from the forcing data
                in_z_idx = 1
                do k = 1, nz
                    ! if the ICAR z level is more than half way to the next forcing z level, then increment the GCM z
                    findz: do while (output_z(i,k,j) > ((input_z(i,in_z_idx,j) + input_z(i,min(nz,in_z_idx+1),j)) / 2))
                        in_z_idx = min(nz, in_z_idx + 1)

                        if (in_z_idx == nz) then
                            exit findz
                        endif
                    end do findz
                    ! make a new copy of the pressure and z data from the closest GCM model level
                    temp_z(i,k,j) = input_z(i,in_z_idx,j)
                    temp_p(i,k,j) = pressure(i,in_z_idx,j)
                    temp_t(i,k,j) = exner_function(pressure(i,in_z_idx,j)) * potential_temperature(i,in_z_idx,j)
                end do
            enddo
        enddo

        ! put the updated pressure data into the pressure variable prior to adjustments
        pressure = temp_p


        ! update pressure for the change in height between the closest GCM model level and each ICAR level.
        call update_pressure(pressure, temp_z, output_z, temp_t)

    end subroutine

    !> -------------------------------
    !! Interpolate one variable by requesting the forcing data from the boundary data structure then
    !! calling the appropriate interpolation routine (2D vs 3D) with the appropriate grid (mass, u, v)
    !!
    !! -------------------------------
    subroutine interpolate_variable(var_data, input_data, forcing, dom, interpolate_agl_in, var_is_u, var_is_v, nsmooth)
        implicit none
        real,            intent(inout) :: var_data(:,:,:)
        type(variable_t),   intent(inout) :: input_data
        type(boundary_t),   intent(in)    :: forcing
        type(domain_t),     intent(in)    :: dom
        logical,            intent(in),   optional :: interpolate_agl_in
        logical,            intent(in),   optional :: var_is_u, var_is_v
        integer,            intent(in),   optional :: nsmooth

        ! note that 3D variables have a different number of vertical levels, so they have to first be interpolated
        ! to the high res horizontal grid, then vertically interpolated to the actual icar domain
        real, allocatable :: temp_3d(:,:,:)
        logical :: interpolate_agl, uvar, vvar
        integer :: nx, ny, nz, ims, ime, jms, jme
        integer :: windowsize, z

        interpolate_agl=.False.
        if (present(interpolate_agl_in)) interpolate_agl = interpolate_agl_in
        uvar = .False.
        if (present(var_is_u)) uvar = var_is_u
        vvar = .False.
        if (present(var_is_v)) vvar = var_is_v
        windowsize = 0
        if (present(nsmooth)) windowsize = nsmooth

        ims = lbound(var_data,1)
        ime = ubound(var_data,1)
        jms = lbound(var_data,3)
        jme = ubound(var_data,3)

        ! allocate a temporary variable to hold the horizontally interpolated data before vertical interpolation
        allocate(temp_3d(ims:ime, size(input_data%data_3d,2), jms:jme ))

        ! Sequence of if statements to test if this variable needs to be interpolated onto the staggared grids
        ! This could all be combined by passing in the geo data to use, along with a smoothing flag.

        ! Interpolate to the Mass grid
        if ((size(var_data,1) == size(forcing%geo%geolut%x,2)).and.(size(var_data,3) == size(forcing%geo%geolut%x,3))) then

            call geo_interp(temp_3d, input_data%data_3d, forcing%geo%geolut)

            if (interpolate_agl) then
                call vinterp(var_data, temp_3d, forcing%geo_agl%vert_lut)
            else
                call vinterp(var_data, temp_3d, forcing%geo%vert_lut)
            endif
            
        ! Interpolate to the u staggered grid
        else if (uvar) then

            ! One grid cell smoothing of original input data
            if (windowsize > 0) call smooth_array(input_data%data_3d, windowsize=1, ydim=3)
            call geo_interp(temp_3d, input_data%data_3d, forcing%geo_u%geolut)

            call vinterp(var_data, temp_3d, forcing%geo_u%vert_lut)
            ! temp_3d = pre_smooth(:,:nz,:) ! no vertical interpolation option
            if (windowsize > 0) call smooth_array(var_data, windowsize=windowsize, ydim=3)
                        
        ! Interpolate to the v staggered grid
        else if (vvar) then

            ! One grid cell smoothing of original input data
            if (windowsize > 0) call smooth_array(input_data%data_3d, windowsize=1, ydim=3)
            call geo_interp(temp_3d, input_data%data_3d, forcing%geo_v%geolut)
            
            call vinterp(var_data, temp_3d, forcing%geo_v%vert_lut)
            ! temp_3d = pre_smooth(:,:nz,:) ! no vertical interpolation option
            if (windowsize > 0) call smooth_array(var_data, windowsize=windowsize, ydim=3)
        endif
        
    end subroutine



end submodule
