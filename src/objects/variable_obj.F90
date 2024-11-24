submodule(variable_interface) variable_implementation
    implicit none


contains

    !> -------------------------------
    !! Initialize a variable object from a given grid
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_grid(this, grid, forcing_var, force_boundaries, dtype)
        implicit none
        class(variable_t),  intent(inout) :: this
        type(grid_t),       intent(in)    :: grid
        character(len=*),   intent(in), optional :: forcing_var
        logical,            intent(in), optional :: force_boundaries
        integer,            intent(in), optional :: dtype

        integer :: err

        this%dtype = kREAL
        if (present(dtype)) this%dtype = dtype

        this%two_d   = grid%is2d
        this%three_d = grid%is3d

        this%dimensions = grid%dimensions
        this%dim_len    = grid%get_dims()
        this%grid       = grid
        
        if (allocated(this%global_dim_len)) deallocate(this%global_dim_len)
        if (this%two_d) then
            allocate(this%global_dim_len(2))
            this%global_dim_len(1) = grid%ide
            this%global_dim_len(2) = grid%jde
        else
            allocate(this%global_dim_len(3))
            this%global_dim_len(1) = grid%ide
            this%global_dim_len(2) = grid%kde
            this%global_dim_len(3) = grid%jde
        endif
        this%forcing_var = ""
        if (present(forcing_var)) this%forcing_var = forcing_var

        this%force_boundaries = .True.
        if (present(force_boundaries)) this%force_boundaries = force_boundaries

        if (grid%is2d) then
            this%n_dimensions = 2
            if (associated(this%data_2d)) deallocate(this%data_2d)
            if (this%dtype == kREAL) then
                allocate(this%data_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2d = 0
            elseif (this%dtype == kDOUBLE) then
                allocate(this%data_2dd(grid%ims:grid%ime,    &
                                       grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2dd = 0
            elseif (this%dtype == kINTEGER) then
                allocate(this%data_2di(grid%ims:grid%ime,    &
                                       grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2di = 0
        
            endif

            if (trim(this%forcing_var) /= "") then
                allocate(this%dqdt_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_2d: Allocation request failed"

                this%dqdt_2d = 0
            endif

        endif

        if (grid%is3d) then
            this%n_dimensions = 3
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(grid%ims:grid%ime,    &
                                  grid%kms:grid%kme,    &
                                  grid%jms:grid%jme), stat=err)
            if (err /= 0) stop "variable:grid:3d: Allocation request failed"

            this%data_3d = 0

            if (trim(this%forcing_var) /= "") then
                allocate(this%dqdt_3d(grid%ims:grid%ime,    &
                                      grid%kms:grid%kme,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_3d: Allocation request failed"

                this%dqdt_3d = 0
            endif

        endif
        this%xstag = grid%nx_e
        this%ystag = grid%ny_e

    end subroutine

    !> -------------------------------
    !! Initialize a variable object from a given array of dimension sizes
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_dims(this, dims, forcing_var, force_boundaries)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: dims(:)
        character(len=*),   intent(in), optional :: forcing_var
        logical,            intent(in), optional :: force_boundaries

        integer :: err

        this%dim_len    = dims

        this%two_d   = size(dims) == 2
        this%three_d = size(dims) == 3

        this%forcing_var = ""
        if (present(forcing_var)) this%forcing_var = forcing_var

        this%force_boundaries = .True.
        if (present(force_boundaries)) this%force_boundaries = force_boundaries

        if (this%two_d) then
            this%n_dimensions = 2
            this%dimensions = ['x','y']
            if (associated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(dims(1), dims(2)), stat=err)
            if (err /= 0) stop "variable:dims:2d: Allocation request denied"
            this%data_2d = 0

            if (trim(this%forcing_var) /= "") then
                if (associated(this%dqdt_2d)) deallocate(this%dqdt_2d)
                allocate(this%dqdt_2d(dims(1), dims(2)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_2d: Allocation request denied"

                this%dqdt_2d = 0
            endif
        endif

        if (this%three_d) then
            this%n_dimensions = 3
            this%dimensions = ['x','y','z']
            if (associated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(dims(1), dims(2), dims(3)), stat=err)
            if (err /= 0) stop "variable:dims:3d: Allocation request denied"

            this%data_3d = 0

            if (trim(this%forcing_var) /= "") then
                if (associated(this%dqdt_3d)) deallocate(this%dqdt_3d)
                allocate(this%dqdt_3d(dims(1), dims(2), dims(3)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_3d: Allocation request denied"

                this%dqdt_3d = 0
            endif
        endif

    end subroutine


end submodule
