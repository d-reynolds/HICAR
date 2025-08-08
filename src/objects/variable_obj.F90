submodule(variable_interface) variable_implementation
    use output_metadata, only : get_varmeta
    implicit none

contains

    !> -------------------------------
    !! Initialize a variable object from a given grid
    !!
    !! Allocates 2d/3d data structure as appropriate
    !!
    !! -------------------------------
    module subroutine init_grid(this, var_idx, grid, forcing_var)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: var_idx
        type(grid_t),       intent(in)    :: grid
        logical,            intent(in), optional :: forcing_var

        integer :: err


        this%id = var_idx
        call this%set_from_metadata(var_idx)

        this%one_d   = grid%is1d
        this%two_d   = grid%is2d
        this%three_d = grid%is3d
        this%four_d  = grid%is4d
        this%dim_len    = grid%get_dims()
        this%grid       = grid
        if (allocated(this%global_dim_len)) deallocate(this%global_dim_len)
        if (this%one_d) then
            allocate(this%global_dim_len(1))
            this%global_dim_len(1) = grid%kde
        else if (this%two_d) then
            allocate(this%global_dim_len(2))
            this%global_dim_len(1) = grid%ide
            this%global_dim_len(2) = grid%jde
        else if (this%three_d) then
            allocate(this%global_dim_len(3))
            this%global_dim_len(1) = grid%ide
            this%global_dim_len(2) = grid%kde
            this%global_dim_len(3) = grid%jde
        else if (this%four_d) then
            allocate(this%global_dim_len(4))
            this%global_dim_len(1) = grid%ide
            this%global_dim_len(2) = grid%kde
            this%global_dim_len(3) = grid%jde
            this%global_dim_len(4) = grid%n_4d
        endif

        if (present(forcing_var)) this%forcing_var = forcing_var

        if (grid%is1d) then
            if (allocated(this%data_1d)) deallocate(this%data_1d)
            allocate(this%data_1d(grid%kms:grid%kme), stat=err)
            if (err /= 0) stop "variable:grid:1d: Allocation request failed"

            this%data_1d = 0

        endif

        if (grid%is2d) then

            if (allocated(this%data_2d)) deallocate(this%data_2d)
            if (allocated(this%data_2di)) deallocate(this%data_2di)
            if (allocated(this%dqdt_2d)) deallocate(this%dqdt_2d)

            if (this%dtype == kREAL) then
                allocate(this%data_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2d = 0
            ! elseif (this%dtype == kDOUBLE) then
            !     allocate(this%data_2dd(grid%ims:grid%ime,    &
            !                            grid%jms:grid%jme), stat=err)
            !     if (err /= 0) stop "variable:grid:2d: Allocation request failed"

            !     this%data_2dd = 0
            elseif (this%dtype == kINTEGER) then
                allocate(this%data_2di(grid%ims:grid%ime,    &
                                       grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:2d: Allocation request failed"

                this%data_2di = 0
        
            endif

            if (this%forcing_var) then
                allocate(this%dqdt_2d(grid%ims:grid%ime,    &
                                      grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_2d: Allocation request failed"

                this%dqdt_2d = 0
            endif

        endif

        if (grid%is3d) then
            if (allocated(this%data_3d)) deallocate(this%data_3d)
            if (allocated(this%dqdt_3d)) deallocate(this%dqdt_3d)
            allocate(this%data_3d(grid%ims:grid%ime,    &
                                  grid%kms:grid%kme,    &
                                  grid%jms:grid%jme), stat=err)
            if (err /= 0) stop "variable:grid:3d: Allocation request failed"

            this%data_3d = 0
        
            ! note w is special cased because it does not have a forcing variable, so it is not necessarily allocated automatically
            if (this%forcing_var .or. this%id==kVARS%w) then
                allocate(this%dqdt_3d(grid%ims:grid%ime,    &
                                        grid%kms:grid%kme,    &
                                        grid%jms:grid%jme), stat=err)
                if (err /= 0) stop "variable:grid:dqdt_3d: Allocation request failed"

                this%dqdt_3d = 0
            endif
        endif
        if (grid%is4d) then
            if (allocated(this%data_4d)) deallocate(this%data_4d)
            allocate(this%data_4d(grid%ims:grid%ime,    &
                                    grid%kms:grid%kme,    &
                                    grid%jms:grid%jme,    &
                                    1:grid%n_4d), stat=err)
            if (err /= 0) stop "variable:grid:4d: Allocation request failed"

            this%data_4d = 0
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
    module subroutine init_dims(this, var_idx, dims, forcing_var)
        implicit none
        class(variable_t),  intent(inout) :: this
        integer,            intent(in)    :: var_idx
        integer,            intent(in)    :: dims(:)
        logical,            intent(in), optional :: forcing_var

        integer :: err


        call this%set_from_metadata(var_idx)

        this%dim_len    = dims

        this%two_d   = size(dims) == 2
        this%three_d = size(dims) == 3

        this%forcing_var = .False.
        if (present(forcing_var)) this%forcing_var = forcing_var

        if (this%two_d) then
            if (allocated(this%data_2d)) deallocate(this%data_2d)
            allocate(this%data_2d(dims(1), dims(2)), stat=err)
            if (err /= 0) stop "variable:dims:2d: Allocation request denied"
            this%data_2d = 0

            if (this%forcing_var) then
                if (allocated(this%dqdt_2d)) deallocate(this%dqdt_2d)
                allocate(this%dqdt_2d(dims(1), dims(2)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_2d: Allocation request denied"

                this%dqdt_2d = 0
            endif
        endif

        if (this%three_d) then
            if (allocated(this%data_3d)) deallocate(this%data_3d)
            allocate(this%data_3d(dims(1), dims(2), dims(3)), stat=err)
            if (err /= 0) stop "variable:dims:3d: Allocation request denied"

            this%data_3d = 0

            if (this%forcing_var) then
                if (allocated(this%dqdt_3d)) deallocate(this%dqdt_3d)
                allocate(this%dqdt_3d(dims(1), dims(2), dims(3)), stat=err)
                if (err /= 0) stop "variable:dims:dqdt_3d: Allocation request denied"

                this%dqdt_3d = 0
            endif
        endif

    end subroutine

    module subroutine assign_variable(dest, src)
        implicit none
        class(variable_t), intent(out) :: dest
        class(variable_t), intent(in)  :: src
                
        ! Copy all allocatable arrays
        if (allocated(src%data_4d)) then
            if (allocated(dest%data_4d)) deallocate(dest%data_4d)
            allocate(dest%data_4d, source=src%data_4d)
        else if (allocated(dest%data_4d)) then
            deallocate(dest%data_4d)
        endif
        
        if (allocated(src%data_3d)) then
            if (allocated(dest%data_3d)) deallocate(dest%data_3d)
            allocate(dest%data_3d, source=src%data_3d)
        else if (allocated(dest%data_3d)) then
            deallocate(dest%data_3d)
        endif
        
        if (allocated(src%data_2d)) then
            if (allocated(dest%data_2d)) deallocate(dest%data_2d)
            allocate(dest%data_2d, source=src%data_2d)
        else if (allocated(dest%data_2d)) then
            deallocate(dest%data_2d)
        endif
        
        if (allocated(src%data_1d)) then
            if (allocated(dest%data_1d)) deallocate(dest%data_1d)
            allocate(dest%data_1d, source=src%data_1d)
        else if (allocated(dest%data_1d)) then
            deallocate(dest%data_1d)
        endif
        
        if (allocated(src%data_2di)) then
            if (allocated(dest%data_2di)) deallocate(dest%data_2di)
            allocate(dest%data_2di, source=src%data_2di)
        else if (allocated(dest%data_2di)) then
            deallocate(dest%data_2di)
        endif
        
        if (allocated(src%dqdt_3d)) then
            if (allocated(dest%dqdt_3d)) deallocate(dest%dqdt_3d)
            allocate(dest%dqdt_3d, source=src%dqdt_3d)
        else if (allocated(dest%dqdt_3d)) then
            deallocate(dest%dqdt_3d)
        endif
        
        if (allocated(src%dqdt_2d)) then
            if (allocated(dest%dqdt_2d)) deallocate(dest%dqdt_2d)
            allocate(dest%dqdt_2d, source=src%dqdt_2d)
        else if (allocated(dest%dqdt_2d)) then
            deallocate(dest%dqdt_2d)
        endif
        
        ! Copy all scalar variables
        dest%one_d = src%one_d
        dest%three_d = src%three_d
        dest%two_d = src%two_d
        dest%four_d = src%four_d
        dest%force_boundaries = src%force_boundaries
        dest%computed = src%computed
        dest%forcing_var = src%forcing_var
        
        dest%dtype = src%dtype
        dest%grid = src%grid
        dest%id = src%id
        
        ! Copy allocatable arrays
        if (allocated(src%dim_len)) then
            if (allocated(dest%dim_len)) deallocate(dest%dim_len)
            allocate(dest%dim_len, source=src%dim_len)
        else if (allocated(dest%dim_len)) then
            deallocate(dest%dim_len)
        endif
        
        if (allocated(src%global_dim_len)) then
            if (allocated(dest%global_dim_len)) deallocate(dest%global_dim_len)
            allocate(dest%global_dim_len, source=src%global_dim_len)
        else if (allocated(dest%global_dim_len)) then
            deallocate(dest%global_dim_len)
        endif
                
        if (allocated(src%dim_ids)) then
            if (allocated(dest%dim_ids)) deallocate(dest%dim_ids)
            allocate(dest%dim_ids, source=src%dim_ids)
        else if (allocated(dest%dim_ids)) then
            deallocate(dest%dim_ids)
        endif
        
        dest%file_var_id = src%file_var_id
        dest%xstag = src%xstag
        dest%ystag = src%ystag
        
    end subroutine assign_variable

    module subroutine set_from_metadata(this,var_id)
        implicit none
        class(variable_t), intent(inout) :: this
        integer, intent(in) :: var_id
        logical :: force_boundaries

        type(meta_data_t) :: var_meta

        var_meta = get_varmeta(var_id,force_boundaries=force_boundaries)

        this%id = var_meta%id
        this%xstag = var_meta%xstag
        this%ystag = var_meta%ystag
        this%one_d = var_meta%one_d
        this%two_d = var_meta%two_d
        this%three_d = var_meta%three_d
        this%four_d = var_meta%four_d
        this%dtype = var_meta%dtype
        this%dim_len = var_meta%dim_len
        ! this%global_dim_len = var_meta%global_dim_len
        ! this%dim_ids = var_meta%dim_ids
        ! this%file_var_id = var_meta%file_var_id
        this%force_boundaries = force_boundaries

    end subroutine set_from_metadata

end submodule
