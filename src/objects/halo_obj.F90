!>------------------------------------------------------------
!!  Implementation of an object for handling halo exchanges on behalf
!!  of the domain object
!!
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
submodule(halo_interface) halo_implementation
use icar_constants
use iso_fortran_env
use, intrinsic :: iso_c_binding

implicit none


contains


!> -------------------------------
!! Initialize the exchange arrays and dimensions
!!
!! -------------------------------
module subroutine init(this, exch_vars, adv_vars, grid, comms)
    class(halo_t), intent(inout) :: this
    type(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(grid_t), intent(in) :: grid
    type(MPI_comm), intent(inout) :: comms

    type(MPI_Group) :: comp_proc, neighbor_group
    type(c_ptr) :: tmp_ptr
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: current, my_index, n_neighbors, nx, nz, ny, ierr
    integer :: i,j,k, real_size

    CALL MPI_Type_size(MPI_REAL, real_size)

    !Some stuff we can just copy right over
    this%grid = grid
    this%halo_size = grid%halo_size

    this%north_boundary = (this%grid%yimg == this%grid%yimages)
    this%south_boundary = (this%grid%yimg == 1)
    this%east_boundary  = (this%grid%ximg == this%grid%ximages)
    this%west_boundary  = (this%grid%ximg == 1)

    this%northwest_boundary = this%north_boundary .or. this%west_boundary
    this%southeast_boundary = this%south_boundary .or. this%east_boundary
    this%southwest_boundary = this%south_boundary .or. this%west_boundary
    this%northeast_boundary = this%north_boundary .or. this%east_boundary

    this%ims = this%grid%ims; this%its = this%grid%its; this%ids = this%grid%ids
    this%ime = this%grid%ime; this%ite = this%grid%ite; this%ide = this%grid%ide
    this%kms = this%grid%kms; this%kts = this%grid%kts; this%kds = this%grid%kds
    this%kme = this%grid%kme; this%kte = this%grid%kte; this%kde = this%grid%kde
    this%jms = this%grid%jms; this%jts = this%grid%jts; this%jds = this%grid%jds
    this%jme = this%grid%jme; this%jte = this%grid%jte; this%jde = this%grid%jde

    call MPI_Comm_Rank(comms,this%halo_rank)
    !Setup the parent-child group used for buffer communication
    call MPI_Comm_Group(comms,comp_proc)

    !The following if statements are structured as follows:
    ! 1. Check if the boundary is set
    ! 2. If not, set the neighbor index using process numbering for the computational group
    ! 3. Create a group for the neighbor
    !Need to make neighbor group
    if (.not.(this%south_boundary)) then
        this%south_neighbor = this%halo_rank - this%grid%ximages
        call MPI_Group_Incl(comp_proc, 1, [this%south_neighbor], this%south_neighbor_grp)
    endif
    if (.not.(this%west_boundary)) then
        this%west_neighbor  = this%halo_rank-1
        call MPI_Group_Incl(comp_proc, 1, [this%west_neighbor], this%west_neighbor_grp)
    endif
    if (.not.(this%east_boundary)) then
        this%east_neighbor  = this%halo_rank+1
        call MPI_Group_Incl(comp_proc, 1, [this%east_neighbor], this%east_neighbor_grp)
    endif
    if (.not.(this%north_boundary)) then
        this%north_neighbor = this%halo_rank + this%grid%ximages
        call MPI_Group_Incl(comp_proc, 1, [this%north_neighbor], this%north_neighbor_grp)
    endif
    if (.not.(this%northwest_boundary)) then
        this%northwest_neighbor = this%halo_rank + this%grid%ximages - 1
        call MPI_Group_Incl(comp_proc, 1, [this%northwest_neighbor], this%northwest_neighbor_grp)
    endif
    if (.not.(this%southeast_boundary)) then
        this%southeast_neighbor = this%halo_rank - this%grid%ximages + 1
        call MPI_Group_Incl(comp_proc, 1, [this%southeast_neighbor], this%southeast_neighbor_grp)
    endif
    if (.not.(this%southwest_boundary)) then
        this%southwest_neighbor = this%halo_rank - this%grid%ximages - 1
        call MPI_Group_Incl(comp_proc, 1, [this%southwest_neighbor], this%southwest_neighbor_grp)
    endif
    if (.not.(this%northeast_boundary)) then
        this%northeast_neighbor = this%halo_rank + this%grid%ximages + 1
        call MPI_Group_Incl(comp_proc, 1, [this%northeast_neighbor], this%northeast_neighbor_grp)
    endif

    ! Detect if neighbors are on shared memory hardware
    call detect_shared_memory(this, comms)

    !Now allocate the actual 3D halo
    !We only want to set up remote windows for domain objects which are part of the actual domain
    if (.not.(comms == MPI_COMM_NULL)) then

        nx = this%grid%ns_halo_nx
        nz = this%grid%halo_nz
        ny = this%halo_size+1
        win_size = nx*nz*ny
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%south_in_win)
        call C_F_POINTER(tmp_ptr, this%south_in_3d, [nx, nz, ny])
        this%south_in_3d = 1
        
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%north_in_win)
        call C_F_POINTER(tmp_ptr, this%north_in_3d, [nx, nz, ny])
        this%north_in_3d = 1

        nx = this%halo_size+1
        nz = this%grid%halo_nz
        ny = this%grid%ew_halo_ny
        win_size = nx*nz*ny
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%east_in_win)
        call C_F_POINTER(tmp_ptr, this%east_in_3d, [nx, nz, ny])
        this%east_in_3d = 2

        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%west_in_win)
        call C_F_POINTER(tmp_ptr, this%west_in_3d, [nx, nz, ny])
        this%west_in_3d = 1

        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
    endif

    !...and the larger 3D halo for batch exchanges
    call setup_batch_exch(this, exch_vars, adv_vars, comms)

end subroutine init

module subroutine finalize(this)
    class(halo_t), intent(inout) :: this

    integer :: ierr
    
    if (this%n_2d > 0) then
        if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_2d_win)
        if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_2d_win)
        if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_2d_win)
        if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_2d_win)
        if (.not.(this%north_boundary)) call MPI_Win_Complete( this%north_2d_win)
        if (.not.(this%south_boundary)) call MPI_Win_Complete( this%south_2d_win)
        if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_2d_win)
        if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_2d_win)
        if (.not.(this%north_boundary)) call MPI_Win_Wait( this%north_2d_win)
        if (.not.(this%south_boundary)) call MPI_Win_Wait( this%south_2d_win)
        if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_2d_win)
        if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_2d_win)

        if (.not.(this%north_boundary)) call MPI_WIN_FREE(this%north_2d_win, ierr)
        if (.not.(this%south_boundary)) call MPI_WIN_FREE(this%south_2d_win, ierr)
        if (.not.(this%east_boundary)) call MPI_WIN_FREE(this%east_2d_win, ierr)
        if (.not.(this%west_boundary)) call MPI_WIN_FREE(this%west_2d_win, ierr)
        
        call MPI_Type_free(this%NS_2d_win_halo_type, ierr)
        call MPI_Type_free(this%EW_2d_win_halo_type, ierr)
    endif


    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%north_boundary)) call MPI_Win_Complete( this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete( this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_3d_win)
    if (.not.(this%north_boundary)) call MPI_Win_Wait( this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait( this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_3d_win)

    if (.not.(this%north_boundary)) call MPI_WIN_FREE(this%north_3d_win, ierr)
    if (.not.(this%south_boundary)) call MPI_WIN_FREE(this%south_3d_win, ierr)
    if (.not.(this%east_boundary)) call MPI_WIN_FREE(this%east_3d_win, ierr)
    if (.not.(this%west_boundary)) call MPI_WIN_FREE(this%west_3d_win, ierr)

    call MPI_WIN_FREE(this%north_in_win, ierr)
    call MPI_WIN_FREE(this%south_in_win, ierr)
    call MPI_WIN_FREE(this%east_in_win, ierr)
    call MPI_WIN_FREE(this%west_in_win, ierr)
    call MPI_Type_free(this%NS_3d_win_halo_type, ierr)
    call MPI_Type_free(this%EW_3d_win_halo_type, ierr)
    call MPI_Type_free(this%corner_3d_win_halo_type, ierr)

end subroutine finalize

!> -------------------------------
!! Exchange a given variable with neighboring processes.
!! This function will determine from var if the variable
!! is 2D or 3D, and if it has x- or y-staggering, and perform
!! the according halo exchange using the pre-allocated halo
!!
!! -------------------------------
module subroutine exch_var(this, var, do_dqdt, corners)
    class(halo_t),     intent(inout) :: this
    type(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt, corners
    
    integer :: xdim, ydim
    logical :: dqdt, do_corners

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    do_corners=.False.
    if (present(corners)) do_corners=corners

    if (do_corners) then
            if (.not. this%north_boundary .and. .not.this%east_boundary) call this%put_northeast(var, dqdt)
            if (.not. this%north_boundary .and. .not.this%west_boundary) call this%put_northwest(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%east_boundary)  call this%put_southeast(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%west_boundary)  call this%put_southwest(var, dqdt)

            call MPI_Win_fence(0,this%south_in_win)
            call MPI_Win_fence(0,this%north_in_win)
            call MPI_Win_fence(0,this%east_in_win)
            call MPI_Win_fence(0,this%west_in_win)

            if (.not. this%north_boundary .and. .not.this%west_boundary) call this%retrieve_northwest_halo(var, dqdt)
            if (.not. this%north_boundary .and. .not.this%east_boundary) call this%retrieve_northeast_halo(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%east_boundary)  call this%retrieve_southeast_halo(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%west_boundary)  call this%retrieve_southwest_halo(var, dqdt)

    endif

    ! if staggered in x direction, we need to carefully call the put and get commands
    if(var%xstag>0) then
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)

        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win) 

        if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
        if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)

        call MPI_Win_fence(0,this%south_in_win) 
        call MPI_Win_fence(0,this%north_in_win)

        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)

        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win) 

        if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
        if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)
    ! if staggered in y direction, we need to carefully call the put and get commands
    elseif(var%ystag>0) then
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)

        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)

        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)

        if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
        if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)

        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)                

        if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
        if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)

    else
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)

        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)
        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)

        if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
        if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)
        if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
        if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)
    endif

end subroutine exch_var


!> -------------------------------
!! Initialize the arrays and co-arrays needed to perform a batch exchange
!!
!! -------------------------------
module subroutine setup_batch_exch(this, exch_vars, adv_vars, comms)
    type(halo_t), intent(inout) :: this
    type(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(MPI_comm), intent(in) :: comms
    type(variable_t) :: var

    integer :: nx, ny, nz = 0
    type(c_ptr) :: tmp_ptr, tmp_ptr_2d
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d, win_size_corner
    integer :: ierr, i, real_size, size_out
    type(MPI_Info) :: info_in 
    type(MPI_Group) :: comp_proc, tmp_MPI_grp
    type(MPI_Comm) :: tmp_MPI_comm

    CALL MPI_Type_size(MPI_REAL, real_size)

    ! Loop over all adv_vars and count how many are 3D
    call adv_vars%reset_iterator()
    
    this%n_2d = 0
    this%n_3d = 0

    do while (adv_vars%has_more_elements())
        var = adv_vars%next()
        if (var%three_d) this%n_3d = this%n_3d + 1
    end do
    
    ! Loop over all exch vars and count how many are 3D
    call exch_vars%reset_iterator()
    
    do while (exch_vars%has_more_elements())
        var = exch_vars%next()
        if (var%three_d) this%n_3d = this%n_3d + 1
    end do
    if (STD_OUT_PE) write(*,*) "In Setup Batch Exch"
    if (STD_OUT_PE) flush(output_unit)

    ! Determine number of 2D and 3D vars present
    this%n_2d = (adv_vars%n_vars+exch_vars%n_vars)-this%n_3d

    if (.not.(comms == MPI_COMM_NULL)) then
        call MPI_Info_Create(info_in,ierr)
        call MPI_INFO_SET(info_in, 'no_locks', '.true.')
        call MPI_INFO_SET(info_in, 'same_size', '.true.')
        call MPI_INFO_SET(info_in, 'same_disp_unit', '.true.')
        ! call MPI_INFO_SET(info_in, 'alloc_shared_noncontig', '.true.')

        !First do NS
        nx = this%grid%ns_halo_nx
        nz = this%grid%halo_nz
        ny = this%halo_size
        win_size = nx*nz*ny*this%n_3d
        win_size_2d = nx*ny*this%n_2d

        call MPI_Type_contiguous(nx*ny*this%n_2d, MPI_REAL, this%NS_2d_win_halo_type)
        call MPI_Type_commit(this%NS_2d_win_halo_type)

        call MPI_Type_contiguous(nx*nz*ny*this%n_3d, MPI_REAL, this%NS_3d_win_halo_type)
        call MPI_Type_commit(this%NS_3d_win_halo_type)

        ! Group processes for creation of communication halos. 
        ! First create MPI communicator for each pair of north-south processes that will exchange data
        ! Group according to this%grid%yimg, where processes with yimg=(0,1) are in the same group, (2,3) in the same group, etc.

        !If this%grid%yimg is odd, then we allocate the north_3d_win. If it is even, we allocate the south_3d_win.
        if ((mod(this%grid%yimg,2) == 1)) then
            ! Create a group for the north-south exchange
            if (.not.(this%north_boundary)) call setup_batch_exch_north_wins(this, comms, info_in)
            if (.not.(this%south_boundary)) call setup_batch_exch_south_wins(this, comms, info_in)
        else if((mod(this%grid%yimg,2) == 0)) then
            ! Create a group for the north-south exchange
            if (.not.(this%south_boundary)) call setup_batch_exch_south_wins(this, comms, info_in)
            if (.not.(this%north_boundary)) call setup_batch_exch_north_wins(this, comms, info_in)
        endif
    
        if (.not.(this%south_boundary)) then
    
            this%south_batch_in_3d = 1
            if (this%n_2d > 0) this%south_batch_in_2d = 1
    
            if (this%south_shared) then
                call MPI_WIN_SHARED_QUERY(this%south_3d_win, 0, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%south_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then 
                    call MPI_WIN_SHARED_QUERY(this%south_2d_win, 0, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%south_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
                allocate(this%south_buffer_3d(this%n_3d,1:this%grid%ns_halo_nx,this%kms:this%kme,1:this%halo_size))
                if (this%n_2d > 0) allocate(this%south_buffer_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size))
            endif
        endif
        if (.not.(this%north_boundary)) then
            this%north_batch_in_3d = 1
            if (this%n_2d > 0) this%north_batch_in_2d = 1

            if (this%north_shared) then
                call MPI_WIN_SHARED_QUERY(this%north_3d_win, 1, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%north_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%north_2d_win, 1, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%north_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
                allocate(this%north_buffer_3d(this%n_3d,1:this%grid%ns_halo_nx,this%kms:this%kme,1:this%halo_size))
                if (this%n_2d > 0) allocate(this%north_buffer_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size))
            endif
        endif

        !Then do EW
        nx = this%halo_size
        nz = this%grid%halo_nz
        ny = this%grid%ew_halo_ny
        win_size = nx*nz*ny*this%n_3d
        win_size_2d = nx*ny*this%n_2d

        call MPI_Type_contiguous(nx*ny*this%n_2d, MPI_REAL, this%EW_2d_win_halo_type)
        call MPI_Type_commit(this%EW_2d_win_halo_type)

        call MPI_Type_contiguous(nx*nz*ny*this%n_3d, MPI_REAL, this%EW_3d_win_halo_type)
        call MPI_Type_commit(this%EW_3d_win_halo_type)

        ! Group processes for creation of communication halos. 
        ! First create MPI communicator for each pair of north-south processes that will exchange data
        ! Group according to this%grid%yimg, where processes with yimg=(0,1) are in the same group, (2,3) in the same group, etc.

        !If this%grid%yimg is odd, then we allocate the north_3d_win. If it is even, we allocate the south_3d_win.
        if ((mod(this%grid%ximg,2) == 1)) then
            ! Create a group for the east-west exchange
            if (.not.(this%east_boundary)) call setup_batch_exch_east_wins(this, comms, info_in)
            if (.not.(this%west_boundary)) call setup_batch_exch_west_wins(this, comms, info_in)
        else if((mod(this%grid%ximg,2) == 0)) then
            ! Create a group for the east-west exchange
            if (.not.(this%west_boundary)) call setup_batch_exch_west_wins(this, comms, info_in)
            if (.not.(this%east_boundary)) call setup_batch_exch_east_wins(this, comms, info_in)
        endif

        if (.not.(this%east_boundary)) then    
            this%east_batch_in_3d = 1
            if (this%n_2d > 0) this%east_batch_in_2d = 1

            if (this%east_shared) then

                call MPI_WIN_SHARED_QUERY(this%east_3d_win, 1, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%east_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%east_2d_win, 1, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%east_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
                allocate(this%east_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%grid%ew_halo_ny))
                if (this%n_2d > 0) allocate(this%east_buffer_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny))
            endif
        endif
        if (.not.(this%west_boundary)) then
            this%west_batch_in_3d = 1
            if (this%n_2d > 0) this%west_batch_in_2d = 1

            if (this%west_shared) then
                call MPI_WIN_SHARED_QUERY(this%west_3d_win, 0, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%west_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%west_2d_win, 0, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%west_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
                allocate(this%west_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%grid%ew_halo_ny))
                if (this%n_2d > 0) allocate(this%west_buffer_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny))
            endif
        endif


        !Then do corners
        win_size_corner = this%halo_size*nz*this%halo_size*this%n_3d

        call MPI_Type_contiguous(this%halo_size*nz*this%halo_size*this%n_3d, MPI_REAL, this%corner_3d_win_halo_type)
        call MPI_Type_commit(this%corner_3d_win_halo_type)

        if ((mod(this%grid%ximg,2) == 1)) then
            ! Create a group for the east-west exchange
            if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
            if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
        else if((mod(this%grid%ximg,2) == 0)) then
            ! Create a group for the east-west exchange
            if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
            if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
        endif

        if ((mod(this%grid%ximg,2) == 1)) then
            ! Create a group for the east-west exchange
            if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
            if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
        else if((mod(this%grid%ximg,2) == 0)) then
            ! Create a group for the east-west exchange
            if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
            if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
        endif

        if (.not.(this%northwest_boundary)) then

            if (this%northwest_shared) then
                call MPI_WIN_SHARED_QUERY(this%northwest_3d_win, 1, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%northwest_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%northwest_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%northwest_batch_in_3d = 1
        endif
        if (.not.(this%southeast_boundary)) then

            if (this%southeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%southeast_3d_win, 0, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%southeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%southeast_batch_in_3d = 1
        endif
        if (.not.(this%southwest_boundary)) then

            if (this%southwest_shared) then
                call MPI_WIN_SHARED_QUERY(this%southwest_3d_win, 0, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southwest_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%southwest_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%southwest_batch_in_3d = 1
        endif
        if (.not.(this%northeast_boundary)) then

            if (this%northeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%northeast_3d_win, 1, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%northeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%northeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%northeast_batch_in_3d = 1
        endif

        if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
        if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
        if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
        if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)

        if (.not.(this%northwest_boundary)) call MPI_Win_Post(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
        if (.not.(this%southeast_boundary)) call MPI_Win_Post(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
        if (.not.(this%southwest_boundary)) call MPI_Win_Post(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
        if (.not.(this%northeast_boundary)) call MPI_Win_Post(this%northeast_neighbor_grp, 0, this%northeast_3d_win)

        if (this%n_2d > 0) then
            if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_2d_win)
            if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_2d_win)
            if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_2d_win)
            if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_2d_win)
        
        endif
    endif

end subroutine setup_batch_exch


module subroutine halo_3d_send_batch(this, exch_vars, adv_vars,exch_var_only)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    logical, optional, intent(in) :: exch_var_only
    
    type(variable_t) :: var
    logical :: exch_v_only
    integer :: n, k_max, msg_size, indx
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    if (this%n_3d <= 0) return

    msg_size = 1
    disp = 0

    exch_v_only = .False.
    if (present(exch_var_only)) exch_v_only=exch_var_only

    call adv_vars%reset_iterator()
    call exch_vars%reset_iterator()
    n = 1


    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Start(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Start(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Start(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Start(this%northeast_neighbor_grp, 0, this%northeast_3d_win)

    ! Now iterate through the dictionary as long as there are more elements present. If two processors are on shared memory
    ! this step will directly copy the data to the other PE
    if (.not.(exch_v_only)) then
        do while (adv_vars%has_more_elements())
            ! get the next variable
            var = adv_vars%next()
            if (var%three_d) then
                if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+1),:,:) = &
                        var%data_3d(this%its:this%ite,:,(this%jte-this%halo_size+1):this%jte)
                if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+1),:,:) = &
                        var%data_3d(this%its:this%ite,:,this%jts:(this%jts+this%halo_size-1))
                if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,:,1:(this%jte-this%jts+1)) = &
                        var%data_3d((this%ite-this%halo_size+1):this%ite,:,this%jts:this%jte)
                if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,:,1:(this%jte-this%jts+1)) = &
                        var%data_3d(this%its:(this%its+this%halo_size-1),:,this%jts:this%jte)

                if (.not.(this%northwest_boundary)) this%northwest_buffer_3d(n,1:this%halo_size,:,1:this%halo_size) = &
                        var%data_3d(this%its:(this%its+this%halo_size-1),:,(this%jte-this%halo_size+1):this%jte)
                if (.not.(this%southeast_boundary)) this%southeast_buffer_3d(n,1:this%halo_size,:,1:this%halo_size) = &
                        var%data_3d((this%ite-this%halo_size+1):this%ite,:,this%jts:(this%jts+this%halo_size-1))
                if (.not.(this%southwest_boundary)) this%southwest_buffer_3d(n,1:this%halo_size,:,1:this%halo_size) = &
                        var%data_3d(this%its:(this%its+this%halo_size-1),:,this%jts:(this%jts+this%halo_size-1))
                if (.not.(this%northeast_boundary)) this%northeast_buffer_3d(n,1:this%halo_size,:,1:this%halo_size) = &
                        var%data_3d((this%ite-this%halo_size+1):this%ite,:,(this%jte-this%halo_size+1):this%jte)
                n = n+1
            endif
        enddo
    endif

    ! Now iterate through the exchange-only objects as long as there are more elements present
    do while (exch_vars%has_more_elements())
        ! get the next variable
        var = exch_vars%next()
        if (var%three_d) then
            k_max = ubound(var%data_3d,2)
            if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+1),1:k_max,:) = &
                    var%data_3d(this%its:this%ite,1:k_max,(this%jte-this%halo_size+1):this%jte)
            if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+1),1:k_max,:) = &
                    var%data_3d(this%its:this%ite,1:k_max,this%jts:(this%jts+this%halo_size-1))
            if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+1)) = &
                    var%data_3d((this%ite-this%halo_size+1):this%ite,1:k_max,this%jts:this%jte)
            if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+1)) = &
                    var%data_3d(this%its:(this%its+this%halo_size)-1,1:k_max,this%jts:this%jte)

            n = n+1
        endif
    enddo

    if (.not.(this%south_boundary)) then
        if (.not.(this%south_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%south_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 0, disp, msg_size, &
                this%NS_3d_win_halo_type, this%south_3d_win)
        endif
    endif

    if (.not.(this%north_boundary)) then
        if (.not.(this%north_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%north_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 1, disp, msg_size, &
                this%NS_3d_win_halo_type, this%north_3d_win)
        endif
    endif

    if (.not.(this%east_boundary)) then
        if (.not.(this%east_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%east_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 1, disp, msg_size, &
                this%EW_3d_win_halo_type, this%east_3d_win)
        endif
    endif

    if (.not.(this%west_boundary)) then
        if (.not.(this%west_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%west_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 0, disp, msg_size, &
                this%EW_3d_win_halo_type, this%west_3d_win)
        endif
    endif

    if (.not.(this%northwest_boundary)) then
        if (.not.(this%northwest_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%northwest_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, 1, disp, msg_size, &
                this%corner_3d_win_halo_type, this%northwest_3d_win)
        endif
    endif
    if (.not.(this%southeast_boundary)) then
        if (.not.(this%southeast_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%southeast_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, 0, disp, msg_size, &
                this%corner_3d_win_halo_type, this%southeast_3d_win)
        endif
    endif
    if (.not.(this%southwest_boundary)) then
        if (.not.(this%southwest_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%southwest_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, 0, disp, msg_size, &
                this%corner_3d_win_halo_type, this%southwest_3d_win)
        endif
    endif
    if (.not.(this%northeast_boundary)) then
        if (.not.(this%northeast_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%northeast_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, 1, disp, msg_size, &
                this%corner_3d_win_halo_type, this%northeast_3d_win)
        endif
    endif

    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Complete(this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Complete(this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Complete(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Complete(this%northeast_3d_win)

end subroutine halo_3d_send_batch

module subroutine halo_3d_retrieve_batch(this,exch_vars, adv_vars,exch_var_only, wait_timer)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    logical, optional, intent(in) :: exch_var_only
    type(timer_t), optional,     intent(inout)   :: wait_timer

    type(variable_t) :: var
    integer :: n, k_max
    logical :: exch_v_only

    if (this%n_3d <= 0) return

    exch_v_only = .False.
    if (present(exch_var_only)) exch_v_only=exch_var_only

    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Wait(this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Wait(this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Wait(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Wait(this%northeast_3d_win)

    call adv_vars%reset_iterator()
    call exch_vars%reset_iterator()
    n = 1
    ! Now iterate through the dictionary as long as there are more elements present
    if (.not.(exch_v_only)) then
        do while (adv_vars%has_more_elements())
            ! get the next variable
            var = adv_vars%next()
            if (var%three_d) then
                if (.not.(this%north_boundary)) var%data_3d(this%its:this%ite,:,(this%jte+1):this%jme) = &
                        this%north_batch_in_3d(n,1:(this%ite-this%its+1),:,1:this%halo_size)
                if (.not.(this%south_boundary)) var%data_3d(this%its:this%ite,:,this%jms:(this%jts-1)) = &
                        this%south_batch_in_3d(n,1:(this%ite-this%its+1),:,1:this%halo_size)
                if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,:,this%jts:this%jte) = &
                        this%east_batch_in_3d(n,:,:,1:(this%jte-this%jts+1))
                if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),:,this%jts:this%jte) = &
                        this%west_batch_in_3d(n,:,:,1:(this%jte-this%jts+1))

                if (.not.(this%northwest_boundary)) var%data_3d(this%ims:(this%its-1),:,(this%jte+1):this%jme) = &
                        this%northwest_batch_in_3d(n,:,:,:)
                if (.not.(this%southeast_boundary)) var%data_3d((this%ite+1):this%ime,:,this%jms:(this%jts-1)) = &
                        this%southeast_batch_in_3d(n,:,:,:)
                if (.not.(this%southwest_boundary)) var%data_3d(this%ims:(this%its-1),:,this%jms:(this%jts-1)) = &
                        this%southwest_batch_in_3d(n,:,:,:)
                if (.not.(this%northeast_boundary)) var%data_3d((this%ite+1):this%ime,:,(this%jte+1):this%jme) = &
                        this%northeast_batch_in_3d(n,:,:,:)
                n = n+1
            endif
        enddo
    endif
    ! Now iterate through the exchange-only objects as long as there are more elements present
    do while (exch_vars%has_more_elements())
        ! get the next variable
        var = exch_vars%next()
        if (var%three_d) then
            k_max = ubound(var%data_3d,2)
            if (.not.(this%north_boundary)) var%data_3d(this%its:this%ite,1:k_max,(this%jte+1):this%jme) = &
                    this%north_batch_in_3d(n,1:(this%ite-this%its+1),1:k_max,:)
            if (.not.(this%south_boundary)) var%data_3d(this%its:this%ite,1:k_max,this%jms:(this%jts-1)) = &
                    this%south_batch_in_3d(n,1:(this%ite-this%its+1),1:k_max,:)
            if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,1:k_max,this%jts:this%jte) = &
                    this%east_batch_in_3d(n,:,1:k_max,1:(this%jte-this%jts+1))
            if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),1:k_max,this%jts:this%jte) = &
                    this%west_batch_in_3d(n,:,1:k_max,1:(this%jte-this%jts+1))
            n = n+1
        endif
    enddo
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Post(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Post(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Post(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Post(this%northeast_neighbor_grp, 0, this%northeast_3d_win)

end subroutine halo_3d_retrieve_batch

module subroutine halo_2d_send_batch(this, exch_vars, adv_vars)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(variable_t) :: var
    integer :: n, msg_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    if (this%n_2d <= 0) return

    msg_size = 1
    disp = 0

    call exch_vars%reset_iterator()


    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_2d_win)

    n = 1
    ! Now iterate through the exchange-only objects as long as there are more elements present
    do while (exch_vars%has_more_elements())
        ! get the next variable
        var = exch_vars%next()
        if (var%two_d) then
            if (.not.(this%north_boundary)) this%north_buffer_2d(n,1:(this%ite-this%its+1),:) = &
                    var%data_2d(this%its:this%ite,(this%jte-this%halo_size+1):this%jte)
            if (.not.(this%south_boundary)) this%south_buffer_2d(n,1:(this%ite-this%its+1),:) = &
                    var%data_2d(this%its:this%ite,this%jts:(this%jts+this%halo_size-1))
            if (.not.(this%east_boundary)) this%east_buffer_2d(n,:,1:(this%jte-this%jts+1)) = &
                    var%data_2d((this%ite-this%halo_size+1):this%ite,this%jts:this%jte)
            if (.not.(this%west_boundary)) this%west_buffer_2d(n,:,1:(this%jte-this%jts+1)) = &
                    var%data_2d(this%its:(this%its+this%halo_size)-1,this%jts:this%jte)

            n = n+1
        endif
    enddo

    if (.not.(this%north_boundary)) then
        if (.not.(this%north_shared)) then
            call MPI_Put(this%north_buffer_2d, size(this%north_buffer_2d), &
                MPI_REAL, 1, disp, size(this%north_buffer_2d), MPI_REAL, this%north_2d_win)
        endif
    endif

    if (.not.(this%south_boundary)) then
        if (.not.(this%south_shared)) then
            call MPI_Put(this%south_buffer_2d, size(this%south_buffer_2d), &
                MPI_REAL, 0, disp, size(this%south_buffer_2d), MPI_REAL, this%south_2d_win)
        endif
    endif

    if (.not.(this%east_boundary)) then
        if (.not.(this%east_shared)) then
            call MPI_Put(this%east_buffer_2d, size(this%east_buffer_2d), &
                MPI_REAL, 1, disp, size(this%east_buffer_2d), MPI_REAL, this%east_2d_win)
        endif
    endif

    if (.not.(this%west_boundary)) then
        if (.not.(this%west_shared)) then
            call MPI_Put(this%west_buffer_2d, size(this%west_buffer_2d), &
                MPI_REAL, 0, disp, size(this%west_buffer_2d), MPI_REAL, this%west_2d_win)
        endif
    endif
    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_2d_win)

end subroutine halo_2d_send_batch

module subroutine halo_2d_retrieve_batch(this, exch_vars, adv_vars)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(variable_t) :: var
    integer :: n

    if (this%n_2d <= 0) return
    ! if (.not.(this%north_boundary)) call MPI_Win_fence(0, this%north_2d_win)
    ! if (.not.(this%south_boundary)) call MPI_Win_fence(0, this%south_2d_win)
    ! if (.not.(this%east_boundary)) call MPI_Win_fence(0, this%east_2d_win)
    ! if (.not.(this%west_boundary)) call MPI_Win_fence(0, this%west_2d_win)

    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_2d_win)

    call exch_vars%reset_iterator()
    n = 1    
    ! Now iterate through the exchange-only objects as long as there are more elements present
    do while (exch_vars%has_more_elements())
        ! get the next variable
        var = exch_vars%next()
        if (var%two_d) then
            if (.not.(this%north_boundary)) var%data_2d(this%its:this%ite,(this%jte+1):this%jme) = this%north_batch_in_2d(n,1:(this%ite-this%its+1),:)
            if (.not.(this%south_boundary)) var%data_2d(this%its:this%ite,this%jms:(this%jts-1)) = this%south_batch_in_2d(n,1:(this%ite-this%its+1),:)
            if (.not.(this%east_boundary)) var%data_2d((this%ite+1):this%ime,this%jts:this%jte) = this%east_batch_in_2d(n,:,1:(this%jte-this%jts+1))
            if (.not.(this%west_boundary)) var%data_2d(this%ims:(this%its-1),this%jts:this%jte) = this%west_batch_in_2d(n,:,1:(this%jte-this%jts+1))
            n = n+1
        endif
    enddo

    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_2d_win)

end subroutine halo_2d_retrieve_batch

!> -------------------------------
!! Send and get the data from all exch+adv objects to/from their neighbors (3D)
!!
!! -------------------------------
module subroutine batch_exch(this, exch_vars, adv_vars, two_d, three_d, exch_var_only)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    logical, optional, intent(in) :: two_d,three_d,exch_var_only
    
    logical :: twod, threed, exch_only
    
    exch_only = .False.
    if (present(exch_var_only)) exch_only = exch_var_only

    twod = .False.
    if(present(two_d)) twod = two_d
    threed = .True.
    if(present(three_d)) threed = three_d

    if (twod) then
        call this%halo_2d_send_batch(exch_vars, adv_vars)
     
        call this%halo_2d_retrieve_batch(exch_vars, adv_vars)
    endif
    if (threed) then
        call this%halo_3d_send_batch(exch_vars, adv_vars, exch_var_only=exch_only)

        call this%halo_3d_retrieve_batch(exch_vars, adv_vars, exch_var_only=exch_only)
    endif

end subroutine


module subroutine put_north(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, nx, offs, msg_size
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%ystag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      n = ubound(var%data_2d,2)
      nx = size(var%data_2d,1)
        call MPI_Put(var%data_2d(var%grid%ims,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
        var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
  else
      n = ubound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%ims,var%grid%kts,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%ims,var%grid%kts,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
      endif
  endif
end subroutine



module subroutine put_south(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, nx, offs, msg_size
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%ystag
  
  disp = 0
  msg_size = 1
  if (var%two_d) then
      start = lbound(var%data_2d,2)
      nx = size(var%data_2d,1)
        call MPI_Put(var%data_2d(var%grid%ims,var%grid%jts), msg_size, &
        var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
  else
      start = lbound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%ims,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%ims,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
      endif
  endif
end subroutine


module subroutine put_east(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, ny, msg_size, offs
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%xstag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      n = ubound(var%data_2d,1)
      ny = size(var%data_2d,2)
        call MPI_Put(var%data_2d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%jms), msg_size, &
        var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
  else
      n = ubound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
      endif
  endif
end subroutine




module subroutine put_west(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, ny, msg_size, offs
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%xstag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      start = lbound(var%data_2d,1)
      ny = size(var%data_2d,2)
        call MPI_Put(var%data_2d(var%grid%its,var%grid%jms), msg_size, &
        var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
  else
      start = lbound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
      endif
  endif
end subroutine

module subroutine retrieve_north_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, nx, offs_x, offs_y
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag
  offs_x=var%xstag
  
  if (var%two_d) then
      n = ubound(var%data_2d,2)
      nx = size(var%data_2d,1)
        var%data_2d(var%grid%its:var%grid%ite-offs_x,n-this%halo_size+1-offs_y:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,1,1:(this%halo_size+offs_y))
  else
      n = ubound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          var%dqdt_3d(var%grid%its:var%grid%ite-offs_x,var%grid%kts:var%grid%kte,n-this%halo_size+1-offs_y:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,var%grid%kts:var%grid%kte,1:(this%halo_size+offs_y))
      else
          var%data_3d(var%grid%its:var%grid%ite-offs_x,var%grid%kts:var%grid%kte,n-this%halo_size+1-offs_y:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,var%grid%kts:var%grid%kte,1:(this%halo_size+offs_y))
      endif
  endif
end subroutine

module subroutine retrieve_south_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, nx, offs_y, offs_x
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag
  offs_x=var%xstag

  
  if (var%two_d) then
      start = lbound(var%data_2d,2)
      nx = size(var%data_2d,1)
        var%data_2d(var%grid%its:var%grid%ite-offs_x,start:start+this%halo_size-1) = this%south_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,1,1:this%halo_size)
  else
      start = lbound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          var%dqdt_3d(var%grid%its:var%grid%ite-offs_x,var%grid%kts:var%grid%kte,start:start+this%halo_size-1) = this%south_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,var%grid%kts:var%grid%kte,1:this%halo_size)
      else
          var%data_3d(var%grid%its:var%grid%ite-offs_x,var%grid%kts:var%grid%kte,start:start+this%halo_size-1) = this%south_in_3d(1+this%halo_size:nx-this%halo_size-offs_x,var%grid%kts:var%grid%kte,1:this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_east_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, ny, offs_x, offs_y
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag
  offs_x=var%xstag

  if (var%two_d) then
      n = ubound(var%data_2d,1)
      ny = size(var%data_2d,2)
        var%data_2d(n-this%halo_size+1-offs_x:n,var%grid%jts:var%grid%jte-offs_y) = this%east_in_3d(1:(this%halo_size+offs_x),1,1+this%halo_size:ny-this%halo_size-offs_y)
  else
      n = ubound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          var%dqdt_3d(n-this%halo_size+1-offs_x:n,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte-offs_y) = this%east_in_3d(1:(this%halo_size+offs_x),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size-offs_y)
      else
          var%data_3d(n-this%halo_size+1-offs_x:n,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte-offs_y) = this%east_in_3d(1:(this%halo_size+offs_x),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size-offs_y)
      endif
  endif
end subroutine

module subroutine retrieve_west_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, ny, offs_x, offs_y
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag
  offs_x=var%xstag
  
  if (var%two_d) then
      start = lbound(var%data_2d,1)
      ny = size(var%data_2d,2)
        var%data_2d(start:start+this%halo_size-1,var%grid%jts:var%grid%jte-offs_y) = this%west_in_3d(1:this%halo_size,1,1+this%halo_size:ny-this%halo_size-offs_y)
  else
      start = lbound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          var%dqdt_3d(start:start+this%halo_size-1,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte-offs_y) = this%west_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size-offs_y)
      else
          var%data_3d(start:start+this%halo_size-1,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte-offs_y) = this%west_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size-offs_y)
      endif
  endif
end subroutine



module subroutine put_northeast(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  logical :: dqdt
  integer :: msg_size, offs_x, offs_y
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_x=var%xstag
  offs_y=var%ystag
  disp = 0
  msg_size = 1

  if (var%two_d) then
        call MPI_Put(var%data_2d(var%grid%ite-this%halo_size+1-offs_x,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
  else
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%ite-this%halo_size+1-offs_x,var%grid%kts,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%ite-this%halo_size+1-offs_x,var%grid%kts,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
      endif
  endif
end subroutine

module subroutine put_northwest(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  logical :: dqdt
  integer :: msg_size, offs_x, offs_y
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_x=var%xstag
  offs_y=var%ystag
  disp = 0
  msg_size = 1

  if (var%two_d) then
        call MPI_Put(var%data_2d(var%grid%its,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
  else
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jte-this%halo_size+1-offs_y), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
      endif
  endif
end subroutine


module subroutine put_southwest(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  logical :: dqdt
  integer :: msg_size
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  disp = 0
  msg_size = 1

  if (var%two_d) then
        call MPI_Put(var%data_2d(var%grid%its,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
  else
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
      endif
  endif
end subroutine

module subroutine put_southeast(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  logical :: dqdt
  integer :: msg_size, offs_x, offs_y
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_x=var%xstag
  offs_y=var%ystag
  disp = 0
  msg_size = 1

  if (var%two_d) then
        call MPI_Put(var%data_2d(var%grid%ite-this%halo_size+1-offs_x,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
  else
      if (dqdt) then
          call MPI_Put(var%dqdt_3d(var%grid%ite-this%halo_size+1-offs_x,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
      else
          call MPI_Put(var%data_3d(var%grid%ite-this%halo_size+1-offs_x,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
      endif
  endif
end subroutine


module subroutine retrieve_northeast_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
        var%data_2d(var%grid%ite+1-offs_x:var%grid%ime,var%grid%jte+1-offs_y:var%grid%jme) = this%east_in_3d(1:(this%halo_size+offs_x),1,1:(this%halo_size+offs_y))
  else
      if (dqdt) then
          var%dqdt_3d(var%grid%ite+1-offs_x:var%grid%ime,this%kts:this%kte,var%grid%jte+1-offs_y:var%grid%jme) = this%east_in_3d(1:(this%halo_size+offs_x),this%kts:this%kte,1:(this%halo_size+offs_y))
      else
          var%data_3d(var%grid%ite+1-offs_x:var%grid%ime,this%kts:this%kte,var%grid%jte+1-offs_y:var%grid%jme) = this%east_in_3d(1:(this%halo_size+offs_x),this%kts:this%kte,1:(this%halo_size+offs_y))
        endif
  endif
end subroutine

module subroutine retrieve_northwest_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
        var%data_2d(var%grid%ims:var%grid%its-1,var%grid%jte+1-offs_y:var%grid%jme) = this%north_in_3d(1:this%halo_size,1,1:(this%halo_size+offs_y))
  else
      if (dqdt) then
          var%dqdt_3d(var%grid%ims:var%grid%its-1,this%kts:this%kte,var%grid%jte+1-offs_y:var%grid%jme) = this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:(this%halo_size+offs_y))
      else
          var%data_3d(var%grid%ims:var%grid%its-1,this%kts:this%kte,var%grid%jte+1-offs_y:var%grid%jme) = this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:(this%halo_size+offs_y))
      endif
  endif
end subroutine

module subroutine retrieve_southwest_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt

  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
        var%data_2d(var%grid%ims:var%grid%its-1,var%grid%jms:var%grid%jts-1) = this%west_in_3d(1:this%halo_size,1,1:this%halo_size)
  else
      if (dqdt) then
          var%dqdt_3d(var%grid%ims:var%grid%its-1,this%kts:this%kte,var%grid%jms:var%grid%jts-1) = this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(var%grid%ims:var%grid%its-1,this%kts:this%kte,var%grid%jms:var%grid%jts-1) = this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_southeast_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  if (var%two_d) then
        var%data_2d(var%grid%ite+1-offs_x:var%grid%ime,var%grid%jms:var%grid%jts-1) = this%south_in_3d(1:(this%halo_size+offs_x),1,1:this%halo_size)
  else
      if (dqdt) then
          var%dqdt_3d(var%grid%ite+1-offs_x:var%grid%ime,this%kts:this%kte,var%grid%jms:var%grid%jts-1) = this%south_in_3d(1:(this%halo_size+offs_x),this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(var%grid%ite+1-offs_x:var%grid%ime,this%kts:this%kte,var%grid%jms:var%grid%jts-1) = this%south_in_3d(1:(this%halo_size+offs_x),this%kts:this%kte,1:this%halo_size)
      endif
  endif
end subroutine

!> -------------------------------
!! Detect if neighbors are on shared memory hardware
!!
!! -------------------------------
module subroutine detect_shared_memory(this, comms)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    
    integer :: ierr, shared_size, shared_rank
    type(MPI_Comm) :: shared_comm
    type(MPI_Group) :: shared_comm_grp
    integer :: neighbor_shared_rank(1), neighbor_rank
    
    ! Create communicator for processes that can create shared memory
    call MPI_Comm_split_type(comms, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
    
    ! Get rank and size in shared memory communicator
    call MPI_Comm_rank(shared_comm, shared_rank, ierr)
    call MPI_Comm_size(shared_comm, shared_size, ierr)
    call MPI_Comm_group(shared_comm, shared_comm_grp, ierr)

    ! Initialize all shared flags to false
    this%north_shared = .false.
    this%south_shared = .false.
    this%east_shared = .false.
    this%west_shared = .false.
    this%northwest_shared = .false.
    this%northeast_shared = .false.
    this%southwest_shared = .false.
    this%southeast_shared = .false.
    
    ! Check if each neighbor is in same shared memory space
    if (.not. this%north_boundary) then
        neighbor_rank = this%north_neighbor
        call MPI_Group_translate_ranks(this%north_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%north_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%south_boundary) then
        neighbor_rank = this%south_neighbor
        call MPI_Group_translate_ranks(this%south_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%south_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%east_boundary) then
        neighbor_rank = this%east_neighbor
        call MPI_Group_translate_ranks(this%east_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%east_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%west_boundary) then
        neighbor_rank = this%west_neighbor
        call MPI_Group_translate_ranks(this%west_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%west_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif

    if (.not. this%northwest_boundary) then
        neighbor_rank = this%northwest_neighbor
        call MPI_Group_translate_ranks(this%northwest_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%northwest_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%northeast_boundary) then
        neighbor_rank = this%northeast_neighbor
        call MPI_Group_translate_ranks(this%northeast_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%northeast_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%southwest_boundary) then
        neighbor_rank = this%southwest_neighbor
        call MPI_Group_translate_ranks(this%southwest_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%southwest_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%southeast_boundary) then
        neighbor_rank = this%southeast_neighbor
        call MPI_Group_translate_ranks(this%southeast_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%southeast_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    ! Free the shared memory communicator
    call MPI_Comm_free(shared_comm, ierr)
    
    ! Set whether to use shared windows based on if any neighbors are shared
    this%use_shared_windows = this%north_shared .or. this%south_shared .or. &
                             this%east_shared .or. this%west_shared
end subroutine detect_shared_memory

module subroutine setup_batch_exch_north_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nx, ny, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr, tmp_ptr_2d

    CALL MPI_Type_size(MPI_REAL, real_size)

    !NS
    nx = this%grid%ns_halo_nx
    ny = this%halo_size
    nz = this%grid%halo_nz
    win_size = nx*nz*ny*this%n_3d
    win_size_2d = nx*ny*this%n_2d

    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%north_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%north_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%north_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%north_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%north_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%north_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%north_batch_in_2d, [this%n_2d, nx, ny])
end subroutine setup_batch_exch_north_wins

module subroutine setup_batch_exch_south_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nx, ny, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr, tmp_ptr_2d

    CALL MPI_Type_size(MPI_REAL, real_size)

    !NS
    nx = this%grid%ns_halo_nx
    ny = this%halo_size
    nz = this%grid%halo_nz
    win_size = nx*nz*ny*this%n_3d
    win_size_2d = nx*ny*this%n_2d

    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%south_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%south_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%south_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%south_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%south_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%south_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%south_batch_in_2d, [this%n_2d, nx, ny])
end subroutine setup_batch_exch_south_wins

module subroutine setup_batch_exch_east_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nx, ny, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr, tmp_ptr_2d

    CALL MPI_Type_size(MPI_REAL, real_size)

    !EW
    nx = this%halo_size
    nz = this%grid%halo_nz
    ny = this%grid%ew_halo_ny
    win_size = nx*nz*ny*this%n_3d
    win_size_2d = nx*ny*this%n_2d

    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%east_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%east_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%east_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%east_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%east_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%east_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%east_batch_in_2d, [this%n_2d, nx, ny])

end subroutine setup_batch_exch_east_wins

module subroutine setup_batch_exch_west_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nx, ny, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr, tmp_ptr_2d

    CALL MPI_Type_size(MPI_REAL, real_size)

    !EW
    nx = this%halo_size
    nz = this%grid%halo_nz
    ny = this%grid%ew_halo_ny
    win_size = nx*nz*ny*this%n_3d
    win_size_2d = nx*ny*this%n_2d

    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%west_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%west_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%west_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%west_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%west_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%west_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%west_batch_in_2d, [this%n_2d, nx, ny])
end subroutine setup_batch_exch_west_wins

module subroutine setup_batch_exch_northwest_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr

    CALL MPI_Type_size(MPI_REAL, real_size)

    !NW
    nz = this%grid%halo_nz
    win_size = this%halo_size*nz*this%halo_size*this%n_3d

    ! Create a group for the northwest-southwest exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages-1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%northwest_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northwest_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northwest_3d_win)
    endif
    call C_F_POINTER(tmp_ptr, this%northwest_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])

end subroutine setup_batch_exch_northwest_wins

module subroutine setup_batch_exch_northeast_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr

    CALL MPI_Type_size(MPI_REAL, real_size)

    !NE
    nz = this%grid%halo_nz
    win_size = this%halo_size*nz*this%halo_size*this%n_3d

    ! Create a group for the northeast-southeast exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages+1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%northeast_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northeast_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northeast_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%northeast_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])

end subroutine setup_batch_exch_northeast_wins

module subroutine setup_batch_exch_southwest_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr

    CALL MPI_Type_size(MPI_REAL, real_size)

    !SW
    nz = this%grid%halo_nz
    win_size = this%halo_size*nz*this%halo_size*this%n_3d

    ! Create a group for the southwest-southeast exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages-1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%southwest_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southwest_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southwest_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%southwest_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])

end subroutine setup_batch_exch_southwest_wins

module subroutine setup_batch_exch_southeast_wins(this, comms, info_in)
    class(halo_t), intent(inout) :: this
    type(MPI_comm), intent(in) :: comms
    type(MPI_Info), intent(in) :: info_in

    integer :: ierr, nz
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: real_size
    type(MPI_Group) :: tmp_MPI_grp, comp_proc
    type(MPI_Comm) :: tmp_MPI_comm
    type(C_PTR) :: tmp_ptr

    CALL MPI_Type_size(MPI_REAL, real_size)

    !SE
    nz = this%grid%halo_nz
    win_size = this%halo_size*nz*this%halo_size*this%n_3d

    ! Create a group for the southeast-southwest exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages+1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

    if (this%southeast_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southeast_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southeast_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%southeast_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])

end subroutine setup_batch_exch_southeast_wins

end submodule
