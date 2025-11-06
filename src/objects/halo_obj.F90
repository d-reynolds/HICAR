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
use output_metadata,            only : get_varmeta, get_varindx
use meta_data_interface,        only : meta_data_t
use, intrinsic :: iso_c_binding

implicit none


contains


!> -------------------------------
!! Initialize the exchange arrays and dimensions
!!
!! -------------------------------
module subroutine init_halo(this, exch_vars, adv_vars, grid, comms)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(in) :: adv_vars(:), exch_vars(:)
    type(grid_t), intent(in) :: grid
    type(MPI_comm), intent(inout) :: comms

    type(MPI_Group) :: comp_proc, neighbor_group
    type(c_ptr) :: tmp_ptr
    type(MPI_Info) :: info_in
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

        call MPI_Info_create(info_in, ierr)
#ifdef _OPENACC
        ! Check if MPI supports CUDA-aware memory
        call MPI_Info_set(info_in, "alloc_shm", ".true.", ierr)
        ! Set memory type hint for GPU accessibility
        call MPI_Info_set(info_in, "alloc_mem", "device", ierr)
        call MPI_Info_set(info_in, "mpi_assert_memory_alloc_kinds", "gpu:device", ierr)
        call MPI_Info_set(info_in, "cuda_aware", ".true.", ierr)
#endif
        nx = this%grid%ns_halo_nx
        nz = this%grid%halo_nz
        ny = this%halo_size+1
        win_size = nx*nz*ny

#ifdef _OPENACC
        allocate(this%south_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%south_in_3d)
        !$acc host_data use_device(this%south_in_3d)
        call MPI_WIN_CREATE(this%south_in_3d, win_size*real_size, real_size, info_in, comms, this%south_in_win, ierr)
        !$acc end host_data
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%south_in_win)
        call C_F_POINTER(tmp_ptr, this%south_in_3d, [nx, nz, ny])
#endif
        allocate(this%south_in_buffer(1:this%grid%ns_halo_nx, 1:nz, 1:this%halo_size+1))
        !$acc enter data copyin(this%south_in_buffer)

#ifdef _OPENACC
        allocate(this%north_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%north_in_3d)
        !$acc host_data use_device(this%north_in_3d)
        call MPI_WIN_CREATE(this%north_in_3d, win_size*real_size, real_size, info_in, comms, this%north_in_win, ierr)
        !$acc end host_data
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%north_in_win)
        call C_F_POINTER(tmp_ptr, this%north_in_3d, [nx, nz, ny])
#endif
        allocate(this%north_in_buffer(1:this%grid%ns_halo_nx, 1:nz, 1:this%halo_size+1))
        !$acc enter data copyin(this%north_in_buffer)

        nx = this%halo_size+1
        nz = this%grid%halo_nz
        ny = this%grid%ew_halo_ny
        win_size = nx*nz*ny

#ifdef _OPENACC
        allocate(this%east_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%east_in_3d)
        !$acc host_data use_device(this%east_in_3d)
        call MPI_WIN_CREATE(this%east_in_3d, win_size*real_size, real_size, info_in, comms, this%east_in_win, ierr)
        !$acc end host_data
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%east_in_win)
        call C_F_POINTER(tmp_ptr, this%east_in_3d, [nx, nz, ny])
#endif
        allocate(this%east_in_buffer(1:this%halo_size+1, 1:nz, 1:this%grid%ew_halo_ny))
        !$acc enter data copyin(this%east_in_buffer)

#ifdef _OPENACC
        allocate(this%west_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%west_in_3d)
        !$acc host_data use_device(this%west_in_3d)
        call MPI_WIN_CREATE(this%west_in_3d, win_size*real_size, real_size, info_in, comms, this%west_in_win, ierr)
        !$acc end host_data
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%west_in_win)
        call C_F_POINTER(tmp_ptr, this%west_in_3d, [nx, nz, ny])
#endif
        allocate(this%west_in_buffer(1:this%halo_size+1, 1:nz, 1:this%grid%ew_halo_ny))
        !$acc enter data copyin(this%west_in_buffer)
        this%north_in_3d = 1
        this%south_in_3d = 1
        this%east_in_3d = 1
        this%west_in_3d = 1

    endif

    !...and the larger 3D halo for batch exchanges
    call setup_batch_exch(this, exch_vars, adv_vars, comms)

end subroutine init_halo

module subroutine finalize(this)
    implicit none
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

    !$acc exit data finalize delete(this%north_in_3d, this%south_in_3d, this%east_in_3d, this%west_in_3d, &
    !$acc                      this%north_in_buffer, this%south_in_buffer, this%east_in_buffer, this%west_in_buffer, &
    !$acc                      this%south_buffer_2d, this%north_buffer_2d, this%east_buffer_2d, this%west_buffer_2d, &
    !$acc                      this%south_buffer_3d, this%north_buffer_3d, this%east_buffer_3d, this%west_buffer_3d, &
    !$acc                      this%northwest_buffer_3d, this%southwest_buffer_3d, this%northeast_buffer_3d, this%southeast_buffer_3d, &
    !$acc                      this%southwest_batch_in_3d, this%northwest_batch_in_3d, this%southeast_batch_in_3d, this%northeast_batch_in_3d, &
    !$acc                      this%south_batch_in_3d, this%north_batch_in_3d, this%east_batch_in_3d, this%west_batch_in_3d, &
    !$acc                      this%south_batch_in_2d, this%north_batch_in_2d, this%east_batch_in_2d, this%west_batch_in_2d)

end subroutine finalize

!> -------------------------------
!! Exchange a given variable with neighboring processes.
!! This function will determine from var if the variable
!! is 2D or 3D, and if it has x- or y-staggering, and perform
!! the according halo exchange using the pre-allocated halo
!!
!! -------------------------------
module subroutine exch_var(this, var, do_dqdt, corners)
    implicit none
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

        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)

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

end subroutine exch_var


!> -------------------------------
!! Initialize the arrays and co-arrays needed to perform a batch exchange
!!
!! -------------------------------
subroutine setup_batch_exch(this, exch_vars, adv_vars, comms)
    implicit none
    type(halo_t), intent(inout) :: this
    type(index_type), intent(in) :: adv_vars(:), exch_vars(:)
    type(MPI_comm), intent(in) :: comms
    type(meta_data_t) :: var

    integer :: nx, ny, nz = 0
    type(c_ptr) :: tmp_ptr, tmp_ptr_2d
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d, win_size_corner
    integer :: ierr, i, real_size, size_out
    type(MPI_Info) :: info_in 
    type(MPI_Group) :: comp_proc, tmp_MPI_grp
    type(MPI_Comm) :: tmp_MPI_comm

    CALL MPI_Type_size(MPI_REAL, real_size)    


    if (STD_OUT_PE) write(*,*) "In Setup Batch Exch"
    if (STD_OUT_PE) flush(output_unit)

    this%n_2d = 0
    this%n_3d = 0

    ! Loop over all adv vars and count how many are 3D
    do i = 1,size(adv_vars)
        var = get_varmeta(adv_vars(i)%id)
        if (var%three_d) then
            this%n_3d = this%n_3d + 1
        end if
    end do

    ! Loop over all exch vars and count how many are 3D
    do i = 1,size(exch_vars)
        var = get_varmeta(exch_vars(i)%id)
        if (var%three_d) then
            this%n_3d = this%n_3d + 1
        end if
    end do

    ! Determine number of 2D and 3D vars present
    this%n_2d = (size(adv_vars)+size(exch_vars))-this%n_3d

    if (.not.(comms == MPI_COMM_NULL)) then
        call MPI_Info_Create(info_in,ierr)
        call MPI_INFO_SET(info_in, 'no_locks', '.true.')
        call MPI_INFO_SET(info_in, 'same_size', '.true.')
        call MPI_INFO_SET(info_in, 'same_disp_unit', '.true.')
        ! call MPI_INFO_SET(info_in, 'alloc_shared_noncontig', '.true.')
#ifdef _OPENACC
        ! Check if MPI supports CUDA-aware memory
        call MPI_Info_set(info_in, "alloc_shm", ".true.", ierr)
        ! Set memory type hint for GPU accessibility
        call MPI_Info_set(info_in, "alloc_mem", "device", ierr)
        call MPI_Info_set(info_in, "mpi_assert_memory_alloc_kinds", "gpu:device", ierr)
        call MPI_Info_set(info_in, "cuda_aware", ".true.", ierr)
#endif

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
            !$acc enter data copyin(this%south_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%south_buffer_2d)
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
            !$acc enter data copyin(this%north_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%north_buffer_2d)
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
            !$acc enter data copyin(this%east_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%east_buffer_2d)
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
            !$acc enter data copyin(this%west_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%west_buffer_2d)
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
            !$acc enter data copyin(this%northwest_buffer_3d)
        endif
        if (.not.(this%southeast_boundary)) then

            if (this%southeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%southeast_3d_win, 0, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%southeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%southeast_batch_in_3d = 1
            !$acc enter data copyin(this%southeast_buffer_3d)
        endif
        if (.not.(this%southwest_boundary)) then

            if (this%southwest_shared) then
                call MPI_WIN_SHARED_QUERY(this%southwest_3d_win, 0, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southwest_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%southwest_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%southwest_batch_in_3d = 1
            !$acc enter data copyin(this%southwest_buffer_3d)
        endif
        if (.not.(this%northeast_boundary)) then

            if (this%northeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%northeast_3d_win, 1, win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%northeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
                allocate(this%northeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
            endif
            this%northeast_batch_in_3d = 1
            !$acc enter data copyin(this%northeast_buffer_3d)
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


module subroutine halo_3d_send_batch(this, exch_vars, adv_vars, var_data, exch_var_only)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: adv_vars(:), exch_vars(:)
    type(variable_t), intent(inout) :: var_data(:)
    logical, optional, intent(in) :: exch_var_only
    
    logical :: exch_v_only
    integer :: n, p, k_max, msg_size, indx, i, j, k, n_vars
    integer :: kms, kme, its, ite, jts, jte, halo_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    if (this%n_3d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    msg_size = 1
    disp = 0

    halo_size = this%halo_size
    kms = this%kms; kme = this%kme
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

    exch_v_only = .False.
    if (present(exch_var_only)) exch_v_only=exch_var_only

    n = 1
    n_vars = size(var_data)

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
        do p = 1, size(adv_vars)
            if (adv_vars(p)%v <= n_vars) then
                !check that the variable indexed matches the name
                if (var_data(adv_vars(p)%v)%id == adv_vars(p)%id) then
                    associate(var => var_data(adv_vars(p)%v)%data_3d)

                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = its,ite
                            buff(n,i-its+1,k,j) = var(i,k,(jte-halo_size+j))
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = its,ite
                            buff(n,i-its+1,k,j) = var(i,k,jts+j-1)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = jts,jte
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j-jts+1) = var((ite-halo_size+i),k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = jts,jte
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j-jts+1) = var(its+i-1,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%northwest_boundary)) then
                        associate(buff => this%northwest_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j) = var(its+i-1,k,(jte-halo_size+j))
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%southeast_boundary)) then
                        associate(buff => this%southeast_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j) = var((ite-halo_size+i),k,jts+j-1)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    
                    if (.not.(this%southwest_boundary)) then
                        associate(buff => this%southwest_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j) = var(its+i-1,k,jts+j-1)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%northeast_boundary)) then
                        associate(buff => this%northeast_buffer_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            buff(n,i,k,j) = var((ite-halo_size+i),k,(jte-halo_size+j))
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    end associate
                    n = n+1
                endif
            endif
        enddo
    endif

    ! Now iterate through the exchange-only objects as long as there are more elements present
     do p = 1, size(exch_vars)
        !check that the variable indexed matches the name
        if (exch_vars(p)%v <= n_vars) then
            if (var_data(exch_vars(p)%v)%id == exch_vars(p)%id) then
                associate(var => var_data(exch_vars(p)%v)%data_3d)
                k_max = ubound(var,2)

                if (.not.(this%north_boundary)) then
                    associate(buff => this%north_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = 1,k_max
                    do i = its,ite
                        buff(n,i-its+1,k,j) = var(i,k,(jte-halo_size+j))
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%south_boundary)) then
                    associate(buff => this%south_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = 1,k_max
                    do i = its,ite
                        buff(n,i-its+1,k,j) = var(i,k,jts+j-1)
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%east_boundary)) then
                    associate(buff => this%east_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = jts,jte
                    do k = 1,k_max
                    do i = 1,halo_size
                        buff(n,i,k,j-jts+1) = var((ite-halo_size+i),k,j)
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%west_boundary)) then
                    associate(buff => this%west_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = jts,jte
                    do k = 1,k_max
                    do i = 1,halo_size
                        buff(n,i,k,j-jts+1) = var(its+i-1,k,j)
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                end associate
                n = n+1
            endif
        endif
    enddo

    !$acc data present(this)


    !$acc host_data use_device(this%south_buffer_3d, this%north_buffer_3d, &
    !$acc this%east_buffer_3d, this%west_buffer_3d, &
    !$acc this%northwest_buffer_3d, this%southeast_buffer_3d, &
    !$acc this%southwest_buffer_3d, this%northeast_buffer_3d)
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
    !$acc end host_data
    !$acc end data

    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Complete(this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Complete(this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Complete(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Complete(this%northeast_3d_win)

end subroutine halo_3d_send_batch

module subroutine halo_3d_retrieve_batch(this,exch_vars, adv_vars, var_data, exch_var_only, wait_timer)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: adv_vars(:), exch_vars(:)
    type(variable_t), intent(inout) :: var_data(:)
    logical, optional, intent(in) :: exch_var_only
    type(timer_t), optional,     intent(inout)   :: wait_timer

    integer :: n, p, k_max, i, j, k, n_vars
    integer :: halo_size, kms, kme, ims, jms, its, ite, jts, jte
    logical :: exch_v_only

    if (this%n_3d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    halo_size = this%halo_size
    kms = this%kms; kme = this%kme
    ims = this%ims; jms = this%jms
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

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

    n = 1
    n_vars = size(var_data)
    if (present(wait_timer)) call wait_timer%start()

    ! Now iterate through the dictionary as long as there are more elements present
    if (.not.(exch_v_only)) then
        do p = 1, size(adv_vars)
            if (adv_vars(p)%v <= n_vars) then
                !check that the variable indexed matches the name
                if (var_data(adv_vars(p)%v)%id == adv_vars(p)%id) then
                    associate(var => var_data(adv_vars(p)%v)%data_3d)

                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = its,ite
                            var(i,k,(jte+j)) = buff(n,i-its+1,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = its,ite
                            var(i,k,jms+j-1) = buff(n,i-its+1,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = jts,jte
                        do k = kms,kme
                        do i = 1,halo_size
                            var((ite+i),k,j) = buff(n,i,k,j-jts+1)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = jts,jte
                        do k = kms,kme
                        do i = 1,halo_size
                            var(ims+i-1,k,j) = buff(n,i,k,j-jts+1)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    
                    if (.not.(this%northwest_boundary)) then
                        associate(buff => this%northwest_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            var(ims+i-1,k,(jte+j)) = buff(n,i,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%southeast_boundary)) then
                        associate(buff => this%southeast_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            var((ite+i),k,jms+j-1) = buff(n,i,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%southwest_boundary)) then
                        associate(buff => this%southwest_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            var(ims+i-1,k,jms+j-1) = buff(n,i,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif

                    if (.not.(this%northeast_boundary)) then
                        associate(buff => this%northeast_batch_in_3d)
                        !$acc parallel loop gang vector collapse(3) present(var, buff)
                        do j = 1,halo_size
                        do k = kms,kme
                        do i = 1,halo_size
                            var((ite+i),k,(jte+j)) = buff(n,i,k,j)
                        enddo
                        enddo
                        enddo
                        end associate
                    endif
                    n = n+1
                    end associate
                endif
            endif
        enddo
    endif
    ! Now iterate through the exchange-only objects as long as there are more elements present
    do p = 1, size(exch_vars)
        if (exch_vars(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(exch_vars(p)%v)%id == exch_vars(p)%id) then
                associate(var => var_data(exch_vars(p)%v)%data_3d)
                k_max = ubound(var,2)

                if (.not.(this%north_boundary)) then
                    associate(buff => this%north_batch_in_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = 1,k_max
                    do i = its,ite
                        var(i,k,(jte+j)) = buff(n,i-its+1,k,j)
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                if (.not.(this%south_boundary)) then
                    associate(buff => this%south_batch_in_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = 1,k_max
                    do i = its,ite
                        var(i,k,jms+j-1) = buff(n,i-its+1,k,j)
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                if (.not.(this%east_boundary)) then
                    associate(buff => this%east_batch_in_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = jts,jte
                    do k = 1,k_max
                    do i = 1,halo_size
                        var((ite+i),k,j) = buff(n,i,k,j-jts+1)
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                if (.not.(this%west_boundary)) then
                    associate(buff => this%west_batch_in_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = jts,jte
                    do k = 1,k_max
                    do i = 1,halo_size
                        var(ims+i-1,k,j) = buff(n,i,k,j-jts+1)
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                n = n+1
                end associate
            endif
        endif
    enddo

    if (present(wait_timer)) call wait_timer%stop()

    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Post(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Post(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Post(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Post(this%northeast_neighbor_grp, 0, this%northeast_3d_win)

end subroutine halo_3d_retrieve_batch

module subroutine halo_2d_send_batch(this, exch_vars, adv_vars, var_data)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: adv_vars(:), exch_vars(:)
    type(variable_t), intent(inout) :: var_data(:)
    integer :: n, p, msg_size, i, j, k, n_vars
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    if (this%n_2d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    msg_size = 1
    disp = 0

    ! call exch_vars%reset_iterator()


    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_2d_win)

    n = 1
    n_vars = size(var_data)

    !$acc data present(this, exch_vars, adv_vars, var_data)

    ! Now iterate through the exchange-only objects as long as there are more elements present
    do p = 1, size(exch_vars)
        if (exch_vars(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(exch_vars(p)%v)%id == exch_vars(p)%id) then
                    if (.not.(this%north_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = 1,this%halo_size
                        do i = this%its,this%ite
                        this%north_buffer_2d(n,i-this%its+1,j) = &
                            var_data(exch_vars(p)%v)%data_2d(i,(this%jte-this%halo_size+j))
                        enddo
                        enddo
                    endif
                    if (.not.(this%south_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = 1,this%halo_size
                        do i = this%its,this%ite
                        this%south_buffer_2d(n,i-this%its+1,j) = &
                                var_data(exch_vars(p)%v)%data_2d(i,this%jts+j-1)
                        enddo
                        enddo
                    endif
                    if (.not.(this%east_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = this%jts,this%jte
                        do i = 1,this%halo_size
                        this%east_buffer_2d(n,i,j-this%jts+1) = &
                            var_data(exch_vars(p)%v)%data_2d((this%ite-this%halo_size+i),j)
                        enddo
                        enddo
                    endif
                    if (.not.(this%west_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = this%jts,this%jte
                        do i = 1,this%halo_size
                        this%west_buffer_2d(n,i,j-this%jts+1) = &
                            var_data(exch_vars(p)%v)%data_2d(this%its+i-1,j)
                        enddo
                        enddo
                    endif
                n = n+1
            endif
        endif
    enddo

    !$acc host_data use_device(this%south_buffer_2d, this%north_buffer_2d, &
    !$acc this%east_buffer_2d, this%west_buffer_2d)
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
    !$acc end host_data
    !$acc end data
    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_2d_win)

end subroutine halo_2d_send_batch

module subroutine halo_2d_retrieve_batch(this, exch_vars, adv_vars, var_data)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: adv_vars(:), exch_vars(:)
    type(variable_t), intent(inout) :: var_data(:)
    integer :: n, p, i, j, k, n_vars

    if (this%n_2d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_2d_win)

    n = 1    
    n_vars = size(var_data)

    !$acc data present(this, exch_vars, adv_vars, var_data)

    ! Now iterate through the exchange-only objects as long as there are more elements present
    do p = 1, size(exch_vars)
        if (exch_vars(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(exch_vars(p)%v)%id == exch_vars(p)%id) then
                    if (.not.(this%north_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = 1,this%halo_size
                        do i = this%its,this%ite
                        var_data(exch_vars(p)%v)%data_2d(i,(this%jte+j)) = &
                                this%north_batch_in_2d(n,i-this%its+1,j)
                        enddo
                        enddo
                    endif
                    if (.not.(this%south_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = 1,this%halo_size
                        do i = this%its,this%ite
                        var_data(exch_vars(p)%v)%data_2d(i,this%jms+j-1) = &
                                this%south_batch_in_2d(n,i-this%its+1,j)
                        enddo
                        enddo
                    endif
                    if (.not.(this%east_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = this%jts,this%jte
                        do i = 1,this%halo_size
                        var_data(exch_vars(p)%v)%data_2d((this%ite+i),j) = &
                                this%east_batch_in_2d(n,i,j-this%jts+1)
                        enddo
                        enddo
                    endif
                    if (.not.(this%west_boundary)) then
                        !$acc parallel loop gang vector collapse(2)
                        do j = this%jts,this%jte
                        do i = 1,this%halo_size
                        var_data(exch_vars(p)%v)%data_2d(this%ims+i-1,j) = &
                                this%west_batch_in_2d(n,i,j-this%jts+1)
                        enddo
                        enddo
                    endif
                n = n+1
            endif
        endif
    enddo
    !$acc end data

    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_2d_win)

end subroutine halo_2d_retrieve_batch

!> -------------------------------
!! Send and get the data from all exch+adv objects to/from their neighbors (3D)
!!
!! -------------------------------
! module subroutine batch_exch(this, exch_vars, adv_vars, two_d, three_d, var_data, exch_var_only)
!     class(halo_t), intent(inout) :: this
!     type(index_type), intent(inout) :: adv_vars(:), exch_vars(:)
!     type(variable_t), intent(inout) :: var_data(:)
!     logical, optional, intent(in) :: two_d,three_d,exch_var_only
    
!     logical :: twod, threed, exch_only
    
!     exch_only = .False.
!     if (present(exch_var_only)) exch_only = exch_var_only

!     twod = .False.
!     if(present(two_d)) twod = two_d
!     threed = .True.
!     if(present(three_d)) threed = three_d

!     if (twod) then
!         call this%halo_2d_send_batch(exch_vars, adv_vars)
     
!         call this%halo_2d_retrieve_batch(exch_vars, adv_vars)
!     endif
!     if (threed) then
!         call this%halo_3d_send_batch(exch_vars, adv_vars, exch_var_only=exch_only)

!         call this%halo_3d_retrieve_batch(exch_vars, adv_vars, exch_var_only=exch_only)
!     endif

! end subroutine


module subroutine put_north(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, nx, offs, msg_size, i, k, j, indx_start
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%ystag
  disp = 0
  msg_size = 1
  indx_start = var%grid%jte-offs-var%grid%halo_size+1

  !$acc data present(this%north_in_buffer)

  if (var%two_d) then
    !$acc parallel loop gang vector collapse(2) present(var%data_2d)
    do j = indx_start, var%grid%jte
        do i = var%grid%its, var%grid%ite
            this%north_in_buffer(i-var%grid%its+1,1,j-indx_start+1) = var%data_2d(i,j)
        enddo
    enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = indx_start, var%grid%jte
            do k = var%grid%kts, var%grid%kte
                do i = var%grid%its, var%grid%ite
                    this%north_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-indx_start+1) = var%dqdt_3d(i,k,j)
                enddo
            enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = indx_start, var%grid%jte
            do k = var%grid%kts, var%grid%kte
                do i = var%grid%its, var%grid%ite
                    this%north_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-indx_start+1) = var%data_3d(i,k,j)
                enddo
            enddo
          enddo
      endif
  endif

  !$acc host_data use_device(this%north_in_buffer)
  call MPI_Put(this%north_in_buffer, msg_size, &
      var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
  !$acc end host_data

  !$acc end data
end subroutine



module subroutine put_south(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, nx, offs, msg_size, i, k, j, indx_end
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%ystag
  disp = 0
  msg_size = 1
  indx_end = var%grid%jts+offs+var%grid%halo_size-1

  !$acc data present(this%south_in_buffer)

  if (var%two_d) then
          !$acc parallel loop gang vector collapse(2) present(var%data_2d)
          do j = var%grid%jts, indx_end
            do i = var%grid%its, var%grid%ite
              this%south_in_buffer(i-var%grid%its+1,1,j-var%grid%jts+1) = var%data_2d(i,j)
            enddo
          enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jts, indx_end
            do k = var%grid%kts, var%grid%kte
              do i = var%grid%its, var%grid%ite
                this%south_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%dqdt_3d(i,k,j)
              enddo
            enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jts, indx_end
            do k = var%grid%kts, var%grid%kte
              do i = var%grid%its, var%grid%ite
                this%south_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%data_3d(i,k,j)
              enddo
            enddo
          enddo
      endif
  endif

  !$acc host_data use_device(this%south_in_buffer)
  call MPI_Put(this%south_in_buffer, msg_size, &
      var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
  !$acc end host_data

  !$acc end data

end subroutine


module subroutine put_east(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, ny, msg_size, offs, i, k, j, indx_start
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%xstag
  disp = 0
  msg_size = 1
  indx_start = var%grid%ite-offs-var%grid%halo_size+1

  !$acc data present(this%east_in_buffer)

  if (var%two_d) then
    !$acc parallel loop gang vector collapse(2) present(var%data_2d)
    do j = var%grid%jts, var%grid%jte
        do i = indx_start, var%grid%ite
            this%east_in_buffer(i-indx_start+1,1,j-var%grid%jts+1) = var%data_2d(i,j)
        enddo
    enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jts, var%grid%jte
            do k = var%grid%kts, var%grid%kte
                do i = indx_start, var%grid%ite
                    this%east_in_buffer(i-indx_start+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%dqdt_3d(i,k,j)
                enddo
            enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jts, var%grid%jte
            do k = var%grid%kts, var%grid%kte
                do i = indx_start, var%grid%ite
                    this%east_in_buffer(i-indx_start+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%data_3d(i,k,j)
                enddo
            enddo
          enddo
      endif
  endif

  !$acc host_data use_device(this%east_in_buffer)
  call MPI_Put(this%east_in_buffer, msg_size, &
  var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
  !$acc end host_data

  !$acc end data
end subroutine




module subroutine put_west(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, ny, msg_size, offs, i, k, j, indx_end
  logical :: dqdt
  INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=var%xstag
  disp = 0
  msg_size = 1
  indx_end = var%grid%its+offs+var%grid%halo_size-1


    !$acc data present(this%west_in_buffer)

  if (var%two_d) then
          !$acc parallel loop gang vector collapse(2) present(var%data_2d)
          do j = var%grid%jts, var%grid%jte
            do i = var%grid%its, indx_end
              this%west_in_buffer(i-var%grid%its+1,1,j-var%grid%jts+1) = var%data_2d(i,j)
            enddo
          enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jts, var%grid%jte
            do k = var%grid%kts, var%grid%kte
              do i = var%grid%its, indx_end
                this%west_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%dqdt_3d(i,k,j)
              enddo
            enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jts, var%grid%jte
            do k = var%grid%kts, var%grid%kte
              do i = var%grid%its, indx_end
                this%west_in_buffer(i-var%grid%its+1,k-var%grid%kts+1,j-var%grid%jts+1) = var%data_3d(i,k,j)
              enddo
            enddo
          enddo
      endif
  endif

  !$acc host_data use_device(this%west_in_buffer)
  call MPI_Put(this%west_in_buffer, msg_size, &
      var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
  !$acc end host_data

  !$acc end data

end subroutine

module subroutine retrieve_north_halo(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(in) :: this
  class(variable_t), intent(inout) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, nx, offs_x, offs_y, i, k, j
  integer :: its, ite, jts, jte, kts, kte
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag

  its = var%grid%its
  ite = var%grid%ite
  jts = var%grid%jts
  jte = var%grid%jte
  kts = var%grid%kts
  kte = var%grid%kte

  !$acc data present(this%north_in_3d)
if (var%two_d) then
        n = ubound(var%data_2d,2)
        nx = size(var%data_2d,1)
        !$acc parallel loop gang vector collapse(2) present(var%data_2d)
        do j = n-this%halo_size+1-offs_y, n
            do i = var%grid%its, var%grid%ite
                var%data_2d(i,j) = this%north_in_3d(i-var%grid%its+1+this%halo_size,1,j-(n-this%halo_size+1-offs_y)+1)
            enddo
        enddo
else
        n = ubound(var%data_3d,3)
        nx = size(var%data_3d,1)
        if (dqdt) then
                !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
                do j = n-this%halo_size+1-offs_y, n
                    do k = var%grid%kts, var%grid%kte
                        do i = var%grid%its, var%grid%ite
                            var%dqdt_3d(i,k,j) = this%north_in_3d(i-var%grid%its+1+this%halo_size,k,j-(n-this%halo_size+1-offs_y)+1)
                        enddo
                    enddo
                enddo
        else
                !$acc parallel loop gang vector collapse(3) present(var%data_3d)
                do j = n-this%halo_size+1-offs_y, n
                    do k = var%grid%kts, var%grid%kte
                        do i = var%grid%its, var%grid%ite
                            var%data_3d(i,k,j) = this%north_in_3d(i-var%grid%its+1+this%halo_size,k,j-(n-this%halo_size+1-offs_y)+1)
                        enddo
                    enddo
                enddo
        endif
endif

    !$acc end data
end subroutine

module subroutine retrieve_south_halo(this,var,do_dqdt)    
    implicit none
  class(halo_t), intent(in) :: this
  class(variable_t), intent(inout) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, nx, offs_y, offs_x, i, j, k
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs_y=var%ystag

  !$acc data present(this%south_in_3d)
if (var%two_d) then
    start = lbound(var%data_2d,2)
    !$acc parallel loop gang vector collapse(2) present(var%data_2d)
    do j = start,start+this%halo_size-1
    do i = var%grid%its,var%grid%ite
        var%data_2d(i,j) = this%south_in_3d(i-var%grid%its+1+this%halo_size,1,j-start+1)
    enddo
    enddo
else
    start = lbound(var%data_3d,3)
    if (dqdt) then
        !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
        do j = start,start+this%halo_size-1
        do k = var%grid%kts,var%grid%kte
        do i = var%grid%its,var%grid%ite
            var%dqdt_3d(i,k,j) = this%south_in_3d(i-var%grid%its+1+this%halo_size,k,j-start+1)
        enddo
        enddo
        enddo
    else
        !$acc parallel loop gang vector collapse(3) present(var%data_3d)
        do j = start,start+this%halo_size-1
        do k = var%grid%kts,var%grid%kte
        do i = var%grid%its,var%grid%ite
            var%data_3d(i,k,j) = this%south_in_3d(i-var%grid%its+1+this%halo_size,k,j-start+1)
        enddo
        enddo
        enddo
    endif
endif

    !$acc end data
end subroutine

module subroutine retrieve_east_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
integer :: n, ny, offs_x, offs_y, i, j, k
logical :: dqdt

dqdt=.False.
if (present(do_dqdt)) dqdt=do_dqdt

offs_x=var%xstag

!$acc data present(this%east_in_3d)
if (var%two_d) then
    n = ubound(var%data_2d,1)
    !$acc parallel loop gang vector collapse(2) present(var%data_2d)
    do j = var%grid%jts,var%grid%jte
    do i = n-this%halo_size+1-offs_x,n
        var%data_2d(i,j) = this%east_in_3d(i-(n-this%halo_size+1-offs_x)+1,1,j-var%grid%jts+1+this%halo_size)
    enddo
    enddo
else
    n = ubound(var%data_3d,1)
    if (dqdt) then
        !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
        do j = var%grid%jts,var%grid%jte
        do k = var%grid%kts,var%grid%kte
        do i = n-this%halo_size+1-offs_x,n
            var%dqdt_3d(i,k,j) = this%east_in_3d(i-(n-this%halo_size+1-offs_x)+1,k,j-var%grid%jts+1+this%halo_size)
        enddo
        enddo
        enddo
    else
        !$acc parallel loop gang vector collapse(3) present(var%data_3d)
        do j = var%grid%jts,var%grid%jte
        do k = var%grid%kts,var%grid%kte
        do i = n-this%halo_size+1-offs_x,n
            var%data_3d(i,k,j) = this%east_in_3d(i-(n-this%halo_size+1-offs_x)+1,k,j-var%grid%jts+1+this%halo_size)
        enddo
        enddo
        enddo
    endif
endif
    !$acc end data
end subroutine

module subroutine retrieve_west_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
integer :: start, ny, offs_x, offs_y, i, j, k
logical :: dqdt

dqdt=.False.
if (present(do_dqdt)) dqdt=do_dqdt

offs_x=var%xstag

!$acc data present(this%west_in_3d)
if (var%two_d) then
    start = lbound(var%data_2d,1)
    !$acc parallel loop gang vector collapse(2) present(var%data_2d)
    do j = var%grid%jts,var%grid%jte
    do i = start,start+this%halo_size-1
        var%data_2d(i,j) = this%west_in_3d(i-start+1,1,j-var%grid%jts+1+this%halo_size)
    enddo
    enddo
else
    start = lbound(var%data_3d,1)
    if (dqdt) then
        !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
        do j = var%grid%jts,var%grid%jte
        do k = var%grid%kts,var%grid%kte
        do i = start,start+this%halo_size-1
            var%dqdt_3d(i,k,j) = this%west_in_3d(i-start+1,k,j-var%grid%jts+1+this%halo_size)
        enddo
        enddo
        enddo
    else
        !$acc parallel loop gang vector collapse(3) present(var%data_3d)
        do j = var%grid%jts,var%grid%jte
        do k = var%grid%kts,var%grid%kte
        do i = start,start+this%halo_size-1
            var%data_3d(i,k,j) = this%west_in_3d(i-start+1,k,j-var%grid%jts+1+this%halo_size)
        enddo
        enddo
        enddo
    endif
endif
    !$acc end data
end subroutine



module subroutine put_northeast(this,var,do_dqdt)    
    implicit none
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

    !$acc data present(var%data_2d, var%data_3d, var%dqdt_3d)
    !$acc host_data use_device(var%data_2d, var%data_3d, var%dqdt_3d)
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
    !$acc end host_data
    !$acc end data
end subroutine

module subroutine put_northwest(this,var,do_dqdt)    
    implicit none
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

    !$acc data present(var%data_2d, var%data_3d, var%dqdt_3d)
    !$acc host_data use_device(var%data_2d, var%data_3d, var%dqdt_3d)
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
    !$acc end host_data
    !$acc end data
end subroutine


module subroutine put_southwest(this,var,do_dqdt)    
    implicit none
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

    !$acc data present(var%data_2d, var%data_3d, var%dqdt_3d)
    !$acc host_data use_device(var%data_2d, var%data_3d, var%dqdt_3d)
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
    !$acc end host_data
    !$acc end data
end subroutine

module subroutine put_southeast(this,var,do_dqdt)    
    implicit none
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

    !$acc data present(var%data_2d, var%data_3d, var%dqdt_3d)
    !$acc host_data use_device(var%data_2d, var%data_3d, var%dqdt_3d)
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
    !$acc end host_data
    !$acc end data
end subroutine


module subroutine retrieve_northeast_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y, n, nx, i, j, k
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
    !$acc data present(this%east_in_3d)
  if (var%two_d) then
        n = ubound(var%data_2d,2)
        nx = size(var%data_2d,1)
        !$acc parallel loop gang vector collapse(2) present(var%data_2d)
        do j = var%grid%jte+1-offs_y, var%grid%jme
        do i = var%grid%ite+1-offs_x, var%grid%ime
            var%data_2d(i,j) = this%east_in_3d(i-(var%grid%ite+1-offs_x)+1,1,j-(var%grid%jte+1-offs_y)+1)
        enddo
        enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jte+1-offs_y, var%grid%jme
          do k = this%kts, this%kte
          do i = var%grid%ite+1-offs_x, var%grid%ime
              var%dqdt_3d(i,k,j) = this%east_in_3d(i-(var%grid%ite+1-offs_x)+1,k,j-(var%grid%jte+1-offs_y)+1)
          enddo
          enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jte+1-offs_y, var%grid%jme
          do k = this%kts, this%kte
          do i = var%grid%ite+1-offs_x, var%grid%ime
              var%data_3d(i,k,j) = this%east_in_3d(i-(var%grid%ite+1-offs_x)+1,k,j-(var%grid%jte+1-offs_y)+1)
          enddo
          enddo
          enddo
      endif
  endif
    !$acc end data
end subroutine

module subroutine retrieve_northwest_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y, n, nx, i, j, k
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
    !$acc data present(this%north_in_3d)
  if (var%two_d) then
        !$acc parallel loop gang vector collapse(2) present(var%data_2d)
        do j = var%grid%jte+1-offs_y, var%grid%jme
        do i = var%grid%ims, var%grid%its-1
            var%data_2d(i,j) = this%north_in_3d(i-var%grid%ims+1,1,j-(var%grid%jte+1-offs_y)+1)
        enddo
        enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jte+1-offs_y, var%grid%jme
          do k = this%kts, this%kte
          do i = var%grid%ims, var%grid%its-1
              var%dqdt_3d(i,k,j) = this%north_in_3d(i-var%grid%ims+1,k,j-(var%grid%jte+1-offs_y)+1)
          enddo
          enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jte+1-offs_y, var%grid%jme
          do k = this%kts, this%kte
          do i = var%grid%ims, var%grid%its-1
              var%data_3d(i,k,j) = this%north_in_3d(i-var%grid%ims+1,k,j-(var%grid%jte+1-offs_y)+1)
          enddo
          enddo
          enddo
      endif
  endif
    !$acc end data
end subroutine

module subroutine retrieve_southwest_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
  integer :: i, j, k
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
    !$acc data present(this%west_in_3d)
  if (var%two_d) then
        !$acc parallel loop gang vector collapse(2) present(var%data_2d)
        do j = var%grid%jms, var%grid%jts-1
        do i = var%grid%ims, var%grid%its-1
            var%data_2d(i,j) = this%west_in_3d(i-var%grid%ims+1,1,j-var%grid%jms+1)
        enddo
        enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jms, var%grid%jts-1
          do k = this%kts, this%kte
          do i = var%grid%ims, var%grid%its-1
              var%dqdt_3d(i,k,j) = this%west_in_3d(i-var%grid%ims+1,k,j-var%grid%jms+1)
          enddo
          enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jms, var%grid%jts-1
          do k = this%kts, this%kte
          do i = var%grid%ims, var%grid%its-1
              var%data_3d(i,k,j) = this%west_in_3d(i-var%grid%ims+1,k,j-var%grid%jms+1)
          enddo
          enddo
          enddo
      endif
  endif
    !$acc end data
end subroutine

module subroutine retrieve_southeast_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
  integer :: offs_x, offs_y, i, j, k
  logical :: dqdt

  offs_x=var%xstag
  offs_y=var%ystag

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt

    !$acc data present(this%south_in_3d)
  if (var%two_d) then
        !$acc parallel loop gang vector collapse(2) present(var%data_2d)
        do j = var%grid%jms, var%grid%jts-1
        do i = var%grid%ite+1-offs_x, var%grid%ime
            var%data_2d(i,j) = this%south_in_3d(i-(var%grid%ite+1-offs_x)+1,1,j-var%grid%jms+1)
        enddo
        enddo
  else
      if (dqdt) then
          !$acc parallel loop gang vector collapse(3) present(var%dqdt_3d)
          do j = var%grid%jms, var%grid%jts-1
          do k = this%kts, this%kte
          do i = var%grid%ite+1-offs_x, var%grid%ime
              var%dqdt_3d(i,k,j) = this%south_in_3d(i-(var%grid%ite+1-offs_x)+1,k,j-var%grid%jms+1)
          enddo
          enddo
          enddo
      else
          !$acc parallel loop gang vector collapse(3) present(var%data_3d)
          do j = var%grid%jms, var%grid%jts-1
          do k = this%kts, this%kte
          do i = var%grid%ite+1-offs_x, var%grid%ime
              var%data_3d(i,k,j) = this%south_in_3d(i-(var%grid%ite+1-offs_x)+1,k,j-var%grid%jms+1)
          enddo
          enddo
          enddo
      endif
  endif
    !$acc end data

end subroutine

!> -------------------------------
!! Detect if neighbors are on shared memory hardware
!!
!! -------------------------------
subroutine detect_shared_memory(this, comms)    
    implicit none
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
    
#ifndef _OPENACC
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
#endif
    ! Free the shared memory communicator
    call MPI_Comm_free(shared_comm, ierr)
    
    ! Set whether to use shared windows based on if any neighbors are shared
    this%use_shared_windows = this%north_shared .or. this%south_shared .or. &
                             this%east_shared .or. this%west_shared
end subroutine detect_shared_memory

subroutine setup_batch_exch_north_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: north win, nx: ",nx
        write(*,*) "ERROR in setup batch_exch: north win, ny: ",ny
        write(*,*) "ERROR in setup batch_exch: north win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: north win, n_3d: ",this%n_3d
    endif
    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC

    allocate(this%north_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%north_batch_in_3d)
    !$acc host_data use_device(this%north_batch_in_3d)
    call MPI_WIN_CREATE(this%north_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%north_3d_win, ierr)
    !$acc end host_data

    if (this%n_2d > 0) then
        allocate(this%north_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%north_batch_in_2d)
        !$acc host_data use_device(this%north_batch_in_2d)
        call MPI_WIN_CREATE(this%north_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%north_2d_win, ierr)
        !$acc end host_data
    endif

#else
    if (this%north_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%north_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%north_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%north_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%north_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%north_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%north_batch_in_2d, [this%n_2d, nx, ny])

#endif
end subroutine setup_batch_exch_north_wins

subroutine setup_batch_exch_south_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: south win, nx: ",nx
        write(*,*) "ERROR in setup batch_exch: south win, ny: ",ny
        write(*,*) "ERROR in setup batch_exch: south win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: south win, n_3d: ",this%n_3d
    endif
    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC

    allocate(this%south_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%south_batch_in_3d)
    !$acc host_data use_device(this%south_batch_in_3d)
    call MPI_WIN_CREATE(this%south_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%south_3d_win, ierr)
    !$acc end host_data

    if (this%n_2d > 0) then
        allocate(this%south_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%south_batch_in_2d)
        !$acc host_data use_device(this%south_batch_in_2d)
        call MPI_WIN_CREATE(this%south_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%south_2d_win, ierr)
        !$acc end host_data
    endif

#else
    if (this%south_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%south_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%south_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%south_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%south_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%south_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%south_batch_in_2d, [this%n_2d, nx, ny])

#endif
end subroutine setup_batch_exch_south_wins

subroutine setup_batch_exch_east_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: east win, nx: ",nx
        write(*,*) "ERROR in setup batch_exch: east win, ny: ",ny
        write(*,*) "ERROR in setup batch_exch: east win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: east win, n_3d: ",this%n_3d
    endif 
    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%east_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%east_batch_in_3d)
    !$acc host_data use_device(this%east_batch_in_3d)
    call MPI_WIN_CREATE(this%east_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%east_3d_win, ierr)
    !$acc end host_data

    if (this%n_2d > 0) then
        allocate(this%east_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%east_batch_in_2d)
        !$acc host_data use_device(this%east_batch_in_2d)
        call MPI_WIN_CREATE(this%east_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%east_2d_win, ierr)
        !$acc end host_data
    endif
#else
    if (this%east_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%east_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%east_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%east_3d_win, ierr)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%east_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%east_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%east_batch_in_2d, [this%n_2d, nx, ny])
#endif

end subroutine setup_batch_exch_east_wins

subroutine setup_batch_exch_west_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: west win, nx: ",nx
        write(*,*) "ERROR in setup batch_exch: west win, ny: ",ny
        write(*,*) "ERROR in setup batch_exch: west win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: west win, n_3d: ",this%n_3d
    endif
    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%west_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%west_batch_in_3d)
    !$acc host_data use_device(this%west_batch_in_3d)
    call MPI_WIN_CREATE(this%west_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%west_3d_win, ierr)
    !$acc end host_data

    if (this%n_2d > 0) then
        allocate(this%west_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%west_batch_in_2d)
        !$acc host_data use_device(this%west_batch_in_2d)
        call MPI_WIN_CREATE(this%west_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%west_2d_win, ierr)
        !$acc end host_data
    endif
#else
    if (this%west_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%west_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE_SHARED(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%west_2d_win, ierr)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%west_3d_win)
        if (this%n_2d > 0) call MPI_WIN_ALLOCATE(win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr_2d, this%west_2d_win, ierr)
    endif
    call C_F_POINTER(tmp_ptr, this%west_batch_in_3d, [this%n_3d, nx, nz, ny])
    if (this%n_2d > 0) call C_F_POINTER(tmp_ptr_2d, this%west_batch_in_2d, [this%n_2d, nx, ny])
#endif
end subroutine setup_batch_exch_west_wins

subroutine setup_batch_exch_northwest_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: northwest win, halo_size: ",this%halo_size
        write(*,*) "ERROR in setup batch_exch: northwest win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: northwest win, n_3d: ",this%n_3d
    endif

    ! Create a group for the northwest-southwest exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages-1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%northwest_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%northwest_batch_in_3d)
    !$acc host_data use_device(this%northwest_batch_in_3d)
    call MPI_WIN_CREATE(this%northwest_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%northwest_3d_win, ierr)
    !$acc end host_data
#else
    if (this%northwest_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northwest_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northwest_3d_win)
    endif
    call C_F_POINTER(tmp_ptr, this%northwest_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])
#endif
end subroutine setup_batch_exch_northwest_wins

subroutine setup_batch_exch_northeast_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: northeast win, halo_size: ",this%halo_size
        write(*,*) "ERROR in setup batch_exch: northeast win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: northeast win, n_3d: ",this%n_3d
    endif

    ! Create a group for the northeast-southeast exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%halo_rank+this%grid%ximages+1/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%northeast_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%northeast_batch_in_3d)
    !$acc host_data use_device(this%northeast_batch_in_3d)
    call MPI_WIN_CREATE(this%northeast_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%northeast_3d_win, ierr)
    !$acc end host_data
#else
    if (this%northeast_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northeast_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%northeast_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%northeast_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])
#endif
end subroutine setup_batch_exch_northeast_wins

subroutine setup_batch_exch_southwest_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: southwest win, halo_size: ",this%halo_size
        write(*,*) "ERROR in setup batch_exch: southwest win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: southwest win, n_3d: ",this%n_3d
    endif

    ! Create a group for the southwest-southeast exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages-1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%southwest_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%southwest_batch_in_3d)
    !$acc host_data use_device(this%southwest_batch_in_3d)
    call MPI_WIN_CREATE(this%southwest_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%southwest_3d_win, ierr)
    !$acc end host_data
#else
    if (this%southwest_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southwest_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southwest_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%southwest_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])
#endif
end subroutine setup_batch_exch_southwest_wins

subroutine setup_batch_exch_southeast_wins(this, comms, info_in)    
    implicit none
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
    if (win_size <= 0) then
        write(*,*) "ERROR in setup batch_exch: southeast win, halo_size: ",this%halo_size
        write(*,*) "ERROR in setup batch_exch: southeast win, nz: ",nz
        write(*,*) "ERROR in setup batch_exch: southeast win, n_3d: ",this%n_3d
    endif

    ! Create a group for the southeast-southwest exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank-this%grid%ximages+1, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)

#ifdef _OPENACC
    allocate(this%southeast_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%southeast_batch_in_3d)
    !$acc host_data use_device(this%southeast_batch_in_3d)
    call MPI_WIN_CREATE(this%southeast_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%southeast_3d_win, ierr)
    !$acc end host_data
#else
    if (this%southeast_shared) then
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southeast_3d_win)
    else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, tmp_MPI_comm, tmp_ptr, this%southeast_3d_win)
    endif
    call C_F_POINTER(tmp_ptr,this%southeast_batch_in_3d,[this%n_3d,this%halo_size,nz,this%halo_size])
#endif
end subroutine setup_batch_exch_southeast_wins

end submodule
