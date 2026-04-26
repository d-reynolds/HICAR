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
#ifdef USE_NCCL
    use nccl_interface
    use openacc, only: acc_get_device_num, acc_device_nvidia, acc_get_cuda_stream, acc_async_sync
#endif
implicit none


contains


!> -------------------------------
!! Initialize the exchange arrays and dimensions
!!
!! -------------------------------
module subroutine init_halo(this, exch_vars, grid, comms)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(in) :: exch_vars(:)
    type(grid_t), intent(in) :: grid
    type(MPI_comm), intent(inout) :: comms

    type(MPI_Group) :: comp_proc, neighbor_group
    type(c_ptr) :: tmp_ptr
    type(MPI_Info) :: info_in
    integer(KIND=MPI_ADDRESS_KIND) :: win_size
    integer :: current, my_index, n_neighbors, nx, nz, ny, ierr
    integer :: i,j,k, real_size, nccl_nranks, nccl_ierr

    CALL MPI_Type_size(MPI_REAL, real_size)

    !Some stuff we can just copy right over
    this%grid = grid
    this%halo_size = grid%halo_size

    this%north_boundary = (this%grid%yimg == this%grid%yimages)
    this%south_boundary = (this%grid%yimg == 1)
    this%east_boundary  = (this%grid%ximg == this%grid%ximages)
    this%west_boundary  = (this%grid%ximg == 1)

    this%northwest_boundary = this%north_boundary .and. this%west_boundary
    this%southeast_boundary = this%south_boundary .and. this%east_boundary
    this%southwest_boundary = this%south_boundary .and. this%west_boundary
    this%northeast_boundary = this%north_boundary .and. this%east_boundary

    this%corner = this%northwest_boundary .or. this%northeast_boundary .or. &
                  this%southwest_boundary .or. this%southeast_boundary

    this%interior = .not.(this%north_boundary .or. this%south_boundary .or. &
                          this%east_boundary  .or. this%west_boundary)

    this%ims = this%grid%ims; this%its = this%grid%its; this%ids = this%grid%ids
    this%ime = this%grid%ime; this%ite = this%grid%ite; this%ide = this%grid%ide
    this%kms = this%grid%kms; this%kts = this%grid%kts; this%kds = this%grid%kds
    this%kme = this%grid%kme; this%kte = this%grid%kte; this%kde = this%grid%kde
    this%jms = this%grid%jms; this%jts = this%grid%jts; this%jds = this%grid%jds
    this%jme = this%grid%jme; this%jte = this%grid%jte; this%jde = this%grid%jde
    this%north_neighbor = -1; this%south_neighbor = -1; this%west_neighbor = -1; this%east_neighbor = -1
    this%northwest_neighbor = -1; this%southwest_neighbor = -1; this%northeast_neighbor = -1; this%southeast_neighbor = -1

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
        if (this%north_boundary) then
            this%northwest_neighbor = this%west_neighbor
        else if (this%west_boundary) then
            this%northwest_neighbor = this%north_neighbor
        else
            this%northwest_neighbor = this%halo_rank + this%grid%ximages - 1
        endif
        call MPI_Group_Incl(comp_proc, 1, [this%northwest_neighbor], this%northwest_neighbor_grp)
    endif
    if (.not.(this%southeast_boundary)) then
        if (this%south_boundary) then
            this%southeast_neighbor = this%east_neighbor
        else if (this%east_boundary) then
            this%southeast_neighbor = this%south_neighbor
        else
            this%southeast_neighbor = this%halo_rank - this%grid%ximages + 1
        endif
        call MPI_Group_Incl(comp_proc, 1, [this%southeast_neighbor], this%southeast_neighbor_grp)
    endif
    if (.not.(this%southwest_boundary)) then
        if (this%south_boundary) then
            this%southwest_neighbor = this%west_neighbor
        else if (this%west_boundary) then
            this%southwest_neighbor = this%south_neighbor
        else
            this%southwest_neighbor = this%halo_rank - this%grid%ximages - 1
        endif
        call MPI_Group_Incl(comp_proc, 1, [this%southwest_neighbor], this%southwest_neighbor_grp)
    endif
    if (.not.(this%northeast_boundary)) then
        if (this%north_boundary) then
            this%northeast_neighbor = this%east_neighbor
        else if (this%east_boundary) then
            this%northeast_neighbor = this%north_neighbor
        else
            this%northeast_neighbor = this%halo_rank + this%grid%ximages + 1
        endif
        call MPI_Group_Incl(comp_proc, 1, [this%northeast_neighbor], this%northeast_neighbor_grp)
    endif

#ifdef USE_NCCL
        call MPI_Comm_size(comms, nccl_nranks)
        nccl_ierr = nccl_comm_init(this%nccl_comm, nccl_nranks, this%halo_rank, comms%MPI_VAL, &
                                   int(acc_get_device_num(acc_device_nvidia), c_int))
        if (nccl_ierr /= 0) then
            write(*,*) "ERROR: NCCL communicator init failed on rank ", this%halo_rank
            call MPI_Abort(comms, nccl_ierr)
        endif
        ! Use OpenACC's synchronous-queue CUDA stream for NCCL operations.
        ! Sharing the stream provides natural within-stream ordering between
        ! pack/unpack OpenACC kernels and NCCL send/recv, eliminating the need
        ! for host-side stream synchronization and silencing racecheck false
        ! positives that arise from cross-stream host-side syncs.
        this%nccl_stream = transfer(acc_get_cuda_stream(acc_async_sync), this%nccl_stream)
        if (STD_OUT_PE) write(*,*) "NCCL communicator initialized for halo exchange"
#endif

    ! Detect if neighbors are on shared memory hardware
    call detect_shared_memory(this, comms)

    !Now allocate the actual 3D halo
    !We only want to set up remote windows for domain objects which are part of the actual domain
    if (.not.(comms == MPI_COMM_NULL)) then

#ifndef USE_NCCL
        call MPI_Info_create(info_in, ierr)
#ifdef _OPENACC
        ! Check if MPI supports CUDA-aware memory
        call MPI_Info_set(info_in, "alloc_shm", ".true.", ierr)
        ! Set memory type hint for GPU accessibility
        call MPI_Info_set(info_in, "alloc_mem", "device", ierr)
        call MPI_Info_set(info_in, "mpi_assert_memory_alloc_kinds", "gpu:device", ierr)
        call MPI_Info_set(info_in, "cuda_aware", ".true.", ierr)
#else
        call MPI_INFO_SET(info_in, 'no_locks', '.true.')
        call MPI_INFO_SET(info_in, 'same_size', '.true.')
        call MPI_INFO_SET(info_in, 'same_disp_unit', '.true.')
#endif
#endif
        nx = this%grid%ns_halo_nx
        nz = this%grid%halo_win_nz   ! single-var exch_var windows sized to max
        ny = this%halo_size+1
        win_size = nx*nz*ny

#ifdef _OPENACC
        allocate(this%south_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%south_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%south_in_3d)
        call MPI_WIN_CREATE(this%south_in_3d, win_size*real_size, real_size, info_in, comms, this%south_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%south_in_win)
        call C_F_POINTER(tmp_ptr, this%south_in_3d, [nx, nz, ny])
#endif
        allocate(this%south_in_buffer(1:this%grid%ns_halo_nx, 1:nz, 1:this%halo_size+1))
        allocate(this%south_in_buffer_2d(1:this%grid%ns_halo_nx, 1:this%halo_size+1))
        !$acc enter data copyin(this%south_in_buffer_2d, this%south_in_buffer)
#ifdef _OPENACC
        allocate(this%north_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%north_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%north_in_3d)
        call MPI_WIN_CREATE(this%north_in_3d, win_size*real_size, real_size, info_in, comms, this%north_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%north_in_win)
        call C_F_POINTER(tmp_ptr, this%north_in_3d, [nx, nz, ny])
#endif
        allocate(this%north_in_buffer(1:this%grid%ns_halo_nx, 1:nz, 1:this%halo_size+1))
        allocate(this%north_in_buffer_2d(1:this%grid%ns_halo_nx, 1:this%halo_size+1))
        !$acc enter data copyin(this%north_in_buffer_2d, this%north_in_buffer)

        nx = this%halo_size+1
        nz = this%grid%halo_win_nz   ! single-var exch_var windows sized to max
        ny = this%grid%ew_halo_ny
        win_size = nx*nz*ny

#ifdef _OPENACC
        allocate(this%east_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%east_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%east_in_3d)
        call MPI_WIN_CREATE(this%east_in_3d, win_size*real_size, real_size, info_in, comms, this%east_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%east_in_win)
        call C_F_POINTER(tmp_ptr, this%east_in_3d, [nx, nz, ny])
#endif
        allocate(this%east_in_buffer(1:this%halo_size+1, 1:nz, 1:this%grid%ew_halo_ny))
        allocate(this%east_in_buffer_2d(1:this%halo_size+1, 1:this%grid%ew_halo_ny))
        !$acc enter data copyin(this%east_in_buffer_2d, this%east_in_buffer)

#ifdef _OPENACC
        allocate(this%west_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%west_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%west_in_3d)
        call MPI_WIN_CREATE(this%west_in_3d, win_size*real_size, real_size, info_in, comms, this%west_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%west_in_win)
        call C_F_POINTER(tmp_ptr, this%west_in_3d, [nx, nz, ny])
#endif
        allocate(this%west_in_buffer(1:this%halo_size+1, 1:nz, 1:this%grid%ew_halo_ny))
        allocate(this%west_in_buffer_2d(1:this%halo_size+1, 1:this%grid%ew_halo_ny))
        !$acc enter data copyin(this%west_in_buffer_2d, this%west_in_buffer)
        this%north_in_3d = 1
        this%south_in_3d = 1
        this%east_in_3d = 1
        this%west_in_3d = 1

        nx = this%halo_size+1
        nz = this%grid%halo_win_nz   ! single-var exch_var windows sized to max
        ny = this%halo_size+1
        win_size = nx*nz*ny

#ifdef _OPENACC
        allocate(this%southwest_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%southwest_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%southwest_in_3d)
        call MPI_WIN_CREATE(this%southwest_in_3d, win_size*real_size, real_size, info_in, comms, this%southwest_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%southwest_in_win)
        call C_F_POINTER(tmp_ptr, this%southwest_in_3d, [nx, nz, ny])
#endif

#ifdef _OPENACC
        allocate(this%northwest_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%northwest_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%northwest_in_3d)
        call MPI_WIN_CREATE(this%northwest_in_3d, win_size*real_size, real_size, info_in, comms, this%northwest_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%northwest_in_win)
        call C_F_POINTER(tmp_ptr, this%northwest_in_3d, [nx, nz, ny])
#endif

#ifdef _OPENACC
        allocate(this%northeast_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%northeast_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%northeast_in_3d)
        call MPI_WIN_CREATE(this%northeast_in_3d, win_size*real_size, real_size, info_in, comms, this%northeast_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%northeast_in_win)
        call C_F_POINTER(tmp_ptr, this%northeast_in_3d, [nx, nz, ny])
#endif

#ifdef _OPENACC
        allocate(this%southeast_in_3d(nx, nz, ny))
        !$acc enter data copyin(this%southeast_in_3d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%southeast_in_3d)
        call MPI_WIN_CREATE(this%southeast_in_3d, win_size*real_size, real_size, info_in, comms, this%southeast_in_win, ierr)
        !$acc end host_data
#endif
#else
        call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%southeast_in_win)
        call C_F_POINTER(tmp_ptr, this%southeast_in_3d, [nx, nz, ny])
#endif

        this%northwest_in_3d = 1
        this%southwest_in_3d = 1
        this%northeast_in_3d = 1
        this%southeast_in_3d = 1

        allocate(this%ne_corner_send(this%halo_size+1, nz, this%halo_size+1))
        allocate(this%nw_corner_send(this%halo_size+1, nz, this%halo_size+1))
        allocate(this%se_corner_send(this%halo_size+1, nz, this%halo_size+1))
        allocate(this%sw_corner_send(this%halo_size+1, nz, this%halo_size+1))
        this%ne_corner_send = 0
        this%nw_corner_send = 0
        this%se_corner_send = 0
        this%sw_corner_send = 0
        !$acc enter data copyin(this%ne_corner_send, this%nw_corner_send, this%se_corner_send, this%sw_corner_send)

    endif

    !...and the larger 3D halo for batch exchanges
    call setup_batch_exch(this, exch_vars, comms)

end subroutine init_halo

module subroutine finalize(this)
    implicit none
    class(halo_t), intent(inout) :: this

    integer :: ierr

#ifdef USE_NCCL
        ! Note: do NOT destroy this%nccl_stream — it's owned by OpenACC, not us.
        call nccl_comm_destroy(this%nccl_comm)
#else
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

    if (.not.(this%northwest_boundary)) call MPI_Win_Start(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Start(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Start(this%northeast_neighbor_grp, 0, this%northeast_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Start(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Complete( this%northwest_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Complete(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Complete(this%northeast_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Complete( this%southeast_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Wait( this%northwest_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Wait(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Wait(this%northeast_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Wait( this%southeast_3d_win)

    if (.not.(this%north_boundary)) call MPI_WIN_FREE(this%north_3d_win, ierr)
    if (.not.(this%south_boundary)) call MPI_WIN_FREE(this%south_3d_win, ierr)
    if (.not.(this%east_boundary)) call MPI_WIN_FREE(this%east_3d_win, ierr)
    if (.not.(this%west_boundary)) call MPI_WIN_FREE(this%west_3d_win, ierr)
    if (.not.(this%northwest_boundary)) call MPI_WIN_FREE(this%northwest_3d_win, ierr)
    if (.not.(this%southwest_boundary)) call MPI_WIN_FREE(this%southwest_3d_win, ierr)
    if (.not.(this%northeast_boundary)) call MPI_WIN_FREE(this%northeast_3d_win, ierr)
    if (.not.(this%southeast_boundary)) call MPI_WIN_FREE(this%southeast_3d_win, ierr)

    call MPI_Type_free(this%NS_3d_win_halo_type, ierr)
    call MPI_Type_free(this%EW_3d_win_halo_type, ierr)
    call MPI_Type_free(this%corner_3d_win_halo_type, ierr)
#endif

#ifndef USE_NCCL
    call MPI_Win_fence(0, this%north_in_win)
    call MPI_Win_fence(0, this%south_in_win)
    call MPI_Win_fence(0, this%east_in_win)
    call MPI_Win_fence(0, this%west_in_win)
    call MPI_Win_fence(0, this%southwest_in_win)
    call MPI_Win_fence(0, this%northwest_in_win)
    call MPI_Win_fence(0, this%southeast_in_win)
    call MPI_Win_fence(0, this%northeast_in_win)

    call MPI_WIN_FREE(this%north_in_win, ierr)
    call MPI_WIN_FREE(this%south_in_win, ierr)
    call MPI_WIN_FREE(this%east_in_win, ierr)
    call MPI_WIN_FREE(this%west_in_win, ierr)
    call MPI_WIN_FREE(this%southwest_in_win, ierr)
    call MPI_WIN_FREE(this%northwest_in_win, ierr)
    call MPI_WIN_FREE(this%northeast_in_win, ierr)
    call MPI_WIN_FREE(this%southeast_in_win, ierr)
#endif

    !$acc exit data finalize delete(this%north_in_3d, this%south_in_3d, this%east_in_3d, this%west_in_3d, &
    !$acc                      this%north_in_buffer, this%south_in_buffer, this%east_in_buffer, this%west_in_buffer, &
    !$acc                      this%south_buffer_2d, this%north_buffer_2d, this%east_buffer_2d, this%west_buffer_2d, &
    !$acc                      this%south_buffer_3d, this%north_buffer_3d, this%east_buffer_3d, this%west_buffer_3d, &
    !$acc                      this%northwest_buffer_3d, this%southwest_buffer_3d, this%northeast_buffer_3d, this%southeast_buffer_3d, &
    !$acc                      this%ne_corner_send, this%nw_corner_send, this%se_corner_send, this%sw_corner_send, &
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
#ifdef USE_NCCL
    integer(c_int) :: ns_count, ew_count, corner_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    do_corners=.False.
    if (present(corners)) do_corners=corners

#ifdef USE_NCCL
    ! NCCL pure path: pack-and-send each direction, post matching recvs into
    ! the corresponding *_in_3d, all fused into one group. NCCL ops share
    ! OpenACC's sync stream, so within-stream ordering serializes the
    ! preceding pack kernels before NCCL reads, and serializes NCCL writes
    ! before the unpack kernels in retrieve_*_halo.
    ns_count     = int(size(this%north_in_buffer), c_int)
    ew_count     = int(size(this%east_in_buffer),  c_int)
    corner_count = int(size(this%ne_corner_send),  c_int)

    call nccl_group_start()

    if (do_corners) then
        if (.not. this%northeast_boundary) call this%put_northeast(var, dqdt)
        if (.not. this%northwest_boundary) call this%put_northwest(var, dqdt)
        if (.not. this%southeast_boundary) call this%put_southeast(var, dqdt)
        if (.not. this%southwest_boundary) call this%put_southwest(var, dqdt)
    endif

    if (.not. this%north_boundary) call this%put_north(var, dqdt)
    if (.not. this%south_boundary) call this%put_south(var, dqdt)
    if (.not. this%east_boundary)  call this%put_east(var, dqdt)
    if (.not. this%west_boundary)  call this%put_west(var, dqdt)

    ! Post matching recvs. Each recv lands at offset 0 of *_in_3d; the
    ! retrieve_*_halo kernels read with halo_off=0 under USE_NCCL.
    if (.not. this%north_boundary) then
        !$acc host_data use_device(this%north_in_3d)
        nccl_ierr = nccl_recv_float(c_loc(this%north_in_3d), ns_count, &
            this%north_neighbor, this%nccl_comm, this%nccl_stream)
        !$acc end host_data
    endif
    if (.not. this%south_boundary) then
        !$acc host_data use_device(this%south_in_3d)
        nccl_ierr = nccl_recv_float(c_loc(this%south_in_3d), ns_count, &
            this%south_neighbor, this%nccl_comm, this%nccl_stream)
        !$acc end host_data
    endif
    if (.not. this%east_boundary) then
        !$acc host_data use_device(this%east_in_3d)
        nccl_ierr = nccl_recv_float(c_loc(this%east_in_3d), ew_count, &
            this%east_neighbor, this%nccl_comm, this%nccl_stream)
        !$acc end host_data
    endif
    if (.not. this%west_boundary) then
        !$acc host_data use_device(this%west_in_3d)
        nccl_ierr = nccl_recv_float(c_loc(this%west_in_3d), ew_count, &
            this%west_neighbor, this%nccl_comm, this%nccl_stream)
        !$acc end host_data
    endif

    if (do_corners) then
        if (.not. this%northeast_boundary) then
            !$acc host_data use_device(this%northeast_in_3d)
            nccl_ierr = nccl_recv_float(c_loc(this%northeast_in_3d), corner_count, &
                this%northeast_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif
        if (.not. this%northwest_boundary) then
            !$acc host_data use_device(this%northwest_in_3d)
            nccl_ierr = nccl_recv_float(c_loc(this%northwest_in_3d), corner_count, &
                this%northwest_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif
        if (.not. this%southeast_boundary) then
            !$acc host_data use_device(this%southeast_in_3d)
            nccl_ierr = nccl_recv_float(c_loc(this%southeast_in_3d), corner_count, &
                this%southeast_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif
        if (.not. this%southwest_boundary) then
            !$acc host_data use_device(this%southwest_in_3d)
            nccl_ierr = nccl_recv_float(c_loc(this%southwest_in_3d), corner_count, &
                this%southwest_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif
    endif

    call nccl_group_end()

    if (do_corners) then
        if (.not. this%northeast_boundary) call this%retrieve_northeast_halo(var, dqdt)
        if (.not. this%northwest_boundary) call this%retrieve_northwest_halo(var, dqdt)
        if (.not. this%southeast_boundary) call this%retrieve_southeast_halo(var, dqdt)
        if (.not. this%southwest_boundary) call this%retrieve_southwest_halo(var, dqdt)
    endif

    if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
    if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)
    if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
    if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)
#else

    if (do_corners) then

        call MPI_Win_fence(0,this%southwest_in_win)
        call MPI_Win_fence(0,this%northwest_in_win)
        call MPI_Win_fence(0,this%southeast_in_win)
        call MPI_Win_fence(0,this%northeast_in_win)

        if (.not. this%northeast_boundary) call this%put_northeast(var, dqdt)
        if (.not. this%northwest_boundary) call this%put_northwest(var, dqdt)
        if (.not. this%southeast_boundary)  call this%put_southeast(var, dqdt)
        if (.not. this%southwest_boundary)  call this%put_southwest(var, dqdt)


        call MPI_Win_fence(0,this%southwest_in_win)
        call MPI_Win_fence(0,this%northwest_in_win)
        call MPI_Win_fence(0,this%southeast_in_win)
        call MPI_Win_fence(0,this%northeast_in_win)

        if (.not. this%northeast_boundary) call this%retrieve_northeast_halo(var, dqdt)
        if (.not. this%northwest_boundary) call this%retrieve_northwest_halo(var, dqdt)
        if (.not. this%southeast_boundary)  call this%retrieve_southeast_halo(var, dqdt)
        if (.not. this%southwest_boundary)  call this%retrieve_southwest_halo(var, dqdt)
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
#endif

end subroutine exch_var


!> -------------------------------
!! Initialize the arrays and co-arrays needed to perform a batch exchange
!!
!! -------------------------------
subroutine setup_batch_exch(this, exch_vars, comms)
    implicit none
    type(halo_t), intent(inout) :: this
    type(index_type), intent(in) :: exch_vars(:)
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

    ! Loop over all exch vars and count how many are 3D
    do i = 1,size(exch_vars)
        var = get_varmeta(exch_vars(i)%id)
        if (var%three_d) then
            this%n_3d = this%n_3d + 1
        end if
    end do

    ! Determine number of 2D and 3D vars present
    this%n_2d = (size(exch_vars))-this%n_3d

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

#ifndef USE_NCCL
        call MPI_Type_contiguous(nx*ny*this%n_2d, MPI_REAL, this%NS_2d_win_halo_type)
        call MPI_Type_commit(this%NS_2d_win_halo_type)

        call MPI_Type_contiguous(nx*nz*ny*this%n_3d, MPI_REAL, this%NS_3d_win_halo_type)
        call MPI_Type_commit(this%NS_3d_win_halo_type)
#endif

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

#ifndef USE_NCCL
            if (this%south_shared) then
                call MPI_WIN_SHARED_QUERY(this%south_3d_win, 0, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%south_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then 
                    call MPI_WIN_SHARED_QUERY(this%south_2d_win, 0, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%south_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
#endif
                allocate(this%south_buffer_3d(this%n_3d,1:this%grid%ns_halo_nx,this%kms:this%kme,1:this%halo_size))
                if (this%n_2d > 0) allocate(this%south_buffer_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
            !$acc enter data copyin(this%south_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%south_buffer_2d)
            endif
        endif
        if (.not.(this%north_boundary)) then
            this%north_batch_in_3d = 1
            if (this%n_2d > 0) this%north_batch_in_2d = 1

#ifndef USE_NCCL
            if (this%north_shared) then
                call MPI_WIN_SHARED_QUERY(this%north_3d_win, 1, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%north_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%north_2d_win, 1, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%north_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
#endif
                allocate(this%north_buffer_3d(this%n_3d,1:this%grid%ns_halo_nx,this%kms:this%kme,1:this%halo_size))
                if (this%n_2d > 0) allocate(this%north_buffer_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
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

#ifndef USE_NCCL
        call MPI_Type_contiguous(nx*ny*this%n_2d, MPI_REAL, this%EW_2d_win_halo_type)
        call MPI_Type_commit(this%EW_2d_win_halo_type)

        call MPI_Type_contiguous(nx*nz*ny*this%n_3d, MPI_REAL, this%EW_3d_win_halo_type)
        call MPI_Type_commit(this%EW_3d_win_halo_type)
#endif

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

#ifndef USE_NCCL
            if (this%east_shared) then

                call MPI_WIN_SHARED_QUERY(this%east_3d_win, 1, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%east_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%east_2d_win, 1, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%east_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
#endif
                allocate(this%east_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%grid%ew_halo_ny))
                if (this%n_2d > 0) allocate(this%east_buffer_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny))
#ifndef USE_NCCL
            endif
#endif
            !$acc enter data copyin(this%east_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%east_buffer_2d)
            endif
        endif
        if (.not.(this%west_boundary)) then
            this%west_batch_in_3d = 1
            if (this%n_2d > 0) this%west_batch_in_2d = 1

#ifndef USE_NCCL
            if (this%west_shared) then
                call MPI_WIN_SHARED_QUERY(this%west_3d_win, 0, win_size, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%west_buffer_3d, [this%n_3d, nx, nz, ny])

                if (this%n_2d > 0) then
                    call MPI_WIN_SHARED_QUERY(this%west_2d_win, 0, win_size_2d, size_out, tmp_ptr)
                    call C_F_POINTER(tmp_ptr, this%west_buffer_2d, [this%n_2d, nx, ny])
                endif
            else
#endif
                allocate(this%west_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%grid%ew_halo_ny))
                if (this%n_2d > 0) allocate(this%west_buffer_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny))
#ifndef USE_NCCL
            endif
#endif
            !$acc enter data copyin(this%west_buffer_3d)
            if (this%n_2d > 0) then
                !$acc enter data copyin(this%west_buffer_2d)
            endif
        endif


        !Then do corners
        win_size_corner = this%halo_size*nz*this%halo_size*this%n_3d

#ifndef USE_NCCL
        call MPI_Type_contiguous(this%halo_size*nz*this%halo_size*this%n_3d, MPI_REAL, this%corner_3d_win_halo_type)
        call MPI_Type_commit(this%corner_3d_win_halo_type)
#endif

        ! Setup for corner exchanges. This will exclude edge processes which have a corner neighbor
        ! off-grid
        if ((mod(this%grid%ximg,2) == 1)) then
            if (.not.(this%north_boundary .or. this%west_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
            if (.not.(this%south_boundary .or. this%east_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
        else if((mod(this%grid%ximg,2) == 0)) then
            if (.not.(this%south_boundary .or. this%east_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
            if (.not.(this%north_boundary .or. this%west_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
        endif

        if ((mod(this%grid%ximg,2) == 1)) then
            if (.not.(this%north_boundary .or. this%east_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
            if (.not.(this%south_boundary .or. this%west_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
        else if((mod(this%grid%ximg,2) == 0)) then
            if (.not.(this%south_boundary .or. this%west_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
            if (.not.(this%north_boundary .or. this%east_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
        endif

        ! Now we just need to hook up the edge processes which were excluded above.
        ! The setup_batch_exch routines will conntect adjacent processes when they do not
        ! have a direct corner neighbor (i.e., they are on the edge of the grid)
        if (this%west_boundary) then
            if ((mod(this%grid%yimg,2) == 1)) then
                if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
                if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
            else if((mod(this%grid%yimg,2) == 0)) then
                if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
                if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
            endif
        endif
        if (this%east_boundary) then
            if ((mod(this%grid%yimg,2) == 1)) then
                if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
                if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
            else if((mod(this%grid%yimg,2) == 0)) then
                if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
                if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
            endif
        endif
        if (this%north_boundary) then
            if ((mod(this%grid%ximg,2) == 1)) then
                if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
                if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
            else if((mod(this%grid%ximg,2) == 0)) then
                if (.not.(this%northeast_boundary)) call setup_batch_exch_northeast_wins(this, comms, info_in)
                if (.not.(this%northwest_boundary)) call setup_batch_exch_northwest_wins(this, comms, info_in)
            endif
        endif
        if (this%south_boundary) then
            if ((mod(this%grid%ximg,2) == 1)) then
                if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
                if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
            else if((mod(this%grid%ximg,2) == 0)) then
                if (.not.(this%southeast_boundary)) call setup_batch_exch_southeast_wins(this, comms, info_in)
                if (.not.(this%southwest_boundary)) call setup_batch_exch_southwest_wins(this, comms, info_in)
            endif
        endif

        if (.not.(this%northwest_boundary)) then
#ifndef USE_NCCL
            if (this%northwest_shared) then
                call MPI_WIN_SHARED_QUERY(this%northwest_3d_win, merge(0, 1, this%halo_rank > this%northwest_neighbor), win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%northwest_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
#endif
                allocate(this%northwest_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
            this%northwest_batch_in_3d = 1
            !$acc enter data copyin(this%northwest_buffer_3d)
        endif
        if (.not.(this%southeast_boundary)) then
#ifndef USE_NCCL
            if (this%southeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%southeast_3d_win, merge(0, 1, this%halo_rank > this%southeast_neighbor), win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
#endif
                allocate(this%southeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
            this%southeast_batch_in_3d = 1
            !$acc enter data copyin(this%southeast_buffer_3d)
        endif
        if (.not.(this%southwest_boundary)) then
#ifndef USE_NCCL
            if (this%southwest_shared) then
                call MPI_WIN_SHARED_QUERY(this%southwest_3d_win, merge(0, 1, this%halo_rank > this%southwest_neighbor), win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%southwest_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
#endif
                allocate(this%southwest_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
            this%southwest_batch_in_3d = 1
            !$acc enter data copyin(this%southwest_buffer_3d)
        endif
        if (.not.(this%northeast_boundary)) then
#ifndef USE_NCCL
            if (this%northeast_shared) then
                call MPI_WIN_SHARED_QUERY(this%northeast_3d_win, merge(0, 1, this%halo_rank > this%northeast_neighbor), win_size_corner, size_out, tmp_ptr)
                call C_F_POINTER(tmp_ptr, this%northeast_buffer_3d, [this%n_3d, this%halo_size, nz, this%halo_size])
            else
#endif
                allocate(this%northeast_buffer_3d(this%n_3d,1:this%halo_size,this%kms:this%kme,1:this%halo_size))
#ifndef USE_NCCL
            endif
#endif
            this%northeast_batch_in_3d = 1
            !$acc enter data copyin(this%northeast_buffer_3d)
        endif

#ifndef USE_NCCL
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
#endif
    endif

end subroutine setup_batch_exch


module subroutine halo_3d_send_batch(this, vars_to_send, var_data)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: vars_to_send(:)
    type(variable_t), intent(inout) :: var_data(:)
    
    integer :: n, p, k_max, msg_size, indx, i, j, k, n_vars
    integer :: kms, kme, its, ite, jts, jte, halo_size
    integer :: i_start, j_start, kms_var, kme_var
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
    integer :: ierr, ns_msg_size, ew_msg_size, corner_msg_size

    if (this%n_3d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    msg_size = 1
    disp = 0

    halo_size = this%halo_size
    kms = this%kms; kme = this%kme
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

    n = 1
    n_vars = size(var_data)

#ifndef USE_NCCL
    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Start(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Start(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Start(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Start(this%northeast_neighbor_grp, 0, this%northeast_3d_win)
#endif

    ! Now iterate through the dictionary as long as there are more elements present. If two processors are on shared memory
    ! this step will directly copy the data to the other PE

    do p = 1, size(vars_to_send)
        if (vars_to_send(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(vars_to_send(p)%v)%id == vars_to_send(p)%id) then
                associate(var => var_data(vars_to_send(p)%v)%data_3d)
                kms_var = var_data(vars_to_send(p)%v)%grid%kts
                kme_var = var_data(vars_to_send(p)%v)%grid%kte

                if (.not.(this%north_boundary)) then
                    associate(buff => this%north_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
                    do i = 1,halo_size
                        buff(n,i,k,j-jts+1) = var(its+i-1,k,j)
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%northwest_boundary)) then
                    j_start = jte-halo_size
                    i_start = its-1 
                    if (this%north_boundary) then
                        j_start = j_start + halo_size
                    elseif (this%west_boundary) then
                        i_start = i_start - halo_size
                    endif
                    associate(buff => this%northwest_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
                    do i = 1,halo_size
                        buff(n,i,k,j) = var((i_start+i),k,(j_start+j))
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%southeast_boundary)) then
                    j_start = jts-1
                    i_start = ite-halo_size
                    if (this%south_boundary) then
                        j_start = j_start - halo_size
                    elseif (this%east_boundary) then
                        i_start = i_start + halo_size
                    endif
                    associate(buff => this%southeast_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
                    do i = 1,halo_size
                        buff(n,i,k,j) = var((i_start+i),k,(j_start+j))
                    enddo
                    enddo
                    enddo
                    end associate
                endif
                
                if (.not.(this%southwest_boundary)) then
                    j_start = jts-1
                    i_start = its-1
                    if (this%south_boundary) then
                        j_start = j_start - halo_size
                    elseif (this%west_boundary) then
                        i_start = i_start - halo_size
                    endif
                    associate(buff => this%southwest_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
                    do i = 1,halo_size
                        buff(n,i,k,j) = var((i_start+i),k,(j_start+j))
                    enddo
                    enddo
                    enddo
                    end associate
                endif

                if (.not.(this%northeast_boundary)) then
                    j_start = jte-halo_size
                    i_start = ite-halo_size
                    if (this%north_boundary) then
                        j_start = j_start + halo_size
                    elseif (this%east_boundary) then
                        i_start = i_start + halo_size
                    endif
                    associate(buff => this%northeast_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
                    do i = 1,halo_size
                        buff(n,i,k,j) = var((i_start+i),k,(j_start+j))
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

#ifdef USE_NCCL

        ns_msg_size = this%n_3d * this%grid%ns_halo_nx * this%grid%halo_nz * this%halo_size
        ew_msg_size = this%n_3d * this%halo_size * this%grid%halo_nz * this%grid%ew_halo_ny
        corner_msg_size = this%n_3d * this%halo_size * this%grid%halo_nz * this%halo_size

        ! NCCL operations share OpenACC's stream, so within-stream ordering
        ! automatically serializes the preceding pack kernels before NCCL reads.
        call nccl_group_start()

        if (.not. this%north_boundary) then
            !$acc host_data use_device(this%north_buffer_3d, this%north_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%north_buffer_3d), ns_msg_size, &
                this%north_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%north_batch_in_3d), ns_msg_size, &
                this%north_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%south_boundary) then
            !$acc host_data use_device(this%south_buffer_3d, this%south_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%south_buffer_3d), ns_msg_size, &
                this%south_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%south_batch_in_3d), ns_msg_size, &
                this%south_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%east_boundary) then
            !$acc host_data use_device(this%east_buffer_3d, this%east_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%east_buffer_3d), ew_msg_size, &
                this%east_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%east_batch_in_3d), ew_msg_size, &
                this%east_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%west_boundary) then
            !$acc host_data use_device(this%west_buffer_3d, this%west_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%west_buffer_3d), ew_msg_size, &
                this%west_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%west_batch_in_3d), ew_msg_size, &
                this%west_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%northwest_boundary) then
            !$acc host_data use_device(this%northwest_buffer_3d, this%northwest_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%northwest_buffer_3d), corner_msg_size, &
                this%northwest_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%northwest_batch_in_3d), corner_msg_size, &
                this%northwest_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%northeast_boundary) then
            !$acc host_data use_device(this%northeast_buffer_3d, this%northeast_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%northeast_buffer_3d), corner_msg_size, &
                this%northeast_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%northeast_batch_in_3d), corner_msg_size, &
                this%northeast_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%southwest_boundary) then
            !$acc host_data use_device(this%southwest_buffer_3d, this%southwest_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%southwest_buffer_3d), corner_msg_size, &
                this%southwest_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%southwest_batch_in_3d), corner_msg_size, &
                this%southwest_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%southeast_boundary) then
            !$acc host_data use_device(this%southeast_buffer_3d, this%southeast_batch_in_3d)
            ierr = nccl_send_float(c_loc(this%southeast_buffer_3d), corner_msg_size, &
                this%southeast_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%southeast_batch_in_3d), corner_msg_size, &
                this%southeast_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        call nccl_group_end()
#else
    !$acc data present(this)

    !$acc host_data use_device(this%south_buffer_3d, this%north_buffer_3d, &
    !$acc this%east_buffer_3d, this%west_buffer_3d, &
    !$acc this%northwest_buffer_3d, this%southeast_buffer_3d, &
    !$acc this%southwest_buffer_3d, this%northeast_buffer_3d)
    if (.not.(this%south_boundary)) then
        if (.not.(this%south_shared)) then
            call MPI_Put(this%south_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 0, disp, msg_size, &
                this%NS_3d_win_halo_type, this%south_3d_win)
        endif
    endif

    if (.not.(this%north_boundary)) then
        if (.not.(this%north_shared)) then
            call MPI_Put(this%north_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 1, disp, msg_size, &
                this%NS_3d_win_halo_type, this%north_3d_win)
        endif
    endif

    if (.not.(this%east_boundary)) then
        if (.not.(this%east_shared)) then
            call MPI_Put(this%east_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 1, disp, msg_size, &
                this%EW_3d_win_halo_type, this%east_3d_win)
        endif
    endif

    if (.not.(this%west_boundary)) then
        if (.not.(this%west_shared)) then
            call MPI_Put(this%west_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 0, disp, msg_size, &
                this%EW_3d_win_halo_type, this%west_3d_win)
        endif
    endif

    if (.not.(this%northwest_boundary)) then
        if (.not.(this%northwest_shared)) then
            call MPI_Put(this%northwest_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, &
                merge(0, 1, this%halo_rank > this%northwest_neighbor), &
                disp, msg_size, this%corner_3d_win_halo_type, this%northwest_3d_win)
        endif
    endif
    if (.not.(this%southeast_boundary)) then
        if (.not.(this%southeast_shared)) then
            call MPI_Put(this%southeast_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, &
                merge(0, 1, this%halo_rank > this%southeast_neighbor), &
                disp, msg_size, this%corner_3d_win_halo_type, this%southeast_3d_win)
        endif
    endif
    if (.not.(this%southwest_boundary)) then
        if (.not.(this%southwest_shared)) then
            call MPI_Put(this%southwest_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, &
                merge(0, 1, this%halo_rank > this%southwest_neighbor), &
                disp, msg_size, this%corner_3d_win_halo_type, this%southwest_3d_win)
        endif
    endif
    if (.not.(this%northeast_boundary)) then
        if (.not.(this%northeast_shared)) then
            call MPI_Put(this%northeast_buffer_3d, msg_size, &
                this%corner_3d_win_halo_type, &
                merge(0, 1, this%halo_rank > this%northeast_neighbor), &
                disp, msg_size, this%corner_3d_win_halo_type, this%northeast_3d_win)
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
#endif

end subroutine halo_3d_send_batch

module subroutine halo_3d_retrieve_batch(this, vars_to_ret, var_data, wait_timer)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: vars_to_ret(:)
    type(variable_t), intent(inout) :: var_data(:)
    type(timer_t), optional,     intent(inout)   :: wait_timer

    integer :: n, p, k_max, i, j, k, n_vars
    integer :: halo_size, kms, kme, ims, jms, its, ite, jts, jte
    integer :: kms_var, kme_var

    if (this%n_3d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    halo_size = this%halo_size
    kms = this%kms; kme = this%kme
    ims = this%ims; jms = this%jms
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

#ifdef USE_NCCL
        ! NCCL shares OpenACC's stream — within-stream ordering ensures NCCL
        ! writes to receive buffers complete before the unpack kernels below
        ! read them. No explicit sync needed.
#else
    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Wait(this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Wait(this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Wait(this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Wait(this%northeast_3d_win)
#endif

    n = 1
    n_vars = size(var_data)
    if (present(wait_timer)) call wait_timer%start()

    ! Now iterate through the dictionary as long as there are more elements present
    do p = 1, size(vars_to_ret)
        if (vars_to_ret(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(vars_to_ret(p)%v)%id == vars_to_ret(p)%id) then
                associate(var => var_data(vars_to_ret(p)%v)%data_3d)
                kms_var = var_data(vars_to_ret(p)%v)%grid%kts
                kme_var = var_data(vars_to_ret(p)%v)%grid%kte
                if (.not.(this%north_boundary)) then
                    associate(buff => this%north_batch_in_3d)
                    !$acc parallel loop gang vector collapse(3) present(var, buff)
                    do j = 1,halo_size
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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
                    do k = kms_var,kme_var
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

    if (present(wait_timer)) call wait_timer%stop()

#ifndef USE_NCCL
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)
    if (.not.(this%northwest_boundary)) call MPI_Win_Post(this%northwest_neighbor_grp, 0, this%northwest_3d_win)
    if (.not.(this%southeast_boundary)) call MPI_Win_Post(this%southeast_neighbor_grp, 0, this%southeast_3d_win)
    if (.not.(this%southwest_boundary)) call MPI_Win_Post(this%southwest_neighbor_grp, 0, this%southwest_3d_win)
    if (.not.(this%northeast_boundary)) call MPI_Win_Post(this%northeast_neighbor_grp, 0, this%northeast_3d_win)
#endif

end subroutine halo_3d_retrieve_batch

module subroutine halo_2d_send_batch(this, vars_to_send, var_data)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: vars_to_send(:)
    type(variable_t), intent(inout) :: var_data(:)
    integer :: n, p, msg_size, i, j, k, n_vars
    integer :: halo_size, its, ite, jts, jte
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
    integer :: ierr, ns_msg_size, ew_msg_size
    if (this%n_2d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    msg_size = 1
    disp = 0

    halo_size = this%halo_size
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

#ifndef USE_NCCL
    if (.not.(this%north_boundary)) call MPI_Win_Start(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Start(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Start(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Start(this%west_neighbor_grp, 0, this%west_2d_win)
#endif

    n = 1
    n_vars = size(var_data)

    ! Now iterate through the exchange-only objects as long as there are more elements present
    do p = 1, size(vars_to_send)
        if (vars_to_send(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(vars_to_send(p)%v)%id == vars_to_send(p)%id) then
                if (var_data(vars_to_send(p)%v)%dtype == kINTEGER) then
                    associate(var => var_data(vars_to_send(p)%v)%data_2di)
                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            buff(n,i-its+1,j) = var(i,(jte-halo_size+j))
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            buff(n,i-its+1,j) = var(i,jts+j-1)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            buff(n,i,j-jts+1) = var((ite-halo_size+i),j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            buff(n,i,j-jts+1) = var(its+i-1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    end associate
                else
                    associate(var => var_data(vars_to_send(p)%v)%data_2d)
                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            buff(n,i-its+1,j) = var(i,(jte-halo_size+j))
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            buff(n,i-its+1,j) = var(i,jts+j-1)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            buff(n,i,j-jts+1) = var((ite-halo_size+i),j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            buff(n,i,j-jts+1) = var(its+i-1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    end associate
                endif
                n = n+1
            endif
        endif
    enddo

#ifdef USE_NCCL

        ns_msg_size = this%n_2d * this%grid%ns_halo_nx * this%halo_size
        ew_msg_size = this%n_2d * this%halo_size * this%grid%ew_halo_ny

        ! NCCL operations share OpenACC's stream, so within-stream ordering
        ! automatically serializes the preceding pack kernels before NCCL reads.
        call nccl_group_start()

        if (.not. this%north_boundary) then
            !$acc host_data use_device(this%north_buffer_2d, this%north_batch_in_2d)
            ierr = nccl_send_float(c_loc(this%north_buffer_2d), ns_msg_size, &
                this%north_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%north_batch_in_2d), ns_msg_size, &
                this%north_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%south_boundary) then
            !$acc host_data use_device(this%south_buffer_2d, this%south_batch_in_2d)
            ierr = nccl_send_float(c_loc(this%south_buffer_2d), ns_msg_size, &
                this%south_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%south_batch_in_2d), ns_msg_size, &
                this%south_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%east_boundary) then
            !$acc host_data use_device(this%east_buffer_2d, this%east_batch_in_2d)
            ierr = nccl_send_float(c_loc(this%east_buffer_2d), ew_msg_size, &
                this%east_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%east_batch_in_2d), ew_msg_size, &
                this%east_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        if (.not. this%west_boundary) then
            !$acc host_data use_device(this%west_buffer_2d, this%west_batch_in_2d)
            ierr = nccl_send_float(c_loc(this%west_buffer_2d), ew_msg_size, &
                this%west_neighbor, this%nccl_comm, this%nccl_stream)
            ierr = nccl_recv_float(c_loc(this%west_batch_in_2d), ew_msg_size, &
                this%west_neighbor, this%nccl_comm, this%nccl_stream)
            !$acc end host_data
        endif

        call nccl_group_end()
#else
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
    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_2d_win)
#endif

end subroutine halo_2d_send_batch

module subroutine halo_2d_retrieve_batch(this, vars_to_ret, var_data)
    implicit none
    class(halo_t), intent(inout) :: this
    type(index_type), intent(inout) :: vars_to_ret(:)
    type(variable_t), intent(inout) :: var_data(:)
    integer :: n, p, i, j, k, n_vars
    integer :: halo_size, ims, jms, its, ite, jts, jte

    if (this%n_2d <= 0 .or. (this%north_boundary.and.this%east_boundary.and.this%south_boundary.and.this%west_boundary)) return

    halo_size = this%halo_size
    ims = this%ims; jms = this%jms
    its = this%its; ite = this%ite
    jts = this%jts; jte = this%jte

#ifdef USE_NCCL
        ! NCCL shares OpenACC's stream — within-stream ordering ensures NCCL
        ! writes to receive buffers complete before the unpack kernels below
        ! read them. No explicit sync needed.
#else
    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_2d_win)
#endif

    n = 1    
    n_vars = size(var_data)

    ! Now iterate through the exchange-only objects as long as there are more elements present
    do p = 1, size(vars_to_ret)
        if (vars_to_ret(p)%v <= n_vars) then
            !check that the variable indexed matches the name
            if (var_data(vars_to_ret(p)%v)%id == vars_to_ret(p)%id) then
                if (var_data(vars_to_ret(p)%v)%dtype==kINTEGER) then
                    associate(var => var_data(vars_to_ret(p)%v)%data_2di)
                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            var(i,(jte+j)) = buff(n,i-its+1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            var(i,jms+j-1) = buff(n,i-its+1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            var((ite+i),j) = buff(n,i,j-jts+1)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            var(ims+i-1,j) = buff(n,i,j-jts+1)
                        enddo
                        enddo
                        end associate
                    endif
                    end associate
                else
                    associate(var => var_data(vars_to_ret(p)%v)%data_2d)
                    if (.not.(this%north_boundary)) then
                        associate(buff => this%north_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            var(i,(jte+j)) = buff(n,i-its+1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%south_boundary)) then
                        associate(buff => this%south_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = 1,halo_size
                        do i = its,ite
                            var(i,jms+j-1) = buff(n,i-its+1,j)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%east_boundary)) then
                        associate(buff => this%east_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            var((ite+i),j) = buff(n,i,j-jts+1)
                        enddo
                        enddo
                        end associate
                    endif
                    if (.not.(this%west_boundary)) then
                        associate(buff => this%west_batch_in_2d)
                        !$acc parallel loop gang vector collapse(2) present(var, buff)
                        do j = jts,jte
                        do i = 1,halo_size
                            var(ims+i-1,j) = buff(n,i,j-jts+1)
                        enddo
                        enddo
                        end associate
                    endif
                    end associate
                endif
                n = n+1
            endif
        endif
    enddo

#ifndef USE_NCCL
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_2d_win)
#endif

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
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs=var%ystag
    disp = 0
    msg_size = 1
    indx_start = var%grid%jte-offs-var%grid%halo_size+1
    associate(its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        north_in_buffer => this%north_in_buffer, north_in_buffer_2d => this%north_in_buffer_2d)
    !$acc data present(north_in_buffer, north_in_buffer_2d)
    if (var%two_d) then
#ifdef USE_NCCL
        ! Pack 2D data into k=1 slice of the 3D buffer so the flat size
        ! matches the 3D receive buffer (north_in_3d). The receiver's
        ! retrieve_*_halo only reads k=1 for 2D vars, so other slices are don't-care.
        if(var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = indx_start, jte
                do i = its, ite
                    north_in_buffer(i-its+1, 1, j-indx_start+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = indx_start, jte
                do i = its, ite
                    north_in_buffer(i-its+1, 1, j-indx_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
#else
        if(var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = indx_start, jte
                do i = its, ite
                    north_in_buffer_2d(i-its+1,j-indx_start+1) = data_2di(i,j)
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = indx_start, jte
                do i = its, ite
                    north_in_buffer_2d(i-its+1,j-indx_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc host_data use_device(north_in_buffer_2d)
        call MPI_Put(north_in_buffer_2d, msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
        !$acc end host_data
#endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = indx_start, jte
                do k = kts, kte
                    do i = its, ite
                        north_in_buffer(i-its+1,k-kts+1,j-indx_start+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = indx_start, jte
                do k = kts, kte
                    do i = its, ite
                        north_in_buffer(i-its+1,k-kts+1,j-indx_start+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(north_in_buffer)
        call MPI_Put(north_in_buffer, msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
        !$acc end host_data
#endif
    endif

#ifdef USE_NCCL
    ! Send the entire 3D buffer; receiver's recv buffer (north_in_3d on the
    ! peer) has identical shape, so element-for-element layout matches.
    nccl_count = int(size(this%north_in_buffer), c_int)
    !$acc host_data use_device(north_in_buffer)
    nccl_ierr = nccl_send_float(c_loc(north_in_buffer), nccl_count, &
        this%north_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif

    !$acc end data
    end associate
end subroutine



module subroutine put_south(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: start, nx, offs, msg_size, i, k, j, indx_end
    logical :: dqdt
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs=var%ystag
    disp = 0
    msg_size = 1
    indx_end = var%grid%jts+offs+var%grid%halo_size-1

    associate(its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        south_in_buffer => this%south_in_buffer, south_in_buffer_2d => this%south_in_buffer_2d)
    !$acc data present(south_in_buffer, south_in_buffer_2d)

    if (var%two_d) then
#ifdef USE_NCCL
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, indx_end
                do i = its, ite
                    south_in_buffer(i-its+1, 1, j-jts+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, indx_end
                do i = its, ite
                    south_in_buffer(i-its+1, 1, j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
#else
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, indx_end
                do i = its, ite
                    south_in_buffer_2d(i-its+1,j-jts+1) = data_2di(i,j)
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, indx_end
                do i = its, ite
                    south_in_buffer_2d(i-its+1,j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc host_data use_device(south_in_buffer_2d)
        call MPI_Put(south_in_buffer_2d, msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
        !$acc end host_data
#endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = jts, indx_end
                do k = kts, kte
                    do i = its, ite
                        south_in_buffer(i-its+1,k-kts+1,j-jts+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = jts, indx_end
                do k = kts, kte
                    do i = its, ite
                        south_in_buffer(i-its+1,k-kts+1,j-jts+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(south_in_buffer)
        call MPI_Put(south_in_buffer, msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
        !$acc end host_data
#endif
    endif

#ifdef USE_NCCL
    nccl_count = int(size(this%south_in_buffer), c_int)
    !$acc host_data use_device(south_in_buffer)
    nccl_ierr = nccl_send_float(c_loc(south_in_buffer), nccl_count, &
        this%south_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif

    !$acc end data
    end associate

end subroutine


module subroutine put_east(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: n, ny, msg_size, offs, i, k, j, indx_start
    logical :: dqdt
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs=var%xstag
    disp = 0
    msg_size = 1
    indx_start = var%grid%ite-offs-var%grid%halo_size+1

    associate(its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        east_in_buffer => this%east_in_buffer, east_in_buffer_2d => this%east_in_buffer_2d)
    !$acc data present(east_in_buffer, east_in_buffer_2d)

    if (var%two_d) then
#ifdef USE_NCCL
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, jte
                do i = indx_start, ite
                    east_in_buffer(i-indx_start+1, 1, j-jts+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, jte
                do i = indx_start, ite
                    east_in_buffer(i-indx_start+1, 1, j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
#else
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, jte
                do i = indx_start, ite
                    east_in_buffer_2d(i-indx_start+1,j-jts+1) = data_2di(i,j)
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, jte
                do i = indx_start, ite
                    east_in_buffer_2d(i-indx_start+1,j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc host_data use_device(east_in_buffer_2d)
        call MPI_Put(east_in_buffer_2d, msg_size, &
        var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
        !$acc end host_data
#endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = jts, jte
                do k = kts, kte
                    do i = indx_start, ite
                        east_in_buffer(i-indx_start+1,k-kts+1,j-jts+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = jts, jte
                do k = kts, kte
                    do i = indx_start, ite
                        east_in_buffer(i-indx_start+1,k-kts+1,j-jts+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(east_in_buffer)
        call MPI_Put(east_in_buffer, msg_size, &
        var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
        !$acc end host_data
#endif
    endif

#ifdef USE_NCCL
    nccl_count = int(size(this%east_in_buffer), c_int)
    !$acc host_data use_device(east_in_buffer)
    nccl_ierr = nccl_send_float(c_loc(east_in_buffer), nccl_count, &
        this%east_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif

    !$acc end data
    end associate

end subroutine




module subroutine put_west(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: start, ny, msg_size, offs, i, k, j, indx_end
    logical :: dqdt
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs=var%xstag
    disp = 0
    msg_size = 1
    indx_end = var%grid%its+offs+var%grid%halo_size-1

    associate(its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        west_in_buffer => this%west_in_buffer, west_in_buffer_2d => this%west_in_buffer_2d)
    !$acc data present(west_in_buffer, west_in_buffer_2d)

    if (var%two_d) then
#ifdef USE_NCCL
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, jte
                do i = its, indx_end
                    west_in_buffer(i-its+1, 1, j-jts+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, jte
                do i = its, indx_end
                    west_in_buffer(i-its+1, 1, j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
#else
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts, jte
                do i = its, indx_end
                    west_in_buffer_2d(i-its+1,j-jts+1) = data_2di(i,j)
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts, jte
                do i = its, indx_end
                    west_in_buffer_2d(i-its+1,j-jts+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc host_data use_device(west_in_buffer_2d)
        call MPI_Put(west_in_buffer_2d, msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
        !$acc end host_data
#endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = jts, jte
                do k = kts, kte
                    do i = its, indx_end
                        west_in_buffer(i-its+1,k-kts+1,j-jts+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = jts, jte
                do k = kts, kte
                    do i = its, indx_end
                        west_in_buffer(i-its+1,k-kts+1,j-jts+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(west_in_buffer)
        call MPI_Put(west_in_buffer, msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
        !$acc end host_data
#endif
    endif

#ifdef USE_NCCL
    nccl_count = int(size(this%west_in_buffer), c_int)
    !$acc host_data use_device(west_in_buffer)
    nccl_ierr = nccl_send_float(c_loc(west_in_buffer), nccl_count, &
        this%west_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif

    !$acc end data
    end associate

end subroutine

module subroutine retrieve_north_halo(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: n, nx, offs_x, offs_y, i, k, j
    integer :: halo_off
    logical :: dqdt

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag

    ! Under USE_NCCL the recv lands at the start of north_in_3d (no derived-type
    ! offset). MPI mode places the data at i-offset halo_size due to NS_win_halo.
#ifdef USE_NCCL
    halo_off = 0
#else
    halo_off = this%halo_size
#endif

    associate(halo_size => this%halo_size, its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        north_in_3d => this%north_in_3d)
    !$acc data present(north_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            n = ubound(var%data_2di,2)
            nx = size(var%data_2di,1)
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = n-halo_size+1-offs_y, n
                do i = its, ite
                    data_2di(i,j) = north_in_3d(i-its+1+halo_off,1,j-(n-halo_size+1-offs_y)+1)
                enddo
            enddo
        else
            n = ubound(var%data_2d,2)
            nx = size(var%data_2d,1)
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = n-halo_size+1-offs_y, n
                do i = its, ite
                    data_2d(i,j) = north_in_3d(i-its+1+halo_off,1,j-(n-halo_size+1-offs_y)+1)
                enddo
            enddo
        endif
    else
        n = ubound(var%data_3d,3)
        nx = size(var%data_3d,1)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = n-halo_size+1-offs_y, n
                do k = kts, kte
                    do i = its, ite
                        dqdt_3d(i,k,j) = north_in_3d(i-its+1+halo_off,k,j-(n-halo_size+1-offs_y)+1)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = n-halo_size+1-offs_y, n
                do k = kts, kte
                    do i = its, ite
                        data_3d(i,k,j) = north_in_3d(i-its+1+halo_off,k,j-(n-halo_size+1-offs_y)+1)
                    enddo
                enddo
            enddo
        endif
    endif

    !$acc end data
    end associate

end subroutine

module subroutine retrieve_south_halo(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: start, nx, offs_y, offs_x, i, j, k
    integer :: halo_off
    logical :: dqdt

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag

#ifdef USE_NCCL
    halo_off = 0
#else
    halo_off = this%halo_size
#endif

    associate(halo_size => this%halo_size, its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
         data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
         south_in_3d => this%south_in_3d)
    !$acc data present(south_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            start = lbound(var%data_2di,2)
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = start,start+halo_size-1
            do i = its,ite
                data_2di(i,j) = south_in_3d(i-its+1+halo_off,1,j-start+1)
            enddo
            enddo
        else
            start = lbound(var%data_2d,2)
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = start,start+halo_size-1
            do i = its,ite
                data_2d(i,j) = south_in_3d(i-its+1+halo_off,1,j-start+1)
            enddo
            enddo
        endif
    else
        start = lbound(var%data_3d,3)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = start,start+halo_size-1
            do k = kts,kte
            do i = its,ite
                dqdt_3d(i,k,j) = south_in_3d(i-its+1+halo_off,k,j-start+1)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = start,start+halo_size-1
            do k = kts,kte
            do i = its,ite
                data_3d(i,k,j) = south_in_3d(i-its+1+halo_off,k,j-start+1)
            enddo
            enddo
            enddo
        endif
    endif

    !$acc end data
    end associate

end subroutine

module subroutine retrieve_east_halo(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: n, ny, offs_x, offs_y, i, j, k
    integer :: halo_off
    logical :: dqdt

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag

#ifdef USE_NCCL
    halo_off = 0
#else
    halo_off = this%halo_size
#endif

    associate(halo_size => this%halo_size, its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        east_in_3d => this%east_in_3d)
    !$acc data present(this%east_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            n = ubound(var%data_2di,1)
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts,jte
            do i = n-halo_size+1-offs_x,n
                data_2di(i,j) = east_in_3d(i-(n-halo_size+1-offs_x)+1,1,j-jts+1+halo_off)
            enddo
            enddo
        else
            n = ubound(var%data_2d,1)
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts,jte
            do i = n-halo_size+1-offs_x,n
                data_2d(i,j) = east_in_3d(i-(n-halo_size+1-offs_x)+1,1,j-jts+1+halo_off)
            enddo
            enddo
        endif
    else
        n = ubound(var%data_3d,1)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = jts,jte
            do k = kts,kte
            do i = n-halo_size+1-offs_x,n
                dqdt_3d(i,k,j) = east_in_3d(i-(n-halo_size+1-offs_x)+1,k,j-jts+1+halo_off)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = jts,jte
            do k = kts,kte
            do i = n-halo_size+1-offs_x,n
                data_3d(i,k,j) = east_in_3d(i-(n-halo_size+1-offs_x)+1,k,j-jts+1+halo_off)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
end subroutine

module subroutine retrieve_west_halo(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: start, ny, offs_x, offs_y, i, j, k
    integer :: halo_off
    logical :: dqdt

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag

#ifdef USE_NCCL
    halo_off = 0
#else
    halo_off = this%halo_size
#endif

    associate(halo_size => this%halo_size, its => var%grid%its, jts => var%grid%jts, kts => var%grid%kts, ite => var%grid%ite, jte => var%grid%jte, kte => var%grid%kte, &
        data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d, &
        west_in_3d => this%west_in_3d)
    !$acc data present(this%west_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            start = lbound(var%data_2di,1)
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = jts,jte
            do i = start,start+halo_size-1
                data_2di(i,j) = west_in_3d(i-start+1,1,j-jts+1+halo_off)
            enddo
            enddo
        else
            start = lbound(var%data_2d,1)
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = jts,jte
            do i = start,start+halo_size-1
                data_2d(i,j) = west_in_3d(i-start+1,1,j-jts+1+halo_off)
            enddo
            enddo
        endif
    else
        start = lbound(var%data_3d,1)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = jts,jte
            do k = kts,kte
            do i = start,start+halo_size-1
                dqdt_3d(i,k,j) = west_in_3d(i-start+1,k,j-jts+1+halo_off)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = jts,jte
            do k = kts,kte
            do i = start,start+halo_size-1
                data_3d(i,k,j) = west_in_3d(i-start+1,k,j-jts+1+halo_off)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
end subroutine



module subroutine put_northeast(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    logical :: dqdt
    integer :: msg_size, offs_x, offs_y, i_start, j_start, i, j, k
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifndef USE_NCCL
    type(MPI_Win) :: dst_win
#endif
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag
    i_start = var%grid%ite - this%halo_size + 1 - offs_x
    j_start = var%grid%jte - this%halo_size + 1 - offs_y
    disp = 0
    msg_size = 1

#ifndef USE_NCCL
    ! select which destination window to use for MPI_Put. This handles
    ! the case for corner exchanges on domain boundaries
    if (this%north_boundary) then
        dst_win = this%northwest_in_win
        j_start = j_start + this%halo_size
    elseif (this%east_boundary) then
        dst_win = this%southeast_in_win
        i_start = i_start + this%halo_size
    else
        dst_win = this%southwest_in_win
    end if
#endif

    associate(data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    if (var%two_d) then
#ifdef USE_NCCL
        ! Pack 2D corner block into k=1 slice of the 3D corner send buffer.
        associate(buf => this%ne_corner_send)
        !$acc data present(buf)
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc end data
        end associate
#else
        if (var%dtype==kINTEGER) then
            !$acc host_data use_device(data_2di)
            call MPI_Put(data_2di(i_start,j_start), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        else
            !$acc host_data use_device(data_2d)
            call MPI_Put(data_2d(i_start,j_start), msg_size, &
                var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        endif
#endif
    else
        associate(kts => var%grid%kts, kte => var%grid%kte, buf => this%ne_corner_send)
        !$acc data present(buf)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(buf)
        call MPI_Put(buf, msg_size, &
            var%grid%corner_win_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
        !$acc end host_data
#endif
        !$acc end data
        end associate
    endif
    end associate

#ifdef USE_NCCL
    nccl_count = int(size(this%ne_corner_send), c_int)
    !$acc host_data use_device(this%ne_corner_send)
    nccl_ierr = nccl_send_float(c_loc(this%ne_corner_send), nccl_count, &
        this%northeast_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif
end subroutine

module subroutine put_northwest(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    logical :: dqdt
    integer :: msg_size, offs_x, offs_y, i_start, j_start, i, j, k
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifndef USE_NCCL
    type(MPI_Win) :: dst_win
#endif
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag
    i_start = var%grid%its
    j_start = var%grid%jte - this%halo_size + 1 - offs_y

    disp = 0
    msg_size = 1

#ifndef USE_NCCL
    ! select which destination window to use for MPI_Put. This handles
    ! the case for corner exchanges on domain boundaries
    if (this%north_boundary) then
        dst_win = this%northeast_in_win
        j_start = j_start + this%halo_size
    elseif (this%west_boundary) then
        dst_win = this%southwest_in_win
        i_start = i_start - this%halo_size
    else
        dst_win = this%southeast_in_win
    end if
#endif

    associate(data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    if (var%two_d) then
#ifdef USE_NCCL
        associate(buf => this%nw_corner_send)
        !$acc data present(buf)
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc end data
        end associate
#else
        if (var%dtype==kINTEGER) then
            !$acc host_data use_device(data_2di)
            call MPI_Put(data_2di(i_start,j_start), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        else
            !$acc host_data use_device(data_2d)
            call MPI_Put(data_2d(i_start,j_start), msg_size, &
                var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        endif
#endif
    else
        associate(kts => var%grid%kts, kte => var%grid%kte, buf => this%nw_corner_send)
        !$acc data present(buf)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(buf)
        call MPI_Put(buf, msg_size, &
            var%grid%corner_win_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
        !$acc end host_data
#endif
        !$acc end data
        end associate
    endif
    end associate

#ifdef USE_NCCL
    nccl_count = int(size(this%nw_corner_send), c_int)
    !$acc host_data use_device(this%nw_corner_send)
    nccl_ierr = nccl_send_float(c_loc(this%nw_corner_send), nccl_count, &
        this%northwest_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif
end subroutine


module subroutine put_southwest(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    logical :: dqdt
    integer :: msg_size, offs_x, offs_y, i_start, j_start, i, j, k
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifndef USE_NCCL
    type(MPI_Win) :: dst_win
#endif
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag
    disp = 0
    msg_size = 1
    i_start = var%grid%its
    j_start = var%grid%jts

#ifndef USE_NCCL
    ! select which destination window to use for MPI_Put. This handles
    ! the case for corner exchanges on domain boundaries
    if (this%south_boundary) then
        dst_win = this%southeast_in_win
        j_start = j_start - this%halo_size
    elseif (this%west_boundary) then
        dst_win = this%northwest_in_win
        i_start = i_start - this%halo_size
    else
        dst_win = this%northeast_in_win
    end if
#endif

    associate(data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    if (var%two_d) then
#ifdef USE_NCCL
        associate(buf => this%sw_corner_send)
        !$acc data present(buf)
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc end data
        end associate
#else
        if (var%dtype==kINTEGER) then
            !$acc host_data use_device(data_2di)
            call MPI_Put(data_2di(i_start,j_start), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        else
            !$acc host_data use_device(data_2d)
            call MPI_Put(data_2d(i_start,j_start), msg_size, &
                var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        endif
#endif
    else
        associate(kts => var%grid%kts, kte => var%grid%kte, buf => this%sw_corner_send)
        !$acc data present(buf)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(buf)
        call MPI_Put(buf, msg_size, &
            var%grid%corner_win_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
        !$acc end host_data
#endif
        !$acc end data
        end associate
    endif
    end associate

#ifdef USE_NCCL
    nccl_count = int(size(this%sw_corner_send), c_int)
    !$acc host_data use_device(this%sw_corner_send)
    nccl_ierr = nccl_send_float(c_loc(this%sw_corner_send), nccl_count, &
        this%southwest_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif
end subroutine

module subroutine put_southeast(this,var,do_dqdt)
    implicit none
    class(halo_t), intent(inout) :: this
    class(variable_t), intent(in) :: var
    logical, optional, intent(in) :: do_dqdt
    logical :: dqdt
    integer :: msg_size, offs_x, offs_y, i_start, j_start, i, j, k
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
#ifndef USE_NCCL
    type(MPI_Win) :: dst_win
#endif
#ifdef USE_NCCL
    integer(c_int) :: nccl_count, nccl_ierr
#endif

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    offs_x=var%xstag
    offs_y=var%ystag
    i_start = var%grid%ite - this%halo_size + 1 - offs_x
    j_start = var%grid%jts
    disp = 0
    msg_size = 1

#ifndef USE_NCCL
    ! select which destination window to use for MPI_Put. This handles
    ! the case for corner exchanges on domain boundaries
    if (this%south_boundary) then
        dst_win = this%southwest_in_win
        j_start = j_start - this%halo_size
    elseif (this%east_boundary) then
        dst_win = this%northeast_in_win
        i_start = i_start + this%halo_size
    else
        dst_win = this%northwest_in_win
    end if
#endif

    associate(data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    if (var%two_d) then
#ifdef USE_NCCL
        associate(buf => this%se_corner_send)
        !$acc data present(buf)
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = real(data_2di(i,j))
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do i = i_start, i_start + this%halo_size + offs_x - 1
                    buf(i-i_start+1, 1, j-j_start+1) = data_2d(i,j)
                enddo
            enddo
        endif
        !$acc end data
        end associate
#else
        if (var%dtype==kINTEGER) then
            !$acc host_data use_device(data_2di)
            call MPI_Put(data_2di(i_start,j_start), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        else
            !$acc host_data use_device(data_2d)
            call MPI_Put(data_2d(i_start,j_start), msg_size, &
                var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
            !$acc end host_data
        endif
#endif
    else
        associate(kts => var%grid%kts, kte => var%grid%kte, buf => this%se_corner_send)
        !$acc data present(buf)
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = dqdt_3d(i,k,j)
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_start + this%halo_size + offs_y - 1
                do k = kts, kte
                    do i = i_start, i_start + this%halo_size + offs_x - 1
                        buf(i-i_start+1, k-kts+1, j-j_start+1) = data_3d(i,k,j)
                    enddo
                enddo
            enddo
        endif
#ifndef USE_NCCL
        !$acc host_data use_device(buf)
        call MPI_Put(buf, msg_size, &
            var%grid%corner_win_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_win_halo, dst_win)
        !$acc end host_data
#endif
        !$acc end data
        end associate
    endif
    end associate

#ifdef USE_NCCL
    nccl_count = int(size(this%se_corner_send), c_int)
    !$acc host_data use_device(this%se_corner_send)
    nccl_ierr = nccl_send_float(c_loc(this%se_corner_send), nccl_count, &
        this%southeast_neighbor, this%nccl_comm, this%nccl_stream)
    !$acc end host_data
#endif
end subroutine


module subroutine retrieve_northeast_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: offs_x, offs_y, i, j, k
    logical :: dqdt
    integer :: i_start, i_end, j_start, j_end


    offs_x=var%xstag
    offs_y=var%ystag
    i_start=var%grid%ite+1-offs_x
    i_end=var%grid%ime
    j_start=var%grid%jte+1-offs_y
    j_end=var%grid%jme

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt
  
    associate(kts => var%grid%kts, kte => var%grid%kte, northeast_in_3d => this%northeast_in_3d, &
            data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    !$acc data present(northeast_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2di(i,j) = northeast_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2d(i,j) = northeast_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                dqdt_3d(i,k,j) = northeast_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                data_3d(i,k,j) = northeast_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
end subroutine

module subroutine retrieve_northwest_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: offs_x, offs_y, n, nx, i, j, k
    logical :: dqdt
    integer :: i_start, i_end, j_start, j_end


    offs_x=var%xstag
    offs_y=var%ystag
    i_start=var%grid%ims
    i_end=var%grid%its-1
    j_start=var%grid%jte+1-offs_y
    j_end=var%grid%jme

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt
  
    associate(kts => var%grid%kts, kte => var%grid%kte, northwest_in_3d => this%northwest_in_3d, &
            data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    !$acc data present(northwest_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2di(i,j) = northwest_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2d(i,j) = northwest_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                dqdt_3d(i,k,j) = northwest_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                data_3d(i,k,j) = northwest_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
end subroutine

module subroutine retrieve_southwest_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: i, j, k
    logical :: dqdt
    integer :: i_start, i_end, j_start, j_end

    i_start=var%grid%ims
    i_end=var%grid%its-1
    j_start=var%grid%jms
    j_end=var%grid%jts-1

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt
  
    associate(kts => var%grid%kts, kte => var%grid%kte, southwest_in_3d => this%southwest_in_3d, &
            data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    !$acc data present(southwest_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2di(i,j) = southwest_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2d(i,j) = southwest_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                dqdt_3d(i,k,j) = southwest_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                data_3d(i,k,j) = southwest_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
end subroutine

module subroutine retrieve_southeast_halo(this,var,do_dqdt)    
    implicit none
    class(halo_t), intent(in) :: this
    class(variable_t), intent(inout) :: var
    logical, optional, intent(in) :: do_dqdt
    integer :: offs_x, offs_y, i, j, k
    logical :: dqdt
    integer :: i_start, i_end, j_start, j_end

    offs_x=var%xstag
    offs_y=var%ystag
    i_start=var%grid%ite+1-offs_x
    i_end=var%grid%ime
    j_start=var%grid%jms
    j_end=var%grid%jts-1

    dqdt=.False.
    if (present(do_dqdt)) dqdt=do_dqdt

    associate(kts => var%grid%kts, kte => var%grid%kte, southeast_in_3d => this%southeast_in_3d, &
            data_3d => var%data_3d, data_2d => var%data_2d, data_2di => var%data_2di, dqdt_3d => var%dqdt_3d)
    !$acc data present(southeast_in_3d)
    if (var%two_d) then
        if (var%dtype==kINTEGER) then
            !$acc parallel loop gang vector collapse(2) present(data_2di)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2di(i,j) = southeast_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(2) present(data_2d)
            do j = j_start, j_end
            do i = i_start, i_end
                data_2d(i,j) = southeast_in_3d(i-i_start+1,1,j-j_start+1)
            enddo
            enddo
        endif
    else
        if (dqdt) then
            !$acc parallel loop gang vector collapse(3) present(dqdt_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                dqdt_3d(i,k,j) = southeast_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3) present(data_3d)
            do j = j_start, j_end
            do k = kts, kte
            do i = i_start, i_end
                data_3d(i,k,j) = southeast_in_3d(i-i_start+1,k,j-j_start+1)
            enddo
            enddo
            enddo
        endif
    endif
    !$acc end data
    end associate
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
    integer :: neighbor_shared_rank(1)
    
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
        call MPI_Group_translate_ranks(this%north_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%north_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%south_boundary) then
        call MPI_Group_translate_ranks(this%south_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%south_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%east_boundary) then
        call MPI_Group_translate_ranks(this%east_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%east_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    
    if (.not. this%west_boundary) then
        call MPI_Group_translate_ranks(this%west_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%west_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif

    if (.not. this%northwest_boundary) then
        call MPI_Group_translate_ranks(this%northwest_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%northwest_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%northeast_boundary) then
        call MPI_Group_translate_ranks(this%northeast_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%northeast_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%southwest_boundary) then
        call MPI_Group_translate_ranks(this%southwest_neighbor_grp, 1, [0], shared_comm_grp, neighbor_shared_rank, ierr)
        this%southwest_shared = (neighbor_shared_rank(1) /= MPI_UNDEFINED)
    endif
    if (.not. this%southeast_boundary) then
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
!#ifndef USE_NCCL
    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%north_neighbor/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC

    allocate(this%north_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%north_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%north_batch_in_3d)
    call MPI_WIN_CREATE(this%north_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%north_3d_win, ierr)
    !$acc end host_data
#endif

    if (this%n_2d > 0) then
        allocate(this%north_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%north_batch_in_2d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%north_batch_in_2d)
        call MPI_WIN_CREATE(this%north_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%north_2d_win, ierr)
        !$acc end host_data
#endif
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
!#ifndef USE_NCCL
    ! Create a group for the north-south exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%south_neighbor, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC

    allocate(this%south_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%south_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%south_batch_in_3d)
    call MPI_WIN_CREATE(this%south_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%south_3d_win, ierr)
    !$acc end host_data
#endif

    if (this%n_2d > 0) then
        allocate(this%south_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%south_batch_in_2d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%south_batch_in_2d)
        call MPI_WIN_CREATE(this%south_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%south_2d_win, ierr)
        !$acc end host_data
#endif
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
!#ifndef USE_NCCL
    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%halo_rank, this%east_neighbor/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%east_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%east_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%east_batch_in_3d)
    call MPI_WIN_CREATE(this%east_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%east_3d_win, ierr)
    !$acc end host_data
#endif

    if (this%n_2d > 0) then
        allocate(this%east_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%east_batch_in_2d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%east_batch_in_2d)
        call MPI_WIN_CREATE(this%east_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%east_2d_win, ierr)
        !$acc end host_data
#endif
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
!#ifndef USE_NCCL
    ! Create a group for the east-west exchange
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/this%west_neighbor, this%halo_rank/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%west_batch_in_3d(this%n_3d, nx, nz, ny))
    !$acc enter data copyin(this%west_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%west_batch_in_3d)
    call MPI_WIN_CREATE(this%west_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%west_3d_win, ierr)
    !$acc end host_data
#endif

    if (this%n_2d > 0) then
        allocate(this%west_batch_in_2d(this%n_2d, nx, ny))
        !$acc enter data copyin(this%west_batch_in_2d)
#ifndef USE_NCCL
        !$acc host_data use_device(this%west_batch_in_2d)
        call MPI_WIN_CREATE(this%west_batch_in_2d, win_size_2d*real_size, real_size, info_in, tmp_MPI_comm, this%west_2d_win, ierr)
        !$acc end host_data
#endif
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
    integer :: real_size, rank1, rank2
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

!#ifndef USE_NCCL
    ! Create a group for the northwest-southwest exchange
    ! corner neighbor rank could be the direct corner neighbor,
    ! or an adjacent process, so we need to find min/max ranks
    rank1 = min(this%halo_rank, this%northwest_neighbor)
    rank2 = max(this%halo_rank, this%northwest_neighbor)
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/rank1, rank2/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%northwest_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%northwest_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%northwest_batch_in_3d)
    call MPI_WIN_CREATE(this%northwest_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%northwest_3d_win, ierr)
    !$acc end host_data
#endif
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
    integer :: real_size, rank1, rank2
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

!#ifndef USE_NCCL
    ! Create a group for the northeast-southeast exchange
    ! corner neighbor rank could be the direct corner neighbor,
    ! or an adjacent process, so we need to find min/max ranks
    rank1 = min(this%halo_rank, this%northeast_neighbor)
    rank2 = max(this%halo_rank, this%northeast_neighbor)
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/rank1, rank2/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%northeast_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%northeast_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%northeast_batch_in_3d)
    call MPI_WIN_CREATE(this%northeast_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%northeast_3d_win, ierr)
    !$acc end host_data
#endif
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
    integer :: real_size, rank1, rank2
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

!#ifndef USE_NCCL
    ! Create a group for the southwest-southeast exchange
    ! corner neighbor rank could be the direct corner neighbor,
    ! or an adjacent process, so we need to find min/max ranks
    rank1 = min(this%halo_rank, this%southwest_neighbor)
    rank2 = max(this%halo_rank, this%southwest_neighbor)
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/rank1, rank2/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%southwest_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%southwest_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%southwest_batch_in_3d)
    call MPI_WIN_CREATE(this%southwest_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%southwest_3d_win, ierr)
    !$acc end host_data
#endif
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
    integer :: real_size, rank1, rank2
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

!#ifndef USE_NCCL
    ! Create a group for the southeast-southwest exchange
    rank1 = min(this%southeast_neighbor, this%halo_rank)
    rank2 = max(this%southeast_neighbor, this%halo_rank)
    ! corner neighbor rank could be the direct corner neighbor,
    ! or an adjacent process, so we need to find min/max ranks
    call MPI_Comm_group(comms, comp_proc)
    call MPI_Group_incl(comp_proc, 2, (/rank1, rank2/), tmp_MPI_grp, ierr)
    call MPI_Comm_create_group(comms, tmp_MPI_grp, 0, tmp_MPI_comm, ierr)
!#endif

#ifdef _OPENACC
    allocate(this%southeast_batch_in_3d(this%n_3d, this%halo_size, nz, this%halo_size))
    !$acc enter data copyin(this%southeast_batch_in_3d)
#ifndef USE_NCCL
    !$acc host_data use_device(this%southeast_batch_in_3d)
    call MPI_WIN_CREATE(this%southeast_batch_in_3d, win_size*real_size, real_size, info_in, tmp_MPI_comm, this%southeast_3d_win, ierr)
    !$acc end host_data
#endif
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
