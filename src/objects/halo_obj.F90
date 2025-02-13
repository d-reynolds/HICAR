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

    n_neighbors = merge(0,1,this%south_boundary)  &
        +merge(0,1,this%north_boundary)  &
        +merge(0,1,this%east_boundary)   &
        +merge(0,1,this%west_boundary)
    n_neighbors = max(1, n_neighbors)
    allocate(this%neighbors(n_neighbors))

    this%ims = this%grid%ims; this%its = this%grid%its; this%ids = this%grid%ids
    this%ime = this%grid%ime; this%ite = this%grid%ite; this%ide = this%grid%ide
    this%kms = this%grid%kms; this%kts = this%grid%kts; this%kds = this%grid%kds
    this%kme = this%grid%kme; this%kte = this%grid%kte; this%kde = this%grid%kde
    this%jms = this%grid%jms; this%jts = this%grid%jts; this%jds = this%grid%jds
    this%jme = this%grid%jme; this%jte = this%grid%jte; this%jde = this%grid%jde

    my_index = FINDLOC(DOM_IMG_INDX,PE_RANK_GLOBAL+1,dim=1)

    !Compute cardinal direction neighbors
#ifdef CRAY_PE        
    if (.not.(this%south_boundary)) this%south_neighbor = DOM_IMG_INDX(my_index - this%grid%ximages)
    if (.not.(this%north_boundary)) this%north_neighbor = DOM_IMG_INDX(my_index + this%grid%ximages)
    if (.not.(this%east_boundary)) this%east_neighbor  = DOM_IMG_INDX(my_index + 1)
    if (.not.(this%west_boundary)) this%west_neighbor  = DOM_IMG_INDX(my_index - 1)

    current = 1
    if (.not. this%south_boundary) then
        this%neighbors(current) = this%south_neighbor
        current = current+1
    endif
    if (.not. this%north_boundary) then
        this%neighbors(current) = this%north_neighbor
        current = current+1
    endif
    if (.not. this%east_boundary) then
        this%neighbors(current) = this%east_neighbor
        current = current+1
    endif
    if (.not. this%west_boundary) then
        this%neighbors(current) = this%west_neighbor
        current = current+1
    endif

    ! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
    if (current == 1) then
        this%neighbors(current) = PE_RANK_GLOBAL+1
    endif

#else
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
#endif


    !Compute diagonal direction neighbors
#ifdef CRAY_PE
    if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) this%southwest_neighbor = DOM_IMG_INDX(this%halo_rank - this%grid%ximages - 1)
    if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) this%northwest_neighbor = DOM_IMG_INDX(this%halo_rank + this%grid%ximages - 1)
    if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) this%southeast_neighbor  = DOM_IMG_INDX(this%halo_rank - this%grid%ximages + 1)
    if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) this%northeast_neighbor  = DOM_IMG_INDX(this%halo_rank + this%grid%ximages + 1)
#else
    if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) this%southwest_neighbor = this%halo_rank - this%grid%ximages - 1
    if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) this%northwest_neighbor = this%halo_rank + this%grid%ximages - 1
    if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) this%southeast_neighbor  = this%halo_rank - this%grid%ximages + 1
    if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) this%northeast_neighbor  = this%halo_rank + this%grid%ximages + 1
#endif
    n_neighbors = merge(1,0,(.not.(this%south_boundary) .and. .not.(this%west_boundary)))  &
                +merge(1,0,(.not.(this%north_boundary) .and. .not.(this%west_boundary)))  &
                +merge(1,0,(.not.(this%south_boundary) .and. .not.(this%east_boundary)))   &
                +merge(1,0,(.not.(this%north_boundary) .and. .not.(this%east_boundary)))
    n_neighbors = max(1, n_neighbors)

    allocate(this%corner_neighbors(n_neighbors))

    current = 1
    if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) then
        this%corner_neighbors(current) = this%southwest_neighbor
        current = current+1
    endif
    if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) then
        this%corner_neighbors(current) = this%northwest_neighbor
        current = current+1
    endif
    if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) then
        this%corner_neighbors(current) = this%southeast_neighbor
        current = current+1
    endif
    if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) then
        this%corner_neighbors(current) = this%northeast_neighbor
        current = current+1
    endif

! if current = 1 then all of the boundaries were set, just store ourself as our "neighbor"
    if (current == 1) then
        this%corner_neighbors(current) = PE_RANK_GLOBAL+1
    endif
    
    ! Detect if neighbors are on shared memory hardware
    call detect_shared_memory(this, comms)

    !Now allocate the actual 3D halo
#ifdef CRAY_PE
allocate( this%south_in_3d( this%grid%ns_halo_nx, this%grid%halo_nz, this%halo_size+1   )[*])
allocate( this%north_in_3d( this%grid%ns_halo_nx, this%grid%halo_nz, this%halo_size        )[*])
allocate( this%east_in_3d(  this%halo_size       ,  this%grid%halo_nz, this%grid%ew_halo_ny  )[*])
allocate( this%west_in_3d(  this%halo_size+1,  this%grid%halo_nz, this%grid%ew_halo_ny)[*])

#else
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
endif
#endif


    !...and the larger 3D halo for batch exchanges
    call setup_batch_exch(this, exch_vars, adv_vars, comms)


end subroutine init

module subroutine finalize(this)
    class(halo_t), intent(inout) :: this

    integer :: ierr
    
    if (this%n_2d > 0) then
#ifdef CRAY_PE
        deallocate(this%north_batch_in_2d)
        deallocate(this%south_batch_in_2d)
        deallocate(this%east_batch_in_2d)
        deallocate(this%west_batch_in_2d)
#else

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
#endif
    endif

#ifdef CRAY_PE
    deallocate(this%north_batch_in_3d)
    deallocate(this%south_batch_in_3d)
    deallocate(this%east_batch_in_3d)
    deallocate(this%west_batch_in_3d)
#else

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
#endif

    deallocate(this%neighbors)
    deallocate(this%corner_neighbors)

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
#ifndef CRAY_PE
            call MPI_Win_fence(0,this%south_in_win)
            call MPI_Win_fence(0,this%north_in_win)
            call MPI_Win_fence(0,this%east_in_win)
            call MPI_Win_fence(0,this%west_in_win)
#endif
            if (.not. this%north_boundary .and. .not.this%west_boundary) call this%put_northeast(var, dqdt)
            if (.not. this%north_boundary .and. .not.this%east_boundary) call this%put_northwest(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%east_boundary)  call this%put_southeast(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%west_boundary)  call this%put_southwest(var, dqdt)

#ifdef CRAY_PE
            sync images( this%corner_neighbors )
#else
            call MPI_Win_fence(0,this%south_in_win)
            call MPI_Win_fence(0,this%north_in_win)
            call MPI_Win_fence(0,this%east_in_win)
            call MPI_Win_fence(0,this%west_in_win)
#endif

            if (.not. this%north_boundary .and. .not.this%west_boundary) call this%retrieve_northwest_halo(var, dqdt)
            if (.not. this%north_boundary .and. .not.this%east_boundary) call this%retrieve_northeast_halo(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%east_boundary)  call this%retrieve_southeast_halo(var, dqdt)
            if (.not. this%south_boundary .and. .not.this%west_boundary)  call this%retrieve_southwest_halo(var, dqdt)

            return
    endif

    ! if staggered in x direction, we need to carefully call the put and get commands
    if(var%xstag>0) then
#ifndef CRAY_PE                
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
#endif
        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

#ifdef CRAY_PE
        sync images( this%neighbors )
#else
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win) 
#endif

        if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
        if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)

#ifndef CRAY_PE                
        call MPI_Win_fence(0,this%south_in_win) 
        call MPI_Win_fence(0,this%north_in_win)
#endif
        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)

#ifdef CRAY_PE
        sync images( this%neighbors )
#else
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win) 
#endif

        if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
        if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)
    ! if staggered in y direction, we need to carefully call the put and get commands
    elseif(var%ystag>0) then
#ifndef CRAY_PE
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
#endif
        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)

#ifdef CRAY_PE
        sync images( this%neighbors )
#else
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
#endif

        if (.not. this%north_boundary) call this%retrieve_north_halo(var, dqdt)
        if (.not. this%south_boundary) call this%retrieve_south_halo(var, dqdt)

#ifndef CRAY_PE   
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
#endif
        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

#ifdef CRAY_PE
        sync images( this%neighbors )
#else
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)                
#endif
        if (.not. this%east_boundary)  call this%retrieve_east_halo(var, dqdt)
        if (.not. this%west_boundary)  call this%retrieve_west_halo(var, dqdt)

    else
#ifndef CRAY_PE
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
#endif
        if (.not. this%north_boundary) call this%put_north(var, dqdt)
        if (.not. this%south_boundary) call this%put_south(var, dqdt)
        if (.not. this%east_boundary)  call this%put_east(var, dqdt)
        if (.not. this%west_boundary)  call this%put_west(var, dqdt)

#ifdef CRAY_PE
        sync images( this%neighbors )
#else
        call MPI_Win_fence(0,this%south_in_win)
        call MPI_Win_fence(0,this%north_in_win)
        call MPI_Win_fence(0,this%east_in_win)
        call MPI_Win_fence(0,this%west_in_win)
#endif
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
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, win_size_2d
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

#ifdef CRAY_PE
    allocate(this%north_batch_in_3d(this%n_3d,1:this%grid%ns_halo_nx,&
                    this%kms:this%kme,1:this%halo_size)[*])
    allocate(this%south_batch_in_3d(this%n_3d,1:this%grid%ns_halo_nx,&
                    this%kms:this%kme,1:this%halo_size)[*])
    allocate(this%east_batch_in_3d(this%n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:this%grid%ew_halo_ny)[*])
    allocate(this%west_batch_in_3d(this%n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:this%grid%ew_halo_ny)[*])
    if (this%n_2d > 0) then
        allocate(this%north_batch_in_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size)[*])
        allocate(this%south_batch_in_2d(this%n_2d,1:this%grid%ns_halo_nx,1:this%halo_size)[*])
        allocate(this%east_batch_in_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny)[*])
        allocate(this%west_batch_in_2d(this%n_2d,1:this%halo_size,1:this%grid%ew_halo_ny)[*])
    endif
#else
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
        if ((mod(this%grid%yimg,2) == 1) .and. .not.(this%north_boundary)) then
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

        else if((mod(this%grid%yimg,2) == 0) .and. .not.(this%south_boundary)) then
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

        endif

        if ((mod(this%grid%yimg,2) == 0) .and. .not.(this%north_boundary)) then
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

        else if((mod(this%grid%yimg,2) == 1) .and. .not.(this%south_boundary)) then
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
        if ((mod(this%grid%ximg,2) == 1) .and. .not.(this%east_boundary)) then
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

        else if((mod(this%grid%ximg,2) == 0) .and. .not.(this%west_boundary)) then
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
        endif

        if ((mod(this%grid%ximg,2) == 0) .and. .not.(this%east_boundary)) then
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

        else if((mod(this%grid%ximg,2) == 1) .and. .not.(this%west_boundary)) then
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
        ! if (.not.(this%north_boundary)) call MPI_Win_fence(0, this%north_3d_win)
        ! if (.not.(this%south_boundary)) call MPI_Win_fence(0, this%south_3d_win)
        ! if (.not.(this%east_boundary)) call MPI_Win_fence(0, this%east_3d_win)
        ! if (.not.(this%west_boundary)) call MPI_Win_fence(0, this%west_3d_win)

        if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
        if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
        if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
        if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)
    
        if (this%n_2d > 0) then
            ! if (.not.(this%north_boundary)) call MPI_Win_fence(0, this%north_2d_win)
            ! if (.not.(this%south_boundary)) call MPI_Win_fence(0, this%south_2d_win)
            ! if (.not.(this%east_boundary)) call MPI_Win_fence(0, this%east_2d_win)
            ! if (.not.(this%west_boundary)) call MPI_Win_fence(0, this%west_2d_win)

            if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_2d_win)
            if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_2d_win)
            if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_2d_win)
            if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_2d_win)
        
        endif
    endif
#endif

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

    ! Now iterate through the dictionary as long as there are more elements present. If two processors are on shared memory
    ! this step will directly copy the data to the other PE
    if (.not.(exch_v_only)) then
        do while (adv_vars%has_more_elements())
            ! get the next variable
            var = adv_vars%next()
            if (var%three_d) then
                if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+3),:,:) = &
                        var%data_3d(this%its-1:this%ite+1,:,(this%jte-this%halo_size+1):this%jte)
                if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+3),:,:) = &
                        var%data_3d(this%its-1:this%ite+1,:,this%jts:(this%jts+this%halo_size-1))
                if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,:,1:(this%jte-this%jts+3)) = &
                        var%data_3d((this%ite-this%halo_size+1):this%ite,:,this%jts-1:this%jte+1)
                if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,:,1:(this%jte-this%jts+3)) = &
                        var%data_3d(this%its:(this%its+this%halo_size-1),:,this%jts-1:this%jte+1)

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
            if (.not.(this%north_boundary)) this%north_buffer_3d(n,1:(this%ite-this%its+3),1:k_max,:) = &
                    var%data_3d(this%its-1:this%ite+1,1:k_max,(this%jte-this%halo_size+1):this%jte)
            if (.not.(this%south_boundary)) this%south_buffer_3d(n,1:(this%ite-this%its+3),1:k_max,:) = &
                    var%data_3d(this%its-1:this%ite+1,1:k_max,this%jts:(this%jts+this%halo_size-1))
            if (.not.(this%east_boundary)) this%east_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+3)) = &
                    var%data_3d((this%ite-this%halo_size+1):this%ite,1:k_max,this%jts-1:this%jte+1)
            if (.not.(this%west_boundary)) this%west_buffer_3d(n,:,1:k_max,1:(this%jte-this%jts+3)) = &
                    var%data_3d(this%its:(this%its+this%halo_size)-1,1:k_max,this%jts-1:this%jte+1)

            n = n+1
        endif
    enddo

    if (.not.(this%south_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%north_batch_in_3d(:,:,:,:)[this%south_neighbor] = this%south_buffer_3d(:,:,:,:)
#else
        if (.not.(this%south_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%south_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 0, disp, msg_size, &
                this%NS_3d_win_halo_type, this%south_3d_win)
        endif
#endif        
    endif

    if (.not.(this%north_boundary)) then
#ifdef CRAY_PE
        !DIR$ PGAS DEFER_SYNC
        this%south_batch_in_3d(:,:,:,:)[this%north_neighbor] = this%north_buffer_3d(:,:,:,:)
#else
        if (.not.(this%north_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%north_buffer_3d, msg_size, &
                this%NS_3d_win_halo_type, 1, disp, msg_size, &
                this%NS_3d_win_halo_type, this%north_3d_win)
        endif
#endif
    endif

    if (.not.(this%east_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%west_batch_in_3d(:,:,:,:)[this%east_neighbor] = this%east_buffer_3d(:,:,:,:)
#else
        if (.not.(this%east_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%east_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 1, disp, msg_size, &
                this%EW_3d_win_halo_type, this%east_3d_win)
        endif
#endif  
    endif

    if (.not.(this%west_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%east_batch_in_3d(:,:,:,:)[this%west_neighbor] = this%west_buffer_3d(:,:,:,:)
#else
        if (.not.(this%west_shared)) then
            ! Use post-start-complete-wait for distributed memory
            call MPI_Put(this%west_buffer_3d, msg_size, &
                this%EW_3d_win_halo_type, 0, disp, msg_size, &
                this%EW_3d_win_halo_type, this%west_3d_win)
        endif
#endif        
    endif

    if (.not.(this%north_boundary)) call MPI_Win_Complete(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Complete(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Complete(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Complete(this%west_3d_win)

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

#ifdef CRAY_PE        
    sync images( this%neighbors )
#else
    ! if (.not.(this%north_boundary)) call MPI_Win_fence(0, this%north_3d_win)
    ! if (.not.(this%south_boundary)) call MPI_Win_fence(0, this%south_3d_win)
    ! if (.not.(this%east_boundary)) call MPI_Win_fence(0, this%east_3d_win)
    ! if (.not.(this%west_boundary)) call MPI_Win_fence(0, this%west_3d_win)

    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_3d_win)

#endif


    call adv_vars%reset_iterator()
    call exch_vars%reset_iterator()
    n = 1
    ! Now iterate through the dictionary as long as there are more elements present
    if (.not.(exch_v_only)) then
        do while (adv_vars%has_more_elements())
            ! get the next variable
            var = adv_vars%next()
            if (var%three_d) then
                if (.not.(this%north_boundary)) var%data_3d(this%its-1:this%ite+1,:,(this%jte+1):this%jme) = &
                        this%north_batch_in_3d(n,1:(this%ite-this%its+3),:,1:this%halo_size)
                if (.not.(this%south_boundary)) var%data_3d(this%its-1:this%ite+1,:,this%jms:(this%jts-1)) = &
                        this%south_batch_in_3d(n,1:(this%ite-this%its+3),:,1:this%halo_size)
                if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,:,this%jts-1:this%jte+1) = &
                        this%east_batch_in_3d(n,:,:,1:(this%jte-this%jts+3))
                if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),:,this%jts-1:this%jte+1) = &
                        this%west_batch_in_3d(n,:,:,1:(this%jte-this%jts+3))
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
            if (.not.(this%north_boundary)) var%data_3d(this%its-1:this%ite+1,1:k_max,(this%jte+1):this%jme) = &
                    this%north_batch_in_3d(n,1:(this%ite-this%its+3),1:k_max,:)
            if (.not.(this%south_boundary)) var%data_3d(this%its-1:this%ite+1,1:k_max,this%jms:(this%jts-1)) = &
                    this%south_batch_in_3d(n,1:(this%ite-this%its+3),1:k_max,:)
            if (.not.(this%east_boundary)) var%data_3d((this%ite+1):this%ime,1:k_max,this%jts-1:this%jte+1) = &
                    this%east_batch_in_3d(n,:,1:k_max,1:(this%jte-this%jts+3))
            if (.not.(this%west_boundary)) var%data_3d(this%ims:(this%its-1),1:k_max,this%jts-1:this%jte+1) = &
                    this%west_batch_in_3d(n,:,1:k_max,1:(this%jte-this%jts+3))
            n = n+1
        endif
    enddo
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp, 0, this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp, 0, this%south_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp, 0, this%east_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp, 0, this%west_3d_win)

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
            if (.not.(this%north_boundary)) this%north_buffer_2d(n,1:(this%ite-this%its+3),:) = &
                    var%data_2d(this%its-1:this%ite+1,(this%jte-this%halo_size+1):this%jte)
            if (.not.(this%south_boundary)) this%south_buffer_2d(n,1:(this%ite-this%its+3),:) = &
                    var%data_2d(this%its-1:this%ite+1,this%jts:(this%jts+this%halo_size-1))
            if (.not.(this%east_boundary)) this%east_buffer_2d(n,:,1:(this%jte-this%jts+3)) = &
                    var%data_2d((this%ite-this%halo_size+1):this%ite,this%jts-1:this%jte+1)
            if (.not.(this%west_boundary)) this%west_buffer_2d(n,:,1:(this%jte-this%jts+3)) = &
                    var%data_2d(this%its:(this%its+this%halo_size)-1,this%jts-1:this%jte+1)

            n = n+1
        endif
    enddo

    if (.not.(this%north_boundary)) then
#ifdef CRAY_PE
        !DIR$ PGAS DEFER_SYNC
        this%south_batch_in_2d(:,:,:)[this%north_neighbor] = this%north_buffer_2d(:,:,:)
#else
        if (.not.(this%north_shared)) then
            call MPI_Put(this%north_buffer_2d, size(this%north_buffer_2d), &
                MPI_REAL, 1, disp, size(this%north_buffer_2d), MPI_REAL, this%north_2d_win)
        endif
#endif
    endif

    if (.not.(this%south_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%north_batch_in_2d(:,:,:)[this%south_neighbor] = this%south_buffer_2d(:,:,:)
#else
        if (.not.(this%south_shared)) then
            call MPI_Put(this%south_buffer_2d, size(this%south_buffer_2d), &
                MPI_REAL, 0, disp, size(this%south_buffer_2d), MPI_REAL, this%south_2d_win)
        endif
#endif        
    endif

    if (.not.(this%east_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%west_batch_in_2d(:,:,:)[this%east_neighbor] = this%east_buffer_2d(:,:,:)
#else
        if (.not.(this%east_shared)) then
            call MPI_Put(this%east_buffer_2d, size(this%east_buffer_2d), &
                MPI_REAL, 1, disp, size(this%east_buffer_2d), MPI_REAL, this%east_2d_win)
        endif
#endif  
    endif

    if (.not.(this%west_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%east_batch_in_2d(:,:,:)[this%west_neighbor] = this%west_buffer_2d(:,:,:)
#else
        if (.not.(this%west_shared)) then
            call MPI_Put(this%west_buffer_2d, size(this%west_buffer_2d), &
                MPI_REAL, 0, disp, size(this%west_buffer_2d), MPI_REAL, this%west_2d_win)
        endif
#endif        
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
#ifdef CRAY_PE        
    sync images( this%neighbors )
#else
    ! if (.not.(this%north_boundary)) call MPI_Win_fence(0, this%north_2d_win)
    ! if (.not.(this%south_boundary)) call MPI_Win_fence(0, this%south_2d_win)
    ! if (.not.(this%east_boundary)) call MPI_Win_fence(0, this%east_2d_win)
    ! if (.not.(this%west_boundary)) call MPI_Win_fence(0, this%west_2d_win)

    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_2d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_2d_win)
#endif

    call exch_vars%reset_iterator()
    n = 1    
    ! Now iterate through the exchange-only objects as long as there are more elements present
    do while (exch_vars%has_more_elements())
        ! get the next variable
        var = exch_vars%next()
        if (var%two_d) then
            if (.not.(this%north_boundary)) var%data_2d(this%its-1:this%ite+1,(this%jte+1):this%jme) = this%north_batch_in_2d(n,1:(this%ite-this%its+3),:)
            if (.not.(this%south_boundary)) var%data_2d(this%its-1:this%ite+1,this%jms:(this%jts-1)) = this%south_batch_in_2d(n,1:(this%ite-this%its+3),:)
            if (.not.(this%east_boundary)) var%data_2d((this%ite+1):this%ime,this%jts-1:this%jte+1) = this%east_batch_in_2d(n,:,1:(this%jte-this%jts+3))
            if (.not.(this%west_boundary)) var%data_2d(this%ims:(this%its-1),this%jts-1:this%jte+1) = this%west_batch_in_2d(n,:,1:(this%jte-this%jts+3))
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
  
  offs=0!var%ystag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      n = ubound(var%data_2d,2)
      nx = size(var%data_2d,1)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1+this%halo_size:nx-this%halo_size,1,1:(this%halo_size+offs))[this%north_neighbor] = var%dqdt_2d(var%grid%its:var%grid%ite,(n-this%halo_size*2+1-offs):(n-this%halo_size))
#else
          call MPI_Put(var%dqdt_2d(var%grid%ims,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1+this%halo_size:nx-this%halo_size,1,1:(this%halo_size+offs))[this%north_neighbor] = var%data_2d(var%grid%its:var%grid%ite,(n-this%halo_size*2+1-offs):(n-this%halo_size))
#else
          call MPI_Put(var%data_2d(var%grid%ims,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
#endif
      endif
  else
      n = ubound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:(this%halo_size+offs))[this%north_neighbor] = &
             var%dqdt_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,(n-this%halo_size*2+1-offs):(n-this%halo_size))
#else
          call MPI_Put(var%dqdt_3d(var%grid%ims,var%grid%kts,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:(this%halo_size+offs))[this%north_neighbor] = &
             var%data_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,(n-this%halo_size*2+1-offs):(n-this%halo_size))
#else
          call MPI_Put(var%data_3d(var%grid%ims,var%grid%kts,var%grid%jte-var%grid%halo_size+1-offs), msg_size, &
            var%grid%NS_halo, this%north_neighbor, disp, msg_size, var%grid%NS_win_halo, this%south_in_win)
#endif
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
  
  offs=0!var%ystag
  
  disp = 0
  msg_size = 1
  if (var%two_d) then
      start = lbound(var%data_2d,2)
      nx = size(var%data_2d,1)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1+this%halo_size:nx-this%halo_size,1,1:this%halo_size)[this%south_neighbor] = var%dqdt_2d(var%grid%its:var%grid%ite,(start+this%halo_size+offs):(start+this%halo_size*2-1+offs))
#else
          call MPI_Put(var%dqdt_2d(var%grid%ims,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1+this%halo_size:nx-this%halo_size,1,1:this%halo_size)[this%south_neighbor] = var%data_2d(var%grid%its:var%grid%ite,(start+this%halo_size+offs):(start+this%halo_size*2-1+offs))
#else
          call MPI_Put(var%data_2d(var%grid%ims,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
#endif
      endif
  else
      start = lbound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:this%halo_size)[this%south_neighbor] = &
            var%dqdt_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,(start+this%halo_size+offs):(start+this%halo_size*2-1+offs))
#else
          call MPI_Put(var%dqdt_3d(var%grid%ims,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:this%halo_size)[this%south_neighbor] = &
            var%data_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,(start+this%halo_size+offs):(start+this%halo_size*2-1+offs))
#else
          call MPI_Put(var%data_3d(var%grid%ims,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%NS_halo, this%south_neighbor, disp, msg_size, var%grid%NS_win_halo, this%north_in_win)
#endif
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
  
  offs=0!var%xstag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      n = ubound(var%data_2d,1)
      ny = size(var%data_2d,2)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:(this%halo_size+offs),1,1+this%halo_size:ny-this%halo_size)[this%east_neighbor] = var%dqdt_2d((n-this%halo_size*2+1-offs):(n-this%halo_size),var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%dqdt_2d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:(this%halo_size+offs),1,1+this%halo_size:ny-this%halo_size)[this%east_neighbor] = var%data_2d((n-this%halo_size*2+1-offs):(n-this%halo_size),var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%data_2d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
#endif
      endif
  else
      n = ubound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:(this%halo_size+offs),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)[this%east_neighbor] = &
            var%dqdt_3d((n-this%halo_size*2+1-offs):(n-this%halo_size),var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%dqdt_3d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:(this%halo_size+offs),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)[this%east_neighbor] =&
            var%data_3d((n-this%halo_size*2+1-offs):(n-this%halo_size),var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%data_3d(var%grid%ite-offs-var%grid%halo_size+1,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%east_neighbor, disp, msg_size, var%grid%EW_win_halo, this%west_in_win)
#endif
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
  
  offs=0!var%xstag
  disp = 0
  msg_size = 1

  if (var%two_d) then
      start = lbound(var%data_2d,1)
      ny = size(var%data_2d,2)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,1,1+this%halo_size:ny-this%halo_size)[this%west_neighbor] = var%dqdt_2d((start+this%halo_size+offs):(start+this%halo_size*2-1+offs),var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%dqdt_2d(var%grid%its,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,1,1+this%halo_size:ny-this%halo_size)[this%west_neighbor] = var%data_2d((start+this%halo_size+offs):(start+this%halo_size*2-1+offs),var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%data_2d(var%grid%its,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
#endif
      endif
  else
      start = lbound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)[this%west_neighbor] = &
            var%dqdt_3d((start+this%halo_size+offs):(start+this%halo_size*2-1+offs),var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)[this%west_neighbor] = &
            var%data_3d((start+this%halo_size+offs):(start+this%halo_size*2-1+offs),var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte)
#else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jms), msg_size, &
            var%grid%EW_halo, this%west_neighbor, disp, msg_size, var%grid%EW_win_halo, this%east_in_win)
#endif
      endif
  endif
end subroutine

module subroutine retrieve_north_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, nx, offs
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=0!var%ystag
  
  if (var%two_d) then
      n = ubound(var%data_2d,2)
      nx = size(var%data_2d,1)
      if (dqdt) then
          var%dqdt_2d(var%grid%its:var%grid%ite,n-this%halo_size+1:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size,1,1:this%halo_size)
      else
          var%data_2d(var%grid%its:var%grid%ite,n-this%halo_size+1:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size,1,1:this%halo_size)
      endif
  else
      n = ubound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          var%dqdt_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,n-this%halo_size+1:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:this%halo_size)
      else
          var%data_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,n-this%halo_size+1:n) = this%north_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_south_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, nx, offs
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=0!var%ystag
 
  
  if (var%two_d) then
      start = lbound(var%data_2d,2)
      nx = size(var%data_2d,1)
      if (dqdt) then
          var%dqdt_2d(var%grid%its:var%grid%ite,start:start+this%halo_size-1+offs) = this%south_in_3d(1+this%halo_size:nx-this%halo_size,1,1:(this%halo_size+offs))
      else
          var%data_2d(var%grid%its:var%grid%ite,start:start+this%halo_size-1+offs) = this%south_in_3d(1+this%halo_size:nx-this%halo_size,1,1:(this%halo_size+offs))
      endif
  else
      start = lbound(var%data_3d,3)
      nx = size(var%data_3d,1)
      if (dqdt) then
          var%dqdt_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,start:start+this%halo_size-1+offs) = this%south_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:(this%halo_size+offs))
      else
          var%data_3d(var%grid%its:var%grid%ite,var%grid%kts:var%grid%kte,start:start+this%halo_size-1+offs) = this%south_in_3d(1+this%halo_size:nx-this%halo_size,var%grid%kts:var%grid%kte,1:(this%halo_size+offs))
      endif
  endif
end subroutine

module subroutine retrieve_east_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: n, ny, offs
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=0!var%xstag

  if (var%two_d) then
      n = ubound(var%data_2d,1)
      ny = size(var%data_2d,2)
      if (dqdt) then
          var%dqdt_2d(n-this%halo_size+1:n,var%grid%jts:var%grid%jte) = this%east_in_3d(1:this%halo_size,1,1+this%halo_size:ny-this%halo_size)
      else
          var%data_2d(n-this%halo_size+1:n,var%grid%jts:var%grid%jte) = this%east_in_3d(1:this%halo_size,1,1+this%halo_size:ny-this%halo_size)
      endif
  else
      n = ubound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          var%dqdt_3d(n-this%halo_size+1:n,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte) = this%east_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)
      else
          var%data_3d(n-this%halo_size+1:n,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte) = this%east_in_3d(1:this%halo_size,var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_west_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt
  integer :: start, ny, offs
  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  offs=0!var%xstag
  
  if (var%two_d) then
      start = lbound(var%data_2d,1)
      ny = size(var%data_2d,2)
      if (dqdt) then
          var%dqdt_2d(start:start+this%halo_size-1+offs,var%grid%jts:var%grid%jte) = this%west_in_3d(1:(this%halo_size+offs),1,1+this%halo_size:ny-this%halo_size)
      else
          var%data_2d(start:start+this%halo_size-1+offs,var%grid%jts:var%grid%jte) = this%west_in_3d(1:(this%halo_size+offs),1,1+this%halo_size:ny-this%halo_size)
      endif
  else
      start = lbound(var%data_3d,1)
      ny = size(var%data_3d,3)
      if (dqdt) then
          var%dqdt_3d(start:start+this%halo_size-1+offs,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte) = this%west_in_3d(1:(this%halo_size+offs),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)
      else
          var%data_3d(start:start+this%halo_size-1+offs,var%grid%kts:var%grid%kte,var%grid%jts:var%grid%jte) = this%west_in_3d(1:(this%halo_size+offs),var%grid%kts:var%grid%kte,1+this%halo_size:ny-this%halo_size)
      endif
  endif
end subroutine



module subroutine put_northeast(this,var,do_dqdt)
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
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:this%halo_size,1,1:this%halo_size)[this%northeast_neighbor] = var%dqdt_2d(this%ite-this%halo_size:this%ite,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%dqdt_2d(var%grid%ite-this%halo_size,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:this%halo_size,1,1:this%halo_size)[this%northeast_neighbor] = var%data_2d(this%ite-this%halo_size:this%ite,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%data_2d(var%grid%ite-this%halo_size,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
#endif
      endif
  else
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%northeast_neighbor] = var%dqdt_3d(this%ite-this%halo_size:this%ite,this%kts:this%kte,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%dqdt_3d(var%grid%ite-this%halo_size,var%grid%kts,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%northeast_neighbor] = var%data_3d(this%ite-this%halo_size:this%ite,this%kts:this%kte,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%data_3d(var%grid%ite-this%halo_size,var%grid%kts,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northeast_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%west_in_win)
#endif
      endif
  endif
end subroutine

module subroutine put_northwest(this,var,do_dqdt)
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
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1:this%halo_size,1,1:this%halo_size)[this%northwest_neighbor] = var%dqdt_2d(this%its:this%its+this%halo_size,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%dqdt_2d(var%grid%its,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1:this%halo_size,1,1:this%halo_size)[this%northwest_neighbor] = var%data_2d(this%its:this%its+this%halo_size,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%data_2d(var%grid%its,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
#endif
      endif
  else
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%northwest_neighbor] = var%dqdt_3d(this%its:this%its+this%halo_size,this%kts:this%kte,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%south_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%northwest_neighbor] = var%data_3d(this%its:this%its+this%halo_size,this%kts:this%kte,this%jte-this%halo_size:this%jte)
#else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jte-this%halo_size), msg_size, &
            var%grid%corner_halo, this%northwest_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%south_in_win)
#endif
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
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,1,1:this%halo_size)[this%southwest_neighbor] = var%dqdt_2d(this%its:this%its+this%halo_size,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%dqdt_2d(var%grid%its,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,1,1:this%halo_size)[this%southwest_neighbor] = var%data_2d(this%its:this%its+this%halo_size,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%data_2d(var%grid%its,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
#endif
      endif
  else
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%southwest_neighbor] = var%dqdt_3d(this%its:this%its+this%halo_size,this%kts:this%kte,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%dqdt_3d(var%grid%its,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%east_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%southwest_neighbor] = var%data_3d(this%its:this%its+this%halo_size,this%kts:this%kte,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%data_3d(var%grid%its,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southwest_neighbor, disp, msg_size, var%grid%corner_EW_win_halo, this%east_in_win)
#endif
      endif
  endif
end subroutine

module subroutine put_southeast(this,var,do_dqdt)
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
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1:this%halo_size,1,1:this%halo_size)[this%southeast_neighbor] = var%dqdt_2d(this%ite-this%halo_size:this%ite,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%dqdt_2d(var%grid%ite-this%halo_size,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1:this%halo_size,1,1:this%halo_size)[this%southeast_neighbor] = var%data_2d(this%ite-this%halo_size:this%ite,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%data_2d(var%grid%ite-this%halo_size,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
#endif
      endif
  else
      if (dqdt) then
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%southeast_neighbor] = var%dqdt_3d(this%ite-this%halo_size:this%ite,this%kts:this%kte,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%dqdt_3d(var%grid%ite-this%halo_size,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
#endif
      else
#ifdef CRAY_PE
          !DIR$ PGAS DEFER_SYNC
          this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)[this%southeast_neighbor] = var%data_3d(this%ite-this%halo_size:this%ite,this%kts:this%kte,this%jts:this%jts+this%halo_size)
#else
          call MPI_Put(var%data_3d(var%grid%ite-this%halo_size,var%grid%kts,var%grid%jts), msg_size, &
            var%grid%corner_halo, this%southeast_neighbor, disp, msg_size, var%grid%corner_NS_win_halo, this%north_in_win)
#endif
      endif
  endif
end subroutine


module subroutine retrieve_northeast_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt

  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
      if (dqdt) then
          var%dqdt_2d(this%ite+1:this%ime,this%jte+1:this%jme) = this%east_in_3d(1:this%halo_size,1,1:this%halo_size)
      else
          var%data_2d(this%ite+1:this%ime,this%jte+1:this%jme) = this%east_in_3d(1:this%halo_size,1,1:this%halo_size)
      endif
  else
      if (dqdt) then
          var%dqdt_3d(this%ite+1:this%ime,this%kts:this%kte,this%jte+1:this%jme) = this%east_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(this%ite+1:this%ime,this%kts:this%kte,this%jte+1:this%jme) = this%east_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_northwest_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt

  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
      if (dqdt) then
          var%dqdt_2d(this%ims:this%its-1,this%jte+1:this%jme) = this%north_in_3d(1:this%halo_size,1,1:this%halo_size)
      else
          var%data_2d(this%ims:this%its-1,this%jte+1:this%jme) = this%north_in_3d(1:this%halo_size,1,1:this%halo_size)
      endif
  else
      if (dqdt) then
          var%dqdt_3d(this%ims:this%its-1,this%kts:this%kte,this%jte+1:this%jme) = this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(this%ims:this%its-1,this%kts:this%kte,this%jte+1:this%jme) = this%north_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
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
      if (dqdt) then
          var%dqdt_2d(this%ims:this%its-1,this%jms:this%jts-1) = this%west_in_3d(1:this%halo_size,1,1:this%halo_size)
      else
          var%data_2d(this%ims:this%its-1,this%jms:this%jts-1) = this%west_in_3d(1:this%halo_size,1,1:this%halo_size)
      endif
  else
      if (dqdt) then
          var%dqdt_3d(this%ims:this%its-1,this%kts:this%kte,this%jms:this%jts-1) = this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(this%ims:this%its-1,this%kts:this%kte,this%jms:this%jts-1) = this%west_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      endif
  endif
end subroutine

module subroutine retrieve_southeast_halo(this,var,do_dqdt)
  class(halo_t), intent(inout) :: this
  class(variable_t), intent(in) :: var
  logical, optional, intent(in) :: do_dqdt

  logical :: dqdt

  dqdt=.False.
  if (present(do_dqdt)) dqdt=do_dqdt
  
  if (var%two_d) then
      if (dqdt) then
          var%dqdt_2d(this%ite+1:this%ime,this%jms:this%jts-1) = this%south_in_3d(1:this%halo_size,1,1:this%halo_size)
      else
          var%data_2d(this%ite+1:this%ime,this%jms:this%jts-1) = this%south_in_3d(1:this%halo_size,1,1:this%halo_size)
      endif
  else
      if (dqdt) then
          var%dqdt_3d(this%ite+1:this%ime,this%kts:this%kte,this%jms:this%jts-1) = this%south_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
      else
          var%data_3d(this%ite+1:this%ime,this%kts:this%kte,this%jms:this%jts-1) = this%south_in_3d(1:this%halo_size,this%kts:this%kte,1:this%halo_size)
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
    
    ! Free the shared memory communicator
    call MPI_Comm_free(shared_comm, ierr)
    
    ! Set whether to use shared windows based on if any neighbors are shared
    this%use_shared_windows = this%north_shared .or. this%south_shared .or. &
                             this%east_shared .or. this%west_shared
end subroutine detect_shared_memory

end submodule
