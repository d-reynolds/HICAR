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
use iso_c_binding

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
    call MPI_Comm_Rank(comms,my_index)
    !Setup the parent-child group used for buffer communication
    call MPI_Comm_Group(comms,comp_proc)

    !The following if statements are structured as follows:
    ! 1. Check if the boundary is set
    ! 2. If not, set the neighbor index using process numbering for the computational group
    ! 3. Create a group for the neighbor
    !Need to make neighbor group
    if (.not.(this%south_boundary)) then
        this%south_neighbor = my_index - this%grid%ximages
        call MPI_Group_Incl(comp_proc, 1, [this%south_neighbor], this%south_neighbor_grp)
    endif
    if (.not.(this%west_boundary)) then
        this%west_neighbor  = my_index-1
        call MPI_Group_Incl(comp_proc, 1, [this%west_neighbor], this%west_neighbor_grp)
    endif
    if (.not.(this%east_boundary)) then
        this%east_neighbor  = my_index+1
        call MPI_Group_Incl(comp_proc, 1, [this%east_neighbor], this%east_neighbor_grp)
    endif
    if (.not.(this%north_boundary)) then
        this%north_neighbor = my_index + this%grid%ximages
        call MPI_Group_Incl(comp_proc, 1, [this%north_neighbor], this%north_neighbor_grp)
    endif
#endif


    !Compute diagonal direction neighbors
#ifdef CRAY_PE
    if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) this%southwest_neighbor = DOM_IMG_INDX(my_index - this%grid%ximages - 1)
    if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) this%northwest_neighbor = DOM_IMG_INDX(my_index + this%grid%ximages - 1)
    if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) this%southeast_neighbor  = DOM_IMG_INDX(my_index - this%grid%ximages + 1)
    if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) this%northeast_neighbor  = DOM_IMG_INDX(my_index + this%grid%ximages + 1)
#else
    if (.not.(this%south_boundary) .and. .not.(this%west_boundary)) this%southwest_neighbor = my_index - this%grid%ximages - 1
    if (.not.(this%north_boundary) .and. .not.(this%west_boundary)) this%northwest_neighbor = my_index + this%grid%ximages - 1
    if (.not.(this%south_boundary) .and. .not.(this%east_boundary)) this%southeast_neighbor  = my_index - this%grid%ximages + 1
    if (.not.(this%north_boundary) .and. .not.(this%east_boundary)) this%northeast_neighbor  = my_index + this%grid%ximages + 1
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
    
    !Now allocate the actual 3D halo
#ifdef CRAY_PE
allocate( this%south_in_3d( this%grid%ns_halo_nx+this%halo_size*2, this%grid%halo_nz, this%halo_size+1   )[*])
allocate( this%north_in_3d( this%grid%ns_halo_nx+this%halo_size*2, this%grid%halo_nz, this%halo_size        )[*])
allocate( this%east_in_3d(  this%halo_size       ,  this%grid%halo_nz, this%grid%ew_halo_ny+this%halo_size*2  )[*])
allocate( this%west_in_3d(  this%halo_size+1,  this%grid%halo_nz, this%grid%ew_halo_ny+this%halo_size*2  )[*])
#else
!We only want to set up remote windows for domain objects which are part of the actual domain
if (.not.(comms == MPI_COMM_NULL)) then

    nx = this%grid%ns_halo_nx+this%halo_size*2
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
    ny = this%grid%ew_halo_ny+this%halo_size*2
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

end subroutine

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
        endif

        ! if staggered in y direction, we need to carefully call the put and get commands
        if(var%ystag>0) then
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
        endif

        if((var%ystag+var%xstag)==0) then
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

    integer :: nx, ny, nz, n_2d, n_3d = 0
    type(c_ptr) :: tmp_ptr
    integer(KIND=MPI_ADDRESS_KIND) :: win_size, tmp_ptr2
    integer :: ierr, i, real_size
    type(MPI_Info) :: info_in 

    CALL MPI_Type_size(MPI_REAL, real_size)

    ! Loop over all adv_vars and count how many are 3D
    call adv_vars%reset_iterator()
    
    do while (adv_vars%has_more_elements())
        var = adv_vars%next()
        if (var%three_d) n_3d = n_3d + 1
    end do
    
    ! Loop over all exch vars and count how many are 3D
    call exch_vars%reset_iterator()
    
    do while (exch_vars%has_more_elements())
        var = exch_vars%next()
        if (var%three_d) n_3d = n_3d + 1
    end do
    if (STD_OUT_PE) write(*,*) "In Setup Batch Exch"
    if (STD_OUT_PE) flush(output_unit)

    ! Determine number of 2D and 3D vars present
    n_2d = (adv_vars%n_vars+exch_vars%n_vars)-n_3d

#ifdef CRAY_PE
    allocate(this%north_batch_in_3d(n_3d,1:(this%grid%ns_halo_nx+2),&
                    this%kms:this%kme,1:this%halo_size)[*])
    allocate(this%south_batch_in_3d(n_3d,1:(this%grid%ns_halo_nx+2),&
                    this%kms:this%kme,1:this%halo_size)[*])
    allocate(this%east_batch_in_3d(n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:(this%grid%ew_halo_ny+2))[*])
    allocate(this%west_batch_in_3d(n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:(this%grid%ew_halo_ny+2))[*])

#else
    if (.not.(comms == MPI_COMM_NULL)) then
            call MPI_Info_Create(info_in,ierr)
            call MPI_INFO_SET(info_in, 'no_locks', '.true.')
            call MPI_INFO_SET(info_in, 'same_size', '.true.')
            call MPI_INFO_SET(info_in, 'same_disp_unit', '.true.')
            call MPI_INFO_SET(info_in, 'alloc_shared_noncontig', '.true.')

            !First do NS
            nx = this%grid%ns_halo_nx+2
            nz = this%grid%halo_nz
            ny = this%halo_size
            win_size = nx*nz*ny*n_3d

            call MPI_Type_contiguous(nx*nz*ny*n_3d, MPI_REAL, this%N_3d_win_halo_type)
            call MPI_Type_commit(this%N_3d_win_halo_type)
            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%north_3d_win, ierr)
            call C_F_POINTER(tmp_ptr, this%north_batch_in_3d, [n_3d, nx, nz, ny])

            call MPI_Type_contiguous(nx*nz*ny*n_3d, MPI_REAL, this%S_3d_win_halo_type)
            call MPI_Type_commit(this%S_3d_win_halo_type)
            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%south_3d_win)
            call C_F_POINTER(tmp_ptr, this%south_batch_in_3d, [n_3d, nx, nz, ny])

            !Then do EW
            nx = this%halo_size
            nz = this%grid%halo_nz
            ny = this%grid%ew_halo_ny+2
            win_size = nx*nz*ny*n_3d
            call MPI_Type_contiguous(nx*nz*ny*n_3d, MPI_REAL, this%EW_3d_win_halo_type)
            call MPI_Type_commit(this%EW_3d_win_halo_type)

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%east_3d_win)
            call C_F_POINTER(tmp_ptr, this%east_batch_in_3d, [n_3d, nx, nz, ny])
            this%east_batch_in_3d = 1

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, info_in, comms, tmp_ptr, this%west_3d_win)
            call C_F_POINTER(tmp_ptr, this%west_batch_in_3d, [n_3d, nx, nz, ny])
            this%west_batch_in_3d = 1

            if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp,0,this%north_3d_win)
            if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp,0,this%south_3d_win)
            if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp,0,this%west_3d_win)
            if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp,0,this%east_3d_win)
    endif
#endif

    if (.not.(this%north_boundary)) allocate(this%north_buffer_3d(n_3d,1:(this%grid%ns_halo_nx+2),&
                    this%kms:this%kme,1:this%halo_size))
    if (.not.(this%south_boundary)) allocate(this%south_buffer_3d(n_3d,1:(this%grid%ns_halo_nx+2),&
                    this%kms:this%kme,1:this%halo_size))
    if (.not.(this%east_boundary)) allocate(this%east_buffer_3d(n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:(this%grid%ew_halo_ny+2)))
    if (.not.(this%west_boundary)) allocate(this%west_buffer_3d(n_3d,1:this%halo_size,&
                    this%kms:this%kme,1:(this%grid%ew_halo_ny+2)))


    ! If no 2D vars present, don't allocate arrays (nothing should be calling exch 2D then)
    if (n_2d > 0) then
#ifdef CRAY_PE
        allocate(this%north_batch_in_2d(n_2d,1:(this%grid%ns_halo_nx+2),1:this%halo_size)[*])
        allocate(this%south_batch_in_2d(n_2d,1:(this%grid%ns_halo_nx+2),1:this%halo_size)[*])
        allocate(this%east_batch_in_2d(n_2d,1:this%halo_size,1:(this%grid%ew_halo_ny+2))[*])
        allocate(this%west_batch_in_2d(n_2d,1:this%halo_size,1:(this%grid%ew_halo_ny+2))[*])
#else
        if (.not.(comms == MPI_COMM_NULL)) then

            nx = this%grid%ns_halo_nx+2
            ny = this%halo_size
            win_size = nx*ny*n_2d
            call MPI_Type_contiguous(nx*ny*n_2d, MPI_REAL, this%NS_2d_win_halo_type)
            call MPI_Type_commit(this%NS_2d_win_halo_type)

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%north_2d_win)
            call C_F_POINTER(tmp_ptr, this%north_batch_in_2d, [n_2d, nx, ny])
            this%north_batch_in_2d = 1

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%south_2d_win)
            call C_F_POINTER(tmp_ptr, this%south_batch_in_2d, [n_2d, nx, ny])
            this%south_batch_in_2d = 1

            nx = this%halo_size
            ny = this%grid%ew_halo_ny+2
            win_size = nx*ny*n_2d
            call MPI_Type_contiguous(nx*ny*n_2d, MPI_REAL, this%EW_2d_win_halo_type)
            call MPI_Type_commit(this%EW_2d_win_halo_type)

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%east_2d_win)
            call C_F_POINTER(tmp_ptr, this%east_batch_in_2d, [n_2d, nx, ny])
            this%east_batch_in_2d = 1

            call MPI_WIN_ALLOCATE(win_size*real_size, real_size, MPI_INFO_NULL, comms, tmp_ptr, this%west_2d_win)
            call C_F_POINTER(tmp_ptr, this%west_batch_in_2d, [n_2d, nx, ny])
            this%west_batch_in_2d = 1

            if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp,0,this%south_2d_win)
            if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp,0,this%north_2d_win)
            if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp,0,this%east_2d_win)
            if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp,0,this%west_2d_win)        
        endif
#endif

        if (.not.(this%north_boundary)) allocate(this%north_buffer_2d(n_2d,1:(this%grid%ns_halo_nx+2),1:this%halo_size))
        if (.not.(this%south_boundary)) allocate(this%south_buffer_2d(n_2d,1:(this%grid%ns_halo_nx+2),1:this%halo_size))
        if (.not.(this%east_boundary)) allocate(this%east_buffer_2d(n_2d,1:this%halo_size,1:(this%grid%ew_halo_ny+2)))
        if (.not.(this%west_boundary)) allocate(this%west_buffer_2d(n_2d,1:this%halo_size,1:(this%grid%ew_halo_ny+2)))
    endif

end subroutine setup_batch_exch



module subroutine halo_3d_send_batch(this, exch_vars, adv_vars,exch_var_only)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    logical, optional, intent(in) :: exch_var_only
    
    type(variable_t) :: var
    logical :: exch_v_only
    integer :: n, k_max, msg_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    msg_size = 1
    disp = 0

    exch_v_only = .False.
    if (present(exch_var_only)) exch_v_only=exch_var_only

    call adv_vars%reset_iterator()
    call exch_vars%reset_iterator()
    n = 1
    ! Now iterate through the dictionary as long as there are more elements present
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
        call MPI_Win_Start(this%south_neighbor_grp,0,this%north_3d_win)
        call MPI_Put(this%south_buffer_3d, size(this%south_buffer_3d), &
            MPI_REAL, this%south_neighbor, disp, size(this%south_buffer_3d), MPI_REAL, this%north_3d_win)
        !call MPI_Put(this%south_buffer_3d, msg_size, &
        !    this%N_3d_win_halo_type, this%south_neighbor, disp, msg_size, this%N_3d_win_halo_type, this%north_3d_win)
        call MPI_Win_Complete(this%north_3d_win)
#endif        
    endif
    if (.not.(this%north_boundary)) then
#ifdef CRAY_PE
        !DIR$ PGAS DEFER_SYNC
        this%south_batch_in_3d(:,:,:,:)[this%north_neighbor] = this%north_buffer_3d(:,:,:,:)
#else
        call MPI_Win_Start(this%north_neighbor_grp,0,this%south_3d_win)
        call MPI_Put(this%north_buffer_3d, size(this%north_buffer_3d), &
            MPI_REAL, this%north_neighbor, disp, size(this%north_buffer_3d), MPI_REAL, this%south_3d_win)
        !call MPI_Put(this%north_buffer_3d, msg_size, &
        !    this%S_3d_win_halo_type, this%north_neighbor, disp, msg_size, this%S_3d_win_halo_type, this%south_3d_win)
        call MPI_Win_Complete(this%south_3d_win)

#endif
    endif
    if (.not.(this%east_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%west_batch_in_3d(:,:,:,:)[this%east_neighbor] = this%east_buffer_3d(:,:,:,:)
#else
        call MPI_Win_Start(this%east_neighbor_grp,0,this%west_3d_win)
        call MPI_Put(this%east_buffer_3d, size(this%east_buffer_3d), &
            MPI_REAL, this%east_neighbor, disp, size(this%east_buffer_3d), MPI_REAL, this%west_3d_win)
        !call MPI_Put(this%east_buffer_3d, msg_size, &
        !    this%EW_3d_win_halo_type, this%east_neighbor, disp, msg_size, this%EW_3d_win_halo_type, this%west_3d_win)
        call MPI_Win_Complete(this%west_3d_win)
#endif  
    endif
    if (.not.(this%west_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%east_batch_in_3d(:,:,:,:)[this%west_neighbor] = this%west_buffer_3d(:,:,:,:)
#else
        call MPI_Win_Start(this%west_neighbor_grp,0,this%east_3d_win)
        call MPI_Put(this%west_buffer_3d, size(this%west_buffer_3d), &
            MPI_REAL, this%west_neighbor, disp, size(this%west_buffer_3d), MPI_REAL, this%east_3d_win)
        !call MPI_Put(this%west_buffer_3d, msg_size, &
        !    this%EW_3d_win_halo_type, this%west_neighbor, disp, msg_size, this%EW_3d_win_halo_type, this%east_3d_win)
        call MPI_Win_Complete(this%east_3d_win)
#endif        
    endif

end subroutine halo_3d_send_batch

module subroutine halo_3d_retrieve_batch(this,exch_vars, adv_vars,exch_var_only, wait_timer)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    logical, optional, intent(in) :: exch_var_only
    type(timer_t), optional,     intent(inout)   :: wait_timer

    type(variable_t) :: var
    integer :: n, k_max
    logical :: exch_v_only

    exch_v_only = .False.
    if (present(exch_var_only)) exch_v_only=exch_var_only
    if (present(wait_timer)) call wait_timer%start()

#ifdef CRAY_PE        
    sync images( this%neighbors )
#else
    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Wait(this%west_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Wait(this%east_3d_win)
#endif
    if (present(wait_timer)) call wait_timer%stop()

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

#ifndef CRAY_PE
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp,MPI_MODE_NOSTORE,this%north_3d_win)
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp,MPI_MODE_NOSTORE,this%south_3d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp,MPI_MODE_NOSTORE,this%west_3d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp,MPI_MODE_NOSTORE,this%east_3d_win)
#endif

end subroutine halo_3d_retrieve_batch

module subroutine halo_2d_send_batch(this, exch_vars, adv_vars)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(variable_t) :: var
    integer :: n, msg_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

    msg_size = 1
    disp = 0

    call exch_vars%reset_iterator()
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
        call MPI_Win_Start(this%north_neighbor_grp,0,this%south_2d_win)
        call MPI_Put(this%north_buffer_2d, size(this%north_buffer_2d), &
            MPI_REAL, this%north_neighbor, disp, size(this%north_buffer_2d), MPI_REAL, this%south_2d_win)
        call MPI_Win_Complete(this%south_2d_win)
#endif
    endif
    if (.not.(this%south_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%north_batch_in_2d(:,:,:)[this%south_neighbor] = this%south_buffer_2d(:,:,:)
#else
        call MPI_Win_Start(this%south_neighbor_grp,0,this%north_2d_win)
        call MPI_Put(this%south_buffer_2d, size(this%south_buffer_2d), &
            MPI_REAL, this%south_neighbor, disp, size(this%south_buffer_2d), MPI_REAL, this%north_2d_win)
        call MPI_Win_Complete(this%north_2d_win)
#endif        
    endif
    if (.not.(this%east_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%west_batch_in_2d(:,:,:)[this%east_neighbor] = this%east_buffer_2d(:,:,:)
#else
        call MPI_Win_Start(this%east_neighbor_grp,0,this%west_2d_win)
        call MPI_Put(this%east_buffer_2d, size(this%east_buffer_2d), &
            MPI_REAL, this%east_neighbor, disp, size(this%east_buffer_2d), MPI_REAL, this%west_2d_win)
        call MPI_Win_Complete(this%west_2d_win)
#endif  
    endif
    if (.not.(this%west_boundary)) then
#ifdef CRAY_PE        
        !DIR$ PGAS DEFER_SYNC
        this%east_batch_in_2d(:,:,:)[this%west_neighbor] = this%west_buffer_2d(:,:,:)
#else
        call MPI_Win_Start(this%west_neighbor_grp,0,this%east_2d_win)
        call MPI_Put(this%west_buffer_2d, size(this%west_buffer_2d), &
            MPI_REAL, this%west_neighbor, disp, size(this%west_buffer_2d), MPI_REAL, this%east_2d_win)
        call MPI_Win_Complete(this%east_2d_win)
#endif        
    endif

end subroutine halo_2d_send_batch

module subroutine halo_2d_retrieve_batch(this, exch_vars, adv_vars)
    class(halo_t), intent(inout) :: this
    class(var_dict_t), intent(inout) :: exch_vars, adv_vars
    type(variable_t) :: var
    integer :: n

#ifdef CRAY_PE        
    sync images( this%neighbors )
#else
    if (.not.(this%south_boundary)) call MPI_Win_Wait(this%south_2d_win)
    if (.not.(this%north_boundary)) call MPI_Win_Wait(this%north_2d_win)
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

#ifndef CRAY_PE
    if (.not.(this%south_boundary)) call MPI_Win_Post(this%south_neighbor_grp,MPI_MODE_NOSTORE,this%south_2d_win)
    if (.not.(this%north_boundary)) call MPI_Win_Post(this%north_neighbor_grp,MPI_MODE_NOSTORE,this%north_2d_win)
    if (.not.(this%east_boundary)) call MPI_Win_Post(this%east_neighbor_grp,MPI_MODE_NOSTORE,this%east_2d_win)
    if (.not.(this%west_boundary)) call MPI_Win_Post(this%west_neighbor_grp,MPI_MODE_NOSTORE,this%west_2d_win)
#endif

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
  
  offs=var%ystag
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
  
  offs=var%xstag
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
  
  offs=var%xstag
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
  
  offs=var%ystag
  
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
  
  offs=var%ystag
 
  
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
  
  offs=var%xstag

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
  
  offs=var%xstag
  
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

end submodule
