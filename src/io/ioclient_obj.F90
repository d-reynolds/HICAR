
!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!----------------------------------------------------------
submodule(ioclient_interface) ioclient_implementation
  use debug_module,             only : check_ncdf, check_var
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use output_metadata,          only : get_varindx, get_varmeta
  use meta_data_interface,      only : meta_data_t
  use geo,                      only : geo_interp2d
  use string,           only  : str

  implicit none

contains        
    
    module subroutine init_ioclient(this, domain, forcing, options, n_indx)
        implicit none
        class(ioclient_t),  intent(inout)  :: this
        type(domain_t),     intent(inout)  :: domain
        type(boundary_t),   intent(in)     :: forcing
        type(options_t),    intent(in)     :: options(:)
        integer,            intent(in)     :: n_indx

        type(variable_t) :: var
        type(MPI_Group) :: family_group

        integer :: my_rank, comm_size, some_child_id
        
        this%i_s_r = forcing%its; this%i_e_r = forcing%ite
        this%k_s_r = forcing%kts; this%k_e_r = forcing%kte
        this%j_s_r = forcing%jts; this%j_e_r = forcing%jte
        this%written = .False.
        this%nest_updated = .False.

        this%i_s_w = domain%its; this%i_e_w = domain%ite
        this%k_s_w = domain%kts; this%k_e_w = domain%kte+1
        this%j_s_w = domain%jts; this%j_e_w = domain%jte

        this%i_s_re = domain%ims; this%i_e_re = domain%ime
        this%j_s_re = domain%jms; this%j_e_re = domain%jme

        this%ide = domain%ide; this%kde = domain%kde; this%jde = domain%jde

        if (domain%ims == domain%ids) this%i_s_w = domain%ids
        if (domain%ime == domain%ide) this%i_e_w = domain%ide
        if (domain%jms == domain%jds) this%j_s_w = domain%jds
        if (domain%jme == domain%jde) this%j_e_w = domain%jde


        ! If this run uses nests, then set the variables which we need to output to a child nest.
        ! Microphysics must be the same for all nests, so any child's forcing_options may be used;
        ! pick the first child nest id.
        if (size(options(n_indx)%general%child_nests) > 0) then
            some_child_id = options(n_indx)%general%child_nests(1)
            this%vars_for_nest = options(some_child_id)%forcing%vars_to_read
            call classify_nest_init_vars(options(some_child_id), this%vars_for_nest, this%send_init_vars)
        endif

        ! If this is a child nest, compute the 2D + 3D restart vars to receive from parent.
        if (options(n_indx)%general%parent_nest > 0) then
            call classify_nest_init_vars(options(n_indx), options(n_indx)%forcing%vars_to_read, this%recv_init_vars)
        endif

        ! Change 2: Store output/restart classification and compute counts
        this%vars_for_output = options(n_indx)%output%vars_for_output
        this%vars_for_restart = options(n_indx)%vars_for_restart
        this%restart_count = options(n_indx)%restart%restart_count
        this%restart_counter = 1
        if (options(n_indx)%restart%restart) this%restart_counter = 2

        call init_with_server(this)

        call setup_MPI_windows(this)

        !Setup the parent-child group used for buffer communication
        call MPI_Comm_Group(this%parent_comms,family_group)
        call MPI_Comm_size(this%parent_comms, comm_size)

        !Our rank is the last process in the group, which is equal to n_children
        call MPI_Group_Incl(family_group, 1, [comm_size-1], this%parent_group)

        ! Do MPI_Win_Post on read_buffer to indicate that we are open for delivery of input data
        call MPI_Win_Post(this%parent_group,0,this%read_win)

    end subroutine init_ioclient

    subroutine init_with_server(this)
        implicit none
        class(ioclient_t),   intent(inout) :: this

        integer :: ierr, PE_parent_comm

        ! The PE of the parent communicator is the last PE in the communicator
        call MPI_Comm_Size(this%parent_comms,PE_parent_comm)
        PE_parent_comm = PE_parent_comm - 1

        call MPI_Gatherv(this%i_s_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%i_e_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%k_s_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%k_e_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_s_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_e_w, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%i_s_re, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%i_e_re, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_s_re, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_e_re, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%i_s_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%i_e_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%k_s_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%k_e_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_s_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)
        call MPI_Gatherv(this%j_e_r, 1, MPI_INTEGER, 0, [0], [0], &
            MPI_INTEGER, PE_parent_comm, this%parent_comms)

        call MPI_Allreduce(MPI_IN_PLACE,this%ide,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%kde,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%jde,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
    
    end subroutine init_with_server


    subroutine setup_MPI_windows(this)
        class(ioclient_t),   intent(inout)  :: this

        type(c_ptr) :: tmp_ptr
        integer(KIND=MPI_ADDRESS_KIND) :: win_size
        integer :: nx_w, nz_w, ny_w, nx_re, ny_re, n_w_2d, n_f, n_f_2d, n_f_3d_init, n_w_3d, ierr
        integer :: nx_r, nz_r, ny_r, n_r, real_size

        CALL MPI_Type_size(MPI_REAL, real_size)

        nx_w = 0
        nz_w = 0
        ny_w = 0
        n_w_3d = 0
        n_w_2d = 0
        n_f = 0
        n_f_2d = 0
        n_f_3d_init = 0
        nx_re = 0
        ny_re = 0

        nx_r = 0
        nz_r = 0
        ny_r = 0
        n_r = 0

       ! Setup MPI windows for inter-process communication        
        call MPI_Allreduce(MPI_IN_PLACE,nx_w,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,ny_w,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,nz_w,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,nx_re,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,ny_re,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,nx_r,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,ny_r,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,nz_r,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,n_w_3d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,n_w_2d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,n_f,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,n_f_2d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,n_f_3d_init,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%nz_init_3d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,n_r,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        ! Shared memory windows for all builds; GPU builds get device copies for kernel packing
        win_size = n_w_3d*nx_re*nz_w*ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%write_win_3d)
        call C_F_POINTER(tmp_ptr, this%write_buffer_3d, [n_w_3d, nx_re, nz_w, ny_re])
        this%write_buffer_3d = kEMPT_BUFF
        !$acc enter data copyin(this%write_buffer_3d)

        win_size = n_w_2d*nx_re*ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%write_win_2d)
        call C_F_POINTER(tmp_ptr, this%write_buffer_2d, [n_w_2d, nx_re, ny_re])
        this%write_buffer_2d = kEMPT_BUFF
        !$acc enter data copyin(this%write_buffer_2d)

        if (n_f > 0) then
            win_size = n_f*nx_w*ny_w*nz_w
            call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%forcing_win)
            call C_F_POINTER(tmp_ptr, this%forcing_buffer, [n_f, nx_w, nz_w, ny_w])
            this%forcing_buffer = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer)
        endif

        if (n_f_2d > 0) then
            win_size = n_f_2d*nx_w*ny_w
            call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%forcing_win_2d)
            call C_F_POINTER(tmp_ptr, this%forcing_buffer_2d, [n_f_2d, nx_w, ny_w])
            this%forcing_buffer_2d = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer_2d)
        endif

        if (n_f_3d_init > 0) then
            win_size = n_f_3d_init*nx_w*nz_w*ny_w
            call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%forcing_win_3d_init)
            call C_F_POINTER(tmp_ptr, this%forcing_buffer_3d_init, [n_f_3d_init, nx_w, nz_w, ny_w])
            this%forcing_buffer_3d_init = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer_3d_init)
        endif

        win_size = n_r*nx_r*nz_r*ny_r
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%read_win)
        call C_F_POINTER(tmp_ptr, this%read_buffer, [n_r, nx_r, nz_r, ny_r])
        this%read_buffer = kEMPT_BUFF
        !$acc enter data copyin(this%read_buffer)

    end subroutine setup_MPI_windows


    ! This subroutine pushes the output fields from the domain object
    ! to the write buffer for the IO processes to use
    module subroutine push(this, domain)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),   intent(inout)    :: domain
        
        type(variable_t) :: var
        type(meta_data_t) :: tmp_var
        integer :: i, n_3d, n_2d, nx, ny, nz_v, i_s_w, i_e_w, j_s_w, j_e_w, idx
        logical :: should_do_restart
        integer :: ii, jj, kk
        integer :: comm_size, my_rank, ierr

        n_3d = 1
        n_2d = 1

        !This is false only when it is the first call to push (i.e. first write call)
        if (this%written) then
            ! Do MPI_Win_Wait on write_buffer to make sure that server process has completed writing of previous data
            call smart_wait(this%write_win_3d, 'Waiting for write_win_3d completion')
            call smart_wait(this%write_win_2d, 'Waiting for write_win_2d completion')
        endif
        this%written = .False.

        ! Change 2: Determine if restart variables should be packed this step
        should_do_restart = (this%restart_counter > this%restart_count) .or. this%first_push

        ! Pass 1: output variables (always packed)
        do i = 1, kMAX_STORAGE_VARS
            ! get the next variable in the structure
            if (domain%vars_to_out(i)%v <= 0) cycle
            if (this%vars_for_output(i) <= 0) cycle  ! skip non-output vars
            idx = domain%vars_to_out(i)%v
            tmp_var = get_varmeta(i)

            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            if (domain%ime == domain%ide) i_e_w = i_e_w+tmp_var%xstag
            if (domain%jme == domain%jde) j_e_w = j_e_w+tmp_var%ystag
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            ! GPU path: pack on device, no GPU->CPU transfer
            if (tmp_var%three_d) then
                nz_v = domain%vars_3d(idx)%dim_len(2)
                associate(src => domain%vars_3d(idx)%data_3d, dst => this%write_buffer_3d)
                !$acc parallel loop gang vector collapse(3) present(dst, src)
                do jj = 1, ny
                  do kk = 1, nz_v
                    do ii = 1, nx
                      dst(n_3d, ii, kk, jj) = src(i_s_w+ii-1, kk, j_s_w+jj-1)
                    enddo
                  enddo
                enddo
                end associate
                n_3d = n_3d + 1
            else if (tmp_var%two_d) then
                if (domain%vars_2d(idx)%dtype == kREAL) then
                    associate(src => domain%vars_2d(idx)%data_2d, dst => this%write_buffer_2d)
                    !$acc parallel loop gang vector collapse(2) present(dst, src)
                    do jj = 1, ny
                      do ii = 1, nx
                        dst(n_2d, ii, jj) = src(i_s_w+ii-1, j_s_w+jj-1)
                      enddo
                    enddo
                    end associate
                elseif (domain%vars_2d(idx)%dtype == kINTEGER) then
                    associate(src => domain%vars_2d(idx)%data_2di, dst => this%write_buffer_2d)
                    !$acc parallel loop gang vector collapse(2) present(dst, src)
                    do jj = 1, ny
                      do ii = 1, nx
                        dst(n_2d, ii, jj) = real(src(i_s_w+ii-1, j_s_w+jj-1))
                      enddo
                    enddo
                    end associate
                endif
                n_2d = n_2d + 1
            else
                write(*,*) 'Error: Variable ', tmp_var%name, ' not found in parent domain: ', domain%nest_indx
                stop
            endif
        enddo

        ! Pass 2: restart-only variables (only on restart steps or first push)
        if (should_do_restart) then
            ! Send dt to ioserver for restart file (rank 0 only; all compute ranks have same dt)
            call MPI_Comm_Rank(this%parent_comms, my_rank)
            if (my_rank == 0) then
                call MPI_Comm_Size(this%parent_comms, comm_size)
                call MPI_Send(domain%dt, 1, MPI_REAL, comm_size-1, 42, this%parent_comms, ierr)
            endif

            do i = 1, kMAX_STORAGE_VARS
                if (domain%vars_to_out(i)%v <= 0) cycle
                if (this%vars_for_output(i) > 0) cycle      ! already packed in pass 1
                if (this%vars_for_restart(i) <= 0) cycle     ! not a restart var
                idx = domain%vars_to_out(i)%v
                tmp_var = get_varmeta(i)

                i_s_w = this%i_s_w; i_e_w = this%i_e_w
                j_s_w = this%j_s_w; j_e_w = this%j_e_w
                if (domain%ime == domain%ide) i_e_w = i_e_w+tmp_var%xstag
                if (domain%jme == domain%jde) j_e_w = j_e_w+tmp_var%ystag
                nx = i_e_w - i_s_w + 1
                ny = j_e_w - j_s_w + 1

                if (tmp_var%three_d) then
                    nz_v = domain%vars_3d(idx)%dim_len(2)
                    associate(src => domain%vars_3d(idx)%data_3d, dst => this%write_buffer_3d)
                    !$acc parallel loop gang vector collapse(3) present(dst, src)
                    do jj = 1, ny
                      do kk = 1, nz_v
                        do ii = 1, nx
                          dst(n_3d, ii, kk, jj) = src(i_s_w+ii-1, kk, j_s_w+jj-1)
                        enddo
                      enddo
                    enddo
                    end associate
                    n_3d = n_3d + 1
                else if (tmp_var%two_d) then
                    if (domain%vars_2d(idx)%dtype == kREAL) then
                        associate(src => domain%vars_2d(idx)%data_2d, dst => this%write_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(dst, src)
                        do jj = 1, ny
                          do ii = 1, nx
                            dst(n_2d, ii, jj) = src(i_s_w+ii-1, j_s_w+jj-1)
                          enddo
                        enddo
                        end associate
                    elseif (domain%vars_2d(idx)%dtype == kINTEGER) then
                        associate(src => domain%vars_2d(idx)%data_2di, dst => this%write_buffer_2d)
                        !$acc parallel loop gang vector collapse(2) present(dst, src)
                        do jj = 1, ny
                          do ii = 1, nx
                            dst(n_2d, ii, jj) = real(src(i_s_w+ii-1, j_s_w+jj-1))
                          enddo
                        enddo
                        end associate
                    endif
                    n_2d = n_2d + 1
                endif
            enddo
        endif

        !$acc wait
        !$acc update host(this%write_buffer_3d, this%write_buffer_2d)

        this%written = .True.
        call domain%increment_output_time()

        ! Update restart counter (mirror server logic)
        if (this%restart_counter > this%restart_count) this%restart_counter = 1
        this%restart_counter = this%restart_counter + 1
        if (this%first_push) this%first_push = .false.

        ! Do MPI_Win_Post on write_win to inform that server process can begin writing of data
        call MPI_Win_Post(this%parent_group,0,this%write_win_3d)
        call MPI_Win_Post(this%parent_group,0,this%write_win_2d)


    end subroutine 

    module subroutine update_nest(this, domain)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),      intent(in)    :: domain

        type(variable_t) :: var
        type(meta_data_t) :: tmp_var
        integer :: i, n_3d, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w, var_indx, idx
        logical :: var_val_check
        character(len=kMAX_NAME_LENGTH) :: err_msg
        integer :: ii, jj, kk, nz_v

        n_3d = 1
        ! Do MPI_Win_Wait on forcing_win. This is mostly unnecesarry, since the server process is what will be waiting on us, which is handeled by the MPI_Win_Start call
        ! in ioserver%gather_foring. Still, MPI_Win_Wait is called here for completeness of PSCW model
        if (this%nest_updated) then
            call smart_wait(this%forcing_win, 'Waiting for forcing_win completion')
        endif

        this%nest_updated = .False.
        !This is false only when it is the first call to push (i.e. first write call)

        do i = 1, kMAX_STORAGE_VARS
            ! get the next variable in the structure
            if (this%vars_for_nest(i) == '') cycle
            var_indx = get_varindx(this%vars_for_nest(i))
            tmp_var = get_varmeta(var_indx)
            var_val_check = (var%maxval /= kUNSET_REAL .and. var%minval /= kUNSET_REAL)

            if (tmp_var%two_d) cycle

            if (domain%var_indx(var_indx)%v <= 0) then
                write(*,*) 'Error: Variable ', this%vars_for_nest(i), ' not found in parent domain: ', domain%nest_indx
                stop
            end if

            idx = domain%var_indx(var_indx)%v

            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            var = domain%vars_3d(idx)
            i_e_w = i_e_w+var%xstag !Add extra to accomodate staggered vars
            j_e_w = j_e_w+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            nz_v = var%dim_len(2)
            associate(src => domain%vars_3d(idx)%data_3d, dst => this%forcing_buffer)
            !$acc parallel loop gang vector collapse(3) present(dst, src)
            do jj = 1, ny
              do kk = 1, nz_v
                do ii = 1, nx
                  dst(n_3d, ii, kk, jj) = src(i_s_w+ii-1, kk, j_s_w+jj-1)
                enddo
              enddo
            enddo
            end associate

            if (var_val_check) then
                err_msg = 'Warning on ioclient_obj::update_nest: Nest level: '// str(domain%nest_indx)
                call check_var(var, trim(err_msg))
            endif

            n_3d = n_3d+1
        enddo

        !$acc wait
        !$acc update host(this%forcing_buffer)

        this%nest_updated = .True.
        ! Do MPI_Win_Post on forcing_win to inform that server process can begin gathering of nest data
        call MPI_Win_Post(this%parent_group,0,this%forcing_win)

        ! One-time init: pack 2D + extra 3D restart vars for child nest initialization
        if (.not. this%initial_nest_done) then
            if (associated(this%forcing_buffer_2d)) then
                call pack_init_vars_2d(this, domain, this%send_init_vars%vars_2d)
                !$acc wait
                !$acc update host(this%forcing_buffer_2d)
                call MPI_Win_Post(this%parent_group, 0, this%forcing_win_2d)
            endif

            if (associated(this%forcing_buffer_3d_init)) then
                call pack_init_vars_3d(this, domain, this%send_init_vars%vars_3d)
                !$acc wait
                !$acc update host(this%forcing_buffer_3d_init)
                call MPI_Win_Post(this%parent_group, 0, this%forcing_win_3d_init)
            endif

            this%initial_nest_done = .true.
        endif

    end subroutine

    ! This subroutine receives the input fields from the IO buffer
    ! for assignment to the forcing object
    module subroutine receive(this, forcing, domain)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(boundary_t), intent(inout)  :: forcing
        type(domain_t),   intent(inout)  :: domain

        integer :: ii, jj, kk, nz_v
        type(variable_t)     :: var
        integer :: i, n, nx, ny, var_id
        logical :: var_val_check
        character(len=kMAX_NAME_LENGTH) :: err_msg
                
        n = 1
        ! Do MPI_Win_Wait on read_buffer to make sure that server process has completed data transfer
        call smart_wait(this%read_win, 'Waiting for read_win completion')

        !$acc update device(this%read_buffer)

        do i = 1, forcing%variables%n_vars
            if (forcing%variables%var_list(i)%var%computed) cycle
            var_val_check = (forcing%variables%var_list(i)%var%maxval /= kUNSET_REAL .and. forcing%variables%var_list(i)%var%minval /= kUNSET_REAL)
            if (forcing%variables%var_list(i)%var%two_d .and. &
                forcing%variables%var_list(i)%var%dtype == kREAL) then
                nx = size(forcing%variables%var_list(i)%var%data_2d, 1)
                ny = size(forcing%variables%var_list(i)%var%data_2d, 2)

                associate(dst => forcing%variables%var_list(i)%var%data_2d, &
                          src => this%read_buffer)
                !$acc parallel loop gang vector collapse(2) present(dst, src)
                do jj = 1, ny
                    do ii = 1, nx
                        dst(ii, jj) = src(n, ii, 1, jj)
                    enddo
                enddo
                end associate
            else if (forcing%variables%var_list(i)%var%three_d) then
                nx = size(forcing%variables%var_list(i)%var%data_3d, 1)
                ny = size(forcing%variables%var_list(i)%var%data_3d, 3)
                nz_v = forcing%variables%var_list(i)%var%dim_len(2)

                associate(dst => forcing%variables%var_list(i)%var%data_3d, &
                          src => this%read_buffer)
                !$acc parallel loop gang vector collapse(3) present(dst, src)
                do jj = 1, ny
                    do kk = 1, nz_v
                        do ii = 1, nx
                            dst(ii, kk, jj) = src(n, ii, kk, jj)
                        enddo
                    enddo
                enddo
                end associate
            endif
            if (var_val_check) then
                err_msg = 'Warning on ioclient_obj::receive: Nest level: '// str(domain%nest_indx)
                call check_var(forcing%variables%var_list(i)%var, trim(err_msg))
            endif

            n = n + 1
        enddo
        !$acc wait
        ! Do MPI_Win_Post on read_buffer to indicate that we are open for delivery of new input data
        call MPI_Win_Post(this%parent_group,0,this%read_win)


    end subroutine 
    
    ! This subroutine assigns the data from the write_buffer to the appropriate fields
    ! of the domain object. It is a "reverse write"
    module subroutine receive_rst(this, domain, options)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(domain_t),   intent(inout)  :: domain
        type(options_t),  intent(in)     :: options

        type(meta_data_t)     :: tmp_var
        type(variable_t)      :: var
        integer :: i, n_2d, n_3d, nx, ny, i_s_re, i_e_re, j_s_re, j_e_re
        character(len=kMAX_NAME_LENGTH) :: varname

        ! Because this is for reading restart data, performance is not critical, and 
        ! we use a simple MPI_fence syncronization
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)

        n_3d = 1
        n_2d = 1
        
        do i = 1, kMAX_STORAGE_VARS

            !See if var is in restart vars
            if (options%vars_for_restart(i) <= 0) then
                cycle
            endif

            ! get the next variable in the structure
            if (domain%vars_to_out(i)%v <= 0) cycle
            tmp_var = get_varmeta(i)
            if (tmp_var%two_d) then
                var = domain%vars_2d(domain%vars_to_out(i)%v)
            else if (tmp_var%three_d) then
                var = domain%vars_3d(domain%vars_to_out(i)%v)
            else
                write(*,*) 'Error: Variable ', tmp_var%name, ' not found in parent domain: ', domain%nest_indx
                stop
            endif


            i_s_re = this%i_s_re; i_e_re = this%i_e_re
            j_s_re = this%j_s_re; j_e_re = this%j_e_re
            i_e_re = i_e_re+var%xstag !Add extra to accomodate staggered vars
            j_e_re = j_e_re+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_re - i_s_re + 1
            ny = j_e_re - j_s_re + 1

            if (var%three_d) then
                domain%vars_3d(domain%vars_to_out(i)%v)%data_3d(i_s_re:i_e_re,1:var%dim_len(2),j_s_re:j_e_re) = &
                    this%write_buffer_3d(n_3d,1:nx,1:var%dim_len(2),1:ny)
                !$acc update device(domain%vars_3d(domain%vars_to_out(i)%v)%data_3d)
                n_3d = n_3d+1
            else
                if (var%dtype == kREAL) then
                    domain%vars_2d(domain%vars_to_out(i)%v)%data_2d(i_s_re:i_e_re,j_s_re:j_e_re) = &
                        this%write_buffer_2d(n_2d,1:nx,1:ny)
                    !$acc update device(domain%vars_2d(domain%vars_to_out(i)%v)%data_2d)
                elseif (var%dtype == kINTEGER) then
                    domain%vars_2d(domain%vars_to_out(i)%v)%data_2di(i_s_re:i_e_re,j_s_re:j_e_re) = &
                        int(this%write_buffer_2d(n_2d,1:nx,1:ny))
                    !$acc update device(domain%vars_2d(domain%vars_to_out(i)%v)%data_2di)
                ! elseif (var%dtype == kDOUBLE) then
                !     var%data_2dd(i_s_re:i_e_re,j_s_re:j_e_re) = &
                !             dble(this%write_buffer_2d(n_2d,1:nx,1:ny))
                endif                
                n_2d = n_2d+1
            endif
        enddo

        !$acc wait  ! Ensure all restart data uploads complete before proceeding

        ! Pack this%write_win_2d and this%write_win_3d with kEMPT_BUFF. This is done so that the write_buffers are
        ! again initialized to empty buffers, which is needed for the output objects on the IOserver processes
        ! to find the "blocks" which are needed to output
        this%write_buffer_2d = kEMPT_BUFF
        this%write_buffer_3d = kEMPT_BUFF

        !$acc update device(this%write_buffer_2d, this%write_buffer_3d)

    end subroutine


    ! Receive initial 2D+3D restart state from parent nest via temporary MPI_Recv.
    ! Each variable is geo-interpolated from forcing grid to domain grid using the
    ! geolut that was set up in get_initial_conditions.
    module subroutine receive_nest_init(this, domain, forcing)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(inout) :: domain
        type(boundary_t),  intent(in)    :: forcing

        integer :: n_init_2d, n_init_3d, comm_size, parent_rank

        n_init_2d = count(this%recv_init_vars%vars_2d /= '')
        n_init_3d = count(this%recv_init_vars%vars_3d /= '')

        if (n_init_2d == 0 .and. n_init_3d == 0) return

        call MPI_Comm_Size(this%parent_comms, comm_size)
        parent_rank = comm_size - 1

        if (n_init_2d > 0) then
            call unpack_init_vars_2d(this, domain, forcing, this%recv_init_vars%vars_2d, n_init_2d, parent_rank)
        endif

        if (n_init_3d > 0) then
            call unpack_init_vars_3d(this, domain, forcing, this%recv_init_vars%vars_3d, n_init_3d, parent_rank)
        endif

    end subroutine


    ! Filter/classify kMAX_STORAGE_VARS into 2D and 3D name lists for the init-only
    ! nest transfer. Used for both the parent-send path (opts = child's options,
    ! forcing_vars = child's forcing var list) and the child-receive path
    ! (opts = this node's options, forcing_vars = this node's forcing var list).
    subroutine classify_nest_init_vars(opts, forcing_vars, out)
        type(options_t),             intent(in)    :: opts
        character(len=*),            intent(in)    :: forcing_vars(:)
        type(nest_init_var_list_t),  intent(out)   :: out

        type(meta_data_t) :: meta
        integer :: v, j_2d, j_3d

        out%vars_2d = ''
        out%vars_3d = ''
        j_2d = 1; j_3d = 1
        do v = 1, kMAX_STORAGE_VARS
            if (opts%vars_to_allocate(v) <= 0) cycle
            if (opts%vars_for_restart(v) <= 0) cycle
            meta = get_varmeta(v)
            if (len_trim(meta%name) == 0) cycle
            if (any(forcing_vars == meta%name)) cycle
            if (meta%two_d) then
                out%vars_2d(j_2d) = meta%name
                j_2d = j_2d + 1
            endif
            if (meta%three_d) then
                out%vars_3d(j_3d) = meta%name
                j_3d = j_3d + 1
            endif
        enddo
    end subroutine classify_nest_init_vars

    ! Pack 2D init-only variables from domain storage into forcing_buffer_2d.
    ! The caller is responsible for the surrounding !$acc wait / update host /
    ! MPI_Win_Post sequence so that window logic stays visible at the call site.
    subroutine pack_init_vars_2d(this, domain, var_names)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(in)    :: domain
        character(len=*),  intent(in)    :: var_names(:)

        type(variable_t) :: var
        integer :: i, v_indx, v_idx, n_2d, nx, ny, ii, jj, i_s_w, i_e_w, j_s_w, j_e_w

        n_2d = 1
        do i = 1, kMAX_STORAGE_VARS
            if (var_names(i) == '') cycle
            v_indx = get_varindx(var_names(i))
            v_idx = domain%var_indx(v_indx)%v
            if (v_idx <= 0) then; n_2d = n_2d + 1; cycle; endif

            var = domain%vars_2d(v_idx)
            i_s_w = this%i_s_w; i_e_w = this%i_e_w + var%xstag
            j_s_w = this%j_s_w; j_e_w = this%j_e_w + var%ystag
            nx = i_e_w - i_s_w + 1; ny = j_e_w - j_s_w + 1

            if (var%dtype == kREAL) then
                associate(src => domain%vars_2d(v_idx)%data_2d, dst => this%forcing_buffer_2d)
                !$acc parallel loop gang vector collapse(2) present(dst, src)
                do jj = 1, ny
                  do ii = 1, nx
                    dst(n_2d, ii, jj) = src(i_s_w+ii-1, j_s_w+jj-1)
                  enddo
                enddo
                end associate
            elseif (var%dtype == kINTEGER) then
                associate(src => domain%vars_2d(v_idx)%data_2di, dst => this%forcing_buffer_2d)
                !$acc parallel loop gang vector collapse(2) present(dst, src)
                do jj = 1, ny
                  do ii = 1, nx
                    dst(n_2d, ii, jj) = real(src(i_s_w+ii-1, j_s_w+jj-1))
                  enddo
                enddo
                end associate
            endif
            n_2d = n_2d + 1
        enddo
    end subroutine pack_init_vars_2d

    ! Pack 3D init-only restart variables from domain storage into forcing_buffer_3d_init.
    ! Caller handles surrounding !$acc wait / update host / MPI_Win_Post.
    subroutine pack_init_vars_3d(this, domain, var_names)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(in)    :: domain
        character(len=*),  intent(in)    :: var_names(:)

        type(variable_t) :: var
        integer :: i, v_indx, v_idx, n_3d_i, nx, ny, nz_v, ii, jj, kk, i_s_w, i_e_w, j_s_w, j_e_w

        n_3d_i = 1
        do i = 1, kMAX_STORAGE_VARS
            if (var_names(i) == '') cycle
            v_indx = get_varindx(var_names(i))
            v_idx = domain%var_indx(v_indx)%v
            if (v_idx <= 0) then; n_3d_i = n_3d_i + 1; cycle; endif

            var = domain%vars_3d(v_idx)
            i_s_w = this%i_s_w; i_e_w = this%i_e_w + var%xstag
            j_s_w = this%j_s_w; j_e_w = this%j_e_w + var%ystag
            nx = i_e_w - i_s_w + 1; ny = j_e_w - j_s_w + 1
            nz_v = var%dim_len(2)

            associate(src => domain%vars_3d(v_idx)%data_3d, dst => this%forcing_buffer_3d_init)
            !$acc parallel loop gang vector collapse(3) present(dst, src)
            do jj = 1, ny
              do kk = 1, nz_v
                do ii = 1, nx
                  dst(n_3d_i, ii, kk, jj) = src(i_s_w+ii-1, kk, j_s_w+jj-1)
                enddo
              enddo
            enddo
            end associate
            n_3d_i = n_3d_i + 1
        enddo
    end subroutine pack_init_vars_3d

    ! Receive 2D init-only variables from parent nest via MPI_Recv, then geo-interpolate
    ! each from forcing grid to domain grid and write into domain storage.
    subroutine unpack_init_vars_2d(this, domain, forcing, var_names, n_init_2d, parent_rank)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(inout) :: domain
        type(boundary_t),  intent(in)    :: forcing
        character(len=*),  intent(in)    :: var_names(:)
        integer,           intent(in)    :: n_init_2d, parent_rank

        integer :: n, v, var_indx, idx, nx_r, ny_r, ierr
        real, allocatable :: recv_2d(:,:,:)
        real, allocatable :: forcing_field(:,:), domain_field(:,:)
        type(MPI_Status)  :: status
        type(meta_data_t) :: meta

        nx_r = this%i_e_r - this%i_s_r + 2
        ny_r = this%j_e_r - this%j_s_r + 2

        allocate(recv_2d(n_init_2d, nx_r, ny_r))
        call MPI_Recv(recv_2d, size(recv_2d), MPI_REAL, parent_rank, 98, &
                      this%parent_comms, status, ierr)

        ! 1-based bounds: the geolut stores 1-based local forcing-tile indices,
        ! not global parent-grid indices.
        allocate(forcing_field(this%i_e_r - this%i_s_r + 1, this%j_e_r - this%j_s_r + 1))
        allocate(domain_field(domain%ims:domain%ime, domain%jms:domain%jme))

        n = 1
        do v = 1, kMAX_STORAGE_VARS
            if (var_names(v) == '') cycle
            var_indx = get_varindx(var_names(v))
            idx = domain%var_indx(var_indx)%v
            if (idx <= 0) then; n = n + 1; cycle; endif

            ! Per-nest static var (geometry / grid metric / static surface
            ! descriptor) — the child reads/derives these from its own input
            ! file. Do NOT overwrite. Slot counter still advances so the next
            ! var lands in the correct recv_2d position.
            meta = get_varmeta(var_indx)
            if (meta%static_data) then
                n = n + 1
                cycle
            endif

            ! Unpack onto forcing grid (without stagger extra)
            forcing_field(:, :) = &
                recv_2d(n, 1:(this%i_e_r-this%i_s_r+1), 1:(this%j_e_r-this%j_s_r+1))

            call geo_interp2d(domain_field, forcing_field, forcing%geo%geolut)

            ! Write to domain
            if (domain%vars_2d(idx)%dtype == kREAL) then
                domain%vars_2d(idx)%data_2d(domain%ims:domain%ime, domain%jms:domain%jme) = domain_field
                !$acc update device(domain%vars_2d(idx)%data_2d)
            elseif (domain%vars_2d(idx)%dtype == kINTEGER) then
                domain%vars_2d(idx)%data_2di(domain%ims:domain%ime, domain%jms:domain%jme) = nint(domain_field)
                !$acc update device(domain%vars_2d(idx)%data_2di)
            endif
            n = n + 1
        enddo
    end subroutine unpack_init_vars_2d

    ! Receive 3D init-only variables from parent nest via MPI_Recv, then layer-by-layer
    ! geo-interpolate each from forcing grid to domain grid and write into domain storage.
    subroutine unpack_init_vars_3d(this, domain, forcing, var_names, n_init_3d, parent_rank)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(inout) :: domain
        type(boundary_t),  intent(in)    :: forcing
        character(len=*),  intent(in)    :: var_names(:)
        integer,           intent(in)    :: n_init_3d, parent_rank

        integer :: n, v, var_indx, idx, nx_r, ny_r, nz_r, nz_v, k, ierr
        integer :: msg_count
        real, allocatable :: recv_3d(:,:,:,:)
        real, allocatable :: forcing_field(:,:), domain_field(:,:)
        type(MPI_Status)  :: status, probe_status
        type(meta_data_t) :: meta

        nx_r = this%i_e_r - this%i_s_r + 2
        ny_r = this%j_e_r - this%j_s_r + 2

        ! Learn the k extent of the incoming init-3D message from its actual
        ! element count.
        call MPI_Probe(parent_rank, 99, this%parent_comms, probe_status, ierr)
        call MPI_Get_count(probe_status, MPI_REAL, msg_count, ierr)
        nz_r = msg_count / (n_init_3d * nx_r * ny_r)

        allocate(recv_3d(n_init_3d, nx_r, nz_r, ny_r))
        call MPI_Recv(recv_3d, size(recv_3d), MPI_REAL, parent_rank, 99, &
                      this%parent_comms, status, ierr)

        ! 1-based bounds: the geolut stores 1-based local forcing-tile indices,
        ! not global parent-grid indices.
        allocate(forcing_field(this%i_e_r - this%i_s_r + 1, this%j_e_r - this%j_s_r + 1))
        allocate(domain_field(domain%ims:domain%ime, domain%jms:domain%jme))

        n = 1
        do v = 1, kMAX_STORAGE_VARS
            if (var_names(v) == '') cycle
            var_indx = get_varindx(var_names(v))
            idx = domain%var_indx(var_indx)%v
            if (idx <= 0) then; n = n + 1; cycle; endif

            ! Per-nest static var (geometry / grid metric / static surface
            ! descriptor) — the child reads/derives these from its own input
            ! file. Do NOT overwrite. Slot counter still advances so the next
            ! var lands in the correct recv_3d position.
            meta = get_varmeta(var_indx)
            if (meta%static_data) then
                n = n + 1
                cycle
            endif

            nz_v = domain%vars_3d(idx)%dim_len(2)

            do k = 1, nz_v
                forcing_field(:, :) = &
                    recv_3d(n, 1:(this%i_e_r-this%i_s_r+1), k, 1:(this%j_e_r-this%j_s_r+1))

                call geo_interp2d(domain_field, forcing_field, forcing%geo%geolut)

                domain%vars_3d(idx)%data_3d(domain%ims:domain%ime, k, domain%jms:domain%jme) = domain_field
            enddo
            !$acc update device(domain%vars_3d(idx)%data_3d)
            n = n + 1
        enddo
    end subroutine unpack_init_vars_3d

    !Necesary? Should communicate relevant post/wait calls and free window
    subroutine close_client(this)
        implicit none
        class(ioclient_t), intent(inout) :: this

    end subroutine close_client

    subroutine smart_wait(window, err_msg)
        implicit none
        type(MPI_Win), intent(in) :: window
        character(len=*), intent(in) :: err_msg

        integer :: ierr
        logical :: flag
        double precision :: t_start, t_elapsed

        flag = .False.
        t_start = MPI_Wtime()

        call MPI_Win_Test(window, flag, ierr)
        do while (.not.(flag))

            t_elapsed = MPI_Wtime() - t_start
            if (t_elapsed > 600.0d0) then
                if (STD_OUT_PE) write(*,*) err_msg
                if (STD_OUT_PE) flush(output_unit)
                stop
            endif

            call MPI_Win_Test(window, flag, ierr)
            if (ierr /= 0) write(*,*) "MPI_Win_Test returned error: ",ierr
        end do
    end subroutine

end submodule
