
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

        integer :: my_rank, comm_size, some_child_id, ierr, my_nx, my_ny
        integer :: c_idx, child_id, p_idx
        logical :: child_needs_recv
        integer :: parent_idx, first_sibling_idx

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
        ! pick the first child nest id. Pass options(n_indx) — this nest IS the parent — as
        ! parent_opts so vars the parent doesn't allocate are filtered out of send_init_vars.
        if (size(options(n_indx)%general%child_nests) > 0) then
            some_child_id = options(n_indx)%general%child_nests(1)
            this%vars_for_nest = options(some_child_id)%forcing%vars_to_read
            call classify_nest_init_vars(options(some_child_id), options(n_indx), &
                                         this%vars_for_nest, this%send_init_vars)
        endif

        ! If this is a child nest, compute the 2D + 3D restart vars to receive
        ! from parent. The list MUST match the parent's send_init_vars exactly
        ! (same options used for classification, same iteration order in
        ! classify_nest_init_vars), or pack and unpack will use different
        ! slot indices and we get aliasing / kEMPT_BUFF leaking into vars.
        ! The parent's send-side (above, line ~67, plus init_ioserver in
        ! ioserver_obj.F90) classifies using its first child's options, so
        ! every child of a given parent must mirror that here — including
        ! children that aren't themselves the first sibling.
        if (options(n_indx)%general%parent_nest > 0) then
            parent_idx        = options(n_indx)%general%parent_nest
            first_sibling_idx = options(parent_idx)%general%child_nests(1)
            call classify_nest_init_vars(options(first_sibling_idx), &
                                            options(parent_idx), &
                                            options(first_sibling_idx)%forcing%vars_to_read, &
                                            this%recv_init_vars)
        endif

        ! Three-condition gate for the parent->child init transfer.
        ! recv side: this nest needs the transfer iff
        !   (.not. restart) .and. (parent_nest > 0)
        !                   .and. (start_time > parent's start_time).
        ! send side: this nest needs to send iff any of its children needs to recv.
        ! See ioclient_h.F90 for the flag semantics.

        this%nest_init_recv_done = .true.
        if (.not. options(n_indx)%restart%restart) then
            p_idx = options(n_indx)%general%parent_nest
            if (p_idx > 0) then
                if (options(n_indx)%general%start_time > options(p_idx)%general%start_time) then
                    this%nest_init_recv_done = .false.
                endif
            endif
        endif

        this%nest_init_send_done = .true.
        do c_idx = 1, size(options(n_indx)%general%child_nests)
            child_id = options(n_indx)%general%child_nests(c_idx)
            child_needs_recv = .false.
            if (.not. options(child_id)%restart%restart) then
                if (options(child_id)%general%parent_nest > 0) then
                    if (options(child_id)%general%start_time > &
                        options(options(child_id)%general%parent_nest)%general%start_time) then
                        child_needs_recv = .true.
                    endif
                endif
            endif
            if (child_needs_recv) then
                this%nest_init_send_done = .false.
                exit
            endif
        enddo

        ! Change 2: Store output/restart classification and compute counts
        this%vars_for_output = options(n_indx)%output%vars_for_output
        this%vars_for_restart = options(n_indx)%vars_for_restart
        this%restart_count = options(n_indx)%restart%restart_count
        this%restart_counter = 1
        if (options(n_indx)%restart%restart) this%restart_counter = 2

        call init_with_server(this)

        call setup_MPI_windows(this)

        ! Record the server's rank in parent_comms — used as destination /
        ! source for all two-sided MPI calls with the IO server.
        call MPI_Comm_Size(this%parent_comms, comm_size)
        this%server = comm_size - 1

        ! Pre-post the first main-forcing Irecv so the server's Isend from
        ! scatter_forcing (during our own wake OR during parent->child
        ! handoff that fires before we wake) does not stall on an
        ! unposted receive. Mirrors the old PSCW MPI_Win_Post(read_win)
        ! that lived at the end of init_ioclient in the RMA design.
        call MPI_Irecv(this%read_buffer, size(this%read_buffer), MPI_REAL, &
                       this%server, kIO_TAG_READ, this%parent_comms, &
                       this%read_req, ierr)

        ! Restart: pre-post Irecvs so the server's Isend in
        ! read_restart_file (which fires during the ioserver's wake, while
        ! the compute ranks are still in init_model_state) does not stall.
        if (options(n_indx)%restart%restart) then
            my_nx = this%i_e_re - this%i_s_re + 2
            my_ny = this%j_e_re - this%j_s_re + 2
            allocate(this%rst_scratch_3d(size(this%write_buffer_3d, 1), my_nx, &
                                          size(this%write_buffer_3d, 3), my_ny))
            allocate(this%rst_scratch_2d(size(this%write_buffer_2d, 1), my_nx, my_ny))
            call MPI_Irecv(this%rst_scratch_3d, size(this%rst_scratch_3d), MPI_REAL, &
                           this%server, kIO_TAG_RST_3D, this%parent_comms, &
                           this%rst_req_3d, ierr)
            call MPI_Irecv(this%rst_scratch_2d, size(this%rst_scratch_2d), MPI_REAL, &
                           this%server, kIO_TAG_RST_2D, this%parent_comms, &
                           this%rst_req_2d, ierr)
            this%rst_posted = .true.
        endif

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

        integer :: nx_w, nz_w, ny_w, nx_re, ny_re, n_w_2d, n_f, n_f_2d, n_f_3d_init, n_w_3d, ierr
        integer :: nx_r, nz_r, ny_r, n_r, n_out_3d, n_out_2d

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

        n_out_3d = 0
        n_out_2d = 0

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

        ! Output-only var counts (populated by the server end of the
        ! Allreduce; clients pass 0). Needed by the vector-type build below.
        call MPI_Allreduce(MPI_IN_PLACE,n_out_3d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,n_out_2d,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        ! Plain allocations — sent to server via Isend in push()
        allocate(this%write_buffer_3d(n_w_3d, nx_re, nz_w, ny_re))
        this%write_buffer_3d = kEMPT_BUFF
        !$acc enter data copyin(this%write_buffer_3d)

        allocate(this%write_buffer_2d(n_w_2d, nx_re, ny_re))
        this%write_buffer_2d = kEMPT_BUFF
        !$acc enter data copyin(this%write_buffer_2d)

        if (n_f > 0) then
            allocate(this%forcing_buffer(n_f, nx_w, nz_w, ny_w))
            this%forcing_buffer = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer)
        endif

        if (n_f_2d > 0) then
            allocate(this%forcing_buffer_2d(n_f_2d, nx_w, ny_w))
            this%forcing_buffer_2d = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer_2d)
        endif

        if (n_f_3d_init > 0) then
            allocate(this%forcing_buffer_3d_init(n_f_3d_init, nx_w, nz_w, ny_w))
            this%forcing_buffer_3d_init = kEMPT_BUFF
            !$acc enter data copyin(this%forcing_buffer_3d_init)
        endif

        ! read_buffer: plain heap allocation, filled by MPI_Irecv in receive()
        allocate(this%read_buffer(n_r, nx_r, nz_r, ny_r))
        this%read_buffer = kEMPT_BUFF
        !$acc enter data copyin(this%read_buffer)

        ! Build MPI vector types for the output-only Isend. On non-restart
        ! steps push() ships only the populated leading-dim slice
        ! write_buffer_*(1:n_out_*, :, ...) rather than the union-sized
        ! buffer; the vector type packs that slice from each (i,k,j) tile
        ! with stride n_w_* between blocks. If there is no restart-only
        ! tail (n_out == n_w) we keep using MPI_REAL — a flat send is
        ! cheaper than going through the type machinery.
        if (n_out_3d > 0 .and. n_out_3d < n_w_3d) then
            call MPI_Type_vector(nx_re*nz_w*ny_re, n_out_3d, n_w_3d, &
                                 MPI_REAL, this%send_type_3d_out, ierr)
            call MPI_Type_commit(this%send_type_3d_out, ierr)
        else
            this%send_type_3d_out = MPI_REAL
        endif
        if (n_out_2d > 0 .and. n_out_2d < n_w_2d) then
            call MPI_Type_vector(nx_re*ny_re, n_out_2d, n_w_2d, &
                                 MPI_REAL, this%send_type_2d_out, ierr)
            call MPI_Type_commit(this%send_type_2d_out, ierr)
        else
            this%send_type_2d_out = MPI_REAL
        endif

    end subroutine setup_MPI_windows


    ! This subroutine pushes the output fields from the domain object
    ! to the write buffer for the IO processes to use
    module subroutine push(this, domain)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),   intent(inout)    :: domain
        
        type(variable_t) :: var
        type(meta_data_t) :: tmp_var
        type(MPI_Request) :: reqs(2)
        integer :: i, n_3d, n_2d, nx, ny, nz_v, i_s_w, i_e_w, j_s_w, j_e_w, idx
        logical :: should_do_restart
        integer :: ii, jj, kk
        integer :: comm_size, my_rank, ierr

        n_3d = 1
        n_2d = 1
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
                call MPI_Send(domain%dt, 1, MPI_REAL, comm_size-1, kIO_TAG_DT_RESTART, this%parent_comms, ierr)
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

        ! Ship the packed buffers to the IO server.  Server posts matching
        ! MPI_Irecv in ioserver_obj::write_file.  Blocking Waitall below
        ! keeps the Isend buffers live until the server has received; they
        ! can then be safely overwritten by the next push().

        if (should_do_restart) then
            call MPI_Isend(this%write_buffer_3d, size(this%write_buffer_3d), MPI_REAL, &
                            this%server, kIO_TAG_WRITE_3D, this%parent_comms, reqs(1), ierr)
            call MPI_Isend(this%write_buffer_2d, size(this%write_buffer_2d), MPI_REAL, &
                            this%server, kIO_TAG_WRITE_2D, this%parent_comms, reqs(2), ierr)
        else
            call MPI_Isend(this%write_buffer_3d, 1, this%send_type_3d_out, &
                            this%server, kIO_TAG_WRITE_3D, this%parent_comms, reqs(1), ierr)
            call MPI_Isend(this%write_buffer_2d, 1, this%send_type_2d_out, &
                            this%server, kIO_TAG_WRITE_2D, this%parent_comms, reqs(2), ierr)
        endif
        call MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE, ierr)

    end subroutine 

    module subroutine update_nest(this, domain, options)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),      intent(in)    :: domain
        type(options_t),     intent(in)    :: options

        type(meta_data_t) :: tmp_var
        type(MPI_Request) :: req
        integer :: i, n_3d, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w, var_indx, idx
        logical :: var_val_check
        character(len=kMAX_NAME_LENGTH) :: err_msg
        integer :: ii, jj, kk, nz_v, ierr

        n_3d = 1
        this%nest_updated = .False.
        !This is false only when it is the first call to push (i.e. first write call)

        do i = 1, kMAX_STORAGE_VARS
            ! get the next variable in the structure
            if (this%vars_for_nest(i) == '') cycle
            var_indx = get_varindx(this%vars_for_nest(i))
            tmp_var = get_varmeta(var_indx)
            var_val_check = (tmp_var%maxval /= kUNSET_REAL .and. tmp_var%minval /= kUNSET_REAL)

            if (tmp_var%two_d) cycle

            if (domain%var_indx(var_indx)%v <= 0) then
                write(*,*) 'Error: Variable ', this%vars_for_nest(i), ' not found in parent domain: ', domain%nest_indx
                stop
            end if

            idx = domain%var_indx(var_indx)%v

            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            i_e_w = i_e_w+tmp_var%xstag !Add extra to accomodate staggered vars
            j_e_w = j_e_w+tmp_var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            associate(src => domain%vars_3d(idx)%data_3d, dst => this%forcing_buffer)
            nz_v = ubound(src, 2)
            
            !$acc parallel loop gang vector collapse(3) present(dst, src)
            do jj = 1, ny
              do kk = 1, nz_v
                do ii = 1, nx
                  dst(n_3d, ii, kk, jj) = src(i_s_w+ii-1, kk, j_s_w+jj-1)
                enddo
              enddo
            enddo
            end associate

            if (options%general%debug .and. var_val_check) then
                err_msg = 'Warning on ioclient_obj::update_nest: Nest level: '// str(domain%nest_indx)
                ! Check only the tile/interior region that is actually packed and sent
                ! to the child nest. Halo cells hold stale intermediate-stage values
                ! (uncorrected high-order RK3 leftovers) that are never used here.
                call check_var(domain%vars_3d(idx), trim(err_msg), is=i_s_w, ie=i_e_w, js=j_s_w, je=j_e_w)
            endif

            n_3d = n_3d+1
        enddo

        !$acc wait
        !$acc update host(this%forcing_buffer)

        this%nest_updated = .True.

        ! Ship packed forcing to the parent IO server.  Matching Irecv in
        ! ioserver_obj::gather_forcing.  Blocking Wait keeps forcing_buffer
        ! live until receive completes before we overwrite it next step.
        call MPI_Isend(this%forcing_buffer, size(this%forcing_buffer), MPI_REAL, &
                        this%server, kIO_TAG_NEST_3D, this%parent_comms, req, ierr)
        call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)

        ! One-time init: pack 2D + extra 3D restart vars for child nest initialization.
        ! Skipped when nest_init_send_done is pre-set in init_ioclient (no child of this
        ! nest needs the transfer per the three-condition rule). Sentinel'd at the end
        ! so subsequent update_nest calls during the run don't repeat the init send.
        if (.not. this%nest_init_send_done) then
            if (associated(this%forcing_buffer_2d)) then
                call pack_init_vars_2d(this, domain, this%send_init_vars%vars_2d)
                !$acc wait
                !$acc update host(this%forcing_buffer_2d)
                call MPI_Isend(this%forcing_buffer_2d, size(this%forcing_buffer_2d), MPI_REAL, &
                                this%server, kIO_TAG_NEST_2D, this%parent_comms, req, ierr)
                call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
            endif

            if (associated(this%forcing_buffer_3d_init)) then
                call pack_init_vars_3d(this, domain, this%send_init_vars%vars_3d)
                !$acc wait
                !$acc update host(this%forcing_buffer_3d_init)
                call MPI_Isend(this%forcing_buffer_3d_init, size(this%forcing_buffer_3d_init), MPI_REAL, &
                                this%server, kIO_TAG_NEST_3D_INIT, this%parent_comms, req, ierr)
                call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
            endif

            this%nest_init_send_done = .true.
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
        integer :: dbg_rank, dbg_ierr, ierr

        n = 1
        ! Wait on the Irecv posted last time (in init_ioclient or the
        ! tail of the previous receive() call). This lets the server's
        ! scatter_forcing Isend complete without blocking on our arrival
        ! in this routine — critical for init-time parent->child handoff.
        call MPI_Wait(this%read_req, MPI_STATUS_IGNORE, ierr)
        !$acc update device(this%read_buffer)
        ! Post the next Irecv now so a future scatter_forcing (next
        ! timestep or next parent->child handoff) has a matching
        ! receive already waiting.
        call MPI_Irecv(this%read_buffer, size(this%read_buffer), MPI_REAL, &
                        this%server, kIO_TAG_READ, this%parent_comms, &
                        this%read_req, ierr)

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
        integer :: my_nx, my_ny, ierr
        character(len=kMAX_NAME_LENGTH) :: varname

        ! Wait on the restart Irecvs pre-posted in init_ioclient (matching
        ! Isend is in ioserver_obj::read_restart_file, tags kIO_TAG_RST_3D
        ! / _2D). Pre-posting prevents the server's Isend from stalling
        ! while the compute rank is still in init_model_state.
        if (this%rst_posted) then
            call MPI_Wait(this%rst_req_3d, MPI_STATUS_IGNORE, ierr)
            call MPI_Wait(this%rst_req_2d, MPI_STATUS_IGNORE, ierr)
            my_nx = this%i_e_re - this%i_s_re + 2
            my_ny = this%j_e_re - this%j_s_re + 2
            this%write_buffer_3d(:, 1:my_nx, :, 1:my_ny) = this%rst_scratch_3d
            this%write_buffer_2d(:, 1:my_nx,    1:my_ny) = this%rst_scratch_2d
            !$acc update device(this%write_buffer_3d, this%write_buffer_2d)
            deallocate(this%rst_scratch_3d)
            deallocate(this%rst_scratch_2d)
            this%rst_posted = .false.
        endif

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

        integer :: n_init_2d, n_init_3d, parent_rank, nx_r, ny_r, ierr
        type(MPI_Status) :: status
        real, allocatable :: recv_2d(:,:,:)

        ! Skip the parent->child init transfer entirely when the flag is
        ! pre-set: the three conditions checked in init_ioclient determined
        ! that this nest does not need a state transfer from its parent
        ! (restart, no parent, or start_time not later than parent's).
        if (this%nest_init_recv_done) return

        n_init_2d = count(this%recv_init_vars%vars_2d /= '')
        n_init_3d = count(this%recv_init_vars%vars_3d /= '')

        if (n_init_2d == 0 .and. n_init_3d == 0) return

        if (STD_OUT_PE) write(*,*) "Receiving 2D+3D state from parent nest"
        if (STD_OUT_PE) flush(output_unit)

        ! init_ioclient already cached this%server (= comm_size-1 in parent_comms).
        parent_rank = this%server

        if (n_init_2d > 0) then
            nx_r = this%i_e_r - this%i_s_r + 2
            ny_r = this%j_e_r - this%j_s_r + 2
            allocate(recv_2d(n_init_2d, nx_r, ny_r))
            call MPI_Recv(recv_2d, size(recv_2d), MPI_REAL, parent_rank, 98, &
                          this%parent_comms, status, ierr)
            call unpack_init_vars_2d(this, domain, forcing, this%recv_init_vars%vars_2d, n_init_2d, recv_2d)
            deallocate(recv_2d)
        endif

        if (n_init_3d > 0) then
            call unpack_init_vars_3d(this, domain, forcing, this%recv_init_vars%vars_3d, n_init_3d, parent_rank)
        endif

        this%nest_init_recv_done = .true.

    end subroutine


    ! Filter/classify kMAX_STORAGE_VARS into 2D and 3D name lists for the
    ! init-only nest transfer. Both the parent-send path (in init_ioclient and
    ! init_ioserver) and the child-receive path (in init_ioclient) call this
    ! with the SAME `opts` (the parent's first child's options — which fixes
    ! slot ordering across siblings) and the SAME `parent_opts` (the parent's
    ! own options — which excludes vars the parent doesn't allocate so we
    ! never end up sending kEMPT_BUFF for a var the parent simply doesn't
    ! carry). The two filters together ensure pack and unpack iterate
    ! identical lists in identical order on every rank.
    subroutine classify_nest_init_vars(opts, parent_opts, forcing_vars, out)
        type(options_t),             intent(in)    :: opts
        type(options_t),             intent(in)    :: parent_opts
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
            ! The parent must actually have this var allocated, otherwise its
            ! pack would skip the slot and leave kEMPT_BUFF in forcing_buffer
            ! for the child to read.
            if (parent_opts%vars_to_allocate(v) <= 0) cycle
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
    ! MPI_Isend sequence so that the transport plumbing stays visible at the call site.
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
    ! Caller handles surrounding !$acc wait / update host / MPI_Isend.
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

    ! Geo-interpolate the 2D init slabs in recv_2d (filled by the caller via
    ! MPI_Recv) directly into domain storage on the device. The geolut
    ! (forcing%geo%geolut%x/y/w) is already enter-data'd on the device by
    ! domain_obj at init, so the only host->device traffic is the bulk
    ! recv_2d transfer issued here.
    subroutine unpack_init_vars_2d(this, domain, forcing, var_names, n_init_2d, recv_2d)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(inout) :: domain
        type(boundary_t),  intent(in)    :: forcing
        character(len=*),  intent(in)    :: var_names(:)
        integer,           intent(in)    :: n_init_2d
        real,              intent(inout) :: recv_2d(:,:,:)

        integer :: n, v, var_indx, idx
        integer :: i, j, l, localx, localy
        integer :: ims, ime, jms, jme
        real :: localw, local_center, tmp_val
        type(meta_data_t) :: meta

        ! Single bulk host->device transfer. The MPI buffer must stay on the
        ! host (parent_comms carries mpi_memory_alloc_kinds=system, set in
        ! init.F90), so we copy here once and then run all per-var interps
        ! on the device.
        !$acc enter data copyin(recv_2d)

        associate( gx => forcing%geo%geolut%x, &
                   gy => forcing%geo%geolut%y, &
                   gw => forcing%geo%geolut%w )

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

            if (domain%vars_2d(idx)%dtype == kREAL) then
                associate(dst => domain%vars_2d(idx)%data_2d)
                ims = lbound(dst, 1); ime = ubound(dst, 1)
                jms = lbound(dst, 2); jme = ubound(dst, 2)
                !$acc parallel loop gang vector collapse(2) &
                !$acc& present(dst, recv_2d, gx, gy, gw) &
                !$acc& private(l, localx, localy, localw, local_center)
                do j = jms, jme
                    do i = ims, ime
                        dst(i, j) = 0
                        local_center = 0
                        do l = 1, 4
                            localx = gx(l, i, j)
                            localy = gy(l, i, j)
                            local_center = local_center + recv_2d(n, localx, localy)
                            if (l < 3) then
                                localw = gw(l, i, j)
                                dst(i, j) = dst(i, j) + recv_2d(n, localx, localy) * localw
                            endif
                        enddo
                        localw = gw(3, i, j)
                        dst(i, j) = dst(i, j) + local_center / 4 * localw
                    enddo
                enddo
                end associate
            elseif (domain%vars_2d(idx)%dtype == kINTEGER) then
                associate(dst => domain%vars_2d(idx)%data_2di)
                ims = lbound(dst, 1); ime = ubound(dst, 1)
                jms = lbound(dst, 2); jme = ubound(dst, 2)
                !$acc parallel loop gang vector collapse(2) &
                !$acc& present(dst, recv_2d, gx, gy, gw) &
                !$acc& private(l, localx, localy, localw, local_center, tmp_val)
                do j = jms, jme
                    do i = ims, ime
                        tmp_val = 0
                        local_center = 0
                        do l = 1, 4
                            localx = gx(l, i, j)
                            localy = gy(l, i, j)
                            local_center = local_center + recv_2d(n, localx, localy)
                            if (l < 3) then
                                localw = gw(l, i, j)
                                tmp_val = tmp_val + recv_2d(n, localx, localy) * localw
                            endif
                        enddo
                        localw = gw(3, i, j)
                        tmp_val = tmp_val + local_center / 4 * localw
                        dst(i, j) = nint(tmp_val)
                    enddo
                enddo
                end associate
            endif
            n = n + 1
        enddo

        end associate

        !$acc exit data delete(recv_2d)
    end subroutine unpack_init_vars_2d

    ! Receive 3D init-only variables from parent nest via MPI_Recv, then geo-interpolate
    ! each from forcing grid directly into domain storage on the device, k-layer kernel
    ! collapsed into one launch per variable. As with the 2D path, the geolut is already
    ! resident on the device so the only host->device traffic is the bulk recv_3d transfer.
    subroutine unpack_init_vars_3d(this, domain, forcing, var_names, n_init_3d, parent_rank)
        class(ioclient_t), intent(inout) :: this
        type(domain_t),    intent(inout) :: domain
        type(boundary_t),  intent(in)    :: forcing
        character(len=*),  intent(in)    :: var_names(:)
        integer,           intent(in)    :: n_init_3d, parent_rank

        integer :: n, v, var_indx, idx, nx_r, ny_r, nz_r, nz_v, ierr
        integer :: msg_count
        integer :: i, j, k, l, localx, localy
        integer :: ims, ime, jms, jme
        real :: localw, local_center
        real, allocatable :: recv_3d(:,:,:,:)
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

        !$acc enter data copyin(recv_3d)

        associate( gx => forcing%geo%geolut%x, &
                   gy => forcing%geo%geolut%y, &
                   gw => forcing%geo%geolut%w )

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

            ! Skip the var entirely when parent and child disagree on its
            ! z-extent. Pack on the parent only writes [1:nz_v_parent]
            ! into the slot, leaving [nz_v_parent+1:nz_w] at the buffer's
            ! kEMPT_BUFF init value. If the cell at this child's nz_v layer
            ! is still the sentinel, the parent didn't write that deep, so
            ! its nz_v < ours and the protocol's slot cannot give us valid
            ! data for this var — child's own init pipeline (read_land,
            ! init_physics, restart) is responsible for it instead. Slot
            ! counter still advances to keep iteration in lockstep with pack.
            if (recv_3d(n, 1, nz_v, 1) == kEMPT_BUFF) then
                n = n + 1
                cycle
            endif

            associate(dst => domain%vars_3d(idx)%data_3d)
            ims = lbound(dst, 1); ime = ubound(dst, 1)
            jms = lbound(dst, 3); jme = ubound(dst, 3)
            !$acc parallel loop gang vector collapse(3) &
            !$acc& present(dst, recv_3d, gx, gy, gw) &
            !$acc& private(l, localx, localy, localw, local_center)
            do j = jms, jme
                do k = 1, nz_v
                    do i = ims, ime
                        dst(i, k, j) = 0
                        local_center = 0
                        do l = 1, 4
                            localx = gx(l, i, j)
                            localy = gy(l, i, j)
                            local_center = local_center + recv_3d(n, localx, k, localy)
                            if (l < 3) then
                                localw = gw(l, i, j)
                                dst(i, k, j) = dst(i, k, j) + recv_3d(n, localx, k, localy) * localw
                            endif
                        enddo
                        localw = gw(3, i, j)
                        dst(i, k, j) = dst(i, k, j) + local_center / 4 * localw
                    enddo
                enddo
            enddo
            end associate

            n = n + 1
        enddo

        end associate

        !$acc exit data delete(recv_3d)
        deallocate(recv_3d)
    end subroutine unpack_init_vars_3d

    !Necesary? Should communicate relevant post/wait calls and free window
    subroutine close_client(this)
        implicit none
        class(ioclient_t), intent(inout) :: this

        integer :: ierr

        ! Cancel the pre-posted main-forcing Irecv left outstanding by the
        ! last receive() call. MPI_Finalize may otherwise complain about
        ! a live request on some implementations.
        call MPI_Cancel(this%read_req, ierr)
        call MPI_Wait(this%read_req, MPI_STATUS_IGNORE, ierr)

        ! If restart was set but receive_rst was never reached (e.g. the
        ! job errored out during init), free the pre-posted restart reqs.
        if (this%rst_posted) then
            call MPI_Cancel(this%rst_req_3d, ierr)
            call MPI_Wait(this%rst_req_3d, MPI_STATUS_IGNORE, ierr)
            call MPI_Cancel(this%rst_req_2d, ierr)
            call MPI_Wait(this%rst_req_2d, MPI_STATUS_IGNORE, ierr)
            if (allocated(this%rst_scratch_3d)) deallocate(this%rst_scratch_3d)
            if (allocated(this%rst_scratch_2d)) deallocate(this%rst_scratch_2d)
            this%rst_posted = .false.
        endif

    end subroutine close_client

end submodule
