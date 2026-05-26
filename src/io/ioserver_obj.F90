
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
submodule(ioserver_interface) ioserver_implementation
  use debug_module,             only : check_ncdf
  use iso_fortran_env
  use iso_c_binding
  use output_metadata,          only : get_varindx, get_varmeta
  use meta_data_interface,      only : meta_data_t
  use time_object,              only : Time_type
  use string,                   only : split_str, as_string
  use io_routines,              only : check_file_exists, io_write
  use string,                   only : str
  use time_io,                  only : find_timestep_in_file
  use mpi_utils_module,         only : get_mpi_global_rank
  implicit none

contains


    module subroutine init_ioserver(this, options, nest_indx)
        class(ioserver_t),  intent(inout)  :: this
        type(options_t), intent(in)     :: options(:)
        integer,            intent(in)     :: nest_indx

        integer ::  n, n_3d, n_2d, var_indx, out_i, rst_i, some_child_id, ierr
        integer :: out_ord_i, rst_ord_i
        logical :: is_output, is_rst_only
        type(meta_data_t) :: tmp_meta
        integer :: v
        integer :: c_idx, child_id, p_idx
        logical :: child_needs_recv

        ! Call the parent type's init procedure
        call this%init_flow_obj(options(nest_indx),nest_indx)

        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        if (STD_OUT_PE_IO) write(*,*) "IOServer: initializing with clients"
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        ! Perform a series of MPI Gather/reduction to communicate the grid dimensions of the children clients to the parent server
        call init_with_clients(this)
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        if (STD_OUT_PE_IO) write(*,*) "IOServer: After init with clients in ioserver init"
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------" 
        !Determine number of IOServers in the this%IO_Comms communicator, i.e. size of communicator
        call MPI_Comm_size(this%IO_Comms, this%n_servers, ierr)

        !Setup reading capability
        if (options(nest_indx)%general%parent_nest == 0) then
            call this%reader%init(this%i_s_r,this%i_e_r,this%k_s_r,this%k_e_r,this%j_s_r,this%j_e_r,options(nest_indx))
            !this%n_children = size(this%children)
            this%n_r = this%reader%n_vars
            this%files_to_read = .not.(this%reader%eof)

            allocate(this%forcing_buffer(this%n_r,this%i_s_r:this%i_e_r+1,this%k_s_r:this%k_e_r, this%j_s_r:this%j_e_r+1))
        else
            this%n_r = count(options(nest_indx)%forcing%vars_to_read /= "")
            this%files_to_read = .False.
        endif

        this%n_child_ioservers = size(options(nest_indx)%general%child_nests)
        this%n_f = 0
        this%n_f_2d = 0
        this%n_f_3d_init = 0

        ! Three-condition gate for the parent->child init transfer.
        ! recv side: this nest needs the transfer iff
        !   (.not. restart) .and. (parent_nest > 0)
        !                   .and. (start_time > parent's start_time).
        ! send side: this nest needs to send iff any of its children needs to recv.
        ! Both flags are queried by distribute_init_forcing (send + per-child).

        this%nest_init_recv_done = .true.
        if (.not. options(nest_indx)%restart%restart) then
            p_idx = options(nest_indx)%general%parent_nest
            if (p_idx > 0) then
                if (options(nest_indx)%general%start_time > options(p_idx)%general%start_time) then
                    this%nest_init_recv_done = .false.
                endif
            endif
        endif

        this%nest_init_send_done = .true.
        do c_idx = 1, this%n_child_ioservers
            child_id = options(nest_indx)%general%child_nests(c_idx)
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

        !determine if we need to increase our k index due to some very large soil field
        this%n_w_3d = 0

        ! Save the number of vertical levels in the domain, since this is what our child ioserver will expect as the z-dimension for forcing data
        this%k_e_f = options(nest_indx)%domain%nz

        if (size(options(nest_indx)%general%child_nests) > 0) then
            some_child_id = options(nest_indx)%general%child_nests(1)
            this%n_f = count(options(some_child_id)%forcing%vars_to_read /= "")
            if (this%n_f > 0) then
                allocate(this%gather_buffer(this%n_f,this%i_s_w:this%i_e_w+1,this%k_s_w:this%k_e_f,this%j_s_w:this%j_e_w+1))
            endif

            allocate(this%send_nest_types(this%n_child_ioservers,this%n_servers))
            allocate(this%buffer_nest_types(this%n_child_ioservers,this%n_servers))
            allocate(this%nest_types_initialized(this%n_child_ioservers))
            this%nest_types_initialized = .false.

            ! Shared per-child geometry cache for the three setup_nest_types_* variants
            allocate(this%nest_geometry(this%n_child_ioservers))

            ! Count init-only restart vars (2D + extra 3D not in atmospheric
            ! forcing). Same filter set as classify_nest_init_vars in
            ! ioclient_obj.F90: must be allocated AND wanted-on-restart by
            ! the first child, AND allocated by the parent (this nest), so
            ! pack never has to skip a slot and leak kEMPT_BUFF.
            do v = 1, kMAX_STORAGE_VARS
                if (options(some_child_id)%vars_to_allocate(v) <= 0) cycle
                if (options(some_child_id)%vars_for_restart(v) <= 0) cycle
                if (options(nest_indx)%vars_to_allocate(v) <= 0) cycle
                tmp_meta = get_varmeta(v)
                if (len_trim(tmp_meta%name) == 0) cycle
                if (any(options(some_child_id)%forcing%vars_to_read == tmp_meta%name)) cycle
                if (tmp_meta%two_d)   this%n_f_2d = this%n_f_2d + 1
                if (tmp_meta%three_d) this%n_f_3d_init = this%n_f_3d_init + 1
            enddo

            if (this%n_f_2d > 0) then
                allocate(this%gather_buffer_2d(this%n_f_2d, this%i_s_w:this%i_e_w+1, this%j_s_w:this%j_e_w+1))
                allocate(this%send_nest_types_2d(this%n_child_ioservers, this%n_servers))
                allocate(this%buffer_nest_types_2d(this%n_child_ioservers, this%n_servers))
                allocate(this%nest_types_2d_initialized(this%n_child_ioservers))
                this%nest_types_2d_initialized = .false.
            endif
            if (this%n_f_3d_init > 0) then
                ! gather_buffer_3d_init is allocated later, after outputter init, so that
                ! its k extent can cover snow/soil vars whose dim_len(2) exceeds atm nz.
                allocate(this%send_nest_types_3d_init(this%n_child_ioservers, this%n_servers))
                allocate(this%buffer_nest_types_3d_init(this%n_child_ioservers, this%n_servers))
                allocate(this%nest_types_3d_init_initialized(this%n_child_ioservers))
                this%nest_types_3d_init_initialized = .false.
            endif
        endif
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        if (STD_OUT_PE_IO) write(*,*) "IOServer: Setting up outputter..."
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------" 

        !Setup writing capability
        call this%outputer%init(options(nest_indx),this%i_s_w,this%i_e_w,this%k_s_w,this%k_e_w,this%j_s_w,this%j_e_w,this%ide,this%kde,this%jde)


        do n = 1,this%outputer%n_vars
            if (this%outputer%var_meta(n)%three_d) then
                this%n_w_3d = this%n_w_3d+1
                if(this%outputer%var_meta(n)%dim_len(2) > this%k_e_w) this%k_e_w = this%outputer%var_meta(n)%dim_len(2)
            endif
        enddo

        this%n_w_2d = this%outputer%n_vars - this%n_w_3d

        ! Init-3D transport buffers: size their k extent to the max dim_len(2) across
        ! output vars (= k_e_w), so snow/soil restart vars (nsnow/nsoil layers) fit
        ! instead of being truncated to the atmospheric k_e_f.
        if (this%n_f_3d_init > 0) then
            this%nz_init_3d = this%k_e_w
            allocate(this%gather_buffer_3d_init(this%n_f_3d_init, this%i_s_w:this%i_e_w+1, &
                     1:this%nz_init_3d, this%j_s_w:this%j_e_w+1))
        endif

        if (this%n_w_2d == 0) write(*,*) 'Warning: No 2D variables set for output, this should never happen'
        if (this%n_w_3d == 0) write(*,*) 'Warning: No 3D variables set for output, this should never happen'

        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        if (STD_OUT_PE_IO) write(*,*) "IOServer: Setting up MPI Windows..."
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------" 

        this%n_out_3d = 0; this%n_out_2d = 0
        this%n_rst_only_3d = 0; this%n_rst_only_2d = 0
        do n = 1, this%outputer%n_vars
            var_indx = get_varindx(this%outputer%var_meta(n)%name)
            is_output = (options(nest_indx)%output%vars_for_output(var_indx) > 0)
            is_rst_only = (.not. is_output) .and. (options(nest_indx)%vars_for_restart(var_indx) > 0)
            if (this%outputer%var_meta(n)%three_d) then
                if (is_output) this%n_out_3d = this%n_out_3d + 1
                if (is_rst_only) this%n_rst_only_3d = this%n_rst_only_3d + 1
            else if (this%outputer%var_meta(n)%two_d) then
                if (is_output) this%n_out_2d = this%n_out_2d + 1
                if (is_rst_only) this%n_rst_only_2d = this%n_rst_only_2d + 1
            endif
        enddo

        call setup_MPI_windows(this)



        !Setup arrays for information about accessing variables from write buffer
        ! Count actual output/restart vars from the outputer to ensure allocation matches fill
        out_i = 0
        rst_i = 0
        do n=1,this%outputer%n_vars
            var_indx = get_varindx(this%outputer%var_meta(n)%name)
            if (options(nest_indx)%output%vars_for_output(var_indx) > 0) out_i = out_i + 1
            if (options(nest_indx)%vars_for_restart(var_indx) > 0) rst_i = rst_i + 1
        enddo

        allocate(this%out_var_indices(out_i))
        allocate(this%rst_var_indices(rst_i))

        this%restart_counter = 1
        this%restart_count = options(nest_indx)%restart%restart_count
        
        out_i = 1
        rst_i = 1
        
        do n=1,this%outputer%n_vars
            var_indx = get_varindx(this%outputer%var_meta(n)%name)
            if (options(nest_indx)%output%vars_for_output(var_indx) > 0) then
                this%out_var_indices(out_i) = n
                out_i = out_i + 1
            endif
            if (options(nest_indx)%vars_for_restart(var_indx) > 0) then
                this%rst_var_indices(rst_i) = n
                rst_i = rst_i + 1
            endif
        enddo


        ! Build ordered index arrays mapping buffer position -> outputer var index
        allocate(this%out_ordered_indices(this%n_out_3d + this%n_out_2d))
        allocate(this%rst_only_ordered_indices(this%n_rst_only_3d + this%n_rst_only_2d))
        out_ord_i = 1
        rst_ord_i = 1
        do n = 1, this%outputer%n_vars
            var_indx = get_varindx(this%outputer%var_meta(n)%name)
            is_output = (options(nest_indx)%output%vars_for_output(var_indx) > 0)
            is_rst_only = (.not. is_output) .and. (options(nest_indx)%vars_for_restart(var_indx) > 0)
            if (is_output) then
                this%out_ordered_indices(out_ord_i) = n
                out_ord_i = out_ord_i + 1
            else if (is_rst_only) then
                this%rst_only_ordered_indices(rst_ord_i) = n
                rst_ord_i = rst_ord_i + 1
            endif
        enddo

        if (options(nest_indx)%restart%restart) then
            call this%outputer%init_restart(options(nest_indx), this%IO_Comms, this%out_var_indices)
            ! if this is a restart run, then there will be no initial call to save_out_file, meaning
            ! that the restart counter should be auto-incremented to 2 (1 + 1) here 
            this%restart_counter = this%restart_counter + 1
        endif
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------"
        if (STD_OUT_PE_IO) write(*,*) "IOServer: Initialization finished"
        if (STD_OUT_PE_IO) write(*,*) "---------------------------------------------------" 

    end subroutine init_ioserver

    ! Build the per-child nest geometry: for each peer ioserver n, the list of
    ! (i,j) cells that participate in sends (this's write buffer -> n) and the
    ! list that participate in the receive (n -> child's read buffer). The
    ! geometry is invariant in the k axis and in n_vars, so it is reused across
    ! the 3D-per-step, 2D-init, and 3D-init datatype variants.
    !
    ! This is the shared first half of the old setup_nest_types routines: the
    ! Allreduces, the test_write_buffer_2d construction, and the per-server
    ! mask broadcasts. Only the cell enumeration is split off from the
    ! displacement arithmetic.
    module subroutine build_nest_geometry(this, child_ioserver, geom)
        class(ioserver_t),     intent(inout) :: this
        type(ioserver_t),      intent(in)    :: child_ioserver
        type(nest_geometry_t), intent(inout) :: geom

        integer :: n, my_rank, i, j, ierr, i_start, i_end, j_start, j_end, mask_size, cnt, c
        integer, allocatable, dimension(:) :: parent_ims, parent_ime, parent_jms, parent_jme
        integer, allocatable, dimension(:) :: child_isr, child_ier, child_jsr, child_jer
        real,    allocatable, dimension(:,:) :: mask, test_write_buffer_2d

        if (geom%built) return

        call MPI_Comm_rank(this%IO_Comms, my_rank, ierr)

        allocate(parent_ims(this%n_servers), parent_ime(this%n_servers), parent_jms(this%n_servers), parent_jme(this%n_servers))
        allocate(child_isr(this%n_servers),  child_ier(this%n_servers),  child_jsr(this%n_servers),  child_jer(this%n_servers))

        parent_ims = 0; parent_ime = 0; parent_jms = 0; parent_jme = 0
        child_isr  = 0; child_ier  = 0; child_jsr  = 0; child_jer  = 0

        parent_ims(my_rank+1) = this%i_s_w
        parent_ime(my_rank+1) = this%i_e_w
        parent_jms(my_rank+1) = this%j_s_w
        parent_jme(my_rank+1) = this%j_e_w

        child_isr(my_rank+1) = child_ioserver%i_s_r
        child_ier(my_rank+1) = child_ioserver%i_e_r
        child_jsr(my_rank+1) = child_ioserver%j_s_r
        child_jer(my_rank+1) = child_ioserver%j_e_r

        call MPI_Allreduce(MPI_IN_PLACE,child_isr,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_ier,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_jsr,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_jer,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,parent_ims,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_ime,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_jms,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_jme,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)

        ! Designate all cells for which this ioserver is responsible for gathering data.
        ! Only at the domain edge do we also include the staggered cells, since no other
        ! ioserver will cover them.
        allocate(test_write_buffer_2d(this%i_s_w:this%i_e_w+1, this%j_s_w:this%j_e_w+1))
        test_write_buffer_2d = kEMPT_BUFF
        do i=1,this%n_children
            if (this%iewc(i)==this%ide .and. this%jewc(i)==this%jde) then
                test_write_buffer_2d(this%iswc(i):this%iewc(i)+1,this%jswc(i):this%jewc(i)+1) = 1.0
            else if (this%iewc(i)==this%ide) then
                test_write_buffer_2d(this%iswc(i):this%iewc(i)+1,this%jswc(i):this%jewc(i)) = 1.0
            else if (this%jewc(i)==this%jde) then
                test_write_buffer_2d(this%iswc(i):this%iewc(i),this%jswc(i):this%jewc(i)+1) = 1.0
            else
                test_write_buffer_2d(this%iswc(i):this%iewc(i),this%jswc(i):this%jewc(i)) = 1.0
            endif
        enddo

        allocate(geom%send_cells(this%n_servers))
        allocate(geom%buff_cells(this%n_servers))

        do n = 1, this%n_servers
            ! ---------- send side: this's write buffer -> peer server n ----------
            i_start = max(child_isr(n),this%i_s_w)
            i_end   = min(child_ier(n),this%i_e_w)+1
            j_start = max(child_jsr(n),this%j_s_w)
            j_end   = min(child_jer(n),this%j_e_w)+1

            if (i_end >= i_start .and. j_end >= j_start) then
                cnt = count(test_write_buffer_2d(i_start:i_end,j_start:j_end) /= kEMPT_BUFF)
            else
                cnt = 0
            endif
            geom%send_cells(n)%n = cnt
            allocate(geom%send_cells(n)%i(cnt))
            allocate(geom%send_cells(n)%j(cnt))
            if (cnt > 0) then
                c = 0
                ! Enumerate in j-outer, i-inner order so that downstream displacement
                ! builders can walk the cell list as a grouped-by-j loop and reproduce
                ! the legacy do j; do k; do i; do v traversal byte-for-byte.
                do j = j_start, j_end
                    do i = i_start, i_end
                        if (test_write_buffer_2d(i,j) /= kEMPT_BUFF) then
                            c = c + 1
                            geom%send_cells(n)%i(c) = i
                            geom%send_cells(n)%j(c) = j
                        endif
                    enddo
                enddo
            endif

            ! ---------- buffer side: peer server n -> child's read buffer ----------
            ! We need the mask of the sending server to know which cells it actually owns.
            ! Each peer broadcasts its own test_write_buffer_2d in turn.
            if (allocated(mask)) deallocate(mask)
            allocate(mask(parent_ims(n):parent_ime(n)+1,parent_jms(n):parent_jme(n)+1))
            mask_size = (parent_ime(n)-parent_ims(n)+2)*(parent_jme(n)-parent_jms(n)+2)
            mask = kEMPT_BUFF
            if (my_rank == n-1) then
                mask(parent_ims(n):parent_ime(n)+1,parent_jms(n):parent_jme(n)+1) = &
                    test_write_buffer_2d(parent_ims(n):parent_ime(n)+1,parent_jms(n):parent_jme(n)+1)
            endif
            call MPI_Bcast(mask, mask_size, MPI_REAL, n-1, this%IO_Comms, ierr)

            i_start = max(parent_ims(n),child_ioserver%i_s_r)
            i_end   = min(parent_ime(n),child_ioserver%i_e_r)+1
            j_start = max(parent_jms(n),child_ioserver%j_s_r)
            j_end   = min(parent_jme(n),child_ioserver%j_e_r)+1

            if (i_end >= i_start .and. j_end >= j_start) then
                cnt = count(mask(i_start:i_end,j_start:j_end) /= kEMPT_BUFF)
            else
                cnt = 0
            endif
            geom%buff_cells(n)%n = cnt
            allocate(geom%buff_cells(n)%i(cnt))
            allocate(geom%buff_cells(n)%j(cnt))
            if (cnt > 0) then
                c = 0
                do j = j_start, j_end
                    do i = i_start, i_end
                        if (mask(i,j) /= kEMPT_BUFF) then
                            c = c + 1
                            geom%buff_cells(n)%i(c) = i
                            geom%buff_cells(n)%j(c) = j
                        endif
                    enddo
                enddo
            endif
        enddo

        geom%built = .true.
    end subroutine build_nest_geometry

    ! Build a committed MPI_Type_Indexed from a cached (i,j) cell list.
    !
    ! Displacement formula (matches the layout of a (v, i, k, j) Fortran array):
    !     disp = (v-1) + ((i - i_origin) + (k - k_s)*ni + (j - j_origin)*ni*nk) * n_vars
    !
    ! The buffer layout has v as the fastest dimension, so for any fixed (i,k,j)
    ! the n_vars values are contiguous in memory, and stepping i->i+1 lands
    ! immediately after the previous (v=1..n_vars) run. We exploit this by
    ! emitting one block per maximal contiguous (v=1..n_vars, i=run) span at
    ! fixed (k, j_cur), with block_length = n_vars * run_length. This is
    ! many orders of magnitude fewer indexed entries than the prior
    ! per-(v,i,k,j) length-1 emission, which is the dominant cost inside
    ! MPI_Alltoallw under MPI implementations that walk the indexed list
    ! sequentially. Both sides of the alltoallw use this same builder over
    ! their corresponding cell lists (built in j-outer/i-inner order by
    ! build_nest_geometry), so sender and receiver agree on the byte order.
    !
    ! For the 2D variant, pass k_s = k_e = 1 (nk = 1). The k-term then collapses
    ! to zero and the formula reduces to the 2D (v, i, j) layout.
    subroutine build_indexed_type_from_cells(cells, n_vars, i_origin, j_origin, ni, k_s, k_e, mpi_type)
        type(nest_cell_list_t), intent(in)  :: cells
        integer,                intent(in)  :: n_vars, i_origin, j_origin, ni, k_s, k_e
        type(MPI_Datatype),     intent(out) :: mpi_type

        integer :: nk, counter, c_start, c_end, cc, k, j_cur, run_start, run_end
        integer, allocatable :: block_lengths(:), displacements(:)

        nk = k_e - k_s + 1

        if (cells%n * nk == 0) then
            mpi_type = MPI_REAL
            call MPI_Type_commit(mpi_type)
            return
        endif

        ! Worst case: every cell is its own i-run (no two consecutive i's),
        ! repeated for every k. Real counts are typically far smaller.
        allocate(block_lengths(cells%n * nk))
        allocate(displacements(cells%n * nk))

        counter = 0
        c_start = 1
        do while (c_start <= cells%n)
            ! Advance c_end to cover the full run of cells sharing this j value.
            j_cur = cells%j(c_start)
            c_end = c_start
            do while (c_end < cells%n)
                if (cells%j(c_end + 1) /= j_cur) exit
                c_end = c_end + 1
            enddo

            do k = k_s, k_e
                ! Within this (j_cur, k), fuse cells with consecutive i values
                ! into a single block. Run [run_start..run_end] becomes one
                ! length-(n_vars*run_length) block at the run-start displacement.
                cc = c_start
                do while (cc <= c_end)
                    run_start = cc
                    do while (cc < c_end)
                        if (cells%i(cc + 1) /= cells%i(cc) + 1) exit
                        cc = cc + 1
                    enddo
                    run_end = cc

                    counter = counter + 1
                    block_lengths(counter) = n_vars * (run_end - run_start + 1)
                    displacements(counter) = ((cells%i(run_start) - i_origin) + &
                                              (k     - k_s)      * ni        + &
                                              (j_cur - j_origin) * ni * nk) * n_vars
                    cc = cc + 1
                enddo
            enddo

            c_start = c_end + 1
        enddo

        call MPI_Type_Indexed(counter, block_lengths(1:counter), displacements(1:counter), MPI_REAL, mpi_type)
        call MPI_Type_commit(mpi_type)
    end subroutine build_indexed_type_from_cells

    ! Release the cached (i,j) cell lists for a given child after all three
    ! datatype variants have been committed. The MPI datatypes themselves live
    ! on in this%send_nest_types* / this%buffer_nest_types*.
    subroutine free_nest_geometry(geom)
        type(nest_geometry_t), intent(inout) :: geom
        integer :: n

        if (allocated(geom%send_cells)) then
            do n = 1, size(geom%send_cells)
                if (allocated(geom%send_cells(n)%i)) deallocate(geom%send_cells(n)%i)
                if (allocated(geom%send_cells(n)%j)) deallocate(geom%send_cells(n)%j)
            enddo
            deallocate(geom%send_cells)
        endif
        if (allocated(geom%buff_cells)) then
            do n = 1, size(geom%buff_cells)
                if (allocated(geom%buff_cells(n)%i)) deallocate(geom%buff_cells(n)%i)
                if (allocated(geom%buff_cells(n)%j)) deallocate(geom%buff_cells(n)%j)
            enddo
            deallocate(geom%buff_cells)
        endif
        geom%built = .false.
    end subroutine free_nest_geometry

    ! Set up the per-timestep 3D nest datatypes. Delegates cell enumeration to
    ! build_nest_geometry (cached) and displacement arithmetic to
    ! build_indexed_type_from_cells.
    module subroutine setup_nest_types(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
        class(ioserver_t), intent(inout) :: this
        type(ioserver_t),  intent(in)    :: child_ioserver
        integer,           intent(in)    :: child_indx
        type(MPI_Datatype),intent(out)   :: send_nest_types(:), buffer_nest_types(:)

        integer :: n, ni_send, ni_buff

        call this%build_nest_geometry(child_ioserver, this%nest_geometry(child_indx))

        ni_send = this%i_e_w - this%i_s_w + 2
        ni_buff = child_ioserver%i_e_r - child_ioserver%i_s_r + 2

        associate (geom => this%nest_geometry(child_indx))
            do n = 1, this%n_servers
                call build_indexed_type_from_cells(geom%send_cells(n), this%n_f,       &
                    this%i_s_w, this%j_s_w, ni_send, this%k_s_w, this%k_e_f,           &
                    send_nest_types(n))
                call build_indexed_type_from_cells(geom%buff_cells(n), this%n_f,       &
                    child_ioserver%i_s_r, child_ioserver%j_s_r, ni_buff,               &
                    child_ioserver%k_s_r, child_ioserver%k_e_r,                        &
                    buffer_nest_types(n))
            enddo
        end associate
    end subroutine setup_nest_types

    ! 2D variant: no k dimension. Buffer is (n_f_2d, i, j). We pass k_s = k_e = 1
    ! so the k-term in the displacement formula collapses to zero.
    module subroutine setup_nest_types_2d(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
        class(ioserver_t), intent(inout) :: this
        type(ioserver_t),  intent(in)    :: child_ioserver
        integer,           intent(in)    :: child_indx
        type(MPI_Datatype),intent(out)   :: send_nest_types(:), buffer_nest_types(:)

        integer :: n, ni_send, ni_buff

        call this%build_nest_geometry(child_ioserver, this%nest_geometry(child_indx))

        ni_send = this%i_e_w - this%i_s_w + 2
        ni_buff = child_ioserver%i_e_r - child_ioserver%i_s_r + 2

        associate (geom => this%nest_geometry(child_indx))
            do n = 1, this%n_servers
                call build_indexed_type_from_cells(geom%send_cells(n), this%n_f_2d,    &
                    this%i_s_w, this%j_s_w, ni_send, 1, 1,                             &
                    send_nest_types(n))
                call build_indexed_type_from_cells(geom%buff_cells(n), this%n_f_2d,    &
                    child_ioserver%i_s_r, child_ioserver%j_s_r, ni_buff, 1, 1,         &
                    buffer_nest_types(n))
            enddo
        end associate
    end subroutine setup_nest_types_2d

    ! 3D init variant: same (v, i, k, j) layout as setup_nest_types but with
    ! n_f_3d_init variables instead of n_f.
    module subroutine setup_nest_types_3d_init(this, child_ioserver, child_indx, send_nest_types, buffer_nest_types)
        class(ioserver_t), intent(inout) :: this
        type(ioserver_t),  intent(in)    :: child_ioserver
        integer,           intent(in)    :: child_indx
        type(MPI_Datatype),intent(out)   :: send_nest_types(:), buffer_nest_types(:)

        integer :: n, ni_send, ni_buff

        call this%build_nest_geometry(child_ioserver, this%nest_geometry(child_indx))

        ni_send = this%i_e_w - this%i_s_w + 2
        ni_buff = child_ioserver%i_e_r - child_ioserver%i_s_r + 2

        associate (geom => this%nest_geometry(child_indx))
            do n = 1, this%n_servers
                ! Both sides share the same k extent after distribute_init_forcing
                ! resolves the per-family nz_init_3d and writes it to both
                ! this%nz_init_3d and each child%nz_init_3d on this rank.
                ! Using this%nz_init_3d for send and child_ioserver%nz_init_3d
                ! for buff makes the "nk matches the physical buffer stride on
                ! its side" invariant obvious at each call site.
                call build_indexed_type_from_cells(geom%send_cells(n), this%n_f_3d_init, &
                    this%i_s_w, this%j_s_w, ni_send, 1, this%nz_init_3d,                 &
                    send_nest_types(n))
                call build_indexed_type_from_cells(geom%buff_cells(n), this%n_f_3d_init, &
                    child_ioserver%i_s_r, child_ioserver%j_s_r, ni_buff,                 &
                    1, child_ioserver%nz_init_3d,                                        &
                    buffer_nest_types(n))
            enddo
        end associate
    end subroutine setup_nest_types_3d_init

    subroutine test_nest_types(ioserver, child_ioserver, child_indx)
        implicit none
        class(ioserver_t), intent(in) :: ioserver
        class(ioserver_t), intent(in) :: child_ioserver
        integer, intent(in) :: child_indx

        integer :: i,j, nx, ny, n, ierr, msg_size, real_size
        integer :: ix, jy, k, v
        real :: forcing_buffer_init, gather_buffer_init
        integer, allocatable :: send_msg_size_alltoall(:), buff_msg_size_alltoall(:), disp_alltoall(:)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lowerbound, extent
        real, allocatable, dimension(:,:,:,:) :: forcing_buffer_test, gather_buffer_test
        logical :: tripped = .False.
        integer :: PE_rank_global

        allocate(send_msg_size_alltoall(ioserver%n_servers))
        allocate(buff_msg_size_alltoall(ioserver%n_servers))
        allocate(disp_alltoall(ioserver%n_servers))

        PE_rank_global = get_mpi_global_rank()

        disp_alltoall = 0
        forcing_buffer_init = -1000.0*PE_rank_global
        gather_buffer_init = 1e30
        allocate(forcing_buffer_test(ioserver%n_f,child_ioserver%i_s_r:child_ioserver%i_e_r+1,child_ioserver%k_s_r:child_ioserver%k_e_r, child_ioserver%j_s_r:child_ioserver%j_e_r+1))
        allocate(gather_buffer_test(ioserver%n_f, ioserver%i_s_w:ioserver%i_e_w+1,ioserver%k_s_w:ioserver%k_e_f,ioserver%j_s_w:ioserver%j_e_w+1))

        forcing_buffer_test = forcing_buffer_init
        gather_buffer_test = gather_buffer_init
        ! Loop through child images and get chunks of buffer array from each one
        do i=1,ioserver%n_children
            do jy = ioserver%jswc(i), ioserver%jewc(i)+1
                do k = ioserver%k_s_w, ioserver%k_e_f
                    do ix = ioserver%iswc(i), ioserver%iewc(i)+1
                        gather_buffer_test(1,ix,k,jy) = ix + (k-1)*(ioserver%ide+1) + (jy-1)*(ioserver%ide+1)*(ioserver%kde)
                    enddo
                enddo
            enddo
            do jy = ioserver%jswc(i), ioserver%jewc(i)
                do k = ioserver%k_s_w, ioserver%k_e_f
                    do ix = ioserver%iswc(i), ioserver%iewc(i)+1
                        gather_buffer_test(2,ix,k,jy) = ix + (k-1)*(ioserver%ide+1) + (jy-1)*(ioserver%ide+1)*(ioserver%kde)
                    enddo
                enddo
            enddo
            do jy = ioserver%jswc(i), ioserver%jewc(i)+1
                do k = ioserver%k_s_w, ioserver%k_e_f
                    do ix = ioserver%iswc(i), ioserver%iewc(i)
                        gather_buffer_test(3,ix,k,jy) = ix + (k-1)*(ioserver%ide+1) + (jy-1)*(ioserver%ide+1)*(ioserver%kde)
                    enddo
                enddo
            enddo
            do jy = ioserver%jswc(i), ioserver%jewc(i)
                do k = ioserver%k_s_w, ioserver%k_e_f
                    do ix = ioserver%iswc(i), ioserver%iewc(i)
                        gather_buffer_test(4,ix,k,jy) = ix + (k-1)*(ioserver%ide+1) + (jy-1)*(ioserver%ide+1)*(ioserver%kde)
                    enddo
                enddo
            enddo
        enddo

        ! Get size of an MPI_REAL
        call MPI_Type_size(MPI_REAL, real_size, ierr)
        do n = 1, ioserver%n_servers
            call MPI_Type_get_extent(ioserver%send_nest_types(child_indx,n), lowerbound, extent)
            if (extent > real_size) then
                send_msg_size_alltoall(n) = 1
            else
                send_msg_size_alltoall(n) = 0
            endif
            call MPI_Type_get_extent(ioserver%buffer_nest_types(child_indx,n), lowerbound, extent)
            if (extent > real_size) then
                buff_msg_size_alltoall(n) = 1
            else
                buff_msg_size_alltoall(n) = 0
            endif
        enddo

        call MPI_Alltoallw(gather_buffer_test,  send_msg_size_alltoall, disp_alltoall, ioserver%send_nest_types(child_indx,:), &
                          forcing_buffer_test,  buff_msg_size_alltoall, disp_alltoall, ioserver%buffer_nest_types(child_indx,:), ioserver%IO_Comms)


        ! check that the whole of the forcing buffer has been filled correctly
        do j = child_ioserver%j_s_r, child_ioserver%j_e_r+1
            do k = child_ioserver%k_s_r, child_ioserver%k_e_r
                do i = child_ioserver%i_s_r, child_ioserver%i_e_r+1
                    if (forcing_buffer_test(1,i,k,j) /= i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)) then
                        write(*,*) 'Error in test_nest_types: Mismatch in forcing buffer at (', 1, ',', i, ',', k, ',', j, ') : ', &
                                    ' expected ', (i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)), &
                                    ' but got ', forcing_buffer_test(1,i,k,j)
                        write(*,*) 'child_ioserver%i_s_r: ', child_ioserver%i_s_r, ' child_ioserver%i_e_r: ', child_ioserver%i_e_r
                        write(*,*) 'child_ioserver%j_s_r: ', child_ioserver%j_s_r, ' child_ioserver%j_e_r: ', child_ioserver%j_e_r
                        
                        if (forcing_buffer_test(1,i,k,j) < 0) then
                            forcing_buffer_test(1,i,k,j) = forcing_buffer_init
                            write(*,*) ' (Buffer value indicates that alltoall was done incorrectly for x and y stagger)'
                        endif
                        if (forcing_buffer_test(1,i,k,j) == gather_buffer_init) then
                            write(*,*) ' (Buffer value indicates that gather operation was done incorrectly)'
                        endif
                        write(*,*) '------------------------------------------------------------------------'
                        tripped = .True.
                    endif
                enddo
            enddo
        enddo
        do j = child_ioserver%j_s_r, child_ioserver%j_e_r
            do k = child_ioserver%k_s_r, child_ioserver%k_e_r
                do i = child_ioserver%i_s_r, child_ioserver%i_e_r+1
                    if (forcing_buffer_test(2,i,k,j) /= i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)) then
                        write(*,*) 'Error in test_nest_types: Mismatch in forcing buffer at (', 2, ',', i, ',', k, ',', j, ') : ', &
                                    ' expected ', (i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)), &
                                    ' but got ', forcing_buffer_test(2,i,k,j)
                        write(*,*) 'child_ioserver%i_s_r: ', child_ioserver%i_s_r, ' child_ioserver%i_e_r: ', child_ioserver%i_e_r
                        write(*,*) 'child_ioserver%j_s_r: ', child_ioserver%j_s_r, ' child_ioserver%j_e_r: ', child_ioserver%j_e_r
                        if (forcing_buffer_test(2,i,k,j) < 0) then
                            forcing_buffer_test(2,i,k,j) = forcing_buffer_init
                            write(*,*) ' (Buffer value indicates that alltoall was done incorrectly for x)'
                        endif
                        if (forcing_buffer_test(2,i,k,j) == gather_buffer_init) then
                            write(*,*) ' (Buffer value indicates that gather operation was done incorrectly)'
                        endif
                        write(*,*) '------------------------------------------------------------------------'
                        tripped = .True.
                    endif
                enddo
            enddo
        enddo
        do j = child_ioserver%j_s_r, child_ioserver%j_e_r+1
            do k = child_ioserver%k_s_r, child_ioserver%k_e_r
                do i = child_ioserver%i_s_r, child_ioserver%i_e_r
                    if (forcing_buffer_test(3,i,k,j) /= i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)) then
                        write(*,*) 'Error in test_nest_types: Mismatch in forcing buffer at (', 3, ',', i, ',', k, ',', j, ') : ', &
                                    ' expected ', (i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)), &
                                    ' but got ', forcing_buffer_test(3,i,k,j)
                        write(*,*) 'child_ioserver%i_s_r: ', child_ioserver%i_s_r, ' child_ioserver%i_e_r: ', child_ioserver%i_e_r
                        write(*,*) 'child_ioserver%j_s_r: ', child_ioserver%j_s_r, ' child_ioserver%j_e_r: ', child_ioserver%j_e_r
                        if (forcing_buffer_test(3,i,k,j) < 0) then
                            forcing_buffer_test(3,i,k,j) = forcing_buffer_init
                            write(*,*) ' (Buffer value indicates that alltoall was done incorrectly for y stagger)'
                        endif
                        if (forcing_buffer_test(3,i,k,j) == gather_buffer_init) then
                            write(*,*) ' (Buffer value indicates that gather operation was done incorrectly)'
                        endif
                        write(*,*) '------------------------------------------------------------------------'
                        tripped = .True.
                    endif
                enddo
            enddo
        enddo
        do j = child_ioserver%j_s_r, child_ioserver%j_e_r
            do k = child_ioserver%k_s_r, child_ioserver%k_e_r
                do i = child_ioserver%i_s_r, child_ioserver%i_e_r
                    if (forcing_buffer_test(4,i,k,j) /= i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)) then
                        write(*,*) 'Error in test_nest_types: Mismatch in forcing buffer at (', 4, ',', i, ',', k, ',', j, ') : ', &
                                    ' expected ', (i + (k-1)*(ioserver%ide+1) + (j-1)*(ioserver%ide+1)*(ioserver%kde)), &
                                    ' but got ', forcing_buffer_test(4,i,k,j)
                        write(*,*) 'child_ioserver%i_s_r: ', child_ioserver%i_s_r, ' child_ioserver%i_e_r: ', child_ioserver%i_e_r
                        write(*,*) 'child_ioserver%j_s_r: ', child_ioserver%j_s_r, ' child_ioserver%j_e_r: ', child_ioserver%j_e_r
                        if (forcing_buffer_test(4,i,k,j) < 0) then
                            forcing_buffer_test(4,i,k,j) = forcing_buffer_init
                            write(*,*) ' (Buffer value indicates that alltoall was done incorrectly)'
                        endif
                        if (forcing_buffer_test(4,i,k,j) == gather_buffer_init) then
                            write(*,*) ' (Buffer value indicates that gather operation was done incorrectly)'
                        endif
                        write(*,*) '------------------------------------------------------------------------'
                        tripped = .True.
                    endif
                enddo
            enddo
        enddo

        if (tripped) then
            !get rank of this MPI process on the global communicator
            call io_write('debug/buffer_'//trim(str(PE_rank_global))//".nc",'buffer',forcing_buffer_test)
        endif
    end subroutine test_nest_types
    
    subroutine setup_MPI_windows(this)
        class(ioserver_t),   intent(inout)  :: this

        integer :: ierr, n

        ! +1 added to handle variables on staggered grids
        this%nx_w = maxval(this%iewc-this%iswc+1)+1
        this%nz_w = this%k_e_w ! this will have been updated in init to reflect the largest z dimension to be output
        this%ny_w = maxval(this%jewc-this%jswc+1)+1

        this%nx_re = maxval(this%ierec-this%isrec+1)+1
        this%ny_re = maxval(this%jerec-this%jsrec+1)+1

        this%nx_r = maxval(this%ierc-this%isrc+1)+1
        this%nz_r = maxval(this%kerc-this%ksrc+1)
        this%ny_r = maxval(this%jerc-this%jsrc+1)+1

       ! Setup MPI windows for inter-process communication        
        call MPI_Allreduce(MPI_IN_PLACE,this%nx_w,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%ny_w,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%nz_w,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,this%nx_re,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%ny_re,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,this%nx_r,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%ny_r,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%nz_r,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,this%n_w_3d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%n_w_2d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%n_f,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%n_f_2d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%n_f_3d_init,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%nz_init_3d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,this%n_r,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        ! Propagate output-only counts to the matching MPI_Allreduce on the
        ! client side (ioclient_obj::setup_MPI_windows). Both ends use the
        ! same client_comms / parent_comms; the server holds the real value
        ! and clients pass 0, so MAX gives both sides a consistent count
        ! for building the output-only Isend/Irecv MPI vector types.
        call MPI_Allreduce(MPI_IN_PLACE,this%n_out_3d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%n_out_2d,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        allocate(this%write_buffer_2d(this%n_children))
        allocate(this%write_buffer_3d(this%n_children))
        allocate(this%read_buffer(this%n_children))
        allocate(this%child_gather_buffers(this%n_children))

        ! Plain per-child receive buffers — filled by Irecv in write_file
        do n = 1,this%n_children
            allocate(this%write_buffer_3d(n)%buff(this%n_w_3d, this%nx_re, this%nz_w, this%ny_re))
            this%write_buffer_3d(n)%buff = kEMPT_BUFF
            allocate(this%write_buffer_2d(n)%buff(this%n_w_2d, this%nx_re, this%ny_re))
            this%write_buffer_2d(n)%buff = kEMPT_BUFF
        enddo

        ! MPI vector types matching the output-only Isends from the client
        ! (ioclient_obj::push, output-only branch). On non-restart steps the
        ! Irecv at write_file uses these to consume just the leading-dim
        ! slice (1:n_out_*) of the write_buffer with stride n_w_*. If there
        ! is no restart-only tail (n_out == n_w), fall back to flat MPI_REAL.
        if (this%n_out_3d > 0 .and. this%n_out_3d < this%n_w_3d) then
            call MPI_Type_vector(this%nx_re*this%nz_w*this%ny_re, this%n_out_3d, this%n_w_3d, &
                                 MPI_REAL, this%recv_type_3d_out, ierr)
            call MPI_Type_commit(this%recv_type_3d_out, ierr)
        else
            this%recv_type_3d_out = MPI_REAL
        endif
        if (this%n_out_2d > 0 .and. this%n_out_2d < this%n_w_2d) then
            call MPI_Type_vector(this%nx_re*this%ny_re, this%n_out_2d, this%n_w_2d, &
                                 MPI_REAL, this%recv_type_2d_out, ierr)
            call MPI_Type_commit(this%recv_type_2d_out, ierr)
        else
            this%recv_type_2d_out = MPI_REAL
        endif

        ! Plain per-child receive buffers for nest forcing — filled by
        ! Irecv in gather_forcing / gather_forcing_2d / gather_forcing_3d_init.
        if (this%n_f > 0) then
            do n = 1,this%n_children
                allocate(this%child_gather_buffers(n)%buff(this%n_f, this%nx_w, this%nz_w, this%ny_w))
                this%child_gather_buffers(n)%buff = kEMPT_BUFF
            enddo
            ! Persistent Irecv requests for gather_forcing pre-posting.
            allocate(this%gather_reqs(this%n_children))
        endif

        if (this%n_f_2d > 0) then
            allocate(this%child_gather_buffers_2d(this%n_children))
            do n = 1,this%n_children
                allocate(this%child_gather_buffers_2d(n)%buff(this%n_f_2d, this%nx_w, this%ny_w))
                this%child_gather_buffers_2d(n)%buff = kEMPT_BUFF
            enddo
        endif

        if (this%n_f_3d_init > 0) then
            allocate(this%child_gather_buffers_3d_init(this%n_children))
            do n = 1,this%n_children
                allocate(this%child_gather_buffers_3d_init(n)%buff(this%n_f_3d_init, this%nx_w, this%nz_w, this%ny_w))
                this%child_gather_buffers_3d_init(n)%buff = kEMPT_BUFF
            enddo
        endif

        ! read_buffer: plain per-child send buffers used by scatter_forcing
        ! (one MPI_Isend per child, kIO_TAG_READ).  No shared-memory window.
        do n = 1,this%n_children
            allocate(this%read_buffer(n)%buff(this%n_r, this%nx_r, this%nz_r, this%ny_r))
            this%read_buffer(n)%buff = kEMPT_BUFF
        enddo

    end subroutine setup_MPI_windows

    subroutine init_with_clients(this)
        class(ioserver_t),   intent(inout)  :: this
        integer :: n, comm_size, ierr, i
        integer, allocatable, dimension(:) :: cnts, disps

        !get number of clients on this communicator
        call MPI_Comm_size(this%client_comms, comm_size)
        ! don't forget about oursevles
        this%n_children = comm_size - 1

        n = 0

        allocate(this%iswc(this%n_children))
        allocate(this%iewc(this%n_children))
        allocate(this%kswc(this%n_children))
        allocate(this%kewc(this%n_children))
        allocate(this%jswc(this%n_children))
        allocate(this%jewc(this%n_children))
        allocate(this%isrec(this%n_children))
        allocate(this%ierec(this%n_children))
        allocate(this%jsrec(this%n_children))
        allocate(this%jerec(this%n_children))
        allocate(this%isrc(this%n_children))
        allocate(this%ierc(this%n_children))
        allocate(this%ksrc(this%n_children))
        allocate(this%kerc(this%n_children))
        allocate(this%jsrc(this%n_children))
        allocate(this%jerc(this%n_children))

        allocate(cnts(this%n_children+1))
        allocate(disps(this%n_children+1))

        cnts = 1
        cnts(this%n_children+1) = 0

        disps = (/(i, i=0,this%n_children, 1)/)
        disps(this%n_children+1) = disps(this%n_children+1)-1

        ! The PE of the parent communicator is the last PE in the communicator, i.e. n_children

        call MPI_Gatherv(n, 0, MPI_INTEGER, this%iswc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%iewc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%kswc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%kewc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jswc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jewc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%isrec, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%ierec, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jsrec, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jerec, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%isrc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%ierc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%ksrc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%kerc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jsrc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)
        call MPI_Gatherv(n, 0, MPI_INTEGER, this%jerc, cnts, disps, MPI_INTEGER, this%n_children, this%client_comms)

        this%ide = 0
        this%kde = 0
        this%jde = 0

        call MPI_Allreduce(MPI_IN_PLACE,this%ide,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%kde,1,MPI_INT,MPI_MAX,this%client_comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,this%jde,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        allocate(this%children_ranks(this%n_children))

        do n = 1,this%n_children
            this%children_ranks(n) = n-1
        enddo

        this%i_s_r = minval(this%isrc)
        this%i_e_r = maxval(this%ierc)
        this%i_s_w = minval(this%iswc)
        this%i_e_w = maxval(this%iewc)
        this%i_s_re = minval(this%isrec)
        this%i_e_re = maxval(this%ierec)

        this%j_s_r = minval(this%jsrc)
        this%j_e_r = maxval(this%jerc)
        this%j_s_w = minval(this%jswc)
        this%j_e_w = maxval(this%jewc)
        this%j_s_re = minval(this%jsrec)
        this%j_e_re = maxval(this%jerec)

        this%k_s_r = minval(this%ksrc)
        this%k_e_r = maxval(this%kerc)
        this%k_s_w = minval(this%kswc)
        this%k_e_w = maxval(this%kewc)

    end subroutine init_with_clients
    
    ! This subroutine gathers the write buffers of its children 
    ! compute processes and then writes them to the output file
    module subroutine write_file(this)
        implicit none
        class(ioserver_t), intent(inout)  :: this

        integer :: i, msg_size, ierr, cc
        integer :: n, n_3d, n_2d, x_stag, y_stag, oi
        logical :: should_write_restart
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
        type(MPI_Request) :: reqs(2 * this%n_children)

        msg_size = 1
        disp = 0

        ! check to see if we should collect and write restart data
        should_write_restart = .False.
        should_write_restart = (this%restart_counter > this%restart_count)

        if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",this%nest_indx," --------------"
        if (STD_OUT_PE_IO) write(*,*) "Fetching data from child images for output"
        if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",this%nest_indx," --------------"

        ! Receive dt from compute rank 0 for restart file
        if (should_write_restart .or. this%first_write) then
            call MPI_Recv(this%dt, 1, MPI_REAL, 0, kIO_TAG_DT_RESTART, this%client_comms, MPI_STATUS_IGNORE, ierr)
        endif

        ! Gather packed buffers from every child.  Matching Isend pairs are
        ! posted in ioclient_obj::push (tags kIO_TAG_WRITE_3D / _2D).
        ! Restart steps receive the full union-sized buffer flat; non-restart
        ! steps receive just the leading-dim output slice via the prebuilt
        ! vector type (recv_type_*_out) — matching the client's conditional
        ! send size.
        do cc = 1, this%n_children
            if (should_write_restart .or. this%first_write) then
                call MPI_Irecv(this%write_buffer_3d(cc)%buff, &
                                size(this%write_buffer_3d(cc)%buff), MPI_REAL, &
                                cc - 1, kIO_TAG_WRITE_3D, this%client_comms, &
                                reqs(2*cc - 1), ierr)
                call MPI_Irecv(this%write_buffer_2d(cc)%buff, &
                                size(this%write_buffer_2d(cc)%buff), MPI_REAL, &
                                cc - 1, kIO_TAG_WRITE_2D, this%client_comms, &
                                reqs(2*cc), ierr)
            else
                call MPI_Irecv(this%write_buffer_3d(cc)%buff, 1, this%recv_type_3d_out, &
                                cc - 1, kIO_TAG_WRITE_3D, this%client_comms, &
                                reqs(2*cc - 1), ierr)
                call MPI_Irecv(this%write_buffer_2d(cc)%buff, 1, this%recv_type_2d_out, &
                                cc - 1, kIO_TAG_WRITE_2D, this%client_comms, &
                                reqs(2*cc), ierr)
            endif
        enddo
        call MPI_Waitall(2 * this%n_children, reqs, MPI_STATUSES_IGNORE, ierr)

        ! Unpack into outputer buffers (two-pass)
        ! Pass 1: output variables
        do i=1,this%n_children
            n_3d = 0; n_2d = 0
            do oi = 1, size(this%out_ordered_indices)
                n = this%out_ordered_indices(oi)
                x_stag = this%outputer%var_meta(n)%xstag
                y_stag = this%outputer%var_meta(n)%ystag
                if (this%outputer%var_meta(n)%two_d) then
                    n_2d = n_2d + 1
                    this%outputer%variables(n)%data_2d((this%iswc(i)-this%i_s_w+1):(this%iewc(i)-this%i_s_w+1+x_stag),(this%jswc(i)-this%j_s_w+1):(this%jewc(i)-this%j_s_w+1+y_stag)) = &
                        this%write_buffer_2d(i)%buff(n_2d,1:(this%iewc(i)-this%iswc(i)+1+x_stag),1:(this%jewc(i)-this%jswc(i)+1+y_stag))
                elseif (this%outputer%var_meta(n)%three_d) then
                    n_3d = n_3d + 1
                    this%outputer%variables(n)%data_3d((this%iswc(i)-this%i_s_w+1):(this%iewc(i)-this%i_s_w+1+x_stag),:,(this%jswc(i)-this%j_s_w+1):(this%jewc(i)-this%j_s_w+1+y_stag)) = &
                        this%write_buffer_3d(i)%buff(n_3d,1:(this%iewc(i)-this%iswc(i)+1+x_stag),1:this%outputer%var_meta(n)%dim_len(2),1:(this%jewc(i)-this%jswc(i)+1+y_stag))
                endif
            enddo

            ! Pass 2: restart-only variables (only if restart step or first write)
            if (should_write_restart .or. this%first_write) then
                do oi = 1, size(this%rst_only_ordered_indices)
                    n = this%rst_only_ordered_indices(oi)
                    x_stag = this%outputer%var_meta(n)%xstag
                    y_stag = this%outputer%var_meta(n)%ystag
                    if (this%outputer%var_meta(n)%two_d) then
                        n_2d = n_2d + 1
                        this%outputer%variables(n)%data_2d((this%iswc(i)-this%i_s_w+1):(this%iewc(i)-this%i_s_w+1+x_stag),(this%jswc(i)-this%j_s_w+1):(this%jewc(i)-this%j_s_w+1+y_stag)) = &
                            this%write_buffer_2d(i)%buff(n_2d,1:(this%iewc(i)-this%iswc(i)+1+x_stag),1:(this%jewc(i)-this%jswc(i)+1+y_stag))
                    elseif (this%outputer%var_meta(n)%three_d) then
                        n_3d = n_3d + 1
                        this%outputer%variables(n)%data_3d((this%iswc(i)-this%i_s_w+1):(this%iewc(i)-this%i_s_w+1+x_stag),:,(this%jswc(i)-this%j_s_w+1):(this%jewc(i)-this%j_s_w+1+y_stag)) = &
                            this%write_buffer_3d(i)%buff(n_3d,1:(this%iewc(i)-this%iswc(i)+1+x_stag),1:this%outputer%var_meta(n)%dim_len(2),1:(this%jewc(i)-this%jswc(i)+1+y_stag))
                    endif
                enddo
            endif
        enddo
        if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",this%nest_indx," --------------"
        if (STD_OUT_PE_IO) write(*,*) "Done data from child images for output"
        if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",this%nest_indx," --------------"

        if (ALL(this%write_buffer_3d(1)%buff(1:this%n_out_3d,:,:,:)==kEMPT_BUFF) .or. &
            ALL(this%write_buffer_2d(1)%buff(1:this%n_out_2d,:,:)==kEMPT_BUFF)) then
            stop 'Error, all of write buffer used for output was still set to empty buffer flag at time of writing.'
        endif

        if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",this%nest_indx," --------------"
        if (STD_OUT_PE_IO) write(*,*) "Writing out file"
        if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",this%nest_indx," --------------"

        call this%outputer%save_out_file(this%sim_time,this%IO_Comms,this%out_var_indices)
        if (STD_OUT_PE_IO) write(*,"(/ A23,I2,A16)") "-------------- IOserver",this%nest_indx," --------------"
        if (STD_OUT_PE_IO) write(*,*) "Done writing out file"
        if (STD_OUT_PE_IO) write(*,"(A23,I2,A16 /)") "-------------- IOserver",this%nest_indx," --------------"


        if (should_write_restart) then
            call this%outputer%save_rst_file(this%sim_time,this%IO_Comms,this%rst_var_indices,this%dt)
            this%restart_counter = 1
        endif

        if (this%first_write) this%first_write = .false.

        call this%increment_output_time()
        this%restart_counter = this%restart_counter + 1
    end subroutine 

    
    ! This subroutine calls the read file function from the input object
    ! and then passes the read-in data to the read buffer
    module subroutine read_file(this)
        class(ioserver_t), intent(inout) :: this

        !See if we even have files to read
        if (this%files_to_read) then
            ! read file into buffer array
            call this%reader%read_next_step(this%forcing_buffer,this%IO_Comms)
            this%files_to_read = .not.(this%reader%eof)

            ! Loop through child images and send chunks of buffer array to each one
            call this%scatter_forcing()
        endif

        ! increment input time whether we have files to read or not
        ! this is because input time is just controlling if we are at an input event -- not necessarily
        ! if we read in data or not
        ! the reader object will handle if we try to read a time step which doesn't exist in the forcing
        ! call this%increment_input_time()
    end subroutine 

    module subroutine gather_forcing(this)
        class(ioserver_t), intent(inout) :: this

        integer :: i, ierr

        if (this%n_f <= 0) then
            if (STD_OUT_PE) write(*,*) 'No forcing fields to gather, but we entered gather_forcing. This is a bug.'
            return
        endif

        ! Persistent-request pattern: by the time we re-enter this routine,
        ! a previous call's tail has already pre-posted the Irecvs we need
        ! to wait on. On the very first call, post them just-in-time so the
        ! bootstrap behaves identically to the original code. Matching
        ! Isend is in ioclient::update_nest.
        if (.not. this%gather_posted) then
            do i = 1, this%n_children
                call MPI_Irecv(this%child_gather_buffers(i)%buff, &
                               size(this%child_gather_buffers(i)%buff), MPI_REAL, &
                               i - 1, kIO_TAG_NEST_3D, this%client_comms, this%gather_reqs(i), ierr)
            enddo
            this%gather_posted = .true.
        endif

        call MPI_Waitall(this%n_children, this%gather_reqs, MPI_STATUSES_IGNORE, ierr)

        ! Unpack received child buffers into the server's gather buffer
        do i=1,this%n_children
            this%gather_buffer(:,this%iswc(i):this%iewc(i)+1,this%k_s_w:this%k_e_f,this%jswc(i):this%jewc(i)+1) = &
                this%child_gather_buffers(i)%buff(:,1:(this%iewc(i)-this%iswc(i)+2),this%k_s_w:this%k_e_f,1:(this%jewc(i)-this%jswc(i)+2))
        enddo

        ! Pre-post the next round of Irecvs now so the parent's next
        ! update_nest Isend finds a matching receive already in flight.
        ! Safe to overlap with whatever the IO server does next because
        ! child_gather_buffers won't be touched again until the matching
        ! Waitall at the top of the next gather_forcing call.
        do i = 1, this%n_children
            call MPI_Irecv(this%child_gather_buffers(i)%buff, &
                           size(this%child_gather_buffers(i)%buff), MPI_REAL, &
                           i - 1, kIO_TAG_NEST_3D, this%client_comms, this%gather_reqs(i), ierr)
        enddo

    end subroutine gather_forcing

    module subroutine gather_forcing_2d(this)
        class(ioserver_t), intent(inout) :: this
        integer :: i, ierr
        type(MPI_Request), allocatable :: reqs(:)

        if (this%n_f_2d <= 0) return

        allocate(reqs(this%n_children))
        do i = 1, this%n_children
            call MPI_Irecv(this%child_gather_buffers_2d(i)%buff, &
                           size(this%child_gather_buffers_2d(i)%buff), MPI_REAL, &
                           i - 1, kIO_TAG_NEST_2D, this%client_comms, reqs(i), ierr)
        enddo
        call MPI_Waitall(this%n_children, reqs, MPI_STATUSES_IGNORE, ierr)
        deallocate(reqs)

        do i = 1, this%n_children
            this%gather_buffer_2d(:, this%iswc(i):this%iewc(i)+1, this%jswc(i):this%jewc(i)+1) = &
                this%child_gather_buffers_2d(i)%buff(:, 1:(this%iewc(i)-this%iswc(i)+2), 1:(this%jewc(i)-this%jswc(i)+2))
        enddo

    end subroutine gather_forcing_2d

    module subroutine gather_forcing_3d_init(this)
        class(ioserver_t), intent(inout) :: this
        integer :: i, ierr
        type(MPI_Request), allocatable :: reqs(:)

        if (this%n_f_3d_init <= 0) return

        allocate(reqs(this%n_children))
        do i = 1, this%n_children
            call MPI_Irecv(this%child_gather_buffers_3d_init(i)%buff, &
                           size(this%child_gather_buffers_3d_init(i)%buff), MPI_REAL, &
                           i - 1, kIO_TAG_NEST_3D_INIT, this%client_comms, reqs(i), ierr)
        enddo
        call MPI_Waitall(this%n_children, reqs, MPI_STATUSES_IGNORE, ierr)
        deallocate(reqs)

        ! The parent's own init-3D data only extends to nz_w k-levels.
        do i = 1, this%n_children
            this%gather_buffer_3d_init(:, this%iswc(i):this%iewc(i)+1, 1:this%nz_w, this%jswc(i):this%jewc(i)+1) = &
                this%child_gather_buffers_3d_init(i)%buff(:, 1:(this%iewc(i)-this%iswc(i)+2), 1:this%nz_w, 1:(this%jewc(i)-this%jswc(i)+2))
        enddo

    end subroutine gather_forcing_3d_init

    module subroutine distribute_forcing(this, child_ioserver, child_indx)
        class(ioserver_t), intent(inout) :: this
        class(ioserver_t), intent(inout)    :: child_ioserver
        integer, intent(in) :: child_indx

        integer :: i, nx, ny, n, ierr, msg_size, real_size
        integer, allocatable :: send_msg_size_alltoall(:), buff_msg_size_alltoall(:), disp_alltoall(:)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lowerbound, extent

        allocate(send_msg_size_alltoall(this%n_servers))
        allocate(buff_msg_size_alltoall(this%n_servers))
        allocate(disp_alltoall(this%n_servers))

        disp_alltoall = 0

        ! If this is the first time calling gather_forcing, we are still in initialization. Call setup_nest_types now, passing in the child ioserver
        if (this%nest_types_initialized(child_indx) .eqv. .False.) then
            if (allocated(child_ioserver%forcing_buffer))  deallocate(child_ioserver%forcing_buffer)

            allocate(child_ioserver%forcing_buffer(this%n_f,child_ioserver%i_s_r:child_ioserver%i_e_r+1,child_ioserver%k_s_r:child_ioserver%k_e_r, child_ioserver%j_s_r:child_ioserver%j_e_r+1))
            call this%setup_nest_types(child_ioserver, child_indx, &
                this%send_nest_types(child_indx,:), this%buffer_nest_types(child_indx,:))
            this%nest_types_initialized(child_indx) = .true.
            call test_nest_types(this, child_ioserver, child_indx)
        endif

        ! What would be smart here, is to send our forcing buffer to only the processes, on which the child ioservers need the data
        ! The forcing array that we need has the extent child_ioserver%i_r_s:child_ioserver%i_r_e, child_ioserver%j_r_s:child_ioserver%j_r_e
        ! As a good parent, we will provide for our child, and gather this data from all ioservers

        ! Get size of an MPI_REAL
        call MPI_Type_size(MPI_REAL, real_size, ierr)
        do n = 1, this%n_servers
            call MPI_Type_get_extent(this%send_nest_types(child_indx,n), lowerbound, extent)
            if (extent > real_size) then
                send_msg_size_alltoall(n) = 1
            else
                send_msg_size_alltoall(n) = 0
            endif
            call MPI_Type_get_extent(this%buffer_nest_types(child_indx,n), lowerbound, extent)
            if (extent > real_size) then
                buff_msg_size_alltoall(n) = 1
            else
                buff_msg_size_alltoall(n) = 0
            endif
        enddo

        call MPI_Alltoallw(this%gather_buffer,  send_msg_size_alltoall, disp_alltoall, this%send_nest_types(child_indx,:), &
                        child_ioserver%forcing_buffer, buff_msg_size_alltoall, disp_alltoall, this%buffer_nest_types(child_indx,:), this%IO_Comms)

        ! This call will scatter the forcing fields to the ioclients of the nest child
        call child_ioserver%scatter_forcing()
    end subroutine

    ! One-time init transfer for a parent nest. Resolves the per-family
    ! nz_init_3d (max across this parent and its children) on this rank,
    ! grows this%gather_buffer_3d_init if needed, drives gather_forcing_2d
    ! and gather_forcing_3d_init, then pushes the init 2D + 3D restart vars
    ! to every child via MPI_Alltoallw + fire-and-forget MPI_Isend (Request_free per send). Called once per parent
    ! from flow_events.component_read.
    !
    ! Gating (set in init_ioserver per the three-condition rule):
    !   nest_init_send_done == .true. -> early return: no child of this nest
    !                                    needs the transfer.
    ! Per-child loop body cycles when the child's own nest_init_recv_done is
    ! set so that we don't Isend to a child whose compute side will skip the
    ! recv (which would leave the parent's MPI_Waitall hanging).
    module subroutine distribute_init_forcing(this, components, child_nests)
        class(ioserver_t), intent(inout) :: this
        type(comp_arr_t),  intent(inout) :: components(:)
        integer,           intent(in)    :: child_nests(:)

        integer :: n, ci, nx_c, ny_c, nz_c, msg_size_2d, msg_size_3d
        integer :: family_nz, ierr, real_size
        integer, allocatable :: send_2d_sizes(:), buff_2d_sizes(:), disp_2d(:)
        integer, allocatable :: send_3d_sizes(:), buff_3d_sizes(:), disp_3d(:)
        type(MPI_Request) :: req
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lowerbound, extent

        if (this%nest_init_send_done) return

        ! 1. Per-family max across this parent and its children.
        family_nz = this%nz_init_3d
        do n = 1, size(child_nests)
            select type (child => components(child_nests(n))%comp)
            type is (ioserver_t)
                if (child%nz_init_3d > family_nz) family_nz = child%nz_init_3d
            end select
        enddo

        ! 2. Grow this%gather_buffer_3d_init to family_nz if needed. Safe:
        ! this runs before the first gather_forcing_3d_init, so no data
        ! would be lost.
        if (family_nz > this%nz_init_3d .and. this%n_f_3d_init > 0) then
            if (allocated(this%gather_buffer_3d_init)) deallocate(this%gather_buffer_3d_init)
            allocate(this%gather_buffer_3d_init(this%n_f_3d_init, &
                     this%i_s_w:this%i_e_w+1, 1:family_nz, this%j_s_w:this%j_e_w+1))
            this%gather_buffer_3d_init = 0.0
        endif

        ! 3. Write family_nz into this parent and every child's ioserver on
        ! this rank. Direct struct access is sufficient: every IO rank runs
        ! this same code over the same components array and arrives at the
        ! same result, so all IO ranks stay consistent without an MPI
        ! reduction. Compute ranks learn the size from the actual incoming
        ! message in unpack_init_vars_3d (MPI_Probe), not from this field.
        this%nz_init_3d = family_nz
        do n = 1, size(child_nests)
            select type (child => components(child_nests(n))%comp)
            type is (ioserver_t)
                child%nz_init_3d = family_nz
            end select
        enddo

        ! 4+5. One-time init gathers.
        call this%gather_forcing_2d()
        call this%gather_forcing_3d_init()

        ! 6. Distribute init 2D + 3D to each child.
        call MPI_Type_size(MPI_REAL, real_size, ierr)

        do n = 1, size(child_nests)
            select type (child_ioserver => components(child_nests(n))%comp)
            type is (ioserver_t)
                ! Skip children that don't need the transfer (e.g. they're in
                ! restart mode, or share start_time with this parent). All
                ! parent IO ranks read the same options(:) so they reach the
                ! same skip decision per child, keeping any inner Alltoallw
                ! collective consistent.
                if (child_ioserver%nest_init_recv_done) cycle

                ! Per-child-rank sends are MPI_Isend + MPI_Request_free
                ! (fire-and-forget). A blocking Waitall here would deadlock
                ! when this parent has multiple children: compute processes
                ! children sequentially through component_loop, and between
                ! two children's wakes compute does push (nest N) which
                ! needs IO's component_write (nest N) -- but component_write
                ! only runs after distribute_init_forcing returns. The
                ! buffers (forcing_buffer_2d/3d_init) are allocated once
                ! per child (gated by nest_types_*_initialized below) and
                ! persist for the run, so the slices stay valid for the
                ! Isends. The per-child cycle above guarantees every freed
                ! request has a matching MPI_Recv coming from the child's
                ! receive_nest_init.

                ! 2D init distribute
                if (this%n_f_2d > 0) then
                    if (.not. this%nest_types_2d_initialized(n)) then
                        if (allocated(child_ioserver%forcing_buffer_2d)) deallocate(child_ioserver%forcing_buffer_2d)
                        allocate(child_ioserver%forcing_buffer_2d(this%n_f_2d, &
                            child_ioserver%i_s_r:child_ioserver%i_e_r+1, child_ioserver%j_s_r:child_ioserver%j_e_r+1))
                        call this%setup_nest_types_2d(child_ioserver, n, &
                            this%send_nest_types_2d(n,:), this%buffer_nest_types_2d(n,:))
                        ! Per-child-rank contiguous send buffers. The slice
                        ! forcing_buffer_2d(:, isrc:ierc+1, jsrc:jerc+1) is
                        ! non-contiguous when this IO server hosts more than
                        ! one compute rank; basic-count Isend over a non-
                        ! contiguous slice is unreliable under mpi_f08+NVHPC.
                        ! Copy into these contiguous buffers, then Isend.
                        if (.not. allocated(child_ioserver%client_send_buffers_2d)) &
                            allocate(child_ioserver%client_send_buffers_2d(child_ioserver%n_children))
                        do ci = 1, child_ioserver%n_children
                            nx_c = child_ioserver%ierc(ci) - child_ioserver%isrc(ci) + 2
                            ny_c = child_ioserver%jerc(ci) - child_ioserver%jsrc(ci) + 2
                            allocate(child_ioserver%client_send_buffers_2d(ci)%buff(this%n_f_2d, nx_c, ny_c))
                        enddo
                        this%nest_types_2d_initialized(n) = .true.
                    endif

                    allocate(send_2d_sizes(this%n_servers), buff_2d_sizes(this%n_servers), disp_2d(this%n_servers))
                    disp_2d = 0
                    do ci = 1, this%n_servers
                        call MPI_Type_get_extent(this%send_nest_types_2d(n,ci), lowerbound, extent)
                        send_2d_sizes(ci) = merge(1, 0, extent > real_size)
                        call MPI_Type_get_extent(this%buffer_nest_types_2d(n,ci), lowerbound, extent)
                        buff_2d_sizes(ci) = merge(1, 0, extent > real_size)
                    enddo

                    call MPI_Alltoallw(this%gather_buffer_2d, send_2d_sizes, disp_2d, this%send_nest_types_2d(n,:), &
                        child_ioserver%forcing_buffer_2d, buff_2d_sizes, disp_2d, this%buffer_nest_types_2d(n,:), this%IO_Comms)
                    deallocate(send_2d_sizes, buff_2d_sizes, disp_2d)

                    ! Fire-and-forget Isend of 2D data to each child ioclient.
                    ! Copy the per-rank slice into a contiguous buffer first
                    ! (Fortran assignment handles strides correctly), then
                    ! Isend the contiguous buffer.
                    do ci = 1, child_ioserver%n_children
                        nx_c = child_ioserver%ierc(ci) - child_ioserver%isrc(ci) + 2
                        ny_c = child_ioserver%jerc(ci) - child_ioserver%jsrc(ci) + 2
                        msg_size_2d = this%n_f_2d * nx_c * ny_c
                        child_ioserver%client_send_buffers_2d(ci)%buff(:, 1:nx_c, 1:ny_c) = &
                            child_ioserver%forcing_buffer_2d(:, &
                                child_ioserver%isrc(ci):child_ioserver%ierc(ci)+1, &
                                child_ioserver%jsrc(ci):child_ioserver%jerc(ci)+1)
                        call MPI_Isend(child_ioserver%client_send_buffers_2d(ci)%buff, &
                            msg_size_2d, MPI_REAL, ci-1, 98, child_ioserver%client_comms, &
                            req, ierr)
                        call MPI_Request_free(req, ierr)
                    enddo
                endif

                ! 3D init distribute. Step 3 above ensures this%nz_init_3d ==
                ! child_ioserver%nz_init_3d, so send/buff derived types and
                ! both physical buffers share the same k extent.
                if (this%n_f_3d_init > 0) then
                    if (.not. this%nest_types_3d_init_initialized(n)) then
                        if (allocated(child_ioserver%forcing_buffer_3d_init)) deallocate(child_ioserver%forcing_buffer_3d_init)
                        allocate(child_ioserver%forcing_buffer_3d_init(this%n_f_3d_init, &
                            child_ioserver%i_s_r:child_ioserver%i_e_r+1, 1:child_ioserver%nz_init_3d, &
                            child_ioserver%j_s_r:child_ioserver%j_e_r+1))
                        call this%setup_nest_types_3d_init(child_ioserver, n, &
                            this%send_nest_types_3d_init(n,:), this%buffer_nest_types_3d_init(n,:))
                        ! Per-child-rank contiguous send buffers (see 2D
                        ! block above for rationale).
                        if (.not. allocated(child_ioserver%client_send_buffers_3d_init)) &
                            allocate(child_ioserver%client_send_buffers_3d_init(child_ioserver%n_children))
                        do ci = 1, child_ioserver%n_children
                            nx_c = child_ioserver%ierc(ci) - child_ioserver%isrc(ci) + 2
                            ny_c = child_ioserver%jerc(ci) - child_ioserver%jsrc(ci) + 2
                            allocate(child_ioserver%client_send_buffers_3d_init(ci)%buff( &
                                this%n_f_3d_init, nx_c, child_ioserver%nz_init_3d, ny_c))
                        enddo
                        this%nest_types_3d_init_initialized(n) = .true.
                    endif

                    allocate(send_3d_sizes(this%n_servers), buff_3d_sizes(this%n_servers), disp_3d(this%n_servers))
                    disp_3d = 0
                    do ci = 1, this%n_servers
                        call MPI_Type_get_extent(this%send_nest_types_3d_init(n,ci), lowerbound, extent)
                        send_3d_sizes(ci) = merge(1, 0, extent > real_size)
                        call MPI_Type_get_extent(this%buffer_nest_types_3d_init(n,ci), lowerbound, extent)
                        buff_3d_sizes(ci) = merge(1, 0, extent > real_size)
                    enddo

                    call MPI_Alltoallw(this%gather_buffer_3d_init, send_3d_sizes, disp_3d, this%send_nest_types_3d_init(n,:), &
                        child_ioserver%forcing_buffer_3d_init, buff_3d_sizes, disp_3d, this%buffer_nest_types_3d_init(n,:), this%IO_Comms)
                    deallocate(send_3d_sizes, buff_3d_sizes, disp_3d)

                    ! Fire-and-forget Isend of 3D init data to each child
                    ! ioclient. Use the full k extent (1:nz_init_3d) so
                    ! variables with per-var nz_v larger than atmospheric
                    ! (snow/soil) are transmitted in their entirety. Copy
                    ! the per-rank slice into a contiguous buffer first.
                    do ci = 1, child_ioserver%n_children
                        nx_c = child_ioserver%ierc(ci) - child_ioserver%isrc(ci) + 2
                        nz_c = child_ioserver%nz_init_3d
                        ny_c = child_ioserver%jerc(ci) - child_ioserver%jsrc(ci) + 2
                        msg_size_3d = this%n_f_3d_init * nx_c * nz_c * ny_c
                        child_ioserver%client_send_buffers_3d_init(ci)%buff(:, 1:nx_c, 1:nz_c, 1:ny_c) = &
                            child_ioserver%forcing_buffer_3d_init(:, &
                                child_ioserver%isrc(ci):child_ioserver%ierc(ci)+1, &
                                1:child_ioserver%nz_init_3d, &
                                child_ioserver%jsrc(ci):child_ioserver%jerc(ci)+1)
                        call MPI_Isend(child_ioserver%client_send_buffers_3d_init(ci)%buff, &
                            msg_size_3d, MPI_REAL, ci-1, 99, child_ioserver%client_comms, &
                            req, ierr)
                        call MPI_Request_free(req, ierr)
                    enddo
                endif

                ! All three datatype variants for this child have now been
                ! committed (main forcing in distribute_forcing before this
                ! call, plus init 2D and init 3D just above). The cached
                ! (i,j) cell lists are no longer needed.
                call free_nest_geometry(this%nest_geometry(n))
            end select
        enddo

        this%nest_init_send_done = .true.
    end subroutine distribute_init_forcing

    module subroutine scatter_forcing(this)
        class(ioserver_t), intent(inout) :: this

        integer :: i, ierr
        type(MPI_Request) :: reqs(this%n_children)

        ! Pack each child's slice into its send buffer and post one Isend
        ! per child on kIO_TAG_READ.  Children post matching Irecv in
        ! ioclient::receive.
        do i = 1, this%n_children
            this%read_buffer(i)%buff(:, 1:(this%ierc(i)-this%isrc(i)+1)+1, &
                                      :, 1:(this%jerc(i)-this%jsrc(i)+1)+1) = &
                this%forcing_buffer(:, this%isrc(i):this%ierc(i)+1, &
                                     :, this%jsrc(i):this%jerc(i)+1)
            call MPI_Isend(this%read_buffer(i)%buff, &
                           size(this%read_buffer(i)%buff), MPI_REAL, &
                           i - 1, kIO_TAG_READ, this%client_comms, reqs(i), ierr)
        enddo
        call MPI_Waitall(this%n_children, reqs, MPI_STATUSES_IGNORE, ierr)

    end subroutine

    ! This subroutine reads in the restart file and then sends the data to the children
    ! Same as above, but for restart file
    module subroutine read_restart_file(this, options)
        class(ioserver_t),   intent(inout) :: this
        type(options_t),     intent(in)    :: options

        integer :: i, n_3d, n_2d, nx, ny, i_s_re, i_e_re, j_s_re, j_e_re
        integer :: ncid, file_var_id, dimid_3d(4), nz, err, varid, start_3d(4), cnt_3d(4), start_2d(3), cnt_2d(3)
        integer :: nx_c, ny_c, nz_c, ierr

        type(MPI_Info) :: IO_Comms_info
        real, allocatable :: data3d(:,:,:,:), restart_data_3d(:,:,:,:), restart_data_2d(:,:,:)
        type(meta_data_t)  :: var_meta
        character(len=kMAX_NAME_LENGTH) :: name
        character(len=kMAX_FILE_LENGTH) :: base_rst_file_name, tmp_str, restart_in_file
        character(len=kMAX_STRING_LENGTH), allocatable :: tokens(:)

        integer :: restart_step                         ! time step relative to the start of the restart file
        type(Time_type) :: time_at_step   ! restart date as a modified julian day
        type(MPI_Request), allocatable :: reqs(:)
        type(rst_pack_t), allocatable :: packs(:)

        allocate(restart_data_3d(this%n_w_3d,this%i_s_re:this%i_e_re+1,this%k_s_w:this%k_e_w,this%j_s_re:this%j_e_re+1))
        allocate(restart_data_2d(this%n_w_2d,this%i_s_re:this%i_e_re+1,                      this%j_s_re:this%j_e_re+1))

        tokens = split_str(trim(options%domain%init_conditions_file), "/")

        ! store the base file name temporarily here.
        tmp_str = trim(tokens(size(tokens)))

        ! Set the output filename to be the location of the output folder, and the name of the initial conditions file, with ".nc" removed
        base_rst_file_name = trim(options%restart%restart_folder) // '/' // tmp_str(1:(len_trim(tmp_str)-3)) // "_"

        write(restart_in_file, '(A,A,".nc")')    &
                trim(base_rst_file_name),   &
                trim(as_string(options%restart%restart_time,'(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'))

        call check_file_exists(restart_in_file, message='Restart file does not exist.')

        ! find the time step that most closely matches the requested restart time (<=)
        restart_step = find_timestep_in_file(restart_in_file, 'time', options%restart%restart_time, time_at_step)

        if (options%general%debug) then
            write(*,*) " ------------------ "
            write(*,*) "RESTART INFORMATION"
            write(*,*) "mjd",         options%restart%restart_time%mjd()
            write(*,*) "date:",       trim(as_string(options%restart%restart_time))
            write(*,*) "file:",   trim(restart_in_file)
            write(*,*) "forcing step",options%restart%restart_step_in_file
            write(*,*) " ------------------ "
        endif

        call MPI_Comm_get_info(this%IO_Comms, IO_Comms_info)

        call check_ncdf(nf90_open(restart_in_file, IOR(nf90_nowrite,NF90_NETCDF4), ncid, &
                comm = this%IO_Comms%MPI_VAL, info = IO_Comms_info%MPI_VAL), " Opening file "//trim(restart_in_file))

        ! Validate that the restart file config matches the current config
        call compare_restart_config(ncid, options)

        ! Read saved dt for restart reproducibility (graceful fallback for old restart files)
        err = nf90_get_att(ncid, NF90_GLOBAL, "dt_seconds", this%restart_dt)
        if (err /= NF90_NOERR) this%restart_dt = 0.0

        ! setup start/count arrays accordingly. k_s_w and k_e_w forseen to always cover the bounds of what k_s_re and k_e_re would be
        ! This is because the domain is only decomposed in 2D.
        start_3d = (/ this%i_s_re,this%j_s_re,this%k_s_w,restart_step /)
        start_2d = (/ this%i_s_re,this%j_s_re,restart_step /)
        cnt_3d = (/ (this%i_e_re-this%i_s_re+1),(this%j_e_re-this%j_s_re+1),(this%k_e_w-this%k_s_w+1),1 /)
        cnt_2d = (/ (this%i_e_re-this%i_s_re+1),(this%j_e_re-this%j_s_re+1),1 /)

        n_2d = 1
        n_3d = 1

        do i = 1,size(this%rst_var_indices)
            var_meta = this%outputer%var_meta(this%rst_var_indices(i))
            name = var_meta%name
            call check_ncdf( nf90_inq_varid(ncid, name, file_var_id), " Getting var ID for "//trim(name))
            call check_ncdf( nf90_var_par_access(ncid, file_var_id, nf90_collective))
            
            nx = cnt_3d(1) + var_meta%xstag
            ny = cnt_3d(2) + var_meta%ystag

            if (var_meta%three_d) then
                ! Get length of z dim
                call check_ncdf( nf90_inquire_variable(ncid, file_var_id, dimids = dimid_3d), " Getting dim IDs for "//trim(name))
                call check_ncdf( nf90_inquire_dimension(ncid, dimid_3d(3), len = nz), " Getting z dim len for "//trim(name))

                !clip nz to be the extent needed, or the full extent of the variable in the file. This clip is a dummy check in case
                ! the user changed the nz value for the domain and didn't delete the old restart file. Speaking from experience.
                nz = min(nz, (this%k_e_w-this%k_s_w+1) )

                
                if (allocated(data3d)) deallocate(data3d)
                allocate(data3d(nx,ny,nz,1))
                call check_ncdf( nf90_get_var(ncid, file_var_id, data3d, start=start_3d, count=(/ nx, ny, nz /)), " Getting 3D var "//trim(name))

                restart_data_3d(n_3d,this%i_s_re:this%i_e_re+var_meta%xstag,1:nz,this%j_s_re:this%j_e_re+var_meta%ystag) = &
                        reshape(data3d(:,:,:,1), shape=[nx,nz,ny], order=[1,3,2])
                n_3d = n_3d+1
            else if (var_meta%two_d) then
                call check_ncdf( nf90_get_var(ncid, file_var_id, restart_data_2d(n_2d,this%i_s_re:this%i_e_re+var_meta%xstag,this%j_s_re:this%j_e_re+var_meta%ystag), &
                        start=start_2d, count=(/ nx, ny /)), " Getting 2D "//trim(name))
                        n_2d = n_2d+1
            endif
        end do
        
        call check_ncdf(nf90_close(ncid), "Closing file "//trim(restart_in_file))

        ! Scatter restart slices to each child via flat Isend. One
        ! per-child pack holds a contiguous copy of that child's
        ! restart_data_{3d,2d} sub-block; clients post matching Irecv in
        ! ioclient_obj::receive_rst (tags kIO_TAG_RST_3D / kIO_TAG_RST_2D).
        allocate(packs(this%n_children))
        allocate(reqs(2 * this%n_children))
        nz_c = this%k_e_w - this%k_s_w + 1
        do i = 1, this%n_children
            nx_c = this%ierec(i) - this%isrec(i) + 2
            ny_c = this%jerec(i) - this%jsrec(i) + 2
            allocate(packs(i)%buff_3d(this%n_w_3d, nx_c, nz_c, ny_c))
            allocate(packs(i)%buff_2d(this%n_w_2d, nx_c,       ny_c))
            packs(i)%buff_3d = restart_data_3d(1:this%n_w_3d, &
                                                this%isrec(i):this%isrec(i) + nx_c - 1, &
                                                this%k_s_w:this%k_e_w, &
                                                this%jsrec(i):this%jsrec(i) + ny_c - 1)
            packs(i)%buff_2d = restart_data_2d(1:this%n_w_2d, &
                                                this%isrec(i):this%isrec(i) + nx_c - 1, &
                                                this%jsrec(i):this%jsrec(i) + ny_c - 1)
            call MPI_Isend(packs(i)%buff_3d, size(packs(i)%buff_3d), MPI_REAL, &
                            i - 1, kIO_TAG_RST_3D, this%client_comms, reqs(2*i - 1), ierr)
            call MPI_Isend(packs(i)%buff_2d, size(packs(i)%buff_2d), MPI_REAL, &
                            i - 1, kIO_TAG_RST_2D, this%client_comms, reqs(2*i),     ierr)
        enddo
        call MPI_Waitall(2 * this%n_children, reqs, MPI_STATUSES_IGNORE, ierr)
        ! packs + reqs auto-deallocate on block exit

    end subroutine 


    ! This function closes all open file handles. Files are left open by default
    ! to minimize I/O calls. When the program exits, this must be called
    module subroutine close_files(this)
        class(ioserver_t), intent(inout) :: this
                
        ! close files
        call this%reader%close_file()
        call this%outputer%close_output_files()

    end subroutine


    ! =========================================================================
    ! Restart config validation
    ! =========================================================================
    subroutine compare_restart_config(ncid, options)
        integer,         intent(in) :: ncid
        type(options_t), intent(in) :: options

        character(len=kMAX_CONFIG_STRING_LENGTH) :: config_str
        character(len=512) :: line
        character(len=256) :: key, current_val, attr_name
        character(len=kMAX_ATTR_LENGTH) :: stored_val
        character(len=2512) :: diff_report(500)
        integer :: i, line_start, str_len, eq_pos, n_diffs, slash_pos
        integer :: status, attr_len
        logical :: first_key
        integer :: j

        ! Generate current config (excluding restart-specific fields)
        call options%generate_config_string(config_str, exclude_restart_fields=.true.)

        ! Parse config string and compare each key against stored attributes
        n_diffs = 0
        line_start = 1
        str_len = len_trim(config_str)
        first_key = .true.

        do i = 1, str_len
            if (config_str(i:i) == char(10) .or. i == str_len) then
                if (config_str(i:i) == char(10)) then
                    line = config_str(line_start:i-1)
                else
                    line = config_str(line_start:i)
                endif
                line_start = i + 1

                if (len_trim(line) == 0) cycle

                eq_pos = index(line, '=')
                if (eq_pos <= 0) cycle

                key = adjustl(line(1:eq_pos-1))
                current_val = adjustl(line(eq_pos+1:))

                ! Build NetCDF-safe attribute name (replace '/' with '.')
                attr_name = key
                slash_pos = index(attr_name, '/')
                if (slash_pos > 0) attr_name(slash_pos:slash_pos) = '.'

                ! Check if this key should be excluded from comparison
                ! do j = 1, N_EXCLUDE
                !     if (index(trim(key), trim(exclude_prefixes(j))) == 1) cycle
                ! end do

                ! On the first key, check if cfg.* attributes exist at all (backward compat)
                if (first_key) then
                    first_key = .false.
                    status = nf90_inquire_attribute(ncid, NF90_GLOBAL, trim(attr_name), len=attr_len)
                    if (status == NF90_ENOTATT) then
                        write(*,*) "  WARNING: Restart file does not contain config attributes."
                        write(*,*) "           Unable to verify configuration consistency. Continuing..."
                        return
                    endif
                endif

                ! Read stored attribute value
                status = nf90_inquire_attribute(ncid, NF90_GLOBAL, trim(attr_name), len=attr_len)
                if (status /= NF90_NOERR) then
                    ! Attribute missing from restart file
                    n_diffs = n_diffs + 1
                    if (n_diffs <= 500) then
                        write(diff_report(n_diffs), '(A,A,A,A,A)') &
                            "  ", trim(key), ": restart=<missing>, current='", &
                            trim(current_val), "'"
                    endif
                    cycle
                endif

                stored_val = ''
                call check_ncdf(nf90_get_att(ncid, NF90_GLOBAL, trim(attr_name), stored_val), &
                        " Reading " // trim(attr_name) // " from restart file")

                if (trim(stored_val) /= trim(current_val)) then
                    n_diffs = n_diffs + 1
                    if (n_diffs <= 500) then
                        write(diff_report(n_diffs), '(A,A,A,A,A,A,A)') &
                            "  ", trim(key), ": restart='", trim(stored_val), &
                            "', current='", trim(current_val), "'"
                    endif
                endif
            endif
        end do

        if ( (n_diffs > 0) .and. (options%restart%override_check .eqv. .False.)) then
            write(*,*) ""
            write(*,*) "ERROR: Configuration mismatch between restart file and current namelist"
            write(*,*) "  Number of differences: ", n_diffs
            do i = 1, min(n_diffs, 500)
                write(*,*) trim(diff_report(i))
            end do
            if (n_diffs > 500) then
                write(*,*) "  ... and ", n_diffs - 500, " more differences"
            endif
            write(*,*) ""
            write(*,*) "-------------------------------------------------------------------------------------------"
            write(*,*) "To ignore and continue, set override_check = .True. in restart section of the namelist"
            write(*,*) "-------------------------------------------------------------------------------------------"
            error stop "Configuration mismatch: restart file was generated with different options than current namelist"
        endif

    end subroutine compare_restart_config


end submodule
