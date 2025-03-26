
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
  use output_metadata,          only : get_varindx
  use string,                   only : split_str
  use io_routines,              only : check_file_exists
  use time_io,                    only : find_timestep_in_file
  implicit none

contains


    module subroutine init(this, options, nest_indx)
        class(ioserver_t),  intent(inout)  :: this
        type(options_t), intent(in)     :: options(:)
        integer,            intent(in)     :: nest_indx

        type(variable_t) :: var
        integer ::  n, n_3d, n_2d, var_indx, out_i, rst_i, some_child_id, ierr
        
        ! Call the parent type's init procedure
        call this%init_flow_obj(options(nest_indx),nest_indx)

        ! Perform a series of MPI Gather/reduction to communicate the grid dimensions of the children clients to the parent server
        call init_with_clients(this)

        !Determine number of IOServers in the this%IO_Comms communicator, i.e. size of communicator
        call MPI_Comm_size(this%IO_Comms, this%n_servers, ierr)

        !Setup reading capability
        if (options(nest_indx)%general%parent_nest == 0) then
            call this%reader%init(this%i_s_r,this%i_e_r,this%k_s_r,this%k_e_r,this%j_s_r,this%j_e_r,options(nest_indx))
            !this%n_children = size(this%children)
            this%n_r = this%reader%n_vars
            this%files_to_read = .not.(this%reader%eof)
        else
            this%n_r = count(options(nest_indx)%forcing%vars_to_read /= "")
            this%files_to_read = .False.
        endif

        this%n_child_ioservers = size(options(nest_indx)%general%child_nests)

        if (size(options(nest_indx)%general%child_nests) > 0) then
            some_child_id = options(nest_indx)%general%child_nests(1)
            this%n_f = count(options(some_child_id)%forcing%vars_to_read /= "")
            if (this%n_f > 0) then
                allocate(this%gather_buffer(this%n_f,this%i_s_w:this%i_e_w+1,this%k_s_w:this%k_e_w,this%j_s_w:this%j_e_w+1))
            endif

            allocate(this%send_nest_types(this%n_child_ioservers,this%n_servers))
            allocate(this%buffer_nest_types(this%n_child_ioservers,this%n_servers))    
        endif

        !Setup writing capability
        call this%outputer%init(options(nest_indx),this%i_s_w,this%i_e_w,this%k_s_w,this%k_e_w,this%j_s_w,this%j_e_w,this%ide,this%kde,this%jde)

        !determine if we need to increase our k index due to some very large soil field
        this%n_w_3d = 0

        ! Save the number of vertical levels in the domain, since this is what our child ioserver will expect as the z-dimension for forcing data
        this%k_e_f = options(nest_indx)%domain%nz

        do n = 1,this%outputer%n_vars
            if (this%outputer%variables(n)%three_d) then
                this%n_w_3d = this%n_w_3d+1
            endif
            if(this%outputer%variables(n)%dim_len(3) > this%k_e_w) this%k_e_w = this%outputer%variables(n)%dim_len(3)
        enddo

        this%n_w_2d = this%outputer%n_vars - this%n_w_3d

        if (this%n_w_2d == 0) write(*,*) 'Warning: No 2D variables set for output, this should never happen'
        if (this%n_w_3d == 0) write(*,*) 'Warning: No 3D variables set for output, this should never happen'

        call setup_MPI_windows(this)
        call setup_MPI_types(this)


        !Link local buffer to the outputer variables
        allocate(this%parent_write_buffer_3d(this%n_w_3d,this%i_s_w:this%i_e_w+1,this%k_s_w:this%k_e_w,this%j_s_w:this%j_e_w+1))
        allocate(this%parent_write_buffer_2d(this%n_w_2d,this%i_s_w:this%i_e_w+1,                      this%j_s_w:this%j_e_w+1))

        n_3d = 1
        n_2d = 1

        do n = 1,this%outputer%n_vars
            if (this%outputer%variables(n)%three_d) then
                this%outputer%variables(n)%data_3d => this%parent_write_buffer_3d(n_3d,:,:,:)
                n_3d = n_3d + 1
            else
                this%outputer%variables(n)%data_2d => this%parent_write_buffer_2d(n_2d,:,:)
                n_2d = n_2d + 1
            endif
        enddo

        !Setup arrays for information about accessing variables from write buffer
        allocate(this%out_var_indices(count(options(nest_indx)%output%vars_for_output > 0)))
        allocate(this%rst_var_indices(count(options(nest_indx)%vars_for_restart > 0)))

        out_i = 1
        rst_i = 1
        
        do n=1,this%outputer%n_vars
            var_indx = get_varindx(this%outputer%variables(n)%name)
            if (options(nest_indx)%output%vars_for_output(var_indx) > 0) then
                this%out_var_indices(out_i) = n
                out_i = out_i + 1
            endif
            if (options(nest_indx)%vars_for_restart(var_indx) > 0) then
                this%rst_var_indices(rst_i) = n
                rst_i = rst_i + 1
            endif
        enddo

        if (options(nest_indx)%restart%restart) call this%outputer%init_restart(options(nest_indx), this%IO_Comms, this%out_var_indices)

    end subroutine

    ! This subroutine creates MPI datatypes which are used to collect domain data from all parent ioserver processes
    ! which our child ioserver process needs
    subroutine setup_nest_types(this, child_ioserver, send_nest_types, buffer_nest_types)
        class(ioserver_t), intent(inout) :: this
        type(ioserver_t), intent(in)    :: child_ioserver
        type(MPI_Datatype), intent(out) :: send_nest_types(:), buffer_nest_types(:)


        integer :: n, my_rank, i, j, k, v, mask_size
        integer :: ierr, i_start, i_end, j_start, j_end, counter

        integer, allocatable, dimension(:) :: parent_ims, parent_ime, parent_jms, parent_jme
        integer, allocatable, dimension(:) :: child_isr, child_ier, child_jsr, child_jer, block_lengths, displacements
        real,    allocatable, dimension(:,:) :: mask
        !Determine our rank in the this%IO_Comms communicator
        call MPI_Comm_rank(this%IO_Comms, my_rank, ierr)
        !Determine number of IOServers in the this%IO_Comms communicator, i.e. size of communicator

        allocate(parent_ims(this%n_servers), parent_ime(this%n_servers), parent_jms(this%n_servers), parent_jme(this%n_servers))
        allocate(child_isr(this%n_servers), child_ier(this%n_servers), child_jsr(this%n_servers), child_jer(this%n_servers))

        parent_ims = 0; parent_ime = 0; parent_jms = 0; parent_jme = 0
        child_isr = 0; child_ier = 0; child_jsr = 0; child_jer = 0

        parent_ims(my_rank+1) = this%i_s_w
        parent_ime(my_rank+1) = this%i_e_w
        parent_jms(my_rank+1) = this%j_s_w
        parent_jme(my_rank+1) = this%j_e_w

        child_isr(my_rank+1) = child_ioserver%i_s_r
        child_ier(my_rank+1) = child_ioserver%i_e_r
        child_jsr(my_rank+1) = child_ioserver%j_s_r
        child_jer(my_rank+1) = child_ioserver%j_e_r

        !Need access to all children ioservers -- then I can construct the send_nest_types based on my domain extent and their domain extents
        call MPI_Allreduce(MPI_IN_PLACE,child_isr,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_ier,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_jsr,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,child_jer,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)

        !Need access to all ioservers -- then I can construct the buffer_nest_types, based on my child ioserver domain extents and all the parent domain extents
        call MPI_Allreduce(MPI_IN_PLACE,parent_ims,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_ime,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_jms,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,parent_jme,this%n_servers,MPI_INT,MPI_MAX,this%IO_Comms,ierr)

        call MPI_Allreduce(MPI_IN_PLACE,parent_ims,this%n_servers,MPI_INT,MPI_MIN,this%IO_Comms,ierr)
        do n = 1,this%n_servers
            ! find where we have a "block" (hole) in the domain, since the child ioclients may not give us a perfect rectangle
            
            i_start = max(child_isr(n),this%i_s_w)
            i_end = min(child_ier(n),this%i_e_w)
            j_start = max(child_jsr(n),this%j_s_w)
            j_end = min(child_jer(n),this%j_e_w)

            ! See if any of the child ioserver domain is within the parent ioserver domain
            if (i_end > i_start .and. j_end > j_start) then
                counter = count(this%parent_write_buffer_3d(1,i_start:i_end,1,j_start:j_end)/=kEMPT_BUFF)*this%n_f*(this%k_e_f-this%k_s_w+1)

                if (allocated(block_lengths)) deallocate(block_lengths)
                if (allocated(displacements)) deallocate(displacements)
                allocate(block_lengths(counter))
                allocate(displacements(counter))

                block_lengths = 1

                counter = 0
                do j = j_start, j_end
                    do k = this%k_s_w, this%k_e_f
                        do i = i_start, i_end
                            do v = 1,this%n_f
                                if (this%parent_write_buffer_3d(1,i,1,j) /= kEMPT_BUFF) then
                                    counter = counter + 1
                                    displacements(counter) = (v-1) + ((i-this%i_s_w) + (k-this%k_s_w) * (this%i_e_w-this%i_s_w+2) + &
                                                                      (j-this%j_s_w) * (this%i_e_w-this%i_s_w+2) * (this%k_e_w-this%k_s_w+1))*this%n_f
                                endif
                            enddo
                        enddo
                    enddo
                enddo
                call MPI_Type_Indexed(counter, block_lengths, displacements, MPI_REAL, send_nest_types(n))
            else
                send_nest_types(n) = MPI_REAL
            endif
            call MPI_Type_commit(send_nest_types(n))

            i_start = max(parent_ims(n),child_ioserver%i_s_r)
            i_end = min(parent_ime(n),child_ioserver%i_e_r)
            j_start = max(parent_jms(n),child_ioserver%j_s_r)
            j_end = min(parent_jme(n),child_ioserver%j_e_r)

            ! To calculate the receive buffer, we need to know the mask of the sending IO process
            ! Do a MPI_Scatter here to get the mask of the sending IO process
            if (allocated(mask)) deallocate(mask)
            allocate(mask(parent_ims(n):parent_ime(n)+1,parent_jms(n):parent_jme(n)+1))
            mask_size = (parent_ime(n)-parent_ims(n)+2)*(parent_jme(n)-parent_jms(n)+2)

            if (my_rank == n-1) then
                mask = this%parent_write_buffer_3d(1,parent_ims(n):parent_ime(n)+1,1,parent_jms(n):parent_jme(n)+1)
            endif
            call MPI_Bcast(mask, mask_size, MPI_REAL, n-1, this%IO_Comms, ierr)

            ! See if any of the child ioserver domain is within the parent ioserver domain
            if (i_end > i_start .and. j_end > j_start) then
                counter = count(mask(i_start:i_end,j_start:j_end)/=kEMPT_BUFF)*this%n_f*(child_ioserver%k_e_r-child_ioserver%k_s_r+1)

                if (allocated(block_lengths)) deallocate(block_lengths)
                if (allocated(displacements)) deallocate(displacements)
                allocate(block_lengths(counter))
                allocate(displacements(counter))

                block_lengths = 1

                counter = 0
                do j = j_start, j_end
                    do k = child_ioserver%k_s_r, child_ioserver%k_e_r
                        do i = i_start, i_end
                            do v = 1,this%n_f
                                if (mask(i,j) /= kEMPT_BUFF) then
                                    counter = counter + 1
                                    displacements(counter) = (v-1) + ((i-child_ioserver%i_s_r) + (k-child_ioserver%k_s_r) * (child_ioserver%i_e_r-child_ioserver%i_s_r+1) + &
                                                                    (j-child_ioserver%j_s_r) * (child_ioserver%i_e_r-child_ioserver%i_s_r+1) * (child_ioserver%k_e_r-child_ioserver%k_s_r+1))*this%n_f
                                endif
                            enddo
                        enddo
                    enddo
                enddo
                call MPI_Type_Indexed(counter, block_lengths, displacements, MPI_REAL, buffer_nest_types(n))
            else
                buffer_nest_types(n) = MPI_REAL
            endif
            call MPI_Type_commit(buffer_nest_types(n))
        enddo

        this%nest_types_initialized = .true.

    end subroutine setup_nest_types
    
    subroutine setup_MPI_windows(this)
        class(ioserver_t),   intent(inout)  :: this

        type(c_ptr) :: tmp_ptr
        integer(KIND=MPI_ADDRESS_KIND) :: win_size
        integer :: ierr, real_size

        CALL MPI_Type_size(MPI_REAL, real_size)

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

        call MPI_Allreduce(MPI_IN_PLACE,this%n_r,1,MPI_INT,MPI_MAX,this%client_comms,ierr)

        win_size = this%n_w_3d*this%nx_re*this%nz_w*this%ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%client_comms, tmp_ptr, this%write_win_3d)

        win_size = this%n_w_2d*this%nx_re*this%ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%client_comms, tmp_ptr, this%write_win_2d)

        if (this%n_f > 0) then
            win_size = this%n_f*this%nx_w*this%ny_w*this%nz_w
            call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%client_comms, tmp_ptr, this%nest_win)
        endif

        win_size = this%n_r*this%nx_r*this%nz_r*this%ny_r
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%client_comms, tmp_ptr, this%read_win)
    
    end subroutine setup_MPI_windows

    subroutine setup_MPI_types(this)
        class(ioserver_t),   intent(inout)  :: this

        integer :: i

        allocate(this%get_types_3d(this%n_children))
        allocate(this%get_types_2d(this%n_children))
        allocate(this%rst_types_3d(this%n_children))
        allocate(this%rst_types_2d(this%n_children))
        allocate(this%put_types(this%n_children))
        allocate(this%force_types(this%n_children))

        allocate(this%child_get_types_3d(this%n_children))
        allocate(this%child_get_types_2d(this%n_children))
        allocate(this%child_rst_types_3d(this%n_children))
        allocate(this%child_rst_types_2d(this%n_children))
        allocate(this%child_put_types(this%n_children))
        allocate(this%child_force_types(this%n_children))

        do i = 1,this%n_children
            ! +2 included to account for staggered grids for output variables
            call MPI_Type_create_subarray(4, [this%n_w_3d, (this%i_e_w-this%i_s_w+2), (this%k_e_w-this%k_s_w+1), (this%j_e_w-this%j_s_w+2)], &
                [this%n_w_3d, (this%iewc(i)-this%iswc(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jewc(i)-this%jswc(i)+2)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%get_types_3d(i))
        
            call MPI_Type_create_subarray(3, [this%n_w_2d, (this%i_e_w-this%i_s_w+2), (this%j_e_w-this%j_s_w+2)], &
                [this%n_w_2d, (this%iewc(i)-this%iswc(i)+2), (this%jewc(i)-this%jswc(i)+2)], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%get_types_2d(i))

            call MPI_Type_create_subarray(4, [this%n_w_3d, (this%i_e_re-this%i_s_re+2), (this%k_e_w-this%k_s_w+1), (this%j_e_re-this%j_s_re+2)], &
                [this%n_w_3d, (this%ierec(i)-this%isrec(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jerec(i)-this%jsrec(i)+2)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%rst_types_3d(i))
        
            call MPI_Type_create_subarray(3, [this%n_w_2d, (this%i_e_re-this%i_s_re+2), (this%j_e_re-this%j_s_re+2)], &
                [this%n_w_2d, (this%ierec(i)-this%isrec(i)+2), (this%jerec(i)-this%jsrec(i)+2)], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%rst_types_2d(i))

            call MPI_Type_create_subarray(4, [this%n_r, (this%i_e_r-this%i_s_r+1), (this%k_e_r-this%k_s_r+1), (this%j_e_r-this%j_s_r+1)], &
                [this%n_r, (this%ierc(i)-this%isrc(i)+1), (this%kerc(i)-this%ksrc(i)+1), (this%jerc(i)-this%jsrc(i)+1)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%put_types(i))

            call MPI_Type_create_subarray(4, [this%n_w_3d, this%nx_re, this%nz_w, this%ny_re], &
                [this%n_w_3d, (this%iewc(i)-this%iswc(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jewc(i)-this%jswc(i)+2)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_get_types_3d(i))

            call MPI_Type_create_subarray(3, [this%n_w_2d, this%nx_re, this%ny_re], &
                [this%n_w_2d, (this%iewc(i)-this%iswc(i)+2), (this%jewc(i)-this%jswc(i)+2)], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_get_types_2d(i))

            call MPI_Type_create_subarray(4, [this%n_w_3d, this%nx_re, this%nz_w, this%ny_re], &
                [this%n_w_3d, (this%ierec(i)-this%isrec(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jerec(i)-this%jsrec(i)+2)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_rst_types_3d(i))

            call MPI_Type_create_subarray(3, [this%n_w_2d, this%nx_re, this%ny_re], &
                [this%n_w_2d, (this%ierec(i)-this%isrec(i)+2), (this%jerec(i)-this%jsrec(i)+2)], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_rst_types_2d(i))

            call MPI_Type_create_subarray(4, [this%n_r, this%nx_r, this%nz_r, this%ny_r], &
                [this%n_r, (this%ierc(i)-this%isrc(i)+1), (this%kerc(i)-this%ksrc(i)+1), (this%jerc(i)-this%jsrc(i)+1)], &
                [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_put_types(i))

            if (this%n_f > 0) then
                call MPI_Type_create_subarray(4, [this%n_f, (this%i_e_w-this%i_s_w+2), (this%k_e_w-this%k_s_w+1), (this%j_e_w-this%j_s_w+2)], &
                    [this%n_f, (this%iewc(i)-this%iswc(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jewc(i)-this%jswc(i)+2)], &
                    [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%force_types(i))
                call MPI_Type_commit(this%force_types(i))

                call MPI_Type_create_subarray(4, [this%n_f, this%nx_w, this%nz_w, this%ny_w], &
                    [this%n_f, (this%iewc(i)-this%iswc(i)+2), (this%kewc(i)-this%kswc(i)+1), (this%jewc(i)-this%jswc(i)+2)], &
                    [0,0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, this%child_force_types(i))
                call MPI_Type_commit(this%child_force_types(i))

            endif
            call MPI_Type_commit(this%get_types_3d(i))
            call MPI_Type_commit(this%get_types_2d(i))
            call MPI_Type_commit(this%rst_types_3d(i))
            call MPI_Type_commit(this%rst_types_2d(i))
            call MPI_Type_commit(this%put_types(i))
            call MPI_Type_commit(this%child_get_types_3d(i))
            call MPI_Type_commit(this%child_get_types_2d(i))
            call MPI_Type_commit(this%child_rst_types_3d(i))
            call MPI_Type_commit(this%child_rst_types_2d(i))
            call MPI_Type_commit(this%child_put_types(i))
        enddo
    end subroutine setup_MPI_types

    subroutine init_with_clients(this)
        class(ioserver_t),   intent(inout)  :: this
        integer :: n, comm_size, ierr, i
        integer, allocatable, dimension(:) :: cnts, disps
                        
        type(MPI_Group) :: family_group

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

        call MPI_Comm_Group(this%client_comms,family_group)
        call MPI_Comm_rank(this%client_comms,n)

        call MPI_Group_Excl(family_group,1,[n],this%children_group)

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

        integer :: i, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w, msg_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
        msg_size = 1
        disp = 0

        this%parent_write_buffer_3d = kEMPT_BUFF
        this%parent_write_buffer_2d = kEMPT_BUFF

        ! Do MPI_Win_Start on write_win to initiate get
        call MPI_Win_Start(this%children_group,0,this%write_win_3d)
        call MPI_Win_Start(this%children_group,0,this%write_win_2d)

        ! Loop through child images and get chunks of buffer array from each one
        do i=1,this%n_children
            call MPI_Get(this%parent_write_buffer_3d(1,this%iswc(i),1,this%jswc(i)), msg_size, &
                this%get_types_3d(i), this%children_ranks(i), disp, msg_size, this%child_get_types_3d(i), this%write_win_3d)

            call MPI_Get(this%parent_write_buffer_2d(1,this%iswc(i),this%jswc(i)), msg_size, &
                this%get_types_2d(i), this%children_ranks(i), disp, msg_size, this%child_get_types_2d(i), this%write_win_2d)
        enddo
        ! Do MPI_Win_Complete on write_win to end get
        call MPI_Win_Complete(this%write_win_3d)
        call MPI_Win_Complete(this%write_win_2d)

        if (ALL(this%parent_write_buffer_3d==kEMPT_BUFF) .or. ALL(this%parent_write_buffer_2d==kEMPT_BUFF))then
            stop 'Error, all of write buffer used for output was still set to empty buffer flag at time of writing.'
        endif

        call this%outputer%save_out_file(this%sim_time,this%IO_Comms,this%out_var_indices,this%rst_var_indices)        
        
        call this%increment_output_time()

    end subroutine 

    
    ! This subroutine calls the read file function from the input object
    ! and then passes the read-in data to the read buffer
    module subroutine read_file(this)
        class(ioserver_t), intent(inout) :: this

        real, allocatable, dimension(:,:,:,:) :: parent_read_buffer

        !See if we even have files to read
        if (this%files_to_read) then
            ! read file into buffer array
            call this%reader%read_next_step(parent_read_buffer,this%IO_Comms)
            this%files_to_read = .not.(this%reader%eof)

            ! Loop through child images and send chunks of buffer array to each one
            call this%scatter_forcing(parent_read_buffer)
        endif

        ! increment input time whether we have files to read or not
        ! this is because input time is just controlling if we are at an input event -- not necessarily
        ! if we read in data or not
        ! the reader object will handle if we try to read a time step which doesn't exist in the forcing
        ! call this%increment_input_time()
    end subroutine 

    module subroutine gather_forcing(this)
        class(ioserver_t), intent(inout) :: this

        integer :: i, msg_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

        msg_size = 1
        disp = 0

        if (this%n_f <= 0) then
            if (STD_OUT_PE) write(*,*) 'No forcing fields to gather, but we entered gather_forcing. This is a bug.'
            return
        endif
        
        ! Do MPI_Win_Start on nest_win to initiate get
        call MPI_Win_Start(this%children_group,0,this%nest_win)

        ! Loop through child images and get chunks of buffer array from each one
        do i=1,this%n_children
            call MPI_Get(this%gather_buffer(1,this%iswc(i),1,this%jswc(i)), msg_size, &
                this%force_types(i), this%children_ranks(i), disp, msg_size, this%child_force_types(i), this%nest_win)
        enddo

        ! Do MPI_Win_Complete on read_win to end put
        call MPI_Win_Complete(this%nest_win)

    end subroutine gather_forcing

    module subroutine distribute_forcing(this, child_ioserver, child_indx)
        class(ioserver_t), intent(inout) :: this
        class(ioserver_t), intent(inout)    :: child_ioserver
        integer, intent(in) :: child_indx

        integer :: i, nx, ny, n, ierr, msg_size, real_size
        integer, allocatable :: send_msg_size_alltoall(:), buff_msg_size_alltoall(:), disp_alltoall(:)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lowerbound, extent
        real, allocatable, dimension(:,:,:,:) :: forcing_buffer

        allocate(send_msg_size_alltoall(this%n_servers))
        allocate(buff_msg_size_alltoall(this%n_servers))
        allocate(disp_alltoall(this%n_servers))

        disp_alltoall = 0
    
        allocate(forcing_buffer(this%n_f,child_ioserver%i_s_r:child_ioserver%i_e_r,child_ioserver%k_s_r:child_ioserver%k_e_r, child_ioserver%j_s_r:child_ioserver%j_e_r))

        ! If this is the first time calling gather_forcing, we are still in initialization. Call setup_nest_types now, passing in the child ioserver
        if (this%nest_types_initialized .eqv. .False.) then
            call this%setup_nest_types(child_ioserver, this%send_nest_types(child_indx,:), this%buffer_nest_types(child_indx,:))
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
                        forcing_buffer, buff_msg_size_alltoall, disp_alltoall, this%buffer_nest_types(child_indx,:), this%IO_Comms)
        ! This call will scatter the forcing fields to the ioclients of the nest child
        call child_ioserver%scatter_forcing(forcing_buffer)
        ! call child_ioserver%increment_input_time()
    end subroutine

    module subroutine scatter_forcing(this, forcing_buffer)
        class(ioserver_t), intent(inout) :: this
        real, dimension(:,:,:,:), intent(in) :: forcing_buffer

        integer :: i, msg_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp

        msg_size = 1
        disp = 0

        ! Do MPI_Win_Start on read_win to initiate put
        call MPI_Win_Start(this%children_group,0,this%read_win)

        ! Loop through child images and send chunks of buffer array to each one
        do i=1,this%n_children
            call MPI_Put(forcing_buffer(1,(this%isrc(i)-this%i_s_r+1),1,(this%jsrc(i)-this%j_s_r+1)), msg_size, &
                this%put_types(i), this%children_ranks(i), disp, msg_size, this%child_put_types(i), this%read_win)
        enddo

        ! Do MPI_Win_Complete on read_win to end put
        call MPI_Win_Complete(this%read_win)

    end subroutine

    ! This subroutine reads in the restart file and then sends the data to the children
    ! Same as above, but for restart file
    module subroutine read_restart_file(this, options)
        class(ioserver_t),   intent(inout) :: this
        type(options_t),     intent(in)    :: options

        integer :: i, n_3d, n_2d, nx, ny, i_s_re, i_e_re, j_s_re, j_e_re
        integer :: ncid, var_id, dimid_3d(4), nz, err, varid, start_3d(4), cnt_3d(4), start_2d(3), cnt_2d(3), msg_size
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
        real, allocatable :: data3d(:,:,:,:)
        type(variable_t)  :: var
        character(len=kMAX_NAME_LENGTH) :: name
        character(len=kMAX_FILE_LENGTH) :: base_rst_file_name, tmp_str, restart_in_file
        character(len=kMAX_STRING_LENGTH), allocatable :: tokens(:)

        integer :: restart_step                         ! time step relative to the start of the restart file
        type(Time_type) :: time_at_step   ! restart date as a modified julian day        

        msg_size = 1
        disp = 0

        tokens = split_str(trim(options%domain%init_conditions_file), "/")

        ! store the base file name temporarily here.
        tmp_str = trim(tokens(size(tokens)))

        ! Set the output filename to be the location of the output folder, and the name of the initial conditions file, with ".nc" removed
        base_rst_file_name = trim(options%restart%restart_folder) // '/' // tmp_str(1:(len_trim(tmp_str)-3)) // "_"

        write(restart_in_file, '(A,A,".nc")')    &
                trim(base_rst_file_name),   &
                trim(options%restart%restart_time%as_string('(I4,"-",I0.2,"-",I0.2,"_",I0.2,"-",I0.2,"-",I0.2)'))

        call check_file_exists(restart_in_file, message='Restart file does not exist.')

        ! find the time step that most closely matches the requested restart time (<=)
        restart_step = find_timestep_in_file(restart_in_file, 'time', options%restart%restart_time, time_at_step)

        if (options%general%debug) then
            write(*,*) " ------------------ "
            write(*,*) "RESTART INFORMATION"
            write(*,*) "mjd",         options%restart%restart_time%mjd()
            write(*,*) "date:",       trim(options%restart%restart_time%as_string())
            write(*,*) "file:",   trim(restart_in_file)
            write(*,*) "forcing step",options%restart%restart_step_in_file
            write(*,*) " ------------------ "
        endif

        err = nf90_open(restart_in_file, IOR(nf90_nowrite,NF90_NETCDF4), ncid, &
                comm = this%IO_Comms%MPI_VAL, info = MPI_INFO_NULL%MPI_VAL)
        
        ! setup start/count arrays accordingly. k_s_w and k_e_w forseen to always cover the bounds of what k_s_re and k_e_re would be
        ! This is because the domain is only decomposed in 2D.
        start_3d = (/ this%i_s_re,this%j_s_re,this%k_s_w,restart_step /)
        start_2d = (/ this%i_s_re,this%j_s_re,restart_step /)
        cnt_3d = (/ (this%i_e_re-this%i_s_re+1),(this%j_e_re-this%j_s_re+1),(this%k_e_w-this%k_s_w+1),1 /)
        cnt_2d = (/ (this%i_e_re-this%i_s_re+1),(this%j_e_re-this%j_s_re+1),1 /)

        this%parent_write_buffer_3d = kEMPT_BUFF
        this%parent_write_buffer_2d = kEMPT_BUFF

        n_2d = 1
        n_3d = 1

        do i = 1,size(this%rst_var_indices)
            var = this%outputer%variables(this%rst_var_indices(i))
            name = var%name
            call check_ncdf( nf90_inq_varid(ncid, name, var_id), " Getting var ID for "//trim(name))
            call check_ncdf( nf90_var_par_access(ncid, var_id, nf90_collective))
            
            nx = cnt_3d(1) + var%xstag
            ny = cnt_3d(2) + var%ystag

            if (var%three_d) then
                ! Get length of z dim
                call check_ncdf( nf90_inquire_variable(ncid, var_id, dimids = dimid_3d), " Getting dim IDs for "//trim(name))
                call check_ncdf( nf90_inquire_dimension(ncid, dimid_3d(3), len = nz), " Getting z dim len for "//trim(name))
                
                if (allocated(data3d)) deallocate(data3d)
                allocate(data3d(nx,ny,nz,1))
                call check_ncdf( nf90_get_var(ncid, var_id, data3d, start=start_3d, count=(/ nx, ny, nz /)), " Getting 3D var "//trim(name))

                this%parent_write_buffer_3d(n_3d,this%i_s_re:this%i_e_re+var%xstag,1:nz,this%j_s_re:this%j_e_re+var%ystag) = &
                        reshape(data3d(:,:,:,1), shape=[nx,nz,ny], order=[1,3,2])
                n_3d = n_3d+1
            else if (var%two_d) then
                call check_ncdf( nf90_get_var(ncid, var_id, this%parent_write_buffer_2d(n_2d,this%i_s_re:this%i_e_re+var%xstag,this%j_s_re:this%j_e_re+var%ystag), &
                        start=start_2d, count=(/ nx, ny /)), " Getting 2D "//trim(name))
                        n_2d = n_2d+1
            endif
        end do
        
        call check_ncdf(nf90_close(ncid), "Closing file "//trim(restart_in_file))
        
        ! Because this is for reading restart data, performance is not critical, and 
        ! we use a simple MPI_fence syncronization
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)

        ! Loop through child images and send chunks of buffer array to each one
        do i=1,this%n_children
            !Note that the MPI datatypes here are reversed since we are working with the write window
            call MPI_Put(this%parent_write_buffer_3d(1,this%isrec(i),1,this%jsrec(i)), msg_size, &
                this%rst_types_3d(i), this%children_ranks(i), disp, msg_size, this%child_rst_types_3d(i), this%write_win_3d)

            call MPI_Put(this%parent_write_buffer_2d(1,this%isrec(i),this%jsrec(i)), msg_size, &
                this%rst_types_2d(i), this%children_ranks(i), disp, msg_size, this%child_rst_types_2d(i), this%write_win_2d)
        enddo
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)

    end subroutine 


    ! This function closes all open file handles. Files are left open by default
    ! to minimize I/O calls. When the program exits, this must be called
    module subroutine close_files(this)
        class(ioserver_t), intent(inout) :: this
                
        ! close files
        call this%reader%close_file()
        call this%outputer%close_files()

    end subroutine 
    
end submodule
