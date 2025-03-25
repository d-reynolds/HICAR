
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
  use debug_module,             only : check_ncdf
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use output_metadata,          only : get_varindx


  implicit none

contains        
    
    module subroutine init(this, domain, forcing, options, n_indx)
        implicit none
        class(ioclient_t),  intent(inout)  :: this
        type(domain_t),     intent(inout)  :: domain
        type(boundary_t),   intent(in)     :: forcing
        type(options_t),    intent(in)     :: options(:)
        integer,            intent(in)     :: n_indx

        type(variable_t) :: var
        type(MPI_Group) :: family_group

        integer :: my_rank, comm_size, some_child_id
        
        if (STD_OUT_PE) write(*,*) 'Initializing I/O Clients'
        if (STD_OUT_PE) flush(output_unit)
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
        if (domain%ime == domain%ide) this%i_e_w = domain%ide !Add extra to accomodate staggered vars
        if (domain%jms == domain%jds) this%j_s_w = domain%jds
        if (domain%jme == domain%jde) this%j_e_w = domain%jde !Add extra to accomodate staggered vars


        ! If this run uses nests, then set the variables which we need to output to a child nest
        ! Because microphysics must be the same for all nests, we can use any of the forcing_options
        ! structures on a child nest, since these are all the same. Use domain #2, since there will always
        ! be at least 1 nest following this condition
        if (size(options(n_indx)%general%child_nests) > 0) then
            some_child_id = options(n_indx)%general%child_nests(1)
            this%vars_for_nest = options(some_child_id)%forcing%vars_to_read
        endif

        call init_with_server(this)

        call setup_MPI_windows(this)

        !Setup the parent-child group used for buffer communication
        call MPI_Comm_Group(this%parent_comms,family_group)
        call MPI_Comm_size(this%parent_comms, comm_size)

        !Our rank is the last process in the group, which is equal to n_children
        call MPI_Group_Incl(family_group, 1, [comm_size-1], this%parent_group)

        ! Do MPI_Win_Post on read_buffer to indicate that we are open for delivery of input data
        call MPI_Win_Post(this%parent_group,0,this%read_win)

    end subroutine init

    module subroutine init_with_server(this)
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
        integer :: nx_w, nz_w, ny_w, nx_re, ny_re, n_w_2d, n_f, n_w_3d, ierr
        integer :: nx_r, nz_r, ny_r, n_r, real_size

        CALL MPI_Type_size(MPI_REAL, real_size)

        nx_w = 0
        nz_w = 0
        ny_w = 0
        n_w_3d = 0
        n_w_2d = 0
        n_f = 0
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

        call MPI_Allreduce(MPI_IN_PLACE,n_r,1,MPI_INT,MPI_MAX,this%parent_comms,ierr)

        ! +1 added to handle variables on staggered grids
        win_size = n_w_3d*nx_re*nz_w*ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%write_win_3d)
        call C_F_POINTER(tmp_ptr, this%write_buffer_3d, [n_w_3d, nx_re, nz_w, ny_re])
        this%write_buffer_3d = kEMPT_BUFF

        ! +1 added to handle variables on staggered grids
        win_size = n_w_2d*nx_re*ny_re
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%write_win_2d)
        call C_F_POINTER(tmp_ptr, this%write_buffer_2d, [n_w_2d, nx_re, ny_re])
        this%write_buffer_2d = kEMPT_BUFF

        ! This is the buffer for the forcing data. Necesarry so that output data can sit around in its own buffer while the ioserver
        ! is busy distributing the forcing data to the child nests 
        if (n_f > 0) then
            win_size = n_f*nx_w*ny_w*nz_w
            call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%forcing_win)
            call C_F_POINTER(tmp_ptr, this%forcing_buffer, [n_f, nx_w, nz_w, ny_w])
            this%forcing_buffer = kEMPT_BUFF
        endif

        win_size = n_r*nx_r*nz_r*ny_r
        call MPI_WIN_ALLOCATE_SHARED(win_size*real_size, real_size, MPI_INFO_NULL, this%parent_comms, tmp_ptr, this%read_win)
        call C_F_POINTER(tmp_ptr, this%read_buffer, [n_r, nx_r, nz_r, ny_r])
        this%read_buffer = kEMPT_BUFF

    end subroutine setup_MPI_windows


    ! This subroutine pushes the output fields from the domain object
    ! to the write buffer for the IO processes to use
    module subroutine push(this, domain)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),   intent(inout)    :: domain
        
        type(variable_t) :: var
        integer :: i, n_3d, n_2d, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w
                        
        n_3d = 1
        n_2d = 1

        !This is false only when it is the first call to push (i.e. first write call)
        if (this%written) then
            ! Do MPI_Win_Wait on write_buffer to make sure that server process has completed writing of previous data
            call smart_wait(this%write_win_3d, 'Error: Timeout waiting for write_win_3d completion')
            call smart_wait(this%write_win_2d, 'Error: Timeout waiting for write_win_2d completion')
        endif
        this%written = .False.

        ! loop through the list of variables that need to be written out
        call domain%vars_to_out%reset_iterator()

        do while (domain%vars_to_out%has_more_elements())
            ! get the next variable in the structure
            var = domain%vars_to_out%next()
            
            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            if (domain%ime == domain%ide) i_e_w = i_e_w+var%xstag !Add extra to accomodate staggered vars
            if (domain%jme == domain%jde) j_e_w = j_e_w+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1
            if (var%two_d) then
                if (var%dtype == kREAL) then
                    this%write_buffer_2d(n_2d,1:nx,1:ny) = &
                        var%data_2d(i_s_w:i_e_w,j_s_w:j_e_w)
                elseif (var%dtype == kDOUBLE) then
                    this%write_buffer_2d(n_2d,1:nx,1:ny) = &
                        real(var%data_2dd(i_s_w:i_e_w,j_s_w:j_e_w))
                endif
                n_2d = n_2d+1

            else
                this%write_buffer_3d(n_3d,1:nx,1:var%dim_len(2),1:ny) = &
                        var%data_3d(i_s_w:i_e_w,1:var%dim_len(2),j_s_w:j_e_w)
                n_3d = n_3d+1
            endif
        enddo
        
        this%written = .True.
        call domain%increment_output_time()

        ! Do MPI_Win_Post on write_win to inform that server process can begin writing of data
        call MPI_Win_Post(this%parent_group,0,this%write_win_3d)
        call MPI_Win_Post(this%parent_group,0,this%write_win_2d)


    end subroutine 

    module subroutine update_nest(this, domain)
        implicit none
        class(ioclient_t),   intent(inout) :: this
        type(domain_t),      intent(in)    :: domain

        type(variable_t) :: var
        integer :: i, n_3d, nx, ny, i_s_w, i_e_w, j_s_w, j_e_w, err

        n_3d = 1

        ! Do MPI_Win_Wait on forcing_win. This is mostly unnecesarry, since the server process is what will be waiting on us, which is handeled by the MPI_Win_Start call
        ! in ioserver%gather_foring. Still, MPI_Win_Wait is called here for completeness of PSCW model
        if (this%nest_updated) then
            call smart_wait(this%forcing_win, 'Error: Timeout waiting for forcing_win completion')
        endif

        this%nest_updated = .False.
        !This is false only when it is the first call to push (i.e. first write call)

        do i = 1, size(this%vars_for_nest)
            ! get the next variable in the structure
            if (this%vars_for_nest(i) == '') cycle
            var = domain%vars_to_out%get_var(this%vars_for_nest(i),err=err)
            if (var%two_d) cycle

            if (err > 0) then
                write(*,*) 'Error: Variable ', this%vars_for_nest(i), ' not found in parent domain: ', domain%nest_indx
                stop
            end if
            i_s_w = this%i_s_w; i_e_w = this%i_e_w
            j_s_w = this%j_s_w; j_e_w = this%j_e_w
            if (domain%ime == domain%ide) i_e_w = i_e_w+var%xstag !Add extra to accomodate staggered vars
            if (domain%jme == domain%jde) j_e_w = j_e_w+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_w - i_s_w + 1
            ny = j_e_w - j_s_w + 1

            this%forcing_buffer(n_3d,1:nx,1:var%dim_len(2),1:ny) = &
                var%data_3d(i_s_w:i_e_w,1:var%dim_len(2),j_s_w:j_e_w)
            n_3d = n_3d+1
        enddo

        this%nest_updated = .True.
        ! Do MPI_Win_Post on forcing_win to inform that server process can begin gathering of nest data
        call MPI_Win_Post(this%parent_group,0,this%forcing_win)

    end subroutine
    
    ! This subroutine receives the input fields from the IO buffer
    ! for assignment to the forcing object
    module subroutine receive(this, forcing, domain)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(boundary_t), intent(inout)  :: forcing
        type(domain_t),   intent(inout)  :: domain

        type(variable_t)     :: var
        integer :: i, n, nx, ny
                
        n = 1
        nx = this%i_e_r - this%i_s_r + 1
        ny = this%j_e_r - this%j_s_r + 1

        ! Do MPI_Win_Wait on read_buffer to make sure that server process has completed data transfer
        call smart_wait(this%read_win, 'Error: Timeout waiting for read_win completion')

        ! loop through the list of variables that need to be read in
        call forcing%variables%reset_iterator()
        
        !If the parent I/O server has not yet written all of our input vars, wait

        do while (forcing%variables%has_more_elements())
            ! get the next variable in the structure
            var = forcing%variables%next()
            if (var%computed) then
                cycle
            else
                if (var%two_d) then
                    if (var%dtype == kREAL) then
                        var%data_2d = this%read_buffer(n,1:nx,1,1:ny)
                    elseif (var%dtype == kDOUBLE) then
                        var%data_2dd = dble(this%read_buffer(n,1:nx,1,1:ny))
                    endif
                else
                    var%data_3d = this%read_buffer(n,1:nx,1:var%dim_len(2),1:ny)
                endif
                n = n+1
            endif
        enddo
    
        ! Do MPI_Win_Post on read_buffer to indicate that we are open for delivery of new input data
        call MPI_Win_Post(this%parent_group,0,this%read_win)
        call domain%increment_input_time()


    end subroutine 
    
    ! This subroutine assigns the data from the write_buffer to the appropriate fields
    ! of the domain object. It is a "reverse write"
    module subroutine receive_rst(this, domain, options)
        implicit none
        class(ioclient_t), intent(inout) :: this
        type(domain_t),   intent(inout)  :: domain
        type(options_t),  intent(in)     :: options

        type(variable_t)     :: var
        integer :: i, n_2d, n_3d, nx, ny, i_s_re, i_e_re, j_s_re, j_e_re

        ! Because this is for reading restart data, performance is not critical, and 
        ! we use a simple MPI_fence syncronization
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)
        call MPI_Win_fence(0,this%write_win_3d)
        call MPI_Win_fence(0,this%write_win_2d)

        n_3d = 1
        n_2d = 1

        ! loop through the list of variables that need to be written out
        call domain%vars_to_out%reset_iterator()
        
        do while (domain%vars_to_out%has_more_elements())
            ! get the next variable in the structure
            var = domain%vars_to_out%next()

            !See if var is in restart vars
            if (options%vars_for_restart(get_varindx(var%name)) <= 0) then
                cycle
            endif

            i_s_re = this%i_s_re; i_e_re = this%i_e_re
            j_s_re = this%j_s_re; j_e_re = this%j_e_re
            i_e_re = i_e_re+var%xstag !Add extra to accomodate staggered vars
            j_e_re = j_e_re+var%ystag !Add extra to accomodate staggered vars
            nx = i_e_re - i_s_re + 1
            ny = j_e_re - j_s_re + 1

            if (var%three_d) then
                var%data_3d(i_s_re:i_e_re,1:var%dim_len(2),j_s_re:j_e_re) = &
                    this%write_buffer_3d(n_3d,1:nx,1:var%dim_len(2),1:ny)
                n_3d = n_3d+1
            else
                if (var%dtype == kREAL) then
                    var%data_2d(i_s_re:i_e_re,j_s_re:j_e_re) = &
                        this%write_buffer_2d(n_2d,1:nx,1:ny)
                elseif (var%dtype == kDOUBLE) then
                    var%data_2dd(i_s_re:i_e_re,j_s_re:j_e_re) = &
                            dble(this%write_buffer_2d(n_2d,1:nx,1:ny))
                endif                
                n_2d = n_2d+1
            endif
        enddo

    end subroutine 


    !Necesary? Should communicate relevant post/wait calls and free window
    module subroutine close_client(this)
        implicit none
        class(ioclient_t), intent(inout) :: this

    end subroutine close_client

    subroutine smart_wait(window, err_msg)
        implicit none
        type(MPI_Win), intent(in) :: window
        character(len=*), intent(in) :: err_msg

        integer :: wait_count, ierr
        logical :: flag
        
        wait_count = 0
        flag = .False.

        call MPI_Win_Test(window, flag, ierr)
        do while (.not.(flag))
            
            ! Check if we've waited too long
            wait_count = wait_count + 1

            if (wait_count > 60.0) then
                write(*,*) err_msg
                stop
            endif
            
            ! Small sleep to avoid busy waiting
            call sleep(1)

            call MPI_Win_Test(window, flag, ierr)
        end do
    end subroutine

end submodule
