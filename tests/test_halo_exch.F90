program test_halo_exch

    use variable_dict_interface, only: var_dict_t
    use variable_interface,      only: variable_t
    use mpi_f08
    use grid_interface, only: grid_t
    use halo_interface, only: halo_t
    use icar_constants

    implicit none
    
    type(halo_t) :: halo
    type(var_dict_t) :: exch_vars, adv_vars
    type(variable_t) :: qv, u, v, temperature
    integer :: my_index
    integer :: ierr
    type(grid_t) :: grid, grid_u, grid_v
    type(MPI_Comm) :: comms
    logical :: init_flag

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
    my_index = my_index + 1

    !set the comm MPI_com type to be MPI_COMM_WORLD
    CALL MPI_Comm_dup( MPI_COMM_WORLD, comms, ierr )

    !initialize grids
    call grid%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3)

    call grid_u%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3, nx_extra = 1)
    call grid_v%set_grid_dimensions( 20, 20, 20, image=my_index, comms=comms, adv_order=3, ny_extra = 1)

    !initialize variables to exchange
    call qv%initialize(grid)
    call temperature%initialize(grid)

    !make a staggered x and y variables to test staggered exchange
    call u%initialize(grid_u)
    call v%initialize(grid_v)

    !populate adv_vars with two test variables
    call adv_vars%add_var('qv', qv) 
    !call adv_vars%add_var('temperature', temperature) 

    call halo%init(exch_vars, adv_vars, grid, comms)
  
    ! Initialize fields with my_index values
    qv%data_3d = my_index
    temperature%data_3d = my_index
    u%data_3d = my_index
    v%data_3d = my_index

    call halo%batch_exch(exch_vars, adv_vars)
    call halo%exch_var(temperature)
    call halo%exch_var(u)
    call halo%exch_var(v)
    
    
    ! Verify exchange worked
    ! check if this image is not on the eastern boundary
    if (.not.(grid%ximg == grid%ximages)) then
        !check that the eastern halo is filled with the value of my_index for the eastern neighbor
        if (.not.(ALL(qv%data_3d((grid%ite+1):grid%ime,:,grid%jts:grid%jte) == my_index+1))) then
            write(*,*) "Failed for eastern halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d((grid_u%ite+1):grid_u%ime,:,grid_u%jts:grid_u%jte) == my_index+1))) then
            write(*,*) "Failed for eastern halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d((grid_v%ite+1):grid_v%ime,:,grid_v%jts:grid_v%jte) == my_index+1))) then
            write(*,*) "Failed for eastern halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! check if this image is not on the western boundary
    if (.not.(grid%ximg == 1)) then
        !check that the western halo is filled with the value of my_index for the western neighbor
        if (.not.(ALL(qv%data_3d(grid%ims:grid%its-1,:,grid%jts:grid%jte) == my_index-1))) then
            write(*,*) "Failed for western halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%ims:grid_u%its-1,:,grid_u%jts:grid_u%jte) == my_index-1))) then
            write(*,*) "Failed for western halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%ims:grid_v%its-1,:,grid_v%jts:grid_v%jte) == my_index-1))) then
            write(*,*) "Failed for western halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! check if this image is not on the southern boundary
    if (.not.(grid%yimg == 1)) then
        !check that the southern halo is filled with the value of my_index for the southern neighbor
        if (.not.(ALL(qv%data_3d(grid%its:grid%ite,:,grid%jms:grid%jts-1) == my_index-grid%ximages))) then
            write(*,*) "Failed for southern halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%its:grid_u%ite,:,grid_u%jms:grid_u%jts-1) == my_index-grid%ximages))) then
            write(*,*) "Failed for southern halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%its:grid_v%ite,:,grid_v%jms:grid_v%jts-1) == my_index-grid%ximages))) then
            write(*,*) "Failed for southern halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! check if this image is not on the northern boundary
    if (.not.(grid%yimg == grid%yimages)) then
        !check that the northern halo is filled with the value of my_index for the northern neighbor
        if (.not.(ALL(qv%data_3d(grid%its:grid%ite,:,grid%jte+1:grid%jme) == my_index+grid%ximages))) then
            write(*,*) "Failed for northern halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%its:grid_u%ite,:,grid_u%jte+1:grid_u%jme) == my_index+grid%ximages))) then
            write(*,*) "Failed for northern halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%its:grid_v%ite,:,grid_v%jte+1:grid_v%jme) == my_index+grid%ximages))) then
            write(*,*) "Failed for northern halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif

    call halo%exch_var(u, corners=.true.)
    call halo%exch_var(v, corners=.true.)
    call halo%exch_var(temperature, corners=.true.)

    ! if this image is in the north eastern corner
    if ((grid%yimg == 1 .and. grid%ximg == 1)) then
        !check that the north eastern corner halo is filled with the value of my_index for the north eastern neighbor
        if (.not.(ALL(temperature%data_3d(grid%ite+1:grid%ime,:,grid%jte+1:grid%jme) == my_index+grid%ximages+1))) then
            write(*,*) "Failed for north eastern corner halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%ite+1:grid_u%ime,:,grid_u%jte+1:grid_u%jme) == my_index+grid%ximages+1))) then
            write(*,*) "Failed for north eastern corner halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%ite+1:grid_v%ime,:,grid_v%jte+1:grid_v%jme) == my_index+grid%ximages+1))) then
            write(*,*) "Failed for north eastern corner halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! if this image is in the north western corner
    if ((grid%yimg == 1 .and. grid%ximg == grid%ximages)) then
        !check that the north western corner halo is filled with the value of my_index for the north western neighbor
        if (.not.(ALL(temperature%data_3d(grid%ims:grid%its-1,:,grid%jte+1:grid%jme) == my_index+grid%ximages-1))) then
            write(*,*) "Failed for north western corner halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%ims:grid_u%its-1,:,grid_u%jte+1:grid_u%jme) == my_index+grid%ximages-1))) then
            write(*,*) "Failed for north western corner halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%ims:grid_v%its-1,:,grid_v%jte+1:grid_v%jme) == my_index+grid%ximages-1))) then
            write(*,*) "Failed for north western corner halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! if this image is in the south eastern corner
    if ((grid%yimg == grid%yimages .and. grid%ximg == 1)) then
        !check that the south eastern corner halo is filled with the value of my_index for the south eastern neighbor
        if (.not.(ALL(temperature%data_3d(grid%ite+1:grid%ime,:,grid%jms:grid%jts-1) == my_index-grid%ximages+1))) then
            write(*,*) "Failed for south eastern corner halo exchange on image ", my_index
            print *, temperature%data_3d(grid%ite+1:grid%ime,1,grid%jms:grid%jts-1)
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%ite+1:grid_u%ime,:,grid_u%jms:grid_u%jts-1) == my_index-grid%ximages+1))) then
            write(*,*) "Failed for south eastern corner halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%ite+1:grid_v%ime,:,grid_v%jms:grid_v%jts-1) == my_index-grid%ximages+1))) then
            write(*,*) "Failed for south eastern corner halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! if this image is in the south western corner
    if ((grid%yimg == grid%yimages .and. grid%ximg == grid%ximages)) then
        !check that the south western corner halo is filled with the value of my_index for the south western neighbor
        if (.not.(ALL(temperature%data_3d(grid%ims:grid%its-1,:,grid%jms:grid%jts-1) == my_index-grid%ximages-1))) then
            write(*,*) "Failed for south western corner halo exchange on image ", my_index
            stop
        endif
        if (.not.(ALL(u%data_3d(grid_u%ims:grid_u%its-1,:,grid_u%jms:grid_u%jts-1) == my_index-grid%ximages-1))) then
            write(*,*) "Failed for south western corner halo exchange, u_grid, on image ", my_index
            stop
        endif
        if (.not.(ALL(v%data_3d(grid_v%ims:grid_v%its-1,:,grid_v%jms:grid_v%jts-1) == my_index-grid%ximages-1))) then
            write(*,*) "Failed for south western corner halo exchange, v_grid, on image ", my_index
            stop
        endif
    endif
    call MPI_Finalize(ierr)
    ! Print message indicating successful completion
    if (my_index == 1) write(*,*) "Halo exchange completed successfully."
end program test_halo_exch
