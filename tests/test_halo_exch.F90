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
    type(variable_t) :: qv, temperature
    integer :: my_index
    integer :: ierr
    type(grid_t) :: grid
    type(MPI_Comm) :: comm
    logical :: init_flag

    write(*,*) "here"

    !Initialize MPI if needed
    init_flag = .False.
    call MPI_initialized(init_flag, ierr)
    if (.not.(init_flag)) then
        call MPI_INIT(ierr)
        init_flag = .True.
    endif
    write(*,*) "here"

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_index,ierr)
    my_index = my_index + 1

    write(*,*) "my_index = ", my_index
    !set the comm MPI_com type to be MPI_COMM_WORLD
    CALL MPI_Comm_dup( MPI_COMM_WORLD, comm, ierr )
    write(*,*) "my_index = ", my_index

    !initialize grid
    call grid%set_grid_dimensions( 100, 100, 20, image=my_index, comms=comm, adv_order=1)
    write(*,*) "my_index = ", my_index

    !initialize variables to exchange
    call qv%initialize(grid)
    call temperature%initialize(grid)
    write(*,*) "my_index = ", my_index

    !populate adv_vars with two test variables
    call adv_vars%add_var('qv', qv) 
    call adv_vars%add_var('temperature', temperature) 
    write(*,*) "my_index = ", my_index

    call halo%init(exch_vars, adv_vars, grid, comm)

    write(*,*) "my_index = ", my_index

    ! Initialize fields with my_index values
    qv%data_3d = my_index

    call halo%batch_exch(exch_vars, adv_vars)
    write(*,*) "my_index = ", my_index

    ! Verify exchange worked
    ! check if this image is not on the eastern boundary
    if (.not.(grid%ximg == grid%ximages)) then
        !check that the eastern halo is filled with the value of my_index for the eastern neighbor
        if (.not.(ALL(qv%data_3d((grid%ime-grid%halo_size+1):grid%ime,:,:) == my_index+1))) then
            print*, "Field qv after exchange:"
            print*, qv%data_3d
            stop
        endif
    endif
    write(*,*) "my_index = ", my_index

    ! check if this image is not on the western boundary
    if (.not.(grid%ximg == 1)) then
        !check that the western halo is filled with the value of my_index for the western neighbor
        if (.not.(ALL(qv%data_3d(grid%ims:(grid%ims+grid%halo_size-1),:,:) == my_index-1))) then
            print*, "Field qv after exchange:"
            print*, qv%data_3d
            stop
        endif
    endif
    write(*,*) "my_index = ", my_index

    ! check if this image is not on the southern boundary
    if (.not.(grid%yimg == 1)) then
        !check that the southern halo is filled with the value of my_index for the southern neighbor
        if (.not.(ALL(qv%data_3d(:,:,grid%jms:(grid%jms+grid%halo_size-1)) == my_index-grid%ximages))) then
            print*, "Field qv after exchange:"
            print*, qv%data_3d
            stop
        endif
    endif
    write(*,*) "my_index = ", my_index

    ! check if this image is not on the northern boundary
    if (.not.(grid%yimg == grid%yimages)) then
        !check that the northern halo is filled with the value of my_index for the northern neighbor
        if (.not.(ALL(qv%data_3d(:,:,(grid%jme-grid%halo_size+1):grid%jme) == my_index+grid%ximages))) then
            print*, "Field qv after exchange:"
            print*, qv%data_3d
            stop
        endif
    endif

    call MPI_Finalize(ierr)
    ! Print message indicating successful completion
    print*, "Halo exchange completed successfully."
end program test_halo_exch
