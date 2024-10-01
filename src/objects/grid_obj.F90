submodule(grid_interface) grid_implementation
    use assertions_mod, only : assert, assertions
    use iso_fortran_env
    implicit none

contains

    !> -------------------------------
    !! Return the dimensions of this grid as an n element array
    !!
    !! -------------------------------
    module function get_dims(this) result(dims)
        class(grid_t), intent(in) :: this
        integer, allocatable :: dims(:)

        if (this%is2d) then
            allocate(dims(2))
            dims(1) = this%ime - this%ims + 1
            dims(2) = this%jme - this%jms + 1
        endif
        if (this%is3d) then
            allocate(dims(3))
            dims(1) = this%ime - this%ims + 1
            dims(2) = this%kme - this%kms + 1
            dims(3) = this%jme - this%jms + 1
        endif
    end function
    

    !> -------------------------------
    !! Decompose the domain into as even a set of tiles as possible in two dimensions
    !!
    !! Searches through possible numbers of x and y tiles that multiple evenly to
    !! give the total number of images requested.
    !!
    !! For each x/y split compute the number of grid cells in both dimensions in each tile
    !! return the split that provides the closest match between the number of x and y grid cells
    !!
    !! -------------------------------
    module subroutine domain_decomposition(this, nx, ny, nimages, image, ratio)
        class(grid_t),  intent(inout) :: this
        integer,        intent(in)    :: nx, ny, nimages
        integer,        intent(in)    :: image
        real,           intent(in), optional :: ratio
        real    :: multiplier
        integer :: ysplit, xsplit, xs, ys, i
        real    :: best, current, x, y

        multiplier=1
        if (present(ratio)) multiplier = ratio

        xsplit = 1
        ysplit = nimages
        xs = xsplit
        ys = ysplit

        x = (nx/real(xsplit))
        y = (ny/real(ysplit))

        if (y > (multiplier*x)) then
            best = abs(1 - ( y / (multiplier*x) ))
        else
            best = abs(1 - ( (multiplier*x) / y ))
        endif

        do i=nimages,1,-1
            if (mod(nimages,i)==0) then
                ysplit = i
                xsplit = nimages / i

                x = (nx/float(xsplit))
                y = (ny/float(ysplit))

                if (y > (multiplier*x)) then
                    current = abs(1 - ( y / (multiplier*x) ))
                else
                    current = abs(1 - ( (multiplier*x) / y ))
                endif

                if (current < best) then
                    best = current
                    xs = xsplit
                    ys = ysplit
                endif
            endif
        enddo

        this%ximages = xs
        this%yimages = ys

        this%ximg = mod(image-1,  this%ximages)+1
        y = 100 * (real(image-1+epsilon) / real(this%ximages))
        this%yimg = (floor(y+1)/100)+1

        x = (nx/float(xs))
        y = (ny/float(ys))

        if (assertions) call assert((xs*ys) == nimages, "Number of tiles does not sum to number of images")
        ! if (image==1) print*, "ximgs=",xs, "yimgs=",ys

    end subroutine domain_decomposition

    !> -------------------------------
    !! Compute the number of grid cells in the current image along a dimension
    !!
    !! This takes care of the fact that generally the number of images will not evenly divide the number of grid cells
    !! In this case the extra grid cells need to be evenly distributed among all images
    !!
    !! n_global should be the full domain size of the dimension
    !! me should be this images image number along this dimension
    !! nimg should be the number of images this dimension will be divided into
    !!
    !! -------------------------------
    function my_n(n_global, me, nimg) result(n_local)
       integer, intent(in) :: n_global, me, nimg
       integer :: n_local

       ! add 1 if this image is less than the remainder that need an extra grid cell
       n_local = n_global / nimg + merge(1,0,me <= mod(n_global,nimg)  )
    end function

    !> -------------------------------
    !! Return the starting coordinate in the global domain coordinate system for a given image (me)
    !!
    !! -------------------------------
    function my_start(n_global, me, nimg) result(memory_start)
        implicit none
        integer, intent(in) :: n_global, me, nimg
        integer :: memory_start
        integer :: base_n

        base_n = n_global / nimg

        memory_start = (me-1)*(base_n) + min(me-1,mod(n_global,nimg)) + 1

    end function my_start
    
    !> -------------------------------
    !! Generate the domain decomposition mapping and compute the indicies for local memory
    !!
    !! -------------------------------
    module subroutine set_grid_dimensions(this, nx, ny, nz, image, comms, global_nz, adv_order, nx_extra, ny_extra)
      class(grid_t),   intent(inout) :: this
      integer,         intent(in)    :: nx, ny, nz, image
      type(MPI_Comm), optional, intent(in)    :: comms
      integer, optional, intent(in)    :: global_nz, adv_order, nx_extra, ny_extra

      integer :: halo_size, ierr

      halo_size = kDEFAULT_HALO_SIZE
      if (present(adv_order)) halo_size = ceiling(adv_order/2.0)

      this%nx_e = 0
      this%ny_e = 0
      if (present(nx_extra)) this%nx_e = nx_extra ! used to add 1 to the u-field staggered grid
      if (present(ny_extra)) this%ny_e = ny_extra ! used to add 1 to the v-field staggered grid

      call this%domain_decomposition(nx, ny, kNUM_COMPUTE, image=image)

      if (nz<1) then
          this%is2d = .True.
          this%is3d = .False.
          if (allocated(this%dimensions)) deallocate(this%dimensions)
          allocate(this%dimensions(2))
          this%dimensions(1) = "lat"
          this%dimensions(2) = "lon"
      else
          this%is2d = .False.
          this%is3d = .True.
          if (allocated(this%dimensions)) deallocate(this%dimensions)
          allocate(this%dimensions(3))
          this%dimensions(1) = "lat"
          this%dimensions(2) = "height"
          this%dimensions(3) = "lon"
      endif

      this%nz         = nz                                            ! note nz is both global and local
      this%nx         = my_n(nx, this%ximg, this%ximages) ! local grid size
      this%ny         = my_n(ny, this%yimg, this%yimages) ! local grid size

      if (this%nx <= halo_size .or. this%ny <= halo_size) then
          write(*,*) 'ERROR: tile size too small for halo size'
          write(*,*) 'ERROR: this usually results from an unfavorable domain decomposition'
          write(*,*) 'ERROR: number of images in x direction: ',this%ximages
          write(*,*) 'ERROR: number of images in y direction: ',this%yimages
          write(*,*) 'ERROR: try changing the number of processes used in combination with this domain'
          stop
        endif

      ! define the bounds needed for memory to store the data local to this image
      this%ims        = my_start(nx, this%ximg, this%ximages)
      this%ime        = this%ims + this%nx + this%nx_e - 1

      this%jms        = my_start(ny, this%yimg, this%yimages)
      this%jme        = this%jms + this%ny + this%ny_e - 1

      this%kms        = 1
      this%kme        = this%nz

      ! Now define the tile of data to process in physics routines
      this%kts = this%kms
      this%kte = this%kme

      ! The entire model domain begins at 1 and ends at nx,y,z
      this%ny_global  = ny + this%ny_e                                     ! global model domain grid size
      this%nx_global  = nx + this%nx_e                                     ! global model domain grid size

      this%ids = 1
      this%jds = 1
      this%kds = 1
      this%ide = this%nx_global
      this%jde = this%ny_global
      this%kde = this%nz

      this%halo_nz    = this%nz

      this%halo_size = halo_size
      call update_with_halos(this, halo_size)

      ! define the halo needed to manage communications between images
      ! perhaps this should be defined in exchangeable instead though?
      !this%ns_halo_nx = this%nx_global / this%ximages + 1 + this%nx_e  ! number of grid cells in x in the ns halo
      !this%ew_halo_ny = this%ny_global / this%yimages + 1 + this%ny_e  ! number of grid cells in y in the ew halo

      if (present(comms)) then
        if (.not.(comms==MPI_COMM_NULL)) then
            call MPI_Allreduce(this%nx,this%ns_halo_nx,1,MPI_INT,MPI_MAX,comms,ierr)
            !this%ns_halo_nx = this%ns_halo_nx !+ this%nx_e  ! number of grid cells in x in the ns halo

            call MPI_Allreduce(this%ny,this%ew_halo_ny,1,MPI_INT,MPI_MAX,comms,ierr)
            !this%ew_halo_ny = this%ew_halo_ny !+ this%ny_e  ! number of grid cells in y in the ew halo

            !If we have been passed the global_nz, it means that this grid is not the global, 3D grid. Thus, pass the global_nz to
            !create_MPI_types so that the window nz is set correctly in accordance with what is done in halo_obj.f90. If this 3D grid's nz
            !is larger than the global_nz, however, we have a problem, and should not try to make an MPI type for communication. There are
            !only a few 3D grids with nz's larger than 8, none of which should need to be exchanged, and the model should almost always be 
            !run with more than 8 z levels. So this should not be an issue.
            if (present(global_nz)) then
                if (this%nz <= global_nz) then
                call create_MPI_types(this, win_nz=global_nz)
                endif
            else
                call create_MPI_types(this)
            endif
        endif
      else
        this%ns_halo_nx = this%nx_global / this%ximages + 1 + this%nx_e  ! number of grid cells in x in the ns halo
        this%ew_halo_ny = this%ny_global / this%yimages + 1 + this%ny_e  ! number of grid cells in y in the ew halo
      endif
      !if (image >= 30 .and. ((this%nx_e+this%ny_e)==0) .and. this%nz==20) then
      !  write(*,*) 'image: ',image
      !  write(*,*) this%ns_halo_nx
      !  write(*,*) this%ew_halo_ny
      !  write(*,*) this%ximages
      !  write(*,*) this%yimages
      !  write(*,*) nx
      !  write(*,*) ny
      !  write(*,*) 'its: ',this%its
      !  write(*,*) 'ite: ',this%ite
      !  write(*,*) 'ims: ',this%ims
      !  write(*,*) 'ime: ',this%ime
      !  write(*,*) 'jts: ',this%jts
      !  write(*,*) 'jte: ',this%jte
      !  write(*,*) 'jms: ',this%jms
      !  write(*,*) 'jme: ',this%jme


      !endif

  end subroutine

  !> -------------------------------
  !! updates the grid memory dimensions with halo sizes if necessary
  !!
  !! -------------------------------
  subroutine update_with_halos(grid, halo_size)
      type(grid_t), intent(inout)   :: grid
      integer,      intent(in)      :: halo_size

      logical :: north_boundary, south_boundary, &
                 east_boundary,  west_boundary

      north_boundary = (grid%yimg == grid%yimages)
      south_boundary = (grid%yimg == 1)
      east_boundary  = (grid%ximg == grid%ximages)
      west_boundary  = (grid%ximg == 1)

      ! if this is on a given boundary, then add 0, if it is not a boundary than add/subtract halo_size
      !grid%ims = grid%its - merge(0, halo_size, west_boundary)
      !grid%ime = grid%ite + merge(0, halo_size, east_boundary)
      !grid%jms = grid%jts - merge(0, halo_size, south_boundary)
      !grid%jme = grid%jte + merge(0, halo_size, north_boundary)

      grid%ims = grid%ims - merge(0, halo_size, west_boundary)
      grid%ime = grid%ime + merge(0, halo_size, east_boundary)
      grid%jms = grid%jms - merge(0, halo_size, south_boundary)
      grid%jme = grid%jme + merge(0, halo_size, north_boundary)

      ! if this is on a boundary, we should skip 1 grid cell (the boundary conditions) else we should skip the halo
      grid%its = grid%ims + halo_size !merge(1, halo_size, west_boundary)
      grid%ite = grid%ime - halo_size !merge(1, halo_size, east_boundary)
      grid%jts = grid%jms + halo_size !merge(1, halo_size, south_boundary)
      grid%jte = grid%jme - halo_size !merge(1, halo_size, north_boundary)


      grid%nx = grid%ime - grid%ims + 1
      grid%ny = grid%jme - grid%jms + 1

  end subroutine

  subroutine create_MPI_types(grid, win_nz)
      type(grid_t), intent(inout)   :: grid
      integer, optional, intent(in) :: win_nz

      integer :: nz_win

      nz_win = grid%nz
      if (present(win_nz)) nz_win = win_nz

      if (grid%is3d) then
        call MPI_Type_create_subarray(3, [grid%nx, grid%nz, grid%ny], [(grid%ite-grid%its+1), grid%nz, grid%halo_size+grid%ny_e], &
                [grid%halo_size,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%NS_halo)
        call MPI_Type_create_subarray(3, [grid%nx, grid%nz, grid%ny], [grid%halo_size+grid%nx_e, grid%nz, (grid%jte-grid%jts+1)], &
                [0,0,grid%halo_size], MPI_ORDER_FORTRAN, MPI_REAL, grid%EW_halo)
        call MPI_Type_create_subarray(3, [grid%nx, grid%nz, grid%ny], [grid%halo_size, grid%nz, grid%halo_size], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_halo)


        call MPI_Type_create_subarray(3, [grid%ns_halo_nx-grid%nx_e+2*grid%halo_size, nz_win, grid%halo_size+1], [(grid%ite-grid%its+1), grid%nz, grid%halo_size+grid%ny_e], &
                [grid%halo_size,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%NS_win_halo)
        call MPI_Type_create_subarray(3, [grid%halo_size+1, nz_win, grid%ew_halo_ny-grid%ny_e+2*grid%halo_size], [grid%halo_size+grid%nx_e, grid%nz, (grid%jte-grid%jts+1)], &
                [0,0,grid%halo_size], MPI_ORDER_FORTRAN, MPI_REAL, grid%EW_win_halo)
        call MPI_Type_create_subarray(3, [grid%ns_halo_nx-grid%nx_e+2*grid%halo_size, nz_win, grid%halo_size+1], [grid%halo_size, grid%nz, grid%halo_size], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_NS_win_halo)
        call MPI_Type_create_subarray(3, [grid%halo_size+1, nz_win, grid%ew_halo_ny-grid%ny_e+2*grid%halo_size], [grid%halo_size, grid%nz, grid%halo_size], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_EW_win_halo)

      else
        call MPI_Type_create_subarray(2, [grid%nx, grid%ny], [(grid%ite-grid%its+1), grid%halo_size], &
                [grid%halo_size,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%NS_halo)
        call MPI_Type_create_subarray(2, [grid%nx, grid%ny], [grid%halo_size, (grid%jte-grid%jts+1)], &
                [0,grid%halo_size], MPI_ORDER_FORTRAN, MPI_REAL, grid%EW_halo)
        call MPI_Type_create_subarray(2, [grid%nx, grid%ny], [grid%halo_size, grid%halo_size], &
                [0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_halo)

        call MPI_Type_create_subarray(3, [grid%ns_halo_nx-grid%nx_e+2*grid%halo_size, nz_win, grid%halo_size+1], [(grid%ite-grid%its+1), 1, grid%halo_size+grid%ny_e], &
                [grid%halo_size,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%NS_win_halo)
        call MPI_Type_create_subarray(3, [grid%halo_size+1, nz_win, grid%ew_halo_ny-grid%ny_e+2*grid%halo_size], [grid%halo_size+grid%nx_e, 1, (grid%jte-grid%jts+1)], &
                [0,0,grid%halo_size], MPI_ORDER_FORTRAN, MPI_REAL, grid%EW_win_halo)
        call MPI_Type_create_subarray(3, [grid%ns_halo_nx-grid%nx_e+2*grid%halo_size, nz_win, grid%halo_size+1], [grid%halo_size, 1, grid%halo_size], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_NS_win_halo)
        call MPI_Type_create_subarray(3, [grid%halo_size+1, nz_win, grid%ew_halo_ny-grid%ny_e+2*grid%halo_size], [grid%halo_size, 1, grid%halo_size], &
                [0,0,0], MPI_ORDER_FORTRAN, MPI_REAL, grid%corner_EW_win_halo)

      endif
      call MPI_Type_commit(grid%NS_halo)
      call MPI_Type_commit(grid%NS_win_halo)
      call MPI_Type_commit(grid%EW_halo)
      call MPI_Type_commit(grid%EW_win_halo)
      call MPI_Type_commit(grid%corner_halo)
      call MPI_Type_commit(grid%corner_NS_win_halo)
      call MPI_Type_commit(grid%corner_EW_win_halo)

  end subroutine

end submodule
