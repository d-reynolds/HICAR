!>------------------------------------------------------------
!! Module to solve for a 3D wind field following mass-conservation
!! using NVIDIA AmgX GPU-native solver library.
!!
!! This is a replacement for wind_iterative_old.F90 that uses AmgX
!! instead of PETSc, providing better GPU memory management and
!! native CUDA integration.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!  GitHub Copilot
!!
!!------------------------------------------------------------

module wind_iterative_amgx
    use iso_c_binding
    use domain_interface,  only : domain_t
    use icar_constants,    only : STD_OUT_PE, kVARS
    use options_interface, only : options_t
    use iso_fortran_env
    use mpi_f08
    use openacc

    implicit none
    private
    public:: init_iter_winds_amgx, calc_iter_winds_amgx, finalize_amgx

    ! AmgX C API types
    type(c_ptr) :: amgx_config = c_null_ptr
    type(c_ptr) :: amgx_resources = c_null_ptr
    type(c_ptr) :: amgx_matrix = c_null_ptr
    type(c_ptr) :: amgx_vector_solution = c_null_ptr
    type(c_ptr) :: amgx_vector_rhs = c_null_ptr
    type(c_ptr) :: amgx_solver = c_null_ptr
    
    logical :: initialized = .False.
    logical :: first_solve = .True.
    ! Grid parameters (same as PETSc version)
    real, parameter::deg2rad=0.017453293
    real, parameter :: rad2deg=57.2957779371
    
    ! Stencil coefficient arrays
    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    real, allocatable, dimension(:,:)    :: dzdx_surf, dzdy_surf
    integer              :: hs, i_s, i_e, k_s, k_e, j_s, j_e
    integer              :: ims, ime, jms, jme, ids, ide, jds, jde
    integer              :: xs, ys, zs, xm, ym, zm, mx, my, mz  ! PETSc-style with boundary conditions layer
    
    ! Sparse matrix storage (CSR format for AmgX)
    integer(c_int), allocatable, target :: row_ptrs(:)      ! Row pointers (n+1 elements)
    integer(c_int), allocatable, target :: col_indices(:)   ! Column indices (GLOBAL indices for distributed mode)
    real(c_double), allocatable, target :: values(:)        ! Matrix values
    real(c_double), allocatable, target :: rhs(:)           ! Right-hand side
    real(c_double), allocatable, target :: solution(:)      ! Solution vector
    integer(c_int), allocatable, target :: partition_vec(:) ! Partition vector for distributed mode (size n+1)
    integer(c_int) :: n_rows, nnz                   ! Matrix dimensions (local)
    integer(c_int) :: n_rows_global                 ! Global number of rows across all ranks
    integer :: local_offset, global_n               ! MPI partitioning
    
    ! AmgX C interface declarations
    interface
        ! Initialize AmgX library
        subroutine AMGX_initialize() bind(C, name="AMGX_initialize")
        end subroutine
        
        ! Finalize AmgX library  
        subroutine AMGX_finalize() bind(C, name="AMGX_finalize")
        end subroutine
        
        ! Install signal handler for debugging
        subroutine AMGX_install_signal_handler() bind(C, name="AMGX_install_signal_handler")
        end subroutine
        
        ! Register print callback for verbose output
        subroutine AMGX_register_print_callback(func) bind(C, name="AMGX_register_print_callback")
            import :: c_funptr
            type(c_funptr), value :: func
        end subroutine
        
        ! Create config from string
        function AMGX_config_create(config, options, cfg_file) bind(C, name="AMGX_config_create") result(rc)
            import :: c_ptr, c_char, c_int
            type(c_ptr) :: config
            character(kind=c_char) :: options(*), cfg_file(*)
            integer(c_int) :: rc
        end function
        
        ! Create config from file
        function AMGX_config_create_from_file(config, filename) bind(C, name="AMGX_config_create_from_file") result(rc)
            import :: c_ptr, c_char, c_int
            type(c_ptr) :: config
            character(kind=c_char) :: filename(*)
            integer(c_int) :: rc
        end function
        
        ! Create resources (simple version - uses MPI_COMM_WORLD)
        function AMGX_resources_create_simple(resources, config) bind(C, name="AMGX_resources_create_simple") result(rc)
            import :: c_ptr, c_int
            type(c_ptr) :: resources
            type(c_ptr), value :: config
            integer(c_int) :: rc
        end function
        
        ! Create resources with custom MPI communicator
        function AMGX_resources_create_wrapper(resources, config, comm, num_devices, device_ids) bind(C, name="AMGX_resources_create_wrapper") result(rc)
            import :: c_ptr, c_int
            type(c_ptr) :: resources
            type(c_ptr), value :: config
            integer :: comm  ! MPI communicator (pass address)
            integer(c_int), value :: num_devices
            integer(c_int) :: device_ids(*)
            integer(c_int) :: rc
        end function
        
        ! Create matrix
        function AMGX_matrix_create(matrix, resources, mode) bind(C, name="AMGX_matrix_create") result(rc)
            import :: c_ptr, c_int
            type(c_ptr) :: matrix
            type(c_ptr), value :: resources
            integer(c_int), value :: mode
            integer(c_int) :: rc
        end function
        
        ! Upload matrix (CSR format) - local version
        function AMGX_matrix_upload_all(matrix, n, nnz, block_dimx, block_dimy, row_ptrs, col_indices, &
                                         data, diag_data) bind(C, name="AMGX_matrix_upload_all") result(rc)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: matrix
            integer(c_int), value :: n, nnz, block_dimx, block_dimy
            integer(c_int) :: row_ptrs(*), col_indices(*)
            real(c_double) :: data(*)
            type(c_ptr), value :: diag_data
            integer(c_int) :: rc
        end function
        
        ! Upload matrix (CSR format) - global/distributed version
        function AMGX_matrix_upload_all_global(matrix, n_global, n, nnz, block_dimx, block_dimy, &
                                                row_ptrs, col_indices_global, data, diag_data, &
                                                allocated_halo_depth, num_import_rings, &
                                                partition_vector) bind(C, name="AMGX_matrix_upload_all_global_32") result(rc)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: matrix
            integer(c_int), value :: n_global  ! Total number of rows across all ranks
            integer(c_int), value :: n         ! Number of local rows
            integer(c_int), value :: allocated_halo_depth
            integer(c_int), value :: num_import_rings
            integer(c_int), value :: nnz, block_dimx, block_dimy
            integer(c_int) :: row_ptrs(*), col_indices_global(*)
            real(c_double) :: data(*)
            type(c_ptr), value :: diag_data
            integer(c_int) :: partition_vector(*)  ! Size n+1, contains global row indices
            integer(c_int) :: rc
        end function
        
        ! Create vector
        function AMGX_vector_create(vec, resources, mode) bind(C, name="AMGX_vector_create") result(rc)
            import :: c_ptr, c_int
            type(c_ptr) :: vec
            type(c_ptr), value :: resources
            integer(c_int), value :: mode
            integer(c_int) :: rc
        end function
        
        ! Upload vector
        function AMGX_vector_upload(vec, n, block_dim, data) bind(C, name="AMGX_vector_upload") result(rc)
            import :: c_ptr, c_int, c_double
            type(c_ptr), value :: vec
            integer(c_int), value :: n, block_dim
            real(c_double) :: data(*)
            integer(c_int) :: rc
        end function
        
        ! Zero out vector
        function AMGX_vector_set_zero(vec, n, block_dim) bind(C, name="AMGX_vector_set_zero") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: vec
            integer(c_int), value :: n, block_dim
            integer(c_int) :: rc
        end function

        ! Download vector
        function AMGX_vector_download(vec, data) bind(C, name="AMGX_vector_download") result(rc)
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: vec
            real(c_double) :: data(*)
            integer(c_int) :: rc
        end function
        
        ! Get vector size (including halo elements for distributed case)
        function AMGX_vector_get_size(vec, n, block_dim) bind(C, name="AMGX_vector_get_size") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: vec
            integer(c_int) :: n, block_dim
            integer(c_int) :: rc
        end function
        
        ! Bind vector to matrix (required before using vector in distributed mode)
        function AMGX_vector_bind(vec, matrix) bind(C, name="AMGX_vector_bind") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: vec, matrix
            integer(c_int) :: rc
        end function
        
        ! Create solver
        function AMGX_solver_create(solver, resources, mode, config) bind(C, name="AMGX_solver_create") result(rc)
            import :: c_ptr, c_int
            type(c_ptr) :: solver
            type(c_ptr), value :: resources, config
            integer(c_int), value :: mode
            integer(c_int) :: rc
        end function
        
        ! Setup solver
        function AMGX_solver_setup(solver, matrix) bind(C, name="AMGX_solver_setup") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: solver, matrix
            integer(c_int) :: rc
        end function
        
        ! Solve
        function AMGX_solver_solve(solver, rhs, solution) bind(C, name="AMGX_solver_solve") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: solver, rhs, solution
            integer(c_int) :: rc
        end function
        
        ! Reset solver (frees preconditioner but keeps setup)
        function AMGX_solver_resetup(solver, matrix) bind(C, name="AMGX_solver_resetup") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: solver, matrix
            integer(c_int) :: rc
        end function
        
        ! Destroy objects
        function AMGX_solver_destroy(solver) bind(C, name="AMGX_solver_destroy") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int) :: rc
        end function
        
        function AMGX_matrix_destroy(matrix) bind(C, name="AMGX_matrix_destroy") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: matrix
            integer(c_int) :: rc
        end function
        
        function AMGX_vector_destroy(vec) bind(C, name="AMGX_vector_destroy") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: vec
            integer(c_int) :: rc
        end function
        
        function AMGX_resources_destroy(resources) bind(C, name="AMGX_resources_destroy") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: resources
            integer(c_int) :: rc
        end function
        
        function AMGX_config_destroy(config) bind(C, name="AMGX_config_destroy") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: config
            integer(c_int) :: rc
        end function
        
        ! Pin/unpin memory for efficient GPU transfers
        function AMGX_pin_memory(ptr, bytes) bind(C, name="AMGX_pin_memory") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            integer(c_int), value :: bytes
            integer(c_int) :: rc
        end function
        
        function AMGX_unpin_memory(ptr) bind(C, name="AMGX_unpin_memory") result(rc)
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            integer(c_int) :: rc
        end function
    end interface
    
    ! AmgX modes
    integer(c_int), parameter :: AMGX_mode_dDDI = 8193  ! Double precision (8193), distributed, device (GPU), index int (floating point precision: 8465)

contains

    !>------------------------------------------------------------
    !! Initialize AmgX solver
    !!------------------------------------------------------------
    subroutine init_amgx(compute_comm)
        implicit none
        type(MPI_Comm), intent(in), optional :: compute_comm
        integer(c_int) :: rc
        integer :: rank, ierr, num_devices, device_num
        integer(c_int), allocatable :: device_ids(:)
        type(c_ptr) :: mpi_comm_ptr  ! C pointer for MPI communicator
        type(MPI_Comm) :: amgx_comm
        integer(c_intptr_t) :: comm_handle  ! Integer handle sized for pointer conversion
        character(len=1024) :: config_string
        character(len=1024) :: solver_file = 'amgx_solver_config.json'
        logical :: file_exists = .false.

        if (initialized) return
        
        ! Determine which MPI communicator to use
        if (present(compute_comm)) then
            amgx_comm = compute_comm
        else
            amgx_comm = MPI_COMM_WORLD
        endif
                
        ! Get the OpenACC device number
        device_num = acc_get_device_num(acc_device_nvidia)
        
        ! Initialize AmgX library - MUST be called only once globally by rank 0
        ! Use WORLD communicator for this, not the compute communicator
        call MPI_Comm_rank(amgx_comm, rank, ierr)
        call AMGX_initialize()
        call MPI_Barrier(amgx_comm, ierr)
        
        ! Install error handler for verbose output
        call AMGX_install_signal_handler()
        call AMGX_register_print_callback(c_funloc(amgx_print_callback))

        ! Create configuration with verbose error reporting
        ! Using PBICGSTAB with single Block Jacobi preconditioner (robust and fast on GPU)
        config_string = &
            "config_version=2, " // &
            "solver(main)=PBICGSTAB, " // &
            "main:preconditioner(prec)=BLOCK_JACOBI, " // &
            "prec:relaxation_factor=1, " // &
            "prec:max_iters=1, " // &
            "main:use_scalar_norm=1, " // &
            "main:max_iters=1000, " // &
            "main:convergence=RELATIVE_INI_CORE, " // &
            "main:tolerance=1e-3, " // &
            "main:monitor_residual=1, " // &
            "main:store_res_history=1, " // &
            "main:obtain_timings=1, " // &
            "main:norm=L2" // c_null_char
        
        !check if the solver configuration file exists
        inquire(file=solver_file,exist=file_exists)
        if (file_exists) then
            if (STD_OUT_PE) print*, "Using AmgX solver configuration file: ", trim(solver_file)
            rc = AMGX_config_create_from_file(amgx_config, trim(solver_file))
        else
            rc = AMGX_config_create(amgx_config, trim(config_string), c_null_char)
        end if

        ! Create resources using the OpenACC device number
        ! This ensures AmgX uses the same GPU that OpenACC is using
        ! In distributed mode, each MPI rank uses exactly 1 GPU
        num_devices = 1
        allocate(device_ids(num_devices))
        device_ids(1) = int(device_num, c_int)
                
        ! Pass pointer to integer MPI communicator handle (not the F08 derived type)
        ! Use amgx_comm which was set to either compute_comm or MPI_COMM_WORLD
        rc = AMGX_resources_create_wrapper(amgx_resources, amgx_config, amgx_comm%MPI_VAL, num_devices, device_ids)
        if (rc /= 0) then
            print*, "ERROR: AMGX_resources_create failed with rc=", rc
        endif
        
        deallocate(device_ids)
        
        ! Create matrix, vectors, and solver
        rc = AMGX_matrix_create(amgx_matrix, amgx_resources, AMGX_mode_dDDI)
        rc = AMGX_vector_create(amgx_vector_rhs, amgx_resources, AMGX_mode_dDDI)
        rc = AMGX_vector_create(amgx_vector_solution, amgx_resources, AMGX_mode_dDDI)
        rc = AMGX_solver_create(amgx_solver, amgx_resources, AMGX_mode_dDDI, amgx_config)
        
        initialized = .True.
    end subroutine
    !>------------------------------------------------------------
    !! Initialize AmgX solver
    !!------------------------------------------------------------
    subroutine init_iter_winds_amgx(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: ierr, rank, nprocs
        integer :: local_n


        call finalize_iter_winds_amgx()  ! Clean up any previous AmgX state

        ! Initialize module variables
        call init_module_vars(domain)
        
        ! Calculate local number of rows (will be set properly later in build_csr_matrix)
        local_n = xm * zm * ym
        
        ! Get MPI information
        call MPI_Comm_rank(domain%compute_comms, rank, ierr)
        
        ! Calculate total number of rows globally
        n_rows_global = mx * my * mz
                
        if (STD_OUT_PE) then
            print*, "AmgX initialization: rank", rank, "local_n =", local_n, &
                    "local_offset =", local_offset, "global_n =", global_n
        endif
        
        ! Initialize AmgX with the compute communicator (only compute ranks participate)
        call init_amgx(domain%compute_comms)
        
        if (STD_OUT_PE) print*, "AmgX solver initialized successfully"
        
    end subroutine init_iter_winds_amgx


    subroutine calc_nnz()
        implicit none
        integer :: i, j, k

        nnz = 0
        ! Loop over local domain INCLUDING boundary cells
        do j = ys, (ys+ym-1)
            do k = zs, (zs+zm-1)
                do i = xs, (xs+xm-1)                    
                    ! Lateral boundary conditions (identity matrix: λ = 0)
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1) then
                        nnz = nnz + 1
                    
                    ! Top boundary condition (k = mz-1): dλ/dz = 0
                    else if (k >= mz-1) then
                        nnz = nnz + 2
                    
                    ! Bottom boundary condition (k = 0): terrain-following BC
                    else if (k <= 0) then
                        nnz = nnz + 10
                    ! Interior points: full 15-point stencil
                    else
                        ! Diagonal corners in k-1 plane
                        nnz = nnz + 15
                    endif
                enddo
            enddo
        enddo
    end subroutine calc_nnz



    !>------------------------------------------------------------
    !! Build sparse matrix in CSR format
    !!------------------------------------------------------------
    subroutine build_csr_matrix(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, row, cnt, global_row
        integer(c_int) :: rc
        real :: denom
        integer :: ierr, rank
        integer, allocatable :: partition_vec_fint(:)

        ! Get MPI rank for partition vector
        call MPI_Comm_rank(domain%compute_comms, rank, ierr)
        
        ! Calculate number of rows INCLUDING boundary conditions layer
        n_rows = xm * zm * ym
        call calc_nnz()

        ! Allocate CSR arrays
        if (allocated(row_ptrs)) deallocate(row_ptrs)
        if (allocated(col_indices)) deallocate(col_indices)
        if (allocated(values)) deallocate(values)
        if (allocated(partition_vec)) deallocate(partition_vec)
        
        allocate(row_ptrs(n_rows + 1))
        allocate(col_indices(nnz))
        allocate(values(nnz))
        allocate(partition_vec(n_rows_global))
        allocate(partition_vec_fint(n_rows_global))

        ! NOTE: Do NOT use AMGX_pin_memory on these arrays!
        ! AMGX uses Thrust internally for H->D transfers during matrix upload,
        ! and pinned memory registration conflicts with Thrust's expectations,
        ! causing "cudaErrorInvalidValue" in trivial_device_copy.
        ! AMGX handles memory transfers efficiently without manual pinning.
        
        ! Build partition vector: assign each global row to its owning rank
        ! Each rank fills in only its portion of the global partition vector
        ! Note: This needs to be gathered across all ranks to form the complete partition vector
        ! For simplicity in distributed mode, we'll use a contiguous block partitioning
        ! where each rank owns a contiguous block of global rows
        
        ! Initialize entire partition vector (will be gathered later)
        partition_vec_fint = 0
        
        ! Each rank marks which global rows it owns with its rank number
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1
                    ! Calculate global row index using same formula as column indices
                    global_row = j * mx * mz + k * mx + i + 1
                    ! if (global_row >= 1 .and. global_row <= n_rows_global) then
                        partition_vec_fint(global_row) = rank
                    ! endif
                enddo
            enddo
        enddo
        
        ! Gather partition vector across all ranks using MPI_Allreduce with MAX
        ! (each rank sets its owned rows, others are 0, so MAX gives the owner)
        call MPI_Allreduce(MPI_IN_PLACE, partition_vec_fint, n_rows_global, MPI_INTEGER, &
                          MPI_MAX, domain%compute_comms, ierr)

        !convert partition_vec to c_int
        partition_vec = int(partition_vec_fint, c_int)
        
        ! Build CSR matrix with proper boundary conditions
        cnt = 0
        row = 1
        row_ptrs(1) = 0  ! C-style indexing (0-based)
        
        ! Loop over local domain INCLUDING boundary cells
        do j = ys, ys+ym-1
            do k = zs, zs+zm-1
                do i = xs, xs+xm-1                    
                    ! Lateral boundary conditions (identity matrix: λ = 0)
                    if (i <= 0 .or. j <= 0 .or. i >= mx-1 .or. j >= my-1) then
                        call add_entry(i, k, j, 1.0, cnt)
                    
                    ! Top boundary condition (k = mz-1): dλ/dz = 0
                    else if (k >= mz-1) then
                        ! λ(k) - λ(k-1) = 0
                        call add_entry(i, k, j, 1.0/dz_if(i,k,j), cnt)
                        call add_entry(i, k-1, j, -1.0/dz_if(i,k,j), cnt)
                    
                    ! Bottom boundary condition (k = 0): terrain-following BC
                    else if (k <= 0) then
                        denom = 2.0*(dzdx_surf(i,j)**2 + dzdy_surf(i,j)**2 + alpha(i,1,j)**2)/jaco(i,1,j)
                        
                        ! Vertical derivative terms
                        call add_entry(i, k+1, j, 1.0/dz_if(i,k+1,j), cnt)
                        call add_entry(i, k, j, -1.0/dz_if(i,k+1,j), cnt)
                        
                        ! Lateral derivative terms (note: signs are opposite as per comment in PETSc version)
                        ! x-derivatives at k and k+1
                        call add_entry(i-1, k, j, dzdx_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i+1, k, j, -dzdx_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i-1, k+1, j, dzdx_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i+1, k+1, j, -dzdx_surf(i,j)/(denom*2.0*dx), cnt)
                        
                        ! y-derivatives at k and k+1
                        call add_entry(i, k, j-1, dzdy_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i, k, j+1, -dzdy_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i, k+1, j-1, dzdy_surf(i,j)/(denom*2.0*dx), cnt)
                        call add_entry(i, k+1, j+1, -dzdy_surf(i,j)/(denom*2.0*dx), cnt)
                    
                    ! Interior points: full 15-point stencil
                    else
                        ! Diagonal corners in k-1 plane
                        call add_entry(i, k-1, j-1, O_coef(i,k,j), cnt)
                        call add_entry(i-1, k-1, j, K_coef(i,k,j), cnt)
                        
                        ! Diagonal corners in k+1 plane
                        call add_entry(i, k+1, j-1, M_coef(i,k,j), cnt)
                        call add_entry(i-1, k+1, j, I_coef(i,k,j), cnt)
                        
                        ! k-1 level
                        call add_entry(i, k-1, j, C_coef(i,k,j), cnt)
                        
                        ! k level (central plane)
                        call add_entry(i, k, j-1, G_coef(i,k,j), cnt)
                        call add_entry(i-1, k, j, E_coef(i,k,j), cnt)
                        call add_entry(i, k, j, A_coef(i,k,j), cnt)  ! Diagonal
                        call add_entry(i+1, k, j, D_coef(i,k,j), cnt)
                        call add_entry(i, k, j+1, F_coef(i,k,j), cnt)
                        
                        ! k+1 level
                        call add_entry(i, k+1, j, B_coef(i,k,j), cnt)
                        
                        ! More diagonal corners
                        call add_entry(i+1, k-1, j, J_coef(i,k,j), cnt)
                        call add_entry(i, k-1, j+1, N_coef(i,k,j), cnt)
                        call add_entry(i+1, k+1, j, H_coef(i,k,j), cnt)
                        call add_entry(i, k+1, j+1, L_coef(i,k,j), cnt)
                    endif
                    row = row + 1
                    row_ptrs(row) = cnt
                enddo
            enddo
        enddo
                
    end subroutine build_csr_matrix
    
    !>------------------------------------------------------------
    !! Helper to add matrix entry
    !!------------------------------------------------------------
    subroutine add_entry(i, k, j, coef_val, cnt)
        implicit none
        integer, intent(in) :: i, k, j
        integer, intent(inout) :: cnt
        real, intent(in) :: coef_val
        integer :: global_col
        
        ! Convert (i,k,j) to GLOBAL column index
        ! The grid is indexed in z-major order: global_idx = j*mx*mz + k*mx + i
        ! where mx, mz are the GLOBAL dimensions (same on all ranks)
        global_col = j * (mx) * (mz) + &
                     k * (mx) + &
                     i

        cnt = cnt + 1
        col_indices(cnt) = int(global_col, c_int)
        values(cnt) = real(coef_val, c_double)
        
    end subroutine add_entry

    !>------------------------------------------------------------
    !! Main solver routine
    !!------------------------------------------------------------
    subroutine calc_iter_winds_amgx(domain, alpha_in, div_in, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, dimension(ims:ime,domain%kms:domain%kme,jms:jme), intent(in) :: alpha_in, div_in
        logical, intent(in) :: adv_den
        
        integer(c_int) :: rc
        integer(c_int) :: n_rows_with_halo, block_dim
        logical :: varying_alpha
        integer :: i, j, k, idx
        
        ! Copy input data (on host - AmgX will manage GPU transfers)
        !$acc parallel loop gang vector collapse(3) present(alpha,div,alpha_in,div_in)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    div(i,k,j) = div_in(i,k,j)
                    alpha(i,k,j) = alpha_in(i,k,j)
                enddo
            enddo
        enddo

        !$acc update host(div, alpha)

        ! Initialize or update coefficients
        varying_alpha = .not.( ALL(alpha==minval(alpha)) )

        ! Build RHS vector
        ! call compute_rhs(domain)

        if (.not. allocated(A_coef) .or. varying_alpha) then
            if (.not. allocated(A_coef)) then
                call initialize_coefs(domain)
            else
                call update_coefs(domain)
            endif

            call build_csr_matrix(domain)
            
            ! Upload matrix to AmgX
            ! rc = AMGX_matrix_upload_all(amgx_matrix, n_rows, nnz, 1, 1, &
            !                              row_ptrs, col_indices, values, c_null_ptr)
            ! Upload matrix to AmgX - use array slicing to pass only valid data
            ! (arrays may be over-allocated from initial estimate)
            ! Use the partition_vec we built in build_csr_matrix

            !print out size information for all arrays and n_rows, nnz
                !print first index where partition_vec == 1
            do i = 1, size(partition_vec)
                if (partition_vec(i) == 1) then
                    print*, "  First index where partition_vec == 1: ", i
                    exit
                endif
            enddo
            print*, "  ----------------------------------------"
            rc = AMGX_matrix_upload_all_global(amgx_matrix, n_rows_global, n_rows, nnz, 1, 1, &
                                                row_ptrs, col_indices, &
                                                values, c_null_ptr, 1, 1, partition_vec)

            if (rc /= 0) then
                if (STD_OUT_PE) print*, "ERROR: AMGX_matrix_upload_all_global failed with rc=", rc
            endif

            
            ! Bind vectors to matrix (required for distributed mode to set up halo exchange)
            rc = AMGX_vector_bind(amgx_vector_rhs, amgx_matrix)
            if (rc /= 0) then
                if (STD_OUT_PE) print*, "ERROR: AMGX_vector_bind(rhs) failed with rc=", rc
            endif
            
            rc = AMGX_vector_bind(amgx_vector_solution, amgx_matrix)
            if (rc /= 0) then
                if (STD_OUT_PE) print*, "ERROR: AMGX_vector_bind(solution) failed with rc=", rc
            endif
            
            ! Allocate work arrays with local size
            if (.not. allocated(rhs)) then
                allocate(rhs(n_rows))
                ! NOTE: Not using AMGX_pin_memory - see note above in build_csr_matrix
            endif
            
            if (.not. allocated(solution)) then
                allocate(solution(n_rows))
                solution = 0.0_c_double
                ! NOTE: Not using AMGX_pin_memory - see note above in build_csr_matrix
            endif
                                                            
        endif

        ! Zero out the work arrays
        rhs = 0.0_c_double
        ! Don't zero solution - reuse from previous timestep as initial guess
        ! solution = 0.0_c_double  ! Commented out for warm start
        
        ! Build RHS vector (fills only first n_rows elements)
        call compute_rhs(domain)

        ! Upload RHS and solution vectors
        ! Upload only n_rows (local), AMGX handles halo exchange
        ! Solution vector retains previous timestep values (warm start)
        rc = AMGX_vector_upload(amgx_vector_rhs, n_rows, 1, rhs)
        if (rc /= 0) then
            if (STD_OUT_PE) print*, "ERROR: AMGX_vector_upload(rhs) failed with rc=", rc
        endif

        rc = AMGX_vector_upload(amgx_vector_solution, n_rows, 1, solution)
        if (rc /= 0) then
            if (STD_OUT_PE) print*, "ERROR: AMGX_vector_upload(solution) failed with rc=", rc
        endif


        if (first_solve) then
            ! call MPI_Barrier(domain%compute_comms, rc)
            ! Setup solver
            rc = AMGX_solver_setup(amgx_solver, amgx_matrix)
            if (rc /= 0) then
                if (STD_OUT_PE) print*, "ERROR: AMGX_solver_setup failed with rc=", rc
            endif
            ! call MPI_Barrier(domain%compute_comms, rc)
            first_solve = .False.
        endif
        
        ! Solve
        if (STD_OUT_PE) print*, "Solving with AmgX..."
        rc = AMGX_solver_solve(amgx_solver, amgx_vector_rhs, amgx_vector_solution)
        if (rc /= 0) then
            if (STD_OUT_PE) print*, "ERROR: AMGX_solver_solve failed with rc=", rc
        endif
        
        ! Download solution
        rc = AMGX_vector_download(amgx_vector_solution, solution)
        if (rc /= 0) then
            if (STD_OUT_PE) print*, "ERROR: AMGX_vector_download failed with rc=", rc
        endif
                
        ! Update wind field
        call calc_updated_winds(domain, solution, adv_den)
        
        if (STD_OUT_PE) print*, "AmgX solve complete"
        
    end subroutine calc_iter_winds_amgx

    !>------------------------------------------------------------
    !! Compute RHS vector (same as PETSc version)
    !!------------------------------------------------------------
    subroutine compute_rhs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, idx
        
        ! rhs is now pre-allocated with global size, we only fill first n_rows elements
        
        idx = 1
        do j=ys,(ys+ym-1)
        do k=zs,(zs+zm-1)
            do i=xs,(xs+xm-1)
                    !For global boundary conditions
                    if (i.le.0 .or. j.le.0 .or. k.eq.0 .or. &
                        i.ge.mx-1 .or. j.ge.my-1 .or. k.eq.mz-1) then
                        rhs(idx) = 0.0
                    else
                        rhs(idx) = -2*div(i,k,j)
                    endif
                    idx = idx + 1
                enddo
            enddo
        enddo
        
    end subroutine compute_rhs

    !>------------------------------------------------------------
    !! Update winds from solution (simplified version)
    !!------------------------------------------------------------
    subroutine calc_updated_winds(domain, lambda, adv_den)
        implicit none
        type(domain_t), intent(inout) :: domain
        real(c_double), intent(in) :: lambda(1:n_rows)
        logical, intent(in) :: adv_den
        
        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, dlambdz, u_temp, v_temp, lambda_3d, rho, rho_u, rho_v
        integer :: i, j, k, i_start, i_end, j_start, j_end, idx

        i_start = i_s
        i_end   = i_e+1
        j_start = j_s
        j_end   = j_e+1

        allocate(u_temp(i_start:i_end,k_s-1:k_e+1,j_s:j_e))
        allocate(v_temp(i_s:i_e,k_s-1:k_e+1,j_start:j_end))
        allocate(lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1))

        allocate(u_dlambdz(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,k_s:k_e,j_start:j_end))
        allocate(dlambdz(i_s-1:i_e+1,k_s:k_e,j_s-1:j_e+1))

        allocate(rho(domain%ims:domain%ime,k_s:k_e,domain%jms:domain%jme))
        allocate(rho_u(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(rho_v(i_s:i_e,k_s:k_e,j_start:j_end))
        
        ! Initialize lambda_3d with zeros for extended bounds
        lambda_3d = 0.0
        
        ! Reshape 1D solution array to 3D (including boundary condition layer)
        idx = 1
        ! do j = j_s-1, j_e+1
        !     do k = k_s-1, k_e+1
        !         do i = i_s-1, i_e+1
        !             lambda_3d(i,k,j) = real(lambda(idx))
        !             idx = idx + 1
        !         end do
        !     end do
        ! end do
        do j = ys, (ys+ym-1)
            do k = zs, (zs+zm-1)
                do i = xs, (xs+xm-1)
                    lambda_3d(i,k,j) = real(lambda(idx))
                    idx = idx + 1
                end do
            end do
        end do
        
        ! Exchange halos for lambda_3d using simple MPI communication
        ! This fills the ghost cells (i_s-1, i_e+1, j_s-1, j_e+1) needed for stencil operations
        call exchange_lambda_halos(lambda_3d, domain)

        associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                    u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                    v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                    w       => domain%vars_3d(domain%var_indx(kVARS%w_real)%v)%data_3d, &
                    alpha  => domain%vars_3d(domain%var_indx(kVARS%wind_alpha)%v)%data_3d, &
                    jaco_u_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                    jaco_v_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                    jaco_domain   => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                    dzdx_u => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                    dzdy_v => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d)
        ! Create GPU data for local arrays

        !$acc enter data create(u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v)
        !$acc data present(density, u, v, jaco_u_domain, jaco_v_domain, jaco_domain, dzdx_u, dzdy_v, u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v, alpha, dz_if) copyin(lambda_3d)

        !$acc kernels
        rho = 1.0
        rho_u = 1.0
        rho_v = 1.0
        !$acc end kernels

        if (adv_den) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms, jme
                do k = k_s, k_e
                    do i = ims, ime
                        rho(i,k,j) = density(i,k,j)
                    end do
                end do
            end do
        endif
        
        ! Calculate rho_u and rho_v with GPU kernels
        !$acc parallel
        if (i_s==ids .and. i_e==ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end-1
                        rho_u(i,k,j) = 0.5*(rho(i,k,j) + rho(i-1,k,j))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j) = rho(i_end-1,k,j)
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                end do
            end do
        else if (i_s==ids) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start+1, i_end
                        rho_u(i,k,j) = 0.5*(rho(i,k,j) + rho(i-1,k,j))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_start,k,j) = rho(i_start,k,j)
                end do
            end do
        else if (i_e==ide) then
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end-1
                        rho_u(i,k,j) = 0.5*(rho(i,k,j) + rho(i-1,k,j))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do j = j_s, j_e
                do k = k_s, k_e
                    rho_u(i_end,k,j) = rho(i_end-1,k,j)
                end do
            end do
        else
            !$acc loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_start, i_end
                        rho_u(i,k,j) = 0.5*(rho(i,k,j) + rho(i-1,k,j))
                    end do
                end do
            end do
        endif
        !$acc end parallel
        
        !$acc parallel
        if (j_s==jds .and. j_e==jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5*(rho(i,k,j) + rho(i,k,j-1))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                    rho_v(i,k,j_end) = rho(i,k,j_end-1)
                end do
            end do
        else if (j_s==jds) then
            !$acc loop gang vector collapse(3)
            do j = j_start+1, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5*(rho(i,k,j) + rho(i,k,j-1))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_start) = rho(i,k,j_start)
                end do
            end do
        else if (j_e==jde) then
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end-1
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5*(rho(i,k,j) + rho(i,k,j-1))
                    end do
                end do
            end do
            !$acc loop gang vector collapse(2)
            do k = k_s, k_e
                do i = i_s, i_e
                    rho_v(i,k,j_end) = rho(i,k,j_end-1)
                end do
            end do
        else
            !$acc loop gang vector collapse(3)
            do j = j_start, j_end
                do k = k_s, k_e
                    do i = i_s, i_e
                        rho_v(i,k,j) = 0.5*(rho(i,k,j) + rho(i,k,j-1))
                    end do
                end do
            end do
        endif
        !$acc end parallel
        
        !stager lambda to u grid - GPU computation
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s-1, k_e+1
                do i = i_start, i_end
                    u_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i-1,k,j)) / 2 
                end do
            end do
        end do

        !stager lambda to v grid - GPU computation
        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s-1, k_e+1
                do i = i_s, i_e
                    v_temp(i,k,j) = (lambda_3d(i,k,j) + lambda_3d(i,k,j-1)) / 2 
                end do
            end do
        end do

        !divide dz differences by dz. Note that dz will be horizontally constant
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k=k_s,k_e
                do i = i_start, i_end
                    u_dlambdz(i,k,j) = u_temp(i,k+1,j) - u_temp(i,k-1,j)
                    u_dlambdz(i,k,j) = u_dlambdz(i,k,j)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
                end do
            end do
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k=k_s,k_e
                do i = i_s, i_e
                    v_dlambdz(i,k,j) = v_temp(i,k+1,j) - v_temp(i,k-1,j)
                    v_dlambdz(i,k,j) = v_dlambdz(i,k,j)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
                end do
            end do
        enddo

        !$acc parallel loop gang vector collapse(3)
        do j = j_s-1, j_e+1
            do k=k_s,k_e
                do i = i_s-1, i_e+1
                    dlambdz(i,k,j) = lambda_3d(i,k+1,j) - lambda_3d(i,k-1,j)
                    dlambdz(i,k,j) = dlambdz(i,k,j)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
                end do
            end do
        enddo

        ! Update wind fields on GPU
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u(i,k,j) = u(i,k,j) + 0.5*((lambda_3d(i,k,j) - lambda_3d(i-1,k,j))/dx - &
                                        dzdx_u(i,k,j)*(u_dlambdz(i,k,j))/jaco_u_domain(i,k,j))/(rho_u(i,k,j))
                end do
            end do
        end do
        

        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v(i,k,j) = v(i,k,j) + 0.5*((lambda_3d(i,k,j) - lambda_3d(i,k,j-1))/dx - dzdy_v(i,k,j)*(v_dlambdz(i,k,j))/jaco_v_domain(i,k,j))/(rho_v(i,k,j))
                end do
            end do
        end do
        
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    w(i,k,j) = w(i,k,j) + 0.5*(alpha(i,k,j)**2)*dlambdz(i,k,j)/jaco_domain(i,k,j)
                end do
            end do
        end do

        !$acc end data
        !$acc exit data delete(u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v)
        end associate     
    end subroutine calc_updated_winds

    !>------------------------------------------------------------
    !! Initialize module variables (same as PETSc version)
    !!------------------------------------------------------------
    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        integer :: i, j, k, i_s_bnd, i_e_bnd, j_s_bnd, j_e_bnd
        
        i_s = domain%its
        i_e = domain%ite
        k_s = domain%kts
        k_e = domain%kte
        j_s = domain%jts
        j_e = domain%jte
        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        
        ! Adjust for global boundaries
        if (ims==ids) i_s = ids
        if (ime==ide) i_e = ide
        if (jms==jds) j_s = jds
        if (jme==jde) j_e = jde
        
        ! Create boundary layer: extend domain by 1 cell on each face
        ! This matches PETSc's (ide+2) × (jde+2) × (kde+2) convention
        xs = i_s - merge(1, 0, i_s==ids)  ! Add left boundary if at global left edge
        ys = j_s - merge(1, 0, j_s==jds)  ! Add south boundary if at global south edge
        zs = 0                             ! Bottom boundary always at k=0
        
        mx = ide + 2  ! Total x-dimension with boundary layers (0 to ide+1)
        my = jde + 2  ! Total y-dimension with boundary layers (0 to jde+1)
        mz = domain%kde + 2  ! Total z-dimension with boundary layers (0 to kde+1)
        
        xm = (i_e + merge(1, 0, i_e==ide)) - xs + 1  ! Local x extent including boundaries
        ym = (j_e + merge(1, 0, j_e==jde)) - ys + 1  ! Local y extent including boundaries
        zm = mz - zs  ! Local z extent (always includes both top and bottom)
        
        ! Calculate interior boundaries for finite differences
        i_s_bnd = i_s
        i_e_bnd = i_e
        j_s_bnd = j_s
        j_e_bnd = j_e
        if (i_s==ids) i_s_bnd = i_s + 1
        if (i_e==ide) i_e_bnd = i_e - 1
        if (j_s==jds) j_s_bnd = j_s + 1
        if (j_e==jde) j_e_bnd = j_e - 1
        
        hs = domain%grid%halo_size
        
        ! Allocate arrays if not already allocated
        if (.not. allocated(dzdx)) then
            allocate(dzdx(i_s:i_e, k_s:k_e, j_s:j_e))
            allocate(dzdy(i_s:i_e, k_s:k_e, j_s:j_e))
            allocate(jaco(i_s:i_e, k_s:k_e, j_s:j_e))
            
            allocate(dzdx_surf(i_s:i_e, j_s:j_e))
            allocate(dzdy_surf(i_s:i_e, j_s:j_e))
            
            allocate(sigma(i_s:i_e, k_s:k_e, j_s:j_e))
            allocate(dz_if(i_s:i_e, k_s:k_e+1, j_s:j_e))
            allocate(alpha(i_s:i_e, k_s:k_e, j_s:j_e))
            allocate(div(i_s:i_e, k_s:k_e, j_s:j_e))
            dx = domain%dx
            
            ! Initialize geometric arrays from domain
            associate(advection_dz_var => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      neighbor_terrain_var => domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d, &
                      dzdx_domain => domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                      dzdy_domain => domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, &
                      jaco_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)
            
            ! Compute dz_if (interface vertical spacing)
            do j = j_s, j_e
                do k = k_s+1, k_e
                    do i = i_s, i_e
                        dz_if(i,k,j) = (advection_dz_var(i,k,j) + advection_dz_var(i,k-1,j))/2
                    end do
                end do
            end do
            
            ! Boundary vertical spacing
            do j = j_s, j_e
                do i = i_s, i_e
                    dz_if(i,k_s,j) = advection_dz_var(i,k_s,j)
                    dz_if(i,k_e+1,j) = advection_dz_var(i,k_e,j)
                    
                    dzdx_surf(i,j) = dzdx_domain(i,k_s,j)
                    dzdy_surf(i,j) = dzdy_domain(i,k_s,j)
                enddo
            end do
            
            ! Copy terrain gradients and jacobian
            do j = j_s, j_e
                do k = k_s, k_e
                    do i = i_s, i_e
                        dzdx(i,k,j) = dzdx_domain(i,k,j)
                        dzdy(i,k,j) = dzdy_domain(i,k,j)
                        jaco(i,k,j) = jaco_domain(i,k,j)
                        sigma(i,k,j) = dz_if(i,k,j)/dz_if(i,k+1,j)
                    end do
                end do
            end do
            
            ! Update surface derivatives for interior points
            do j = j_s, j_e
                do i = i_s_bnd, i_e_bnd
                    dzdx_surf(i,j) = (neighbor_terrain_var(i+1,j) - neighbor_terrain_var(i-1,j))/(2*dx)
                end do
            end do
            
            do j = j_s_bnd, j_e_bnd
                do i = i_s, i_e
                    dzdy_surf(i,j) = (neighbor_terrain_var(i,j+1) - neighbor_terrain_var(i,j-1))/(2*dx)
                end do
            end do
            
            end associate
        endif

        !$acc enter data copyin(dz_if, div, alpha)
        
    end subroutine init_module_vars

    !>------------------------------------------------------------
    !! Initialize coefficient arrays (placeholder)
    !!------------------------------------------------------------
    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        integer :: k
        real :: mixed_denom

        ! Allocate coefficient arrays
        allocate(A_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(B_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(C_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(D_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(E_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(F_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(G_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(H_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(I_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(J_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(K_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(L_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(N_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(O_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        
        
        ! Initialize diagonal coefficient
        A_coef = 1
        B_coef = 0
        C_coef = 0
        D_coef = 0
        E_coef = 0
        F_coef = 0
        G_coef = 0
        H_coef = 0
        I_coef = 0
        J_coef = 0
        K_coef = 0
        L_coef = 0
        M_coef = 0
        N_coef = 0
        O_coef = 0

        ! Second derivative coefficients in x and y
        D_coef = 1.0/(domain%dx**2)
        E_coef = 1.0/(domain%dx**2)
        F_coef = 1.0/(domain%dx**2)
        G_coef = 1.0/(domain%dx**2)

        ! Coefficients for mixed derivatives
        do k=k_s,k_e
            mixed_denom = 2*domain%dx*(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
            
            H_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            I_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            L_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            M_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom

            J_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            K_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
            N_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            O_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:)+domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,k,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,k,j_s:j_e))/mixed_denom
        enddo
        
        ! Update time-varying coefficients
        call update_coefs(domain)
        
    end subroutine initialize_coefs

    !>------------------------------------------------------------
    !! Update coefficient arrays (placeholder)
    !!------------------------------------------------------------
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        real, allocatable, dimension(:,:,:) :: mixed_denom, X_coef
        real, allocatable, dimension(:,:) :: M_up, M_dwn
        integer :: k
        
        allocate(mixed_denom(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(X_coef(i_s:i_e, k_s:k_e, j_s:j_e))
        allocate(M_up(i_s:i_e, j_s:j_e))
        allocate(M_dwn(i_s:i_e, j_s:j_e))

        ! Compute denominators for B and C coefficients
        mixed_denom = (dz_if(:,k_s+1:k_e+1,:)+dz_if(:,k_s:k_e,:))*2*domain%dx
        
        do k = k_s,k_e
            if (k == k_s) then
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*dzdx_surf!/jaco(i_s:i_e,k,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*dzdy_surf!/jaco(i_s:i_e,k,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + alpha(:,k,:)**2!/jaco(i_s:i_e,k,j_s:j_e)
            else if (k == k_e) then
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s) + &
                        alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            else
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k+1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s))/ &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k-1,j_s))/domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i_s:i_e,k-1,j_s:j_e)
            endif
            B_coef(:,k,:) = M_up/(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k+1,j_s))
            C_coef(:,k,:) = M_dwn/(jaco(:,k,:)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i_s,k,j_s)*dz_if(i_s,k,j_s))
        enddo
                
        A_coef = -4/(domain%dx**2) - B_coef - C_coef
        
        X_coef = -((domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e) - domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,:,j_s:j_e)) + &
                   (domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1) - domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s:j_e)))/mixed_denom
                   
        B_coef = B_coef + X_coef 
        C_coef = C_coef - X_coef
        
        deallocate(mixed_denom, X_coef, M_up, M_dwn)
        
    end subroutine update_coefs

    !>------------------------------------------------------------
    !! Cleanup
    !!------------------------------------------------------------
    subroutine finalize_iter_winds_amgx()
        implicit none
        integer(c_int) :: rc
        
        if (.not. initialized) return
        
        ! Destroy AmgX objects
        rc = AMGX_solver_destroy(amgx_solver)
        rc = AMGX_vector_destroy(amgx_vector_solution)
        rc = AMGX_vector_destroy(amgx_vector_rhs)
        rc = AMGX_matrix_destroy(amgx_matrix)
        rc = AMGX_resources_destroy(amgx_resources)
        rc = AMGX_config_destroy(amgx_config)
        
        
        ! Unpin and deallocate arrays
        ! Deallocate arrays (no unpinning needed since we don't pin them)
        if (allocated(row_ptrs)) deallocate(row_ptrs)
        if (allocated(col_indices)) deallocate(col_indices)
        if (allocated(values)) deallocate(values)
        if (allocated(rhs)) deallocate(rhs)
        if (allocated(solution)) deallocate(solution)
        if (allocated(partition_vec)) deallocate(partition_vec)
        
        initialized = .False.
        
    end subroutine finalize_iter_winds_amgx

    subroutine finalize_amgx()
        implicit none

        if (initialized) then
            ! ! Finalize library
            call AMGX_finalize()
        endif

    end subroutine finalize_amgx
    
    !>------------------------------------------------------------
    !! AmgX print callback for verbose error output
    !!------------------------------------------------------------
    subroutine amgx_print_callback(msg, length) bind(C)
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=c_char), intent(in) :: msg(*)
        integer(c_int), value, intent(in) :: length
        integer :: i
        
        ! Print the message character by character
        do i = 1, length
            write(*, '(A)', advance='no') msg(i)
        end do
        
    end subroutine amgx_print_callback

    !>------------------------------------------------------------
    !! Exchange halos for lambda_3d array using MPI
    !!------------------------------------------------------------
    subroutine exchange_lambda_halos(lambda_3d, domain)
        implicit none
        real, dimension(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1:j_e+1), intent(inout) :: lambda_3d
        type(domain_t), intent(in) :: domain
        integer :: ierr
        type(MPI_Status) :: status
        integer :: nz, nx, ny, east_rank, west_rank, north_rank, south_rank
        integer :: tag_ew, tag_we, tag_ns, tag_sn, my_rank
        real, allocatable, dimension(:,:) :: send_buff, rec_buff
        nz = (k_e+1) - (k_s-1) + 1
        nx = (i_e+1) - (i_s-1) + 1
        ny = (j_e+1) - (j_s-1) + 1
        ! Determine neighbor ranks from domain decomposition
        ! For a 2D domain decomposition, neighbors are determined by ximg and yimg

        !get my MPI rank
        call MPI_Comm_rank(domain%compute_comms, my_rank, ierr)
        east_rank = my_rank + 1
        west_rank = my_rank - 1
        north_rank = my_rank + domain%grid%ximages
        south_rank = my_rank - domain%grid%ximages
        
        tag_ew = 100
        tag_we = 101
        tag_ns = 102
        tag_sn = 103
        
        ! Exchange in X direction (East-West)
        if (.not. domain%east_boundary) then
            ! Send east face, receive west halo

            !package east face of lambda_3d
            if (allocated(send_buff)) deallocate(send_buff)
            if (allocated(rec_buff)) deallocate(rec_buff)

            allocate(send_buff(nz,ny))
            allocate(rec_buff(nz,ny))

            send_buff = lambda_3d(i_e, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, &
                             east_rank, tag_ew, &
                             rec_buff, nz*ny, MPI_REAL, &
                             east_rank, tag_we, &
                             domain%compute_comms, status, ierr)
            lambda_3d(i_e+1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
        endif

        if (.not. domain%west_boundary) then
            ! Send west face, receive east halo
            if (allocated(send_buff)) deallocate(send_buff)
            if (allocated(rec_buff)) deallocate(rec_buff)

            allocate(send_buff(nz,ny))
            allocate(rec_buff(nz,ny))
            send_buff = lambda_3d(i_s, k_s-1:k_e+1, j_s-1:j_e+1)
            call MPI_Sendrecv(send_buff, nz*ny, MPI_REAL, &
                             west_rank, tag_we, &
                             rec_buff, nz*ny, MPI_REAL, &
                             west_rank, tag_ew, &
                             domain%compute_comms, status, ierr)
            lambda_3d(i_s-1, k_s-1:k_e+1, j_s-1:j_e+1) = rec_buff
        endif

        ! Exchange in Y direction (North-South)
        if (.not. domain%north_boundary) then
            ! Send north face, receive south halo
            if (allocated(send_buff)) deallocate(send_buff)
            if (allocated(rec_buff)) deallocate(rec_buff)

            allocate(send_buff(nx,nz))
            allocate(rec_buff(nx,nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, &
                             north_rank, tag_ns, &
                             rec_buff, nz*nx, MPI_REAL, &
                             north_rank, tag_sn, &
                             domain%compute_comms, status, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_e+1) = rec_buff
        endif

        if (.not. domain%south_boundary) then
            ! Send south face, receive north halo
            if (allocated(send_buff)) deallocate(send_buff)
            if (allocated(rec_buff)) deallocate(rec_buff)

            allocate(send_buff(nx,nz))
            allocate(rec_buff(nx,nz))
            send_buff = lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s)
            call MPI_Sendrecv(send_buff, nz*nx, MPI_REAL, &
                             south_rank, tag_sn, &
                             rec_buff, nz*nx, MPI_REAL, &
                             south_rank, tag_ns, &
                             domain%compute_comms, status, ierr)
            lambda_3d(i_s-1:i_e+1, k_s-1:k_e+1, j_s-1) = rec_buff
        endif
        
    end subroutine exchange_lambda_halos

end module wind_iterative_amgx
