!>----------------------------------------------------------
!! This module provides the linear wind theory calculations
!!
!!  Code is based primarily off equations in Barstad and Gronas (2006)
!!  @see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
!!
!! The main entry point to the code is:
!!      linear_perturb(domain, options, vsmooth, reverse, useDensity)
!!
!! <pre>
!! Call tree graph :
!!  linear_perturb -> spatial_winds
!!  setup_linwinds -> [ add_buffer_topo, initialize_spatial_winds ]
!!
!! High level routine descriptions / purpose
!!   linear_perturbation      - computes FFT-based linear wind perturbation for a single (U,V,Nsq) state
!!   add_buffer_topo          - generates a topo grid with a surrounding buffer for the FFT
!!   initialize_spatial_winds - generates the look up tables for spatially varying linear winds
!!   spatial_winds            - per-timestep LUT interpolation and perturbation application
!!   setup_linwinds           - sets up module level variables and calls add_buffer_topo
!!   linear_perturb           - main entry point, calls setup on first entry for a given domain
!!
!! Inputs: domain, options, vsmooth, reverse, useDensity
!!      domain,options  = as defined in data_structures
!!      vsmooth         = number of vertical levels to smooth winds over
!!      reverse         = remove linear perturbation instead of adding it
!!      useDensity      = create a linear field that attempts to mitigate the
!!                          boussinesq approx that is embedded in the linear theory
!!                          so it advection can properly incorporate density.
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module linear_theory_winds
    use, intrinsic :: iso_c_binding
    use iso_fortran_env, only: output_unit
    use mpi_f08
    use fft,                        only: fftw_execute_dft,                    &
                                          fftw_plan_dft_2d, fftw_destroy_plan, &
                                          FFTW_FORWARD, FFTW_MEASURE,          &
                                          FFTW_BACKWARD, FFTW_ESTIMATE ! note fft module is defined in fftshift.f90
    use fftshifter,                 only: ifftshift, fftshift
#ifdef _OPENACC
    use openacc, only: acc_get_cuda_stream, acc_async_sync
    use cufft_interface, only: cufftExecZ2Z, CUFFT_INVERSE, cufftPlan2d, CUFFT_Z2Z, cufftSetStream, cufftDestroy
    use fftshifter,                 only: ifftshift2cc_gpu
#endif
    use domain_interface,           only: domain_t
    use string,                     only: str
    use grid_interface,             only: grid_t
    use linear_theory_lut_disk_io,  only: read_LUT, write_LUT
    use mod_atm_utilities,         only : calc_stability, calc_u, calc_v, &
                                          calc_speed, calc_direction,     &
                                          options_t
    use array_utilities,            only: smooth_array, calc_weight, &
                                          linear_space
    use icar_constants,             only: kMAX_FILE_LENGTH, STD_OUT_PE, kREAL, kVARS
    use data_structures,            only: linear_theory_type
    use mod_wrf_constants,          only: piconst

    implicit none

    private
    public :: linear_perturb, setup_linwinds

    logical :: module_initialized = .False.
    logical :: smooth_nsq
    complex(C_DOUBLE_COMPLEX),  allocatable :: terrain_frequency(:,:) ! FFT(terrain)

    ! Module-scope lt_data to prevent stack issues with OpenMP/ifort
    type(linear_theory_type) :: lt_data_m
    !$omp threadprivate(lt_data_m)

    !! Linear wind look up table values and tables
    real, allocatable,          dimension(:)           :: dir_values, nsq_values, spd_values
    ! Look Up Tables for linear perturbation are nspd x n_dir_values x n_nsq_values x nx x nz x ny
    real, allocatable,          dimension(:,:,:,:,:,:) :: hi_u_LUT, hi_v_LUT !, rev_u_LUT, rev_v_LUT

    ! store the linear perturbation so we can update it slightly each time step
    ! this permits the linear field to build up over time.
    real, allocatable, target,  dimension(:,:,:) :: hi_u_perturbation, hi_v_perturbation, lo_u_perturbation, lo_v_perturbation
    real, pointer,              dimension(:,:,:) :: u_perturbation, v_perturbation

    integer :: buffer, original_buffer ! number of grid cells to buffer around the domain MUST be >=1
    integer :: stability_window_size

    real :: max_stability ! limits on the calculated Brunt Vaisala Frequency
    real :: min_stability ! these may need to be a little narrower.
    real :: linear_contribution = 1.0 ! multiplier on uhat,vhat before adding to u,v
    ! controls the rate at which the linearfield updates (should be calculated as f(in_dt))
    ! new/current perturbation is multiplited by linear_update and added to (1-linear_update) * the previous combined perturbation
    real :: linear_update_fraction = 1.0

    real :: dirmax ! =2*pi
    real :: dirmin ! =0
    real :: spdmax ! =30
    real :: spdmin ! =0
    real :: nsqmax ! =log(max_stability)
    real :: nsqmin ! =log(min_stability)
    real :: minimum_layer_size

    integer :: n_dir_values=36
    integer :: n_nsq_values=10
    integer :: n_spd_values=10

    complex,parameter :: imaginary_number = (0,1)

    real, parameter :: SMALL_VALUE = 1e-15

contains

    !>----------------------------------------------------------
    !! Compute linear wind perturbations given background U, V, Nsq
    !!
    !! see Appendix A of Barstad and Gronas (2006) Tellus,58A,2-18
    !!
    !!----------------------------------------------------------
    !>----------------------------------------------------------
    !! Precompute z-independent Fourier-space quantities for a given (U, V, Nsq).
    !! Must be called once before one or more calls to linear_perturbation_apply_z.
    !!----------------------------------------------------------
    subroutine linear_perturbation_setup(U, V, Nsq, lt_data)
        implicit none
        real,                     intent(in)    :: U, V
        real,                     intent(in)    :: Nsq
        type(linear_theory_type), intent(inout) :: lt_data
#ifdef _OPENACC
        integer :: i, j, nx, ny

        nx = size(lt_data%k, 1)
        ny = size(lt_data%k, 2)

        associate(mimag => lt_data%mimag, &
                  m => lt_data%m, &
                  sig => lt_data%sig, &
                  kl => lt_data%kl, &
                  k => lt_data%k, &
                  l => lt_data%l, &
                  msq => lt_data%msq, &
                  denom => lt_data%denom &
        )

        !$acc parallel loop gang vector collapse(2) present(sig, k, l, &
        !$acc   denom, msq, kl, mimag, m)
        do j = 1, ny
            do i = 1, nx
                sig(i,j) = U * k(i,j) + V * l(i,j)
                if (sig(i,j) == 0.0) sig(i,j) = SMALL_VALUE

                denom(i,j) = sig(i,j)**2
                msq(i,j) = Nsq / denom(i,j) * kl(i,j)
                mimag(i,j) = real(sqrt(-msq(i,j))) * imaginary_number

                m(i,j) = sqrt(msq(i,j))
                if (sig(i,j) < 0.0) m(i,j) = m(i,j) * (-1)
                if (real(msq(i,j)) < 0.0) m(i,j) = mimag(i,j)
            end do
        end do

        end associate
#else
        lt_data%sig  = U * lt_data%k + V * lt_data%l
        where(lt_data%sig == 0) lt_data%sig = SMALL_VALUE

        lt_data%denom = lt_data%sig**2

        lt_data%msq   = Nsq / lt_data%denom * lt_data%kl
        lt_data%mimag = (0, 0)
        lt_data%mimag = lt_data%mimag + (real(sqrt(-lt_data%msq)) * imaginary_number)

        lt_data%m = sqrt(lt_data%msq)
        where(lt_data%sig < 0) lt_data%m = lt_data%m * (-1)
        where(real(lt_data%msq) < 0) lt_data%m = lt_data%mimag
#endif

    end subroutine linear_perturbation_setup

    !>----------------------------------------------------------
    !! Apply the z-dependent phase factor and inverse FFT.
    !! Requires linear_perturbation_setup to have been called first for this (U,V,Nsq).
    !!----------------------------------------------------------
    subroutine linear_perturbation_apply_z(z, fourier_terrain, lt_data)
        implicit none
        real,                     intent(in)    :: z
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:)
        type(linear_theory_type), intent(inout) :: lt_data
#ifdef _OPENACC
        integer :: i, j, nx, ny, ierr

        nx = size(fourier_terrain, 1)
        ny = size(fourier_terrain, 2)

        associate(ineta => lt_data%ineta, &
                  m => lt_data%m, &
                  sig => lt_data%sig, &
                  kl => lt_data%kl, &
                  k => lt_data%k, &
                  l => lt_data%l, &
                  uhat => lt_data%uhat, &
                  vhat => lt_data%vhat  &
        )


        ! Element-wise complex arithmetic on GPU
        !$acc parallel loop gang vector collapse(2) present(ineta, m, sig, &
        !$acc   kl, k, l, uhat, vhat, fourier_terrain)
        do j = 1, ny
            do i = 1, nx
                ineta(i,j) = imaginary_number * fourier_terrain(i,j) * exp(imaginary_number * m(i,j) * z)
                ineta(i,j) = ineta(i,j) / (kl(i,j) / ((0 - m(i,j)) * sig(i,j)))
                uhat(i,j) = k(i,j) * ineta(i,j)
                vhat(i,j) = l(i,j) * ineta(i,j)
            end do
        end do

        end associate

        ! GPU ifftshift using pre-allocated workspace
        call ifftshift2cc_gpu(lt_data%uhat, lt_data%fftshift_tmp, nx, ny)
        call ifftshift2cc_gpu(lt_data%vhat, lt_data%fftshift_tmp, nx, ny)

        ! Execute inverse FFT on GPU via cuFFT
        !$acc host_data use_device(lt_data%uhat, lt_data%u_perturb)
        ierr = cufftExecZ2Z(lt_data%cufft_plan, c_loc(lt_data%uhat), c_loc(lt_data%u_perturb), CUFFT_INVERSE)
        !$acc end host_data
        !$acc host_data use_device(lt_data%vhat, lt_data%v_perturb)
        ierr = cufftExecZ2Z(lt_data%cufft_plan, c_loc(lt_data%vhat), c_loc(lt_data%v_perturb), CUFFT_INVERSE)
        !$acc end host_data
#else
        lt_data%ineta = imaginary_number * fourier_terrain * exp(imaginary_number * lt_data%m * z)

        ! uhat = -m * sig * k * ineta / kl,  vhat = -m * sig * l * ineta / kl  (no coriolis)
        lt_data%ineta = lt_data%ineta / (lt_data%kl / ((0 - lt_data%m) * lt_data%sig))
        lt_data%uhat  = lt_data%k * lt_data%ineta
        lt_data%vhat  = lt_data%l * lt_data%ineta

        call ifftshift(lt_data%uhat)
        call ifftshift(lt_data%vhat)
        call fftw_execute_dft(lt_data%uplan, lt_data%uhat, lt_data%u_perturb)
        call fftw_execute_dft(lt_data%vplan, lt_data%vhat, lt_data%v_perturb)
#endif

    end subroutine linear_perturbation_apply_z

    !>----------------------------------------------------------
    !! Convenience wrapper: compute linear wind perturbation at a single height.
    !! Calls setup + apply_z internally.
    !!----------------------------------------------------------
    subroutine linear_perturbation_at_height(U, V, Nsq, z, fourier_terrain, lt_data)
        implicit none
        real,                     intent(in)    :: U, V
        real,                     intent(in)    :: Nsq
        real,                     intent(in)    :: z
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:)
        type(linear_theory_type), intent(inout) :: lt_data

        if ((U==0).and.(V==0)) then
            lt_data%u_perturb = 0
            lt_data%v_perturb = 0
            return
        endif

        call linear_perturbation_setup(U, V, Nsq, lt_data)
        call linear_perturbation_apply_z(z, fourier_terrain, lt_data)

    end subroutine linear_perturbation_at_height


    subroutine linear_perturbation(U, V, Nsq, z_bottom, z_top, minimum_step, fourier_terrain, lt_data, skip_setup)
        implicit none
        real,                     intent(in)    :: U, V             ! U and V components of background wind
        real,                     intent(in)    :: Nsq              ! Brunt-Vaisalla frequency (N squared)
        real,                     intent(in)    :: z_top(:,:), z_bottom(:,:)  ! elevation for the top and bottom bound of a layer
        real,                     intent(in)    :: minimum_step     ! minimum layer step size to compute LT for
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:) ! FFT(terrain)
        type(linear_theory_type), intent(inout) :: lt_data          ! intermediate arrays needed for LT calc
        logical,                  intent(in), optional :: skip_setup ! if .true., assume setup was already called

        integer :: n_steps, i
        real    :: step_size, current_z, start_z, end_z

        ! Handle the trivial case quickly
        if ((U==0).and.(V==0)) then
            lt_data%u_perturb = 0
            lt_data%v_perturb = 0
            return
        endif

        start_z = minval(z_bottom)
        end_z = maxval(z_top)

        ! Fill pre-allocated work arrays with buffer-extended z boundaries
        lt_data%internal_z_top = maxval(z_top)
        lt_data%internal_z_top(buffer:buffer+size(z_top,1)-1, buffer:buffer+size(z_top,2)-1) = z_top

        lt_data%internal_z_bottom = minval(z_bottom)
        lt_data%internal_z_bottom(buffer:buffer+size(z_bottom,1)-1, buffer:buffer+size(z_bottom,2)-1) = z_bottom

        lt_data%layer_count = 0

        step_size = min(minimum_step, minval(z_top - z_bottom))
        if (step_size < 5) then
            write(*,*) "WARNING: Very small vertical step size in linear winds. z_top min: ",minval(z_top)," z_top max: ",maxval(z_top)
            flush(output_unit)
        endif

        current_z = start_z + step_size/2 ! we want the value in the middle of each theoretical layer

        ! Compute z-independent quantities (can be skipped if caller already called setup)
        if (.not.(present(skip_setup).and.skip_setup)) then
            call linear_perturbation_setup(U, V, Nsq, lt_data)
        endif

        ! sum up u/v perturbations over all layers evaluated
        lt_data%u_accumulator = 0
        lt_data%v_accumulator = 0

        do while (current_z < end_z)
            call linear_perturbation_apply_z(current_z, fourier_terrain, lt_data)

            lt_data%layer_fraction = max(0.0,                                                                                                   &
                                min(step_size/2, current_z - lt_data%internal_z_bottom) + min(0.0, lt_data%internal_z_top - current_z)          &
                              + min(step_size/2, lt_data%internal_z_top - current_z)    + min(0.0, current_z - lt_data%internal_z_bottom)       ) / step_size

            lt_data%layer_count = lt_data%layer_count + lt_data%layer_fraction
            lt_data%u_accumulator = lt_data%u_accumulator + lt_data%u_perturb * lt_data%layer_fraction
            lt_data%v_accumulator = lt_data%v_accumulator + lt_data%v_perturb * lt_data%layer_fraction

            current_z = current_z + step_size
        enddo

        lt_data%u_perturb = lt_data%u_accumulator / lt_data%layer_count
        lt_data%v_perturb = lt_data%v_accumulator / lt_data%layer_count
    end subroutine linear_perturbation

#ifdef _OPENACC
    !>----------------------------------------------------------
    !! GPU variant of linear_perturbation that takes precomputed scalars
    !! and buffer-extended z arrays to avoid on-device reductions.
    !! All data must be present on GPU before calling.
    !!----------------------------------------------------------
    subroutine linear_perturbation_gpu(start_z, end_z, step_size, int_z_top, int_z_bot, &
                                       fourier_terrain, lt_data)
        implicit none
        real,                     intent(in)    :: start_z, end_z, step_size
        real,                     intent(in)    :: int_z_top(:,:), int_z_bot(:,:)
        complex(C_DOUBLE_COMPLEX),intent(in)    :: fourier_terrain(:,:)
        type(linear_theory_type), intent(inout) :: lt_data

        real    :: current_z, half_step
        integer :: i, j, nx, ny

        nx = size(fourier_terrain, 1)
        ny = size(fourier_terrain, 2)
        half_step = step_size / 2.0

        associate(u_accumulator => lt_data%u_accumulator, &
                  v_accumulator => lt_data%v_accumulator, &
                  layer_count => lt_data%layer_count, &
                  layer_fraction => lt_data%layer_fraction, &
                  u_perturb => lt_data%u_perturb, &
                  v_perturb => lt_data%v_perturb &
        )

        ! Single present lookup for all arrays — avoids repeated per-component
        ! present-table lookups on each inner parallel loop launch
        !$acc data present(u_accumulator, v_accumulator, layer_count, &
        !$acc   layer_fraction, u_perturb, v_perturb, int_z_top, int_z_bot)

        ! Zero accumulators on GPU
        !$acc parallel loop gang vector collapse(2)
        do j = 1, ny
            do i = 1, nx
                u_accumulator(i,j) = (0.0d0, 0.0d0)
                v_accumulator(i,j) = (0.0d0, 0.0d0)
                layer_count(i,j) = 0.0
            end do
        end do

        current_z = start_z + half_step

        do while (current_z < end_z)
            call linear_perturbation_apply_z(current_z, fourier_terrain, lt_data)

            ! Layer fraction and accumulation on GPU
            !$acc parallel loop gang vector collapse(2)
            do j = 1, ny
                do i = 1, nx
                    layer_fraction(i,j) = max(0.0,                                                     &
                        min(half_step, current_z - int_z_bot(i,j)) + min(0.0, int_z_top(i,j) - current_z)     &
                      + min(half_step, int_z_top(i,j) - current_z) + min(0.0, current_z - int_z_bot(i,j))     &
                    ) / step_size

                    layer_count(i,j) = layer_count(i,j) + layer_fraction(i,j)
                    u_accumulator(i,j) = u_accumulator(i,j) + u_perturb(i,j) * layer_fraction(i,j)
                    v_accumulator(i,j) = v_accumulator(i,j) + v_perturb(i,j) * layer_fraction(i,j)
                end do
            end do

            current_z = current_z + step_size
        enddo

        ! Final normalization on GPU
        !$acc parallel loop gang vector collapse(2)
        do j = 1, ny
            do i = 1, nx
                u_perturb(i,j) = u_accumulator(i,j) / layer_count(i,j)
                v_perturb(i,j) = v_accumulator(i,j) / layer_count(i,j)
            end do
        end do

        !$acc end data
        end associate
    end subroutine linear_perturbation_gpu
#endif

    !>----------------------------------------------------------
    !! Add a smoothed buffer around the edge of the terrain to prevent crazy wrap around effects
    !! in the FFT due to discontinuities between the left and right (or top and bottom) edges of the domain
    !!
    !!----------------------------------------------------------
    subroutine add_buffer_topo(terrain, buffer_topo, smooth_window, buffer, debug)
        implicit none
        real, dimension(:,:), intent(in) :: terrain
        complex(C_DOUBLE_COMPLEX),allocatable,dimension(:,:), intent(inout) :: buffer_topo
        integer, intent(in) :: smooth_window
        integer, intent(in) :: buffer
        logical, intent(in), optional :: debug
        real, dimension(:,:),allocatable :: real_terrain
        integer::nx,ny,i,j,pos, xs,xe,ys,ye, window
        real::weight

        nx=size(terrain,1)+buffer*2
        ny=size(terrain,2)+buffer*2
        allocate(buffer_topo(nx,ny))
        buffer_topo=minval(terrain)

        buffer_topo(1+buffer:nx-buffer,1+buffer:ny-buffer)=terrain
        do i=1,buffer
            weight=i/(real(buffer)*2)
            pos=buffer-i
            buffer_topo(pos+1,1+buffer:ny-buffer)  =terrain(1,1:ny-buffer*2)*(1-weight)+terrain(nx-buffer*2,1:ny-buffer*2)*   weight
            buffer_topo(nx-pos,1+buffer:ny-buffer) =terrain(1,1:ny-buffer*2)*(  weight)+terrain(nx-buffer*2,1:ny-buffer*2)*(1-weight)
        enddo
        do i=1,buffer
            weight=i/(real(buffer)*2)
            pos=buffer-i
            buffer_topo(:,pos+1)  =buffer_topo(:,buffer+1)*(1-weight) + buffer_topo(:,ny-buffer)*   weight
            buffer_topo(:,ny-pos) =buffer_topo(:,buffer+1)*(  weight) + buffer_topo(:,ny-buffer)*(1-weight)
        enddo

        ! smooth the outer most grid cells in all directions to minimize artifacts at the borders of real terrain
        ! smoothing effectively increases as it gets further from the real terrain border (window=min(j,smooth_window))
        if (smooth_window>0) then
            do j=1,buffer
                window=min(j,smooth_window)
                ! smooth the top and bottom borders
                do i=1,nx
                    xs=max(1, i-window)
                    xe=min(nx,i+window)

                    ys=max(1, buffer-j+1-window)
                    ye=min(ny,buffer-j+1+window)

                    buffer_topo(i,buffer-j+1)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

                    ys=max(1, ny-(buffer-j)-window)
                    ye=min(ny,ny-(buffer-j)+window)

                    buffer_topo(i,ny-(buffer-j))=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
                end do
                ! smooth the left and right borders
                do i=1,ny
                    xs=max(1, buffer-j+1-window)
                    xe=min(nx,buffer-j+1+window)
                    ys=max(1, i-window)
                    ye=min(ny,i+window)

                    buffer_topo(buffer-j+1,i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))

                    xs=max(1, nx-(buffer-j)-window)
                    xe=min(nx,nx-(buffer-j)+window)

                    buffer_topo(nx-(buffer-j),i)=sum(buffer_topo(xs:xe,ys:ye))/((xe-xs+1) * (ye-ys+1))
                end do
            end do
        endif

    end subroutine add_buffer_topo

    !>----------------------------------------------------------
    !! Allocate and initialize arrays in lt_data structure
    !!
    !! Initialize constant arrays in lt_data, e.g. k, l, kl wave number arrays, fftw plans
    !!
    !!----------------------------------------------------------
    subroutine initialize_linear_theory_data(lt_data, nx, ny, dx)
        implicit none
        type(linear_theory_type), intent(inout) :: lt_data
        integer, intent(in) :: nx, ny
        real,    intent(in) :: dx

        real :: gain, offset
        integer :: i


        allocate(lt_data%k(nx,ny))
        allocate(lt_data%l(nx,ny))
        allocate(lt_data%kl(nx,ny))
        allocate(lt_data%sig(nx,ny))
        allocate(lt_data%denom(nx,ny))
        allocate(lt_data%m(nx,ny))
        allocate(lt_data%msq(nx,ny))
        allocate(lt_data%mimag(nx,ny))
        allocate(lt_data%ineta(nx,ny))

        ! Compute 2D k and l wavenumber fields
        ! Spacing must match DFT frequencies: dk = 2*pi/(nx*dx), range [-pi/dx, pi/dx)
        offset = piconst / dx
        gain = 2 * offset / nx

        ! compute k wave numbers for the first row
        do i=1, nx
            lt_data%k(i,1) = ((i-1) * gain - offset)
        end do
        ! copy k wave numbers to all other rows
        do i=2, ny
            lt_data%k(:,i) = lt_data%k(:,1)
        enddo

        gain = 2 * offset / ny
        ! compute l wave numbers for the first column
        do i=1,ny
            lt_data%l(1,i) = ((i-1)*gain-offset)
        end do
        ! copy l wave numbers to all other columns
        do i=2,nx
            lt_data%l(i,:) = lt_data%l(1,:)
        enddo

        ! finally compute the kl combination array
        lt_data%kl = lt_data%k**2 + lt_data%l**2
        WHERE (lt_data%kl == 0.0) lt_data%kl = SMALL_VALUE

        allocate(lt_data%uhat(nx,ny))
        allocate(lt_data%u_perturb(nx,ny))
        allocate(lt_data%u_accumulator(nx,ny))
        allocate(lt_data%vhat(nx,ny))
        allocate(lt_data%v_perturb(nx,ny))
        allocate(lt_data%v_accumulator(nx,ny))
        allocate(lt_data%fftshift_tmp(nx,ny))

        ! Pre-allocate work arrays for linear_perturbation (reused across calls)
        allocate(lt_data%layer_count(nx,ny))
        allocate(lt_data%layer_fraction(nx,ny))
        allocate(lt_data%internal_z_top(nx,ny))
        allocate(lt_data%internal_z_bottom(nx,ny))

#ifdef _OPENACC
        ! Create cuFFT plan and bind to the OpenACC default CUDA stream
        ! so cuFFT executions are properly ordered with OpenACC kernels
        i = cufftPlan2d(lt_data%cufft_plan, ny, nx, CUFFT_Z2Z)
        i = cufftSetStream(lt_data%cufft_plan, transfer(acc_get_cuda_stream(acc_async_sync), C_NULL_PTR))

        ! Put all lt_data arrays on the GPU
        !$acc enter data copyin(lt_data%k, lt_data%l, lt_data%kl)
        !$acc enter data create(lt_data%sig, lt_data%denom, lt_data%m, lt_data%msq, lt_data%mimag)
        !$acc enter data create(lt_data%ineta, lt_data%uhat, lt_data%vhat)
        !$acc enter data create(lt_data%u_perturb, lt_data%v_perturb)
        !$acc enter data create(lt_data%u_accumulator, lt_data%v_accumulator)
        !$acc enter data create(lt_data%layer_count, lt_data%layer_fraction)
        !$acc enter data create(lt_data%internal_z_top, lt_data%internal_z_bottom)
        !$acc enter data create(lt_data%fftshift_tmp)
#else
        ! FFTW plan creation (not threadsafe)
        !$omp critical (fftw_lock)
        lt_data%uplan = fftw_plan_dft_2d(ny,nx, lt_data%uhat, lt_data%u_perturb, FFTW_BACKWARD, FFTW_MEASURE)
        lt_data%vplan = fftw_plan_dft_2d(ny,nx, lt_data%vhat, lt_data%v_perturb, FFTW_BACKWARD, FFTW_MEASURE)
        !$omp end critical (fftw_lock)
#endif

    end subroutine initialize_linear_theory_data


    !>----------------------------------------------------------
    !! Deallocate lt_data arrays if allocated, includes fftw_* calls where necessary
    !!
    !!----------------------------------------------------------
    subroutine destroy_linear_theory_data(lt_data)
        implicit none
        type(linear_theory_type), intent(inout) :: lt_data
        integer :: ierr

#ifdef _OPENACC
        ! Remove lt_data arrays from GPU before deallocation
        !$acc exit data delete(lt_data%k, lt_data%l, lt_data%kl)
        !$acc exit data delete(lt_data%sig, lt_data%denom, lt_data%m, lt_data%msq, lt_data%mimag)
        !$acc exit data delete(lt_data%ineta, lt_data%uhat, lt_data%vhat)
        !$acc exit data delete(lt_data%u_perturb, lt_data%v_perturb)
        !$acc exit data delete(lt_data%u_accumulator, lt_data%v_accumulator)
        !$acc exit data delete(lt_data%layer_count, lt_data%layer_fraction)
        !$acc exit data delete(lt_data%internal_z_top, lt_data%internal_z_bottom)
        !$acc exit data delete(lt_data%fftshift_tmp)

        ierr = cufftDestroy(lt_data%cufft_plan)
#else
        ! FFTW plan destruction (not threadsafe)
        !$omp critical (fftw_lock)
        call fftw_destroy_plan(lt_data%uplan)
        call fftw_destroy_plan(lt_data%vplan)
        !$omp end critical (fftw_lock)
#endif

        if (allocated(lt_data%k))       deallocate(lt_data%k)
        if (allocated(lt_data%l))       deallocate(lt_data%l)
        if (allocated(lt_data%kl))      deallocate(lt_data%kl)
        if (allocated(lt_data%sig))     deallocate(lt_data%sig)
        if (allocated(lt_data%denom))   deallocate(lt_data%denom)
        if (allocated(lt_data%m))       deallocate(lt_data%m)
        if (allocated(lt_data%msq))     deallocate(lt_data%msq)
        if (allocated(lt_data%mimag))   deallocate(lt_data%mimag)
        if (allocated(lt_data%ineta))   deallocate(lt_data%ineta)

        if (allocated(lt_data%layer_count))      deallocate(lt_data%layer_count)
        if (allocated(lt_data%layer_fraction))   deallocate(lt_data%layer_fraction)
        if (allocated(lt_data%internal_z_top))   deallocate(lt_data%internal_z_top)
        if (allocated(lt_data%internal_z_bottom))deallocate(lt_data%internal_z_bottom)

        if (allocated(lt_data%uhat))          deallocate(lt_data%uhat)
        if (allocated(lt_data%u_perturb))     deallocate(lt_data%u_perturb)
        if (allocated(lt_data%u_accumulator)) deallocate(lt_data%u_accumulator)
        if (allocated(lt_data%vhat))          deallocate(lt_data%vhat)
        if (allocated(lt_data%v_perturb))     deallocate(lt_data%v_perturb)
        if (allocated(lt_data%v_accumulator)) deallocate(lt_data%v_accumulator)
        if (allocated(lt_data%fftshift_tmp))  deallocate(lt_data%fftshift_tmp)

    end subroutine destroy_linear_theory_data


    subroutine setup_remote_grids(u_grids, v_grids, terrain, nz, halo_size, comms)
        implicit none
        type(grid_t), intent(inout), allocatable :: u_grids(:), v_grids(:)
        real,         intent(in)    :: terrain(:,:)
        integer,      intent(in)    :: nz, halo_size
        type(MPI_Comm), intent(in) :: comms
        
        integer :: nx, ny, i, NUM_COMPUTE

        if (allocated(u_grids)) deallocate(u_grids)
        if (allocated(v_grids)) deallocate(v_grids)

        call MPI_Comm_Size(comms,NUM_COMPUTE)

        allocate(u_grids(NUM_COMPUTE))
        allocate(v_grids(NUM_COMPUTE))

        nx = size(terrain, 1)
        ny = size(terrain, 2)

        do i=1,NUM_COMPUTE
            call u_grids(i)%set_grid_dimensions(nx, ny, nz, image=i, comms=comms, adv_order=halo_size*2, nx_extra=1)
            call v_grids(i)%set_grid_dimensions(nx, ny, nz, image=i, comms=comms, adv_order=halo_size*2, ny_extra=1)
        enddo

    end subroutine setup_remote_grids

    subroutine copy_data_to_remote(wind, grids, LUT_win, i, j, k, z, LUT_nx, LUT_nz, LUT_ny)
        implicit none
        real,           intent(in)  :: wind(:,:)
        type(grid_t),   intent(in)  :: grids(:)
        !real,           intent(inout):: LUT(:,:,:,:,:,:)
        type(MPI_Win),  intent(in)  :: LUT_win
        integer,        intent(in)  :: i, j, k, z, LUT_nx, LUT_nz, LUT_ny

        integer :: img, msg_size, wind_nx, wind_ny, NUM_COMPUTE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
        type(MPI_Datatype) :: LUT_type, send_type
        type(MPI_Group) :: win_group

        msg_size = 1
        disp = 0

        wind_nx = size(wind,1)
        wind_ny = size(wind,2)

        !get group using MPI window
        call MPI_Win_get_group(LUT_win, win_group)
        call MPI_Group_size(win_group, NUM_COMPUTE)

        do img = 1, NUM_COMPUTE
            associate(ims => grids(img)%ims, &
                      jms => grids(img)%jms, &
                      nx  => grids(img)%nx,  &
                      ny  => grids(img)%ny   &
                )

            !$omp critical

            ! Create subarray types describing source and destination memory layouts
            call MPI_Type_create_subarray(6, [n_spd_values, n_dir_values, n_nsq_values, LUT_nx, LUT_nz, LUT_ny], [1, 1, 1, nx, 1, ny], &
                [(k-1),(i-1),(j-1),0,(z-1),0], MPI_ORDER_FORTRAN, MPI_REAL, LUT_type)

            call MPI_Type_create_subarray(2, [wind_nx, wind_ny], [nx, ny], &
                [0,0], MPI_ORDER_FORTRAN, MPI_REAL, send_type)

            call MPI_Type_commit(send_type)
            call MPI_Type_commit(LUT_type)

            ! Passive target: lock target rank, put, unlock — no collective barrier needed
            call MPI_Win_lock(MPI_LOCK_SHARED, img-1, 0, LUT_win)
            call MPI_Put(wind(ims,jms), msg_size, send_type, (img-1), disp, msg_size, LUT_type, LUT_win)
            call MPI_Win_unlock(img-1, LUT_win)

            call MPI_Type_free(send_type)
            call MPI_Type_free(LUT_type)

            !$omp end critical

            end associate
        enddo

    end subroutine copy_data_to_remote

    !>----------------------------------------------------------
    !! Compute look up tables for all combinations of U, V, and Nsq
    !!
    !!----------------------------------------------------------
    subroutine initialize_spatial_winds(domain,options,reverse)
        implicit none
        type(domain_t),  intent(inout)::domain
        type(options_t), intent(in) :: options
        logical, intent(in) :: reverse

        ! local variables used to calculate the LUT
        real :: u,v, layer_height, layer_height_bottom, layer_height_top
        integer :: nx,ny,nz, nxu,nyv, i,j,k,z,ik, n, error
        integer :: fftnx, fftny
        integer, dimension(3,2) :: LUT_dims
        integer :: loops_completed ! this is just used to measure progress in the LUT creation
        integer :: total_LUT_entries, ijk, start_pos, stop_pos
        integer :: ims, jms, my_index, NUM_COMPUTE
        real, allocatable :: temporary_u(:,:,:), temporary_v(:,:,:)
        real, allocatable :: z_above_terrain_bot(:,:,:), z_above_terrain_top(:,:,:)
        integer :: global_nx, global_ny
        character(len=kMAX_FILE_LENGTH) :: LUT_file

        type(c_ptr) :: tmp_ptr
        integer(KIND=MPI_ADDRESS_KIND) :: win_size
        type(MPI_Win) :: U_LUT_win, V_LUT_win
        real, pointer, dimension(:,:,:,:,:,:) :: U_LUT_p, V_LUT_p
        real, allocatable :: start_z_arr(:), end_z_arr(:), step_size_arr(:)
        real, allocatable :: int_z_top(:,:,:), int_z_bot(:,:,:)
        integer :: ii, jj

        type(grid_t), allocatable :: u_grids(:), v_grids(:)

        ! append total number of images and the current image number to the LUT filename
        call MPI_Comm_Size(MPI_COMM_WORLD,NUM_COMPUTE)
        call MPI_Comm_rank(domain%compute_comms, my_index)
        ! MPI returns rank, which is 0-indexed
        my_index = my_index + 1

        LUT_file = trim(options%lt%u_LUT_Filename) // "_" // trim(str(NUM_COMPUTE)) // "_" // trim(str(my_index)) // ".nc"

        ims = lbound(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d,1)
        jms = lbound(domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d,3)

        ! the domain to work over
        nz = size(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,  2)

        fftnx = size(terrain_frequency, 1)
        fftny = size(terrain_frequency, 2)

        call setup_remote_grids(u_grids, v_grids, domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, nz, &
        domain%grid%halo_size, domain%compute_comms)

        ! ensure these are at their required size for all images
        nx = maxval(v_grids%nx)
        ny = maxval(u_grids%ny)
        nxu = maxval(u_grids%nx)
        nyv = maxval(v_grids%ny)

        ! store to make it easy to check dim sizes in read_LUT
        LUT_dims(:,1) = [nxu,nz,ny]
        LUT_dims(:,2) = [nx,nz,nyv]

        ! note:
        ! buffer = (fftnx - nx)/2

        ! default assumes no errors in reading the LUT
        error = 0

        total_LUT_entries = n_dir_values * n_spd_values * n_nsq_values
        
        ! create the array of spd, dir, and nsq values to create LUTs for
        ! generates the values for each look up table dimension
        ! generate table of wind directions to be used
        call linear_space(dir_values,dirmin,dirmax,n_dir_values)
        ! generate table of wind speeds to be used
        call linear_space(nsq_values,nsqmin,nsqmax,n_nsq_values)
        ! generate table of Brunt-Vaisalla frequencies (Nsq) to be used
        call linear_space(spd_values,spdmin,spdmax,n_spd_values)

        ! Allocate the (LARGE) look up tables for both U and V
        if (options%lt%read_LUT .or. module_initialized) then
            if (STD_OUT_PE) write(*,*) "    Reading LUT from file: ", trim(LUT_file)

            error=1
            error = read_LUT(LUT_file, hi_u_LUT, hi_v_LUT, options%domain%dz_levels(:nz), options%lt, LUT_dims)
            
            if (error/=0) then
                if (STD_OUT_PE) write(*,*) "WARNING: LUT on disk does not match that specified in the namelist or does not exist."
                if (STD_OUT_PE) write(*,*) "    LUT will be recreated"
            endif
        endif
        
        start_pos = nint((real(my_index-1) / NUM_COMPUTE) * total_LUT_entries)
        if (my_index==NUM_COMPUTE) then
            stop_pos = total_LUT_entries - 1
        else
            stop_pos  = nint((real(my_index) / NUM_COMPUTE) * total_LUT_entries) - 1
        endif

        ! if (options%general%debug) then
        if (STD_OUT_PE) write(*,*) "Local Look up Table size:", 4*product(shape(hi_u_LUT))/real(2**20), "MB"
        if (STD_OUT_PE) write(*,*) "Wind Speeds:",spd_values
        if (STD_OUT_PE) write(*,*) "Directions:",360*dir_values/(2*piconst)
        if (STD_OUT_PE) write(*,*) "Stabilities:",exp(nsq_values)
        ! endif

        if ( (reverse.or.(.not.((options%lt%read_LUT).and.(error==0)))) .and. (.not.module_initialized) ) then

    
            if (allocated(hi_u_LUT)) deallocate(hi_u_LUT)
            allocate(hi_u_LUT(n_spd_values, n_dir_values, n_nsq_values, nxu, nz, ny), source=0.0)
            if (allocated(hi_v_LUT)) deallocate(hi_v_LUT)
            allocate(hi_v_LUT(n_spd_values, n_dir_values, n_nsq_values, nx,  nz, nyv), source=0.0)

            ! loop over combinations of U, V, and Nsq values
            loops_completed = 0
            if (STD_OUT_PE) write(*,*) "    Initializing linear theory"
            if (STD_OUT_PE) flush(output_unit)

            ! Create MPI windows for one-sided communication of LUT entries across ranks
            win_size = n_spd_values*n_dir_values*n_nsq_values*nxu*nz*ny
            call MPI_WIN_ALLOCATE(win_size*kREAL, kREAL, MPI_INFO_NULL, domain%compute_comms, tmp_ptr, U_LUT_win)
            call C_F_POINTER(tmp_ptr, U_LUT_p, [n_spd_values, n_dir_values, n_nsq_values, nxu, nz, ny])
            win_size = n_spd_values*n_dir_values*n_nsq_values*nx*nz*nyv
            call MPI_WIN_ALLOCATE(win_size*kREAL, kREAL, MPI_INFO_NULL, domain%compute_comms, tmp_ptr, V_LUT_win)
            call C_F_POINTER(tmp_ptr, V_LUT_p, [n_spd_values, n_dir_values, n_nsq_values, nx, nz, nyv])
            U_LUT_p = 1.0
            V_LUT_p = 1.0

            call initialize_linear_theory_data(lt_data_m, fftnx, fftny, domain%dx)
            allocate(temporary_u(fftnx - buffer*2+1, fftny - buffer*2,   nz), source=0.0)
            allocate(temporary_v(fftnx - buffer*2,   fftny - buffer*2+1, nz), source=0.0)

            ! Precompute z_interface - terrain arrays (constant across all dir/spd/nsq combinations)
            global_nx = size(domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, 1)
            global_ny = size(domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, 2)
            allocate(z_above_terrain_bot(global_nx, nz, global_ny))
            allocate(z_above_terrain_top(global_nx, nz, global_ny))
            do z = 1, nz
                z_above_terrain_bot(:,z,:) = domain%vars_3d(domain%var_indx(kVARS%global_z_interface)%v)%data_3d(:,z,:) &
                                           - domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d
                z_above_terrain_top(:,z,:) = z_above_terrain_bot(:,z,:) &
                                           + domain%vars_3d(domain%var_indx(kVARS%global_dz_interface)%v)%data_3d(:,z,:)
            end do

#ifdef _OPENACC
            ! GPU path: precompute per-z-level scalars and buffer-extended z arrays on host,
            ! then run the entire LUT computation on GPU

            allocate(start_z_arr(nz), end_z_arr(nz), step_size_arr(nz))
            allocate(int_z_top(fftnx, fftny, nz), int_z_bot(fftnx, fftny, nz))

            ! Precompute reductions and buffer-extended z arrays on host (avoids GPU reductions)
            do z = 1, nz
                start_z_arr(z) = minval(z_above_terrain_bot(:,z,:))
                end_z_arr(z)   = maxval(z_above_terrain_top(:,z,:))
                step_size_arr(z) = min(minimum_layer_size, minval(z_above_terrain_top(:,z,:) - z_above_terrain_bot(:,z,:)))

                int_z_top(:,:,z) = maxval(z_above_terrain_top(:,z,:))
                int_z_top(buffer:buffer+global_nx-1, buffer:buffer+global_ny-1, z) = z_above_terrain_top(:,z,:)
                int_z_bot(:,:,z) = minval(z_above_terrain_bot(:,z,:))
                int_z_bot(buffer:buffer+global_nx-1, buffer:buffer+global_ny-1, z) = z_above_terrain_bot(:,z,:)
            end do

            ! Put read-only and work arrays on GPU
            !$acc enter data copyin(terrain_frequency, int_z_top, int_z_bot)
            !$acc enter data create(temporary_u, temporary_v)

            do ijk = start_pos, stop_pos
                ik = ijk / n_nsq_values
                j = mod(ijk,n_nsq_values) + 1
                i = ik/n_spd_values + 1
                k = mod(ik,n_spd_values) + 1

                u = calc_u( dir_values(i), spd_values(k) )
                v = calc_v( dir_values(i), spd_values(k) )

                if ((u==0.0).and.(v==0.0)) then
                    loops_completed = loops_completed + nz
                    cycle
                endif

                ! Setup z-independent Fourier quantities on GPU
                call linear_perturbation_setup(u, v, exp(nsq_values(j)), lt_data_m)

                do z=1,nz
                    ! Compute perturbation on GPU using precomputed z arrays
                    call linear_perturbation_gpu(start_z_arr(z), end_z_arr(z), step_size_arr(z), &
                                                    int_z_top(:,:,z), int_z_bot(:,:,z),              &
                                                    terrain_frequency, lt_data_m)

                    ! Extract buffer region into temporary arrays on GPU
                    if (nxu /= nx) then
                        associate(u_perturb => lt_data_m%u_perturb, v_perturb => lt_data_m%v_perturb)
                        !$acc parallel present(temporary_u, u_perturb)
                        !$acc loop gang vector collapse(2)
                        do jj = 1, fftny - buffer*2
                            do ii = 1, fftnx - buffer*2 + 1
                                temporary_u(ii,jj,z) = real( real(                                 &
                                    ( u_perturb(buffer+ii-1, buffer+jj)                  &
                                    + u_perturb(buffer+ii,   buffer+jj) ))) / 2.0
                            end do
                        end do
                        !$acc loop gang vector collapse(2)
                        do jj = 1, fftny - buffer*2 + 1
                            do ii = 1, fftnx - buffer*2
                                temporary_v(ii,jj,z) = real( real(                                 &
                                    ( v_perturb(buffer+ii, buffer+jj-1)                  &
                                    + v_perturb(buffer+ii, buffer+jj) ))) / 2.0
                            end do
                        end do
                        !$acc end parallel
                        end associate
                    else
                        stop "ERROR: linear wind LUT creation not set up for non-staggered grids yet"
                    endif

                    loops_completed = loops_completed+1
                enddo

                ! Copy results to host for MPI distribution
                !$acc update host(temporary_u, temporary_v)
                do z=1,nz
                    call copy_data_to_remote(temporary_u(:,:,z), u_grids, U_LUT_win, i,j,k, z, nxu, nz, ny)
                    call copy_data_to_remote(temporary_v(:,:,z), v_grids, V_LUT_win, i,j,k, z, nx, nz, nyv)
                end do

                if (options%general%interactive) then
                    if (STD_OUT_PE) write(*,"(f5.1,A)") loops_completed/real(nz*(stop_pos-start_pos+1))*100," %"
                    if (STD_OUT_PE) flush(output_unit)
                endif
            end do

            !$acc exit data delete(terrain_frequency, int_z_top, int_z_bot)
            !$acc exit data delete(temporary_u, temporary_v)
            deallocate(start_z_arr, end_z_arr, step_size_arr, int_z_top, int_z_bot)
#else
            do ijk = start_pos, stop_pos
                ! Decode (dir, nsq, spd) indices from combined ijk index
                ik = ijk / n_nsq_values
                j = mod(ijk,n_nsq_values) + 1
                i = ik/n_spd_values + 1
                k = mod(ik,n_spd_values) + 1

                u = calc_u( dir_values(i), spd_values(k) )
                v = calc_v( dir_values(i), spd_values(k) )

                ! Setup z-independent Fourier quantities once for this (dir, spd, nsq)
                call linear_perturbation_setup(u, v, exp(nsq_values(j)), lt_data_m)

                ! Compute all z-levels, accumulate into 3D temporary arrays
                do z=1,nz
                    call linear_perturbation(u, v, exp(nsq_values(j)),        &
                                             z_above_terrain_bot(:,z,:),       &
                                             z_above_terrain_top(:,z,:),       &
                                             minimum_layer_size, terrain_frequency, lt_data_m, skip_setup=.true.)

                    if (nxu /= nx) then
                        temporary_u(:,:,z) = real( real(                                              &
                                ( lt_data_m%u_perturb(buffer:fftnx-buffer,     1+buffer:fftny-buffer)           &
                                + lt_data_m%u_perturb(1+buffer:fftnx-buffer+1,     1+buffer:fftny-buffer)) )) / 2

                        temporary_v(:,:,z) = real( real(                                              &
                                ( lt_data_m%v_perturb(1+buffer:fftnx-buffer,     buffer:fftny-buffer)         &
                                + lt_data_m%v_perturb(1+buffer:fftnx-buffer,     1+buffer:fftny-buffer+1)) )) / 2
                    else
                        stop "ERROR: linear wind LUT creation not set up for non-staggered grids yet"
                    endif

                    loops_completed = loops_completed+1
                enddo

                ! Distribute this entry to all ranks via passive target one-sided MPI
                ! (lock/unlock inside copy_data_to_remote — no collective barriers)
                do z=1,nz
                    call copy_data_to_remote(temporary_u(:,:,z), u_grids, U_LUT_win, i,j,k, z, nxu, nz, ny)
                    call copy_data_to_remote(temporary_v(:,:,z), v_grids, V_LUT_win, i,j,k, z, nx, nz, nyv)
                end do

                if (options%general%interactive) then
                    if (STD_OUT_PE) write(*,"(f5.1,A)") loops_completed/real(nz*(stop_pos-start_pos+1))*100," %"
                    if (STD_OUT_PE) flush(output_unit)
                endif
            end do
#endif

            deallocate(z_above_terrain_bot, z_above_terrain_top)

            ! Ensure all ranks have completed their puts before reading from the windows
            call MPI_Barrier(domain%compute_comms)
            call MPI_Win_sync(U_LUT_win)
            call MPI_Win_sync(V_LUT_win)
            hi_u_LUT = U_LUT_p(:,:,:,:,:,:)
            hi_v_LUT = V_LUT_p(:,:,:,:,:,:)
            call MPI_Win_free(U_LUT_win)
            call MPI_Win_free(V_LUT_win)

            ! memory needs to be freed so this structure can be used again when removing linear winds
            call destroy_linear_theory_data(lt_data_m)

            if (STD_OUT_PE) write(*,*) "All images: 100 % Complete"
            if (STD_OUT_PE) write(*,*) char(10),"--------  Linear wind look up table generation complete ---------"
            if (STD_OUT_PE) flush(output_unit)
        endif

        if ((options%lt%write_LUT).and.(.not.reverse).and.(.not.module_initialized)) then
            if ((options%lt%read_LUT) .and. (error == 0)) then
                if (STD_OUT_PE) write(*,*) "    Not writing Linear Theory LUT to file because LUT was read from file"
            else
                if (STD_OUT_PE) write(*,*) "    Writing Linear Theory LUT to file: ", trim(LUT_file)
                error = write_LUT(LUT_file, hi_u_LUT, hi_v_LUT, options%domain%dz_levels(:nz), options%lt)
            endif
        endif

    end subroutine initialize_spatial_winds


    !>----------------------------------------------------------
    !! Compute a spatially variable linear wind perturbation
    !! based off of look uptables computed in via setup
    !! for each grid point, find the closest LUT data in U and V space
    !! then bilinearly interpolate the nearest LUT values for that points linear wind field
    !!
    !!----------------------------------------------------------
    subroutine spatial_winds(domain,reverse, vsmooth, winsz, update)
        implicit none
        type(domain_t), intent(inout):: domain
        logical,        intent(in)   :: reverse
        integer,        intent(in)   :: vsmooth
        integer,        intent(in)   :: winsz
        logical,        intent(in)   :: update

        integer :: nx,nxu, ny,nyv, nz, i,j,k, smoothz
        integer :: uk, vi
        integer :: step, dpos, npos, spos, nexts, nextd, nextn
        integer :: south, east, west, top, bottom
        real :: u, v
        real,allocatable :: u1d(:),v1d(:)
        integer :: ims, ime, jms, jme, ims_u, ime_u, jms_v, jme_v, kms, kme
        real :: dweight, nweight, sweight, curspd, curdir, curnsq, wind_first, wind_second
        real :: hydrometeors

        associate(                                                                                          &
            u_data    => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,                                &
            v_data    => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d,                                &
            u_dqdt    => domain%vars_3d(domain%var_indx(kVARS%u)%v)%dqdt_3d,                                &
            v_dqdt    => domain%vars_3d(domain%var_indx(kVARS%v)%v)%dqdt_3d,                                &
            theta     => domain%vars_3d(domain%var_indx(kVARS%potential_temperature)%v)%data_3d,             &
            exner     => domain%vars_3d(domain%var_indx(kVARS%exner)%v)%data_3d,                            &
            z_data    => domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d,                                &
            qv        => domain%vars_3d(domain%var_indx(kVARS%water_vapor)%v)%data_3d,                      &
            nsquared  => domain%vars_3d(domain%var_indx(kVARS%nsquared)%v)%data_3d                          &
        )

        ! Copy only the fields we need from device to host
        !$acc update host(u_data, v_data, theta, exner, z_data, qv, nsquared)
        if (update) then
            !$acc update host(u_dqdt, v_dqdt)
        endif
        if (domain%var_indx(kVARS%cloud_water_mass)%v > 0) then
            !$acc update host(domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d)
        endif
        if (domain%var_indx(kVARS%ice_mass)%v > 0) then
            !$acc update host(domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d)
        endif
        if (domain%var_indx(kVARS%rain_mass)%v > 0) then
            !$acc update host(domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d)
        endif
        if (domain%var_indx(kVARS%snow_mass)%v > 0) then
            !$acc update host(domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d)
        endif

        ims_u = lbound(u_data,1)
        ime_u = ubound(u_data,1)
        jms = lbound(u_data,3)
        jme = ubound(u_data,3)

        ims = lbound(v_data,1)
        ime = ubound(v_data,1)
        jms_v = lbound(v_data,3)
        jme_v = ubound(v_data,3)

        kms = lbound(u_data,2)
        kme = ubound(u_data,2)

        nx  = size(domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,1)
        ny  = size(domain%vars_2d(domain%var_indx(kVARS%latitude)%v)%data_2d,2)
        nz  = size(u_data,2)
        nxu = size(u_data,1)
        nyv = size(v_data,3)

        if (reverse) then
            ! TODO: lo_*_perturbation arrays are allocated but no low-res LUT is generated yet
            u_perturbation=>lo_u_perturbation
            v_perturbation=>lo_v_perturbation
        else
            u_perturbation=>hi_u_perturbation
            v_perturbation=>hi_v_perturbation
        endif

        ! ---- Compute Brunt-Vaisala frequency (N²) ----
        do k=1,ny
            do j=1,nz
                do i=1,nx
                    top = min(j+vsmooth, nz)
                    bottom = max(1, j - (vsmooth - (top-j)))

                    hydrometeors = 0
                    if (domain%var_indx(kVARS%cloud_water_mass)%v > 0) &
                        hydrometeors = hydrometeors + domain%vars_3d(domain%var_indx(kVARS%cloud_water_mass)%v)%data_3d(i+ims-1,j+kms-1,k+jms-1)
                    if (domain%var_indx(kVARS%ice_mass)%v > 0) &
                        hydrometeors = hydrometeors + domain%vars_3d(domain%var_indx(kVARS%ice_mass)%v)%data_3d(i+ims-1,j+kms-1,k+jms-1)
                    if (domain%var_indx(kVARS%rain_mass)%v > 0) &
                        hydrometeors = hydrometeors + domain%vars_3d(domain%var_indx(kVARS%rain_mass)%v)%data_3d(i+ims-1,j+kms-1,k+jms-1)
                    if (domain%var_indx(kVARS%snow_mass)%v > 0) &
                        hydrometeors = hydrometeors + domain%vars_3d(domain%var_indx(kVARS%snow_mass)%v)%data_3d(i+ims-1,j+kms-1,k+jms-1)

                    if (.not.reverse) then
                        nsquared(i+ims-1,j+kms-1,k+jms-1) = calc_stability(  &
                            theta(i+ims-1,bottom+kms-1,k+jms-1), theta(i+ims-1,top+kms-1,k+jms-1),      &
                            exner(i+ims-1,bottom+kms-1,k+jms-1), exner(i+ims-1,top+kms-1,k+jms-1),      &
                            z_data(i+ims-1,bottom+kms-1,k+jms-1), z_data(i+ims-1,top+kms-1,k+jms-1),    &
                            qv(i+ims-1,bottom+kms-1,k+jms-1), qv(i+ims-1,top+kms-1,k+jms-1),            &
                            hydrometeors)

                        nsquared(i+ims-1,j+kms-1,k+jms-1) = max(min_stability, min(max_stability, &
                                                nsquared(i+ims-1,j+kms-1,k+jms-1)))
                    else
                        nsquared(i+ims-1,j+kms-1,k+jms-1) = 3e-6
                    endif
                end do
                ! LUT is computed in log space
                nsquared(:,j+kms-1,k+jms-1) = log(nsquared(:,j+kms-1,k+jms-1))
            end do

            if (smooth_nsq) then
                do j=1,nz
                    top = min(j+vsmooth,nz)
                    bottom = max(1, j - (vsmooth - (top-j)) )

                    do smoothz = bottom, j-1
                        nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1) + nsquared(:,smoothz+kms-1,k+jms-1)
                    end do
                    do smoothz = j+1, top
                        nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1) + nsquared(:,smoothz+kms-1,k+jms-1)
                    end do
                    nsquared(:,j+kms-1,k+jms-1) = nsquared(:,j+kms-1,k+jms-1)/(top-bottom+1)
                end do
            endif
        end do

        if (smooth_nsq) then
            call smooth_array(nsquared, winsz, ydim=3)
        endif

        ! ---- Interpolate LUT and apply perturbations ----
        allocate(u1d(nxu), v1d(nxu))
        do k=1, nyv

            uk = min(k,ny)
            do i=1,nxu
                vi = min(i,nx)
                u1d(i) = sum(u_data(i+ims_u-1,:,uk+jms-1)) / nz
                v1d(i) = sum(v_data(vi+ims-1, :,k+jms_v-1)) / nz
            enddo

            do j=1, nz
                do i=1, nxu
                        uk = min(k,ny)
                        vi = min(i,nx)

                        bottom= max(j - winsz, 1)
                        top   = min(j + winsz,nz)

                        if (reverse) then
                            u = u_data(i+ims_u-1,j,uk+jms-1)
                            v = v_data(vi+ims-1, j,k+jms_v-1)
                        else
                            u = u1d(i)
                            v = v1d(i)
                        endif

                        ! Find LUT position for wind direction
                        dpos = 1
                        curdir = calc_direction( u, v )
                        do step=1, n_dir_values
                            if (curdir > dir_values(step)) dpos = step
                        end do

                        ! Find LUT position for wind speed
                        spos = 1
                        curspd = calc_speed( u, v )
                        do step=1, n_spd_values
                            if (curspd > spd_values(step)) spos = step
                        end do

                        ! Find LUT position for N² (mean in log space)
                        curnsq = sum(nsquared(vi+ims-1,bottom+kms-1:top+kms-1,uk+jms-1)) / (top - bottom + 1)
                        npos = 1
                        do step=1, n_nsq_values
                            if (curnsq > nsq_values(step)) npos = step
                        end do

                        ! Interpolation weights
                        dweight = calc_weight(dir_values, dpos, nextd, curdir)
                        sweight = calc_weight(spd_values, spos, nexts, curspd)
                        nweight = calc_weight(nsq_values, npos, nextn, curnsq)

                        ! Trilinear interpolation of LUT and apply perturbation to u
                        if (k<=ny) then
                            wind_first =      nweight  * (dweight * hi_u_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * hi_u_LUT(spos, nextd,npos, i,j,k))   &
                                        +  (1-nweight) * (dweight * hi_u_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * hi_u_LUT(spos, nextd,nextn,i,j,k))

                            wind_second=      nweight  * (dweight * hi_u_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * hi_u_LUT(nexts,nextd,npos, i,j,k))   &
                                        +  (1-nweight) * (dweight * hi_u_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * hi_u_LUT(nexts,nextd,nextn,i,j,k))

                            u_perturbation(i,j,k) = u_perturbation(i,j,k) * (1-linear_update_fraction) &
                                        + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                            if (update) then
                                u_dqdt(i+ims_u-1,j,k+jms-1) = u_dqdt(i+ims_u-1,j,k+jms-1) + u_perturbation(i,j,k) * linear_contribution
                            else
                                u_data(i+ims_u-1,j,k+jms-1) = u_data(i+ims_u-1,j,k+jms-1) + u_perturbation(i,j,k) * linear_contribution
                            endif
                        endif

                        ! Trilinear interpolation of LUT and apply perturbation to v
                        if (i<=nx) then
                            wind_first =      nweight  * (dweight * hi_v_LUT(spos, dpos,npos, i,j,k) + (1-dweight) * hi_v_LUT(spos, nextd,npos, i,j,k))    &
                                        +  (1-nweight) * (dweight * hi_v_LUT(spos, dpos,nextn,i,j,k) + (1-dweight) * hi_v_LUT(spos, nextd,nextn,i,j,k))

                            wind_second=      nweight  * (dweight * hi_v_LUT(nexts,dpos,npos, i,j,k) + (1-dweight) * hi_v_LUT(nexts,nextd,npos, i,j,k))    &
                                        +  (1-nweight) * (dweight * hi_v_LUT(nexts,dpos,nextn,i,j,k) + (1-dweight) * hi_v_LUT(nexts,nextd,nextn,i,j,k))

                            v_perturbation(i,j,k) = v_perturbation(i,j,k) * (1-linear_update_fraction) &
                                        + linear_update_fraction * (sweight*wind_first + (1-sweight)*wind_second)

                            if (update) then
                                v_dqdt(i+ims-1,j,k+jms_v-1) = v_dqdt(i+ims-1,j,k+jms_v-1) + v_perturbation(i,j,k) * linear_contribution
                            else
                                v_data(i+ims-1,j,k+jms_v-1) = v_data(i+ims-1,j,k+jms_v-1) + v_perturbation(i,j,k) * linear_contribution
                            endif
                        endif
                end do
            end do
        end do
        deallocate(u1d, v1d)

        ! Convert nsquared back from log-space to linear (LUT interpolation works in log-space)
        nsquared = exp(nsquared)

        ! Copy modified fields back to device
        if (update) then
            !$acc update device(u_dqdt, v_dqdt)
        else
            !$acc update device(u_data, v_data)
        endif
        !$acc update device(nsquared)

        end associate

    end subroutine spatial_winds


    !>----------------------------------------------------------
    !! Copy options from the options data structure into module level variables
    !!
    !!----------------------------------------------------------
    subroutine set_module_options(options)
        implicit none
        type(options_t), intent(in) :: options

        original_buffer       = options%lt%buffer
        smooth_nsq            = options%lt%smooth_nsq

        stability_window_size = options%lt%stability_window_size
        max_stability         = options%lt%max_stability
        min_stability         = options%lt%min_stability
        linear_contribution   = options%lt%linear_contribution
        linear_update_fraction = options%lt%linear_update_fraction

        dirmax = options%lt%dirmax
        dirmin = options%lt%dirmin
        spdmax = options%lt%spdmax
        spdmin = options%lt%spdmin
        nsqmax = options%lt%nsqmax
        nsqmin = options%lt%nsqmin
        n_dir_values = options%lt%n_dir_values
        n_nsq_values = options%lt%n_nsq_values
        n_spd_values = options%lt%n_spd_values
        minimum_layer_size = options%lt%minimum_layer_size

    end subroutine set_module_options

    !>----------------------------------------------------------
    !! Called from linear_perturb the first time perturb is called
    !! compute FFT(terrain), and dzdx,dzdy components
    !!
    !!----------------------------------------------------------
    subroutine setup_linwinds(domain, options, reverse, useDensity)
        implicit none
        type(domain_t),     intent(inout)  :: domain
        type(options_t),    intent(in)     :: options
        logical,            intent(in)     :: reverse, useDensity
        ! locals
        complex(C_DOUBLE_COMPLEX), allocatable  :: complex_terrain_firstpass(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable  :: complex_terrain(:,:)
        type(C_PTR) :: plan
        integer     :: nx, ny, nz

        ! store module level variables so we don't have to pass options through everytime
        ! lots of these things probably need to be moved to the linearizable class so they
        ! can be separated for the domain and bc fields
        ! this is a little tricky, because we don't want to have to calculate the LUTs
        ! twice, once for domain and once for bc%next_domain
        call set_module_options(options)

        if (STD_OUT_PE) write(*,*) "Initializing linear winds"
        if (allocated(terrain_frequency)) deallocate(terrain_frequency)

        ! Create a buffer zone around the topography to smooth the edges
        buffer = original_buffer
        ! first create it including a 5 grid cell smoothing function
        call add_buffer_topo(domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, complex_terrain_firstpass, 5, buffer)
        buffer = 2
        ! then further add a small (~2) grid cell buffer where all cells have the same value
        call add_buffer_topo(real(real(complex_terrain_firstpass)), complex_terrain, 0, buffer, debug=options%general%debug)
        buffer = buffer + original_buffer

        nx = size(complex_terrain, 1)
        ny = size(complex_terrain, 2)

        allocate(terrain_frequency(nx,ny))

        ! calculate the fourier transform of the terrain for use in linear winds
        plan = fftw_plan_dft_2d(ny, nx, complex_terrain, terrain_frequency, FFTW_FORWARD, FFTW_ESTIMATE)
        call fftw_execute_dft(plan, complex_terrain, terrain_frequency)
        call fftw_destroy_plan(plan)
        ! normalize FFT by N - grid cells
        terrain_frequency = terrain_frequency / (nx * ny)
        ! shift the grid cell quadrants
        ! need to test what effect all of the related shifts actually have...
        call fftshift(terrain_frequency)

        if (linear_contribution/=1) then
            if (STD_OUT_PE) write(*,*) "  Using a fraction of the linear perturbation:",linear_contribution
        endif

        nx = size(domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, 1)
        nz = size(domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d,       2)
        ny = size(domain%vars_2d(domain%var_indx(kVARS%global_terrain)%v)%data_2d, 2)

        ! allocate the fields that will hold the perturbation only so we can update it
        ! slowly and add the total to the domain%u,v
        if (reverse) then
            if (allocated(lo_u_perturbation)) deallocate(lo_u_perturbation)
            if (allocated(lo_v_perturbation)) deallocate(lo_v_perturbation)
            allocate(lo_u_perturbation(nx+1,nz,ny))
            lo_u_perturbation=0
            allocate(lo_v_perturbation(nx,nz,ny+1))
            lo_v_perturbation=0
            u_perturbation=>lo_u_perturbation
            v_perturbation=>lo_v_perturbation
        else
            if (allocated(hi_u_perturbation)) deallocate(hi_u_perturbation)
            if (allocated(hi_v_perturbation)) deallocate(hi_v_perturbation)
            allocate(hi_u_perturbation(nx+1,nz,ny))
            hi_u_perturbation=0
            allocate(hi_v_perturbation(nx,nz,ny+1))
            hi_v_perturbation=0
            u_perturbation=>hi_u_perturbation
            v_perturbation=>hi_v_perturbation
        endif

        if (STD_OUT_PE) write(*,*) "  Generating a spatially variable linear perturbation look up table"
        call initialize_spatial_winds(domain, options, reverse)

        module_initialized = .True.

    end subroutine setup_linwinds

    !>----------------------------------------------------------
    !! Initialize and/or apply linear wind solution.
    !!
    !! Called from ICAR to update the U and V wind fields based on linear theory (W is calculated to balance U/V)
    !!
    !!----------------------------------------------------------
    subroutine linear_perturb(domain,options,vsmooth,reverse,useDensity, update)
        implicit none
        type(domain_t),     intent(inout):: domain
        type(options_t),    intent(in)   :: options
        integer,            intent(in)   :: vsmooth
        logical,            intent(in),  optional :: reverse, useDensity, update

        logical :: rev, useD, updt

        rev = .False.
        if (present(reverse)) rev = reverse

        useD = .False.
        if (present(useDensity)) useD = useDensity

        updt = .False.
        if (present(update)) updt = update

        if (rev) then
            linear_contribution = options%lt%rm_linear_contribution
        else
            linear_contribution = options%lt%linear_contribution
        endif

        call spatial_winds(domain, rev, vsmooth, stability_window_size, update=updt)

    end subroutine linear_perturb
end module linear_theory_winds
