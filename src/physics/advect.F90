!> ----------------------------------------------------------------------------
!!  Standard advection scheme with variable order
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_std
    use mpi_f08, only : MPI_Wtime
    use icar_constants
    use options_interface, only: options_t
    use domain_interface,  only: domain_t
    use adv_fluxcorr,      only: WRF_flux_corr
    use timer_interface,   only: timer_t
    implicit none
    private

    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte
    integer :: i_s_w, i_e_w, j_s_w, j_e_w, i_s, i_e, j_s, j_e, horder, vorder
    !type(timer_t) :: flux_time, flux_up_time, flux_corr_time, sum_time
    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v
    real, dimension(:,:,:), allocatable   :: flux_x, flux_y, flux_z

    ! Fine-mesh flux arrays (used by flux3_fm / sum_kernel_fm)
    real, dimension(:,:,:), allocatable   :: flux_x_fm, flux_y_fm, flux_z_fm

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind, adv_std_clean_wind_arrays
    public :: adv_std_compute_wind_3d_fm, flux3_fm, sum_kernel_fm, adv_std_clean_wind_arrays_fm
    public :: flux_x_fm, flux_y_fm, flux_z_fm

contains

    subroutine adv_std_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options


        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        its = domain%its
        ite = domain%ite
        jts = domain%jts
        jte = domain%jte
        
        !set order of advection
        horder = options%adv%h_order
        vorder = options%adv%v_order
        
        !Define bounds of advection computation. If using monotonic flux-limiter, it is necesarry to increase
        !advection bounds by 1. The necesarry extension of the halo is handeled in domain_object
        i_s = its
        i_e = ite
        j_s = jts
        j_e = jte

        i_s_w = i_s 
        i_e_w = i_e 
        j_s_w = j_s 
        j_e_w = j_e 

        if (options%adv%flux_corr==kFLUXCOR_MONO) then
            i_s = its - 1
            i_e = ite + 1
            j_s = jts - 1
            j_e = jte + 1
            
            ! wind arrays need to be extended by 1 in each direction for monotonic flux correction
            ! this allows for 2 upwind advection steps to be computed without needing a halo exchange
            i_s_w = i_s - 1
            i_e_w = i_e + 1
            j_s_w = j_s - 1
            j_e_w = j_e + 1
        endif
        
    end subroutine

    subroutine adv_std_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for adv4 advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

        ! List the variables that are required for restarts with adv4 advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

    end subroutine

    subroutine flux3(q,U_m,V_m,W_m,flux_x,flux_z,flux_y,t_factor)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1), intent(in)       :: U_m, V_m, W_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1),   intent(out)    :: flux_x
        real, dimension(i_s:i_e+1,  kms:kme,j_s:j_e+1), intent(out)    :: flux_y
        real, dimension(i_s:i_e+1,  kms:kme+1,j_s:j_e+1), intent(out)    :: flux_z
        real, intent(in) :: t_factor
        integer :: i, j, k
        real :: tmp, coef, u, v, w, q0, q1, q2, qn1, qn2, qn3, abs_u, t_factor_compact
        real :: qin1, qin2, qi1, qk1, qj1, qkn1, qkn2, qjn1, qjn2
        integer :: top, bot, bot2, nor, sou, sou2, eas, wes, wes2
        ! GPU optimization variables
        real :: q_cache(7), u_val, v_val, w_val, abs_u_val
                      
        !$acc data present(q,U_m,V_m,W_m,flux_x,flux_y,flux_z)

        if (horder==1) then
            t_factor_compact = 0.5  * t_factor
            !$acc parallel async(1)
            ! Using tile(32,2) for better cache locality in 3D loop with moderate k dimension
            !$acc loop gang vector tile(32,2,1) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
            do j = j_s, j_e
                do k = kms+1, kme
                    do i = i_s, i_e
                        ! Cache frequently accessed values to reduce memory traffic
                        q0 = q(i,k,j)
                        qn1 = q(i-1,k,j)
                        qn3 = q(i,k,j-1)

                        ! X-direction flux computation with cached values
                        u_val = U_m(i,k,j)
                        abs_u_val = ABS(u_val)
                        flux_x(i, k, j) = (u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0
                        flux_x(i, k, j) = flux_x(i, k, j) * t_factor_compact

                        ! Y-direction flux computation with cached values
                        v_val = V_m(i,k,j)
                        abs_u_val = ABS(v_val)  ! reuse abs_u_val variable
                        flux_y(i, k, j) = (v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0
                        flux_y(i, k, j) = flux_y(i, k, j) * t_factor_compact
                    enddo
                enddo
            enddo

            ! Using tile(32) for 2D loops - good balance for horizontal dimensions
            !$acc loop gang vector tile(32,2) private(u_val, abs_u_val, q0, qn1)
            do j = j_s, j_e
                do k = kms, kme
                    qn1 = q(i_e,k,j)
                    q0 = q(i_e+1,k,j)

                    u_val = U_m(i_e+1,k,j)
                    abs_u_val = ABS(u_val)
                    flux_x(i_e+1, k, j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector tile(32,2) private(v_val, abs_u_val, q0, qn3)
            do k = kms, kme
                do i = i_s, i_e
                    qn3 = q(i,k,j_e)
                    q0 = q(i,k,j_e+1)

                    v_val = V_m(i,k,j_e+1)
                    abs_u_val = ABS(v_val)  ! reuse variable
                    flux_y(i, k, j_e+1) = ((v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector tile(64,1) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
            do j = j_s, j_e
                do i = i_s, i_e
                    q0 = q(i,kms,j)
                    qn1 = q(i-1,kms,j)
                    qn3 = q(i,kms,j-1)

                    u_val = U_m(i,kms,j)
                    abs_u_val = ABS(u_val)
                    flux_x(i, kms, j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact

                    v_val = V_m(i,kms,j)
                    abs_u_val = ABS(v_val)  ! reuse variable
                    flux_y(i, kms, j) = ((v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo
            !$acc end parallel

        else if (horder==3) then
            coef = (1./12) * t_factor

            ! Using tile(32,2) for 3D loop with complex stencil operations
            !$acc parallel loop gang vector tile(32,2,1) async(1) private(u_val, abs_u_val, q_cache, tmp)
            do j = j_s,j_e+1
            do k = kms,kme
            do i = i_s,i_e+1
                ! Cache q values to reduce memory accesses
                q_cache(1) = q(i,k,j-2)  ! qn2
                q_cache(2) = q(i,k,j-1)  ! qn1
                q_cache(3) = q(i-2,k,j)  ! qn2
                q_cache(4) = q(i-1,k,j)  ! qn1
                q_cache(5) = q(i,k,j)    ! q0
                q_cache(6) = q(i+1,k,j)  ! q1
                q_cache(7) = q(i,k,j+1)  ! q1

                u_val = U_m(i,k,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 4th order flux calculation using cached values
                tmp = 7.0 * (q_cache(5) + q_cache(4)) - (q_cache(6) + q_cache(3))
                tmp = u_val * tmp
                
                ! Optimized 3rd order diffusive terms using cached values
                tmp = tmp - abs_u_val * (3.0 * (q_cache(5) - q_cache(4)) - (q_cache(6) - q_cache(3)))
                
                ! Optimized upwind flux calculation using cached values
                flux_x(i,k,j) = tmp * coef
                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 4th order flux calculation using cached values
                tmp = 7.0 * (q_cache(5) + q_cache(2)) - (q_cache(7) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 3rd order diffusive terms using cached values
                tmp = tmp - abs_u_val * (3.0 * (q_cache(5) - q_cache(2)) - (q_cache(7) - q_cache(1)))

                ! Optimized upwind flux calculation using cached values
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo

        else if (horder==5) then
            coef = (1./60)*t_factor
            !$acc parallel async(1)
            ! Using tile(32,2) for complex 5th order stencil computations
            !$acc loop gang vector tile(32,2,1)  private(u_val, abs_u_val, q_cache, tmp)
            do j = j_s,j_e
            do k = kms,kme
            do i = i_s,i_e+1
                ! Cache all required q values in registers to minimize memory accesses
                q_cache(1) = q(i-3,k,j)  ! qn3
                q_cache(2) = q(i-2,k,j)  ! qn2
                q_cache(3) = q(i-1,k,j)  ! qn1
                q_cache(4) = q(i,k,j)    ! q0
                q_cache(5) = q(i+1,k,j)  ! q1
                ! Note: q2 = q(i+2,k,j) accessed directly when needed
                
                u_val = U_m(i,k,j)
                abs_u_val = ABS(u_val)

                ! Optimized 6th order flux calculation using cached values
                ! Original: tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = 37.0 * (q_cache(4) + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (q(i+2,k,j) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 5th order diffusive terms using cached values
                ! Original: tmp = tmp - abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                tmp = tmp - abs_u_val * (10.0 * (q_cache(4) - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (q(i+2,k,j) - q_cache(1)))
                flux_x(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            
            !$acc loop gang vector tile(32,2,1)  private(u_val, abs_u_val, q_cache, tmp)
            do j = j_s,j_e+1
            do k = kms,kme
            do i = i_s,i_e
                ! Cache all required q values in registers to minimize memory accesses
                q_cache(1) = q(i,k,j-3)  ! qn3
                q_cache(2) = q(i,k,j-2)  ! qn2
                q_cache(3) = q(i,k,j-1)  ! qn1
                q_cache(4) = q(i,k,j)    ! q0
                q_cache(5) = q(i,k,j+1)  ! q1
                ! Note: q2 = q(i,k,j+2) accessed directly when needed
                
                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 6th order flux calculation using cached values
                tmp = 37.0 * (q_cache(4) + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (q(i,k,j+2) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 5th order diffusive terms using cached values
                tmp = tmp - abs_u_val * (10.0 * (q_cache(4) - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (q(i,k,j+2) - q_cache(1)))
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            !$acc end parallel
        endif
        
        if (vorder==1) then
            t_factor_compact = 0.5  * t_factor
            !$acc parallel async(2)
            ! Using tile(32,2) for vertical flux computation
            !$acc loop gang vector tile(32,2,1) private(w_val, abs_u_val, q0, qn1)
            do j = j_s,j_e
               do k = kms+1,kme
                   do i = i_s,i_e
                       ! Cache values to reduce memory accesses
                       w_val = W_m(i,k-1,j)
                       abs_u_val = ABS(w_val)
                       q0 = q(i,k,j)
                       qn1 = q(i,k-1,j)
                       
                       flux_z(i,k,j) = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                   enddo
               enddo
            enddo
            !$acc loop gang vector tile(64,1) private(w_val, q0)
            do j = j_s,j_e
                do i = i_s,i_e
                    ! flux_z(i,kms,j) = 0
                    q0 = q(i,kme,j)
                    w_val = W_m(i,kme,j)
                    flux_z(i,kme+1,j) = q0 * w_val * t_factor
                enddo
            enddo
            !$acc end parallel
        else if (vorder==3) then
            coef = (1./12)*t_factor
            ! Using tile(32,2) for vertical direction with more k levels
            !$acc parallel loop gang vector tile(32,2,1) async(2) private(u_val, abs_u_val, tmp)
            do j = j_s,j_e
            do k = kms+1,kme+1
            do i = i_s,i_e
                u = W_m(i,k-1,j)
                qn1 = q(i,k-1,j); 
                if (k==kms+1 .or. k==kme) then
                    q0  = q(i,k,j);
                    flux_z(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * t_factor
                elseif(k==kme+1) then
                    flux_z(i,k,j) = qn1 * u * t_factor
                else 
                    q0  = q(i,k,j);
                    q1  = q(i,k+1,j)
                    qn2 = q(i,k-2,j)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                    flux_z(i,k,j) = tmp*coef         
                end if
            enddo
            enddo
            enddo
        else if (vorder==5) then
            coef = (1./60)*t_factor
            !$acc parallel async(2)
            ! Using tile(32,2) for complex 5th order vertical computations
            !$acc loop gang vector tile(32,2,1) private(u_val, abs_u_val, q_cache, tmp)
            do j = j_s,j_e
            do k = kms+3,kme-2
            do i = i_s,i_e
                ! Cache all required q values to minimize memory accesses
                q_cache(1) = q(i,k-3,j)  ! qn3
                q_cache(2) = q(i,k-2,j)  ! qn2  
                q_cache(3) = q(i,k-1,j)  ! qn1
                q_cache(4) = q(i,k,j)    ! q0
                q_cache(5) = q(i,k+1,j)  ! q1
                ! Note: q2 = q(i,k+2,j) accessed directly when needed
                
                u_val = W_m(i,k-1,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 6th order flux calculation using cached values
                tmp = 37.0 * (q_cache(4) + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (q(i,k+2,j) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 5th order diffusive terms using cached values
                tmp = tmp - abs_u_val * (10.0 * (q_cache(4) - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (q(i,k+2,j) - q_cache(1)))
                flux_z(i,k,j) = tmp * coef
            enddo
            enddo
            enddo
            
            coef = (1./12)*t_factor
            !$acc loop gang vector tile(64,1) private(u_val, abs_u_val, q_cache, tmp, q0, qn1)
            do j = j_s,j_e
                do i = i_s,i_e
                    
                    ! flux_z(i,kms,j) = 0

                    ! Boundary treatment for cells with insufficient stencil - use simpler upwind
                    u_val = W_m(i,kms,j)
                    abs_u_val = ABS(u_val)
                    q0 = q(i,kms+1,j)
                    qn1 = q(i,kms,j)
                    flux_z(i,kms+1,j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * 0.5 * t_factor

                    ! Use 3rd order for kms+2 cell  
                    u_val = W_m(i,kms+1,j)
                    abs_u_val = ABS(u_val)
                    q_cache(1) = q(i,kms,j)    ! qn2
                    q_cache(2) = q(i,kms+1,j)  ! qn1
                    q_cache(3) = q(i,kms+2,j)  ! q0
                    q_cache(4) = q(i,kms+3,j)  ! q1
                    
                    tmp = 7.0 * (q_cache(3) + q_cache(2)) - (q_cache(4) + q_cache(1))
                    tmp = u_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (q_cache(3) - q_cache(2)) - (q_cache(4) - q_cache(1)))
                    flux_z(i,kms+2,j) = tmp * coef               
                    
                    ! Use 3rd order for kme-1 cell
                    u_val = W_m(i,kme-2,j)
                    abs_u_val = ABS(u_val)
                    q_cache(1) = q(i,kme-3,j)  ! qn2
                    q_cache(2) = q(i,kme-2,j)  ! qn1
                    q_cache(3) = q(i,kme-1,j)  ! q0
                    q_cache(4) = q(i,kme,j)    ! q1

                    tmp = 7.0 * (q_cache(3) + q_cache(2)) - (q_cache(4) + q_cache(1))
                    tmp = u_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (q_cache(3) - q_cache(2)) - (q_cache(4) - q_cache(1)))
                    flux_z(i,kme-1,j) = tmp * coef               

                    u_val = W_m(i,kme-1,j)
                    abs_u_val = ABS(u_val)
                    q0 = q(i,kme,j)
                    qn1 = q(i,kme-1,j)
                    flux_z(i,kme,j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * 0.5 * t_factor

                    flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j) * t_factor
                enddo
            enddo
            !$acc end parallel
        endif

        !$acc wait(1,2)
        !$acc end data
    end subroutine flux3

    !>------------------------------------------------------------
    !! FUSED FLUX COMPUTATION AND ADVECTION
    !! 
    !! Computes fluxes on-the-fly and immediately applies them
    !! Saves 6N reads + 3N writes compared to separate flux3 + summation
    !! 
    !! Only for cases WITHOUT flux correction (flux_corr == 0)
    !!------------------------------------------------------------
    subroutine flux_and_advect_fused(qfluxes,qold,U_m,V_m,W_m,denom,dz,t_factor,q_id)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(out)      :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)       :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)       :: dz, denom
        real, dimension(i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1), intent(in)       :: U_m, V_m, W_m
        real, intent(in) :: t_factor
        integer, intent(in) :: q_id
        
        integer :: i, j, k
        real :: t_factor_compact, coef
        real :: flux_x_i, flux_x_i1, flux_y_j, flux_y_j1, flux_z_k, flux_z_k1
        real :: q_cache(7), u_val, v_val, w_val, abs_u_val, tmp
        real :: q0, qn1, qn2, qn3, q1, q2
        
        !$acc data present(qold,U_m,V_m,W_m,qfluxes,denom,dz)
        t_factor_compact = 0.5 * t_factor

        if (horder==1 .and. vorder==1) then
            
            !$acc parallel loop gang vector tile(64,2,1) async(q_id) private(flux_x_i, flux_x_i1, flux_y_j, flux_y_j1, flux_z_k, flux_z_k1, u_val, v_val, w_val, abs_u_val, q0, qn1)
            do j = jts, jte
                do k = kms, kme
                    do i = its, ite
                        q0 = qold(i,k,j)
                        
                        ! Compute X-flux at i
                        qn1 = qold(i-1,k,j)
                        u_val = U_m(i,k,j)
                        abs_u_val = abs(u_val)
                        flux_x_i = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact
                        
                        ! Compute X-flux at i+1
                        q1 = qold(i+1,k,j)
                        u_val = U_m(i+1,k,j)
                        abs_u_val = abs(u_val)
                        flux_x_i1 = ((u_val + abs_u_val) * q0 + (u_val - abs_u_val) * q1) * t_factor_compact
                        
                        ! Compute Y-flux at j
                        qn1 = qold(i,k,j-1)
                        v_val = V_m(i,k,j)
                        abs_u_val = abs(v_val)
                        flux_y_j = ((v_val + abs_u_val) * qn1 + (v_val - abs_u_val) * q0) * t_factor_compact
                        
                        ! Compute Y-flux at j+1
                        q1 = qold(i,k,j+1)
                        v_val = V_m(i,k,j+1)
                        abs_u_val = abs(v_val)
                        flux_y_j1 = ((v_val + abs_u_val) * q0 + (v_val - abs_u_val) * q1) * t_factor_compact
                        
                        ! Compute Z-flux at k
                        if (k == kms) then
                            flux_z_k = 0.0
                        else
                            qn1 = qold(i,k-1,j)
                            w_val = W_m(i,k-1,j)
                            abs_u_val = abs(w_val)
                            flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                        endif
                        
                        ! Compute Z-flux at k+1
                        if (k == kme) then
                            w_val = W_m(i,kme,j)
                            flux_z_k1 = q0 * w_val * t_factor
                        else
                            q1 = qold(i,k+1,j)
                            w_val = W_m(i,k,j)
                            abs_u_val = abs(w_val)
                            flux_z_k1 = ((w_val + abs_u_val) * q0 + (w_val - abs_u_val) * q1) * t_factor_compact
                        endif
                        
                        ! Apply flux divergence
                        qfluxes(i,k,j) = q0 - ((flux_x_i1 - flux_x_i + &
                                                flux_y_j1 - flux_y_j + &
                                                (flux_z_k1 - flux_z_k) / dz(i,k,j)) * denom(i,k,j))
                    enddo
                enddo
            enddo
            
        else if (horder==3 .and. vorder==3) then
            coef = (1./12) * t_factor
            
            !$acc parallel loop gang vector tile(64,2,1) async(q_id) private(q_cache, tmp, flux_x_i, flux_x_i1, flux_y_j, flux_y_j1, flux_z_k, flux_z_k1, u_val, v_val, w_val, abs_u_val, q0, qn1, qn2, q1)
            do j = jts, jte
                do k = kms, kme
                    do i = its, ite
                        q0 = qold(i,k,j)
                        
                        ! Compute X-flux at i (3rd order)
                        q_cache(3) = qold(i-2,k,j)  ! qn2
                        q_cache(4) = qold(i-1,k,j)  ! qn1
                        q_cache(6) = qold(i+1,k,j)  ! q1
                        
                        u_val = U_m(i,k,j)
                        abs_u_val = abs(u_val)
                        tmp = 7.0 * (q0 + q_cache(4)) - (q_cache(6) + q_cache(3))
                        tmp = u_val * tmp - abs_u_val * (3.0 * (q0 - q_cache(4)) - (q_cache(6) - q_cache(3)))
                        flux_x_i = tmp * coef
                        
                        ! Compute X-flux at i+1 (3rd order)
                        q_cache(3) = q_cache(4)  ! old qn1 becomes new qn2
                        q_cache(4) = q_cache(6)          ! old q0 becomes new qn1
                        q_cache(6) = qold(i+2,k,j)  ! new q1
                        
                        u_val = U_m(i+1,k,j)
                        abs_u_val = abs(u_val)
                        tmp = 7.0 * (q_cache(4) + q0) - (q_cache(6) + q_cache(3))
                        tmp = u_val * tmp - abs_u_val * (3.0 * (q_cache(4) - q0) - (q_cache(6) - q_cache(3)))
                        flux_x_i1 = tmp * coef
                        
                        ! Compute Y-flux at j (3rd order)
                        q_cache(1) = qold(i,k,j-2)  ! qn2
                        q_cache(2) = qold(i,k,j-1)  ! qn1
                        q_cache(7) = qold(i,k,j+1)  ! q1
                        
                        v_val = V_m(i,k,j)
                        abs_u_val = abs(v_val)
                        tmp = 7.0 * (q0 + q_cache(2)) - (q_cache(7) + q_cache(1))
                        tmp = v_val * tmp - abs_u_val * (3.0 * (q0 - q_cache(2)) - (q_cache(7) - q_cache(1)))
                        flux_y_j = tmp * coef
                        
                        ! Compute Y-flux at j+1 (3rd order)
                        q_cache(1) = q_cache(2)  ! old qn1 becomes new qn2
                        q_cache(2) = q_cache(7)          ! old q0 becomes new qn1
                        q_cache(7) = qold(i,k,j+2)  ! new q1
                        
                        v_val = V_m(i,k,j+1)
                        abs_u_val = abs(v_val)
                        tmp = 7.0 * (q_cache(2) + q0) - (q_cache(7) + q_cache(1))
                        tmp = v_val * tmp - abs_u_val * (3.0 * (q_cache(2) - q0) - (q_cache(7) - q_cache(1)))
                        flux_y_j1 = tmp * coef
                        
                        ! Compute Z-fluxes (3rd order or upwind at boundaries)
                        if (k == kms) then
                            flux_z_k = 0.0
                            
                            !Use upwind for upward interface
                            qn1 = q0
                            q1 = qold(i,k+1,j)
                            w_val = W_m(i,k,j)
                            abs_u_val = abs(w_val)
                            flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact
                        
                        else if (k == kms+1) then
                            !Use upwind for downward interface
                            qn1 = qold(i,k-1,j)
                            w_val = W_m(i,k-1,j)
                            abs_u_val = abs(w_val)
                            flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                            
                            !Can already use 3rd order for upwind interface
                            qn2 = qn1
                            qn1 = q0
                            q1 = qold(i,k+2,j)
                            
                            w_val = W_m(i,k,j)
                            abs_u_val = abs(w_val)
                            tmp = 7.0 * (qold(i,k+1,j) + q0) - (q1 + qn2)
                            tmp = w_val * tmp - abs_u_val * (3.0 * (qold(i,k+1,j) - q0) - (q1 - qn2))
                            flux_z_k1 = tmp * (1./12) * t_factor

                        else if (k== kme-1) then
                            !Use 3rd order for downward interface
                            qn2 = qold(i,k-2,j)
                            qn1 = qold(i,k-1,j)
                            q1 = qold(i,k+1,j)
                            
                            w_val = W_m(i,k-1,j)
                            abs_u_val = abs(w_val)
                            tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                            tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                            flux_z_k = tmp * (1./12) * t_factor
                            
                            !Use upwind for upward interface
                            qn1 = q0
                            q1 = qold(i,k+1,j)
                            w_val = W_m(i,k,j)
                            abs_u_val = abs(w_val)
                            flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact

                        else if (k == kme) then
                            !Use upwind for downward interface
                            qn1 = qold(i,k-1,j)
                            w_val = W_m(i,k-1,j)
                            abs_u_val = abs(w_val)
                            flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                            
                            w_val = W_m(i,kme,j)
                            flux_z_k1 = q0 * w_val * t_factor
                        else
                            ! Both k and k+1 use 3rd order
                            qn2 = qold(i,k-2,j)
                            qn1 = qold(i,k-1,j)
                            q1 = qold(i,k+1,j)
                            q2 = qold(i,k+2,j)

                            w_val = W_m(i,k-1,j)
                            abs_u_val = abs(w_val)
                            tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                            tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                            flux_z_k = tmp * coef
                            
                            
                            w_val = W_m(i,k,j)
                            abs_u_val = abs(w_val)
                            tmp = 7.0 * (q1 + q0) - (q2 + qn1)
                            tmp = w_val * tmp - abs_u_val * (3.0 * (q1 - q0) - (q2 - qn1))
                            flux_z_k1 = tmp * coef
                        endif
                        
                        ! Apply flux divergence
                        qfluxes(i,k,j) = q0 - ((flux_x_i1 - flux_x_i + &
                                                flux_y_j1 - flux_y_j + &
                                                (flux_z_k1 - flux_z_k) / dz(i,k,j)) * denom(i,k,j))
                    enddo
                enddo
            enddo
            
        else if (horder==5 .and. vorder==3) then
            coef = (1./60) * t_factor
            
            !$acc parallel loop gang vector tile(64,2,1) async(q_id) private(q_cache, tmp, flux_x_i, flux_x_i1, flux_y_j, flux_y_j1, flux_z_k, flux_z_k1, u_val, v_val, w_val, abs_u_val, q0, qn1, qn2, qn3, q1)
            do j = jts, jte
            do k = kms, kme
                do i = its, ite
                q0 = qold(i,k,j)
                
                ! Compute X-flux at i (5th order)
                q_cache(1) = qold(i-3,k,j)  ! qn3
                q_cache(2) = qold(i-2,k,j)  ! qn2
                q_cache(3) = qold(i-1,k,j)  ! qn1
                q_cache(5) = qold(i+1,k,j)  ! q1
                
                u_val = U_m(i,k,j)
                abs_u_val = abs(u_val)
                tmp = 37.0 * (q0 + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i+2,k,j) + q_cache(1))
                tmp = u_val * tmp - abs_u_val * (10.0 * (q0 - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i+2,k,j) - q_cache(1)))
                flux_x_i = tmp * coef
                
                ! Compute X-flux at i+1 (5th order)
                q_cache(1) = q_cache(2)  ! shift cache
                q_cache(2) = q_cache(3)
                q_cache(3) = q0
                q_cache(5) = qold(i+2,k,j)
                
                u_val = U_m(i+1,k,j)
                abs_u_val = abs(u_val)
                tmp = 37.0 * (qold(i+1,k,j) + q0) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i+3,k,j) + q_cache(1))
                tmp = u_val * tmp - abs_u_val * (10.0 * (qold(i+1,k,j) - q0) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i+3,k,j) - q_cache(1)))
                flux_x_i1 = tmp * coef
                
                ! Compute Y-flux at j (5th order)
                q_cache(1) = qold(i,k,j-3)  ! qn3
                q_cache(2) = qold(i,k,j-2)  ! qn2
                q_cache(3) = qold(i,k,j-1)  ! qn1
                q_cache(5) = qold(i,k,j+1)  ! q1
                
                v_val = V_m(i,k,j)
                abs_u_val = abs(v_val)
                tmp = 37.0 * (q0 + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i,k,j+2) + q_cache(1))
                tmp = v_val * tmp - abs_u_val * (10.0 * (q0 - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i,k,j+2) - q_cache(1)))
                flux_y_j = tmp * coef
                
                ! Compute Y-flux at j+1 (5th order)
                q_cache(1) = q_cache(2)
                q_cache(2) = q_cache(3)
                q_cache(3) = q0
                q_cache(5) = qold(i,k,j+2)
                
                v_val = V_m(i,k,j+1)
                abs_u_val = abs(v_val)
                tmp = 37.0 * (qold(i,k,j+1) + q0) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i,k,j+3) + q_cache(1))
                tmp = v_val * tmp - abs_u_val * (10.0 * (qold(i,k,j+1) - q0) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i,k,j+3) - q_cache(1)))
                flux_y_j1 = tmp * coef
                
                ! Compute Z-fluxes (3rd order or upwind at boundaries)
                if (k == kms) then
                    flux_z_k = 0.0
                    
                    !Use upwind for upward interface
                    qn1 = q0
                    q1 = qold(i,k+1,j)
                    w_val = W_m(i,k,j)
                    abs_u_val = abs(w_val)
                    flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact
                
                else if (k == kms+1) then
                    !Use upwind for downward interface
                    qn1 = qold(i,k-1,j)
                    w_val = W_m(i,k-1,j)
                    abs_u_val = abs(w_val)
                    flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                    
                    !Can already use 3rd order for upwind interface
                    qn2 = qn1
                    qn1 = q0
                    q1 = qold(i,k+2,j)
                    
                    w_val = W_m(i,k,j)
                    abs_u_val = abs(w_val)
                    tmp = 7.0 * (qold(i,k+1,j) + q0) - (q1 + qn2)
                    tmp = w_val * tmp - abs_u_val * (3.0 * (qold(i,k+1,j) - q0) - (q1 - qn2))
                    flux_z_k1 = tmp * (1./12) * t_factor

                else if (k== kme-1) then
                    !Use 3rd order for downward interface
                    qn2 = qold(i,k-2,j)
                    qn1 = qold(i,k-1,j)
                    q1 = qold(i,k+1,j)
                    
                    w_val = W_m(i,k-1,j)
                    abs_u_val = abs(w_val)
                    tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                    tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                    flux_z_k = tmp * (1./12) * t_factor
                    
                    !Use upwind for upward interface
                    qn1 = q0
                    q1 = qold(i,k+1,j)
                    w_val = W_m(i,k,j)
                    abs_u_val = abs(w_val)
                    flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact

                else if (k == kme) then
                    !Use upwind for downward interface
                    qn1 = qold(i,k-1,j)
                    w_val = W_m(i,k-1,j)
                    abs_u_val = abs(w_val)
                    flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                    
                    w_val = W_m(i,kme,j)
                    flux_z_k1 = q0 * w_val * t_factor
                    
                else
                    qn2 = qold(i,k-2,j)
                    qn1 = qold(i,k-1,j)
                    q1 = qold(i,k+1,j)
                    
                    w_val = W_m(i,k-1,j)
                    abs_u_val = abs(w_val)
                    tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                    tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                    flux_z_k = tmp * (1./12) * t_factor
                    
                    qn2 = qn1
                    qn1 = q0
                    q1 = qold(i,k+2,j)
                    
                    w_val = W_m(i,k,j)
                    abs_u_val = abs(w_val)
                    tmp = 7.0 * (qold(i,k+1,j) + q0) - (q1 + qn2)
                    tmp = w_val * tmp - abs_u_val * (3.0 * (qold(i,k+1,j) - q0) - (q1 - qn2))
                    flux_z_k1 = tmp * (1./12) * t_factor
                endif
                
                ! Apply flux divergence
                qfluxes(i,k,j) = q0 - ((flux_x_i1 - flux_x_i + &
                            flux_y_j1 - flux_y_j + &
                            (flux_z_k1 - flux_z_k) / dz(i,k,j)) * denom(i,k,j))
                enddo
            enddo
            enddo
            
        else if (horder==5 .and. vorder==5) then
            coef = (1./60) * t_factor
            
            !$acc parallel loop gang vector tile(64,2,1) async(q_id) private(q_cache, tmp, flux_x_i, flux_x_i1, flux_y_j, flux_y_j1, flux_z_k, flux_z_k1, u_val, v_val, w_val, abs_u_val, q0, qn1, qn2, qn3, q1)
            do j = jts, jte
            do k = kms, kme
                do i = its, ite
                q0 = qold(i,k,j)
                
                ! Compute X-flux at i (5th order)
                q_cache(1) = qold(i-3,k,j)  ! qn3
                q_cache(2) = qold(i-2,k,j)  ! qn2
                q_cache(3) = qold(i-1,k,j)  ! qn1
                q_cache(5) = qold(i+1,k,j)  ! q1
                
                u_val = U_m(i,k,j)
                abs_u_val = abs(u_val)
                tmp = 37.0 * (q0 + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i+2,k,j) + q_cache(1))
                tmp = u_val * tmp - abs_u_val * (10.0 * (q0 - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i+2,k,j) - q_cache(1)))
                flux_x_i = tmp * coef
                
                ! Compute X-flux at i+1 (5th order)
                q_cache(1) = q_cache(2)
                q_cache(2) = q_cache(3)
                q_cache(3) = q0
                q_cache(5) = qold(i+2,k,j)
                
                u_val = U_m(i+1,k,j)
                abs_u_val = abs(u_val)
                tmp = 37.0 * (qold(i+1,k,j) + q0) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i+3,k,j) + q_cache(1))
                tmp = u_val * tmp - abs_u_val * (10.0 * (qold(i+1,k,j) - q0) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i+3,k,j) - q_cache(1)))
                flux_x_i1 = tmp * coef
                
                ! Compute Y-flux at j (5th order)
                q_cache(1) = qold(i,k,j-3)  ! qn3
                q_cache(2) = qold(i,k,j-2)  ! qn2
                q_cache(3) = qold(i,k,j-1)  ! qn1
                q_cache(5) = qold(i,k,j+1)  ! q1
                
                v_val = V_m(i,k,j)
                abs_u_val = abs(v_val)
                tmp = 37.0 * (q0 + q_cache(3)) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i,k,j+2) + q_cache(1))
                tmp = v_val * tmp - abs_u_val * (10.0 * (q0 - q_cache(3)) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i,k,j+2) - q_cache(1)))
                flux_y_j = tmp * coef
                
                ! Compute Y-flux at j+1 (5th order)
                q_cache(1) = q_cache(2)
                q_cache(2) = q_cache(3)
                q_cache(3) = q0
                q_cache(5) = qold(i,k,j+2)
                
                v_val = V_m(i,k,j+1)
                abs_u_val = abs(v_val)
                tmp = 37.0 * (qold(i,k,j+1) + q0) - 8.0 * (q_cache(5) + q_cache(2)) + (qold(i,k,j+3) + q_cache(1))
                tmp = v_val * tmp - abs_u_val * (10.0 * (qold(i,k,j+1) - q0) - 5.0 * (q_cache(5) - q_cache(2)) + (qold(i,k,j+3) - q_cache(1)))
                flux_y_j1 = tmp * coef
                
                ! Compute Z-fluxes (5th order or lower at boundaries)
                if (k <= kms+2 .or. k >= kme-2) then
                    ! Use upwind or 3rd order at boundaries
                    if (k == kms) then
                        flux_z_k = 0.0
                        
                        !Use upwind for upward interface
                        qn1 = q0
                        q1 = qold(i,k+1,j)
                        w_val = W_m(i,k,j)
                        abs_u_val = abs(w_val)
                        flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact
                    
                    else if (k == kms+1) then
                        !Use upwind for downward interface
                        qn1 = qold(i,k-1,j)
                        w_val = W_m(i,k-1,j)
                        abs_u_val = abs(w_val)
                        flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                        
                        !Can already use 3rd order for upwind interface
                        qn2 = qn1
                        qn1 = q0
                        q1 = qold(i,k+2,j)
                        
                        w_val = W_m(i,k,j)
                        abs_u_val = abs(w_val)
                        tmp = 7.0 * (qold(i,k+1,j) + q0) - (q1 + qn2)
                        tmp = w_val * tmp - abs_u_val * (3.0 * (qold(i,k+1,j) - q0) - (q1 - qn2))
                        flux_z_k1 = tmp * (1./12) * t_factor

                    else if (k == kms+2) then
                        !Use 3rd order for downward interface
                        qn2 = qold(i,k-2,j)
                        qn1 = qold(i,k-1,j)
                        q1 = qold(i,k+1,j)

                        w_val = W_m(i,k-1,j)
                        abs_u_val = abs(w_val)
                        tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                        tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                        flux_z_k = tmp * (1./12) * t_factor

                        !Use 5th order for upward interface
                        qn3 = qn2
                        qn2 = qn1
                        qn1 = q0
                        q1 = qold(i,k+2,j)
                        q2 = qold(i,k+3,j)
                        w_val = W_m(i,k,j)
                        abs_u_val = abs(w_val)
                        tmp = 37.0 * (qold(i,k+1,j) + q0) - 8.0 * (q1 + qn2) + (q2 + qn3)
                        tmp = w_val * tmp - abs_u_val * (10.0 * (qold(i,k+1,j) - q0) - 5.0 * (q1 - qn2) + (q2 - qn3))
                        flux_z_k1 = tmp * coef

                    else if (k== kme-2) then
                        !Use 5th order for downward interface
                        qn3 = qold(i,k-3,j)
                        qn2 = qold(i,k-2,j)
                        qn1 = qold(i,k-1,j)
                        q1 = qold(i,k+1,j)
                        w_val = W_m(i,k-1,j)
                        abs_u_val = abs(w_val)
                        tmp = 37.0 * (q0 + qn1) - 8.0 * (q1 + qn2) + (qold(i,k+2,j) + qn3)
                        tmp = w_val * tmp - abs_u_val * (10.0 * (q0 - qn1) - 5.0 * (q1 - qn2) + (qold(i,k+2,j) - qn3))
                        flux_z_k = tmp * coef

                        !Use 3rd order for upward interface
                        qn2 = qn1
                        qn1 = q0
                        q1 = qold(i,k+2,j)
                        w_val = W_m(i,k,j)
                        abs_u_val = abs(w_val)
                        tmp = 7.0 * (qold(i,k+1,j) + q0) - (q1 + qn2)
                        tmp = w_val * tmp - abs_u_val * (3.0 * (qold(i,k+1,j) - q0) - (q1 - qn2))
                        flux_z_k1 = tmp * (1./12) * t_factor
                        
                    else if (k== kme-1) then
                        !Use 3rd order for downward interface
                        qn2 = qold(i,k-2,j)
                        qn1 = qold(i,k-1,j)
                        q1 = qold(i,k+1,j)
                        
                        w_val = W_m(i,k-1,j)
                        abs_u_val = abs(w_val)
                        tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                        tmp = w_val * tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                        flux_z_k = tmp * (1./12) * t_factor
                        
                        !Use upwind for upward interface
                        qn1 = q0
                        q1 = qold(i,k+1,j)
                        w_val = W_m(i,k,j)
                        abs_u_val = abs(w_val)
                        flux_z_k1 = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q1) * t_factor_compact

                    else if (k == kme) then
                        !Use upwind for downward interface
                        qn1 = qold(i,k-1,j)
                        w_val = W_m(i,k-1,j)
                        abs_u_val = abs(w_val)
                        flux_z_k = ((w_val + abs_u_val) * qn1 + (w_val - abs_u_val) * q0) * t_factor_compact
                        
                        w_val = W_m(i,kme,j)
                        flux_z_k1 = q0 * w_val * t_factor
                    endif
                else
                    ! Use 5th order in interior
                    qn3 = qold(i,k-3,j)
                    qn2 = qold(i,k-2,j)
                    qn1 = qold(i,k-1,j)
                    q1 = qold(i,k+1,j)
                    
                    w_val = W_m(i,k-1,j)
                    abs_u_val = abs(w_val)
                    tmp = 37.0 * (q0 + qn1) - 8.0 * (q1 + qn2) + (qold(i,k+2,j) + qn3)
                    tmp = w_val * tmp - abs_u_val * (10.0 * (q0 - qn1) - 5.0 * (q1 - qn2) + (qold(i,k+2,j) - qn3))
                    flux_z_k = tmp * coef
                    
                    qn3 = qn2
                    qn2 = qn1
                    qn1 = q0
                    q1 = qold(i,k+2,j)
                    
                    w_val = W_m(i,k,j)
                    abs_u_val = abs(w_val)
                    tmp = 37.0 * (qold(i,k+1,j) + q0) - 8.0 * (q1 + qn2) + (qold(i,k+3,j) + qn3)
                    tmp = w_val * tmp - abs_u_val * (10.0 * (qold(i,k+1,j) - q0) - 5.0 * (q1 - qn2) + (qold(i,k+3,j) - qn3))
                    flux_z_k1 = tmp * coef
                endif
                
                ! Apply flux divergence
                qfluxes(i,k,j) = q0 - ((flux_x_i1 - flux_x_i + &
                            flux_y_j1 - flux_y_j + &
                            (flux_z_k1 - flux_z_k) / dz(i,k,j)) * denom(i,k,j))
                enddo
            enddo
            enddo
            
        else
            ! For other orders, fall back to computing fluxes then applying
            ! (Could add more fused kernels for other order combinations)

        endif
        
        !$acc end data
        !$acc wait(q_id)
    end subroutine flux_and_advect_fused


    subroutine adv_std_advect3d(qfluxes,qold,U_m,V_m,W_m,denom,dz, flux_time, flux_corr_time, sum_time, t_factor_in,flux_corr_in,q_id_in)
        ! !DIR$ INLINEALWAYS adv_std_advect3d
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz, denom
        real, dimension(i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1), intent(in)       :: U_m, V_m, W_m
        type(timer_t), optional, intent(inout) :: flux_time, flux_corr_time, sum_time
        real, optional,                              intent(in)      :: t_factor_in
        integer, optional,                           intent(in)      :: flux_corr_in, q_id_in

        ! interal parameters
        real    :: t_factor, flux_diff_x, flux_diff_y, denom_val, dz_val
        integer :: i, k, j, flux_corr, q_id

        
        q_id = 1
        if (present(q_id_in)) q_id = q_id_in
        !Initialize t_factor, which is used during RK time stepping to scale the time step
        t_factor = 1.0
        if (present(t_factor_in)) t_factor = t_factor_in
        
        flux_corr = 0
        if (present(flux_corr_in)) flux_corr = flux_corr_in

        ! Choose optimized path based on whether flux correction is needed
        if (flux_corr == 0) then
            ! if ((horder==vorder) .or. (horder==5 .and. vorder==3)) then
            !     ! FUSED PATH: No flux correction needed
            !     ! Compute fluxes and directly update qfluxes

            !     call flux_time%start()
            !     !$acc wait(q_id) !qold is needed for following function
            !     call flux_and_advect_fused(qfluxes,qold,U_m,V_m,W_m,denom,dz,t_factor,q_id)
            !     call flux_time%stop()
            ! else
                if(present(flux_time)) call flux_time%start()
                !$acc wait(q_id) !qold is needed for following function
                call flux3(qfluxes,U_m,V_m,W_m,flux_x,flux_z,flux_y,t_factor)
                if(present(flux_time)) call flux_time%stop()

                if(present(sum_time)) call sum_time%start()
                call sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
                if(present(sum_time)) call sum_time%stop()
            ! endif
            
        else
            ! STANDARD PATH: Flux correction enabled
            ! Must store fluxes for correction step
            if(present(flux_time)) call flux_time%start()
            call flux3(qfluxes,U_m, V_m, W_m, flux_x,flux_z,flux_y,t_factor)
            !$acc wait(q_id) !qold is not needed until after this point
            if(present(flux_time)) call flux_time%stop()
            
            if(present(flux_corr_time)) call flux_corr_time%start()
            ! Use async version which waits on q_id+100 then applies corrections
            call WRF_flux_corr(qold,U_m, V_m, W_m, flux_x,flux_z,flux_y,dz,denom,(q_id+100))
            if(present(flux_corr_time)) call flux_corr_time%stop()

            if(present(sum_time)) call sum_time%start()
            call sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
            if(present(sum_time)) call sum_time%stop()
        endif

    end subroutine adv_std_advect3d

    subroutine sum_kernel(flux_x, flux_y, flux_z, qold, qfluxes, denom, dz, q_id)
        implicit none
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in) :: flux_x, flux_y
        real, dimension(i_s:i_e+1,kms:kme+1,j_s:j_e+1), intent(in) :: flux_z
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: qold, denom, dz
        real, dimension(ims:ime,kms:kme,jms:jme), intent(out) :: qfluxes
        integer, intent(in) :: q_id

        integer :: i, j, k
        real :: flux_diff_z

        !$acc parallel loop gang vector async(q_id) tile(64,2,1) present(qfluxes, qold, flux_x, flux_y, flux_z, denom, dz) private(flux_diff_z)
        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    flux_diff_z = flux_z(i,k+1,j) - merge(0.0, flux_z(i,k,j), k==kms)
                    qfluxes(i,k,j) = qold(i,k,j) - ((flux_x(i+1,k,j) - flux_x(i,k,j) + &
                                        flux_y(i,k,j+1) - flux_y(i,k,j) + &
                                        flux_diff_z / dz(i,k,j)) * denom(i,k,j))
                enddo
            enddo
        enddo
        !$acc wait(q_id)

    end subroutine sum_kernel
    

    ! subroutine test_divergence(dz)
    !     implicit none
    !     real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)

    !     real, allocatable :: du(:,:), dv(:,:), dw(:,:)
    !     integer :: i,j,k

    !     allocate(du(i_s:i_e,j_s:j_e))
    !     allocate(dv(i_s:i_e,j_s:j_e))
    !     allocate(dw(i_s:i_e,j_s:j_e))

    !     do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e)

    !         du(i,j) = (U_m(i+1,k,j)-U_m(i,k,j))
    !         dv(i,j) = (V_m(i,k,j+1)-V_m(i,k,j))
    !         if (k==kms) then
    !             dw(i,j) = (W_m(i,k,j))/dz(i,k,j)
    !         else
    !             dw(i,j) = (W_m(i,k,j)-W_m(i,k-1,j))/dz(i,k,j)
    !         endif
    !         if (abs(du(i,j) + dv(i,j) + dw(i,j)) > 1e-3) then
    !             print*,  i,k,j , abs(du(i,j) + dv(i,j) + dw(i,j))
    !             print*, "Winds are not balanced on entry to advect"
    !             !error stop
    !         endif
    !     enddo

    ! end subroutine test_divergence

    subroutine adv_std_compute_wind(u,v,w,density,jaco,jaco_u,jaco_v,jaco_w,dz, dx, options, dt, U_m, V_m, W_m, denom)
        implicit none

        real, dimension(ims:ime+1,kms:kme,jms:jme), intent(in) :: u, jaco_u
        real, dimension(ims:ime,kms:kme,jms:jme+1), intent(in) :: v, jaco_v
        real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: w, jaco_w, density, dz, jaco
        type(options_t),    intent(in)  :: options
        real, allocatable, intent(out) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        real,intent(in)::dt, dx
        
        integer :: i, j, k
        real, dimension(ims:ime,kms:kme,jms:jme) :: rho
        logical :: advect_density

        ! if arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)
        if (allocated(flux_x)) deallocate(flux_x)
        if (allocated(flux_y)) deallocate(flux_y)
        if (allocated(flux_z)) deallocate(flux_z)

        ! allocate the arrays
        allocate(U_m     (i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1))
        allocate(V_m     (i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1))
        allocate(W_m     (i_s_w:i_e_w+1,kms:kme,j_s_w:j_e_w+1))
        allocate(denom   (ims:ime,  kms:kme,jms:jme  ))
        allocate(flux_x(i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(flux_y(i_s:i_e+1,  kms:kme,j_s:j_e+1))
        allocate(flux_z(i_s:i_e+1,  kms:kme+1,j_s:j_e+1))

        advect_density = options%adv%advect_density
        !$acc enter data create(U_m,V_m,W_m,denom, flux_x, flux_y, flux_z)


        !$acc data present(u,v,w,density,jaco,jaco_u,jaco_v,jaco_w,dz, dx, U_m, V_m, W_m, denom) create(rho)
        if (advect_density) then
            !$acc parallel loop gang vector collapse(3)
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        rho(i,k,j) = density(i,k,j)  
                    enddo
                enddo
            enddo
        else
            !$acc parallel loop gang vector collapse(3)
            do j = jms,jme
                do k = kms,kme
                    do i = ims,ime
                        rho(i,k,j) = 1
                    enddo
                enddo
            enddo
        endif

        !Compute the denomenator for all of the flux summation terms here once

        !$acc parallel loop gang vector collapse(3) async(1)
        do j = jms,jme
            do k = kms,kme
                do i = ims,ime
                    denom(i,k,j) = 1/(rho(i,k,j)*jaco(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector tile(32,2,1) async(2)
        do j = j_s_w,j_e_w+1
            do k = kms,kme
                do i = i_s_w,i_e_w+1
                    U_m(i,k,j) = u(i,k,j) * dt * (rho(i,k,j)+rho(i-1,k,j))*0.5 * &
                        jaco_u(i,k,j) / dx
                
                    V_m(i,k,j) = v(i,k,j) * dt * (rho(i,k,j)+rho(i,k,j-1))*0.5 * &
                        jaco_v(i,k,j) / dx
                enddo
            enddo
        enddo

        !$acc parallel loop gang vector tile(32,2,1) async(3)
        do j = j_s_w,j_e_w+1
            do k = kms,kme-1
                do i = i_s_w,i_e_w+1
                    W_m(i,k,j) = w(i,k,j) * dt * jaco_w(i,k,j) * &
                        ( rho(i,k,j)*dz(i,k+1,j) + &
                        rho(i,k+1,j)*dz(i,k,j) ) / &
                        (dz(i,k,j)+dz(i,k+1,j))
                enddo
            enddo
        enddo
        
        !$acc parallel loop gang vector collapse(2) async(4)
        do j = j_s_w,j_e_w+1
            do i = i_s_w,i_e_w+1
                W_m(i,kme,j) = w(i,kme,j) * dt * jaco_w(i,kme,j) * rho(i,kme,j)
            enddo
        enddo
        !$acc end data

        !$acc wait(1,2,3,4)
    end subroutine adv_std_compute_wind

    subroutine adv_std_clean_wind_arrays(U_m,V_m,W_m,denom)
        implicit none
        real, allocatable, dimension(:,:,:), intent(inout) :: U_m, V_m, W_m, denom

        !$acc exit data delete(U_m,V_m,W_m,denom, flux_x, flux_y, flux_z)
        
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)
        if (allocated(flux_x)) deallocate(flux_x)
        if (allocated(flux_y)) deallocate(flux_y)
        if (allocated(flux_z)) deallocate(flux_z)

    end subroutine adv_std_clean_wind_arrays


    !>------------------------------------------------------------
    !! Compute 3D wind Courant numbers for fine-mesh advection.
    !! Extends adv_std_compute_wind_horiz with vertical wind + Jacobian.
    !! All bounds are explicit (no module-level state).
    !!------------------------------------------------------------
    subroutine adv_std_compute_wind_3d_fm(u_cell, v_cell, w_cell, rho, &
        jaco, jaco_u, jaco_v, dz, dx, dt, &
        U_m, V_m, W_m, denom, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: u_cell, v_cell, w_cell, rho
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: jaco, jaco_u, jaco_v
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: dz
        real, intent(in) :: dx, dt
        real, allocatable, intent(out) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)

        integer :: i, j, k
        integer :: i_s_w_l, i_e_w_l, j_s_w_l, j_e_w_l
        integer :: i_s_l, i_e_l, j_s_l, j_e_l
        real :: dz_sum

        ! Extended bounds for flux correction:
        ! Flux computation range: its-1:ite+1  (extend interior by 1)
        ! Wind array range: its-2:ite+2        (extend flux range by 1)
        i_s_l   = its_l - 1
        i_e_l   = ite_l + 1
        j_s_l   = jts_l - 1
        j_e_l   = jte_l + 1
        i_s_w_l = its_l - 2
        i_e_w_l = ite_l + 2
        j_s_w_l = jts_l - 2
        j_e_w_l = jte_l + 2

        ! Deallocate if already allocated
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)
        if (allocated(flux_x_fm)) deallocate(flux_x_fm)
        if (allocated(flux_y_fm)) deallocate(flux_y_fm)
        if (allocated(flux_z_fm)) deallocate(flux_z_fm)

        ! Allocate wind arrays (extended bounds for flux correction)
        allocate(U_m    (i_s_w_l:i_e_w_l+1, ks:ke, j_s_w_l:j_e_w_l+1))
        allocate(V_m    (i_s_w_l:i_e_w_l+1, ks:ke, j_s_w_l:j_e_w_l+1))
        allocate(W_m    (i_s_w_l:i_e_w_l+1, ks:ke, j_s_w_l:j_e_w_l+1))
        allocate(denom  (ims_l:ime_l, ks:ke, jms_l:jme_l))

        ! Allocate module-level fine-mesh flux arrays
        allocate(flux_x_fm(i_s_l:i_e_l+1, ks:ke,   j_s_l:j_e_l+1))
        allocate(flux_y_fm(i_s_l:i_e_l+1, ks:ke,   j_s_l:j_e_l+1))
        allocate(flux_z_fm(i_s_l:i_e_l+1, ks:ke+1, j_s_l:j_e_l+1))

        ! Compute 1/(rho * jaco) denominator
        do j = jms_l, jme_l
            do k = ks, ke
                do i = ims_l, ime_l
                    denom(i,k,j) = 1.0 / (rho(i,k,j) * jaco(i,k,j))
                enddo
            enddo
        enddo

        ! U_m: face-staggered in x (mirrors adv_std_compute_wind line 1203)
        do j = j_s_w_l, j_e_w_l+1
            do k = ks, ke
                do i = i_s_w_l, i_e_w_l+1
                    U_m(i,k,j) = u_cell(i,k,j) * dt * &
                        0.5 * (rho(i-1,k,j) + rho(i,k,j)) * jaco_u(i,k,j) / dx
                enddo
            enddo
        enddo

        ! V_m: face-staggered in y (mirrors line 1206)
        do j = j_s_w_l, j_e_w_l+1
            do k = ks, ke
                do i = i_s_w_l, i_e_w_l+1
                    V_m(i,k,j) = v_cell(i,k,j) * dt * &
                        0.5 * (rho(i,k,j-1) + rho(i,k,j)) * jaco_v(i,k,j) / dx
                enddo
            enddo
        enddo

        ! W_m at interfaces: face-staggered in z
        ! jaco_w = 1.0 for AGL coordinates, so omitted
        do j = j_s_w_l, j_e_w_l+1
            do k = ks, ke-1
                do i = i_s_w_l, i_e_w_l+1
                    dz_sum = dz(i,k,j) + dz(i,k+1,j)
                    W_m(i,k,j) = w_cell(i,k,j) * dt * &
                        ( rho(i,k,j)*dz(i,k+1,j) + &
                        rho(i,k+1,j)*dz(i,k,j) ) / &
                        dz_sum

                enddo
            enddo
        enddo

        ! Top interface: W_m(ke) represents interface above cell ke
        do j = j_s_w_l, j_e_w_l+1
            do i = i_s_w_l, i_e_w_l+1
                W_m(i,ke,j) = w_cell(i,ke,j) * dt * rho(i,ke,j)
            enddo
        enddo

    end subroutine adv_std_compute_wind_3d_fm


    !>------------------------------------------------------------
    !! 3rd-order flux computation for fine-mesh advection.
    !! Stores fluxes in module-level flux_x_fm, flux_y_fm, flux_z_fm.
    !! Mirrors flux3 h_order=3 / v_order=3 with explicit bounds
    !! and fine-mesh vertical BCs (zero flux bottom, zero-gradient top).
    !!------------------------------------------------------------
    subroutine flux3_fm(q, U_m, V_m, W_m, t_factor, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in) :: q
        real, dimension(its_l-2:ite_l+3, ks:ke, jts_l-2:jte_l+3), intent(in) :: U_m, V_m, W_m
        real, intent(in) :: t_factor

        integer :: i, j, k
        integer :: i_s_l, i_e_l, j_s_l, j_e_l
        real :: coef, t_factor_up
        real :: u_val, v_val, w_val, abs_u_val, tmp
        real :: q0, qn1, qn2, q1

        ! Flux computation bounds (extended for flux correction)
        i_s_l = its_l - 1
        i_e_l = ite_l + 1
        j_s_l = jts_l - 1
        j_e_l = jte_l + 1

        coef = (1.0/12.0) * t_factor
        t_factor_up = 0.5 * t_factor

        ! ==========================================
        ! Horizontal fluxes (3rd order)
        ! ==========================================
        do j = j_s_l, j_e_l+1
            do k = ks, ke
                do i = i_s_l, i_e_l+1
                    ! X-direction flux
                    u_val = U_m(i,k,j)
                    abs_u_val = ABS(u_val)
                    tmp = 7.0 * (q(i,k,j) + q(i-1,k,j)) - (q(i+1,k,j) + q(i-2,k,j))
                    tmp = u_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (q(i,k,j) - q(i-1,k,j)) - (q(i+1,k,j) - q(i-2,k,j)))
                    flux_x_fm(i,k,j) = tmp * coef

                    ! Y-direction flux
                    v_val = V_m(i,k,j)
                    abs_u_val = ABS(v_val)
                    tmp = 7.0 * (q(i,k,j) + q(i,k,j-1)) - (q(i,k,j+1) + q(i,k,j-2))
                    tmp = v_val * tmp
                    tmp = tmp - abs_u_val * (3.0 * (q(i,k,j) - q(i,k,j-1)) - (q(i,k,j+1) - q(i,k,j-2)))
                    flux_y_fm(i,k,j) = tmp * coef
                enddo
            enddo
        enddo

        ! ==========================================
        ! Vertical fluxes (3rd order with fine-mesh BCs)
        ! ==========================================

        do j = j_s_l, j_e_l
            do i = i_s_l, i_e_l
                ! k=ks: zero flux at bottom (ground)
                flux_z_fm(i,ks,j) = 0.0

                ! k=ks+1: 1st-order upwind (insufficient stencil below)
                w_val = W_m(i,ks,j)
                abs_u_val = ABS(w_val)
                flux_z_fm(i,ks+1,j) = ((w_val + abs_u_val) * q(i,ks,j) + &
                    (w_val - abs_u_val) * q(i,ks+1,j)) * t_factor_up
            enddo
        enddo

        ! Interior vertical fluxes (3rd order): ks+2 to ke-1
        if (ke >= ks+3) then
            do j = j_s_l, j_e_l
                do k = ks+2, ke-1
                    do i = i_s_l, i_e_l
                        w_val = W_m(i,k-1,j)
                        abs_u_val = ABS(w_val)
                        q0  = q(i,k,j)
                        qn1 = q(i,k-1,j)
                        q1  = q(i,k+1,j)
                        qn2 = q(i,k-2,j)
                        tmp = 7.0 * (q0 + qn1) - (q1 + qn2)
                        tmp = w_val * tmp
                        tmp = tmp - abs_u_val * (3.0 * (q0 - qn1) - (q1 - qn2))
                        flux_z_fm(i,k,j) = tmp * coef
                    enddo
                enddo
            enddo
        endif

        do j = j_s_l, j_e_l
            do i = i_s_l, i_e_l
                ! k=ke: 1st-order upwind (insufficient stencil above)
                if (ke > ks+1) then
                    w_val = W_m(i,ke-1,j)
                    abs_u_val = ABS(w_val)
                    flux_z_fm(i,ke,j) = ((w_val + abs_u_val) * q(i,ke-1,j) + &
                        (w_val - abs_u_val) * q(i,ke,j)) * t_factor_up
                else
                    flux_z_fm(i,ke,j) = 0.0
                endif

                ! k=ke+1: zero-gradient outflow at top
                ! Use W_m(ke) (top boundary interface)
                flux_z_fm(i,ke+1,j) = q(i,ke,j) * W_m(i,ke,j) * t_factor
            enddo
        enddo

    end subroutine flux3_fm


    !>------------------------------------------------------------
    !! Apply flux divergence for fine-mesh advection.
    !! Mirrors sum_kernel with explicit bounds and fine-mesh BCs.
    !!------------------------------------------------------------
    subroutine sum_kernel_fm(qold, qfluxes, denom, dz, &
        ims_l, ime_l, ks, ke, jms_l, jme_l, its_l, ite_l, jts_l, jte_l)
        implicit none
        integer, intent(in) :: ims_l, ime_l, ks, ke, jms_l, jme_l
        integer, intent(in) :: its_l, ite_l, jts_l, jte_l
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(in)  :: qold, denom, dz
        real, dimension(ims_l:ime_l, ks:ke, jms_l:jme_l), intent(out) :: qfluxes

        integer :: i, j, k
        real :: flux_diff_z

        do j = jts_l, jte_l
            do k = ks, ke
                do i = its_l, ite_l
                    flux_diff_z = flux_z_fm(i,k+1,j) - merge(0.0, flux_z_fm(i,k,j), k==ks)
                    qfluxes(i,k,j) = qold(i,k,j) - ((flux_x_fm(i+1,k,j) - flux_x_fm(i,k,j) + &
                                        flux_y_fm(i,k,j+1) - flux_y_fm(i,k,j) + &
                                        flux_diff_z / dz(i,k,j)) * denom(i,k,j))
                enddo
            enddo
        enddo

    end subroutine sum_kernel_fm


    !>------------------------------------------------------------
    !! Deallocate fine-mesh wind and flux arrays.
    !! Mirrors adv_std_clean_wind_arrays for the fine mesh.
    !!------------------------------------------------------------
    subroutine adv_std_clean_wind_arrays_fm(U_m, V_m, W_m, denom)
        implicit none
        real, allocatable, dimension(:,:,:), intent(inout) :: U_m, V_m, W_m, denom

        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)
        if (allocated(flux_x_fm)) deallocate(flux_x_fm)
        if (allocated(flux_y_fm)) deallocate(flux_y_fm)
        if (allocated(flux_z_fm)) deallocate(flux_z_fm)

    end subroutine adv_std_clean_wind_arrays_fm


end module adv_std
