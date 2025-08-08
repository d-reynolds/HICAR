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

    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte, i_s, i_e, j_s, j_e, horder, vorder
    !type(timer_t) :: flux_time, flux_up_time, flux_corr_time, sum_time
    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v
    real, dimension(:,:,:), allocatable   :: flux_x, flux_x_up, flux_y, flux_y_up, flux_z, flux_z_up

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind, adv_std_clean_wind_arrays

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
        
        if (options%adv%flux_corr==kFLUXCOR_MONO) then
            i_s = its - 1
            i_e = ite + 1
            j_s = jts - 1
            j_e = jte + 1
        endif
        
    end subroutine

    subroutine adv_std_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for adv4 advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

        call options%advect_vars([kVARS%water_vapor, kVARS%potential_temperature])

        ! List the variables that are required for restarts with adv4 advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

    end subroutine

    subroutine flux3(q,U_m,V_m,W_m,flux_x,flux_z,flux_y,t_factor)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: U_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: V_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: W_m
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
            !$acc loop gang vector  collapse(3) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
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

            !$acc loop gang vector collapse(2) private(u_val, abs_u_val, q0, qn1)
            do j = j_s, j_e
                do k = kms, kme
                    qn1 = q(i_e,k,j)
                    q0 = q(i_e+1,k,j)

                    u_val = U_m(i_e+1,k,j)
                    abs_u_val = ABS(u_val)
                    flux_x(i_e+1, k, j) = ((u_val + abs_u_val) * qn1 + (u_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector collapse(2) private(v_val, abs_u_val, q0, qn3)
            do k = kms, kme
                do i = i_s, i_e
                    qn3 = q(i,k,j_e)
                    q0 = q(i,k,j_e+1)

                    v_val = V_m(i,k,j_e+1)
                    abs_u_val = ABS(v_val)  ! reuse variable
                    flux_y(i, k, j_e+1) = ((v_val + abs_u_val) * qn3 + (v_val - abs_u_val) * q0) * t_factor_compact
                enddo
            enddo

            !$acc loop gang vector collapse(2) private(u_val, v_val, abs_u_val, q0, qn1, qn3)
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

            !$acc parallel loop gang vector collapse(3) async(1) private(u_val, abs_u_val, q_cache, tmp)
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
                flux_x_up(i,k,j) = ((u_val + abs_u_val) * q_cache(4) + (u_val - abs_u_val) * q_cache(5)) * 0.25
                flux_x(i,k,j) = tmp * coef
                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 4th order flux calculation using cached values
                tmp = 7.0 * (q_cache(5) + q_cache(2)) - (q_cache(7) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 3rd order diffusive terms using cached values
                tmp = tmp - abs_u_val * (3.0 * (q_cache(5) - q_cache(2)) - (q_cache(7) - q_cache(1)))

                ! Optimized upwind flux calculation using cached values
                flux_y_up(i,k,j) = ((u_val + abs_u_val) * q_cache(2) + (u_val - abs_u_val) * q_cache(5)) * 0.25
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo

        else if (horder==5) then
            coef = (1./60)*t_factor
            !$acc parallel async(1)
            !$acc loop gang vector collapse(3) private(u_val, abs_u_val, q_cache, tmp)
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
            
            !$acc loop gang vector collapse(3) private(u_val, abs_u_val, q_cache, tmp)
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
            !$acc loop gang vector collapse(3) private(w_val, abs_u_val, q0, qn1)
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
            !$acc loop gang vector collapse(2) private(w_val, q0)
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
            !$acc parallel loop gang vector collapse(3) async(2) private(u_val, abs_u_val, tmp)
            do j = j_s,j_e
            do k = kms+1,kme+1
            do i = i_s,i_e
                u = W_m(i,k-1,j)
                qn1 = q(i,k-1,j); 
                q0  = q(i,k,j);
                if (k==kms+1 .or. k==kme) then
                    flux_z(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * t_factor
                elseif(k==kme+1) then
                    flux_z(i,k,j) = qn1 * u * t_factor
                else 
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
            !$acc loop gang vector collapse(3) private(u_val, abs_u_val, q_cache, tmp)
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
            !$acc loop gang vector collapse(2) private(u_val, abs_u_val, q_cache, tmp, q0, qn1)
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

    subroutine flux3_w_up(q,U_m,V_m,W_m,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),   intent(in)       :: q
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: U_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: V_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: W_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1),   intent(out)    :: flux_x, flux_x_up
        real, dimension(i_s:i_e+1,  kms:kme,j_s:j_e+1), intent(out)    :: flux_y, flux_y_up
        real, dimension(i_s:i_e+1,  kms:kme+1,j_s:j_e+1), intent(out)    :: flux_z, flux_z_up
        integer :: i, j, k, top, bot, bot2, nor, sou, sou2, eas, wes, wes2
        real :: tmp, coef, u, q0, q1, q2, qn1, qn2, qn3, abs_u, qin1, qin2, qi1, qk1, qj1, qkn1, qkn2, qjn1, qjn2
        ! GPU optimization variables
        real :: q_cache(7),  u_val, abs_u_val

        !$acc data present(q,U_m,V_m,W_m,flux_x,flux_x_up,flux_y,flux_y_up,flux_z,flux_z_up)

        if (horder==3) then
            coef = (1./12)
            !$acc parallel loop gang vector collapse(3) async(1) private(u_val, abs_u_val, q_cache, tmp)
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
                flux_x_up(i,k,j) = ((u_val + abs_u_val) * q_cache(4) + (u_val - abs_u_val) * q_cache(5)) * 0.25
                flux_x(i,k,j) = tmp * coef

                u_val = V_m(i,k,j)
                abs_u_val = ABS(u_val)
                
                ! Optimized 4th order flux calculation using cached values
                tmp = 7.0 * (q_cache(5) + q_cache(2)) - (q_cache(7) + q_cache(1))
                tmp = u_val * tmp
                
                ! Optimized 3rd order diffusive terms using cached values
                tmp = tmp - abs_u_val * (3.0 * (q_cache(5) - q_cache(2)) - (q_cache(7) - q_cache(1)))

                ! Optimized upwind flux calculation using cached values
                flux_y_up(i,k,j) = ((u_val + abs_u_val) * q_cache(2) + (u_val - abs_u_val) * q_cache(5)) * 0.25
                flux_y(i,k,j) = tmp * coef
            enddo
            enddo
            enddo

        else if (horder==5) then
            coef = (1./60)
            !$acc parallel async(1)
            !$acc loop gang vector collapse(3)
            do j = j_s,j_e
            do k = kms,kme
            do i = i_s,i_e+1
                u = U_m(i,k,j)
                q0  = q(i,k,j);   q1  = q(i+1,k,j); q2  = q(i+2,k,j)
                qn1 = q(i-1,k,j); qn2 = q(i-2,k,j); qn3 = q(i-3,k,j)

                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                !Calculation of Upwind fluxes
                flux_x_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                flux_x(i,k,j) = tmp*coef
            enddo
            enddo
            enddo
            !$acc loop gang vector collapse(3)
            do j = j_s,j_e+1
            do k = kms,kme
            do i = i_s,i_e
                u = V_m(i,k,j)
                q0  = q(i,k,j);   q1  = q(i,k,j+1); q2  = q(i,k,j+2)
                qn1 = q(i,k,j-1); qn2 = q(i,k,j-2); qn3 = q(i,k,j-3)
                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                !Calculation of Upwind fluxes
                flux_y_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                flux_y(i,k,j) = tmp*coef
            enddo
            enddo
            enddo
            !$acc end parallel
        endif

        if (vorder==1) then
            !$acc parallel async(2)
            !$acc loop gang vector collapse(3)
            do j = j_s,j_e
               do k = kms+1,kme
                   do i = i_s,i_e
                       flux_z(i,k,j) = ((W_m(i,k-1,j) + ABS(W_m(i,k-1,j))) * q(i,k-1,j) + &
                                    (W_m(i,k-1,j) - ABS(W_m(i,k-1,j))) * q(i,k,j))  * 0.5
                   enddo
               enddo
            enddo
            !$acc loop gang vector collapse(2)
            do j = j_s,j_e
                do i = i_s,i_e
                    ! flux_z(i,kms,j) = 0
                    ! flux_z_up(i,kms,j) = 0
                    flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j)
                    flux_z_up(i,kme+1,j) = flux_z(i,kme+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                enddo
            enddo
            !$acc end parallel

        else if (vorder==3) then
            coef = (1./12)
            !$acc parallel loop gang vector collapse(3) async(2) private(u_val, abs_u_val, tmp)
            do j = j_s,j_e
            do k = kms+1,kme+1
            do i = i_s,i_e
                u = W_m(i,k-1,j)
                qn1 = q(i,k-1,j); 
                q0  = q(i,k,j);
                if (k==kms+1 .or. k==kme) then
                    flux_z(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5
                    flux_z_up(i,k,j) = flux_z(i,k,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                elseif(k==kme+1) then
                    flux_z(i,k,j) = qn1 * u
                    flux_z_up(i,k,j) = flux_z(i,k,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                else 
                    q1  = q(i,k+1,j)
                    qn2 = q(i,k-2,j)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    !Calculation of Upwind fluxes
                    flux_z_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                    !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                    flux_z(i,k,j) = tmp*coef         
                end if
            enddo
            enddo
            enddo

        else if (vorder==5) then
            coef = (1./60)
            !$acc parallel
            !$acc loop gang vector collapse(3)
            do j = j_s,j_e
            do k = kms+3,kme-2
            do i = i_s,i_e
                u = W_m(i,k-1,j)
                q0  = q(i,k,j);   q1  = q(i,k+1,j);  q2 = q(i,k+2,j)
                qn1 = q(i,k-1,j); qn2 = q(i,k-2,j); qn3 = q(i,k-3,j)
                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                !Calculation of Upwind fluxes
                flux_z_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                flux_z(i,k,j) = tmp*coef
            enddo
            enddo
            enddo
            coef = (1./12)
            !$acc loop gang vector collapse(2)
            do j = j_s,j_e
                do i = i_s,i_e

                    ! flux_z(i,kms,j) = 0
                    ! flux_z_up(i,kms,j) = 0

                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + &
                                         (u - ABS(u)) * q(i,kms+1,j))  * 0.5
                    flux_z_up(i,kms+1,j) = flux_z(i,kms+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kms+1,j)
                    q0  = q(i,kms+2,j);   q1  = q(i,kms+3,j)
                    qn1 = q(i,kms+1,j);   qn2 = q(i,kms,j)
                    
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kms+2,j) = tmp*coef               
                    flux_z_up(i,kms+2,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kme-2,j)
                    q0  = q(i,kme-1,j);   q1  = q(i,kme,j)
                    qn1 = q(i,kme-2,j);   qn2 = q(i,kme-3,j)

                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kme-1,j) = tmp*coef               
                    flux_z_up(i,kme-1,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step


                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + &
                                       (u - ABS(u)) * q(i,kme,j))  * 0.5
                    flux_z_up(i,kme,j) = flux_z(i,kme,j) * 0.5 ! additional "0.5" since we only want half of the upwind step

                    flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j)
                    flux_z_up(i,kme+1,j) = flux_z(i,kme+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step

                enddo
            enddo
            !$acc end parallel
        endif
                                                          
        !$acc wait(1,2)
        !$acc end data
    end subroutine flux3_w_up

    subroutine adv_std_advect3d(qfluxes,qold,U_m,V_m,W_m,denom,dz, flux_time, flux_up_time, flux_corr_time, sum_time, t_factor_in,flux_corr_in)
        ! !DIR$ INLINEALWAYS adv_std_advect3d
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz, denom
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: U_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: V_m
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1), intent(in)       :: W_m
        type(timer_t), intent(inout) :: flux_time, flux_up_time, flux_corr_time, sum_time
        real, optional,                              intent(in)      :: t_factor_in
        integer, optional,                           intent(in)      :: flux_corr_in

        ! interal parameters
        real    :: t_factor, flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val
        integer :: i, k, j, flux_corr

        
        !Initialize t_factor, which is used during RK time stepping to scale the time step
        t_factor = 1.0
        if (present(t_factor_in)) t_factor = t_factor_in
        
        flux_corr = 0
        if (present(flux_corr_in)) flux_corr = flux_corr_in

        if (flux_corr > 0) then
            call flux_up_time%start()
            call flux3_w_up(qfluxes,U_m, V_m, W_m, flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up)
            call flux_up_time%stop()

            call flux_corr_time%start()
            call WRF_flux_corr(qold,U_m, V_m, W_m, flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up,dz,denom)
            call flux_corr_time%stop()
        else
            call flux_time%start()
            call flux3(qfluxes,U_m, V_m, W_m, flux_x,flux_z,flux_y,t_factor)
            call flux_time%stop()
        endif
        call sum_time%start()
        
        !$acc parallel loop gang vector collapse(3) present(qfluxes,qold,flux_x,flux_y,flux_z,denom,dz) private(flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val)
        do j = jts, jte
            do k = kms, kme
            do i = its, ite
                ! Cache frequently accessed values to reduce memory traffic
                flux_diff_x = flux_x(i+1,k,j) - flux_x(i,k,j)
                flux_diff_y = flux_y(i,k,j+1) - flux_y(i,k,j)

                !if k is the bottom level, we expect no flux in/out of bottom boundary. In fact,
                !advection code is optimized such that no flux is ever calculated...
                if (k==kms) then
                    flux_diff_z = flux_z(i,k+1,j)
                else
                    flux_diff_z = flux_z(i,k+1,j) - flux_z(i,k,j)
                endif
                denom_val = denom(i,k,j)
                dz_val = dz(i,k,j)
                
                ! Perform advection calculation with cached values
                qfluxes(i,k,j) = qold(i,k,j) - (flux_diff_x + flux_diff_y + flux_diff_z / dz_val) * denom_val
            enddo
            enddo
        enddo
        
        call sum_time%stop()

    end subroutine adv_std_advect3d
    

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
        if (allocated(flux_x_up)) deallocate(flux_x_up)
        if (allocated(flux_y)) deallocate(flux_y)
        if (allocated(flux_y_up)) deallocate(flux_y_up)
        if (allocated(flux_z)) deallocate(flux_z)
        if (allocated(flux_z_up)) deallocate(flux_z_up)

        ! allocate the arrays
        allocate(U_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(V_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(W_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(denom   (ims:ime,  kms:kme,jms:jme  ))
        allocate(flux_x(i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(flux_x_up(i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(flux_y(i_s:i_e+1,  kms:kme,j_s:j_e+1))
        allocate(flux_y_up(i_s:i_e+1,  kms:kme,j_s:j_e+1))
        allocate(flux_z(i_s:i_e+1,  kms:kme+1,j_s:j_e+1))
        allocate(flux_z_up(i_s:i_e+1,  kms:kme+1,j_s:j_e+1))

        advect_density = options%adv%advect_density
        !$acc enter data create(U_m,V_m,W_m,denom, flux_x, flux_y, flux_z, flux_x_up, flux_y_up, flux_z_up)


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

        !$acc parallel loop gang vector collapse(3)
        do j = jms,jme
            do k = kms,kme
                do i = ims,ime
                    denom(i,k,j) = 1/(rho(i,k,j)*jaco(i,k,j))
                enddo
            enddo
        enddo

        !$acc parallel
        !$acc loop gang vector collapse(3)
        do j = j_s,j_e+1
            do k = kms,kme
                do i = i_s,i_e+1
                    U_m(i,k,j) = u(i,k,j) * dt * (rho(i,k,j)+rho(i-1,k,j))*0.5 * &
                        jaco_u(i,k,j) / dx
                    V_m(i,k,j) = v(i,k,j) * dt * (rho(i,k,j)+rho(i,k,j-1))*0.5 * &
                        jaco_v(i,k,j) / dx
                enddo
            enddo
        enddo

        !$acc loop gang vector collapse(3)
        do j = j_s,j_e+1
            do k = kms,kme-1
                do i = i_s,i_e+1
                    W_m(i,k,j) = w(i,k,j) * dt * jaco_w(i,k,j) * &
                        ( rho(i,k,j)*dz(i,k+1,j) + &
                        rho(i,k+1,j)*dz(i,k,j) ) / &
                        (dz(i,k,j)+dz(i,k+1,j))
                enddo
            enddo
        enddo
        
        !$acc loop gang vector collapse(2)
        do j = j_s,j_e+1
            do i = i_s,i_e+1
                W_m(i,kme,j) = w(i,kme,j) * dt * jaco_w(i,kme,j) * rho(i,kme,j)
            enddo
        enddo
        !$acc end parallel
        !$acc end data

    end subroutine adv_std_compute_wind

    subroutine adv_std_clean_wind_arrays(U_m,V_m,W_m,denom)
        implicit none
        real, allocatable, dimension(:,:,:), intent(inout) :: U_m, V_m, W_m, denom

        !$acc exit data delete(U_m,V_m,W_m,denom, flux_x, flux_y, flux_z, flux_x_up, flux_y_up, flux_z_up)
        
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)
        if (allocated(flux_x)) deallocate(flux_x)
        if (allocated(flux_y)) deallocate(flux_y)
        if (allocated(flux_z)) deallocate(flux_z)
        if (allocated(flux_x_up)) deallocate(flux_x_up)
        if (allocated(flux_y_up)) deallocate(flux_y_up)
        if (allocated(flux_z_up)) deallocate(flux_z_up)

    end subroutine adv_std_clean_wind_arrays

end module adv_std
