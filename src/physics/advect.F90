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
    use adv_fluxcorr,      only: WRF_flux_corr, set_sign_arrays
    use timer_interface,   only: timer_t
    implicit none
    private

    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte, i_s, i_e, j_s, j_e, horder, vorder
    real    :: dx
    !type(timer_t) :: flux_time, flux_up_time, flux_corr_time, sum_time
    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind

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
                      
        if (horder==1) then
            t_factor_compact = 0.5  * t_factor
            do j = j_s, j_e
                do k = kms+1, kme
                    do i = i_s, i_e

                        qn3 = q(i,k,j-1)
                        qn2 = q(i,k-1,j)
                        qn1 = q(i-1,k,j)
                        q0 = q(i,k,j)

                        u = U_m(i,k,j)
                        flux_x(i, k, j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0) * t_factor_compact

                        v = V_m(i,k,j)
                        flux_y(i, k, j) = ((v + ABS(v)) * qn3 + (v - ABS(v)) * q0) * t_factor_compact

                        w = W_m(i,k-1,j)
                        flux_z(i,k,j) = ((w + ABS(w)) * qn2 +  (w - ABS(w)) * q0)  * t_factor_compact
                    enddo
                enddo
            enddo


            do j = j_s, j_e
                do k = kms, kme
                    qn1 = q(i_e,k,j)
                    q0 = q(i_e+1,k,j)

                    u = U_m(i_e+1,k,j)
                    flux_x(i_e+1, k, j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0) * t_factor_compact
                enddo
            enddo


            do k = kms, kme
                do i = i_s, i_e
                    qn3 = q(i,k,j_e)
                    q0 = q(i,k,j_e+1)

                    v = V_m(i,k,j_e+1)
                    flux_y(i, k, j_e+1) = ((v + ABS(v)) * qn3 + (v - ABS(v)) * q0) * t_factor_compact
                enddo
            enddo

            do j = j_s, j_e
                do i = i_s, i_e
                    qn3 = q(i,kms,j-1)
                    qn1 = q(i-1,kms,j)
                    q0 = q(i,kms,j)

                    u = U_m(i,kms,j)
                    flux_x(i, kms, j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0) * t_factor_compact

                    v = V_m(i,kms,j)
                    flux_y(i, kms, j) = ((v + ABS(v)) * qn3 + (v - ABS(v)) * q0) * t_factor_compact
                enddo
            enddo
        else if (horder==3) then
            coef = (1./12) * t_factor
            if (.False.) then
                do j = j_s, j_e+1
                    do k = kms, kme
                        do i = i_s, i_e+1
                            bot = max(k-1,kms)
                            bot2 = max(k-2,kms)
                            top = min(k+1,kme)
                            eas = min(i+1,ime)
                            wes = max(i-1,ims)
                            wes2 = max(i-2,ims)
                            nor = min(j+1,jme)
                            sou = max(j-1,jms)
                            sou2 = max(j-2,jms)

                            qjn2 = q(i,k,sou2); qjn1 = q(i,k,sou)
                            qkn2 = q(i,bot2,j); qkn1 = q(i,bot,j)
                            qin2 = q(wes2,k,j); qin1 = q(wes,k,j)
                            q0  = q(i,k,j)
                            qi1  = q(eas,k,j)
                            qk1  = q(i,top,j)
                            qj1  = q(i,k,nor)

                            u = U_m(i,k,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qin1) - (qi1+qin2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - abs_u * (3*(q0-qin1) - (qi1-qin2))
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_x(i,k,j) = tmp*coef
            
                            u = V_m(i,k,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qjn1) - (qj1+qjn2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - (abs_u) * (3*(q0-qjn1) - (qj1-qjn2))
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_y(i,k,j) = tmp*coef         

                            u = W_m(i,bot,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qkn1) - (qk1+qkn2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - abs_u * (3 * (q0 - qkn1) - (qk1 - qkn2))
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_z(i,k,j) = tmp*coef         
                        enddo
                    enddo
                enddo

                do j = j_s,j_e
                    do i = i_s,i_e
                        u = W_m(i,kms,j)
                        !Do simple upwind for the cells who's stencil does not allow higher-order
                        flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  * t_factor
                        
                        u = W_m(i,kme-1,j)
                        flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  * t_factor
                    enddo
                enddo
            else
                !DIR$ UNROLL 5
                do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e+1)
                    u = U_m(i,k,j)
                    q0  = q(i,k,j);   q1  = q(i+1,k,j)
                    qn1 = q(i-1,k,j); qn2 = q(i-2,k,j)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                    flux_x(i,k,j) = tmp*coef                
                enddo
                !DIR$ UNROLL 5
                do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e)
                    u = V_m(i,k,j)
                    q0  = q(i,k,j);   q1  = q(i,k,j+1)
                    qn1 = q(i,k,j-1); qn2 = q(i,k,j-2)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                    flux_y(i,k,j) = tmp*coef                
                enddo
            endif
        else if (horder==5) then
            coef = (1./60)*t_factor
            !DIR$ UNROLL 5
            do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e+1)
                u = U_m(i,k,j)
                q0  = q(i,k,j);   q1  = q(i+1,k,j); q2  = q(i+2,k,j)
                qn1 = q(i-1,k,j); qn2 = q(i-2,k,j); qn3 = q(i-3,k,j)

                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) - 5*(q1-qn2) + (q2-qn3))
                flux_x(i,k,j) = tmp*coef
            enddo
            !DIR$ UNROLL 5
            do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e)
                u = V_m(i,k,j)
                q0  = q(i,k,j);   q1  = q(i,k,j+1); q2  = q(i,k,j+2)
                qn1 = q(i,k,j-1); qn2 = q(i,k,j-2); qn3 = q(i,k,j-3)
                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                flux_y(i,k,j) = tmp*coef
            enddo
        endif
        
        if (vorder==1) then
            t_factor_compact = 0.5  * t_factor
            do j = j_s,j_e
               do k = kms+1,kme
                   do i = i_s,i_e
                       flux_z(i,k,j) = ((W_m(i,k-1,j) + ABS(W_m(i,k-1,j))) * q(i,k-1,j) + &
                                    (W_m(i,k-1,j) - ABS(W_m(i,k-1,j))) * q(i,k,j))  * t_factor_compact
                   enddo
               enddo
            enddo
        else if (vorder==3) then
            coef = (1./12)*t_factor
            do concurrent (j = j_s:j_e, k = kms+2:kme-1, i = i_s:i_e)
                u = W_m(i,k-1,j)
                q0  = q(i,k,j);   q1  = q(i,k+1,j)
                qn1 = q(i,k-1,j); qn2 = q(i,k-2,j)
                !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                tmp = 7*(q0+qn1) - (q1+qn2)
                tmp = u*tmp
                !Application of 3rd order diffusive terms
                tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                flux_z(i,k,j) = tmp*coef               
            enddo
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  * t_factor
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  * t_factor
                enddo
            enddo
        else if (vorder==5) then
            coef = (1./60)*t_factor
            do concurrent (j = j_s:j_e, k = kms+3:kme-2, i = i_s:i_e)
                u = W_m(i,k-1,j)
                q0  = q(i,k,j);   q1  = q(i,k+1,j);  q2 = q(i,k+2,j)
                qn1 = q(i,k-1,j); qn2 = q(i,k-2,j); qn3 = q(i,k-3,j)
                !Calculation of 6th order fluxes for later application of 5th order diffusive terms
                tmp = 37*(q0+qn1) - 8*(q1+qn2) + (q2+qn3)
                tmp = u*tmp
                !Application of 5th order diffusive terms
                tmp = tmp - abs(u) * (10*(q0-qn1) -  5*(q1-qn2) + (q2-qn3))
                flux_z(i,k,j) = tmp*coef
            enddo
            coef = (1./12)*t_factor
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms+1,j)
                    q0  = q(i,kms+2,j);   q1  = q(i,kms+3,j)
                    qn1 = q(i,kms+1,j);   qn2 = q(i,kms,j)
                    
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kms+2,j) = tmp*coef               
                    
                    u = W_m(i,kme-2,j)
                    q0  = q(i,kme-1,j);   q1  = q(i,kme,j)
                    qn1 = q(i,kme-2,j);   qn2 = q(i,kme-3,j)

                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                    flux_z(i,kme-1,j) = tmp*coef               

                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + &
                                         (u - ABS(u)) * q(i,kms+1,j))  * 0.5 * t_factor
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + &
                                       (u - ABS(u)) * q(i,kme,j))  * 0.5 * t_factor
                enddo
            enddo
        endif
                                                          
        !Handle top and bottom boundaries for z here
        do j = j_s,j_e
            do i = i_s,i_e
                flux_z(i,kms,j) = 0
                flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j) * t_factor
            enddo
        enddo
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
        
        if (horder==3) then
            coef = (1./12)
            if (.False.) then
                do j = j_s, j_e+1
                    do k = kms, kme
                        do i = i_s, i_e+1
                            bot =  max(k-1,kms)
                            bot2 = max(k-2,kms)
                            top =  min(k+1,kme)
                            eas =  min(i+1,ime)
                            wes =  max(i-1,ims)
                            wes2 = max(i-2,ims)
                            nor =  min(j+1,jme)
                            sou =  max(j-1,jms)
                            sou2 = max(j-2,jms)
                            
                            qjn2 = q(i,k,sou2); qjn1 = q(i,k,sou)
                            qkn2 = q(i,bot2,j); qkn1 = q(i,bot,j)
                            qin2 = q(wes2,k,j); qin1 = q(wes,k,j)
                            q0  = q(i,k,j)
                            qi1  = q(eas,k,j)
                            qk1  = q(i,top,j)
                            qj1  = q(i,k,nor)

                            u = U_m(i,k,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qin1) - (qi1+qin2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - abs_u * (3*(q0-qin1) - (qi1-qin2))
                            !Calculation of Upwind fluxes
                            flux_x_up(i,k,j) = ((u + abs_u) * qin1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_x(i,k,j) = tmp*coef
            
                            u = V_m(i,k,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qjn1) - (qj1+qjn2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - (abs_u) * (3*(q0-qjn1) - (qj1-qjn2))
                            !Calculation of Upwind fluxes
                            flux_y_up(i,k,j) = ((u + abs_u) * qjn1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_y(i,k,j) = tmp*coef         

                            u = W_m(i,bot,j)
                            abs_u = ABS(u)
                            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                            tmp = 7*(q0+qkn1) - (qk1+qkn2)
                            tmp = u*tmp
                            !Application of 3rd order diffusive terms
                            tmp = tmp - abs_u * (3 * (q0 - qkn1) - (qk1 - qkn2))
                            !Calculation of Upwind fluxes
                            flux_z_up(i,k,j) = ((u + abs_u) * qkn1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                            !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                            flux_z(i,k,j) = tmp*coef         
                        enddo
                    enddo
                enddo

                do j = j_s,j_e
                    do i = i_s,i_e
                        !Do simple upwind for the cells who's stencil does not allow higher-order
                        flux_z(i,kms+1,j) = flux_z_up(i,kms+1,j) * 2.0 ! ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  
                        
                        flux_z(i,kme,j) = flux_z_up(i,kme,j) * 2.0 ! ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  
                    enddo
                enddo
            else
            ! !DIR$ UNROLL 5
            !     do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e)
            !         u = U_m(i,k,j)
            !         q0  = q(i,k,j);   q1  = q(i+1,k,j)
            !         qn1 = q(i-1,k,j); qn2 = q(i-2,k,j)
            !         !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            !         tmp = 7*(q0+qn1) - (q1+qn2)
            !         tmp = u*tmp
            !         !Application of 3rd order diffusive terms
            !         tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
            !         !Calculation of Upwind fluxes
            !         flux_x_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
            !         !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
            !         flux_x(i,k,j) = tmp*coef
    
            !         u = V_m(i,k,j)
            !         q1  = q(i,k,j+1)
            !         qn1 = q(i,k,j-1); qn2 = q(i,k,j-2)
            !         !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            !         tmp = 7*(q0+qn1) - (q1+qn2)
            !         tmp = u*tmp
            !         !Application of 3rd order diffusive terms
            !         tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
            !         !Calculation of Upwind fluxes
            !         flux_y_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
            !         !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
            !         flux_y(i,k,j) = tmp*coef         
    
            !     enddo
    
            !     do concurrent (j = j_s:j_e, k = kms:kme)
            !         u = U_m(i_e+1,k,j)
            !         q0  = q(i_e+1,k,j);   q1  = q(i_e+2,k,j)
            !         qn1 = q(i_e,k,j); qn2 = q(i_e-1,k,j)
            !         !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            !         tmp = 7*(q0+qn1) - (q1+qn2)
            !         tmp = u*tmp
            !         !Application of 3rd order diffusive terms
            !         tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
            !         !Calculation of Upwind fluxes
            !         flux_x_up(i_e+1,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
            !         !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
            !         flux_x(i_e+1,k,j) = tmp*coef
            !     enddo
    
            !     do concurrent (k = kms:kme, i = i_s:i_e)
            !         u = V_m(i,k,j_e+1)
            !         q0  = q(i,k,j_e+1);   q1  = q(i,k,j_e+2)
            !         qn1 = q(i,k,j_e); qn2 = q(i,k,j_e-1)
            !         !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            !         tmp = 7*(q0+qn1) - (q1+qn2)
            !         tmp = u*tmp
            !         !Application of 3rd order diffusive terms
            !         tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
            !         !Calculation of Upwind fluxes
            !         flux_y_up(i,k,j_e+1) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
            !         !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
            !         flux_y(i,k,j_e+1) = tmp*coef         
            !     enddo    

                !DIR$ UNROLL 5
                do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e+1)
                    u = U_m(i,k,j)
                    q0  = q(i,k,j);   q1  = q(i+1,k,j)
                    qn1 = q(i-1,k,j); qn2 = q(i-2,k,j)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                    !Calculation of Upwind fluxes
                    flux_x_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                    !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                    flux_x(i,k,j) = tmp*coef                
                enddo
                !DIR$ UNROLL 5
                do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e)
                    u = V_m(i,k,j)
                    q0  = q(i,k,j);   q1  = q(i,k,j+1)
                    qn1 = q(i,k,j-1); qn2 = q(i,k,j-2)
                    !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                    tmp = 7*(q0+qn1) - (q1+qn2)
                    tmp = u*tmp
                    !Application of 3rd order diffusive terms
                    tmp = tmp - (abs(u)) * (3*(q0-qn1) - (q1-qn2))
                    !Calculation of Upwind fluxes
                    flux_y_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                    !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                    flux_y(i,k,j) = tmp*coef
                enddo
            endif
        else if (horder==5) then
            coef = (1./60)
            !DIR$ UNROLL 5
            do concurrent (j = j_s:j_e, k = kms:kme, i = i_s:i_e+1)
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
            !DIR$ UNROLL 5
            do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e)
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
        endif

        if (vorder==1) then
            do j = j_s,j_e
               do k = kms+1,kme
                   do i = i_s,i_e
                       flux_z(i,k,j) = ((W_m(i,k-1,j) + ABS(W_m(i,k-1,j))) * q(i,k-1,j) + &
                                    (W_m(i,k-1,j) - ABS(W_m(i,k-1,j))) * q(i,k,j))  * 0.5
                   enddo
               enddo
            enddo
        else if (vorder==3) then
            coef = (1./12)
            do concurrent (j = j_s:j_e, k = kms+2:kme-1, i = i_s:i_e)
                u = W_m(i,k-1,j)
                q0  = q(i,k,j);   q1  = q(i,k+1,j)
                qn1 = q(i,k-1,j); qn2 = q(i,k-2,j)
                !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
                tmp = 7*(q0+qn1) - (q1+qn2)
                tmp = u*tmp
                !Application of 3rd order diffusive terms
                tmp = tmp - abs(u) * (3 * (q0 - qn1) - (q1 - qn2))
                !Calculation of Upwind fluxes
                flux_z_up(i,k,j) = ((u + ABS(u)) * qn1 + (u - ABS(u)) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
                !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
                flux_z(i,k,j) = tmp*coef         
            enddo
            do j = j_s,j_e
                do i = i_s,i_e
                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + (u - ABS(u)) * q(i,kms+1,j))  * 0.5  
                    flux_z_up(i,kms+1,j) = flux_z(i,kms+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                    
                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + (u - ABS(u)) * q(i,kme,j))  * 0.5  
                    flux_z_up(i,kme,j) = flux_z(i,kme,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                enddo
            enddo
        else if (vorder==5) then
            coef = (1./60)
            do concurrent (j = j_s:j_e, k=kms+3:kme-2, i = i_s:i_e)
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
            coef = (1./12)
            do j = j_s,j_e
                do i = i_s,i_e
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

                    u = W_m(i,kms,j)
                    !Do simple upwind for the cells who's stencil does not allow higher-order
                    flux_z(i,kms+1,j) = ((u + ABS(u)) * q(i,kms,j) + &
                                         (u - ABS(u)) * q(i,kms+1,j))  * 0.5
                    flux_z_up(i,kms+1,j) = flux_z(i,kms+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step

                    u = W_m(i,kme-1,j)
                    flux_z(i,kme,j) = ((u + ABS(u)) * q(i,kme-1,j) + &
                                       (u - ABS(u)) * q(i,kme,j))  * 0.5
                    flux_z_up(i,kme,j) = flux_z(i,kme,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
                enddo
            enddo
        endif
                                                          
        !Handle top and bottom boundaries for z here
        do j = j_s,j_e
            do i = i_s,i_e
                flux_z(i,kms,j) = 0
                flux_z_up(i,kms,j) = 0
                flux_z(i,kme+1,j) = q(i,kme,j) * W_m(i,kme,j)
                flux_z_up(i,kme+1,j) = flux_z(i,kme+1,j) * 0.5 ! additional "0.5" since we only want half of the upwind step
            enddo
        enddo
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
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e+1)   :: flux_x, flux_x_up
        real, dimension(i_s:i_e+1,  kms:kme,j_s:j_e+1) :: flux_y, flux_y_up
        real, dimension(i_s:i_e+1,  kms:kme+1,j_s:j_e+1) :: flux_z, flux_z_up
        real    :: t_factor
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
            call WRF_flux_corr(qold,U_m,V_m,W_m,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up,dz,denom)
            call flux_corr_time%stop()
        else
            call flux_time%start()
            call flux3(qfluxes,U_m, V_m, W_m, flux_x,flux_z,flux_y,t_factor)
            call flux_time%stop()
        endif

        ! !$omp parallel default(shared) &
        ! !$omp private(i,j,k) &
        ! !$omp firstprivate(its, ite, jts, jte, kms,kme)
        ! !$omp do
        call sum_time%start()

        do j = jts, jte
            do k = kms, kme
                do i = its, ite
                    ! perform advection, from difference terms
                    qfluxes(i,k,j)  = qold(i,k,j)  - &
                        ((flux_x(i+1,k,j) - flux_x(i,k,j))  + &
                        (flux_y(i,k,j+1) - flux_y(i,k,j))  + &
                        (flux_z(i,k+1,j) - flux_z(i,k,j))  / &
                        dz(i,k,j))*denom(i,k,j)
                enddo
            enddo
        enddo
        call sum_time%stop()

        ! do j = jts, jte
        !     do k = kms, kme
        !         do i = its, ite

        !             bot =  max(k-1,kms)
        !             bot2 = max(k-2,kms)
        !             top =  min(k+1,kme)
        !             eas =  min(i+1,ime)
        !             wes =  max(i-1,ims)
        !             wes2 = max(i-2,ims)
        !             nor =  min(j+1,jme)
        !             sou =  max(j-1,jms)
        !             sou2 = max(j-2,jms)

        !             qjn2 = q(i,k,sou2); qjn1 = q(i,k,sou)
        !             qkn2 = q(i,bot2,j); qkn1 = q(i,bot,j); 
        !             qin2 = q(wes2,k,j); qin1 = q(wes,k,j)
        !             q0  = q(i,k,j)
        !             qi1  = q(eas,k,j)
        !             qk1  = q(i,top,j)
        !             qj1  = q(i,k,nor)

        !             u = U_m(i,k,j)
        !             abs_u = ABS(u)
        !             !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
        !             tmp = 7*(q0+qin1) - (qi1+qin2)
        !             tmp = u*tmp
        !             !Application of 3rd order diffusive terms
        !             tmp = tmp - abs_u * (3*(q0-qin1) - (qi1-qin2))
        !             !Calculation of Upwind fluxes
        !             flux_x_up(i,k,j) = ((u + abs_u) * qin1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
        !             !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
        !             flux_x(i,k,j) = tmp*coef
    
        !             u = V_m(i,k,j)
        !             abs_u = ABS(u)
        !             !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
        !             tmp = 7*(q0+qjn1) - (qj1+qjn2)
        !             tmp = u*tmp
        !             !Application of 3rd order diffusive terms
        !             tmp = tmp - (abs_u) * (3*(q0-qjn1) - (qj1-qjn2))
        !             !Calculation of Upwind fluxes
        !             flux_y_up(i,k,j) = ((u + abs_u) * qjn1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
        !             !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
        !             flux_y(i,k,j) = tmp*coef         

        !             u = W_m(i,bot,j)
        !             abs_u = ABS(u)
        !             !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
        !             tmp = 7*(q0+qkn1) - (qk1+qkn2)
        !             tmp = u*tmp
        !             !Application of 3rd order diffusive terms
        !             tmp = tmp - abs_u * (3 * (q0 - qkn1) - (qk1 - qkn2))
        !             !Calculation of Upwind fluxes
        !             flux_z_up(i,k,j) = ((u + abs_u) * qkn1 + (u - abs_u) * q0)  * 0.5 * 0.5 ! additional "0.5" since we only want half of the upwind step
        !             !Application of Upwind fluxes to higher order fluxes -- needed for flux correction step
        !             flux_z(i,k,j) = tmp*coef         

        !             ! perform advection, from difference terms
        !             qfluxes(i,k,j)  = qold(i,k,j)  - &
        !                 ((flux_x(i+1,k,j) - flux_x(i,k,j))  + &
        !                 (flux_y(i,k,j+1) - flux_y(i,k,j))  + &
        !                 (flux_z(i,k+1,j) - flux_z(i,k,j))  / &
        !                 dz(i,k,j))*denom(i,k,j)

        !         enddo
        !     enddo
        ! enddo


        ! !$omp end do nowait
        ! !$omp end parallel


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

    subroutine adv_std_compute_wind(domain, options, dt, U_m, V_m, W_m, denom)
        implicit none

        type(options_t),    intent(in)  :: options
        type(domain_t),  intent(inout) :: domain
        real, allocatable, intent(out) :: U_m(:,:,:), V_m(:,:,:), W_m(:,:,:), denom(:,:,:)
        real,intent(in)::dt
        
        integer :: i, j, k
        real, dimension(domain%ims:domain%ime,domain%kms:domain%kme,domain%jms:domain%jme) :: rho
        
        ! if arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        if (allocated(denom)) deallocate(denom)

        ! allocate the arrays
        allocate(U_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(V_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(W_m     (i_s:i_e+1,kms:kme,j_s:j_e+1))
        allocate(denom   (ims:ime,  kms:kme,jms:jme  ))

        dx = domain%dx
        
        if (options%adv%advect_density) then
            do i = ims,ime
                do j = jms,jme
                    do k = kms,kme
                        rho(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d(i,k,j)  
                    enddo
                enddo
            enddo
        else
            do i = ims,ime
                do j = jms,jme
                    do k = kms,kme
                        rho(i,k,j) = 1
                    enddo
                enddo
            enddo
        endif

        !Compute the denomenator for all of the flux summation terms here once
        denom = 1/(rho*domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)

        do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e+1)
            U_m(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d(i,k,j) * dt * (rho(i,k,j)+rho(i-1,k,j))*0.5 * &
                    domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i,k,j) / domain%dx
        enddo

        do concurrent (j = j_s:j_e+1, k = kms:kme, i = i_s:i_e+1)
            V_m(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d(i,k,j) * dt * (rho(i,k,j)+rho(i,k,j-1))*0.5 * &
                    domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i,k,j) / domain%dx
        enddo

        do concurrent (j = j_s:j_e+1, k = kms:kme-1, i = i_s:i_e+1)
            W_m(i,k,j) = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,k,j) * dt * domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i,k,j) * &
                    ( rho(i,k,j)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k+1,j) + &
                        rho(i,k+1,j)*domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k,j) ) / &
                        (domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k,j)+domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d(i,k+1,j))
        enddo

        do j = j_s,j_e+1
            do i = i_s,i_e+1
                W_m(i,kme,j) = domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d(i,kme,j) * dt * domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d(i,kme,j) * rho(i,kme,j)
            enddo
        enddo

        if (options%adv%flux_corr==kFLUXCOR_MONO) then
            call set_sign_arrays(U_m,V_m,W_m)
        endif
    end subroutine adv_std_compute_wind

end module adv_std
