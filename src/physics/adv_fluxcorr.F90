!> ----------------------------------------------------------------------------
!!  A collection of flux correction schemes for advection
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_fluxcorr
    use data_structures
    use domain_interface,  only: domain_t
    use mpi_f08, only: MPI_Wtime
    implicit none
    private
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte

    public :: WRF_flux_corr, init_fluxcorr, set_sign_arrays
    integer, allocatable, dimension(:,:,:)   :: usign, vsign, wsign

contains

    subroutine init_fluxcorr(domain)
        implicit none
        type(domain_t), intent(in) :: domain
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

        if (allocated(usign)) deallocate(usign)
        if (allocated(vsign)) deallocate(vsign)
        if (allocated(wsign)) deallocate(wsign)
        
        allocate(usign(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(vsign(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(wsign(its-1:ite+1,kms:kme,jts-1:jte+1))
    end subroutine init_fluxcorr

    subroutine set_sign_arrays(u,v,w)
        implicit none
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+1),  intent(in) :: u
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+2),  intent(in) :: v

        integer :: i, j, k

        do i = its-1,ite+1
            do j = jts-1,jte+1
                do k = kms,kme
                    if (u(i,k,j) > 0) then
                        usign(i,k,j) = -1
                    elseif (u(i+1,k,j) < 0) then
                        usign(i,k,j) = 1
                    else
                        usign(i,k,j) = 0
                    end if
                    if (v(i,k,j) > 0) then
                        vsign(i,k,j) = -1
                    elseif (v(i,k,j+1) < 0) then
                        vsign(i,k,j) = 1
                    else
                        vsign(i,k,j) = 0
                    end if
                    if (k > kms) then
                        if (w(i,k,j) < 0) then
                            if (k < kme) then
                                wsign(i,k,j) = 1
                            else
                                wsign(i,k,j) = 0
                            end if
		        elseif (w(i,k-1,j) > 0) then
                            wsign(i,k,j) = -1
                        else
                            wsign(i,k,j) = 0
                        end if
                    else
                        if (w(i,k,j) < 0) then
                            wsign(i,k,j) = 1
                        else
                            wsign(i,k,j) = 0
                        end if
                    end if
                end do
            end do
        end do
    end subroutine set_sign_arrays

    subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,flux_x_up,flux_z_up,flux_y_up,dz,denom)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, dz, denom
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: u
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: v
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_x, flux_x_up
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_y, flux_y_up
        real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout)    :: flux_z, flux_z_up
        
        
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1)   :: scale_in, scale_out, q3, q4
        real :: dz_t_i, jaco_rho_t_i, fx, fx1, fy, fy1, fz, fz1, qmax, qmin, q_i, q_j, q_k, q0, temp, flux_in, flux_out
        integer :: i, j ,k
        real :: scale
        real :: flux_x_up_0, flux_y_up_0, flux_z_up_0, flux_x_up_1, flux_y_up_1, flux_z_up_1
        integer :: bot_in, bot_out, sou_in, sou_out, wes_in, wes_out, j_block, k_block, usign0, vsign0, wsign0

        !Initialize some internal variables

        ! Get upwind fluxes

        call upwind_flux3(q,u,v,w,flux_x_up,flux_z_up,flux_y_up,dz,denom)

        ! Next compute max and min possible fluxes
        ! $omp parallel default(none) &
        ! $omp shared(q,flux_x_up,flux_y_up,flux_z_up,flux_x,flux_y,flux_z) &
        ! $omp shared(scale_in,scale_out,flux_x_cpy,flux_y_cpy,flux_z_cpy) &
        ! $omp shared(dz,denom,usign,vsign,wsign) &
        ! $omp private(i,j,k, j_block, k_block, q0, q_i, q_j, q_k, qmin, qmax, dz_t_i) &
        ! $omp private(fx, fx1, fy, fy1, fz, fz1, temp, flux_in, flux_out) &
        ! $omp private(flux_x_up_0, flux_y_up_0, flux_z_up_0, flux_x_up_1, flux_y_up_1, flux_z_up_1, bot, wes, sou) &
        ! $omp firstprivate(its, ite, jts, jte, kms,kme)
        ! $omp do schedule(static) collapse(2)
        ! $omp simd collapse(3)
        do j = jts-1, jte+1 
            do k = kms, kme
                do i = its-1, ite+1

                    usign0 = usign(i,k,j)
                    vsign0 = vsign(i,k,j)
                    wsign0 = wsign(i,k,j)

                    q0 = q(i,k,j)
                    q_i = q(i+usign0,k,j)
                    q_j = q(i,k,j+vsign0)
                    q_k = q(i,k+wsign0,j)

                    qmax = max(q0, q_i, q_j, q_k)
                    qmin = min(q0, q_i, q_j, q_k)

                    !This is the original code, which is may be slower than the above
                    !included code, but is more readable
                    ! if (u(i,k,j) > 0) then
                    !     qmax(i,k,j) = max(q(i-1,k,j),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i-1,k,j),qmin(i,k,j))
                    ! else if (u(i+1,k,j) < 0) then
                    !     qmax(i,k,j) = max(q(i+1,k,j),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i+1,k,j),qmin(i,k,j))
                    ! endif

                    ! if (v(i,k,j) > 0) then
                    !     qmax(i,k,j) = max(q(i,k,j-1),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i,k,j-1),qmin(i,k,j))
                    ! else if (v(i,k,j+1) < 0) then
                    !     qmax(i,k,j) = max(q(i,k,j+1),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i,k,j+1),qmin(i,k,j))
                    ! endif
                    
                    ! if (w(i,k,j) < 0 .and. k < kme) then
                    !     qmax(i,k,j) = max(q(i,k+1,j),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i,k+1,j),qmin(i,k,j))
                    ! else if (w(i,k,j) > 0 .and. k > kms) then
                    !     qmax(i,k,j) = max(q(i,k-1,j),qmax(i,k,j))
                    !     qmin(i,k,j) = min(q(i,k-1,j),qmin(i,k,j))
                    ! endif

                    !Store reused variables to minimize memory accesses
                    dz_t_i   = 1./dz(i,k,j)
                    ! jaco_rho_t_i = denom(i,k,j)
                    flux_x_up_1 = flux_x_up(i+1,k,j)
                    flux_y_up_1 = flux_y_up(i,k,j+1)
                    flux_z_up_1 = flux_z_up(i,k+1,j)
                    flux_x_up_0 = flux_x_up(i,k,j)
                    flux_y_up_0 = flux_y_up(i,k,j)
                    flux_z_up_0 = flux_z_up(i,k,j)
                    fx = flux_x(i,k,j)-flux_x_up_0; fx1 = flux_x(i+1,k,j)-flux_x_up_1
                    fy = flux_y(i,k,j)-flux_y_up_0; fy1 = flux_y(i,k,j+1)-flux_y_up_1
                    fz = flux_z(i,k,j)-flux_z_up_0; fz1 = flux_z(i,k+1,j)-flux_z_up_1
                                                
                    ! flux_x(i,k,j) = flux_x_up_0 - flux_x_up_1
                    ! flux_y(i,k,j) = flux_y_up_0 - flux_y_up_1
                    ! flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up(i,k,j)
                    !Compute concentration if upwind only was used
                    temp  = q0 - ((flux_x_up_1 - flux_x_up_0) + &
                                        (flux_y_up_1 - flux_y_up_0) + &
                                        (flux_z_up_1 - flux_z_up_0) * &
                                            dz_t_i)*denom(i,k,j)

                    flux_in = -(( (abs(fx1)-fx1) - (abs(fx)+fx) ) + &
                                ( (abs(fy1)-fy1) - (abs(fy)+fy) ) + &
                                ( (abs(fz1)-fz1) - (abs(fz)+fz) ) * &
                                    dz_t_i)*0.5*denom(i,k,j)
            
                    flux_out = (( (abs(fx1)+fx1) - (abs(fx)-fx) ) + &
                                ( (abs(fy1)+fy1) - (abs(fy)-fy) ) + &
                                ( (abs(fz1)+fz1) - (abs(fz)-fz) ) * &
                                    dz_t_i)*0.5*denom(i,k,j)
    
                    scale_in(i,k,j) = (qmax-temp)/(flux_in  + 0.000000001)
                    scale_out(i,k,j) = (temp-qmin)/(flux_out+ 0.000000001)
                enddo
            enddo
        enddo
        do concurrent (j = jts:jte+1, k = kms:kme, i = its:ite+1)

                    if (i >= its) then
                        flux_x(i,k,j) = flux_x(i,k,j) - flux_x_up(i,k,j)
                        if (flux_x(i,k,j) > 0) then
                            flux_x(i,k,j) = max(0.0,min(scale_in(i,k,j),scale_out(i-1,k,j),1.0))*flux_x(i,k,j) + flux_x_up(i,k,j)
                        else
                            flux_x(i,k,j) = max(0.0,min(scale_out(i,k,j),scale_in(i-1,k,j),1.0))*flux_x(i,k,j) + flux_x_up(i,k,j)
                        endif
                    end if
                    if (j >= jts) then
                        flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up(i,k,j)
                        if (flux_y(i,k,j) > 0) then
                            flux_y(i,k,j) = max(0.0,min(scale_in(i,k,j),scale_out(i,k,j-1),1.0))*flux_y(i,k,j) + flux_y_up(i,k,j)
                        else
                            flux_y(i,k,j) = max(0.0,min(scale_out(i,k,j),scale_in(i,k,j-1),1.0))*flux_y(i,k,j) + flux_y_up(i,k,j)
                        endif
                    endif
                    if (k > kms) then
                        flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up(i,k,j)
                        if (flux_z(i,k,j) > 0) then
                            flux_z(i,k,j) = max(0.0,min(scale_in(i,k,j),scale_out(i,k-1,j),1.0))*flux_z(i,k,j) + flux_z_up(i,k,j)
                        else
                            flux_z(i,k,j) = max(0.0,min(scale_out(i,k,j),scale_in(i,k-1,j),1.0))*flux_z(i,k,j) + flux_z_up(i,k,j)
                        endif
                    endif
            !     enddo
            ! enddo
        enddo

        ! do concurrent (j = jts:jte+1, k = kms:kme, i = its:ite+1)
        !     flux_x (i,k,j) = flux_x (i,k,j) - flux_x_up(i,k,j)
        !     if (flux_x(i,k,j) > 0) then
        !         flux_x(i,k,j) = min(scale_in(i,k,j),scale_out(i-1,k,j))*flux_x(i,k,j) + flux_x_up(i,k,j)
        !     else
        !         flux_x(i,k,j) = min(scale_out(i,k,j),scale_in(i-1,k,j))*flux_x(i,k,j) + flux_x_up(i,k,j)
        !     endif
        ! enddo
        ! do concurrent (j = jts:jte+1, k = kms:kme, i = its:ite+1)
        !     flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up(i,k,j)

        !     if (flux_y(i,k,j) > 0) then
        !         flux_y(i,k,j) = min(scale_in(i,k,j),scale_out(i,k,j-1))*flux_y(i,k,j) + flux_y_up(i,k,j)
        !     else
        !         flux_y(i,k,j) = min(scale_out(i,k,j),scale_in(i,k,j-1))*flux_y(i,k,j) + flux_y_up(i,k,j)
        !     endif
        ! enddo
        ! do concurrent (j = jts:jte+1, k = kms+1:kme, i = its:ite+1)
        !     flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up(i,k,j)
        !     if (flux_z(i,k,j) > 0) then
        !         flux_z(i,k,j) = min(scale_in(i,k,j),scale_out(i,k-1,j))*flux_z(i,k,j) + flux_z_up(i,k,j)
        !     else
        !         flux_z(i,k,j) = min(scale_out(i,k,j),scale_in(i,k-1,j))*flux_z(i,k,j) + flux_z_up(i,k,j)
        !     endif
        ! enddo

    end subroutine WRF_flux_corr

    subroutine upwind_flux3(q,u,v,w,flux_x,flux_z,flux_y,dz,denom)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in)      :: q, dz, denom
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in)    :: u
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in)    :: v
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w

        real, dimension(its-1:ite+2,kms:kme,jts-1:jte+2),intent(inout)          :: flux_x
        real, dimension(its-1:ite+2,kms:kme,jts-1:jte+2),intent(inout)          :: flux_y
        real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout)    :: flux_z
        
        real, dimension(ims:ime,  kms:kme,jms:jme) :: dumb_q
        integer :: i, j, k, bot, wes, sou
        real    :: tmp, abs_tmp
        
        !When using RK3, we may have a time step derived using a CFL constraint larger than 1
        !This means that our upwind advection here may be in violation of the CFL criterion,
        !since this does not take place within the RK3 scheme. Since the possible CFL constraints
        !under RK3 with the available advection orders are all < 2, we can simply do 2 upwind steps

        !Upwind fluxes for first (half of a )step already calculated in advection flux function -- apply them first here

        !Update intermediate concentration
        dumb_q = q
        ! $omp parallel default(none) &
        ! $omp shared(dumb_q, flux_x, flux_y, flux_z) &
        ! $omp shared(q, u, v, w, dz, denom) &
        ! $omp private(i,j,k, j_block, k_block, bot, wes, sou, tmp, abs_tmp) &
        ! $omp firstprivate(its, ite, jts, jte, kms,kme, ims, ime, jms, jme)
        ! $omp do schedule(static) collapse(2)
        do j = jts-1, jte+1
            do k = kms, kme
                do i = its-1, ite+1
                    dumb_q(i,k,j)  = q(i,k,j) - ((flux_x(i+1,k,j) - flux_x(i,k,j)) + &
                                                (flux_y(i,k,j+1) - flux_y(i,k,j)) + &
                                                (flux_z(i,k+1,j) - flux_z(i,k,j)) / &
                                                dz(i,k,j))*denom(i,k,j)
                enddo
            enddo
        enddo
        ! $omp end do

        !Now compute upwind fluxes after second step
        ! $omp do collapse(2)
        do j = jts-1, jte+2
            do k = kms, kme
                do i = its-1, ite+2
                    bot = max(k-1,kms)
                    wes = max(i-1,ims)
                    sou = max(j-1,jms)

                    tmp = u(i,k,j)
                    abs_tmp = ABS(tmp)
                    flux_x(i,k,j) = flux_x(i,k,j) + 0.5*((tmp + abs_tmp) * dumb_q(wes,k,j) + (tmp - abs_tmp) * dumb_q(i,k,j)) * 0.5
                    tmp = v(i,k,j)
                    abs_tmp = ABS(tmp)
                    flux_y(i,k,j) = flux_y(i,k,j) + 0.5*((tmp + abs_tmp) * dumb_q(i,k,sou) + (tmp - abs_tmp) * dumb_q(i,k,j)) * 0.5
                    tmp = w(i,bot,j)
                    abs_tmp = ABS(tmp)
                    flux_z(i,k,j) = flux_z(i,k,j) + 0.5*((tmp + abs_tmp) * dumb_q(i,bot,j) + (tmp - abs_tmp) * dumb_q(i,k,j)) * 0.5
                enddo
            enddo
        enddo
        ! $omp end do nowait
        ! $omp end parallel
        !Handle top and bottom boundaries for z here
        do j = jts-1,jte+1
            do i = its-1,ite+1
                flux_z(i,kme+1,j) = flux_z(i,kme+1,j) + 0.5*dumb_q(i,kme,j) * w(i,kme,j)
                flux_z(i,kms,j) = 0.0
            enddo
        enddo
                                        
    end subroutine upwind_flux3
end module adv_fluxcorr
