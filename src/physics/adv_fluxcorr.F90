!> ----------------------------------------------------------------------------
!!  A collection of flux correction schemes for advection
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_fluxcorr
    use domain_interface,  only: domain_t
    use mpi_f08, only: MPI_Wtime
    implicit none
    private
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte

    public :: WRF_flux_corr, init_fluxcorr, set_sign_arrays, clear_flux_sign_arrays
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

        !$acc enter data create(usign, vsign, wsign)
        !$acc data present(u, v, w, usign, vsign, wsign)
        !$acc parallel loop gang vector collapse(3)
        do j = jts-1,jte+1
            do k = kms,kme
                do i = its-1,ite+1
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

        !$acc end data

    end subroutine set_sign_arrays

    subroutine clear_flux_sign_arrays()
        implicit none

        !$acc exit data delete(usign, vsign, wsign)

    end subroutine clear_flux_sign_arrays

    subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,dz,denom)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, dz, denom
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: u
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: v
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_x
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_y
        real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout)          :: flux_z
        
        
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1)   :: scale_in, scale_out
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2)   :: flux_x_up
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2)   :: flux_y_up
        real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2)   :: flux_z_up

        real :: dz_t_i, fx, fx1, fy, fy1, fz, fz1, qmax, qmin, q_i, q_j, q_k, q0, temp, flux_in, flux_out, scale_in_cur, scale_out_cur
        integer :: i, j ,k
        real :: scale
        real :: flux_x_up_0, flux_y_up_0, flux_z_up_0, flux_x_up_1, flux_y_up_1, flux_z_up_1
        integer :: bot_in, bot_out, sou_in, sou_out, wes_in, wes_out, j_block, k_block, usign0, vsign0, wsign0

        !Initialize some internal variables
        scale = 1.0

        ! Get upwind fluxes

        !$acc data present(q,flux_x,flux_y,flux_z,dz,denom,usign,vsign,wsign) &
        !$acc create(scale_in,scale_out,flux_x_up,flux_y_up,flux_z_up)

        call upwind_flux3(q,u,v,w,flux_x_up,flux_z_up,flux_y_up,dz,denom)

        ! Next compute max and min possible fluxes        
        !$acc parallel loop gang vector collapse(3) async(1)
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

                    ! If the min and max for this point are 0, there is no flux limiter to calculate...
                    if ((abs(qmax)+abs(qmin)) == 0) then
                        scale_in(i,k,j) = 0.0
                        scale_out(i,k,j) = 0.0
                        cycle
                    endif
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
                    
                    ! Optimize flux calculations to reduce redundant memory accesses
                    flux_x_up_0 = flux_x_up(i,k,j)
                    flux_x_up_1 = flux_x_up(i+1,k,j)
                    flux_y_up_0 = flux_y_up(i,k,j)
                    flux_y_up_1 = flux_y_up(i,k,j+1)
                    flux_z_up_1 = flux_z_up(i,k+1,j)

                    fx = flux_x(i,k,j)-flux_x_up_0
                    fx1 = flux_x(i+1,k,j)-flux_x_up_1
                    fy = flux_y(i,k,j)-flux_y_up_0
                    fy1 = flux_y(i,k,j+1)-flux_y_up_1
                    fz1 = flux_z(i,k+1,j)-flux_z_up_1

                    if (k == kms) then
                        flux_z_up_0 = 0.0
                        fz = 0.0
                    else
                        flux_z_up_0 = flux_z_up(i,k,j)
                        fz = flux_z(i,k,j)-flux_z_up_0
                    endif

                    !Compute concentration if upwind only was used
                    temp  = q0 - ((flux_x_up_1 - flux_x_up_0) + &
                                        (flux_y_up_1 - flux_y_up_0) + &
                                        (flux_z_up_1 - flux_z_up_0) * &
                                            dz_t_i)*denom(i,k,j)
                    ! Optimize flux_in and flux_out calculations using fewer operations
                    ! flux_in = positive inflow, flux_out = positive outflow
                    flux_in =  ((max(0.0,-fx1) + max(0.0,fx)) + &
                                (max(0.0,-fy1) + max(0.0,fy)) + &
                                (max(0.0,-fz1) + max(0.0,fz)) * &
                                    dz_t_i)*denom(i,k,j)
            
                    flux_out = ((max(0.0,fx1) + max(0.0,-fx)) + &
                                (max(0.0,fy1) + max(0.0,-fy)) + &
                                (max(0.0,fz1) + max(0.0,-fz)) * &
                                    dz_t_i)*denom(i,k,j)
    
                    scale_in(i,k,j) = (qmax-temp)/(flux_in  + 0.000000001)
                    scale_out(i,k,j) = (temp-qmin)/(flux_out+ 0.000000001)
                enddo
            enddo
        enddo

        !$acc parallel async(2) wait(1)
        !$acc loop gang vector collapse(3)
        do j = jts, jte+1 
            do k = kms, kme
                do i = its, ite+1
                    scale_in_cur = scale_in(i,k,j)
                    scale_out_cur = scale_out(i,k,j)
                    
                    ! X-direction flux correction
                    flux_x(i,k,j) = flux_x(i,k,j) - flux_x_up(i,k,j)
                    scale = merge(max(0.0,min(scale_in_cur,scale_out(i-1,k,j),1.0)), &
                                  max(0.0,min(scale_out_cur,scale_in(i-1,k,j),1.0)), &
                                  flux_x(i,k,j) > 0.0)
                    flux_x(i,k,j) = scale*flux_x(i,k,j) + flux_x_up(i,k,j)

                    ! Y-direction flux correction
                    flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up(i,k,j)
                    scale = merge(max(0.0,min(scale_in_cur,scale_out(i,k,j-1),1.0)), &
                                  max(0.0,min(scale_out_cur,scale_in(i,k,j-1),1.0)), &
                                  flux_y(i,k,j) > 0.0)
                    flux_y(i,k,j) = scale*flux_y(i,k,j) + flux_y_up(i,k,j)

                    ! Z-direction flux correction (only for k > kms)
                    if (k > kms) then
                        flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up(i,k,j)
                        scale = merge(max(0.0,min(scale_in_cur,scale_out(i,k-1,j),1.0)), &
                                      max(0.0,min(scale_out_cur,scale_in(i,k-1,j),1.0)), &
                                      flux_z(i,k,j) > 0.0)
                        flux_z(i,k,j) = scale*flux_z(i,k,j) + flux_z_up(i,k,j)
                    endif
                enddo
            enddo
        enddo

        !$acc loop gang vector collapse(2)
        do j = jts, jte+1 
            do i = its, ite+1
                flux_z(i,kme+1,j) = flux_z(i,kme+1,j) - flux_z_up(i,kme+1,j)
                scale = merge(max(0.0,min(scale_in(i,kme,j),scale_out(i,kme,j),1.0)), &
                              1.0, abs(flux_z(i,kme+1,j)) > 0.0)
                flux_z(i,kme+1,j) = scale*flux_z(i,kme+1,j) + flux_z_up(i,kme+1,j)
            enddo
        enddo
        !$acc end parallel
        !$acc end data
        !$acc wait(2)
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
        real    :: tmp, abs_tmp, flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val, q_cache(4)
        
        !When using RK3, we may have a time step derived using a CFL constraint larger than 1
        !This means that our upwind advection here may be in violation of the CFL criterion,
        !since this does not take place within the RK3 scheme. Since the possible CFL constraints
        !under RK3 with the available advection orders are all < 2, we can simply do 2 upwind steps

        !Upwind fluxes for first (half of a )step already calculated in advection flux function -- apply them first here

        !Update intermediate concentration

        !$acc data present(q,u,v,w,flux_x,flux_y,flux_z,dz,denom) &
        !$acc create(dumb_q)


        !Now compute upwind fluxes after second step
        !!$acc parallel 
        !$acc parallel loop gang vector collapse(3) private(q_cache, tmp) async(2) wait(1)
        do j = jts-1, jte+2
            do k = kms, kme+1
                do i = its-1, ite+2

                    if (k <= kme) then
                        q_cache(1) = q(i,k,j-1)
                        if (k > kms) q_cache(2) = q(i,k-1,j)
                        q_cache(3) = q(i-1,k,j)
                        q_cache(4) = q(i,k,j)

                        ! X-direction flux - optimize using single calculation
                        tmp = u(i,k,j)
                        flux_x(i,k,j) = 0.25 * (tmp * (q_cache(3) + q_cache(4)) + &
                                                                    abs(tmp) * (q_cache(3) - q_cache(4)))

                        ! Y-direction flux - optimize using single calculation
                        tmp = v(i,k,j)
                        flux_y(i,k,j) = 0.25 * (tmp * (q_cache(1) + q_cache(4)) + &
                                                                    abs(tmp) * (q_cache(1) - q_cache(4)))

                        if (k > kms) then
                            ! Z-direction flux - optimize using single calculation
                            tmp = w(i,k-1,j)
                            flux_z(i,k,j) = 0.25 * (tmp * (q_cache(2) + q_cache(4)) + &
                                                                        abs(tmp) * (q_cache(2) - q_cache(4)))
                        endif
                    else
                        flux_z(i,k,j) = 0.5*q(i,k-1,j) * w(i,k-1,j)
                    endif
                enddo
            enddo
        enddo

        
        !$acc parallel loop gang vector collapse(3) private(flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val) async(1)
        do j = jms, jme
            do k = kms, kme
            do i = ims, ime

                if (i > ite+1 .or. i < its-1 .or. j > jte+1 .or. j < jts-1) then
                    dumb_q(i,k,j) = q(i,k,j)
                else
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
                    dumb_q(i,k,j) = q(i,k,j) - (flux_diff_x + flux_diff_y + flux_diff_z / dz_val) * denom_val
                endif
            enddo
            enddo
        enddo


        !Now compute upwind fluxes after second step
        !!$acc parallel 
        !$acc parallel loop gang vector collapse(3) private(q_cache, tmp) async(2) wait(1)
        do j = jts-1, jte+2
            do k = kms, kme+1
                do i = its-1, ite+2

                    if (k <= kme) then
                        q_cache(1) = dumb_q(i,k,j-1)
                        if (k > kms) q_cache(2) = dumb_q(i,k-1,j)
                        q_cache(3) = dumb_q(i-1,k,j)
                        q_cache(4) = dumb_q(i,k,j)

                        ! X-direction flux - optimize using single calculation
                        tmp = u(i,k,j)
                        flux_x(i,k,j) = flux_x(i,k,j) + 0.25 * (tmp * (q_cache(3) + q_cache(4)) + &
                                                                    abs(tmp) * (q_cache(3) - q_cache(4)))

                        ! Y-direction flux - optimize using single calculation
                        tmp = v(i,k,j)
                        flux_y(i,k,j) = flux_y(i,k,j) + 0.25 * (tmp * (q_cache(1) + q_cache(4)) + &
                                                                    abs(tmp) * (q_cache(1) - q_cache(4)))

                        if (k > kms) then
                            ! Z-direction flux - optimize using single calculation
                            tmp = w(i,k-1,j)
                            flux_z(i,k,j) = flux_z(i,k,j) + 0.25 * (tmp * (q_cache(2) + q_cache(4)) + &
                                                                        abs(tmp) * (q_cache(2) - q_cache(4)))
                        endif
                    else
                        flux_z(i,k,j) = flux_z(i,k,j) + 0.5*dumb_q(i,k-1,j) * w(i,k-1,j)
                    endif
                enddo
            enddo
        enddo

        ! !Handle top and bottom boundaries for z here
        ! !$acc parallel loop gang vector collapse(2)
        ! do j = jts-1,jte+1
        !     do i = its-1,ite+1
        !         flux_z(i,kme+1,j) = flux_z(i,kme+1,j) + 0.5*dumb_q(i,kme,j) * w(i,kme,j)
        !         ! flux_z(i,kms,j) = 0.0
        !     enddo
        ! enddo

        !!$acc end parallel
        !$acc end data
        !$acc wait(2)
    end subroutine upwind_flux3
end module adv_fluxcorr
