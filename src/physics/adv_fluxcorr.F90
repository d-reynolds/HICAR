!> ----------------------------------------------------------------------------
!!  A collection of flux correction schemes for advection
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_fluxcorr
    use domain_interface,  only: domain_t
    implicit none
    private
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte
    integer, parameter :: init_async_id = 63853 !random value that hopefully never comes up elsewhere

    public :: WRF_flux_corr, init_fluxcorr, set_sign_arrays, compute_upwind_fluxes_async
    public :: init_fluxcorr_fm, set_sign_arrays_fm, WRF_flux_corr_fm

    integer, allocatable, dimension(:,:,:)   :: usign, vsign, wsign
    real,    allocatable, dimension(:,:,:)   :: scale_in, scale_out, qmax, qmin
    real,    allocatable, dimension(:,:,:)   :: flux_x_up
    real,    allocatable, dimension(:,:,:)   :: flux_y_up
    real,    allocatable, dimension(:,:,:)   :: flux_z_up
    real,    allocatable, dimension(:,:,:)   :: dumb_q

    ! Fine-mesh flux correction arrays
    integer, allocatable, dimension(:,:,:)   :: usign_fm, vsign_fm
    real,    allocatable, dimension(:,:,:)   :: scale_in_fm, scale_out_fm
    real,    allocatable, dimension(:,:,:)   :: flux_x_up_fm, flux_y_up_fm
    real,    allocatable, dimension(:,:,:)   :: dumb_q_fm

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

        if (allocated(usign)) then
            !$acc exit data delete(usign, vsign, wsign, scale_in, scale_out, qmax, qmin, flux_x_up, flux_y_up, flux_z_up, dumb_q)
            deallocate(usign)
        endif
        if (allocated(vsign)) deallocate(vsign)
        if (allocated(wsign)) deallocate(wsign)
        if (allocated(scale_in)) deallocate(scale_in)
        if (allocated(scale_out)) deallocate(scale_out)
        if (allocated(qmax)) deallocate(qmax)
        if (allocated(qmin)) deallocate(qmin)
        if (allocated(flux_x_up)) deallocate(flux_x_up)
        if (allocated(flux_y_up)) deallocate(flux_y_up)
        if (allocated(flux_z_up)) deallocate(flux_z_up)
        if (allocated(dumb_q)) deallocate(dumb_q)

        allocate(usign(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(vsign(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(wsign(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(scale_in(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(scale_out(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(qmax(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(qmin(its-1:ite+1,kms:kme,jts-1:jte+1))
        allocate(flux_x_up(ims:ime,kms:kme,jms:jme))
        allocate(flux_y_up(ims:ime,kms:kme,jms:jme))
        allocate(flux_z_up(ims:ime,kms:kme+1,jms:jme))
        allocate(dumb_q(ims:ime,kms:kme,jms:jme))

        !$acc enter data create(usign, vsign, wsign, scale_in, scale_out, qmax, qmin, flux_x_up, flux_y_up, flux_z_up, dumb_q)
    end subroutine init_fluxcorr

    subroutine set_sign_arrays(u,v,w)
        implicit none
        real, dimension(its-2:ite+3,  kms:kme,jts-2:jte+3),  intent(in) :: w
        real, dimension(its-2:ite+3,  kms:kme,jts-2:jte+3),  intent(in) :: u
        real, dimension(its-2:ite+3,  kms:kme,jts-2:jte+3),  intent(in) :: v

        integer :: i, j, k

        !$acc data present(u, v, w, usign, vsign, wsign)
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256)
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

    ! subroutine clear_flux_sign_arrays()
    !     implicit none

    !     !$acc exit data delete(usign, vsign, wsign, scale_in, scale_out, qmax, qmin, flux_x_up, flux_y_up, flux_z_up, dumb_q)

    ! end subroutine clear_flux_sign_arrays

    ! subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,dz,denom)
    !     implicit none
    !     real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, dz, denom
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: u
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: v
        
    !     real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_x
    !     real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_y
    !     real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout)          :: flux_z
        
        

    !     real :: dz_t_i, fx, fx1, fy, fy1, fz, fz1, q_i, q_j, q_k, q0, temp, flux_in, flux_out, scale_in_cur, scale_out_cur
    !     integer :: i, j ,k
    !     real :: scale
    !     real :: flux_x_up_0, flux_y_up_0, flux_z_up_0, flux_x_up_1, flux_y_up_1, flux_z_up_1
    !     integer :: bot_in, bot_out, sou_in, sou_out, wes_in, wes_out, j_block, k_block, usign0, vsign0, wsign0

    !     !Initialize some internal variables
    !     scale = 1.0

    !     ! Get upwind fluxes

    !     call upwind_flux3(q,u,v,w,flux_x_up,flux_z_up,flux_y_up,dz,denom)

    !     !$acc data present(usign, vsign, wsign, flux_x, flux_y, flux_z, scale_in, scale_out, qmax, qmin, q, denom, dz, flux_x_up, flux_y_up, flux_z_up) async(1)

    !     ! Next compute max and min possible fluxes        
    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(1)
    !     do j = jts-1, jte+1 
    !         do k = kms, kme
    !             do i = its-1, ite+1
    !                 usign0 = usign(i,k,j)
    !                 vsign0 = vsign(i,k,j)
    !                 wsign0 = wsign(i,k,j)

    !                 q0 = q(i,k,j)
    !                 q_i = q(i+usign0,k,j)
    !                 q_j = q(i,k,j+vsign0)
    !                 q_k = q(i,k+wsign0,j)

    !                 qmax(i,k,j) = max(q0, q_i, q_j, q_k)
    !                 qmin(i,k,j) = min(q0, q_i, q_j, q_k)
    !             enddo
    !         enddo
    !     enddo
    !     !$acc wait(10)

    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(1)
    !     do j = jts, jte+1 
    !         do k = kms, kme
    !             do i = its, ite+1
    !                 ! If the min and max for this point are 0, there is no flux limiter to calculate...
    !                 if ((abs(qmax(i,k,j))+abs(qmin(i,k,j))) == 0) then
    !                     scale_in(i,k,j) = 0.0
    !                     scale_out(i,k,j) = 0.0
    !                     cycle
    !                 endif

    !                 !Store reused variables to minimize memory accesses
    !                 dz_t_i   = 1./dz(i,k,j)
                    
    !                 ! Optimize flux calculations to reduce redundant memory accesses
    !                 flux_x_up_0 = flux_x_up(i,k,j)
    !                 flux_x_up_1 = flux_x_up(i+1,k,j)
    !                 flux_y_up_0 = flux_y_up(i,k,j)
    !                 flux_y_up_1 = flux_y_up(i,k,j+1)
    !                 flux_z_up_1 = flux_z_up(i,k+1,j)

    !                 fx = flux_x(i,k,j)-flux_x_up_0
    !                 fx1 = flux_x(i+1,k,j)-flux_x_up_1
    !                 fy = flux_y(i,k,j)-flux_y_up_0
    !                 fy1 = flux_y(i,k,j+1)-flux_y_up_1
    !                 fz1 = flux_z(i,k+1,j)-flux_z_up_1

    !                 if (k == kms) then
    !                     flux_z_up_0 = 0.0
    !                     fz = 0.0
    !                 else
    !                     flux_z_up_0 = flux_z_up(i,k,j)
    !                     fz = flux_z(i,k,j)-flux_z_up_0
    !                 endif

    !                 !Compute concentration if upwind only was used
    !                 temp  = q(i,k,j) - ((flux_x_up_1 - flux_x_up_0) + &
    !                                     (flux_y_up_1 - flux_y_up_0) + &
    !                                     (flux_z_up_1 - flux_z_up_0) * &
    !                                         dz_t_i)*denom(i,k,j)
    !                 ! Optimize flux_in and flux_out calculations using fewer operations
    !                 ! flux_in = positive inflow, flux_out = positive outflow
    !                 flux_in =  ((max(0.0,-fx1) + max(0.0,fx)) + &
    !                             (max(0.0,-fy1) + max(0.0,fy)) + &
    !                             (max(0.0,-fz1) + max(0.0,fz)) * &
    !                                 dz_t_i)*denom(i,k,j)
            
    !                 flux_out = ((max(0.0,fx1) + max(0.0,-fx)) + &
    !                             (max(0.0,fy1) + max(0.0,-fy)) + &
    !                             (max(0.0,fz1) + max(0.0,-fz)) * &
    !                                 dz_t_i)*denom(i,k,j)
    
    !                 scale_in(i,k,j) = (qmax(i,k,j)-temp)/(flux_in  + 0.000000001)
    !                 scale_out(i,k,j) = (temp-qmin(i,k,j))/(flux_out+ 0.000000001)
    !             enddo
    !         enddo
    !     enddo

    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(2) wait(1)
    !     do j = jts, jte+1 
    !         do k = kms, kme
    !             do i = its, ite+1
    !                 scale_in_cur = scale_in(i,k,j)
    !                 scale_out_cur = scale_out(i,k,j)
                    
    !                 ! X-direction flux correction
    !                 flux_x(i,k,j) = flux_x(i,k,j) - flux_x_up(i,k,j)
    !                 scale = merge(max(0.0,min(scale_in_cur,scale_out(i-1,k,j),1.0)), &
    !                               max(0.0,min(scale_out_cur,scale_in(i-1,k,j),1.0)), &
    !                               flux_x(i,k,j) > 0.0)
    !                 flux_x(i,k,j) = scale*flux_x(i,k,j) + flux_x_up(i,k,j)

    !                 ! Y-direction flux correction
    !                 flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up(i,k,j)
    !                 scale = merge(max(0.0,min(scale_in_cur,scale_out(i,k,j-1),1.0)), &
    !                               max(0.0,min(scale_out_cur,scale_in(i,k,j-1),1.0)), &
    !                               flux_y(i,k,j) > 0.0)
    !                 flux_y(i,k,j) = scale*flux_y(i,k,j) + flux_y_up(i,k,j)

    !                 ! Z-direction flux correction (only for k > kms)
    !                 if (k > kms) then
    !                     flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up(i,k,j)
    !                     scale = merge(max(0.0,min(scale_in_cur,scale_out(i,k-1,j),1.0)), &
    !                                   max(0.0,min(scale_out_cur,scale_in(i,k-1,j),1.0)), &
    !                                   flux_z(i,k,j) > 0.0)
    !                     flux_z(i,k,j) = scale*flux_z(i,k,j) + flux_z_up(i,k,j)
    !                 endif
    !             enddo
    !         enddo
    !     enddo

    !     !$acc parallel loop gang vector collapse(2) async(3) wait(1)
    !     do j = jts, jte+1 
    !         do i = its, ite+1
    !             flux_z(i,kme+1,j) = flux_z(i,kme+1,j) - flux_z_up(i,kme+1,j)
    !             scale = merge(max(0.0,min(scale_in(i,kme,j),scale_out(i,kme,j),1.0)), &
    !                           1.0, flux_z(i,kme+1,j) > 0.0)
    !             flux_z(i,kme+1,j) = scale*flux_z(i,kme+1,j) + flux_z_up(i,kme+1,j)
    !         enddo
    !     enddo
    !     !$acc end data
    !     !$acc wait(2,3)

    ! end subroutine WRF_flux_corr

    ! subroutine upwind_flux3(q,u,v,w,flux_x,flux_z,flux_y,dz,denom)
    !     implicit none
    !     real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in)      :: q, dz, denom
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in)    :: u
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in)    :: v
    !     real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+2),  intent(in) :: w

    !     real, dimension(its-1:ite+2,kms:kme,jts-1:jte+2),intent(inout)          :: flux_x
    !     real, dimension(its-1:ite+2,kms:kme,jts-1:jte+2),intent(inout)          :: flux_y
    !     real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout)    :: flux_z
        
    !     real, dimension(ims:ime,  kms:kme,jms:jme) :: dumb_q
    !     integer :: i, j, k, bot, wes, sou
    !     real    :: tmp, abs_tmp, flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val, q_cache(4)
        
    !     !When using RK3, we may have a time step derived using a CFL constraint larger than 1
    !     !This means that our upwind advection here may be in violation of the CFL criterion,
    !     !since this does not take place within the RK3 scheme. Since the possible CFL constraints
    !     !under RK3 with the available advection orders are all < 2, we can simply do 2 upwind steps

    !     !Upwind fluxes for first (half of a )step already calculated in advection flux function -- apply them first here

    !     !Update intermediate concentration

    !     !$acc data present(q,u,v,w,flux_x,flux_y,flux_z,dz,denom) &
    !     !$acc create(dumb_q) async(10)


    !     !Now compute upwind fluxes after second step

    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) private(q_cache, tmp) async(10)
    !     do j = jts-1, jte+2
    !         do k = kms, kme+1
    !             do i = its-1, ite+2

    !                 if (k <= kme) then
    !                     q_cache(1) = q(i,k,j-1)
    !                     if (k > kms) q_cache(2) = q(i,k-1,j)
    !                     q_cache(3) = q(i-1,k,j)
    !                     q_cache(4) = q(i,k,j)

    !                     ! X-direction flux - optimize using single calculation
    !                     tmp = u(i,k,j)
    !                     flux_x(i,k,j) = 0.25 * (tmp * (q_cache(3) + q_cache(4)) + &
    !                                                                 abs(tmp) * (q_cache(3) - q_cache(4)))

    !                     ! Y-direction flux - optimize using single calculation
    !                     tmp = v(i,k,j)
    !                     flux_y(i,k,j) = 0.25 * (tmp * (q_cache(1) + q_cache(4)) + &
    !                                                                 abs(tmp) * (q_cache(1) - q_cache(4)))

    !                     if (k > kms) then
    !                         ! Z-direction flux - optimize using single calculation
    !                         tmp = w(i,k-1,j)
    !                         flux_z(i,k,j) = 0.25 * (tmp * (q_cache(2) + q_cache(4)) + &
    !                                                                     abs(tmp) * (q_cache(2) - q_cache(4)))
    !                     endif
    !                 else
    !                     flux_z(i,k,j) = 0.5*q(i,k-1,j) * w(i,k-1,j)
    !                 endif
    !             enddo
    !         enddo
    !     enddo

        
    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) private(flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val) async(10)
    !     do j = jms, jme
    !         do k = kms, kme
    !         do i = ims, ime

    !             if (i > ite+1 .or. i < its-1 .or. j > jte+1 .or. j < jts-1) then
    !                 dumb_q(i,k,j) = q(i,k,j)
    !             else
    !                 ! Cache frequently accessed values to reduce memory traffic
    !                 flux_diff_x = flux_x(i+1,k,j) - flux_x(i,k,j)
    !                 flux_diff_y = flux_y(i,k,j+1) - flux_y(i,k,j)

    !                 !if k is the bottom level, we expect no flux in/out of bottom boundary. In fact,
    !                 !advection code is optimized such that no flux is ever calculated...
    !                 if (k==kms) then
    !                     flux_diff_z = flux_z(i,k+1,j)
    !                 else
    !                     flux_diff_z = flux_z(i,k+1,j) - flux_z(i,k,j)
    !                 endif
    !                 denom_val = denom(i,k,j)
    !                 dz_val = dz(i,k,j)
                    
    !                 ! Perform advection calculation with cached values
    !                 dumb_q(i,k,j) = q(i,k,j) - (flux_diff_x + flux_diff_y + flux_diff_z / dz_val) * denom_val
    !             endif
    !         enddo
    !         enddo
    !     enddo


    !     !Now compute upwind fluxes after second step

    !     !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) private(q_cache, tmp) async(10)
    !     do j = jts-1, jte+2
    !         do k = kms, kme+1
    !             do i = its-1, ite+2

    !                 if (k <= kme) then
    !                     q_cache(1) = dumb_q(i,k,j-1)
    !                     if (k > kms) q_cache(2) = dumb_q(i,k-1,j)
    !                     q_cache(3) = dumb_q(i-1,k,j)
    !                     q_cache(4) = dumb_q(i,k,j)

    !                     ! X-direction flux - optimize using single calculation
    !                     tmp = u(i,k,j)
    !                     flux_x(i,k,j) = flux_x(i,k,j) + 0.25 * (tmp * (q_cache(3) + q_cache(4)) + &
    !                                                                 abs(tmp) * (q_cache(3) - q_cache(4)))

    !                     ! Y-direction flux - optimize using single calculation
    !                     tmp = v(i,k,j)
    !                     flux_y(i,k,j) = flux_y(i,k,j) + 0.25 * (tmp * (q_cache(1) + q_cache(4)) + &
    !                                                                 abs(tmp) * (q_cache(1) - q_cache(4)))

    !                     if (k > kms) then
    !                         ! Z-direction flux - optimize using single calculation
    !                         tmp = w(i,k-1,j)
    !                         flux_z(i,k,j) = flux_z(i,k,j) + 0.25 * (tmp * (q_cache(2) + q_cache(4)) + &
    !                                                                     abs(tmp) * (q_cache(2) - q_cache(4)))
    !                     endif
    !                 else
    !                     flux_z(i,k,j) = flux_z(i,k,j) + 0.5*dumb_q(i,k-1,j) * w(i,k-1,j)
    !                 endif
    !             enddo
    !         enddo
    !     enddo

    !     !$acc end data
    !     ! !$acc wait(10)
    ! end subroutine upwind_flux3

    !>------------------------------------------------------------
    !! ASYNC UPWIND FLUX COMPUTATION
    !! 
    !! Computes two-step upwind fluxes asynchronously on specified stream
    !! Stores results in module-level flux_x_up, flux_y_up, flux_z_up arrays
    !! 
    !! This can be called BEFORE flux3 to overlap computation
    !!------------------------------------------------------------
    subroutine compute_upwind_fluxes_async(q,u,v,w,dz,denom,async_id)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, dz, denom
        real, dimension(its-2:ite+3,  kms:kme,jts-2:jte+3),  intent(in) :: w, u, v
        integer, intent(in) :: async_id
        
        real :: q0, u_val, v_val, w_val, abs_u, abs_v, abs_w
        real :: flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val
        integer :: i, j, k

        !$acc data present(u, v, w, q, denom, dz, flux_x_up, flux_y_up, flux_z_up, dumb_q) async(async_id)

        ! ==========================================
        ! STEP 1: First half-step upwind fluxes
        ! Split into 3 branchless kernels to keep the dominant interior loop vectorisable.
        ! 1a: interior k=kms+1..kme (x, y, z fluxes — all three directions)
        ! 1b: bottom k=kms (x, y only; flux_z_up at kms left untouched per original)
        ! 1c: top k=kme+1 (only flux_z_up = 0.5*q(k-1)*w(k-1))
        ! ==========================================

        ! 1a: interior — branchless 3D loop, all three directions
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(async_id) &
        !$acc private(q0, u_val, v_val, w_val, abs_u, abs_v, abs_w)
        do j = jts-2, jte+3
            do k = kms+1, kme
                do i = its-2, ite+3
                    q0 = q(i,k,j)
                    u_val = u(i,k,j); abs_u = abs(u_val)
                    v_val = v(i,k,j); abs_v = abs(v_val)
                    w_val = w(i,k-1,j); abs_w = abs(w_val)
                    flux_x_up(i,k,j) = 0.25 * (u_val * (q(i-1,k,j) + q0) + abs_u * (q(i-1,k,j) - q0))
                    flux_y_up(i,k,j) = 0.25 * (v_val * (q(i,k,j-1) + q0) + abs_v * (q(i,k,j-1) - q0))
                    flux_z_up(i,k,j) = 0.25 * (w_val * (q(i,k-1,j) + q0) + abs_w * (q(i,k-1,j) - q0))
                enddo
            enddo
        enddo

        ! 1b: bottom k=kms — only x and y; flux_z_up(i,kms,j) left at its initial value
        !$acc parallel loop gang vector tile(32, 8) async(async_id) &
        !$acc private(q0, u_val, v_val, abs_u, abs_v)
        do j = jts-2, jte+3
            do i = its-2, ite+3
                q0 = q(i,kms,j)
                u_val = u(i,kms,j); abs_u = abs(u_val)
                v_val = v(i,kms,j); abs_v = abs(v_val)
                flux_x_up(i,kms,j) = 0.25 * (u_val * (q(i-1,kms,j) + q0) + abs_u * (q(i-1,kms,j) - q0))
                flux_y_up(i,kms,j) = 0.25 * (v_val * (q(i,kms,j-1) + q0) + abs_v * (q(i,kms,j-1) - q0))
            enddo
        enddo

        ! 1c: top k=kme+1 — only flux_z_up = 0.5*q(kme)*w(kme)
        !$acc parallel loop gang vector tile(32, 8) async(async_id)
        do j = jts-2, jte+3
            do i = its-2, ite+3
                flux_z_up(i,kme+1,j) = 0.5 * q(i,kme,j) * w(i,kme,j)
            enddo
        enddo

        ! ==========================================
        ! STEP 2: Intermediate concentration. Split bottom k=kms slab (different flux_diff_z formula).
        ! ==========================================

        ! 2a: bottom k=kms — flux_diff_z = flux_z_up(kms+1,j) (no flux_z_up at kms in the original)
        !$acc parallel loop gang vector tile(32, 8) async(async_id) &
        !$acc private(flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val)
        do j = jts-2, jte+2
            do i = its-2, ite+2
                flux_diff_x = flux_x_up(i+1,kms,j) - flux_x_up(i,kms,j)
                flux_diff_y = flux_y_up(i,kms,j+1) - flux_y_up(i,kms,j)
                flux_diff_z = flux_z_up(i,kms+1,j)
                denom_val = denom(i,kms,j); dz_val = dz(i,kms,j)
                dumb_q(i,kms,j) = q(i,kms,j) - (flux_diff_x + flux_diff_y + flux_diff_z / dz_val) * denom_val
            enddo
        enddo

        ! 2b: interior k=kms+1..kme — branchless 3D
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(async_id) &
        !$acc private(flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val)
        do j = jts-2, jte+2
            do k = kms+1, kme
                do i = its-2, ite+2
                    flux_diff_x = flux_x_up(i+1,k,j) - flux_x_up(i,k,j)
                    flux_diff_y = flux_y_up(i,k,j+1) - flux_y_up(i,k,j)
                    flux_diff_z = flux_z_up(i,k+1,j) - flux_z_up(i,k,j)
                    denom_val = denom(i,k,j); dz_val = dz(i,k,j)
                    dumb_q(i,k,j) = q(i,k,j) - (flux_diff_x + flux_diff_y + flux_diff_z / dz_val) * denom_val
                enddo
            enddo
        enddo

        ! ==========================================
        ! STEP 3: Second half-step upwind fluxes (accumulate). Mirror of step 1.
        ! ==========================================

        ! 3a: interior k=kms+1..kme
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(async_id) &
        !$acc private(q0, u_val, v_val, w_val, abs_u, abs_v, abs_w)
        do j = jts-1, jte+2
            do k = kms+1, kme
                do i = its-1, ite+2
                    q0 = dumb_q(i,k,j)
                    u_val = u(i,k,j); abs_u = abs(u_val)
                    v_val = v(i,k,j); abs_v = abs(v_val)
                    w_val = w(i,k-1,j); abs_w = abs(w_val)
                    flux_x_up(i,k,j) = flux_x_up(i,k,j) + 0.25 * (u_val * (dumb_q(i-1,k,j) + q0) + abs_u * (dumb_q(i-1,k,j) - q0))
                    flux_y_up(i,k,j) = flux_y_up(i,k,j) + 0.25 * (v_val * (dumb_q(i,k,j-1) + q0) + abs_v * (dumb_q(i,k,j-1) - q0))
                    flux_z_up(i,k,j) = flux_z_up(i,k,j) + 0.25 * (w_val * (dumb_q(i,k-1,j) + q0) + abs_w * (dumb_q(i,k-1,j) - q0))
                enddo
            enddo
        enddo

        ! 3b: bottom k=kms — x and y only
        !$acc parallel loop gang vector tile(32, 8) async(async_id) &
        !$acc private(q0, u_val, v_val, abs_u, abs_v)
        do j = jts-1, jte+2
            do i = its-1, ite+2
                q0 = dumb_q(i,kms,j)
                u_val = u(i,kms,j); abs_u = abs(u_val)
                v_val = v(i,kms,j); abs_v = abs(v_val)
                flux_x_up(i,kms,j) = flux_x_up(i,kms,j) + 0.25 * (u_val * (dumb_q(i-1,kms,j) + q0) + abs_u * (dumb_q(i-1,kms,j) - q0))
                flux_y_up(i,kms,j) = flux_y_up(i,kms,j) + 0.25 * (v_val * (dumb_q(i,kms,j-1) + q0) + abs_v * (dumb_q(i,kms,j-1) - q0))
            enddo
        enddo

        ! 3c: top k=kme+1 — z only, accumulate
        !$acc parallel loop gang vector tile(32, 8) async(async_id)
        do j = jts-1, jte+2
            do i = its-1, ite+2
                flux_z_up(i,kme+1,j) = flux_z_up(i,kme+1,j) + 0.5 * dumb_q(i,kme,j) * w(i,kme,j)
            enddo
        enddo

        !$acc end data
        ! Note: No wait here - async computation continues

    end subroutine compute_upwind_fluxes_async

    !>------------------------------------------------------------
    !! ASYNC VERSION: Uses pre-computed module-level upwind fluxes
    !! 
    !! REQUIRES: compute_upwind_fluxes_async() called beforehand
    !! Waits on async_id, then applies flux corrections
    !!------------------------------------------------------------
    subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,dz,denom,async_id)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q, dz, denom
        real, dimension(its-2:ite+3,  kms:kme,jts-2:jte+3),  intent(in) :: w, u, v
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+2),intent(inout) :: flux_x, flux_y
        real, dimension(its-1:ite+2,kms:kme+1,jts-1:jte+2),intent(inout) :: flux_z
        integer, intent(in) :: async_id
        
        real :: dz_t_i, q0, q_i, q_j, q_k, temp, flux_in, flux_out
        real :: scale, scale_in_val, scale_out_val, qmax_local, qmin_local
        real :: fx, fx1, fy, fy1, fz, fz1
        real :: flux_x_up_0, flux_x_up_1, flux_y_up_0, flux_y_up_1, flux_z_up_0, flux_z_up_1
        integer :: i, j, k, usign_val, vsign_val, wsign_val

        ! Wait for async upwind flux computation to complete
        !$acc wait(async_id)

        !$acc data present(u, v, w, flux_x, flux_y, flux_z, q, denom, dz, scale_in, scale_out, flux_x_up, flux_y_up, flux_z_up, usign, vsign, wsign) async(1000)
        
        ! ==========================================
        ! STEP 1: scale_in/scale_out — split into 2 kernels (bottom k=kms vs interior)
        ! to keep the dominant interior branchless. The k==kms branch sets flux_z_up_0=0
        ! and fz=0; interior reads flux_z_up(i,k,j) and computes fz normally.
        ! Note: the early `cycle` when |qmax|+|qmin|==0 is preserved as is — it's a
        ! per-cell short-circuit, not a row-category branch, so it doesn't inhibit
        ! vectorisation more than absolutely necessary.
        ! ==========================================

        ! 1a: bottom k=kms (flux_z_up_0 = 0, fz = 0)
        !$acc parallel loop gang vector tile(32, 8) async(1000) &
        !$acc private(q0, q_i, q_j, q_k, qmax_local, qmin_local, dz_t_i, temp, &
        !$acc         flux_x_up_0, flux_x_up_1, flux_y_up_0, flux_y_up_1, flux_z_up_1, &
        !$acc         fx, fx1, fy, fy1, fz, fz1, flux_in, flux_out)
        do j = jts-1, jte+1
            do i = its-1, ite+1
                q0 = q(i,kms,j)
                q_i = q(i+usign(i,kms,j),kms,j)
                q_j = q(i,kms,j+vsign(i,kms,j))
                q_k = q(i,kms+wsign(i,kms,j),j)
                qmax_local = max(q0, q_i, q_j, q_k)
                qmin_local = min(q0, q_i, q_j, q_k)
                if ((abs(qmax_local) + abs(qmin_local)) == 0.0) then
                    scale_in(i,kms,j) = 0.0
                    scale_out(i,kms,j) = 0.0
                    cycle
                endif
                dz_t_i = 1.0 / dz(i,kms,j)
                flux_x_up_0 = flux_x_up(i,kms,j);   flux_x_up_1 = flux_x_up(i+1,kms,j)
                flux_y_up_0 = flux_y_up(i,kms,j);   flux_y_up_1 = flux_y_up(i,kms,j+1)
                flux_z_up_1 = flux_z_up(i,kms+1,j)
                fz = 0.0
                fx  = flux_x(i,kms,j)   - flux_x_up_0
                fx1 = flux_x(i+1,kms,j) - flux_x_up_1
                fy  = flux_y(i,kms,j)   - flux_y_up_0
                fy1 = flux_y(i,kms,j+1) - flux_y_up_1
                fz1 = flux_z(i,kms+1,j) - flux_z_up_1
                temp = q0 - ((flux_x_up_1 - flux_x_up_0) + (flux_y_up_1 - flux_y_up_0) + &
                             flux_z_up_1 * dz_t_i) * denom(i,kms,j)
                flux_in  = ((max(0.0,-fx1) + max(0.0,fx)) + (max(0.0,-fy1) + max(0.0,fy)) + &
                            max(0.0,-fz1) * dz_t_i) * denom(i,kms,j)
                flux_out = ((max(0.0,fx1)  + max(0.0,-fx)) + (max(0.0,fy1)  + max(0.0,-fy)) + &
                            max(0.0,fz1)  * dz_t_i) * denom(i,kms,j)
                scale_in(i,kms,j)  = (qmax_local - temp) / (flux_in  + 1.0e-9)
                scale_out(i,kms,j) = (temp - qmin_local) / (flux_out + 1.0e-9)
            enddo
        enddo

        ! 1b: interior k=kms+1..kme — branchless 3D
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(1000) &
        !$acc private(q0, q_i, q_j, q_k, qmax_local, qmin_local, dz_t_i, temp, &
        !$acc         flux_x_up_0, flux_x_up_1, flux_y_up_0, flux_y_up_1, flux_z_up_0, flux_z_up_1, &
        !$acc         fx, fx1, fy, fy1, fz, fz1, flux_in, flux_out)
        do j = jts-1, jte+1
            do k = kms+1, kme
                do i = its-1, ite+1
                    q0 = q(i,k,j)
                    q_i = q(i+usign(i,k,j),k,j)
                    q_j = q(i,k,j+vsign(i,k,j))
                    q_k = q(i,k+wsign(i,k,j),j)
                    qmax_local = max(q0, q_i, q_j, q_k)
                    qmin_local = min(q0, q_i, q_j, q_k)
                    if ((abs(qmax_local) + abs(qmin_local)) == 0.0) then
                        scale_in(i,k,j) = 0.0
                        scale_out(i,k,j) = 0.0
                        cycle
                    endif
                    dz_t_i = 1.0 / dz(i,k,j)
                    flux_x_up_0 = flux_x_up(i,k,j);   flux_x_up_1 = flux_x_up(i+1,k,j)
                    flux_y_up_0 = flux_y_up(i,k,j);   flux_y_up_1 = flux_y_up(i,k,j+1)
                    flux_z_up_0 = flux_z_up(i,k,j);   flux_z_up_1 = flux_z_up(i,k+1,j)
                    fz  = flux_z(i,k,j)   - flux_z_up_0
                    fx  = flux_x(i,k,j)   - flux_x_up_0
                    fx1 = flux_x(i+1,k,j) - flux_x_up_1
                    fy  = flux_y(i,k,j)   - flux_y_up_0
                    fy1 = flux_y(i,k,j+1) - flux_y_up_1
                    fz1 = flux_z(i,k+1,j) - flux_z_up_1
                    temp = q0 - ((flux_x_up_1 - flux_x_up_0) + (flux_y_up_1 - flux_y_up_0) + &
                                 (flux_z_up_1 - flux_z_up_0) * dz_t_i) * denom(i,k,j)
                    flux_in  = ((max(0.0,-fx1) + max(0.0,fx)) + (max(0.0,-fy1) + max(0.0,fy)) + &
                                (max(0.0,-fz1) + max(0.0,fz)) * dz_t_i) * denom(i,k,j)
                    flux_out = ((max(0.0,fx1)  + max(0.0,-fx)) + (max(0.0,fy1)  + max(0.0,-fy)) + &
                                (max(0.0,fz1)  + max(0.0,-fz)) * dz_t_i) * denom(i,k,j)
                    scale_in(i,k,j)  = (qmax_local - temp) / (flux_in  + 1.0e-9)
                    scale_out(i,k,j) = (temp - qmin_local) / (flux_out + 1.0e-9)
                enddo
            enddo
        enddo

        ! ==========================================
        ! STEP 2: Apply flux corrections. Split off bottom k=kms slab (no z-direction work).
        ! ==========================================

        ! 2a: bottom k=kms — only x and y corrections (z skipped per original)
        !$acc parallel loop gang vector tile(32, 8) async(1001) wait(1000) &
        !$acc private(scale_in_val, scale_out_val, flux_x_up_0, flux_y_up_0, scale)
        do j = jts, jte+1
            do i = its, ite+1
                scale_in_val  = scale_in(i,kms,j)
                scale_out_val = scale_out(i,kms,j)
                flux_x_up_0 = flux_x_up(i,kms,j)
                flux_x(i,kms,j) = flux_x(i,kms,j) - flux_x_up_0
                scale = merge(max(0.0, min(scale_in_val,  scale_out(i-1,kms,j), 1.0)), &
                              max(0.0, min(scale_out_val, scale_in(i-1,kms,j),  1.0)), &
                              flux_x(i,kms,j) > 0.0)
                flux_x(i,kms,j) = scale * flux_x(i,kms,j) + flux_x_up_0

                flux_y_up_0 = flux_y_up(i,kms,j)
                flux_y(i,kms,j) = flux_y(i,kms,j) - flux_y_up_0
                scale = merge(max(0.0, min(scale_in_val,  scale_out(i,kms,j-1), 1.0)), &
                              max(0.0, min(scale_out_val, scale_in(i,kms,j-1),  1.0)), &
                              flux_y(i,kms,j) > 0.0)
                flux_y(i,kms,j) = scale * flux_y(i,kms,j) + flux_y_up_0
            enddo
        enddo

        ! 2b: interior k=kms+1..kme — branchless 3D, includes z-direction correction
        !$acc parallel loop gang vector tile(32, 8, 1) vector_length(256) async(1001) wait(1000) &
        !$acc private(scale_in_val, scale_out_val, flux_x_up_0, flux_y_up_0, flux_z_up_0, scale)
        do j = jts, jte+1
            do k = kms+1, kme
                do i = its, ite+1
                    scale_in_val  = scale_in(i,k,j)
                    scale_out_val = scale_out(i,k,j)

                    flux_x_up_0 = flux_x_up(i,k,j)
                    flux_x(i,k,j) = flux_x(i,k,j) - flux_x_up_0
                    scale = merge(max(0.0, min(scale_in_val,  scale_out(i-1,k,j), 1.0)), &
                                  max(0.0, min(scale_out_val, scale_in(i-1,k,j),  1.0)), &
                                  flux_x(i,k,j) > 0.0)
                    flux_x(i,k,j) = scale * flux_x(i,k,j) + flux_x_up_0

                    flux_y_up_0 = flux_y_up(i,k,j)
                    flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up_0
                    scale = merge(max(0.0, min(scale_in_val,  scale_out(i,k,j-1), 1.0)), &
                                  max(0.0, min(scale_out_val, scale_in(i,k,j-1),  1.0)), &
                                  flux_y(i,k,j) > 0.0)
                    flux_y(i,k,j) = scale * flux_y(i,k,j) + flux_y_up_0

                    flux_z_up_0 = flux_z_up(i,k,j)
                    flux_z(i,k,j) = flux_z(i,k,j) - flux_z_up_0
                    scale = merge(max(0.0, min(scale_in_val,  scale_out(i,k-1,j), 1.0)), &
                                  max(0.0, min(scale_out_val, scale_in(i,k-1,j),  1.0)), &
                                  flux_z(i,k,j) > 0.0)
                    flux_z(i,k,j) = scale * flux_z(i,k,j) + flux_z_up_0
                enddo
            enddo
        enddo
        
        ! Top boundary
        !$acc parallel loop gang vector collapse(2) async(1002) wait(1000)
        do j = jts, jte+1
            do i = its, ite+1
                flux_z_up_1 = flux_z_up(i,kme+1,j)
                flux_z(i,kme+1,j) = flux_z(i,kme+1,j) - flux_z_up_1
                scale = merge(max(0.0, min(scale_in(i,kme,j), scale_out(i,kme,j), 1.0)), &
                              1.0, flux_z(i,kme+1,j) > 0.0)
                flux_z(i,kme+1,j) = scale * flux_z(i,kme+1,j) + flux_z_up_1
            enddo
        enddo
        
        !$acc end data
        !$acc wait(1001, 1002)
        
    end subroutine WRF_flux_corr



    !>------------------------------------------------------------
    !! Fine-mesh flux correction routines
    !! All use explicit bounds (no module-level domain state)
    !!------------------------------------------------------------

    subroutine init_fluxcorr_fm(ks, ke)
        implicit none
        integer, intent(in) :: ks, ke

        if (allocated(usign_fm))      then
            !$acc exit data delete(usign_fm, vsign_fm, scale_in_fm, scale_out_fm, flux_x_up_fm, flux_y_up_fm, dumb_q_fm)
            deallocate(usign_fm)
        endif
        if (allocated(vsign_fm))      deallocate(vsign_fm)
        if (allocated(scale_in_fm))   deallocate(scale_in_fm)
        if (allocated(scale_out_fm))  deallocate(scale_out_fm)
        if (allocated(flux_x_up_fm))  deallocate(flux_x_up_fm)
        if (allocated(flux_y_up_fm))  deallocate(flux_y_up_fm)
        if (allocated(dumb_q_fm))     deallocate(dumb_q_fm)

        allocate(usign_fm(its-1:ite+1, ks:ke, jts-1:jte+1))
        allocate(vsign_fm(its-1:ite+1, ks:ke, jts-1:jte+1))

        allocate(scale_in_fm(its-1:ite+1, ks:ke, jts-1:jte+1))
        allocate(scale_out_fm(its-1:ite+1, ks:ke, jts-1:jte+1))
        allocate(flux_x_up_fm(ims:ime, ks:ke, jms:jme))
        allocate(flux_y_up_fm(ims:ime, ks:ke, jms:jme))
        allocate(dumb_q_fm(ims:ime, ks:ke, jms:jme))

        !$acc enter data create(usign_fm, vsign_fm, scale_in_fm, scale_out_fm, flux_x_up_fm, flux_y_up_fm, dumb_q_fm)
    end subroutine init_fluxcorr_fm


    subroutine set_sign_arrays_fm(u, v, ks, ke)
        implicit none
        integer, intent(in) :: ks, ke
        real, dimension(its-2:ite+3, ks:ke, jts-2:jte+3), intent(in) :: u, v

        integer :: i, j, k

        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts-1, jte+1
            do k = ks, ke
                do i = its-1, ite+1
                    if (u(i,k,j) > 0) then
                        usign_fm(i,k,j) = -1
                    elseif (u(i+1,k,j) < 0) then
                        usign_fm(i,k,j) = 1
                    else
                        usign_fm(i,k,j) = 0
                    end if

                    if (v(i,k,j) > 0) then
                        vsign_fm(i,k,j) = -1
                    elseif (v(i,k,j+1) < 0) then
                        vsign_fm(i,k,j) = 1
                    else
                        vsign_fm(i,k,j) = 0
                    end if
                end do
            end do
        end do
    end subroutine set_sign_arrays_fm


    subroutine compute_upwind_fluxes_fm(q, u, v, denom, ks, ke)
        implicit none
        integer, intent(in) :: ks, ke
        real, dimension(ims:ime, ks:ke, jms:jme), intent(in) :: q, denom
        real, dimension(its-2:ite+3, ks:ke, jts-2:jte+3), intent(in) :: u, v

        real :: q0, u_val, v_val, w_val, abs_u, abs_v, abs_w
        real :: flux_diff_x, flux_diff_y, flux_diff_z, denom_val, dz_val
        integer :: i, j, k

        ! STEP 1: First half-step upwind fluxes
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts-2, jte+3
            do k = ks, ke
                do i = its-2, ite+3
                    q0 = q(i,k,j)

                    u_val = u(i,k,j)
                    abs_u = abs(u_val)
                    flux_x_up_fm(i,k,j) = 0.25 * (u_val * (q(i-1,k,j) + q0) + &
                                                    abs_u * (q(i-1,k,j) - q0))

                    v_val = v(i,k,j)
                    abs_v = abs(v_val)
                    flux_y_up_fm(i,k,j) = 0.25 * (v_val * (q(i,k,j-1) + q0) + &
                                                    abs_v * (q(i,k,j-1) - q0))

                enddo
            enddo
        enddo

        ! return ! at the moment, the fine-mesh flux limiter is only stable when computing the upwind flux
        ! ! off of one large upwind step, instead of decomposing it into two steps to keep it under the CFL for upwind
        ! ! This seems to stem from flux_z being calculated differently. 

        ! STEP 2: Compute intermediate concentration after first half-step
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts-2, jte+2
            do k = ks, ke
                do i = its-2, ite+2
                    flux_diff_x = flux_x_up_fm(i+1,k,j) - flux_x_up_fm(i,k,j)
                    flux_diff_y = flux_y_up_fm(i,k,j+1) - flux_y_up_fm(i,k,j)


                    denom_val = denom(i,k,j)

                    dumb_q_fm(i,k,j) = q(i,k,j) - (flux_diff_x + flux_diff_y) * denom_val
                enddo
            enddo
        enddo

        ! STEP 3: Second half-step upwind fluxes (accumulate)
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts-1, jte+2
            do k = ks, ke
                do i = its-1, ite+2
                    q0 = dumb_q_fm(i,k,j)

                    u_val = u(i,k,j)
                    abs_u = abs(u_val)
                    flux_x_up_fm(i,k,j) = flux_x_up_fm(i,k,j) + 0.25 * (u_val * (dumb_q_fm(i-1,k,j) + q0) + &
                                                                            abs_u * (dumb_q_fm(i-1,k,j) - q0))

                    v_val = v(i,k,j)
                    abs_v = abs(v_val)
                    flux_y_up_fm(i,k,j) = flux_y_up_fm(i,k,j) + 0.25 * (v_val * (dumb_q_fm(i,k,j-1) + q0) + &
                                                                            abs_v * (dumb_q_fm(i,k,j-1) - q0))

                enddo
            enddo
        enddo
    end subroutine compute_upwind_fluxes_fm


    subroutine WRF_flux_corr_fm(q, u, v, flux_x, flux_y, denom, ks, ke)
        implicit none
        integer, intent(in) :: ks, ke
        real, dimension(ims:ime, ks:ke, jms:jme), intent(in) :: q, denom
        real, dimension(its-2:ite+3, ks:ke, jts-2:jte+3), intent(in) :: u, v
        real, dimension(its-1:ite+2, ks:ke, jts-1:jte+2), intent(inout) :: flux_x, flux_y

        real :: dz_t_i, q0, q_i, q_j, temp, flux_in, flux_out
        real :: scale, scale_in_val, scale_out_val, qmax_local, qmin_local
        real :: fx, fx1, fy, fy1, fz, fz1
        real :: flux_x_up_0, flux_x_up_1, flux_y_up_0, flux_y_up_1
        integer :: i, j, k

        ! Compute two-step upwind fluxes
        call compute_upwind_fluxes_fm(q, u, v, denom, ks, ke)

        ! STEP 1: Compute scale_in/scale_out using upwind fluxes
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts-1, jte+1
            do k = ks, ke
                do i = its-1, ite+1
                    q0 = q(i,k,j)
                    q_i = q(i+usign_fm(i,k,j),k,j)
                    q_j = q(i,k,j+vsign_fm(i,k,j))

                    qmax_local = max(q0, q_i, q_j)
                    qmin_local = min(q0, q_i, q_j)

                    if ((abs(qmax_local) + abs(qmin_local)) == 0.0) then
                        scale_in_fm(i,k,j) = 0.0
                        scale_out_fm(i,k,j) = 0.0
                        cycle
                    endif

                    flux_x_up_0 = flux_x_up_fm(i,k,j)
                    flux_x_up_1 = flux_x_up_fm(i+1,k,j)
                    flux_y_up_0 = flux_y_up_fm(i,k,j)
                    flux_y_up_1 = flux_y_up_fm(i,k,j+1)

                    fx = flux_x(i,k,j) - flux_x_up_0
                    fx1 = flux_x(i+1,k,j) - flux_x_up_1
                    fy = flux_y(i,k,j) - flux_y_up_0
                    fy1 = flux_y(i,k,j+1) - flux_y_up_1

                    temp = q0 - ((flux_x_up_1 - flux_x_up_0) + &
                                (flux_y_up_1 - flux_y_up_0) ) * denom(i,k,j)

                    flux_in = ((max(0.0,-fx1) + max(0.0,fx)) + &
                               (max(0.0,-fy1) + max(0.0,fy)) ) * denom(i,k,j)

                    flux_out = ((max(0.0,fx1) + max(0.0,-fx)) + &
                                (max(0.0,fy1) + max(0.0,-fy)) ) * denom(i,k,j)

                    scale_in_fm(i,k,j) = (qmax_local - temp) / (flux_in + 1.0e-9)
                    scale_out_fm(i,k,j) = (temp - qmin_local) / (flux_out + 1.0e-9)
                enddo
            enddo
        enddo

        ! STEP 2: Apply flux corrections
        !$acc parallel loop gang vector collapse(3) default(present)
        do j = jts, jte+1
            do k = ks, ke
                do i = its, ite+1
                    scale_in_val = scale_in_fm(i,k,j)
                    scale_out_val = scale_out_fm(i,k,j)

                    ! X-direction
                    flux_x_up_0 = flux_x_up_fm(i,k,j)
                    flux_x(i,k,j) = flux_x(i,k,j) - flux_x_up_0
                    scale = merge(max(0.0, min(scale_in_val, scale_out_fm(i-1,k,j), 1.0)), &
                                  max(0.0, min(scale_out_val, scale_in_fm(i-1,k,j), 1.0)), &
                                  flux_x(i,k,j) > 0.0)
                    flux_x(i,k,j) = scale * flux_x(i,k,j) + flux_x_up_0

                    ! Y-direction
                    flux_y_up_0 = flux_y_up_fm(i,k,j)
                    flux_y(i,k,j) = flux_y(i,k,j) - flux_y_up_0
                    scale = merge(max(0.0, min(scale_in_val, scale_out_fm(i,k,j-1), 1.0)), &
                                  max(0.0, min(scale_out_val, scale_in_fm(i,k,j-1), 1.0)), &
                                  flux_y(i,k,j) > 0.0)
                    flux_y(i,k,j) = scale * flux_y(i,k,j) + flux_y_up_0

                enddo
            enddo
        enddo

    end subroutine WRF_flux_corr_fm


end module adv_fluxcorr
