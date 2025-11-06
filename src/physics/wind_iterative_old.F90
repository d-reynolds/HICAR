!>------------------------------------------------------------
!! Module to solve for a 3D wind field following mass-conservation
!! and reducing differennces between initial and final wind field.
!! Solver requires use of PETSc software package, and so is 
!! separated from the rest of wind core here to allow for dependancy
!! on PETSc library.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
!#include <petsc/finclude/petscdmda.h90>

module wind_iterative_old    
#include "petscversion.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscdmda.h"

    use domain_interface,  only : domain_t
    !use options_interface, only : options_t
    !use grid_interface,    only : grid_t
    use petscksp
    use petscdm
    use petscdmda
    use icar_constants,    only : STD_OUT_PE, kVARS
    use options_interface, only : options_t
    use iso_fortran_env

    implicit none
    private
    public:: init_iter_winds_old, calc_iter_winds_old, finalize_iter_winds_old, finalize_petsc
    real, parameter::deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    logical :: initialized = .False.

    real, allocatable, dimension(:,:,:) :: A_coef, B_coef, C_coef, D_coef, E_coef, F_coef, G_coef, H_coef, I_coef, &
                                           J_coef, K_coef, L_coef, M_coef, N_coef, O_coef
    real    :: dx
    real, allocatable, dimension(:,:,:)  :: div, dz_if, jaco, dzdx, dzdy, sigma, alpha
    real, allocatable, dimension(:,:)    :: dzdx_surf, dzdy_surf
    integer, allocatable :: xl(:), yl(:)
    integer              :: hs, i_s, i_e, k_s, k_e, j_s, j_e


    type(tKSP), allocatable, dimension(:) ::    ksp
    type(tDM)             da
    type(tVec)            localX, b
    type(tMat)            arr_A

contains



    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx + dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine calc_iter_winds_old(domain,alpha_in,div_in,adv_den,update_in)
        implicit none
        type(domain_t), intent(inout) :: domain
        real, intent(in), dimension(domain%grid%ims:domain%grid%ime, &
                                    domain%grid%kms:domain%grid%kme, &
                                    domain%grid%jms:domain%grid%jme) :: alpha_in, div_in
        logical, intent(in) :: adv_den
        logical, optional, intent(in) :: update_in

        PetscScalar,pointer :: lambda(:,:,:)
        logical             :: update

        integer :: i, j, k 
        logical :: varying_alpha

        PetscErrorCode ierr

        PetscInt       one, x_size, iteration, maxits
        PetscReal      rtol, abstol, dtol

        KSPConvergedReason reason

        update=.False.
        if (present(update_in)) update=update_in
                
        !Initialize div to be the initial divergence of the input wind field
        !$acc data present(div_in, alpha_in, div, alpha)
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_s, i_e
                    div(i,k,j) = div_in(i,k,j)
                    alpha(i,k,j) = alpha_in(i,k,j)
                end do
            end do
        end do
        !$acc update host(alpha)
        !$acc end data 

        ! set minimum number of iterations for ksp solver
        rtol = 1.0e-5
        abstol = 1.0e-5
        dtol = 1000.0
        maxits = 1000
        call KSPSetTolerances(ksp(domain%nest_indx),PETSC_DETERMINE_REAL,PETSC_DETERMINE_REAL,dtol, maxits,ierr)

        varying_alpha = .not.( ALL(alpha==minval(alpha)) )
        if (.not.(allocated(A_coef))) then
            ! Can't be called in module-level init function, since we first need alpha
            call initialize_coefs(domain)
            call ComputeMatrix(ksp(domain%nest_indx),arr_A,prealloc=.True.)
        elseif (varying_alpha) then
            call update_coefs(domain)
            call ComputeMatrix(ksp(domain%nest_indx),arr_A)
        endif

        call KSPSetOperators(ksp(domain%nest_indx),arr_A,arr_A,ierr)
        ! if (update) then
        !     call KSPSetReusePreconditioner(ksp(domain%nest_indx),PETSC_TRUE,ierr)
        ! endif

        call ComputeRHS(ksp(domain%nest_indx),b,ierr)
        
        call KSPSetDMActive(ksp(domain%nest_indx),PETSC_FALSE,ierr)
        call KSPSolve(ksp(domain%nest_indx),b,b,ierr)

        call KSPGetConvergedReason(ksp(domain%nest_indx), reason, ierr)
        
        if (reason > 0) then
            call KSPGetIterationNumber(ksp(domain%nest_indx), iteration, ierr)
            if(STD_OUT_PE) write(*,*) 'Solved PETSc after ',iteration,' iterations'
            if (STD_OUT_PE) flush(output_unit)
            !Subset global solution to local grid so that we can access ghost-points
            call DMGlobalToLocal(da,b,INSERT_VALUES,localX,ierr)

            call DMDAVecGetArrayReadF90(da,localX,lambda, ierr)
            call calc_updated_winds(domain, lambda, adv_den)
            call DMDAVecRestoreArrayReadF90(da,localX,lambda, ierr)

          else
            if (STD_OUT_PE) write(*,*) 'WARNING: PETSc did not converge, PETSc reason was: ', reason
            if (STD_OUT_PE) flush(output_unit)
        endif
    end subroutine calc_iter_winds_old

    subroutine calc_updated_winds(domain,lambda,adv_den) !u, v, w, jaco_u,jaco_v,jaco_w,u_dzdx,v_dzdy,lambda, ids, ide, jds, jde)
        type(domain_t), intent(inout) :: domain
        PetscScalar, intent(in), pointer       :: lambda(:,:,:)
        logical,     intent(in)                :: adv_den


        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, dlambdz, u_temp, v_temp, lambda_too, rho, rho_u, rho_v
        integer i, j, k, i_start, i_end, j_start, j_end 

        i_start = i_s
        i_end   = i_e+1
        j_start = j_s
        j_end   = j_e+1

        allocate(u_temp(i_start:i_end,k_s-1:k_e+1,j_s:j_e))
        allocate(v_temp(i_s:i_e,k_s-1:k_e+1,j_start:j_end))

        allocate(u_dlambdz(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,k_s:k_e,j_start:j_end))
        allocate(dlambdz(i_s-1:i_e+1,k_s:k_e,j_s-1:j_e+1))

        allocate(rho(domain%ims:domain%ime,k_s:k_e,domain%jms:domain%jme))
        allocate(rho_u(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(rho_v(i_s:i_e,k_s:k_e,j_start:j_end))


        associate(density => domain%vars_3d(domain%var_indx(kVARS%density)%v)%data_3d, &
                    u       => domain%vars_3d(domain%var_indx(kVARS%u)%v)%data_3d, &
                    v       => domain%vars_3d(domain%var_indx(kVARS%v)%v)%data_3d, &
                    w       => domain%vars_3d(domain%var_indx(kVARS%w)%v)%data_3d, &
                    jaco_u_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d, &
                    jaco_v_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d, &
                    jaco_domain   => domain%vars_3d(domain%var_indx(kVARS%jacobian_w)%v)%data_3d, &
                    dzdx_u => domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d, &
                    dzdy_v => domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d)
        ! Create GPU data for local arrays
        !$acc enter data create(u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v)
        !$acc data present(density, u, v, jaco_u_domain, jaco_v_domain, jaco_domain, dzdx_u, dzdy_v, lambda, u_temp, v_temp, u_dlambdz, v_dlambdz, dlambdz, rho, rho_u, rho_v, alpha, dz_if) copyin(lambda)

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
                    u_temp(i,k,j) = (lambda(i,k,j) + lambda(i-1,k,j)) / 2 
                end do
            end do
        end do

        !stager lambda to v grid - GPU computation
        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s-1, k_e+1
                do i = i_s, i_e
                    v_temp(i,k,j) = (lambda(i,k,j) + lambda(i,k,j-1)) / 2 
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
                    dlambdz(i,k,j) = lambda(i,k+1,j) - lambda(i,k-1,j)
                    dlambdz(i,k,j) = dlambdz(i,k,j)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
                end do
            end do
        enddo

        ! Update wind fields on GPU
        !$acc parallel loop gang vector collapse(3)
        do j = j_s, j_e
            do k = k_s, k_e
                do i = i_start, i_end
                    u(i,k,j) = u(i,k,j) + 0.5*((lambda(i,k,j) - lambda(i-1,k,j))/dx - &
                                        dzdx_u(i,k,j)*(u_dlambdz(i,k,j))/jaco_u_domain(i,k,j))/(rho_u(i,k,j))
                end do
            end do
        end do
        
        !$acc parallel loop gang vector collapse(3)
        do j = j_start, j_end
            do k = k_s, k_e
                do i = i_s, i_e
                    v(i,k,j) = v(i,k,j) + 0.5*((lambda(i,k,j) - lambda(i,k,j-1))/dx - dzdy_v(i,k,j)*(v_dlambdz(i,k,j))/jaco_v_domain(i,k,j))/(rho_v(i,k,j))
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

    subroutine ComputeRHS(ksp,vec_b,ierr)
        implicit none
        
        PetscErrorCode  ierr
        type(tKSP) ksp
        type(tVec) vec_b, x
        type(tDM)             dm

        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar,pointer :: barray(:,:,:)
        PetscBool :: has_cuda

        call KSPGetDM(ksp,dm,ierr)
        
        call DMDAGetLocalInfoF90(dm,info,ierr)
        !swap z and y in info calls since dims 2/3 are transposed for PETSc w.r.t. HICAR indexing
        mx = info(DMDA_LOCAL_INFO_MX)
        my = info(DMDA_LOCAL_INFO_MZ)
        mz = info(DMDA_LOCAL_INFO_MY)
        xm = info(DMDA_LOCAL_INFO_XM)
        ym = info(DMDA_LOCAL_INFO_ZM)
        zm = info(DMDA_LOCAL_INFO_YM)
        xs = info(DMDA_LOCAL_INFO_XS)
        ys = info(DMDA_LOCAL_INFO_ZS)
        zs = info(DMDA_LOCAL_INFO_YS)

        call DMDAVecGetArrayF90(dm,vec_b,barray, ierr)
        
        !$acc parallel present(div) copy(barray)
        !$acc loop gang vector collapse(3)
        do j=ys,(ys+ym-1)
            do k=zs,(zs+zm-1)
                do i=xs,(xs+xm-1)
                    !For global boundary conditions
                    if (i.le.0 .or. j.le.0 .or. k.eq.0 .or. &
                        i.ge.mx-1 .or. j.ge.my-1 .or. k.eq.mz-1) then
                        barray(i,k,j) = 0.0
                    else
                        barray(i,k,j) = -2*div(i,k,j)
                    endif
                enddo
            enddo
        enddo
        !$acc end parallel

        call DMDAVecRestoreArrayF90(dm,vec_b,barray, ierr)
    end subroutine ComputeRHS

    subroutine ComputeMatrix(ksp,array_A,prealloc)
        implicit none
        type(tKSP) ksp
        type(tMat) array_A!,array_B
        logical, intent(in), optional :: prealloc

        PetscErrorCode  ierr
        type(tDM)             da
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i2,i6,i15
        DMDALocalInfo  info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar, allocatable, dimension(:) ::    v
        PetscInt, allocatable, dimension(:) ::       row_indx, col_indx
        MatStencil     row(4),col(4,15),gnd_col(4,10),top_col(4,2)

        real :: denom
        integer :: cnt, num_entries
        logical :: prealloc_flag

        i1 = 1
        i2 = 2
        i6 = 10
        i15 = 15

        prealloc_flag = .False.
        if (present(prealloc)) prealloc_flag = prealloc

        call KSPGetDM(ksp,da,ierr)

        call DMDAGetLocalInfoF90(da,info,ierr)
        !swap z and y in info calls since dims 2/3 are transposed for PETSc w.r.t. HICAR indexing
        mx = info(DMDA_LOCAL_INFO_MX)
        my = info(DMDA_LOCAL_INFO_MZ)
        mz = info(DMDA_LOCAL_INFO_MY)
        xm = info(DMDA_LOCAL_INFO_XM)
        ym = info(DMDA_LOCAL_INFO_ZM)
        zm = info(DMDA_LOCAL_INFO_YM)
        xs = info(DMDA_LOCAL_INFO_XS)
        ys = info(DMDA_LOCAL_INFO_ZS)
        zs = info(DMDA_LOCAL_INFO_YS)
        
        num_entries = (xm)*(ym)*(zm)*15
        allocate(row_indx((num_entries)), col_indx((num_entries)), v(num_entries))

        row_indx = -1
        col_indx = -1
        cnt = 1
        do j=ys,(ys+ym-1)
        do k=zs,(zs+zm-1)
            do i=xs,(xs+xm-1)
                row(MatStencil_i) = i
                row(MatStencil_j) = k
                row(MatStencil_k) = j
                call DMDAMapMatStencilToGlobal(da,1,row,row_indx(cnt),ierr)

                if (i.le.0 .or. j.le.0 .or. &
                    i.ge.mx-1 .or. j.ge.my-1) then
                    v(cnt) = 1.0
                    call DMDAMapMatStencilToGlobal(da, 1, row, col_indx(cnt), ierr)
                    cnt = cnt + 1

                else if (k.eq.mz-1) then
                    !k
                    v(cnt) = 1/dz_if(i,k,j)
                    top_col(MatStencil_i,1) = i
                    top_col(MatStencil_j,1) = k
                    top_col(MatStencil_k,1) = j
                    !k - 1
                    v(cnt+1) = -1/dz_if(i,k,j)
                    top_col(MatStencil_i,2) = i
                    top_col(MatStencil_j,2) = k-1
                    top_col(MatStencil_k,2) = j
                    call DMDAMapMatStencilToGlobal(da, 2, top_col, col_indx(cnt:cnt+1), ierr)
                    row_indx(cnt:cnt+1) = row_indx(cnt)
                    cnt = cnt + 2
                else if (k.eq.0) then
                    denom = 2*(dzdx_surf(i,j)**2 + dzdy_surf(i,j)**2 + alpha(i,1,j)**2)/(jaco(i,1,j))
                    !k
                    v(cnt+1) = -1/dz_if(i,k+1,j)
                    gnd_col(MatStencil_i,2) = i
                    gnd_col(MatStencil_j,2) = k
                    gnd_col(MatStencil_k,2) = j
                    !k + 1
                    v(cnt) = 1/dz_if(i,k+1,j)
                    gnd_col(MatStencil_i,1) = i
                    gnd_col(MatStencil_j,1) = k+1
                    gnd_col(MatStencil_k,1) = j

                    !Reminder: The equation for the lower BC has all of the lateral derivatives (dlambdadx, etc) NEGATIVE, hence opposite signs below
                    !If we are on left most border
                    !i - 1, k + 1
                    v(cnt+2) = dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,3) = i-1
                    gnd_col(MatStencil_j,3) = k
                    gnd_col(MatStencil_k,3) = j
                    !i + 1, k + 1
                    v(cnt+3) = -dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,4) = i+1
                    gnd_col(MatStencil_j,4) = k
                    gnd_col(MatStencil_k,4) = j
                    !i - 1, k + 1
                    v(cnt+4) = dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,5) = i-1
                    gnd_col(MatStencil_j,5) = k+1
                    gnd_col(MatStencil_k,5) = j
                    !i + 1, k + 1
                    v(cnt+5) = -dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,6) = i+1
                    gnd_col(MatStencil_j,6) = k+1
                    gnd_col(MatStencil_k,6) = j
                    !j - 1, k + 1
                    v(cnt+6) = dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,7) = i
                    gnd_col(MatStencil_j,7) = k
                    gnd_col(MatStencil_k,7) = j-1
                    !j + 1, k + 1
                    v(cnt+7) = -dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,8) = i
                    gnd_col(MatStencil_j,8) = k
                    gnd_col(MatStencil_k,8) = j+1
                    !j - 1, k + 1
                    v(cnt+8) = dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,9) = i
                    gnd_col(MatStencil_j,9) = k+1
                    gnd_col(MatStencil_k,9) = j-1
                    !j + 1, k + 1
                    v(cnt+9) = -dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,10) = i
                    gnd_col(MatStencil_j,10) = k+1
                    gnd_col(MatStencil_k,10) = j+1
                    call DMDAMapMatStencilToGlobal(da, 10, gnd_col, col_indx(cnt:cnt+9), ierr)
                    row_indx(cnt:cnt+9) = row_indx(cnt)
                    cnt = cnt + 10
                else
                    !j - 1, k - 1
                    v(cnt) = O_coef(i,k,j)
                    col(MatStencil_i,1) = i
                    col(MatStencil_j,1) = k-1
                    col(MatStencil_k,1) = j-1
                    !i - 1, k - 1
                    v(cnt+1) = K_coef(i,k,j)
                    col(MatStencil_i,2) = i-1
                    col(MatStencil_j,2) = k-1
                    col(MatStencil_k,2) = j
                    !j - 1, k + 1
                    v(cnt+2) = M_coef(i,k,j)
                    col(MatStencil_i,3) = i
                    col(MatStencil_j,3) = k+1
                    col(MatStencil_k,3) = j-1
                    !i - 1, k + 1
                    v(cnt+3) = I_coef(i,k,j)
                    col(MatStencil_i,4) = i-1
                    col(MatStencil_j,4) = k+1
                    col(MatStencil_k,4) = j
                    !k - 1
                    v(cnt+4) = C_coef(i,k,j)
                    col(MatStencil_i,5) = i
                    col(MatStencil_j,5) = k-1
                    col(MatStencil_k,5) = j
                    !j - 1
                    v(cnt+5) = G_coef(i,k,j)
                    col(MatStencil_i,6) = i
                    col(MatStencil_j,6) = k
                    col(MatStencil_k,6) = j-1
                    !i - 1
                    v(cnt+6) = E_coef(i,k,j)
                    col(MatStencil_i,7) = i-1
                    col(MatStencil_j,7) = k
                    col(MatStencil_k,7) = j
                    !Center
                    v(cnt+7) = A_coef(i,k,j)
                    col(MatStencil_i,8) = i
                    col(MatStencil_j,8) = k
                    col(MatStencil_k,8) = j
                    !i + 1
                    v(cnt+8) = D_coef(i,k,j)
                    col(MatStencil_i,9) = i+1
                    col(MatStencil_j,9) = k
                    col(MatStencil_k,9) = j
                    !j + 1
                    v(cnt+9) = F_coef(i,k,j)
                    col(MatStencil_i,10) = i
                    col(MatStencil_j,10) = k
                    col(MatStencil_k,10) = j+1
                    !k + 1
                    v(cnt+10) = B_coef(i,k,j)
                    col(MatStencil_i,11) = i
                    col(MatStencil_j,11) = k+1
                    col(MatStencil_k,11) = j
                    !i + 1, k + 1
                    v(cnt+11) = H_coef(i,k,j)
                    col(MatStencil_i,12) = i+1
                    col(MatStencil_j,12) = k+1
                    col(MatStencil_k,12) = j
                    !j + 1, k + 1
                    v(cnt+12) = L_coef(i,k,j)
                    col(MatStencil_i,13) = i
                    col(MatStencil_j,13) = k+1
                    col(MatStencil_k,13) = j+1
                    !i + 1, k - 1
                    v(cnt+13) = J_coef(i,k,j)
                    col(MatStencil_i,14) = i+1
                    col(MatStencil_j,14) = k-1
                    col(MatStencil_k,14) = j
                    !j + 1, k - 1
                    v(cnt+14) = N_coef(i,k,j)
                    col(MatStencil_i,15) = i
                    col(MatStencil_j,15) = k-1
                    col(MatStencil_k,15) = j+1
                    call DMDAMapMatStencilToGlobal(da, 15, col, col_indx(cnt:cnt+14), ierr)
                    row_indx(cnt:cnt+14) = row_indx(cnt)
                    cnt = cnt + 15

                endif
            enddo
        enddo
        enddo

        if (prealloc_flag) call MatSetPreallocationCOO(array_A,num_entries,row_indx,col_indx, ierr)
        call MatSetValuesCOO(array_A,v,INSERT_VALUES, ierr)


    end subroutine ComputeMatrix

    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        real :: mixed_denom
        integer :: k

        allocate(A_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(B_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(C_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(D_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(E_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(F_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(G_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(H_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(I_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(J_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(K_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(L_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(M_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(N_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(O_coef(i_s:i_e,k_s:k_e,j_s:j_e))
                
        A_coef = 1 !Because this corresponds to the centered node, must be non-zero, otherwise there is no solution
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


        !This function sets A, B, annd C coefficients
        call update_coefs(domain)
                
        !First terms -- d2lam/dx2
                
        !i+1
        D_coef = 1/(domain%dx**2)
        !i-1
        E_coef = 1/(domain%dx**2)
        !j+1
        F_coef = 1/(domain%dx**2)
        !j-1
        G_coef = 1/(domain%dx**2)
        
        !Now consider mixed derivatives -- d2lam/dxdz
        do k=k_s,k_e
            !sigma and dz_if are both horizontally constant, so can be set to constants for each layer
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

        ! --------------------------------------------------------
        ! Alternative method of calculating coefficients, save for posterity
        ! --------------------------------------------------------
        ! mixed_denom = 2*domain%dx*dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)*jaco

        ! H_coef = -(sigma**2)*2*(dzdx)/mixed_denom
        ! I_coef = (sigma**2)*2*(dzdx)/mixed_denom
        ! J_coef = 2*(dzdx)/mixed_denom
        ! K_coef = -2*(dzdx)/mixed_denom

        ! L_coef = -(sigma**2)*2*(dzdy)/mixed_denom
        ! M_coef = (sigma**2)*2*(dzdy)/mixed_denom
        ! N_coef = 2*(dzdy)/mixed_denom
        ! O_coef = -2*(dzdy)/mixed_denom
        ! --------------------------------------------------------
        ! End alternative method of calculating coefficients
        ! --------------------------------------------------------

    end subroutine initialize_coefs
    
    
    !Update the coefs which change with time, i.e. those which depend on alpha
    subroutine update_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain        
        real, allocatable, dimension(:,:,:) :: mixed_denom, X_coef
        real, allocatable, dimension(:,:)   :: M_up, M_dwn
        integer :: k
        real :: denom_tmp, sig_i
        
        allocate(mixed_denom(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(X_coef(i_s:i_e,k_s:k_e,j_s:j_e))
        allocate(M_up(i_s:i_e,j_s:j_e))
        allocate(M_dwn(i_s:i_e,j_s:j_e))
        
        ! --------------------------------------------------------
        ! Alternative method of calculating coefficients, save for posterity
        ! --------------------------------------------------------
        ! allocate(dijacodx(i_s:i_e,k_s:k_e,j_s:j_e))
        ! allocate(dijacody(i_s:i_e,k_s:k_e,j_s:j_e))
        ! allocate(dzdxz(i_s:i_e,k_s:k_e,j_s:j_e))
        ! allocate(dzdyz(i_s:i_e,k_s:k_e,j_s:j_e))
        ! allocate(d2zdx(i_s:i_e,k_s:k_e,j_s:j_e))
        ! allocate(d2zdy(i_s:i_e,k_s:k_e,j_s:j_e))

        ! i_s_bnd = i_s
        ! i_e_bnd = i_e
        ! j_s_bnd = j_s
        ! j_e_bnd = j_e
        ! if (i_s==domain%grid%ids) then
        !     i_s_bnd = i_s + 1
        ! endif
        ! if (i_e==domain%grid%ide) then
        !     i_e_bnd = i_e - 1
        ! endif
        ! if (j_s==domain%grid%jds) then
        !     j_s_bnd = j_s + 1
        ! endif
        ! if (j_e==domain%grid%jde) then
        !     j_e_bnd = j_e - 1
        ! endif
        ! mixed_denom = dz_if(:,k_s+1:k_e+1,:)*(sigma+sigma**2)

        ! ! dijaco -- Derivative of Inverse JACObian
        ! dijacodx = (1/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e) - &
        !           1/domain%vars_3d(domain%var_indx(kVARS%jacobian_u)%v)%data_3d(i_s:i_e,:,j_s:j_e))/dx
        ! dijacody = (1/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1) - &
        !           1/domain%vars_3d(domain%var_indx(kVARS%jacobian_v)%v)%data_3d(i_s:i_e,:,j_s:j_e))/dx

        ! dijacodx(i_s_bnd:i_e_bnd,:,:) = (1/domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d(i_s_bnd+1:i_e_bnd+1,:,j_s:j_e) - &
        !           1/domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d(i_s_bnd-1:i_e_bnd-1,:,j_s:j_e))/(2*dx)
        ! dijacody(:,:,j_s_bnd:j_e_bnd) = (1/domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d(i_s:i_e,:,j_s_bnd+1:j_e_bnd+1) - &
        !           1/domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d(i_s:i_e,:,j_s_bnd-1:j_e_bnd-1))/(2*dx)

        ! d2zdx = (domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s+1:i_e+1,:,j_s:j_e) - &
        !           domain%vars_3d(domain%var_indx(kVARS%dzdx_u)%v)%data_3d(i_s:i_e,:,j_s:j_e))/dx
        ! d2zdy = (domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s+1:j_e+1) - &
        !           domain%vars_3d(domain%var_indx(kVARS%dzdy_v)%v)%data_3d(i_s:i_e,:,j_s:j_e))/dx

        ! d2zdx(i_s_bnd:i_e_bnd,:,:) = (domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s_bnd+1:i_e_bnd+1,:,j_s:j_e) - &
        !        2*domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s_bnd:i_e_bnd,:,j_s:j_e)     + &
        !          domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s_bnd-1:i_e_bnd-1,:,j_s:j_e))/(dx**2)

        ! d2zdy(:,:,j_s_bnd:j_e_bnd) = (domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s:i_e,:,j_s_bnd+1:j_e_bnd+1) - &
        !        2*domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s:i_e,:,j_s_bnd:j_e_bnd)     + &
        !          domain%vars_3d(domain%var_indx(kVARS%z)%v)%data_3d(i_s:i_e,:,j_s_bnd-1:j_e_bnd-1))/(dx**2)

        ! dzdxz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2)*dzdx(:,k_s+2:k_e,:) - (sigma(:,k_s+1:k_e-1,:)**2-1)*dzdx(:,k_s+1:k_e-1,:) - dzdx(:,k_s:k_e-2,:)
        ! dzdxz = dzdxz/mixed_denom
        ! dzdxz(:,k_s,:) = -(sigma(:,k_s+1,:)**2)*dzdx(:,k_s+2,:)/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
        !                   (sigma(:,k_s+1,:)+1)*dzdx(:,k_s+1,:)/dz_if(:,k_s+1,:) - &
        !                   (2*sigma(:,k_s+1,:)+1)*dzdx(:,k_s,:)/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))

        ! dzdxz(:,k_e,:) = (dzdx(:,k_e-2,:))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)+dz_if(:,k_e,:)*sigma(:,k_e-1,:)**2) - &
        !                    (sigma(:,k_e-1,:)+1)*(dzdx(:,k_e-1,:))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
        !                     (2+sigma(:,k_e-1,:))*(dzdx(:,k_e,:))/(dz_if(:,k_e,:)*(sigma(:,k_e-1,:)+1))


        ! dzdyz(:,k_s+1:k_e-1,:) = (sigma(:,k_s+1:k_e-1,:)**2)*dzdy(:,k_s+2:k_e,:) - (sigma(:,k_s+1:k_e-1,:)**2-1)*dzdy(:,k_s+1:k_e-1,:) - dzdy(:,k_s:k_e-2,:)
        ! dzdyz = dzdyz/mixed_denom
        ! dzdyz(:,k_s,:) = -(sigma(:,k_s+1,:)**2)*dzdy(:,k_s+2,:)/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1)) + &
        !                   (sigma(:,k_s+1,:)+1)*dzdy(:,k_s+1,:)/dz_if(:,k_s+1,:) - &
        !                   (2*sigma(:,k_s+1,:)+1)*dzdy(:,k_s,:)/(dz_if(:,k_s+1,:)*(sigma(:,k_s+1,:)+1))
        ! dzdyz(:,k_e,:) = (dzdy(:,k_e-2,:))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)+dz_if(:,k_e,:)*sigma(:,k_e-1,:)**2) - &
        !                   (sigma(:,k_e-1,:)+1)*(dzdy(:,k_e-1,:))/(dz_if(:,k_e,:)*sigma(:,k_e-1,:)) + &
        !                    (2+sigma(:,k_e-1,:))*(dzdy(:,k_e,:))/(dz_if(:,k_e,:)*(sigma(:,k_e-1,:)+1))

        ! !Numerical solution
        ! X_coef = (-(d2zdx/jaco + dijacodx*dzdx + d2zdy/jaco + dijacody*dzdy) + &
        !             (dzdxz*dzdx+dzdyz*dzdy)/(jaco**2)) /mixed_denom

        ! B_coef = 2*sigma * ( (alpha**2 + dzdy**2 + dzdx**2)) / ((sigma+sigma**2)*dz_if(:,k_s+1:k_e+1,:)**2)
        ! C_coef = 2*( (alpha**2 + dzdy**2 + dzdx**2)) / ((sigma+sigma**2)*dz_if(:,k_s+1:k_e+1,:)**2)
                
        ! B_coef = B_coef / (jaco**2)
        ! C_coef = C_coef / (jaco**2)
        ! --------------------------------------------------------
        ! End alternative method of calculating coefficients
        ! --------------------------------------------------------

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
        
                                                       
    end subroutine
    

    subroutine finalize_iter_winds_old()
        implicit none

        PetscErrorCode ierr

        if (allocated(A_coef)) deallocate(A_coef)
        if (allocated(B_coef)) deallocate(B_coef)
        if (allocated(C_coef)) deallocate(C_coef)
        if (allocated(D_coef)) deallocate(D_coef)
        if (allocated(E_coef)) deallocate(E_coef)
        if (allocated(F_coef)) deallocate(F_coef)
        if (allocated(G_coef)) deallocate(G_coef)
        if (allocated(H_coef)) deallocate(H_coef)
        if (allocated(I_coef)) deallocate(I_coef)
        if (allocated(J_coef)) deallocate(J_coef)
        if (allocated(K_coef)) deallocate(K_coef)
        if (allocated(L_coef)) deallocate(L_coef)
        if (allocated(M_coef)) deallocate(M_coef)
        if (allocated(N_coef)) deallocate(N_coef)
        if (allocated(O_coef)) deallocate(O_coef)
        
        ! Clean up GPU data before deallocating
        if (allocated(div)) then
            !$acc exit data delete(dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma, dz_if, alpha, div)
            deallocate(div)
        endif
        if (allocated(dz_if)) deallocate(dz_if)
        if (allocated(jaco)) deallocate(jaco)
        if (allocated(dzdx)) deallocate(dzdx)
        if (allocated(dzdy)) deallocate(dzdy)
        if (allocated(dzdx_surf)) deallocate(dzdx_surf)
        if (allocated(dzdy_surf)) deallocate(dzdy_surf)
        if (allocated(sigma)) deallocate(sigma)
        if (allocated(alpha)) deallocate(alpha)
        if (allocated(xl)) deallocate(xl)

        ! because if yl is allocated, then the rest of the da variables would be allocated
        if (allocated(yl)) then
             deallocate(yl)
        
            call VecDestroy(localX,ierr)
            call VecDestroy(b, ierr)
            call MatDestroy(arr_A,ierr)
            call DMDestroy(da,ierr)
        endif

    end subroutine finalize_iter_winds_old

    subroutine init_petsc_comms(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options
        PetscErrorCode ierr
        type(tPC) precond
        integer :: i

        PETSC_COMM_WORLD = domain%compute_comms%MPI_VAL
        call PetscInitialize(ierr)

        allocate(ksp(options%general%nests))

        do i=1,options%general%nests
            call KSPCreate(domain%compute_comms%MPI_VAL,ksp(i),ierr)
            call KSPSetType(ksp(i),KSPPIPEGCR,ierr) !KSPPIPEGCR <-- this one tested to give fastest convergence...
                                            !KSPFBCGS <-- this one used previously...
            call KSPSetFromOptions(ksp(i),ierr)
            call KSPSetReusePreconditioner(ksp(domain%nest_indx),PETSC_FALSE,ierr)
            call KSPGetPC(ksp(i),precond,ierr)
            call PCFactorSetUseInPlace(precond,PETSC_TRUE,ierr)
            call PCFactorSetReuseOrdering(precond,PETSC_TRUE,ierr)

        enddo
        initialized = .True.

        if (ierr .ne. 0) then
            print*,'Unable to initialize PETSc'
            stop
        endif

    end subroutine init_petsc_comms

    subroutine finalize_petsc()
        implicit none
        PetscErrorCode ierr

        if (initialized) then
            call PetscFinalize(ierr)
        endif

    end subroutine finalize_petsc

    subroutine init_iter_winds_old(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options

        PetscInt       one, iter
        PetscErrorCode ierr
        PetscBool :: has_cuda
        ISLocalToGlobalMapping isltog

        ! call finalize routine to deallocate any arrays that are already allocated. 
        ! This would only occur if another nest was using this module previously. 
        call finalize_iter_winds_old()

        if (.not.(initialized)) then
            call init_petsc_comms(domain, options)
        endif


        call init_module_vars(domain)

        one = 1

        call DMDACreate3d(domain%compute_comms%MPI_VAL,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          (domain%ide+2),(domain%kde+2),(domain%jde+2),domain%grid%ximages,one,domain%grid%yimages,one,one, &
                          xl, PETSC_NULL_INTEGER_ARRAY ,yl,da,ierr)

        
        call DMSetMatType(da,MATIS,ierr)

#ifdef _OPENACC
        has_cuda = PETSC_FALSE

        call PetscHasExternalPackage('cuda', has_cuda, ierr)

        if (has_cuda .eqv. PETSC_TRUE) then
            call DMSetMatType(da,MATMPIAIJCUSPARSE ,ierr)
            call DMSetVecType(da,VECMPICUDA,ierr)
        endif
#endif
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)

        call KSPSetDM(ksp(domain%nest_indx),da,ierr)

        call DMCreateLocalVector(da,localX,ierr)
        call DMCreateGlobalVector(da,b,ierr)
        call DMCreateMatrix(da,arr_A,ierr)

        call DMGetLocalToGlobalMapping(da,isltog,ierr)
        call MatSetLocalToGlobalMapping(arr_A,isltog,isltog,ierr)

    end subroutine

    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain

        integer :: ierr, i, j, k, i_s_bnd, i_e_bnd, j_s_bnd, j_e_bnd

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

        !i_s+hs, unless we are on global boundary, then i_s
        if (ims==ids) i_s = ids

        !i_e, unless we are on global boundary, then i_e+1
        if (ime==ide) i_e = ide

        !j_s+hs, unless we are on global boundary, then j_s
        if (jms==jds) j_s = jds
        
        !j_e, unless we are on global boundary, then j_e+1
        if (jme==jde) j_e = jde

        i_s_bnd = i_s
        i_e_bnd = i_e
        j_s_bnd = j_s
        j_e_bnd = j_e
        if (i_s==ids) then
            i_s_bnd = i_s + 1
        endif
        if (i_e==ide) then
            i_e_bnd = i_e - 1
        endif
        if (j_s==jds) then
            j_s_bnd = j_s + 1
        endif
        if (j_e==jde) then
            j_e_bnd = j_e - 1
        endif

        hs = domain%grid%halo_size
        if (.not.(allocated(dzdx))) then
            !call MPI_INIT(ierr)
            allocate(dzdx(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dzdy(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(jaco(i_s:i_e,k_s:k_e,j_s:j_e))

            allocate(dzdx_surf(i_s:i_e,j_s:j_e))
            allocate(dzdy_surf(i_s:i_e,j_s:j_e))

            allocate(sigma(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(dz_if(i_s:i_e,k_s:k_e+1,j_s:j_e))
            allocate(alpha(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate(div(i_s:i_e,k_s:k_e,j_s:j_e))
            allocate( xl( 1:domain%grid%ximages ))
            allocate( yl( 1:domain%grid%yimages ))
            xl = 0
            yl = 0

            dx = domain%dx
            
            ! Create GPU data for persistent arrays
            associate(advection_dz_var => domain%vars_3d(domain%var_indx(kVARS%advection_dz)%v)%data_3d, &
                      neighbor_terrain_var => domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d, &
                      dzdx_domain => domain%vars_3d(domain%var_indx(kVARS%dzdx)%v)%data_3d, &
                      dzdy_domain => domain%vars_3d(domain%var_indx(kVARS%dzdy)%v)%data_3d, &
                      jaco_domain => domain%vars_3d(domain%var_indx(kVARS%jacobian)%v)%data_3d)
            !$acc enter data create(dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, sigma, dz_if, alpha, div)

            !$acc data present(advection_dz_var, neighbor_terrain_var, dzdx_domain, dzdy_domain, jaco_domain, domain%vars_2d, domain%var_indx, dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, dz_if, sigma)

            !$acc parallel loop gang vector collapse(3)
            do j = j_s, j_e
                do k = k_s+1, k_e
                    do i = i_s, i_e
                        dz_if(i,k,j) = (advection_dz_var(i,k,j) + advection_dz_var(i,k-1,j))/2
                    end do
                end do
            end do            

            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s, i_e
                    dz_if(i,k_s,j) = advection_dz_var(i,k_s,j)
                    dz_if(i,k_e+1,j) = advection_dz_var(i,k_e,j)


                    dzdx_surf(i,j) = dzdx_domain(i,k_s,j)
                    dzdy_surf(i,j) = dzdy_domain(i,k_s,j)
                enddo
            end do
            
            !$acc parallel loop gang vector collapse(3)
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
            !$acc parallel loop gang vector collapse(2)
            do j = j_s, j_e
                do i = i_s_bnd, i_e_bnd
                    dzdx_surf(i,j) = (neighbor_terrain_var(i+1,j) - neighbor_terrain_var(i-1,j))/(2*dx)
                end do
            end do
            !$acc parallel loop gang vector collapse(2)
            do j = j_s_bnd, j_e_bnd
                do i = i_s, i_e
                    dzdy_surf(i,j) = (neighbor_terrain_var(i,j+1) - neighbor_terrain_var(i,j-1))/(2*dx)
                end do
            end do
            !$acc update host(dzdx, dzdy, jaco, dzdx_surf, dzdy_surf, dz_if, sigma)

            !$acc end data
            end associate

            !Calculate how global grid is decomposed for DMDA
            !subtract halo size from boundaries of each cell to get x/y extent
            if (domain%grid%yimg == 1) xl(domain%grid%ximg) = domain%grid%nx-hs*2
            if (domain%grid%ximg == 1) yl(domain%grid%yimg) = domain%grid%ny-hs*2
        
            !Wait for all images to contribute their dimension            
            call MPI_Allreduce(MPI_IN_PLACE,xl,domain%grid%ximages,MPI_INT,MPI_MAX,domain%compute_comms%MPI_VAL, ierr)
            call MPI_Allreduce(MPI_IN_PLACE,yl,domain%grid%yimages,MPI_INT,MPI_MAX,domain%compute_comms%MPI_VAL, ierr)

            !Add points to xy-edges to accomodate ghost-points of DMDA grid
            !cells at boundaries have 1 extra for ghost-point, and should also be corrected
            !to have the hs which was falsely removed added back on
            xl(1) = xl(1)+1+hs
            xl(domain%grid%ximages) = xl(domain%grid%ximages)+1+hs

            yl(1) = yl(1)+1+hs
            yl(domain%grid%yimages) = yl(domain%grid%yimages)+1+hs
            
        endif

    
    end subroutine

end module wind_iterative_old
