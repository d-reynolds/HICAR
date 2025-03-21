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
    !include 'petsc/finclude/petscksp.h'
    !include 'petsc/finclude/petscdm.h'
    !include 'petsc/finclude/petscdmda.h'
    
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    !use exchangeable_interface,   only : exchangeable_t
    use domain_interface,  only : domain_t
    !use options_interface, only : options_t
    !use grid_interface,    only : grid_t
    use petscksp
    use petscdm
    use petscdmda
    use array_utilities,      only : smooth_array_3d
    use icar_constants,       only : STD_OUT_PE
    use options_interface, only : options_t

    implicit none
    private
    public:: init_iter_winds_old, calc_iter_winds_old, finalize_iter_winds_old
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
    DM             da
    Vec            localX, b
    Mat            arr_A,arr_B

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

        integer k !, i_s, i_e, k_s, k_e, j_s, j_e
        
        PetscErrorCode ierr
        PetscInt       one, x_size, iteration
        PetscReal      norm, conv_tol
        PetscBool      reuse_bool
        Vec   x
        KSPConvergedReason reason

        update=.False.
        if (present(update_in)) update=update_in
                
        !Initialize div to be the initial divergence of the input wind field
        div = div_in(i_s:i_e,k_s:k_e,j_s:j_e) 
        one = 1
        alpha = alpha_in(i_s:i_e,k_s:k_e,j_s:j_e) 


        if (.not.(allocated(A_coef))) then
            ! Can't be called in module-level init function, since we first need alpha
            call initialize_coefs(domain)
        elseif (.not.( ALL(alpha==minval(alpha)) )) then
            call update_coefs(domain)
        endif

        call KSPSetComputeOperators(ksp(domain%nest_indx),ComputeMatrix,0,ierr)
        if (update_in) then
            call KSPSetReusePreconditioner(ksp(domain%nest_indx),PETSC_TRUE,ierr)
        endif
        ! call KSPGetOperators(ksp, A, PETSC_NULL_MAT, ierr) ! The second parameter is for the right-hand matrix, which can be NULL if not needed
        ! call MatIsSymmetric(A, PETSC_SMALL, isSymmetric, ierr);
        ! if (isSymmetric) then
        !     if(STD_OUT_PE) write(*,*) 'Matrix is symmetric'
        !     call KSPSetType(ksp,KSPPIPEGCR,ierr) !KSPIBCGS <-- this one tested to give fastest convergence...
        ! else 
        !     if(STD_OUT_PE) write(*,*) 'Matrix is not symmetric'
        !     call KSPSetType(ksp,KSPPIPEGCR,ierr) !KSPIBCGS <-- this one tested to give fastest convergence...
        ! endif

        call KSPSetComputeRHS(ksp(domain%nest_indx),ComputeRHS,0,ierr)

        call KSPSolve(ksp(domain%nest_indx),PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)
        call KSPGetSolution(ksp(domain%nest_indx),x,ierr)

        call KSPGetIterationNumber(ksp(domain%nest_indx), iteration, ierr)
        if(STD_OUT_PE) write(*,*) 'Solved PETSc after ',iteration,' iterations'
        
        !Subset global solution x to local grid so that we can access ghost-points
        call DMGlobalToLocal(da,x,INSERT_VALUES,localX,ierr)

        call DMDAVecGetArrayF90(da,localX,lambda, ierr)
        call calc_updated_winds(domain, lambda, update, adv_den)
        call DMDAVecRestoreArrayF90(da,localX,lambda, ierr)
        
        !Exchange u and v, since the outer points are not updated in above function
        call domain%halo%exch_var(domain%u,do_dqdt=.True.,corners=.True.)
        call domain%halo%exch_var(domain%v,do_dqdt=.True.,corners=.True.)
                                 
    end subroutine calc_iter_winds_old

    subroutine calc_updated_winds(domain,lambda,update,adv_den) !u, v, w, jaco_u,jaco_v,jaco_w,u_dzdx,v_dzdy,lambda, ids, ide, jds, jde)
        type(domain_t), intent(inout) :: domain
        !real, intent(inout), dimension(:,:,:)  :: u,v,w
        !real, intent(in), dimension(:,:,:)     :: jaco_u,jaco_v,jaco_w, u_dzdx, v_dzdy
        PetscScalar, intent(in), pointer       :: lambda(:,:,:)
        logical,     intent(in)                :: update
        logical,     intent(in)                :: adv_den


        real, allocatable, dimension(:,:,:)    :: u_dlambdz, v_dlambdz, u_temp, v_temp, lambda_too, rho, rho_u, rho_v
        integer k, i_start, i_end, j_start, j_end 

        i_start = i_s
        i_end   = i_e+1
        j_start = j_s
        j_end   = j_e+1

        allocate(u_temp(i_start:i_end,k_s-1:k_e+1,j_s:j_e))
        allocate(v_temp(i_s:i_e,k_s-1:k_e+1,j_start:j_end))

        allocate(u_dlambdz(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(v_dlambdz(i_s:i_e,k_s:k_e,j_start:j_end))
        
        allocate(rho(domain%ims:domain%ime,k_s:k_e,domain%jms:domain%jme))
        allocate(rho_u(i_start:i_end,k_s:k_e,j_s:j_e))
        allocate(rho_v(i_s:i_e,k_s:k_e,j_start:j_end))

        rho = 1.0
        rho_u = 1.0
        rho_v = 1.0
        
        if (adv_den) rho(domain%ims:domain%ime,:,domain%jms:domain%jme)=domain%density%data_3d(domain%ims:domain%ime,:,domain%jms:domain%jme)
        
        if (i_s==domain%grid%ids .and. i_e==domain%grid%ide) then
            rho_u(i_start+1:i_end-1,:,j_s:j_e) = 0.5*(rho(i_start+1:i_end-1,:,j_s:j_e) + rho(i_start:i_end-2,:,j_s:j_e))
            rho_u(i_end,:,j_s:j_e) = rho(i_end-1,:,j_s:j_e)
            rho_u(i_start,:,j_s:j_e) = rho(i_start,:,j_s:j_e)
        else if (i_s==domain%grid%ids) then
            rho_u(i_start+1:i_end,:,j_s:j_e) = 0.5*(rho(i_start+1:i_end,:,j_s:j_e) + rho(i_start:i_end-1,:,j_s:j_e))
            rho_u(i_start,:,j_s:j_e) = rho(i_start,:,j_s:j_e)
        else if (i_e==domain%grid%ide) then
            rho_u(i_start:i_end-1,:,j_s:j_e) = 0.5*(rho(i_start:i_end-1,:,j_s:j_e) + rho(i_start-1:i_end-2,:,j_s:j_e))
            rho_u(i_end,:,j_s:j_e) = rho(i_end-1,:,j_s:j_e)
        else
            rho_u(i_start:i_end,:,j_s:j_e) = 0.5*(rho(i_start:i_end,:,j_s:j_e) + rho(i_start-1:i_end-1,:,j_s:j_e))
        endif
        
        if (j_s==domain%grid%jds .and. j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start+1:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end-1) + rho(i_s:i_e,:,j_start:j_end-2))
            rho_v(i_s:i_e,:,j_start) = rho(i_s:i_e,:,j_start)
            rho_v(i_s:i_e,:,j_end) = rho(i_s:i_e,:,j_end-1)
        else if (j_s==domain%grid%jds) then
            rho_v(i_s:i_e,:,j_start+1:j_end) = 0.5*(rho(i_s:i_e,:,j_start+1:j_end) + rho(i_s:i_e,:,j_start:j_end-1))
            rho_v(i_s:i_e,:,j_start) = rho(i_s:i_e,:,j_start)
        else if (j_e==domain%grid%jde) then
            rho_v(i_s:i_e,:,j_start:j_end-1) = 0.5*(rho(i_s:i_e,:,j_start:j_end-1) + rho(i_s:i_e,:,j_start-1:j_end-2))
            rho_v(i_s:i_e,:,j_end) = rho(i_s:i_e,:,j_end-1)
        else
            rho_v(i_s:i_e,:,j_start:j_end) = 0.5*(rho(i_s:i_e,:,j_start:j_end) + rho(i_s:i_e,:,j_start-1:j_end-1))
        endif
        
        !stager lambda to u grid
        u_temp = (lambda(i_start:i_end,k_s-1:k_e+1,j_s:j_e) + lambda(i_start-1:i_end-1,k_s-1:k_e+1,j_s:j_e)) / 2 

        !stager lambda to v grid
        v_temp = (lambda(i_s:i_e,k_s-1:k_e+1,j_start:j_end) + lambda(i_s:i_e,k_s-1:k_e+1,j_start-1:j_end-1)) / 2 

        !divide dz differennces by dz. Note that dz will be horizontally constant
        !STRUCT_DIFF: u_dlambdz and v_dlambdz are are calculated using "normal" numerics, and not the derivative for
        !             a stretched vertical grid.
        do k=k_s,k_e
            u_dlambdz(:,k,:) = u_temp(:,k+1,:) - u_temp(:,k-1,:)
            v_dlambdz(:,k,:) = v_temp(:,k+1,:) - v_temp(:,k-1,:)
        
            u_dlambdz(:,k,:) = u_dlambdz(:,k,:)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
            v_dlambdz(:,k,:) = v_dlambdz(:,k,:)/(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))
        enddo
        
        !u_dlambdz(:,k_s,:) = -(u_temp(:,k_s+2,:)*sigma(i_s,k_s+1,j_s)**2) + &
        !                    u_temp(:,k_s+1,:)*(sigma(i_s,k_s+1,j_s)+1)**2 - u_temp(:,k_s,:)*(2*sigma(i_s,k_s+1,j_s)+1)
        !v_dlambdz(:,k_s,:) = -(v_temp(:,k_s+2,:)*sigma(i_s,k_s+1,j_s)**2) + &
        !                    v_temp(:,k_s+1,:)*(sigma(i_s,k_s+1,j_s)+1)**2 - v_temp(:,k_s,:)*(2*sigma(i_s,k_s+1,j_s)+1)
        !
        !u_dlambdz(:,k_s,:) = u_dlambdz(:,k_s,:)/(dz_if(i_s,k_s+1,j_s)*(sigma(i_s,k_s+1,j_s)+1))
        !v_dlambdz(:,k_s,:) = v_dlambdz(:,k_s,:)/(dz_if(i_s,k_s+1,j_s)*(sigma(i_s,k_s+1,j_s)+1))
        
        !PETSc arrays are zero-indexed
        
        ! STRUCT_DIFF: divide by jacobian again when it is an "update" call
        domain%u%dqdt_3d(i_start:i_end,:,j_s:j_e) = domain%u%dqdt_3d(i_start:i_end,:,j_s:j_e) + &
                                                        0.5*((lambda(i_start:i_end,k_s:k_e,j_s:j_e) - &
                                                        lambda(i_start-1:i_end-1,k_s:k_e,j_s:j_e))/dx - &
        (1/domain%jacobian_u(i_start:i_end,:,j_s:j_e))*domain%dzdx_u(i_start:i_end,:,j_s:j_e)*(u_dlambdz))/(rho_u(i_start:i_end,:,j_s:j_e)*domain%jacobian_u(i_start:i_end,:,j_s:j_e))
        
        domain%v%dqdt_3d(i_s:i_e,:,j_start:j_end) = domain%v%dqdt_3d(i_s:i_e,:,j_start:j_end) + &
                                                        0.5*((lambda(i_s:i_e,k_s:k_e,j_start:j_end) - &
                                                        lambda(i_s:i_e,k_s:k_e,j_start-1:j_end-1))/dx - &
        (1/domain%jacobian_v(i_s:i_e,:,j_start:j_end))*domain%dzdy_v(i_s:i_e,:,j_start:j_end)*(v_dlambdz))/(rho_v(i_s:i_e,:,j_start:j_end)*domain%jacobian_v(i_s:i_e,:,j_start:j_end))
        

    end subroutine calc_updated_winds

    subroutine ComputeRHS(ksp,vec_b,dummy,ierr)
        implicit none
        
        PetscErrorCode  ierr
        KSP ksp
        Vec vec_b
        integer dummy(*)
        DM      ::       dm

        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs
        DMDALocalInfo       :: info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar,pointer :: barray(:,:,:)

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

        call DMDAVecRestoreArrayF90(dm,vec_b,barray, ierr)
    end subroutine ComputeRHS

    subroutine ComputeInitialGuess(ksp,vec_b,ctx,ierr)
        implicit none
        PetscErrorCode  ierr
        KSP ksp
        PetscInt ctx(*)
        Vec vec_b
        PetscScalar  i_guess

        i_guess = 0.0

        call VecSet(vec_b,i_guess,ierr)
    end subroutine ComputeInitialGuess

    subroutine ComputeMatrix(ksp,array_A,array_B)
        implicit none
        KSP ksp
        Mat array_A,array_B
        ! integer dummy(*)

        PetscErrorCode  ierr
        DM             da
        PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,i1,i2,i6,i15
        DMDALocalInfo  info(DMDA_LOCAL_INFO_SIZE)
        PetscScalar    v(15)
        MatStencil     row(4),col(4,15),gnd_col(4,10),top_col(4,2)

        real :: denom

        i1 = 1
        i2 = 2
        i6 = 10
        i15 = 15

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
        
        do j=ys,(ys+ym-1)
        do k=zs,(zs+zm-1)
            do i=xs,(xs+xm-1)
                row(MatStencil_i) = i
                row(MatStencil_j) = k
                row(MatStencil_k) = j
                if (i.le.0 .or. j.le.0 .or. &
                    i.ge.mx-1 .or. j.ge.my-1) then
                    v(1) = 1.0
                    call MatSetValuesStencil(array_B,i1,row,i1,row,v,INSERT_VALUES, ierr)
                else if (k.eq.mz-1) then
                    !k
                    v(1) = 1/dz_if(i,k,j)
                    top_col(MatStencil_i,1) = i
                    top_col(MatStencil_j,1) = k
                    top_col(MatStencil_k,1) = j
                    !k - 1
                    v(2) = -1/dz_if(i,k,j)
                    top_col(MatStencil_i,2) = i
                    top_col(MatStencil_j,2) = k-1
                    top_col(MatStencil_k,2) = j
                    call MatSetValuesStencil(array_B,i1,row,i2,top_col,v,INSERT_VALUES, ierr)

                else if (k.eq.0) then
                    denom = 2*(dzdx_surf(i,j)**2 + dzdy_surf(i,j)**2 + alpha(i,1,j)**2)/(jaco(i,1,j))
                    !k
                    v(2) = -1/dz_if(i,k+1,j)
                    gnd_col(MatStencil_i,2) = i
                    gnd_col(MatStencil_j,2) = k
                    gnd_col(MatStencil_k,2) = j
                    !k + 1
                    v(1) = 1/dz_if(i,k+1,j)
                    gnd_col(MatStencil_i,1) = i
                    gnd_col(MatStencil_j,1) = k+1
                    gnd_col(MatStencil_k,1) = j

                    !Reminder: The equation for the lower BC has all of the lateral derivatives (dlambdadx, etc) NEGATIVE, hence opposite signs below
                    !If we are on left most border
                    !i - 1, k + 1
                    v(3) = dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,3) = i-1
                    gnd_col(MatStencil_j,3) = k
                    gnd_col(MatStencil_k,3) = j
                    !i + 1, k + 1
                    v(4) = -dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,4) = i+1
                    gnd_col(MatStencil_j,4) = k
                    gnd_col(MatStencil_k,4) = j
                    !i - 1, k + 1
                    v(5) = dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,5) = i-1
                    gnd_col(MatStencil_j,5) = k+1
                    gnd_col(MatStencil_k,5) = j
                    !i + 1, k + 1
                    v(6) = -dzdx_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,6) = i+1
                    gnd_col(MatStencil_j,6) = k+1
                    gnd_col(MatStencil_k,6) = j
                    !j - 1, k + 1
                    v(7) = dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,7) = i
                    gnd_col(MatStencil_j,7) = k
                    gnd_col(MatStencil_k,7) = j-1
                    !j + 1, k + 1
                    v(8) = -dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,8) = i
                    gnd_col(MatStencil_j,8) = k
                    gnd_col(MatStencil_k,8) = j+1
                    !j - 1, k + 1
                    v(9) = dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,9) = i
                    gnd_col(MatStencil_j,9) = k+1
                    gnd_col(MatStencil_k,9) = j-1
                    !j + 1, k + 1
                    v(10) = -dzdy_surf(i,j)/(denom*2*dx)
                    gnd_col(MatStencil_i,10) = i
                    gnd_col(MatStencil_j,10) = k+1
                    gnd_col(MatStencil_k,10) = j+1
                    call MatSetValuesStencil(array_B,i1,row,i6,gnd_col,v,INSERT_VALUES, ierr)
                else
                    !j - 1, k - 1
                    v(1) = O_coef(i,k,j)
                    col(MatStencil_i,1) = i
                    col(MatStencil_j,1) = k-1
                    col(MatStencil_k,1) = j-1
                    !i - 1, k - 1
                    v(2) = K_coef(i,k,j)
                    col(MatStencil_i,2) = i-1
                    col(MatStencil_j,2) = k-1
                    col(MatStencil_k,2) = j
                    !j - 1, k + 1
                    v(3) = M_coef(i,k,j)
                    col(MatStencil_i,3) = i
                    col(MatStencil_j,3) = k+1
                    col(MatStencil_k,3) = j-1
                    !i - 1, k + 1
                    v(4) = I_coef(i,k,j)
                    col(MatStencil_i,4) = i-1
                    col(MatStencil_j,4) = k+1
                    col(MatStencil_k,4) = j
                    !k - 1
                    v(5) = C_coef(i,k,j)
                    col(MatStencil_i,5) = i
                    col(MatStencil_j,5) = k-1
                    col(MatStencil_k,5) = j
                    !j - 1
                    v(6) = G_coef(i,k,j)
                    col(MatStencil_i,6) = i
                    col(MatStencil_j,6) = k
                    col(MatStencil_k,6) = j-1
                    !i - 1
                    v(7) = E_coef(i,k,j)
                    col(MatStencil_i,7) = i-1
                    col(MatStencil_j,7) = k
                    col(MatStencil_k,7) = j
                    !Center
                    v(8) = A_coef(i,k,j)
                    col(MatStencil_i,8) = i
                    col(MatStencil_j,8) = k
                    col(MatStencil_k,8) = j
                    !i + 1
                    v(9) = D_coef(i,k,j)
                    col(MatStencil_i,9) = i+1
                    col(MatStencil_j,9) = k
                    col(MatStencil_k,9) = j
                    !j + 1
                    v(10) = F_coef(i,k,j)
                    col(MatStencil_i,10) = i
                    col(MatStencil_j,10) = k
                    col(MatStencil_k,10) = j+1
                    !k + 1
                    v(11) = B_coef(i,k,j)
                    col(MatStencil_i,11) = i
                    col(MatStencil_j,11) = k+1
                    col(MatStencil_k,11) = j
                    !i + 1, k + 1
                    v(12) = H_coef(i,k,j)
                    col(MatStencil_i,12) = i+1
                    col(MatStencil_j,12) = k+1
                    col(MatStencil_k,12) = j
                    !j + 1, k + 1
                    v(13) = L_coef(i,k,j)
                    col(MatStencil_i,13) = i
                    col(MatStencil_j,13) = k+1
                    col(MatStencil_k,13) = j+1
                    !i + 1, k - 1
                    v(14) = J_coef(i,k,j)
                    col(MatStencil_i,14) = i+1
                    col(MatStencil_j,14) = k-1
                    col(MatStencil_k,14) = j
                    !j + 1, k - 1
                    v(15) = N_coef(i,k,j)
                    col(MatStencil_i,15) = i
                    col(MatStencil_j,15) = k-1
                    col(MatStencil_k,15) = j+1
                    call MatSetValuesStencil(array_B,i1,row,i15,col,v,INSERT_VALUES, ierr)
                endif
            enddo
        enddo
        enddo

        call MatAssemblyBegin(array_B,MAT_FINAL_ASSEMBLY,ierr)
        call MatAssemblyEnd(array_B,MAT_FINAL_ASSEMBLY,ierr)
        if (array_A .ne. array_B) then
         call MatAssemblyBegin(array_A,MAT_FINAL_ASSEMBLY,ierr)
         call MatAssemblyEnd(array_A,MAT_FINAL_ASSEMBLY,ierr)
        endif

    end subroutine ComputeMatrix
    
    

    subroutine initialize_coefs(domain)
        implicit none
        type(domain_t), intent(in) :: domain
        
        real :: sig_i, mixed_denom
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
            mixed_denom = 2*domain%dx*(dz_if(i_s,k+1,j_s)+dz_if(i_s,k,j_s))!*(sig_i+sig_i**2)
            
            H_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:)+domain%dzdx_u(i_s+1:i_e+1,k,j_s:j_e)/domain%jacobian_u(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            I_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:)+domain%dzdx_u(i_s:i_e,k,j_s:j_e)/domain%jacobian_u(i_s:i_e,k,j_s:j_e))/mixed_denom
            L_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:)+domain%dzdy_v(i_s:i_e,k,j_s+1:j_e+1)/domain%jacobian_v(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            M_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:)+domain%dzdy_v(i_s:i_e,k,j_s:j_e)/domain%jacobian_v(i_s:i_e,k,j_s:j_e))/mixed_denom

            J_coef(:,k,:) =  (dzdx(:,k,:)/jaco(:,k,:)+domain%dzdx_u(i_s+1:i_e+1,k,j_s:j_e)/domain%jacobian_u(i_s+1:i_e+1,k,j_s:j_e))/mixed_denom
            K_coef(:,k,:) = -(dzdx(:,k,:)/jaco(:,k,:)+domain%dzdx_u(i_s:i_e,k,j_s:j_e)/domain%jacobian_u(i_s:i_e,k,j_s:j_e))/mixed_denom
            N_coef(:,k,:) =  (dzdy(:,k,:)/jaco(:,k,:)+domain%dzdy_v(i_s:i_e,k,j_s+1:j_e+1)/domain%jacobian_v(i_s:i_e,k,j_s+1:j_e+1))/mixed_denom
            O_coef(:,k,:) = -(dzdy(:,k,:)/jaco(:,k,:)+domain%dzdy_v(i_s:i_e,k,j_s:j_e)/domain%jacobian_v(i_s:i_e,k,j_s:j_e))/mixed_denom

        enddo


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
        
        mixed_denom = (dz_if(:,k_s+1:k_e+1,:)+dz_if(:,k_s:k_e,:))*2*domain%dx

        do k = k_s,k_e
            if (k == k_s) then
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%advection_dz(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*dzdx_surf/jaco(i_s:i_e,k,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%advection_dz(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*dzdy_surf/jaco(i_s:i_e,k,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%advection_dz(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + alpha(:,k,:)**2/jaco(i_s:i_e,k,j_s:j_e)
            else if (k == k_e) then
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%advection_dz(i_s,k,j_s) + &
                        dzdx(:,k,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k,j_s))
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%advection_dz(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%advection_dz(i_s,k,j_s) + &
                        dzdy(:,k,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k,j_s))
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%advection_dz(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%advection_dz(i_s,k,j_s) + &
                        alpha(:,k,:)**2*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k,j_s))
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%advection_dz(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
            else
                !Terms for dzdx
                M_up = dzdx(:,k,:)*(dzdx(:,k,:)*domain%advection_dz(i_s,k+1,j_s) + &
                        dzdx(:,k+1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = dzdx(:,k,:)*(dzdx(:,k,:)*domain%advection_dz(i_s,k-1,j_s) + &
                        dzdx(:,k-1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
                !Terms for dzdy
                M_up = M_up + dzdy(:,k,:)*(dzdy(:,k,:)*domain%advection_dz(i_s,k+1,j_s) + &
                        dzdy(:,k+1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + dzdy(:,k,:)*(dzdy(:,k,:)*domain%advection_dz(i_s,k-1,j_s) + &
                        dzdy(:,k-1,:)*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
                !Terms for alpha
                M_up = M_up + (alpha(:,k,:)**2*domain%advection_dz(i_s,k+1,j_s) + &
                        alpha(:,k+1,:)**2*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k+1,j_s))/domain%jacobian_w(i_s:i_e,k,j_s:j_e)
                M_dwn = M_dwn + (alpha(:,k,:)**2*domain%advection_dz(i_s,k-1,j_s) + &
                        alpha(:,k-1,:)**2*domain%advection_dz(i_s,k,j_s))/ &
                        (domain%advection_dz(i_s,k,j_s)+domain%advection_dz(i_s,k-1,j_s))/domain%jacobian_w(i_s:i_e,k-1,j_s:j_e)
            endif
            B_coef(:,k,:) = M_up/(jaco(:,k,:)*domain%advection_dz(i_s,k,j_s)*dz_if(i_s,k+1,j_s))
            C_coef(:,k,:) = M_dwn/(jaco(:,k,:)*domain%advection_dz(i_s,k,j_s)*dz_if(i_s,k,j_s))
        enddo
                
        A_coef = -4/(domain%dx**2) - B_coef - C_coef
        
        X_coef = -((domain%dzdx_u(i_s+1:i_e+1,:,j_s:j_e)/domain%jacobian_u(i_s+1:i_e+1,:,j_s:j_e) - domain%dzdx_u(i_s:i_e,:,j_s:j_e)/domain%jacobian_u(i_s:i_e,:,j_s:j_e)) + &
                   (domain%dzdy_v(i_s:i_e,:,j_s+1:j_e+1)/domain%jacobian_v(i_s:i_e,:,j_s+1:j_e+1) - domain%dzdy_v(i_s:i_e,:,j_s:j_e)/domain%jacobian_v(i_s:i_e,:,j_s:j_e)))/mixed_denom
                   
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
        if (allocated(div)) deallocate(div)
        if (allocated(dz_if)) deallocate(dz_if)
        if (allocated(jaco)) deallocate(jaco)
        if (allocated(dzdx)) deallocate(dzdx)
        if (allocated(dzdy)) deallocate(dzdy)
        if (allocated(dzdx_surf)) deallocate(dzdx_surf)
        if (allocated(dzdy_surf)) deallocate(dzdy_surf)
        if (allocated(sigma)) deallocate(sigma)
        if (allocated(alpha)) deallocate(alpha)
        if (allocated(xl)) deallocate(xl)
        if (allocated(yl)) deallocate(yl)
        
        call VecDestroy(localX,ierr)
        call VecDestroy(b, ierr)
        call MatDestroy(arr_A,ierr)
        call MatDestroy(arr_B,ierr)
        call DMDestroy(da,ierr)
        ! call KSPDestroy(ksp,ierr)
        ! call PetscFinalize(ierr)

    end subroutine

    subroutine init_petsc_comms(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options
        PetscErrorCode ierr

        integer :: i

        PETSC_COMM_WORLD = domain%compute_comms%MPI_VAL
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

        allocate(ksp(options%general%nests))

        do i=1,options%general%nests
            ! if (options%physics%windtype == kITERATIVE_WINDS .or. &
            !     options%physics%windtype == kLINEAR_ITERATIVE_WINDS) then
                call KSPCreate(domain%compute_comms%MPI_VAL,ksp(i),ierr)
                call KSPSetFromOptions(ksp(i),ierr)
        
                call KSPSetType(ksp(i),KSPPIPEGCR,ierr) !KSPPIPEFCG <-- this one does not converge.
                                                !KSPPIPEGCR <-- this one tested to give fastest convergence...
                call KSPSetComputeInitialGuess(ksp(i),ComputeInitialGuess,0,ierr)
            ! endif
        enddo
        initialized = .True.

        if (ierr .ne. 0) then
            print*,'Unable to initialize PETSc'
            stop
        endif

    end subroutine

    subroutine init_iter_winds_old(domain, options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t), intent(in) :: options

        PetscInt       one, iter
        PetscErrorCode ierr

        if (.not.(initialized)) then
            call init_petsc_comms(domain, options)
        endif


        ! call KSPSetReusePreconditioner(ksp,PETSC_TRUE,ierr)
        ! call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)

        call init_module_vars(domain)

        one = 1

        call DMDACreate3d(domain%compute_comms%MPI_VAL,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX, &
                          (domain%ide+2),(domain%kde+2),(domain%jde+2),domain%grid%ximages,one,domain%grid%yimages,one,one, &
                          xl, PETSC_NULL_INTEGER_ARRAY ,yl,da,ierr)
        
        call DMSetFromOptions(da,ierr)
        call DMSetUp(da,ierr)

        call KSPSetDM(ksp(domain%nest_indx),da,ierr)
        ! call KSPSetDMActive(ksp,PETSC_FALSE, ierr)

        call DMCreateLocalVector(da,localX,ierr)
        call DMCreateGlobalVector(da,b,ierr)
        call DMCreateMatrix(da,arr_A,ierr)
        call DMCreateMatrix(da,arr_B,ierr)
    end subroutine

    subroutine init_module_vars(domain)
        implicit none
        type(domain_t), intent(in) :: domain

        integer :: ierr

        i_s = domain%its
        i_e = domain%ite
        k_s = domain%kts  
        k_e = domain%kte  
        j_s = domain%jts
        j_e = domain%jte

        !i_s+hs, unless we are on global boundary, then i_s
        if (domain%grid%ims==domain%grid%ids) i_s = domain%grid%ids
        
        !i_e, unless we are on global boundary, then i_e+1
        if (domain%grid%ime==domain%grid%ide) i_e = domain%grid%ide
        
        !j_s+hs, unless we are on global boundary, then j_s
        if (domain%grid%jms==domain%grid%jds) j_s = domain%grid%jds
        
        !j_e, unless we are on global boundary, then j_e+1
        if (domain%grid%jme==domain%grid%jde) j_e = domain%grid%jde

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
            dzdx  = domain%dzdx%data_3d(i_s:i_e,k_s:k_e,j_s:j_e) 
            dzdy  = domain%dzdy%data_3d(i_s:i_e,k_s:k_e,j_s:j_e)
            jaco = domain%jacobian(i_s:i_e,k_s:k_e,j_s:j_e)
            
            dzdx_surf = 0.1
            dzdy_surf = 0.1
            dzdx_surf(i_s+1:i_e-1,j_s:j_e) = (domain%neighbor_terrain(i_s+2:i_e,j_s:j_e)-domain%neighbor_terrain(i_s:i_e-2,j_s:j_e))/(2*dx)
            dzdy_surf(i_s:i_e,j_s+1:j_e-1) = (domain%neighbor_terrain(i_s:i_e,j_s+2:j_e)-domain%neighbor_terrain(i_s:i_e,j_s:j_e-2))/(2*dx)
                          
            dzdx_surf(i_s,j_s:j_e) = dzdx_surf(i_s+1,j_s:j_e)
            dzdx_surf(i_e,j_s:j_e) = dzdx_surf(i_e-1,j_s:j_e)
            dzdy_surf(i_s:i_e,j_s) = dzdy_surf(i_s:i_e,j_s+1)
            dzdy_surf(i_s:i_e,j_e) = dzdy_surf(i_s:i_e,j_e-1)
            
            dzdy_surf(i_s:i_e,j_s:j_e) = dzdy(i_s:i_e,1,j_s:j_e)
            dzdx_surf(i_s:i_e,j_s:j_e) = dzdx(i_s:i_e,1,j_s:j_e)

            
            dz_if(:,k_s+1:k_e,:) = (domain%advection_dz(i_s:i_e,k_s+1:k_e,j_s:j_e) + &
                                   domain%advection_dz(i_s:i_e,k_s:k_e-1,j_s:j_e))/2
            dz_if(:,k_s,:) = domain%advection_dz(i_s:i_e,k_s,j_s:j_e)
            dz_if(:,k_e+1,:) = domain%advection_dz(i_s:i_e,k_e,j_s:j_e)
            sigma = dz_if(:,k_s:k_e,:)/dz_if(:,k_s+1:k_e+1,:)
            
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
