!=================================================================================================================
!module_bl_ysu.F was originally copied from ./phys/module_bl_ysu.F from WRF version 3.8.1.
!Laura D. Fowler (laura@ucar.edu) / 2016-10-26.

!modifications to sourcecode for MPAS:
!   * calculated the dry hydrostatic pressure using the dry air density.
!   * added outputs of the vertical diffusivity coefficients.
!     Laura D. Fowler (laura@ucar.edu) / 2016-10-26.

!modifications to code for HICAR integration 2025-09-30 by Dylan Reynolds (dylan.reynolds@epfl.ch)
!     * refactoring of code to run over full 2d/3d arrays instead of iterating over j-loops
!     * fusing of various loops, introducing asynchronous OpenACC operations
!
!=================================================================================================================
!WRF:model_layer:physics
!
!
!
!
!
!
!
module module_bl_ysu_gpu
contains
!
!
!-------------------------------------------------------------------------------
!
   subroutine ysu_gpu(u3d,v3d,th3d,t3d,qv3d,qc3d,qi3d,p3d,p3di,pi3d,               &
                  rthblten,                                    &
                  rqvblten,rqcblten,rqiblten,flag_qi,                          &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
                  dz8w,z8w,psfc,                                                   &
                  znu,znw,mut,p_top,                                           &
                  znt,ust,hpbl,psim,psih,                                      &
                  xland,hfx,qfx,wspd,br,                                       &
                  dt,kpbl2d,                                                   &
                  exch_h,                                                      &
                  ! wstar,delta,                                                 &
                  ! hgamt,hgamq,entpbl,                                          &
                  u10,v10,                                                     &
                  ! uoce,voce, tke_pbl,                                                   &
                  rthraten,ysu_topdown_pblmix,                         &
                  ! ctopo,ctopo2,                                                &
                  ids,ide, jds,jde, kds,kde,                                   &
                  ims,ime, jms,jme, kms,kme,                                   &
                  its,ite, jts,jte, kts,kte,                                   &
                !optional
                  regime,rho,rublten,rvblten                                   &
#if defined(mpas)
                !MPAS specific optional arguments for additional diagnostics:
                  ,kzhout,kzmout,kzqout                                    &
#endif
                 )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!               (note: if P_QI<PARAM_FIRST_SCALAR this should be zero filled)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- rr3d        3d dry air density (kg/m^3)
!-- rublten     u tendency due to
!               pbl parameterization (m/s/s)
!-- rvblten     v tendency due to
!               pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to
!               pbl parameterization (K/s)
!-- rqvblten    qv tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqcblten    qc tendency due to
!               pbl parameterization (kg/kg/s)
!-- rqiblten    qi tendency due to
!               pbl parameterization (kg/kg/s)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rovcp       r/cp
!-- rd          gas constant for dry air (j/kg/k)
!-- rovg        r/g
!-- dz8w        dz between full levels (m)
!-- z8w         height of full levels (m)
!-- xlv         latent heat of vaporization (j/kg)
!-- rv          gas constant for water vapor (j/kg/k)
!-- psfc        pressure at the surface (pa)
!-- znt         roughness length (m)
!-- ust         u* in similarity theory (m/s)
!-- hpbl        pbl height (m)
!-- psim        similarity stability function for momentum
!-- psih        similarity stability function for heat
!-- xland       land mask (1 for land, 2 for water)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- wspd        wind speed at lowest model level (m/s)
!-- u10         u-wind speed at 10 m (m/s)
!-- v10         v-wind speed at 10 m (m/s)
!-- uoce        sea surface zonal currents (m s-1)
!-- voce        sea surface meridional currents (m s-1)
!-- br          bulk richardson number in surface layer
!-- dt          time step (s)
!-- rvovrd      r_v divided by r_d (dimensionless)
!-- ep1         constant for virtual temperature (r_v/r_d - 1)
!-- ep2         constant for specific humidity calculation
!-- karman      von karman constant
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------------------
!
   integer,parameter ::  ndiff = 3
   real,parameter    ::  rcl = 1.0
!
   integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte

   integer,  intent(in)      ::      ysu_topdown_pblmix
!
   real,     intent(in   )   ::      dt,cp,g,rovcp,rovg,rd,xlv,rv
!
   real,     intent(in )     ::      ep1,ep2,karman
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                          qv3d, &
                                                                         qc3d, &
                                                                         qi3d, &
                                                                          p3d, &
                                                                         pi3d, &
                                                                         th3d, &
                                                                          t3d, &
                                                                         dz8w, &
                                                                     rthraten

   real,     dimension( ims:ime, kms:kme+1, jms:jme )                          , &
             intent(in   )   ::                                          z8w

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                          p3di
!
   real,  optional,   dimension( ims:ime, kms:kme, jms:jme )                 , &
             intent(inout)   ::                                       rublten, &
                                                                      rvblten

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(inout)   ::                                      rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten

  !  real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
  !            intent(inout)   ::                                       tke_pbl
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(inout)   ::                                        exch_h
  !  real,     dimension( ims:ime, jms:jme )                                   , &
  !            intent(inout)   ::                                         wstar
  !  real,     dimension( ims:ime, jms:jme )                                   , &
  !            intent(inout)   ::                                         delta
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(inout)   ::                                           u10, &
                                                                          v10
  !  real,     dimension( ims:ime, jms:jme )                                   , &
  !            intent(in   )   ::                                          uoce, &
  !                                                                        voce
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                           hfx, &
                                                                          qfx, &
                                                                           br, &
                                                                         psfc
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(in   )   ::                                                &
                                                                         psim, &
                                                                         psih
   real,     dimension( ims:ime, jms:jme )                                   , &
             intent(inout)   ::                                           znt, &
                                                                          ust, &
                                                                        ! hgamt, &
                                                                        ! hgamq, &
                                                                      !  entpbl, &
                                                                         hpbl, &
                                                                          wspd
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in   )   ::                                           u3d, &
                                                                          v3d
!
   integer,  dimension( ims:ime, jms:jme )                                   , &
             intent(out  )   ::                                        kpbl2d

   integer,  dimension( ims:ime, jms:jme )                                   , &
             intent(in  )   ::                                        xland

   logical,  intent(in)      ::                                       flag_qi
!
!optional
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             optional                                                        , &
             intent(inout)   ::                                        regime
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             optional                                                        , &
             intent(inout)   ::                                       rqiblten
!
   real,     dimension( kms:kme )                                            , &
             optional                                                        , &
             intent(in   )   ::                                           znu, &
                                                                          znw
!
   real,     dimension( ims:ime, jms:jme )                                   , &
             optional                                                        , &
             intent(in   )   ::                                           mut
!
   real,     optional, intent(in   )   ::                               p_top
!
  !  real,     dimension( ims:ime, jms:jme )                                   , &
  !            optional                                                        , &
  !            intent(in   )   ::                                         ctopo, &
  !                                                                      ctopo2
!local
   integer ::  i,j,k
   real,     dimension( its:ite, kts:kte, jts:jte, ndiff )  ::                 qtnp

   real,     dimension( its:ite, kts:kte, jts:jte )  ::                            pdh, dpdh
   real,     dimension( its:ite, kts:kte+1, jts:jte )  ::                         pdhi
   real,     dimension( its:ite, jts:jte )  ::                                          &
                                                                        dusfc, &
                                                                        dvsfc, &
                                                                        dtsfc, &
                                                                        dqsfc
   real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: rho
   real:: rho_d

   real,     dimension( ims:ime, jms:jme ) ::                                  hgamt, &
                                                                 hgamq, &
                                                                 entpbl, wstar, delta, ctopo, ctopo2, uoce, voce
   real,     dimension( ims:ime, kms:kme, jms:jme ) ::                                       tke_pbl

   real,     dimension( its:ite, jts:jte )            ::                           hol
   real,     dimension( its:ite, kts:kte+1, jts:jte ) ::                            zq
!
   real,     dimension( its:ite, kts:kte, jts:jte )   ::                                &
                                                               thx,thvx,thlix, &
                                                                          del, &
                                                                          dza, &
                                                                          dzq, &
                                                                         xkzo, &
                                                                           za
!
   real,    dimension( its:ite, jts:jte )             ::                                &
                                                                         rhox, &
                                                                       govrth, &
                                                                  zl1,thermal, &
                                                                       wscale, &
                                                                    brdn,brup, &
                                                                    phim,phih, &
                                                                        prpbl, &
                                                              wspd1,thermalli
!
   real,    dimension( its:ite, kts:kte, jts:jte )    ::                     xkzm,xkzh, &
                                                                        f1,f2, &
                                                                        r1,r2, &
                                                                        ad,au, &
                                                                           cu, &
                                                                           al, &
                                                                         xkzq, &
                                                                         zfac, &
                                                                        rhox2, &
                                                                       hgamt2
!
  !  real,    dimension( ims:ime )                                             , &
  !           intent(inout)    ::                                           u10, &
  !                                                                         v10
   real,    dimension( its:ite, jts:jte )    ::                                         &
                                                                         brcr, &
                                                                        sflux, &
                                                                         zol1, &
                                                                    brcr_sbro
!
   real,    dimension( its:ite, kts:kte, jts:jte, ndiff)  ::                     r3,f3
   integer, dimension( its:ite, jts:jte )             ::                  kpbl,kpblold
!
   logical, dimension( its:ite, jts:jte )             ::                        pblflg, &
                                                                       sfcflg, &
                                                                       stable, &
                                                                     cloudflg

   logical                                   ::                     definebrup
!
   integer ::  n,p,l,ic,is,kk
   integer ::  klpbl, ktrace1, ktrace2, ktrace3
!
!
   real    ::  dt2,rdt,spdk2,fm,fh,hol1,gamfac,vpert,prnum,prnum0
   real    ::  ss,ri,qmean,tmean,alph,chi,zk,rl2,dk,sri
   real    ::  brint,dtodsd,dtodsu,rdz,dsdzt,dsdzq,dsdz2,rlamdz
   real    ::  utend,vtend,ttend,qtend
   real    ::  dtstep,govrthv,c
   real    ::  cont, conq, conw, conwrc
   real,parameter    ::  xkzminm = 0.1,xkzminh = 0.01
   real,parameter    ::  xkzmin = 0.01,xkzmax = 1000.,rimin = -100.
   real,parameter    ::  rlam = 30.,prmin = 0.25,prmax = 4.
   real,parameter    ::  brcr_ub = 0.0,brcr_sb = 0.25,cori = 1.e-4
   real,parameter    ::  afac = 6.8,bfac = 6.8,pfac = 2.0,pfac_q = 2.0
   real,parameter    ::  phifac = 8.,sfcfrac = 0.1
   real,parameter    ::  d1 = 0.02, d2 = 0.05, d3 = 0.001
   real,parameter    ::  h1 = 0.33333335, h2 = 0.6666667
   real,parameter    ::  ckz = 0.001,zfmin = 1.e-8,aphi5 = 5.,aphi16 = 16.
   real,parameter    ::  tmin=1.e-2
   real,parameter    ::  gamcrt = 3.,gamcrq = 2.e-3
   real,parameter    ::  xka = 2.4e-5
   integer,parameter ::  imvdif = 1


   real, dimension( its:ite, kts:kte, jts:jte )     ::                wscalek,wscalek2, brup_precalc
  !  real, dimension( ims:ime )              ::                           delta
   real, dimension( its:ite, kts:kte, jts:jte )     ::                     xkzml,xkzhl, &
                                                               zfacent,entfac
   real, dimension( its:ite, jts:jte )              ::                            ust3, &
                                                                       wstar3, &
                                                                     wstar3_2, &
                                                                  hgamu,hgamv, &
                                                                      wm2, we, &
                                                                       bfxpbl, &
                                                                hfxpbl,qfxpbl, &
                                                                ufxpbl,vfxpbl, &
                                                                        dthvx
   real    ::  prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,                           &
               dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,tvcon,conpr,     &
               prfac,prfac2,phim8z,radsum,tmp1,templ,rvls,temps,ent_eff,    &
               rcldb,bruptmp,radflux

#if defined(mpas)
!MPAS specific optional arguments for additional diagnostics:
   real,intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: kzhout,kzmout,kzqout
#endif

!$acc data create(hol,zq,thx,thvx,thlix,qtnp,del,dza,dzq,xkzo,za, pdh, dpdh, pdhi, &
!$acc hgamt,hgamq,entpbl,rhox,govrth,zl1,thermal,wscale,brdn,brup,phim,phih, &
!$acc dusfc,dvsfc,dtsfc,dqsfc,prpbl,wspd1,thermalli,xkzm,xkzh,f1,f2, brup_precalc, &
!$acc r1,r2,ad,au,cu,al,xkzq,zfac,rhox2,hgamt2,brcr,sflux,zol1,brcr_sbro, &
!$acc r3,f3,kpbl,kpblold,pblflg,sfcflg,stable,cloudflg,wscalek,wscalek2, &
#if defined(mpas)
!$acc kzhout,kzmout,kzqout, &
#endif
!$acc delta,wstar,xkzml,xkzhl,zfacent,entfac,ust3,wstar3,wstar3_2,uoce,voce, &
!$acc hgamu,hgamv,wm2, we,bfxpbl,hfxpbl,qfxpbl,ufxpbl,vfxpbl,dthvx, tke_pbl) &
!$acc present(u3d,v3d,t3d,qv3d,qc3d,qi3d,p3d,p3di,pi3d,                               &
!$acc         rublten,rvblten,rthblten,                                   &
!$acc         rqvblten,rqcblten,rqiblten,                         &
!$acc         dz8w,z8w,psfc,                                               &
!$acc         znu,znw,mut,p_top,                                        &
!$acc         znt,ust,hpbl,psim,psih,                   &
!$acc         hfx,qfx,wspd,br,                                       &
!$acc         kpbl2d,                                               &
!$acc         exch_h,                                                     &
!$acc         u10,v10,                                                     &
!$acc         rthraten, rho,                                            &
!$acc         regime, xland)

#if defined(mpas)
!$acc parallel loop gang vector collapse(3)
   do j = jts,jte
   do k = kts,kte
   do i = its,ite
      kzhout(i,k,j) = 0.
      kzmout(i,k,j) = 0.
      kzqout(i,k,j) = 0.
   enddo
   enddo
   enddo
!MPAS specific end.
#endif

  if(present(mut))then
!
! For ARW we will replace p and p8w with dry hydrostatic pressure
!
      !$acc parallel loop gang vector collapse(3) async(1)
      do j = jts,jte
        do k = kts,kte+1
          do i = its,ite
             if(k.le.kte)pdh(i,k,j) = mut(i,j)*znu(k) + p_top
             pdhi(i,k,j) = mut(i,j)*znw(k) + p_top
          enddo
        enddo
      enddo

  elseif(present(rho)) then
  ! 203 format(1x,i4,1x,i2,10(1x,e15.8))
      !For MPAS, we replace the hydrostatic pressures defined at theta and w points by
      !the dry hydrostatic pressures (Laura D. Fowler):
      !!!        k = kte+1

      !$acc parallel loop gang vector collapse(2) async(1)
      do j = jts,jte 
        do i = its,ite
           pdhi(i,kte+1,j) = p3di(i,kte+1,j)
        enddo
      enddo

      !$acc parallel loop gang vector private(rho_d) collapse(2) async(1)
      do j = jts,jte
        do i = its,ite
        !$acc loop seq
        do k = kte,kts,-1
           rho_d = rho(i,k,j) / (1. + qv3d(i,k,j))
           if(k.le.kte) pdhi(i,k,j) = pdhi(i,k+1,j) + g*rho_d*dz8w(i,k,j)
        enddo
        enddo
      enddo
      !$acc end parallel

        !$acc parallel async(1)
        !$acc loop gang vector collapse(3)
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           pdh(i,k,j) = 0.5*(pdhi(i,k,j) + pdhi(i,k+1,j))
        enddo
        enddo
        enddo
      !$acc end parallel

  else
        !$acc parallel loop gang vector collapse(3) async(1)
        do j = jts,jte
        do k = kts,kte+1
          do i = its,ite
            if(k.le.kte)pdh(i,k,j) = p3d(i,k,j)
            pdhi(i,k,j) = p3di(i,k,j)
          enddo
        enddo
        enddo

  endif
  !$acc wait(1)

!
!       call ysu2d(J=j,ux=u3d(ims,kms,j),vx=v3d(ims,kms,j)                       &
!               ,tx=t3d(ims,kms,j)                                               &
!               ,qx=qv2d(its,kts)                                                &
!               ,p2d=pdh(its,kts),p2di=pdhi(its,kts)                             &
!               ,pi2d=pi3d(ims,kms,j)                                            &
!               ,utnp=rublten(ims,kms,j),vtnp=rvblten(ims,kms,j)                 &
!               ,ttnp=rthblten(ims,kms,j),qtnp=qtnp(its,kts),ndiff=ndiff     &
!               ,cp=cp,g=g,rovcp=rovcp,rd=rd,rovg=rovg                           &    
!               ,xlv=xlv,rv=rv                                                   &
!               ,ep1=ep1,ep2=ep2,karman=karman                                   &
!               ,dz8w2d=dz8w(ims,kms,j)                                          &
!               ,psfcpa=psfc(ims,j),znt=znt(ims,j),ust=ust(ims,j)                &
!               ,hgamt=hgamt(ims,j),hgamq=hgamq(ims,j)                           &
!               ,entpbl=entpbl(ims,j),hpbl=hpbl(ims,j)                           &
!               ,regime=regime(ims,j),psim=psim(ims,j)                           &
!               ,psih=psih(ims,j),xland=xland(ims,j)                             &
!               ,hfx=hfx(ims,j),qfx=qfx(ims,j)                                   &
!               ,wspd=wspd(ims,j),br=br(ims,j)                                   &
!               ,dusfc=dusfc,dvsfc=dvsfc,dtsfc=dtsfc,dqsfc=dqsfc                 &
!               ,dt=dt,rcl=1.0,kpbl1d=kpbl2d(ims,j)                              &
!               ,exch_hx=exch_h(ims,kms,j)                                       &
!               ,wstar=wstar(ims,j)                                              &
!               ,delta=delta(ims,j)                                              &
!               ,u10=u10(ims,j),v10=v10(ims,j)                                   &
!               ,uox=uoce(ims,j),vox=voce(ims,j)                                 &
!               ,rthraten=rthraten(ims,kms,j),p2diORG=p3di(ims,kms,j)            &
!               ,tke_pbl=tke_pbl(ims,kms,j)                                      &
!               ,ysu_topdown_pblmix=ysu_topdown_pblmix                           &
!               ,ctopo=ctopo(ims,j),ctopo2=ctopo2(ims,j)                         &
! #if defined(mpas)
! !MPAS specific optional arguments for additional diagnostics:
!               ,kzh=kzhout(ims,kms,j)                                           &
!               ,kzm=kzmout(ims,kms,j)                                           &
!               ,kzq=kzqout(ims,kms,j)                                           &
! #endif
!               ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde               &
!               ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme               &
!               ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte   )
!
!
! !$acc end data
!
  !  end subroutine ysu_gpu
!
!-------------------------------------------------------------------------------
!
!    subroutine ysu2d(j,ux,vx,tx,qx,p2d,p2di,pi2d,                               &
!                   utnp,vtnp,ttnp,qtnp,ndiff,                                   &
!                   cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,                    &
!                   dz8w2d,psfcpa,                                               &
!                   znt,ust,hgamt,hgamq,entpbl,hpbl,psim,psih,                   &
!                   xland,hfx,qfx,wspd,br,                                       &
!                   dusfc,dvsfc,dtsfc,dqsfc,                                     &
!                   dt,rcl,kpbl1d,                                               &
!                   exch_hx,                                                     &
!                   wstar,delta,                                                 &
!                   u10,v10,                                                     &
!                   uox,vox,                                                     &
!                   rthraten,p2diORG,                                            &
!                   tke_pbl,                                                     &
!                   ysu_topdown_pblmix,                                          &
!                   ctopo,ctopo2,                                                &
!                   ids,ide, jds,jde, kds,kde,                                   &
!                   ims,ime, jms,jme, kms,kme,                                   &
!                   its,ite, jts,jte, kts,kte,                                   &
!                 !optional
!                   regime                                                       &
!                    )
! !-------------------------------------------------------------------------------
!    implicit none
!-------------------------------------------------------------------------------
!
!     this code is a revised vertical diffusion package ("ysupbl")
!     with a nonlocal turbulent mixing in the pbl after "mrfpbl".
!     the ysupbl (hong et al. 2006) is based on the study of noh
!     et al.(2003) and accumulated realism of the behavior of the
!     troen and mahrt (1986) concept implemented by hong and pan(1996).
!     the major ingredient of the ysupbl is the inclusion of an explicit
!     treatment of the entrainment processes at the entrainment layer.
!     this routine uses an implicit approach for vertical flux
!     divergence and does not require "miter" timesteps.
!     it includes vertical diffusion in the stable atmosphere
!     and moist vertical diffusion in clouds.
!
!     mrfpbl:
!     coded by song-you hong (ncep), implemented by jimy dudhia (ncar)
!              fall 1996
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     further modifications :
!              an enhanced stable layer mixing, april 2008
!              ==> increase pbl height when sfc is stable (hong 2010)
!              pressure-level diffusion, april 2009
!              ==> negligible differences
!              implicit forcing for momentum with clean up, july 2009
!              ==> prevents model blowup when sfc layer is too low
!              incresea of lamda, maximum (30, 0.1 x del z) feb 2010
!              ==> prevents model blowup when delz is extremely large
!              revised prandtl number at surface, peggy lemone, feb 2010
!              ==> increase kh, decrease mixing due to counter-gradient term
!              revised thermal, shin et al. mon. wea. rev. , songyou hong, aug 2011
!              ==> reduce the thermal strength when z1 < 0.1 h
!              revised prandtl number for free convection, dudhia, mar 2012
!              ==> pr0 = 1 + bke (=0.272) when newtral, kh is reduced
!              minimum kzo = 0.01, lo = min (30m,delz), hong, mar 2012
!              ==> weaker mixing when stable, and les resolution in vertical
!              gz1oz0 is removed, and phim phih are ln(z1/z0)-phim,h, hong, mar 2012
!              ==> consider thermal z0 when differs from mechanical z0
!              a bug fix in wscale computation in stable bl, sukanta basu, jun 2012
!              ==> wscale becomes small with height, and less mixing in stable bl
!              revision in background diffusion (kzo), jan 2016
!              ==> kzo = 0.1 for momentum and = 0.01 for mass to account for
!                  internal wave mixing of large et al. (1994), songyou hong, feb 2016
!              ==> alleviate superious excessive mixing when delz is large
!
!     references:
!
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev.
!        hong and pan (1996), mon. wea. rev.
!        noh, chun, hong, and raasch (2003), boundary layer met.
!        troen and mahrt (1986), boundary layer met.
!
!-------------------------------------------------------------------------------
!
!
!    integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,                 &
!                                     ims,ime, jms,jme, kms,kme,                 &
!                                     its,ite, jts,jte, kts,kte,                 &
!                                     j,ndiff

!    integer,  intent(in)      ::     ysu_topdown_pblmix 
! !
!    real,     intent(in   )   ::     dt,rcl,cp,g,rovcp,rovg,rd,xlv,rv
! !
!    real,     intent(in )     ::     ep1,ep2,karman
! !
!    real,     dimension( ims:ime, kms:kme ),                                    &
!              intent(in)      ::                                        pi2d
! !
!    real,     dimension( ims:ime, kms:kme )                                   , &
!              intent(in   )   ::                                            tx
!    real,     dimension( its:ite, kts:kte*ndiff )                             , &
!              intent(in   )   ::                                            qx
! !
!    real,     dimension( ims:ime, kms:kme )                                   , &
!              intent(inout)   ::                                          utnp, &
!                                                                          vtnp, &
!                                                                          ttnp
!    real,     dimension( its:ite, kts:kte, jts:jte, ndiff )                             , &
!              intent(inout)   ::                                          qtnp
! !
!    real,     dimension( its:ite, kts:kte+1 )                                 , &
!              intent(in   )   ::                                          p2di
! !
!    real,     dimension( its:ite, kts:kte )                                   , &
!              intent(in   )   ::                                           p2d

!    real,     dimension( its:ite, kts:kte )                                   , &
!              intent(inout)   ::                                       tke_pbl
! !
!    real,     dimension( ims:ime )                                            , &
!              intent(inout)   ::                                           ust, &
!                                                                         wstar, &
!                                                                         hgamt, &
!                                                                         hgamq, &
!                                                                        entpbl, &
!                                                                          hpbl, &
!                                                                           znt
!    real,     dimension( ims:ime )                                            , &
!              intent(in   )   ::                                         xland, &
!                                                                           hfx, &
!                                                                           qfx
! !
!    real,     dimension( ims:ime ), intent(inout)   ::                    wspd
!    real,     dimension( ims:ime ), intent(in  )    ::                      br
! !
!    real,     dimension( ims:ime ), intent(in   )   ::                    psim, &
!                                                                          psih
! !
! !
!    real,     dimension( ims:ime, kms:kme )                                   , &
!              intent(in   )   ::                                            ux, &
!                                                                            vx, &
!                                                                       rthraten
!    real,     dimension( ims:ime )                                            , &
!              optional                                                        , &
!              intent(in   )   ::                                         ctopo, &
!                                                                        ctopo2
!    real,     dimension( ims:ime )                                            , &
!              optional                                                        , &
!              intent(inout)   ::                                        regime
!
! local vars
!



!
!-------------------------------------------------------------------------------
!
! !$acc parallel num_gangs(1) num_workers(1)
   klpbl = kte
!
   cont=cp/g
   conq=xlv/g
   conw=1./g
   conwrc = conw*sqrt(rcl)
   conpr = bfac*karman*sfcfrac

   dtstep = dt
   dt2 = 2.*dtstep
   rdt = 1./dt2

!
! !$acc end parallel
!

!$acc parallel loop gang vector collapse(3) private(tvcon) async(10)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       thx(i,k,j) = t3d(i,k,j)/pi3d(i,k,j)
       thlix(i,k,j) = (t3d(i,k,j)-xlv*qc3d(i,k,j)/cp-2.834E6*qi3d(i,k,j)/cp)/pi3d(i,k,j)
!!!     enddo
!!!   enddo
!
!!!!$acc loop gang vector collapse(3) private(tvcon)
!!!   do k = kts,kte
!!!     do i = its,ite
!!!       tvcon = (1.+ep1*qv3d(i,k,j))
       thvx(i,k,j) = thx(i,k,j)*(1.+ep1*qv3d(i,k,j))
       tvcon = (1.+ep1*qv3d(i,k,j))
       rhox2(i,k,j) = pdh(i,k,j)/(rd*t3d(i,k,j)*tvcon)

     enddo
   enddo
   enddo

!$acc parallel loop gang vector collapse(2) async(10)
do j = jts,jte
   do i = its,ite
!!!     tvcon = (1.+ep1*qv3d(i,k,j))
!!!     rhox(i,j) = psfcpa(i,j)/(rd*t3d(i,1,j)*tvcon)
     rhox(i,j) = psfc(i,j)/(rd*t3d(i,1,j)*(1.+ep1*qv3d(i,1,j)))
     govrth(i,j) = g/thx(i,1,j)
!!!   enddo
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
!!!!$acc loop gang vector
!!!   do i = its,ite
     zq(i,1,j) = 0.
     uoce(i,j) = 0.
     voce(i,j) = 0.
   enddo
enddo

!$acc parallel loop gang vector collapse(3) wait(10) async(11)
  do j = jts,jte
   do k = kts+1,kte+1
     do i = its,ite
       zq(i,k,j) = z8w(i,k,j) - z8w(i,1,j)
     enddo
   enddo
   enddo

!$acc parallel async(11)
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       za(i,k,j) = 0.5*(zq(i,k,j)+zq(i,k+1,j))
       dzq(i,k,j) = zq(i,k+1,j)-zq(i,k,j)
       del(i,k,j) = pdhi(i,k,j)-pdhi(i,k+1,j)
       if (k<kte) dpdh(i,k,j) = pdh(i,k,j)-pdh(i,k+1,j)
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc parallel async(11)
!$acc loop gang vector collapse(2)
do j = jts,jte
   do i = its,ite
     dza(i,1,j) = za(i,1,j)
   enddo
  enddo
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts+1,kte
     do i = its,ite
       dza(i,k,j) = za(i,k,j)-za(i,k-1,j)
     enddo
   enddo
  enddo
!$acc end parallel
!
!-----initialize vertical tendencies and
!
!$acc parallel async(10)
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kms, kme
   do i = ims, ime
     rthblten(i,k,j) = 0.
     tke_pbl(i,k,j) = 0.
   enddo
   enddo
  enddo
!$acc loop gang vector collapse(4)
   do p = 1, ndiff
   do j = jts, jte
   do k = kts, kte
   do i = its, ite
     qtnp(i,k,j,p) = 0.
   enddo
   enddo
  enddo
  enddo
!$acc loop gang vector collapse(2)
do j = jts,jte
   do i = its,ite
     wspd1(i,j) = sqrt( (u3d(i,1,j)-uoce(i,j))*(u3d(i,1,j)-uoce(i,j)) + (v3d(i,1,j)-voce(i,j))*(v3d(i,1,j)-voce(i,j)) )+1.e-9
! !---- compute vertical diffusion
! !
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! !     compute preliminary variables
! !
     bfxpbl(i,j) = 0.0
     hfxpbl(i,j) = 0.0
     qfxpbl(i,j) = 0.0
     ufxpbl(i,j) = 0.0
     vfxpbl(i,j) = 0.0
     dthvx(i,j)  = 0.0
     hgamu(i,j)  = 0.0
     hgamv(i,j)  = 0.0
     wm2(i,j)    = 0.0
     we(i,j)     = 0.0
     delta(i,j)  = 0.0
     ust3(i,j)   = 0.0
     wstar3(i,j) = 0.0
     wstar3_2(i,j) = 0.0
     dusfc(i,j) = 0.
     dvsfc(i,j) = 0.
     dtsfc(i,j) = 0.
     dqsfc(i,j) = 0.

     hgamt(i,j)  = 0.
     hgamq(i,j)  = 0.
     wscale(i,j) = 0.
     kpbl(i,j)   = 1
     hpbl(i,j)   = 0.0
     zl1(i,j)    = dz8w(i,1,j)*0.5
     thermal(i,j)= thvx(i,1,j)
     thermalli(i,j) = thlix(i,1,j)
     pblflg(i,j) = .true.
     sfcflg(i,j) = .true.
     sflux(i,j) = hfx(i,j)/rhox(i,j)/cp + qfx(i,j)/rhox(i,j)*ep1*thx(i,1,j)
     if(br(i,j).gt.0.0) sfcflg(i,j) = .false.

     stable(i,j) = .false.
     brup(i,j) = br(i,j)
     brcr(i,j) = brcr_ub
   enddo
  enddo

!$acc loop gang vector collapse(3)
do j = jts,jte
    do k = kts,kte
    do i = its,ite
        xkzh(i,k,j)  = 0.0
        xkzm(i,k,j)  = 0.0
        xkzhl(i,k,j) = 0.0
        xkzml(i,k,j) = 0.0
    enddo
    enddo
    enddo

!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,klpbl
     do i = its,ite
       wscalek(i,k,j) = 0.0
       wscalek2(i,k,j) = 0.0
       zfac(i,k,j) = 0.0
     enddo
   enddo
  enddo

!$acc end parallel

!$acc parallel loop gang vector collapse(3) wait(11) async(12)
do j = jts,jte
   do k = kts,klpbl-1
     do i = its,ite
       xkzo(i,k,j) = ckz*dza(i,k+1,j)
     enddo
   enddo
enddo
!!
!     compute the first guess of pbl height
!
!

!$acc parallel loop gang vector collapse(3) wait(11) async(10)
do j = jts,jte
   do k = kts, kte
   do i = its, ite
     brup_precalc(i,k,j) = (thvx(i,k,j)-thermal(i,j))*(g*za(i,k,j)/thvx(i,1,j))/ &
            max(u3d(i,k,j)**2+v3d(i,k,j)**2,1.)

   enddo
   enddo
  enddo

!$acc parallel async(10)
!$acc loop gang vector tile(32,16)
  do j = jts,jte
   do i = its,ite
!$acc loop seq
     do k = 2,klpbl
       if(.not.stable(i,j))then
         brdn(i,j) = brup(i,j)
         brup(i,j) = brup_precalc(i,k,j)
         kpbl(i,j) = k
         stable(i,j) = brup(i,j).gt.brcr(i,j)
       else
          exit
       endif
     enddo
     k = kpbl(i,j)
     if(brdn(i,j).ge.brcr(i,j))then
       brint = 0.
     elseif(brup(i,j).le.brcr(i,j))then
       brint = 1.
     else
       brint = (brcr(i,j)-brdn(i,j))/(brup(i,j)-brdn(i,j))
     endif
     hpbl(i,j) = za(i,max(1,k-1),j)+brint*(za(i,k,j)-za(i,max(1,k-1),j))
     if(hpbl(i,j).lt.zq(i,2,j)) kpbl(i,j) = 1
     if(kpbl(i,j).le.1) pblflg(i,j) = .false.
     fm = psim(i,j)
     fh = psih(i,j)
     zol1(i,j) = max(br(i,j)*fm*fm/fh,rimin)
     if(sfcflg(i,j))then
       zol1(i,j) = min(zol1(i,j),-zfmin)
     else
       zol1(i,j) = max(zol1(i,j),zfmin)
     endif
     hol1 = zol1(i,j)*hpbl(i,j)/zl1(i,j)*sfcfrac
     if(sfcflg(i,j))then
       phim(i,j) = (1.-aphi16*hol1)**(-1./4.)
       phih(i,j) = (1.-aphi16*hol1)**(-1./2.)
       bfx0 = max(sflux(i,j),0.)
       hfx0 = max(hfx(i,j)/rhox(i,j)/cp,0.)
       qfx0 = max(ep1*thx(i,1,j)*qfx(i,j)/rhox(i,j),0.)
       wstar3(i,j) = (govrth(i,j)*bfx0*hpbl(i,j))
       wstar(i,j) = (wstar3(i,j))**h1
     else
       phim(i,j) = (1.+aphi5*hol1)
       phih(i,j) = phim(i,j)
       wstar(i,j)  = 0.
       wstar3(i,j) = 0.
     endif
     ust3(i,j)   = ust(i,j)**3.
     wscale(i,j) = (ust3(i,j)+phifac*karman*wstar3(i,j)*0.5)**h1
     wscale(i,j) = min(wscale(i,j),ust(i,j)*aphi16)
     wscale(i,j) = max(wscale(i,j),ust(i,j)/aphi5)
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
     if(sfcflg(i,j).and.sflux(i,j).gt.0.0)then
       gamfac   = bfac/rhox(i,j)/wscale(i,j)
       hgamt(i,j) = min(gamfac*hfx(i,j)/cp,gamcrt)
       hgamq(i,j) = min(gamfac*qfx(i,j),gamcrq)
       vpert = (hgamt(i,j)+ep1*thx(i,1,j)*hgamq(i,j))/bfac*afac
       thermal(i,j) = thermal(i,j)+max(vpert,0.)*min(za(i,1,j)/(sfcfrac*hpbl(i,j)),1.0)
       thermalli(i,j)= thermalli(i,j)+max(vpert,0.)*min(za(i,1,j)/(sfcfrac*hpbl(i,j)),1.0)
       hgamt(i,j) = max(hgamt(i,j),0.0)
       hgamq(i,j) = max(hgamq(i,j),0.0)
       brint    = -15.9*ust(i,j)*ust(i,j)/wspd(i,j)*wstar3(i,j)/(wscale(i,j)**4.)
       hgamu(i,j) = brint*u3d(i,1,j)
       hgamv(i,j) = brint*v3d(i,1,j)
     else
       pblflg(i,j) = .false.
     endif
   enddo
   enddo
!$acc end parallel

!$acc parallel loop gang vector collapse(3) async(10)
  do j = jts,jte
   do k = kts,kte
      do i = its,ite
        brup_precalc(i,k,j) = (thvx(i,k,j)-thermal(i,j))*(g*za(i,k,j)/thvx(i,1,j))/ &
            max(u3d(i,k,j)**2+v3d(i,k,j)**2,1.)
      enddo
    enddo
  enddo
!     enhance the pbl height by considering the thermal
!
!$acc parallel loop gang vector collapse(2) async(10)
  do j = jts,jte
   do i = its,ite
     if(pblflg(i,j))then
       kpbl(i,j) = 1
       hpbl(i,j) = zq(i,1,j)
       stable(i,j) = .false.
       brup(i,j) = br(i,j)
       brcr(i,j) = brcr_ub
     endif
   enddo
   enddo
!
!$acc parallel async(10)
!$acc loop gang vector tile(32,16)
  do j = jts,jte
   do i = its,ite
!$acc loop seq
     do k = 2,klpbl
       if(.not.stable(i,j).and.pblflg(i,j))then
         brdn(i,j) = brup(i,j)
         brup(i,j) = brup_precalc(i,k,j)
         kpbl(i,j) = k
         stable(i,j) = brup(i,j).gt.brcr(i,j)
       else
          exit
        endif
     enddo
   enddo
   enddo
!$acc end parallel
!
!     enhance pbl by theta-li
!
   if (ysu_topdown_pblmix.eq.1)then
!$acc parallel loop gang vector private(definebrup,spdk2,bruptmp) tile(32,16) async(10) 
    do j = jts,jte
     do i = its,ite
        kpblold(i,j) = kpbl(i,j)
        definebrup=.false.
!$acc loop seq
        do k = kpblold(i,j), kte-1
           spdk2   = max(u3d(i,k,j)**2+v3d(i,k,j)**2,1.)
           bruptmp = (thlix(i,k,j)-thermalli(i,j))*(g*za(i,k,j)/thlix(i,1,j))/spdk2
           stable(i,j) = bruptmp.ge.brcr(i,j)
           if (definebrup) then
            kpbl(i,j) = k
            brup(i,j) = bruptmp
            definebrup=.false.
           endif
           if (.not.stable(i,j)) then !overwrite brup brdn values
            brdn(i,j)=bruptmp
            definebrup=.true.
            pblflg(i,j)=.true.
           endif
        enddo
     enddo
     enddo
   endif

!$acc parallel async(10)
!$acc loop gang vector private(k,brint) tile(32,16)
  do j = jts,jte
   do i = its,ite
     if(pblflg(i,j)) then
       k = kpbl(i,j)
       if(brdn(i,j).ge.brcr(i,j))then
         brint = 0.
       elseif(brup(i,j).le.brcr(i,j))then
         brint = 1.
       else
         brint = (brcr(i,j)-brdn(i,j))/(brup(i,j)-brdn(i,j))
       endif
       hpbl(i,j) = za(i,max(1,k-1),j)+brint*(za(i,k,j)-za(i,max(1,k-1),j))
       if(hpbl(i,j).lt.zq(i,2,j)) kpbl(i,j) = 1
       if(kpbl(i,j).le.1) pblflg(i,j) = .false.
     endif
     if((.not.sfcflg(i,j)).and.hpbl(i,j).lt.zq(i,2,j)) then
       brup(i,j) = br(i,j)
       stable(i,j) = .false.
     else
       stable(i,j) = .true.
     endif
     if((.not.stable(i,j)).and.((xland(i,j)-1.5).ge.0))then
       wspd10 = u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j)
       wspd10 = sqrt(wspd10)
       ross = wspd10 / (cori*znt(i,j))
       brcr_sbro(i,j) = min(0.16*(1.e-7*ross)**(-0.18),.3)
     endif
     if(.not.stable(i,j))then
       if((xland(i,j)-1.5).ge.0)then
         brcr(i,j) = brcr_sbro(i,j)
       else
         brcr(i,j) = brcr_sb
       endif
     endif
!$acc loop seq
     do k = 2,klpbl
       if(.not.stable(i,j))then
         brdn(i,j) = brup(i,j)
         brup(i,j) = brup_precalc(i,k,j)
         kpbl(i,j) = k
         stable(i,j) = brup(i,j).gt.brcr(i,j)
       else
          exit
       endif
     enddo
     if((.not.sfcflg(i,j)).and.hpbl(i,j).lt.zq(i,2,j)) then
       k = kpbl(i,j)
       if(brdn(i,j).ge.brcr(i,j))then
         brint = 0.
       elseif(brup(i,j).le.brcr(i,j))then
         brint = 1.
       else
         brint = (brcr(i,j)-brdn(i,j))/(brup(i,j)-brdn(i,j))
       endif
       hpbl(i,j) = za(i,max(1,k-1),j)+brint*(za(i,k,j)-za(i,max(1,k-1),j))
       if(hpbl(i,j).lt.zq(i,2,j)) kpbl(i,j) = 1
       if(kpbl(i,j).le.1) pblflg(i,j) = .false.
     endif
!
!     estimate the entrainment parameters
!
     cloudflg(i,j)=.false. 
     if(pblflg(i,j)) then
       k = kpbl(i,j) - 1
       wm3       = wstar3(i,j) + 5. * ust3(i,j)
       wm2(i,j)    = wm3**h2
       bfxpbl(i,j) = -0.15*thvx(i,1,j)/g*wm3/hpbl(i,j)
       dthvx(i,j)  = max(thvx(i,k+1,j)-thvx(i,k,j),tmin)
       we(i,j) = max(bfxpbl(i,j)/dthvx(i,j),-sqrt(wm2(i,j)))
       if((qc3d(i,k,j)+qi3d(i,k,j)).gt.0.01e-3.and.ysu_topdown_pblmix.eq.1)then
           if ( kpbl(i,j) .ge. 2) then
                cloudflg(i,j)=.true. 
                templ=thlix(i,k,j)*(pdhi(i,k+1,j)/100000)**rovcp
                !rvls is ws at full level
                rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep2/pdhi(i,k+1,j))
                temps=templ + ((qv3d(i,k,j)+qc3d(i,k,j))-rvls)/(cp/xlv  + &
                ep2*xlv*rvls/(rd*templ**2))
                rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep2/pdhi(i,k+1,j))
                rcldb=max((qv3d(i,k,j)+qc3d(i,k,j))-rvls,0.)
                !entrainment efficiency
                dthvx(i,j)  = (thlix(i,k+2,j)+thx(i,k+2,j)*ep1*(qv3d(i,k+2,j)+qc3d(i,k+2,j))) &
                          - (thlix(i,k,j) + thx(i,k,j)  *ep1*(qv3d(i,k,j)  +qc3d(i,k,j)))
                dthvx(i,j)  = max(dthvx(i,j),0.1)
                tmp1      = xlv/cp * rcldb/(pi3d(i,k,j)*dthvx(i,j))
                ent_eff   = 0.2 * 8. * tmp1 +0.2

                radsum=0.
!$acc loop seq
                do kk = 1,kpbl(i,j)-1
                   radflux=rthraten(i,kk,j)*pi3d(i,kk,j) !converts theta/s to temp/s
                   radflux=radflux*cp/g*(p3di(i,kk,j)-p3di(i,kk+1,j)) ! converts temp/s to W/m^2
                   if (radflux < 0.0 ) radsum=abs(radflux)+radsum
                enddo
                radsum=max(radsum,0.0)

                !recompute entrainment from sfc thermals
                bfx0 = max(max(sflux(i,j),0.0)-radsum/rhox2(i,k,j)/cp,0.)
                bfx0 = max(sflux(i,j),0.0)
                wm3 = (govrth(i,j)*bfx0*hpbl(i,j))+5. * ust3(i,j)
                wm2(i,j)    = wm3**h2
                bfxpbl(i,j) = -0.15*thvx(i,1,j)/g*wm3/hpbl(i,j)
                dthvx(i,j)  = max(thvx(i,k+1,j)-thvx(i,k,j),tmin)
                we(i,j) = max(bfxpbl(i,j)/dthvx(i,j),-sqrt(wm2(i,j)))

                !entrainment from PBL top thermals
                bfx0 = max(radsum/rhox2(i,k,j)/cp-max(sflux(i,j),0.0),0.)
                bfx0 = max(radsum/rhox2(i,k,j)/cp,0.)
                wm3       = (g/thvx(i,k,j)*bfx0*hpbl(i,j)) ! this is wstar3(i,j)
                wm2(i,j)    = wm2(i,j)+wm3**h2
                bfxpbl(i,j) = - ent_eff * bfx0
                dthvx(i,j)  = max(thvx(i,k+1,j)-thvx(i,k,j),0.1)
                we(i,j) = we(i,j) + max(bfxpbl(i,j)/dthvx(i,j),-sqrt(wm3**h2))

                !wstar3_2
                bfx0 = max(radsum/rhox2(i,k,j)/cp,0.)
                wstar3_2(i,j) =  (g/thvx(i,k,j)*bfx0*hpbl(i,j))
                !recompute hgamt 
                wscale(i,j) = (ust3(i,j)+phifac*karman*(wstar3(i,j)+wstar3_2(i,j))*0.5)**h1
                wscale(i,j) = min(wscale(i,j),ust(i,j)*aphi16)
                wscale(i,j) = max(wscale(i,j),ust(i,j)/aphi5)
                gamfac   = bfac/rhox(i,j)/wscale(i,j)
                hgamt(i,j) = min(gamfac*hfx(i,j)/cp,gamcrt)
                hgamq(i,j) = min(gamfac*qfx(i,j),gamcrq)
                gamfac   = bfac/rhox2(i,k,j)/wscale(i,j)
                hgamt2(i,k,j) = min(gamfac*radsum/cp,gamcrt)
                hgamt(i,j) = max(hgamt(i,j),0.0) + max(hgamt2(i,k,j),0.0)
                brint    = -15.9*ust(i,j)*ust(i,j)/wspd(i,j)*(wstar3(i,j)+wstar3_2(i,j))/(wscale(i,j)**4.)
                hgamu(i,j) = brint*u3d(i,1,j)
                hgamv(i,j) = brint*v3d(i,1,j)
           endif
       endif
       prpbl(i,j) = 1.0
       dthx  = max(thx(i,k+1,j)-thx(i,k,j),tmin)
       dqx   = min(qv3d(i,k+1,j)-qv3d(i,k,j),0.0)
       hfxpbl(i,j) = we(i,j)*dthx
       qfxpbl(i,j) = we(i,j)*dqx
!
       dux = u3d(i,k+1,j)-u3d(i,k,j)
       dvx = v3d(i,k+1,j)-v3d(i,k,j)
       if(dux.gt.tmin) then
         ufxpbl(i,j) = max(prpbl(i,j)*we(i,j)*dux,-ust(i,j)*ust(i,j))
       elseif(dux.lt.-tmin) then
         ufxpbl(i,j) = min(prpbl(i,j)*we(i,j)*dux,ust(i,j)*ust(i,j))
       else
         ufxpbl(i,j) = 0.0
       endif
       if(dvx.gt.tmin) then
         vfxpbl(i,j) = max(prpbl(i,j)*we(i,j)*dvx,-ust(i,j)*ust(i,j))
       elseif(dvx.lt.-tmin) then
         vfxpbl(i,j) = min(prpbl(i,j)*we(i,j)*dvx,ust(i,j)*ust(i,j))
       else
         vfxpbl(i,j) = 0.0
       endif
       delb  = govrth(i,j)*d3*hpbl(i,j)
       delta(i,j) = min(d1*hpbl(i,j) + d2*wm2(i,j)/delb,100.)
     endif
   enddo
  enddo
!$acc end parallel
!
!$acc parallel wait(11) async(10)
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       if(pblflg(i,j).and.k.ge.kpbl(i,j))then
         entfac(i,k,j) = ((zq(i,k+1,j)-hpbl(i,j))/delta(i,j))**2.
       else
         entfac(i,k,j) = 1.e30
       endif
!
!     compute diffusion coefficients below pbl
!
       if(k.lt.kpbl(i,j)) then
         zfac(i,k,j) = min(max((1.-(zq(i,k+1,j)-zl1(i,j))/(hpbl(i,j)-zl1(i,j))),zfmin),1.)
         zfacent(i,k,j) = (1.-zfac(i,k,j))**3.
         wscalek(i,k,j) = (ust3(i,j)+phifac*karman*wstar3(i,j)*(1.-zfac(i,k,j)))**h1
         wscalek2(i,k,j) = (phifac*karman*wstar3_2(i,j)*(zfac(i,k,j)))**h1
         if(sfcflg(i,j)) then
           prfac = conpr
           prfac2 = 15.9*(wstar3(i,j)+wstar3_2(i,j))/ust3(i,j)/(1.+4.*karman*(wstar3(i,j)+wstar3_2(i,j))/ust3(i,j))
           prnumfac = -3.*(max(zq(i,k+1,j)-sfcfrac*hpbl(i,j),0.))**2./hpbl(i,j)**2.
         else
           prfac = 0.
           prfac2 = 0.
           prnumfac = 0.
           phim8z = 1.+aphi5*zol1(i,j)*zq(i,k+1,j)/zl1(i,j)
           wscalek(i,k,j) = ust(i,j)/phim8z
           wscalek(i,k,j) = max(wscalek(i,k,j),0.001)
         endif
         prnum0 = (phih(i,j)/phim(i,j)+prfac)
         prnum0 = max(min(prnum0,prmax),prmin)
           xkzm(i,k,j) = wscalek(i,k,j) *karman*    zq(i,k+1,j)      *    zfac(i,k,j)**pfac+ &
                       wscalek2(i,k,j)*karman*(hpbl(i,j)-zq(i,k+1,j))*(1-zfac(i,k,j))**pfac
         !Do not include xkzm at kpbl-1 since it changes entrainment
         if (k.eq.kpbl(i,j)-1.and.cloudflg(i,j).and.we(i,j).lt.0.0) then
           xkzm(i,k,j) = 0.0
         endif
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzq(i,k,j) = xkzm(i,k,j)/prnum*zfac(i,k,j)**(pfac_q-pfac)
         prnum0 = prnum0/(1.+prfac2*karman*sfcfrac)
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzh(i,k,j) = xkzm(i,k,j)/prnum
         xkzm(i,k,j) = max(xkzm(i,k,j),xkzo(i,k,j))
         xkzh(i,k,j) = max(xkzh(i,k,j),xkzo(i,k,j))
         xkzq(i,k,j) = max(xkzq(i,k,j),xkzo(i,k,j))
         xkzm(i,k,j) = min(xkzm(i,k,j),xkzmax)
         xkzh(i,k,j) = min(xkzh(i,k,j),xkzmax)
         xkzq(i,k,j) = min(xkzq(i,k,j),xkzmax)

! implicit tke
         c = zfac(i,k,j)**pfac*0.425 + (1. - zfac(i,k,j)**pfac)
         tke_pbl(i,k,j)=(xkzm(i,k,j)/max(1.,c*zq(i,k+1,j)*zfac(i,k,j)**pfac))**2

       endif
     enddo
   enddo
  enddo
!$acc end parallel
!
!     compute diffusion coefficients over pbl (free atmosphere)
!
!$acc parallel async(10)
!$acc loop gang vector collapse(3) private(ss,govrthv,ri,qmean,tmean,alph,chi,zk,rlamdz,rl2,dk,sri,prnum)
  do j = jts,jte
   do k = kts,kte-1
     do i = its,ite
       if(k.ge.kpbl(i,j)) then
         ss = ((u3d(i,k+1,j)-u3d(i,k,j))*(u3d(i,k+1,j)-u3d(i,k,j))                         &
              +(v3d(i,k+1,j)-v3d(i,k,j))*(v3d(i,k+1,j)-v3d(i,k,j)))                        &
              /(dza(i,k+1,j)*dza(i,k+1,j))+1.e-9
         govrthv = g/(0.5*(thvx(i,k+1,j)+thvx(i,k,j)))
         ri = govrthv*(thvx(i,k+1,j)-thvx(i,k,j))/(ss*dza(i,k+1,j))
         if(imvdif.eq.1.and.ndiff.ge.3)then
           if((qc3d(i,k,j)+qi3d(i,k,j)).gt.0.01e-3.and.(qc3d(i,k+1,j)+qi3d(i,k+1,j)).gt.0.01e-3)then
!      in cloud
             qmean = 0.5*(qv3d(i,k,j)+qv3d(i,k+1,j))
             tmean = 0.5*(t3d(i,k,j)+t3d(i,k+1,j))
             alph  = xlv*qmean/rd/tmean
             chi   = xlv*xlv*qmean/cp/rv/tmean/tmean
             ri    = (1.+alph)*(ri-g*g/ss/tmean/cp*((chi-alph)/(1.+chi)))
           endif
         endif
         zk = karman*zq(i,k+1,j)
         rlamdz = min(max(0.1*dza(i,k+1,j),rlam),300.)
         rlamdz = min(dza(i,k+1,j),rlamdz)
         rl2 = (zk*rlamdz/(rlamdz+zk))**2
         dk = rl2*sqrt(ss)
         if(ri.lt.0.)then
! unstable regime
           ri = max(ri, rimin)
           sri = sqrt(-ri)
           xkzm(i,k,j) = dk*(1+8.*(-ri)/(1+1.746*sri))
           xkzh(i,k,j) = dk*(1+8.*(-ri)/(1+1.286*sri))
         else
! stable regime
           xkzh(i,k,j) = dk/(1+5.*ri)**2
           prnum = 1.0+2.1*ri
           prnum = min(prnum,prmax)
           xkzm(i,k,j) = xkzh(i,k,j)*prnum
         endif
!
         xkzm(i,k,j) = max(xkzm(i,k,j),xkzo(i,k,j))
         xkzh(i,k,j) = max(xkzh(i,k,j),xkzo(i,k,j))
         xkzm(i,k,j) = min(xkzm(i,k,j),xkzmax)
         xkzh(i,k,j) = min(xkzh(i,k,j),xkzmax)
         xkzml(i,k,j) = xkzm(i,k,j)
         xkzhl(i,k,j) = xkzh(i,k,j)

! implicit tke
         tke_pbl(i,k,j)=(xkzm(i,k,j)/max(1.,(zk*rlamdz/(rlamdz+zk))))**2

        if(pblflg(i,j).and.entfac(i,k,j).lt.4.6) then
          xkzh(i,k,j) = -we(i,j)*dza(i,kpbl(i,j),j)*exp(-entfac(i,k,j))
          xkzh(i,k,j) = sqrt(xkzh(i,k,j)*xkzhl(i,k,j))
          xkzh(i,k,j) = max(xkzh(i,k,j),xkzo(i,k,j))
          xkzh(i,k,j) = min(xkzh(i,k,j),xkzmax)
          xkzq(i,k,j) = xkzh(i,k,j)
          xkzq(i,k,j) = sqrt(xkzq(i,k,j)*xkzhl(i,k,j))
          xkzq(i,k,j) = max(xkzq(i,k,j),xkzo(i,k,j))
          xkzq(i,k,j) = min(xkzq(i,k,j),xkzmax)
          xkzm(i,k,j) = prpbl(i,j)*xkzh(i,k,j)
          xkzm(i,k,j) = sqrt(xkzm(i,k,j)*xkzml(i,k,j))
          xkzm(i,k,j) = max(xkzm(i,k,j),xkzo(i,k,j))
          xkzm(i,k,j) = min(xkzm(i,k,j),xkzmax)
        else
          xkzq(i,k,j) = xkzh(i,k,j)
        endif

       endif
     enddo
   enddo
  enddo
!$acc end parallel

!
!     compute tridiagonal matrix elements for heat
!
!$acc parallel loop gang vector collapse(2) async(13)
do j = jts,jte
     do i = its,ite
       au(i,kte,j) = 0.
       al(i,kte,j) = 0.
     enddo
  enddo
!$acc parallel loop gang vector tile(32,16) async(13) wait(10)
  do j = jts,jte
   do i = its,ite
   ! implicit tke lower boundary condition
     tke_pbl(i,kts,j)=5.20*ust(i,j)*ust(i,j)*0.5
     entpbl(i,j) = max(0.,-1.0*we(i,j))

     f1(i,1,j) = thx(i,1,j)-300.+hfx(i,j)/cont/del(i,1,j)*dt2

    enddo
  enddo
!$acc parallel async(13)
!$acc loop gang vector collapse(3) 
  do j = jts,jte
  do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k,j)
       dtodsu = dt2/del(i,k+1,j)
       rdz    = 1./dza(i,k+1,j)
       tem1   = dpdh(i,k,j)*xkzh(i,k,j)*rdz
       dsdz2     = tem1*rdz
       exch_h(i,k+1,j) = xkzh(i,k,j)
       au(i,k,j)   = -dtodsd*dsdz2
       al(i,k,j)   = -dtodsu*dsdz2
       if (k==kts) then
        ad(i,k,j) = 1.0-au(i,k,j)
       else
        ad(i,k,j)   = 1.-(-dtodsd*dpdh(i,k-1,j)*xkzh(i,k-1,j)*(1./dza(i,k,j))**2)-au(i,k,j)
        if(k==kte-1) ad(i,k+1,j) = 1.-al(i,k,j)
       endif


      if(pblflg(i,j).and.k.lt.kpbl(i,j)) then
          dsdzt = tem1*(-hgamt(i,j)/hpbl(i,j)-hfxpbl(i,j)*zfacent(i,k,j)/xkzh(i,k,j))

          if (k==kpbl(i,j)-1) f1(i,k+1,j) = thx(i,k+1,j)-300.-dtodsu*dsdzt

          if (k==kts) then
            f1(i,k,j) = f1(i,k,j) + dtodsd*dsdzt
          else
            f1(i,k,j) = dtodsd*dsdzt

            rdz    = 1./dza(i,k,j)
            tem1   = dpdh(i,k-1,j)*xkzh(i,k-1,j)*rdz
            dsdzt = tem1*(-hgamt(i,j)/hpbl(i,j)-hfxpbl(i,j)*zfacent(i,k-1,j)/xkzh(i,k-1,j))

            f1(i,k,j)   = (thx(i,k,j)-300.-dtodsd*dsdzt) + f1(i,k,j)
          endif

      else
         f1(i,k+1,j) = thx(i,k+1,j)-300.
      endif
     enddo
   enddo
  enddo
!$acc end parallel
!
! copies here to avoid duplicate input args for tridin
!
!$acc parallel async(13)
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       cu(i,k,j) = au(i,k,j)
       r1(i,k,j) = f1(i,k,j)
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc wait(13)
   call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,jts,jte,kts,kte,1)
!
!     recover tendencies of heat
!
!$acc parallel async(16)
!$acc loop gang vector collapse(3)
  do j = jts,jte
    do k = kts,kte
     do i = its,ite
! !
! !     compute tridiagonal matrix elements for moisture, clouds, and gases
! !
       au(i,k,j) = 0.
       al(i,k,j) = 0.
       ad(i,k,j) = 0.
     enddo
   enddo
  enddo

!$acc end parallel

!$acc parallel loop gang vector tile(32,16) async(16)
  do j = jts,jte
   do i = its,ite
     ad(i,1,j) = 1.
     f3(i,1,j,1) = qv3d(i,1,j)+qfx(i,j)*g/del(i,1,j)*dt2
     f3(i,1,j,2) = qc3d(i,1,j)
     f3(i,1,j,3) = qi3d(i,1,j)

    enddo
  enddo
!$acc parallel async(16)
!$acc loop gang vector collapse(3) 
  do j = jts,jte
  do k = kts,kte-1
     do i = its,ite
       dtodsd = dt2/del(i,k,j)
       dtodsu = dt2/del(i,k+1,j)
       rdz    = 1./dza(i,k+1,j)
       tem1   = dpdh(i,k,j)*xkzq(i,k,j)*rdz
       dsdz2     = tem1*rdz
       au(i,k,j)   = -dtodsd*dsdz2
       al(i,k,j)   = -dtodsu*dsdz2
       if (k==kts) then
        ad(i,k,j) = 1.0-au(i,k,j)
       else
        ad(i,k,j)   = 1.-(-dtodsd*dpdh(i,k-1,j)*xkzq(i,k-1,j)*(1./dza(i,k,j))**2)-au(i,k,j)
        if(k==kte-1) ad(i,k+1,j) = 1.-al(i,k,j)
       endif


      if(pblflg(i,j).and.k.lt.kpbl(i,j)) then
          dsdzq = tem1*(-qfxpbl(i,j)*zfacent(i,k,j)/xkzq(i,k,j))

          if (k==kpbl(i,j)-1) f3(i,k+1,j,1) = qv3d(i,k+1,j)-dtodsu*dsdzq

          if (k==kts) then
            f3(i,k,j,1) = f3(i,k,j,1) + dtodsd*dsdzq
          else
            f3(i,k,j,1) = dtodsd*dsdzq

            rdz    = 1./dza(i,k,j)
            tem1   = dpdh(i,k-1,j)*xkzq(i,k-1,j)*rdz
            dsdzq = tem1*(-qfxpbl(i,j)*zfacent(i,k-1,j)/xkzq(i,k-1,j))

            f3(i,k,j,1)   = (qv3d(i,k,j)-dtodsd*dsdzq) + f3(i,k,j,1)
          endif

      else
         f3(i,k+1,j,1) = qv3d(i,k+1,j)
      endif
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc parallel async(16)
!$acc loop gang vector collapse(3)
    do j = jts,jte  
       do k = kts+1,kte
         do i = its,ite
           f3(i,k,j,2) = qc3d(i,k,j)
           f3(i,k,j,3) = qi3d(i,k,j)
         enddo
       enddo
      enddo

!
! copies here to avoid duplicate input args for tridin
!
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       cu(i,k,j) = au(i,k,j)
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc parallel async(16)
!$acc loop gang vector collapse(4)
   do ic = 1,ndiff
    do j = jts,jte
     do k = kts,kte
       do i = its,ite
         r3(i,k,j,ic) = f3(i,k,j,ic)
       enddo
     enddo
   enddo
  enddo
!$acc end parallel
!
!     solve tridiagonal problem for moisture, clouds, and gases
!
!$acc wait(16)
   call tridin_ysu(al,ad,cu,r3,au,f3,its,ite,jts,jte,kts,kte,ndiff)
!
!     recover tendencies of heat and moisture
!
!$acc parallel async(17)
!$acc loop gang vector private(qtend,ttend) collapse(3)
  do j = jts,jte
    do k = kts,kte
     do i = its,ite
! !$acc loop seq
       qtend = (f3(i,k,j,1)-qv3d(i,k,j))*rdt
       qtnp(i,k,j,1) = qtnp(i,k,j,1)+qtend
       qtnp(i,k,j,2) = qtnp(i,k,j,2)+((f3(i,k,j,2)-qc3d(i,k,j))*rdt)
       qtnp(i,k,j,3) = qtnp(i,k,j,3)+((f3(i,k,j,3)-qi3d(i,k,j))*rdt)

       ttend = (f1(i,k,j)-thx(i,k,j)+300.)*rdt*pi3d(i,k,j)
       rthblten(i,k,j) = rthblten(i,k,j)+ttend

      !  dqsfc(i,j) = dqsfc(i,j)+qtend*conq*del(i,k,j) ! DR: Commented out Sept27 2025, not used anywhere
                                                       ! and necesitates sequential loop
      !  dtsfc(i,j) = dtsfc(i,j)+ttend*cont*del(i,k,j)/pi3d(i,k,j)  ! DR: Commented out Sept27 2025, not used anywhere
                                                                    ! and necesitates sequential loop

     enddo
   enddo
  enddo
!$acc end parallel

if (present(rublten) .and. present(rvblten)) then
!
!     compute tridiagonal matrix elements for momentum
!
!$acc parallel loop gang vector collapse(3) async(17)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       au(i,k,j) = 0.
       al(i,k,j) = 0.
       ad(i,k,j) = 0.
       f1(i,k,j) = 0.
       f2(i,k,j) = 0.
     enddo
   enddo
  enddo

!$acc parallel loop gang vector tile(32,16) async(17)
  do j = jts,jte
   do i = its,ite
! paj: ctopo=1 if topo_wind=0 (default)
! mchen add this line to make sure NMM can still work with YSU PBL
      ad(i,1,j) = 1.+ust(i,j)**2/wspd1(i,j)*rhox(i,j)*g/del(i,1,j)*dt2                  &
      *(wspd1(i,j)/wspd(i,j))**2
     f1(i,1,j) = u3d(i,1,j)+uoce(i,j)*ust(i,j)**2*g/del(i,1,j)*dt2/wspd1(i,j)
     f2(i,1,j) = v3d(i,1,j)+voce(i,j)*ust(i,j)**2*g/del(i,1,j)*dt2/wspd1(i,j)
!$acc loop seq
     do k = kts,kte-1
       dtodsd = dt2/del(i,k,j)
       dtodsu = dt2/del(i,k+1,j)
       rdz    = 1./dza(i,k+1,j)
       tem1   = dpdh(i,k,j)*xkzm(i,k,j)*rdz
       if(pblflg(i,j).and.k.lt.kpbl(i,j))then
         dsdzu     = tem1*(-hgamu(i,j)/hpbl(i,j)-ufxpbl(i,j)*zfacent(i,k,j)/xkzm(i,k,j))
         dsdzv     = tem1*(-hgamv(i,j)/hpbl(i,j)-vfxpbl(i,j)*zfacent(i,k,j)/xkzm(i,k,j))
         f1(i,k,j)   = f1(i,k,j)+dtodsd*dsdzu
         f1(i,k+1,j) = u3d(i,k+1,j)-dtodsu*dsdzu
         f2(i,k,j)   = f2(i,k,j)+dtodsd*dsdzv
         f2(i,k+1,j) = v3d(i,k+1,j)-dtodsu*dsdzv
       else
         f1(i,k+1,j) = u3d(i,k+1,j)
         f2(i,k+1,j) = v3d(i,k+1,j)
       endif
       dsdz2     = tem1*rdz
       au(i,k,j)   = -dtodsd*dsdz2
       al(i,k,j)   = -dtodsu*dsdz2
       ad(i,k,j)   = ad(i,k,j)-au(i,k,j)
       ad(i,k+1,j) = 1.-al(i,k,j)
     enddo
   enddo
   enddo
!$acc end parallel
!
! copies here to avoid duplicate input args for tridi1n
!
!$acc parallel async(17)
!$acc loop gang vector collapse(3)
do j = jts,jte
   do k = kts,kte
     do i = its,ite
       cu(i,k,j) = au(i,k,j)
       r1(i,k,j) = f1(i,k,j)
       r2(i,k,j) = f2(i,k,j)
     enddo
   enddo
  enddo
!$acc end parallel
!$acc wait(17)
!
!     solve tridiagonal problem for momentum
!
   call tridi1n(al,ad,cu,r1,r2,au,f1,f2,its,ite,jts,jte,kts,kte)
!
!     recover tendencies of momentum

!$acc parallel loop gang vector collapse(3)
  do j = jts,jte
     do k = kts,kte
      do i = its,ite
! !$acc loop seq
       rublten(i,k,j) = (f1(i,k,j)-u3d(i,k,j))*rdt
       rvblten(i,k,j) = (f2(i,k,j)-v3d(i,k,j))*rdt
      !  dusfc(i,j) = dusfc(i,j) + rublten(i,k,j)*conwrc*del(i,k,j)! DR: Commented out Sept27 2025, not used anywhere
                                                          ! and necesitates sequential loop
      !  dvsfc(i,j) = dvsfc(i,j) + rvblten(i,k,j)*conwrc*del(i,k,j)! DR: Commented out Sept27 2025, not used anywhere
                                                          ! and necesitates sequential loop
     enddo
   enddo
   enddo
endif !end if momentum tendency present


!$acc wait(17)

!$acc parallel
!$acc loop gang vector collapse(3)
  do j = jts,jte
     do k = kts,kte
      do i = its,ite
         rthblten(i,k,j) = rthblten(i,k,j)/pi3d(i,k,j)
         rqvblten(i,k,j) = qtnp(i,k,j,1)
         rqcblten(i,k,j) = qtnp(i,k,j,2)
         rqiblten(i,k,j) = qtnp(i,k,j,3)
     enddo
   enddo
   enddo
!
!---- end of vertical diffusion
!
!$acc loop gang vector collapse(2)
  do j = jts,jte
   do i = its,ite
     kpbl2d(i,j) = kpbl(i,j)
   enddo
  enddo
!
! #if defined(mpas)
! !$acc loop gang vector collapse(3)
! do j = jts,jte
!       do k = kts,kte
!       do i = its,ite
!          kzhout(i,k,j) = xkzh(i,k,j)
!          kzmout(i,k,j) = xkzm(i,k,j)
!          kzqout(i,k,j) = xkzq(i,k,j)
!       enddo
!       enddo
!       enddo
! #endif
!

! !$acc loop gang vector collapse(3)
!     do j = jts,jte
!      do k = kts,kte
!        do i = its,ite
!          rthblten(i,k,j) = rthblten(i,k,j)/pi3d(i,k,j)
!          rqvblten(i,k,j) = qtnp(i,k,j,1)
!          rqcblten(i,k,j) = qtnp(i,k,j,2)
!          if(present(rqiblten)) rqiblten(i,k,j) = qtnp(i,k,j,3)
!        enddo
!      enddo
!      enddo
!$acc end parallel
!$acc end data
   end subroutine ysu_gpu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridi1n(cl,cm,cu,r1,r2,au,f1,f2,its,ite,jts,jte,kts,kte)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, jts, jte, kts,kte
!
   real, dimension( its:ite, kts+1:kte+1, jts:jte )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte, jts:jte )                                       , &
         intent(in   )  ::                                                 cm, &
                                                                           r1, &
                                                                           r2, &
                                                                           cu
                                                                           
  real, dimension( its:ite, kts:kte, jts:jte )                                       , &
         intent(inout   )  ::                                              f1, f2, au

   real    :: fk
   integer :: i,k,j,it
!
!-------------------------------------------------------------------------------
!!
!$acc data present(cl,cm,cu,r1,r2,au,f1,f2)
!$acc parallel async(100)
!$acc loop gang vector private(fk) collapse(2)
  do j = jts,jte
   do i = its,ite
     fk = 1./cm(i,1,j)
     au(i,1,j) = fk*cu(i,1,j)
     f1(i,1,j) = fk*r1(i,1,j)
     f2(i,1,j) = fk*r2(i,1,j)
   enddo
  enddo
!$acc end parallel
!
!$acc parallel wait(100) async(101)
!$acc loop gang vector private(fk) collapse(2)
  do j = jts,jte
   do i = its,ite
!$acc loop seq
     do k = kts+1,kte-1
       fk = 1./(cm(i,k,j)-cl(i,k,j)*au(i,k-1,j))
       au(i,k,j) = fk*cu(i,k,j)
       f1(i,k,j) = fk*(r1(i,k,j)-cl(i,k,j)*f1(i,k-1,j))
       f2(i,k,j) = fk*(r2(i,k,j)-cl(i,k,j)*f2(i,k-1,j))

     enddo
   enddo
  enddo
!$acc end parallel
!
!
!$acc parallel wait(101) async(102)
!$acc loop gang vector private(fk) collapse(2)
  do j = jts,jte
   do i = its,ite
     fk = 1./(cm(i,kte,j)-cl(i,kte,j)*au(i,kte-1,j))
     f1(i,kte,j) = fk*(r1(i,kte,j)-cl(i,kte,j)*f1(i,kte-1,j))
     f2(i,kte,j) = fk*(r2(i,kte,j)-cl(i,kte,j)*f2(i,kte-1,j))
   enddo
  enddo
!$acc end parallel
!
!
!$acc parallel wait(102)
!$acc loop gang vector collapse(2)
  do j = jts,jte
   do i = its,ite
!$acc loop seq
     do k = kte-1,kts,-1
       f1(i,k,j) = f1(i,k,j)-au(i,k,j)*f1(i,k+1,j)
       f2(i,k,j) = f2(i,k,j)-au(i,k,j)*f2(i,k+1,j)
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc end data
!
   end subroutine tridi1n
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,jts,jte,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, jts, jte, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1, jts:jte )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte, jts:jte )                                       , &
         intent(in   )  ::                                                 cm
   real, dimension( its:ite, kts:kte, jts:jte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte, jts:jte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu
   real, dimension( its:ite, kts:kte, jts:jte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,j,k,it
!
!-------------------------------------------------------------------------------
!
!
!$acc data present(cl,cm,cu,r2,au,f2)
!$acc parallel async(200)
!$acc loop gang vector collapse(3) private(fk)
   do it = 1,nt
     do j = jts,jte
       do i = its,ite
         fk = 1./cm(i,1,j)
         au(i,1,j) = fk*cu(i,1,j)
         f2(i,1,j,it) = fk*r2(i,1,j,it)
       enddo
     enddo
   enddo
!$acc end parallel
!
!$acc parallel async(200)
!$acc loop gang vector collapse(3) private(fk)
   do it = 1,nt
    do j = jts,jte
     do i = its,ite
!$acc loop seq
          do k = kts+1,kte-1
            fk = 1./(cm(i,k,j)-cl(i,k,j)*au(i,k-1,j))
            au(i,k,j) = fk*cu(i,k,j)
            f2(i,k,j,it) = fk*(r2(i,k,j,it)-cl(i,k,j)*f2(i,k-1,j,it))
       enddo
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc parallel async(200)
!$acc loop gang vector collapse(3) private(fk)
   do it = 1,nt
     do j = jts,jte
     do i = its,ite
       fk = 1./(cm(i,kte,j)-cl(i,kte,j)*au(i,kte-1,j))
       f2(i,kte,j,it) = fk*(r2(i,kte,j,it)-cl(i,kte,j)*f2(i,kte-1,j,it))
     enddo
   enddo
  enddo
!$acc end parallel
!
!$acc parallel async(200)
!$acc loop gang vector collapse(3)
   do it = 1,nt
     do j = jts,jte
     do i = its,ite
!$acc loop seq
         do k = kte-1,kts,-1
           f2(i,k,j,it) = f2(i,k,j,it)-au(i,k,j)*f2(i,k+1,j,it)
         enddo
       enddo
     enddo
   enddo
!$acc end parallel
!$acc end data
!$acc wait(200)
!
   end subroutine tridin_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine ysuinit_gpu(rublten,rvblten,rthblten,rqvblten,                       &
                      rqcblten,rqiblten,p_qi,p_first_scalar,                   &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   logical , intent(in)          :: restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real , dimension( ims:ime , kms:kme , jms:jme ), intent(out) ::             &
                                                                      rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   integer :: i, j, k, itf, jtf, ktf
!
   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)
!
   if(.not.restart)then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
            rublten(i,k,j) = 0.
            rvblten(i,k,j) = 0.
            rthblten(i,k,j) = 0.
            rqvblten(i,k,j) = 0.
            rqcblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   if (p_qi .ge. p_first_scalar .and. .not.restart) then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
           rqiblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   end subroutine ysuinit_gpu
!-------------------------------------------------------------------------------
end module module_bl_ysu_gpu
!-------------------------------------------------------------------------------
