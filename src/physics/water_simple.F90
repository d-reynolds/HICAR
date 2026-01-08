!>----------------------------------------------------------
!!  Simple open water flux calculations
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_water_simple
    use data_structures
    use options_interface,   only : options_t
    use mod_wrf_constants,   only : XLV, gravity, KARMAN
    use mod_atm_utilities,   only : sat_mr
    implicit none

contains

    subroutine water_simple(options, sst, psfc, wind, qv, temperature,  &
                            sensible_heat, latent_heat, landmask, &
                            qv_surf, evap_flux, tskin, coef_heat_exch, vegtype, ims, ime, kms, kme, jms, jme)
        implicit none
        type(options_t),intent(in)    :: options        
        real,    dimension(ims:ime,jms:jme),  intent(inout) :: sensible_heat, latent_heat, qv_surf, evap_flux, tskin
        real,    dimension(ims:ime,jms:jme),  intent(in)    :: sst, psfc, wind, coef_heat_exch
        real,    dimension(ims:ime,jms:jme),  intent(in)    :: landmask
        real,    dimension(ims:ime,kms:kme,jms:jme),  intent(in)    :: qv, temperature
        integer, dimension(ims:ime,jms:jme),  optional, intent(in)    :: vegtype
        integer, intent(in)                     :: ims, ime, kms, kme, jms, jme

        integer :: i, j, options_water, options_water_cat, options_lake_cat

        options_water = options%physics%watersurface
        options_water_cat = options%lsm%water_category
        options_lake_cat  = options%lsm%lake_category

        !$acc data present(sst,psfc,wind,qv,temperature,sensible_heat,latent_heat,landmask,qv_surf,evap_flux,tskin,coef_heat_exch,vegtype) copyin(ims,ime,jms,jme)
        !$acc parallel loop gang vector collapse(2)
        do j=jms,jme
            do i=ims,ime
                if(                                                                         &
                    ( (options_water==kWATER_SIMPLE) .AND.                   &   ! If lakemodel is not selected, use this
                      (landmask(i,j)==kLC_WATER)                                            & !(n.b. in case noah (mp or lsm) is not used, landmask may not be set correctly!)
                    )                                                                       &
                    .OR.                                                                    &
                    ( (options_water==kWATER_LAKE) .AND.                     &    ! if lake model is selected, and
                      (vegtype(i,j).eq.options_water_cat) .AND.            &    !   gridcell is water,
                      (vegtype(i,j).ne.options_lake_cat)                   &    !   but not lake (i.e ocean)
                    )                                                                         &
                )then
                    qv_surf(i,j) = 0.98 * sat_mr(sst(i,j),psfc(i,j)) ! multiply by 0.98 to account for salinity
                    qv_surf(i,j) = max(min(qv_surf(i,j),qv(i,j)), 1.0e-6)

                    sensible_heat(i,j) = coef_heat_exch(i,j) * wind(i,j) * (sst(i,j)-temperature(i,1,j))
                    evap_flux(i,j)     = coef_heat_exch(i,j) * wind(i,j) * (qv_surf(i,j)-qv(i,1,j))
                    
                    latent_heat(i,j)   = evap_flux(i,j) * XLV
                    tskin(i,j)   = sst(i,j)

                endif
            end do
        end do
        !$acc end data

    end subroutine water_simple

end module module_water_simple
