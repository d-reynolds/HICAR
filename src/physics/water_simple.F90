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
                            qv_surf, evap_flux, tskin, coef_heat_exch, vegtype, its, ite, kts, kte, jts, jte)
        implicit none
        type(options_t),intent(in)    :: options        
        real,    dimension(its:ite,jts:jte),  intent(inout) :: sensible_heat, latent_heat, qv_surf, evap_flux, tskin
        real,    dimension(its:ite,jts:jte),  intent(in)    :: sst, psfc, wind, qv, temperature, coef_heat_exch
        integer, dimension(its:ite,jts:jte),  intent(in)    :: landmask
        integer, dimension(its:ite,jts:jte),  optional, intent(in)    :: vegtype
        integer, intent(in)                     :: its, ite, kts, kte, jts, jte

        integer :: i, j


        do j=jts,jte
            do i=its,ite
                if(                                                                         &
                    ( (options%physics%watersurface==kWATER_SIMPLE) .AND.                   &   ! If lakemodel is not selected, use this
                      (landmask(i,j)==kLC_WATER)                                            & !(n.b. in case noah (mp or lsm) is not used, landmask may not be set correctly!)
                    )                                                                       &
                    .OR.                                                                    &
                    ( (options%physics%watersurface==kWATER_LAKE) .AND.                     &    ! if lake model is selected, and
                      (vegtype(i,j).eq.options%lsm%water_category) .AND.            &    !   gridcell is water,
                      (vegtype(i,j).ne.options%lsm%lake_category)                   &    !   but not lake (i.e ocean)
                    )                                                                         &
                )then
                    qv_surf(i,j) = 0.98 * sat_mr(sst(i,j),psfc(i,j)) ! multiply by 0.98 to account for salinity
                    qv_surf(i,j) = max(min(qv_surf(i,j),qv(i,j)), 1.0e-6)
                    sensible_heat(i,j) = coef_heat_exch(i,j) * wind(i,j) * (sst(i,j)-temperature(i,j))
                    evap_flux(i,j)     = coef_heat_exch(i,j) * wind(i,j) * (qv_surf(i,j)-qv(i,j))
                    latent_heat(i,j)   = evap_flux(i,j) * XLV
                    tskin(i,j)   = sst(i,j)

                endif
            end do
        end do

    end subroutine water_simple

end module module_water_simple
