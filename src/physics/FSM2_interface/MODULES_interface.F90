!-----------------------------------------------------------------------
! output to HICAR variables which are visible to the HICAR model
!-----------------------------------------------------------------------
module MODULES_interface

real, allocatable :: &
  Esrf_(:,:),       &! Moisture flux from the surface (kg/m^2/s)
  Gsoil_(:,:),      &! Heat flux into soil (W/m^2)
  H_(:,:),           &! Sensible heat flux to the atmosphere (W/m^2)
  LE_(:,:),          &! Latent heat flux to the atmosphere (W/m^2)
  Melt_(:,:),       &! Surface melt rate (kg/m^2/s)
  Rnet_(:,:),       &! Net radiation (W/m^2)
  Roff_(:,:),       &! Total runoff (kg/m^2)
  snowdepth_(:,:),  &! Snow depth (m)
  SWE_(:,:),        & ! Snow water equivalent (kg/m^2)
  KH_(:,:),         &! Eddy diffusivity for heat to the atmosphere (m/s)  
  meltflux_out_(:,:), &! Runoff from snowmelt at base of snow (kg/m^2)
  Sliq_out_(:,:),   & ! Total LWC (kg/m^2)
  dSWE_salt_(:,:),    &! SWE change due to saltation (kg/m^2)
  dSWE_susp_(:,:),    &! SWE change due to suspension (kg/m^2)
  dSWE_subl_(:,:),    &! SWE change due to sublimation (kg/m^2)
  dSWE_slide_(:,:)     ! SWE change due to snow slides (kg/m^2)

integer :: &
  firstit           ! Flag for first iteration

end module MODULES_interface
