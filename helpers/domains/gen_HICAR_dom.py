import xarray as xr
import os
import sys
sys.path.append(os.getcwd())
import HICAR_Domain as hd

###########################################################################
#######################  User options, to be edited  ######################
###########################################################################

#The resolution of the domain
res = 250

# The target domain, including lat and lon variables named as "lat" and "lon", and
# a DEM labeled as "topo". Optionally, landuse and landmask variables should be specified here.
target_domain_fn = '/capstor/store/cscs/userlab/s1308/dreynold/HIMA/domains/Tajiki_250m_filled_nested.nc'
# A domain with extent ~20km beyond the borders of the above target domain.
# Only lat,lon, and topo are required variables.
large_domain_fn = '/capstor/store/cscs/userlab/s1308/dreynold/HIMA/domains/Tajiki_1km.nc'
# Name of output file
output_domain_fn = '/capstor/store/cscs/userlab/s1308/dreynold/HIMA/domains/Tajiki_250m_filled_nested_rad.nc'

topo_var = 'HGT_M'
lat_var = 'XLAT_M'
lon_var = 'XLONG_M'

# classification system for land use categories. Used to create land mask
# based on what the water type is for the land use classification. Currently
# only USGS is supported
LU_Category = 'USGS'
# These are used in the calculation of ridelines, and can be tuned if the user
# is not satisfied with the deliniation of ridgelines in the output file
# terr_filter = 10
# TPI_thresh = 100

###########################################################################
############################  End user options  ###########################
###########################################################################

dom = xr.open_dataset(target_domain_fn)
if (not(large_domain_fn=="")):
    dom_rad = xr.open_dataset(large_domain_fn)
else:
    dom_rad = 0

# For all 2D variables dom, rename the first dimension to "x" and the second to "y"
for var in dom.data_vars:
    if len(dom[var].dims) == 2:
        # get the first dimension name
        dim1 = dom[var].dims[0]
        dim2 = dom[var].dims[1]
        #dom[var] = dom[var].rename({dim1:'y',dim2:'x'})
for var in dom.coords:
    if len(dom[var].dims) == 2:
        # get the first dimension name
        dim1 = dom[var].dims[0]
        dim2 = dom[var].dims[1]
        #dom[var] = dom[var].rename({dim1:'y',dim2:'x'})

dom_out = hd.wholeShebang(dom,dom_rad,res=res,LU_Category=LU_Category,topo_var=topo_var,lat_var=lat_var,lon_var=lon_var)

dom_out.to_netcdf(output_domain_fn)
