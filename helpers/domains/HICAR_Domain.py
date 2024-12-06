# Load modules
import os
import numpy as np
import xarray as xr
import pyproj
import horayzon as hray

import warnings
import os
import sys
sys.path.append(os.getcwd())
import ProjHelpers as ph

# from pysheds.grid import Grid
# import skimage.filters
# import skimage.morphology
# from scipy.spatial import KDTree
# import rasterio as rio
# import rioxarray

# GLOBAL definitions
WGS64_crs = pyproj.CRS('EPSG:4326')

def wholeShebang(ds_in,ds_in_rad=0,res=50,terr_filter=10,TPI_thresh=100,valley_thresh=1500,LL_border=0.05,LU_Category='USGS'):
    
    if ('lat' in ds_in_rad):
        print('Adding Radiation Shading variable...')
        ds_in = addHorayzonParms(ds_in,ds_in_rad)
    #ds_in = addRidgeValleyDists(ds_in,ds_in_rad,res=res,terr_filter=terr_filter,TPI_thresh=TPI_thresh,\
    #                            valley_thresh=valley_thresh,LL_border=LL_border)
    
    ds_in = addSlopeAspect(ds_in,res)
    ds_in = addSnowHoldingDepth(ds_in)

    if (not('landmask' in ds_in.variables)): 
        ds_in['landmask'] = ds_in['landuse'].copy()
        WATER_LANDUSE = -1
        if (LU_Category=='USGS'): WATER_LANDUSE = 16
        if (WATER_LANDUSE > 0):
            ds_in['landmask'].values = np.where(ds_in.landuse.values==WATER_LANDUSE,2,1)

    
    return ds_in

def addHorayzonParms(ds_in,ds_in_rad):
    # Description: Compute gridded topographic parameters (slope angle and aspect,
    #              horizon and sky view factor) for an input DEM Ignore Earth's surface
    #              curvature.
    #
    # Copyright (c) 2022 ETH Zurich, Christian R. Steger
    # Edited        2023 SLF, Dylan Reynolds
    # MIT License

    
    # -----------------------------------------------------------------------------
    # Settings
    # -----------------------------------------------------------------------------

    # Domain size and computation settings
    dist_search = 20.0  # search distance for horizon [kilometre]
    azim_num = 90  # number of azimuth sampling directions [-]
    ellps = "WGS84"  # Earth's surface approximation (sphere, GRS80 or WGS84)

    domain = {"lon_min": np.amin(ds_in.lon.values), "lon_max": np.amax(ds_in.lon.values),
              "lat_min": np.amin(ds_in.lat.values), "lat_max": np.amax(ds_in.lat.values)}
    file_hori_temp = "./horayzon_temp_out.nc"

    # -----------------------------------------------------------------------------
    # Compute and save topographic parameters
    # -----------------------------------------------------------------------------

    # Load required DEM data (including outer boundary zone)
    lon = ds_in_rad.lon.values.astype('float64')[0,:]
    lat = np.flipud(ds_in_rad.lat.values.astype('float64'))[:,0]
    elevation = np.flipud(ds_in_rad.topo.values.astype('float32'))

    # -> GeoTIFF can also be read with GDAL if available (-> faster)

    # Compute ellipsoidal heights
    elevation += hray.geoid.undulation(lon, lat, geoid="EGM96")  # [m]

    # Compute indices of inner domain
    slice_in = (slice(np.where(lat >= domain["lat_max"])[0][-1],
                      np.where(lat <= domain["lat_min"])[0][0] + 1),
                slice(np.where(lon <= domain["lon_min"])[0][-1],
                      np.where(lon >= domain["lon_max"])[0][0] + 1))
    offset_0 = slice_in[0].start
    offset_1 = slice_in[1].start

    # Compute ECEF coordinates
    x_ecef, y_ecef, z_ecef = hray.transform.lonlat2ecef(*np.meshgrid(lon, lat),
                                                        elevation, ellps=ellps)
    dem_dim_0, dem_dim_1 = elevation.shape

    # Compute ENU coordinates
    trans_ecef2enu = hray.transform.TransformerEcef2enu(
        lon_or=lon[int(len(lon) / 2)], lat_or=lat[int(len(lat) / 2)], ellps=ellps)
    x_enu, y_enu, z_enu = hray.transform.ecef2enu(x_ecef, y_ecef, z_ecef,
                                                  trans_ecef2enu)

    # Compute unit vectors (up and north) in ENU coordinates for inner domain
    vec_norm_ecef = hray.direction.surf_norm(*np.meshgrid(lon[slice_in[1]],
                                                          lat[slice_in[0]]))
    vec_north_ecef = hray.direction.north_dir(x_ecef[slice_in], y_ecef[slice_in],
                                              z_ecef[slice_in], vec_norm_ecef,
                                              ellps=ellps)
    del x_ecef, y_ecef, z_ecef
    vec_norm_enu = hray.transform.ecef2enu_vector(vec_norm_ecef, trans_ecef2enu)
    vec_north_enu = hray.transform.ecef2enu_vector(vec_north_ecef, trans_ecef2enu)
    del vec_norm_ecef, vec_north_ecef

    # Merge vertex coordinates and pad geometry buffer
    vert_grid = hray.auxiliary.rearrange_pad_buffer(x_enu, y_enu, z_enu)

    # Compute horizon
    hori, azim = hray.horizon.horizon_gridded(vert_grid, dem_dim_0, dem_dim_1,
                                 vec_norm_enu, vec_north_enu,
                                 offset_0, offset_1,
                                 #x_axis_val=lon[slice_in[1]].astype(np.float32),
                                 #y_axis_val=lat[slice_in[0]].astype(np.float32),
                                 #x_axis_name="lon", y_axis_name="lat",
                                 dist_search=dist_search,
                                 azim_num=azim_num)


    # Compute rotation matrix (global ENU -> local ENU)
    rot_mat_glob2loc = hray.transform.rotation_matrix_glob2loc(vec_north_enu,
                                                               vec_norm_enu)
    del vec_north_enu, vec_norm_enu

    # Compute slope
    slice_in_a1 = (slice(slice_in[0].start - 1, slice_in[0].stop + 1),
                   slice(slice_in[1].start - 1, slice_in[1].stop + 1))
    vec_tilt = hray.topo_param.slope_plane_meth(x_enu[slice_in_a1],
                                                y_enu[slice_in_a1],
                                                z_enu[slice_in_a1],
                                                rot_mat=rot_mat_glob2loc,
                                                output_rot=True)[1:-1, 1:-1]
    del rot_mat_glob2loc
    del x_enu, y_enu, z_enu

    # Load horizon data
    lon_2d, lat_2d = np.meshgrid(lon[slice_in[1]],lat[slice_in[0]])

    # Compute Sky View Factor
    svf = hray.topo_param.sky_view_factor(azim, hori, vec_tilt)
    rad_ds = xr.Dataset(data_vars=dict(
            hlm=(["azim", "y", "x"], 90.0 - np.maximum(0.0,np.rad2deg(np.transpose(hori,[2,0,1])))),
            svf=(["y", "x"], svf),
        ),
        coords=dict(
            lon=(["y", "x"], lon_2d),
            lat=(["y", "x"], lat_2d),
            azim=(["azim"],np.maximum(0.0,np.rad2deg(azim)))))
    
    rad_ds.hlm.attrs['units'] = "degrees"
    rad_ds.azim.attrs['units'] = "degrees"
    rad_ds.lon.attrs['units'] = "degrees"
    rad_ds.lat.attrs['units'] = "degrees"

    rad_ds.lat.values = np.flipud(rad_ds.lat.values)
    rad_ds.svf.values = np.flipud(rad_ds.svf.values)
    rad_ds.hlm.values = np.flip(rad_ds.hlm.values,axis=1)

    hlm_ds, __ = ph.regridToDs(rad_ds,ds_in,WGS64_crs,ds_in_xvar='lon',ds_in_yvar='lat',\
                        ds_in_xdim='x',ds_in_ydim='y',to_ds_xvar='lon',to_ds_yvar='lat',\
                        to_ds_xdim='x',to_ds_ydim='y',method='bilinear')
    
    ds_in["hlm"]=(("a", "y", "x"), hlm_ds.hlm.values)
    ds_in["svf"]=(("y", "x"), hlm_ds.svf.values)
    
    return ds_in


def makeGrid(src_crs,resolution,x_min,x_max,y_min,y_max):
    if (src_crs == WGS64_crs):
        LoLat = 111165.92953253252
        
        mean_lat = (y_min+y_max)/2.0
        LoLat, LoLon = get_length_of_degree(mean_lat)

        n_lat = int((abs(y_max-y_min)*LoLat)/resolution)
        n_lon = int((abs(x_max-x_min)*LoLon)/resolution)

        lats = np.linspace(y_max,y_min,n_lat)
        lons = np.linspace(x_min,x_max,n_lon)
    
        lons_2d, lats_2d = np.meshgrid(lons,lats)

    else:
        tformer = pyproj.Transformer.from_crs(src_crs, WGS64_crs,always_xy=True)

        n_y = int((abs(y_max-y_min))/resolution)
        n_x = int((abs(x_max-x_min))/resolution)

        ys = np.linspace(y_max,y_min,n_y)
        xs = np.linspace(x_min,x_max,n_x)
    
        x_grid, y_grid = np.meshgrid(xs,ys)

        xs, ys = tformer.transform(x_grid.flat, y_grid.flat) 

        #out output will be as 1D arrays, so reshape back to the input x-y shape
        lons_2d = np.reshape(xs,(n_y,n_x))
        lats_2d = np.reshape(ys,(n_y,n_x))

    to_ds = xr.Dataset({"lat": (["y", "x"], np.flipud(lats_2d)),\
                        "lon": (["y", "x"], np.flipud(lons_2d))})      
        
    return to_ds

def addSlopeAspect(ds_in,res=50):
    # we easily calculate the sf/dx, df/dy by giving the resolution as below
    slpx=np.gradient(ds_in.topo.data,res,res)[1]
    slpy=np.gradient(ds_in.topo.data,res,res)[0]

    # first, we calculate everything in dgree then conert it to radian
    slope=np.hypot(slpy,slpx)
    slope_angle = np.arctan(slope)*180/np.pi              # slp is the slope angle in degree, i.e. flat is slp = 0; a vertical face is slp = 90
    aspect_angle = np.arctan2(-slpx,-slpy)*180/np.pi
    aspect_angle[aspect_angle<0] = 360 + aspect_angle[aspect_angle<0]   # asp is the aspect angle in degree, so that 0 deg is North, 90 deg is East, 180 deg is South, and 270 deg is East                                      

    #
    slope_angle = slope_angle*np.pi/180    # slp is the slope angle in degree, i.e. flat is slp = 0; a vertical face is slp = 90
    aspect_angle = aspect_angle*np.pi/180


    # adding slope, slope_angle and aspect_angle
    ds_in["slope"]=(("y", "x"), slope)
    ds_in["slope_rad"]=(("y", "x"), slope_angle)
    ds_in["aspect_rad"]=(("y", "x"), aspect_angle)
    
    return ds_in

def addSnowHoldingDepth(ds_in):
    #The following code comes from the SLF OSHD pre-processing code, written by Louis Quéno 2022
    #Shd (Snow holding depth) is needed by the Snowslide parameterization (Bernhardt and Schulz 2010)
    #Values for the parameterization of Shd come from (Marsh) et al., 2022
    
    slope_thresh = np.maximum(np.rad2deg(ds_in["slope_rad"].values),10.0)
    
    shd_norm = 3178.4 * slope_thresh**(-1.998)
    cos_slope_thresh = np.cos(np.deg2rad(slope_thresh))
    cos_slope_thresh = np.maximum(cos_slope_thresh,0.001)
    Shd = shd_norm * cos_slope_thresh
    
    # adding slope, slope_angle and aspect_angle
    ds_in["Shd"]=(("y", "x"), Shd)
    
    return ds_in

def get_length_of_degree(latitude):
    #Courtesy of CHAT-GPT!!!
    # Constants for the WGS-84 ellipsoid
    a = 6378137.0
    b = 6356752.3142
    e_squared = 0.006694379990197

    # Convert latitude to radians
    lat_radians = np.deg2rad(latitude)

    # Calculate the length of a degree of latitude
    numerator = a*(1-e_squared)
    denominator = np.power(1-e_squared*np.power(np.sin(lat_radians),2), 1.5)
    lat_length = (np.pi/180)*numerator/denominator

    # Calculate the length of a degree of longitude
    lon_length = (np.pi/180)*a*np.cos(lat_radians)/np.sqrt(1-e_squared*np.power(np.sin(lat_radians),2))

    return lat_length, lon_length


# this function takes in a xarray dataset ds with coordinates lat and lon, and a larger
# xarray dataset large_dom with coordinates lat and lon. It also takes in the resolution of
# the large_dom dataset and the ds dataset, as well as the window size in km.
# it first finds the inndices in large_dom for the lower left and upper right corners of the ds dataset
# it then extracts the data from large_dom corresponding to a window 15 km larger than the ds dataset
# it returns the extracted data in a new xarray dataset

def extract_rad_window(ds, large_dom, large_dom_res=50, window_km=15):

    # get the coordinates of the ds dataset
    lat_min = ds.lat.min().values
    lat_max = ds.lat.max().values
    lon_min = ds.lon.min().values
    lon_max = ds.lon.max().values

    # get the coordinates of the swiss_dom dataset
    lat_min_large = large_dom.lat.min().values
    lat_max_large = large_dom.lat.max().values
    lon_min_large = large_dom.lon.min().values
    lon_max_large = large_dom.lon.max().values

    mean_lat = (lat_min+lat_max)/2.0
    LoLat, LoLon = get_length_of_degree(mean_lat)

    # calculate the indices of the lower left and upper right corners of the ds dataset
    distances = np.sqrt((large_dom.lat - lat_min)**2 + 
                        (large_dom.lon - lon_min)**2)
    
    # Step 3: Find the index of the minimum distance
    closest_index = distances.argmin()
    
    # Step 4: Convert the 1D index to 2D coordinates
    j_min, i_min = np.unravel_index(closest_index, distances.shape)   
    
    distances = np.sqrt((large_dom.lat - lat_max)**2 + 
                        (large_dom.lon - lon_max)**2)
    
    # Step 3: Find the index of the minimum distance
    closest_index = distances.argmin()
    
    # Step 4: Convert the 1D index to 2D coordinates
    j_max, i_max = np.unravel_index(closest_index, distances.shape)  

    # calculate the indices of the lower left and upper right corners of the window
    i_min_window = i_min - int(window_km*1000 / large_dom_res)
    i_max_window = i_max + int(window_km*1000 / large_dom_res)
    j_min_window = j_min - int(window_km*1000 / large_dom_res)
    j_max_window = j_max + int(window_km*1000 / large_dom_res)

    # extract the data from the swiss_dom dataset
    dom_window = large_dom.isel(x=slice(i_min_window, i_max_window), y=slice(j_min_window, j_max_window))

    return dom_window

# This function was once used to compute the ridge and valley distances for a thermal wind parameterization
# It is not used in the current version of the HICAR model, but the code may be useful to future users
# so it is kept in here. Uncomment and awake from the dead at your owk risk
# def addRidgeValleyDists(ds_in,ds_rad,res=250,terr_filter=10,TPI_thresh=200,valley_thresh=1500,LL_border=0.05):
    
#     #Regridding usually triggers some annoying warnings, so silence these
#     warnings.filterwarnings("ignore")
    
#     x_min = np.amin(ds_in.lon.values)-LL_border
#     x_max = np.amax(ds_in.lon.values)+LL_border
#     y_min = np.amin(ds_in.lat.values)-LL_border
#     y_max = np.amax(ds_in.lat.values)+LL_border
    
#     #Length of degree latitude/longitude
#     LoLat, LoLon = get_length_of_degree((y_min+y_max)*0.5)
    
#     x_ind_offset = int(LL_border*LoLon/res)-5
#     y_ind_offset = int(LL_border*LoLat/res)-5
    
#     print('Cutting dataset for ridge+valley computations')
    
#     to_ds = makeGrid(WGS64_crs,resolution=res,x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max)

#     ds_rad_cut = ds_rad.where((ds_rad.lon < np.amax(to_ds.lon.values)) & \
#                               (ds_rad.lon > np.amin(to_ds.lon.values)) & \
#                               (ds_rad.lat < np.amax(to_ds.lat.values)) & \
#                               (ds_rad.lat > np.amin(to_ds.lat.values)),drop=True)


#     large_ds = ph.regridToDs(ds_rad_cut,to_ds,WGS64_crs,ds_in_xvar='lon',ds_in_yvar='lat',\
#                         ds_in_xdim='x',ds_in_ydim='y',to_ds_xvar='lon',to_ds_yvar='lat',\
#                         to_ds_xdim='x',to_ds_ydim='y',method='bilinear')[0]

    
#     dem_large = large_ds['topo'].values
    
#     #First compute the ridge mask by finding exposed elements
#     TPI = dem_large - skimage.filters.gaussian(dem_large, sigma=terr_filter)
#     ridges = (TPI > TPI_thresh)

#     ridges_skel = skimage.morphology.skeletonize(ridges)

#     #Skeletonization of ridge can shift the mask off of the ridge line, so now correct the mask
#     # Get the coordinates of the skel lines
#     y_points, x_points = np.nonzero(ridges_skel)
#     skel_line_coords = np.column_stack((x_points, y_points))

#     # Get the coordinates of the ridges lines
#     y_lines, x_lines = np.nonzero(ridges)
#     line_coords = np.column_stack((x_lines, y_lines))

#     # Build a KDTree from the point coordinates
#     skel_tree = KDTree(skel_line_coords)

#     # Find the nearest point to each line coordinate
#     distances, indices = skel_tree.query(line_coords)

#     for i, (x, y) in enumerate(line_coords):
#         if (ridges[y,x] and not(ridges_skel[y,x])):
#                 # If our elevation is higher than closest point
#                 if (dem_large[y,x] > dem_large[y_points[indices[i]], x_points[indices[i]]]):
#                     # Add ourselves to res_skel
#                     ridges_skel[y,x] = 1.0
#                     # Delete closest point from res_skel
#                     ridges_skel[y_points[indices[i]], x_points[indices[i]]] = np.nan

                    
#     # Now compute valley mask using flow lines and concentrations: 
#     large_ds.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
#     large_ds.rio.write_crs(4326, inplace=True)
    
#     DEM_tiff_fn = "temp.tif"

#     large_ds['topo'].rio.to_raster(DEM_tiff_fn)
#     grid = Grid.from_raster(DEM_tiff_fn)
#     dem = grid.read_raster(DEM_tiff_fn)

#     # (taken from example code on pyshed documentation)
#     # Condition DEM
#     # ----------------------
#     # Fill pits in DEM
#     pit_filled_dem = grid.fill_pits(dem)

#     # Fill depressions in DEM
#     flooded_dem = grid.fill_depressions(pit_filled_dem)

#     # Resolve flats in DEM
#     dem_fixed = grid.resolve_flats(flooded_dem)
                    
#     #Compute flow direction and accumulation 
#     # Determine D8 flow directions from DEM
#     # ----------------------
#     # Specify directional mapping
#     dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

#     # Compute flow directions
#     # -------------------------------------
#     fdir = grid.flowdir(dem_fixed, dirmap=dirmap)

#     # Calculate flow accumulation
#     # --------------------------
#     acc = grid.accumulation(fdir, dirmap=dirmap)

#     valleys = acc > valley_thresh
#     valleys_skel = valleys #skimage.morphology.skeletonize(valleys)
        
        
#     print('Done masking ridges and valleys')
#     i_indx,j_indx = np.indices(ridges_skel.shape, sparse=True)

#     r_dists = np.zeros(ridges_skel.shape)
#     r_drop  = np.zeros(ridges_skel.shape)
#     v_dists = np.zeros(valleys_skel.shape)
#     v_drop  = np.zeros(valleys_skel.shape)

#     max_dist = res*len(i_indx)*len(j_indx)
#     print('Computing ridge/valley distances...')
#     for i in range(y_ind_offset,ridges_skel.shape[0]-y_ind_offset):
#         if (((i-y_ind_offset)%10)==0): 
#             print('In row '+str(i-y_ind_offset)+' of '+str(ridges_skel.shape[0]-y_ind_offset*2))
#         for j in range(x_ind_offset,ridges_skel.shape[1]-x_ind_offset):
#             if (ridges_skel[i,j] < 1 or valleys_skel[i,j] < 1):
#                 temp_dists = np.sqrt((i_indx-i)**2 + (j_indx-j)**2)*res
#                 is_higher = dem_large[i,j] < dem_large
#             # Only do this for points not on a ridge
#             if (ridges_skel[i,j] < 1):
#                 # Mask distances using ridge-mask and all cells higher than current cell
#                 possible_dists = np.where((is_higher==True) & (ridges_skel==True),temp_dists,max_dist)
#                 # Get location of the minimum value to find height drop
#                 min_ind = np.argmin(possible_dists)
#                 # Distance returned is minimum of masked
#                 r_dists[i,j] = possible_dists.flat[min_ind]
#                 r_drop[i,j] = abs(dem_large[i,j]-dem_large.flat[min_ind])

#             # Only do this for points not in the valley bottom
#             if (valleys_skel[i,j] < 1):
#                 # Mask distances using valley-mask and all cells lower than current cell
#                 possible_dists = np.where((is_higher==False) & (valleys_skel==True),temp_dists,max_dist)
#                 # Get location of the minimum value to find height drop
#                 min_ind = np.argmin(possible_dists)
#                 # Distance returned is minimum of masked
#                 v_dists[i,j] = possible_dists.flat[min_ind]
#                 v_drop[i,j] = abs(dem_fixed[i,j]-dem_fixed.flat[min_ind])

#     valley_mask = xr.DataArray(valleys_skel, coords=large_ds.coords, dims=large_ds.topo.dims, name='valley_mask')
#     valley_dists = xr.DataArray(v_dists, coords=large_ds.coords, dims=large_ds.topo.dims, name='valley_dists')
#     valley_drop = xr.DataArray(v_drop, coords=large_ds.coords, dims=large_ds.topo.dims, name='valley_drop')

#     ridge_mask = xr.DataArray(ridges_skel, coords=large_ds.coords, dims=large_ds.topo.dims, name='ridge_mask')
#     ridge_dists = xr.DataArray(r_dists, coords=large_ds.coords, dims=large_ds.topo.dims, name='ridge_dists')
#     ridge_drop = xr.DataArray(r_drop, coords=large_ds.coords, dims=large_ds.topo.dims, name='ridge_drop')
#     TPI = xr.DataArray(TPI, coords=large_ds.coords, dims=large_ds.topo.dims, name='TPI')

#     large_ds['TPI'] = TPI
#     large_ds['valley_mask'] = valley_mask
#     large_ds['valley_dists'] = valley_dists
#     large_ds['valley_drop'] = valley_drop

#     large_ds['ridge_mask'] = ridge_mask
#     large_ds['ridge_dists'] = ridge_dists
#     large_ds['ridge_drop'] = ridge_drop
    
#     cut_ds = ph.regridToDs(large_ds,ds_in,pyproj.CRS.from_epsg(4326),ds_in_xvar='lon',ds_in_yvar='lat',\
#                         ds_in_xdim='x',ds_in_ydim='y',to_ds_xvar='lon',to_ds_yvar='lat',\
#                         to_ds_xdim='x',to_ds_ydim='y',method='bilinear')[0]

#     cut_ds.valley_mask.values = cut_ds.valley_mask.values > 0
#     cut_ds.ridge_mask.values  = cut_ds.ridge_mask.values > 0

#     ds_in['TPI'] = cut_ds.TPI
#     ds_in['valley_mask'] = cut_ds.valley_mask
#     ds_in['valley_dists'] = cut_ds.valley_dists
#     ds_in['valley_drop'] = cut_ds.valley_drop

#     ds_in['ridge_mask'] = cut_ds.ridge_mask
#     ds_in['ridge_dists'] = cut_ds.ridge_dists
#     ds_in['ridge_drop'] = cut_ds.ridge_drop
    
#     warnings.filterwarnings("ignore")
    
#     return ds_in

