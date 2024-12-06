import xarray as xr
import numpy as np
import sys

# Function to stagger the values of an array along the x-axis
def stagger_x(arr, dir):
    arr_staggered = 0.5 * (arr[:,:-1] + arr[:,1:])
    if dir == 'forward':
        back_end = arr[:,-1] + (arr[:,-1] - arr[:,-2])*0.5
        arr_staggered = xr.concat([arr_staggered, back_end], dim='x')
    elif dir == 'backward':
        front_end = arr[:,0] + (arr[:,0] - arr[:,1])*0.5
        arr_staggered = xr.concat([front_end, arr_staggered], dim='x')
    return arr_staggered

# Function to stagger the values of an array along the y-axis
def stagger_y(arr, dir):
    arr_staggered = 0.5 * (arr[:-1,:] + arr[1:,:])
    if dir == 'forward':
        back_end = arr[-1,:] + (arr[-1,:] - arr[-2,:])*0.5
        arr_staggered = xr.concat([arr_staggered, back_end.expand_dims("y", axis=0)], dim='y')
    elif dir == 'backward':
        front_end = arr[0,:] + (arr[0,:] - arr[1,:])*0.5
        arr_staggered = xr.concat([front_end.expand_dims("y", axis=0), arr_staggered], dim='y')
    return arr_staggered

def enforceESMFConformance(src_ds, dst_ds):
    #For this to be true, they must contain a lat and lon variable
    if 'lat' not in src_ds.variables:
        raise ValueError("The src topographic dataset does not contain a lat variable")
    if 'lon' not in src_ds.variables:
        raise ValueError("The src topographic dataset does not contain a lon variable")
    if 'lat' not in dst_ds.variables:
        raise ValueError("The dst dataset does not contain a lat variable")
    if 'lon' not in dst_ds.variables:
        raise ValueError("The dst dataset does not contain a lon variable")

    # Next, the lat and lon variables must contain a units attribute
    if 'units' not in src_ds['lat'].attrs:
        # Set the units attribute to be 'degrees_north'
        src_ds['lat'].attrs['units'] = 'degrees_north'
    if 'units' not in src_ds['lon'].attrs:
        # Set the units attribute to be 'degrees_east'
        src_ds['lon'].attrs['units'] = 'degrees_east'
    if 'units' not in dst_ds['lat'].attrs:
        # Set the units attribute to be 'degrees_north'
        dst_ds['lat'].attrs['units'] = 'degrees_north'
    if 'units' not in dst_ds['lon'].attrs:
        # Set the units attribute to be 'degrees_east'
        dst_ds['lon'].attrs['units'] = 'degrees_east'

    # Next, the lat and lon variables must be coordinates for the datset
    if 'lat' not in src_ds.coords:
        src_ds = src_ds.set_coords('lat')
    if 'lon' not in src_ds.coords:
        src_ds = src_ds.set_coords('lon')
    if 'lat' not in dst_ds.coords:
        dst_ds = dst_ds.set_coords('lat')
    if 'lon' not in dst_ds.coords:
        dst_ds = dst_ds.set_coords('lon')

    return src_ds, dst_ds

# Open the datasets
src_ds = xr.open_mfdataset(sys.argv[1])
dst_ds = xr.open_mfdataset(sys.argv[2])

print('Testing that datasets are ESMF_Regrid conformant...\n')
src_ds, dst_ds = enforceESMFConformance(src_ds, dst_ds)

# Extract lat and lon values
src_lat = src_ds['lat'].values
src_lon = src_ds['lon'].values
dst_lat = dst_ds['lat'].values
dst_lon = dst_ds['lon'].values

# Check if dst_fn dataset fits within src_fn dataset
if ((min(dst_lat.max(),src_lat.max()) > max(dst_lat.min(),src_lat.min()))  and
    (min(dst_lon.max(),src_lon.max()) > max(dst_lon.min(),src_lon.min()))):

    print('Cropping src topographic dataset...\n')

    # Find min and max values of lat and lon in dst_fn dataset
    min_lat = dst_lat.min() - 0.05
    max_lat = dst_lat.max() + 0.05
    min_lon = dst_lon.min() - 0.05
    max_lon = dst_lon.max() + 0.05
    
    # Use `.where` with `.compute()` to handle Dask-backed variables efficiently
    lat_mask = (src_ds["lat"] >= min_lat) & (src_ds["lat"] <= max_lat)
    lon_mask = (src_ds["lon"] >= min_lon) & (src_ds["lon"] <= max_lon)

    # Combine masks and apply them
    mask = lat_mask & lon_mask
    cropped_ds = src_ds.where(mask.compute(), drop=True)

    # Create new variables, lat_b and lon_b, to store the corner coordinates for each
    # grid cell of lat/lon on the cropped_ds, and on the dst_ds. lat_b and lon_b are
    # x_b by y_b arrays, where the last dimension stores the corner coordinates.
    # lat_b and lon_b are used to create the grid for regridding.

    # having x and y dimensions defined with values can screw up the concatenation here
    if 'x' in cropped_ds.coords and 'y' in cropped_ds.coords:
        cropped_ds = cropped_ds.drop_vars(['x','y'])

    lat_b_1 = stagger_x(stagger_y(cropped_ds['lat'],dir='backward'),dir='backward')
    lat_b_2 = stagger_x(stagger_y(cropped_ds['lat'],dir='backward'),dir='forward')
    lat_b_3 = stagger_x(stagger_y(cropped_ds['lat'],dir='forward'),dir='forward')
    lat_b_4 = stagger_x(stagger_y(cropped_ds['lat'],dir='forward'),dir='backward')

    lon_b_1 = stagger_x(stagger_y(cropped_ds['lon'],dir='backward'),dir='backward')
    lon_b_2 = stagger_x(stagger_y(cropped_ds['lon'],dir='backward'),dir='forward')
    lon_b_3 = stagger_x(stagger_y(cropped_ds['lon'],dir='forward'),dir='forward')
    lon_b_4 = stagger_x(stagger_y(cropped_ds['lon'],dir='forward'),dir='backward')

    #Add lat_b and lon_b to the cropped_ds
    lat_b = xr.concat([lat_b_1,lat_b_2,lat_b_3,lat_b_4],dim='corner')
    lat_b = lat_b.transpose('y', 'x', 'corner')
    lon_b = xr.concat([lon_b_1,lon_b_2,lon_b_3,lon_b_4],dim='corner')
    lon_b = lon_b.transpose('y', 'x', 'corner')

    #Drop coordinates form lat/lon_b, since they themselves will be coordinates
    lon_b = lon_b.drop_vars(['lat','lon'])
    lat_b = lat_b.drop_vars(['lat','lon'])
    cropped_ds = cropped_ds.assign_coords(lat_b=lat_b, lon_b=lon_b)
    cropped_ds['lon'].attrs['bounds'] = 'lon_b'
    cropped_ds['lat'].attrs['bounds'] = 'lat_b'

    # having x and y dimensions defined with values can screw up the concatenation here
    if 'x' in dst_ds.coords and 'y' in dst_ds.coords:
        dst_ds = dst_ds.drop_vars(['x','y'])

    lat_b_1 = stagger_x(stagger_y(dst_ds['lat'],dir='backward'),dir='backward')
    lat_b_2 = stagger_x(stagger_y(dst_ds['lat'],dir='backward'),dir='forward')
    lat_b_3 = stagger_x(stagger_y(dst_ds['lat'],dir='forward'),dir='forward')
    lat_b_4 = stagger_x(stagger_y(dst_ds['lat'],dir='forward'),dir='backward')

    lon_b_1 = stagger_x(stagger_y(dst_ds['lon'],dir='backward'),dir='backward')
    lon_b_2 = stagger_x(stagger_y(dst_ds['lon'],dir='backward'),dir='forward')
    lon_b_3 = stagger_x(stagger_y(dst_ds['lon'],dir='forward'),dir='forward')
    lon_b_4 = stagger_x(stagger_y(dst_ds['lon'],dir='forward'),dir='backward')

    #Add lat_b and lon_b to the cropped_ds
    lat_b = xr.concat([lat_b_1,lat_b_2,lat_b_3,lat_b_4],dim='corner')
    lat_b = lat_b.transpose('y', 'x', 'corner')
    lon_b = xr.concat([lon_b_1,lon_b_2,lon_b_3,lon_b_4],dim='corner')
    lon_b = lon_b.transpose('y', 'x', 'corner')
    #Drop coordinates form lat/lon_b, since they themselves will be coordinates
    lon_b = lon_b.drop_vars(['lat','lon'])
    lat_b = lat_b.drop_vars(['lat','lon'])
    dst_ds = dst_ds.assign_coords(lat_b=lat_b, lon_b=lon_b)
    dst_ds['lon'].attrs['bounds'] = 'lon_b'
    dst_ds['lat'].attrs['bounds'] = 'lat_b'

    # Get the size of the cropped dataset in memory
    size = cropped_ds.lat.nbytes

    # Convert to a number of files to save, using the cutoff size
    # which was passed in from the Regrid_script
    cutoff = int(sys.argv[3])
    num_files = int(np.ceil(size / cutoff))

    if (num_files > 1):
        print('Cropped src topographic dataset is too large to regrid in one file')
    print('Saving into:' + str(num_files) + ' files')
    # Split the cropped dataset into a series of smaller datasets
    max_y = cropped_ds.lat.shape[0]
    increment = (6.0/5)*(max_y/(num_files))
    y_start = 0
    y_end = int(increment)
    # Save the cropped dataset to a series of NetCDF files
    for i in range(1,int(num_files+1)):
        cropped_ds.isel(y=slice(y_start,y_end)).to_netcdf('.tmpsrc_mosaic'+str(i)+'.nc')

        #Include a little overlap so that regridded variables requirig surrounding data are OK
        y_start = int(y_end - (increment/5))
        y_end = int(min(y_start + increment, max_y))


    # Save the filtered dataset to a new NetCDF file
    dst_ds.to_netcdf('.tmpdst1234.nc')

    print('Cropping completed successfully')
else:
    raise ValueError("The dst dataset does not fit within the src topographic dataset")