import sys
import numpy as np
import xarray as xr

def merge_fields(dataset_path, dataset_out_path, base_names, num_fields):
    ds = xr.open_dataset(dataset_path)
    base_names = base_names.split(',')
    num_fields = int(num_fields)
    
    for base_name in base_names:
        merged_field = None
        for i in range(1, num_fields + 1):
            field_name = f"{base_name}{i}"
            if field_name in ds:
                field_data = ds[field_name]
                if merged_field is None:
                    merged_field = field_data.where(~np.isnan(field_data), other=0)
                else:
                    merged_field = field_data.where(~np.isnan(field_data), other=merged_field)
                ds = ds.drop_vars(field_name)
        
        if merged_field is not None:
            ds[base_name] = merged_field
    
    if "lat_b" in ds:
        ds = ds.drop_vars(["lat_b"])
    if "lon_b" in ds:
        ds = ds.drop_vars(["lon_b"])
    ds.to_netcdf(dataset_out_path)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python merge_fields.py <dataset_path> <dataset_out_path> <base_names> <num_fields>")
    else:
        merge_fields(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])