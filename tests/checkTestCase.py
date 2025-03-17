#!/usr/bin/env python3
import xarray as xr
import numpy as np
import sys
import os


def calculate_percent_diff(var1, var2):
    """Calculate the absolute percent difference between two arrays."""
    # Avoid division by zero by adding a small value where var2 is zero
    denominator = np.where(np.abs(var2) > 1e-10, var2, 1e-10)
    percent_diff = np.abs((var1 - var2) / denominator) * 100.0
    return percent_diff

def main():
    # Define file paths
    output_file = "output.nc"
    reference_file = "reference_data.nc"
    
    # Check if files exist
    if not os.path.exists(output_file):
        print(f"Error: {output_file} does not exist")
        return 1
    
    if not os.path.exists(reference_file):
        print(f"Error: {reference_file} does not exist")
        return 1
    
    # Open NetCDF files with xarray
    try:
        ds_output = xr.open_dataset(output_file)
        ds_reference = xr.open_dataset(reference_file)
    except Exception as e:
        print(f"Error opening NetCDF files: {e}")
        return 1
    
    # Variables to check
    variables = ['u', 'v', 'w', 'temperature', 'qv']
    
    # Flag to track if any variable exceeds threshold
    error_flag = False
    threshold = 5.0  # 5% threshold
    
    print("Comparing variables between output.nc and reference_data.nc:")
    print("-" * 60)
    
    # Check each variable
    for var_name in variables:
        try:
            if var_name not in ds_output or var_name not in ds_reference:
                print(f"Warning: Variable {var_name} not found in one or both files")
                continue
            
            var_output = ds_output[var_name].values
            var_reference = ds_reference[var_name].values
            
            # Calculate percent difference
            percent_diff = calculate_percent_diff(var_output, var_reference)
            max_diff = np.max(percent_diff)
            
            # Print results
            print(f"{var_name:12s}: Max absolute percent difference = {max_diff:.4f}%")
            
            # Check if difference exceeds threshold
            if max_diff > threshold:
                print(f"  ERROR: {var_name} exceeds the {threshold}% threshold")
                error_flag = True
                
        except Exception as e:
            print(f"Error comparing {var_name}: {e}")
            error_flag = True
    
    # Close files (xarray automatically closes files when datasets go out of scope,
    # but it's good practice to close them explicitly)
    ds_output.close()
    ds_reference.close()
    
    print("-" * 60)
    if error_flag:
        print("Test FAILED: One or more variables exceed the allowed difference threshold!")
        return 1
    else:
        print("Test PASSED: All variables within allowed difference threshold.")
        return 0

if __name__ == "__main__":
    sys.exit(main())