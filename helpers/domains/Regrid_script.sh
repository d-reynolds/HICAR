#
# Use as: 
# ./Regrid_script.sh dst_file.nc var1,var2 100
#
# Note that ESMF_Regrid has a limit to grids with a number of coordinates
# less than 2^31. If you get an error, try to reduce the number of
# coordinates in the grid by using a coarser grid, or by using a smaller
# grid (i.e. chopping up the grid into smaller pieces, regridding, then
# stitching the pieces back together). This is a limitation of the ESMF
# library, where 32-bit integers are used to store array indices. Future
# versions of ESMF may address this issue.
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# If debugging, remove the "--no_log" flag below
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##################################################################
##################################################################
# Arguments:
# $1: dst_file.nc -- destination file, must have a 2D lat and lon grid
# $2: var1,var2   -- variables to regrid, separated by commas. Variables must
#                    be present in the master files, defined below 
# $3: nprocs      -- number of processors to use for regridding
##################################################################
##################################################################
# Arbitrary cutoff for the size of the lat or lon array in the src file,
# expressed in bytes. If the size of the lat or lon array is larger than
# this value, the src file will be cut into smaller pieces for regridding.
# This value may need to be adjusted depending on the machine memory available.
kFILE_CUTOFF=1000000000
# The relative path to the master file for landuse data
TOPO_MASTER_FN="../MasterDomains/SwissMaster.nc"
# The relative path to the master file for topographic data
LU_MASTER_FN="../MasterDomains/CORINEMaster.nc"
##################################################################
##################################################################


#Delete tmp files, if they still exist from a prior failed run
rm -f .tmpsrc_mosaic*.nc .tmpdst1234.nc

echo "Attempting to crop the src topographic dataset for a faster regrid"

#Cut src domain for faster regridding
result=$(python ds_cutter.py "$TOPO_MASTER_FN" "$1" $kFILE_CUTOFF)

# Extract the value after "into:" which will be our n_iters
n_iters=$(echo "$result" | grep -oP 'into:\K[0-9]+')
echo "$result"
echo " "

# A constant is used above, and passed to the python script, to ensure
# that both this script and the python script use the same value for the
# file cutoff. This is important, as the python script will cut the file
# into n_iters pieces

if [ $n_iters -eq 0 ]; then
    n_iters=1
fi

# # Parse the variable list provided in $2
# # If "landuse" is in the list, set the
# # landuse flag to true and remove it from
# # the list
landuse_flag=0
IFS=','  # Set comma as the delimiter
for name in $2; do
    if [ "$name" == "landuse" ]; then
        landuse_flag=1
    else
        var_list+="${name},"
    fi
    # Remove the trailing comma
    var_list=${var_list%,}
done
unset IFS
# Check that the landuse master file exists
if [ $landuse_flag -eq 1 ]; then
    if [ ! -f "$LU_MASTER_FN" ]; then
        echo "$LU_MASTER_FN does not exist. Exiting."
        exit 1
    else
        echo "Regridding landuse data..."
        mpirun -np $3 ESMF_Regrid -s $LU_MASTER_FN -d .tmpdst1234.nc -m neareststod --src_var landuse --dst_var landuse -r -i --no_log
    fi
fi

#Regrid
echo "Regridding topographic data..."
for i in $(seq 1 $n_iters); do
    if [ $n_iters -gt 1 ]; then
        echo "Regridding piece $i of $n_iters"
    fi
    src_file=".tmpsrc_mosaic${i}.nc"

    # If we are cutting up the fields, change the naming convention
    if [ $n_iters -gt 1 ]; then
        # Initialize an empty result variable
        dst_var_list=""
        IFS=','  # Set comma as the delimiter
        for name in "$var_list"; do
            # Append 'i' to each name
            dst_var_list+="${name}${i},"
        done
        # Remove the trailing comma
        dst_var_list=${dst_var_list%,}
    else
        dst_var_list=$var_list
    fi
    mpirun -np $3 ESMF_Regrid -s $src_file -d .tmpdst1234.nc -m conserve --src_var $var_list --dst_var $dst_var_list -r -i --no_log
done

#Clean up the dst file if we saved multiple fields
if [ $n_iters -gt 1 ]; then
    echo "Cleaning up the dst file, merging the multiple fields calculated"
    python merge_fields.py .tmpdst1234.nc "${1%.nc}_final.nc" "$var_list" $n_iters
    rm .tmpdst1234.nc
else
    mv .tmpdst1234.nc "${1%.nc}_final.nc"
fi

#Clean up
rm -f .tmpsrc_mosaic*.nc
