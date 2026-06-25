#!/bin/bash

#Call like: ./gen_HICAR_dir.sh path/to/desired/parent/directory path/to/HICAR/repo
# i.e. :
#./gen_HICAR_dir.sh ./Model_runs/ /home/user/HICAR/

# check if the user supplied the correct number of arguments
if [ "$#" -ne 2 ]; then
	echo
	echo "Usage: $0 <parent_directory> <HICAR_repo_directory>"
	echo
	exit 1
fi

# check if the directory at path $1 exists
if [ ! -d "$1" ]; then
	mkdir -p $1
fi

# check if the directory at path $2 exists
if [ ! -d "$2" ]; then
	echo
	echo "Error: HICAR repo directory $2 does not exist."
	echo
	exit 1
else [ -d "$2/bin" ]; then
	echo
	echo "Error: HICAR repo directory $2 does not contain a bin/ directory."
	echo
	exit 1
fi

# check for any files in $2/bin starting with HICAR
found_hicar=false
for _f in "$2"/bin/HICAR*; do
	if [ -e "$_f" ]; then
		found_hicar=true
		break
	fi
done

if [ "$found_hicar" = false ]; then
	echo
	echo "Error: HICAR repo directory $2 does not contain a HICAR executable in bin/"
	echo
	exit 1
fi

# Use readlink -f or realpath to get absolute paths
parent_dir=$(readlink -f "$1")
HICAR_dir=$(readlink -f "$2")

echo
echo '#######################################################'
echo '############# Setting up HICAR directory ##############'
echo '#######################################################'
echo
echo '   This script will create the following directories:'
echo '        '$parent_dir'/HICAR'
echo '        '$parent_dir'/HICAR/input'
echo '        '$parent_dir'/HICAR/output'
echo '        '$parent_dir'/HICAR/restart'
echo '        '$parent_dir'/HICAR/forcing'
echo '        '$parent_dir'/HICAR/domains'
echo
echo '  If directories in the file structure already exist,  '
echo '             they will not be overwritten.  '
echo
echo '#######################################################'
echo
#Create the parent directory if it doesn't exist and enter it

cd $parent_dir
if [ ! -d ./HICAR ]; then
	mkdir HICAR
else
	echo
	echo 'HICAR directory already exists, continuing...'
	echo
fi
cd HICAR
########## INPUT ####################
if [ ! -d ./input ]; then
	echo 'Creating Input Directory (./input)'
	mkdir input
fi
cd input

if [ ! -f ./NoahmpTable.TBL ]; then
	echo 'Copying .TBL files needed by NoahMP, which are found in'
	echo $HICAR_dir/run
	echo 'to ./input'

	# check if the .TBL files exist in the HICAR repo run directory, if not,
	# ask the user to re-run the configure script to download them
	if [ ! -f $HICAR_dir/run/NoahmpTable.TBL ]; then
		echo
		echo 'Error: .TBL files not found in the HICAR repo run directory.'
		echo 'Please re-run the cmake configuration step to download them.'
		echo
		exit 1
	fi
	cp $HICAR_dir/run/*.TBL ./
fi

# See if the uesr has already cloned the supporting files
if [ ! -d ./rrtmg_support -o ! -d ./mp_support ]; then

	# Fetch HICAR supporting files from run to the input dir:
	echo
	echo 'Copying supporting files from the HICAR repo /run'
	echo 'directory to ./input'
	echo
	echo 'These files were downloaded at configure time from the repo:'
	echo 'https://github.com/NCAR/icar_supporting_files.git'
	echo 

	# check if the physics support folders exist in the HICAR repo run directory, if not,
	# ask the user to re-run the configure script to download them
	if [ ! -d $HICAR_dir/run/rrtmg_support ]; then
		echo
		echo 'Error: rrtmg_support directory not found in the HICAR repo run directory.'
		echo 'Please re-run the cmake configuration step to download them.'
		echo
		exit 1
	fi

	if [ ! -d $HICAR_dir/run/rrtmgp_support ]; then
		echo
		echo 'Error: rrtmgp_support directory not found in the HICAR repo run directory.'
		echo 'Please re-run the cmake configuration step to download them.'
		echo
		exit 1
	fi

	if [ ! -d $HICAR_dir/run/mp_support ]; then
		echo
		echo 'Error: mp_support directory not found in the HICAR repo run directory.'
		echo 'Please re-run the cmake configuration step to download them.'
		echo
		exit 1
	fi

	cd $HICAR_dir
	cp -r run/rrtmg_support $parent_dir/HICAR/input
	cp -r run/rrtmgp_support $parent_dir/HICAR/input
	cp -r run/mp_support $parent_dir/HICAR/input

fi
####################################
cd $parent_dir/HICAR

if [ ! -d ./output ]; then
echo 'Creating Output Directory (./output)'
mkdir output
fi

if [ ! -d ./restart ]; then
echo 'Creating Restart Files Directory (./restart)'
mkdir restart
fi

if [ ! -d ./forcing ]; then
echo 'Creating Forcing Data Directory (./forcing)'
mkdir forcing
fi

if [ ! -d ./domains ]; then
	echo 'Creating Domains Directory  (./domains)'
	mkdir domains
fi

echo Setup of HICAR directory complete
