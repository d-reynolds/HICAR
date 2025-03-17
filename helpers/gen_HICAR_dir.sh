#!/bin/bash

#Call like: ./gen_HICAR_dir.sh path/to/desired/parent/directory path/to/HICAR/repo
# i.e. :
#./gen_HICAR_dir.sh ./Model_runs/ /home/user/HICAR/

parent_dir=$(realpath $1)
HICAR_dir=$(realpath $2)

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

if [ ! -f ./VEGPARM.TBL ]; then
	echo 'Copying .TBL files needed by NoahMP, which are found in'
	echo $HICAR_dir/run
	echo 'to ./input'

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

	cd $HICAR_dir
	cp -r run/rrtmg_support $parent_dir/HICAR/input
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
