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

# Fetch HICAR supporting files from git repo:
# Ask user if they want us to fetch supporting files
echo
echo Would you like to install supporting files needed 
echo for the RRTMG radiaiton scheme
echo and the ISHMAEL microphysics scheme?
echo 
echo This will clone the icar_supporting_files git
echo repo to the parent directory of the HICAR repo
echo
read -p '>>(y/n):' ans

if [ "$ans" = "y" ];then
	cd $HICAR_dir
	cd ..
	if [ ! -d ./icar_supporting_files ]; then
		git clone https://github.com/NCAR/icar_supporting_files.git
	fi
	cp -r icar_supporting_files/rrtmg_support $parent_dir/HICAR/input
	cp -r icar_supporting_files/mp_support $parent_dir/HICAR/input

	echo
	echo Supporting directories installed to ./input/mp_support and
	echo ./input/rrtmg_support. 
	echo They need to be in the same directory as the namelist 
	echo used for a HICAR run to be found.
else
	#If the user decided not to...
	echo
	echo OK, you can always find them later at: 
	echo https://github.com/NCAR/icar_supporting_files.git
fi
echo
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
echo Should we install the default test case as well?
echo
echo This will clone the HICAR_Test_Data git
echo repo to the parent directory of the HICAR repo.
echo
read -p '>>(y/n):' ans
echo

if [ "$ans" = "y" ];then

	if [ ! -d ./output/TestCase ]; then
	echo 'Creating Output/TestCase Directory (./output/TestCase)'
	mkdir output/TestCase
	fi
	if [ ! -d ./restart/TestCase ]; then
	echo 'Creating Restart/TestCase Directory (./restart/TestCase)'
	mkdir restart/TestCase
	fi
	if [ ! -d ./forcing/COSMO_2017 ]; then
	echo 'Creating Test Case Forcing Data Directory (./forcing/COSMO_2017)'
	mkdir forcing/COSMO_2017
	fi

	cd $HICAR_dir
	cd ..
	if [ ! -d ./HICAR_Test_Data ]; then
		git clone https://github.com/HICAR-Model/Test-Data.git HICAR_Test_Data
	fi

	echo 'Copying domain file (Gaudergrat_250m.nc) to ./domains'
	cp HICAR_Test_Data/static_data/Gaudergrat_250m.nc $parent_dir/HICAR/domains
	
	echo 'Generating namelist (HICAR_Test_Case.nml) to ./input'
	$HICAR_dir/bin/HICAR --gen-nml $HICAR_dir/run/default.nml
	./HICAR_Test_Data/gen_TestCase.sh $HICAR_dir/run/default.nml $parent_dir/HICAR/input/HICAR_Test_Case.nml

	echo 'Copying forcing data to ./forcing/COSMO_2017'
	cp HICAR_Test_Data/input_data/laf*.nc $parent_dir/HICAR/forcing/COSMO_2017

	cd $parent_dir/HICAR

	echo
	echo Creating forcing file list.
	echo You can do this in the future
	echo by running the following command:
	echo
	echo $HICAR_dir'/helpers/filelist_script.sh "forcing/COSMO_2017/laf*.nc" input/file_list_TestCase.txt'
	echo
	if [ -f ./input/file_list_TestCase.txt ]; then
		rm ./input/file_list_TestCase.txt
	fi

	$HICAR_dir/helpers/filelist_script.sh "forcing/COSMO_2017/laf*.nc" input/file_list_TestCase.txt

	echo
	echo Test case installed, for templates of job scripts, see:
	echo $HICAR_dir'/helpers'
	echo
fi

echo Setup of HICAR directory complete
