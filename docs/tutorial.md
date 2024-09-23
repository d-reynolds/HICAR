
# Tutorial

This tutorial is designed to get you accustomed to launching a run with HICAR.

## Compiling

First, consult the information in [compiling.md](compiling.md) to compile the model.

## Create working directories

Once the model has been succesfully compiled, from the root repository type:

`./helpers/gen_HICAR_dir.sh path/to/desired/parent/directory path/to/HICAR/repo`

Replacing the two paths with their actual values in your filesystem. Importantly, the parent directory cannot be the parent directory of the HICAR code repo. Answer "y" to all prompts posed by the script. This script will generate a directory tree at `path/to/desired/working/directory` which organises the input and output data used when running HICAR. The directories generated are:

<pre>

        PARENT_DIR/HICAR         <-- root working directory
        PARENT_DIR/HICAR/input   <-- folder to contain namelists, forcing file lists, and supporting files
        PARENT_DIR/HICAR/output  <-- folder to store model output
        PARENT_DIR/HICAR/restart <-- folder to store model restart files
        PARENT_DIR/HICAR/forcing <-- folder to store forcing files
        PARENT_DIR/HICAR/domains <-- folder to store domain files (static input)
</pre>

This script will automatically populate the input folder with the supporting files needed by HICAR, and download the HICAR Test Data repo (https://github.com/HICAR-Model/Test-Data/tree/main), populating the `forcing` and `domains` folders with forcing data and a domain file.

## Running the model

Once the `gen_HICAR_dir.sh` script has finished, all of the input data and file structure should be in place for launching a test run. Navigate to the `PARENT_DIR/HICAR/input` folder generated in the last step and type:

`mpirun -np 2 /absolute/path/to/bin/HICAR HICAR_Test_Case.nml`

For more information on how to run the model, including a template Slurm script, see [running.md](running.md)

## Graduating

Congratulations, if you made it this far you've likely launched your first HICAR run. To continue with the model and setup your own custom runs, you can consult the other documentation here:
<pre>
- domain_generation.md    <-- for creating your own domains to use with HICAR
- forcing_data.md         <-- information on the requirements for forcing data, and how to generate forcing file lists
- namelist_options.md     <-- for creating custom namelists and getting information on namelist options
</pre>
