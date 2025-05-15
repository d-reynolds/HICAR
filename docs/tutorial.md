
# Tutorial

This tutorial is designed to get you accustomed to launching a run with HICAR.

## Compiling

First, consult the information in [compiling.md](compiling.md) to compile the model.

## Testing

To test the installation, run

```bash
cd build/
make test_cases
```

This will launch a small test case to ensure that libraries have been succesfully installed and linked. If you are installing on a system with the SLURM job scheduler, follow the steps mentioned in [testing.md](testing.md).

## Create working directories

Once the model has been succesfully compiled, from the root repository type:

`./helpers/gen_HICAR_dir.sh path/to/desired/parent/directory path/to/HICAR/repo`

Replacing the two paths with their actual values in your filesystem. Importantly, the parent directory cannot be the parent directory of the HICAR git repo. This script will generate a directory tree at `path/to/desired/working/directory` which organises the input and output data used when running HICAR. The directories generated are:

<pre>

        PARENT_DIR/HICAR         <-- root working directory
        PARENT_DIR/HICAR/input   <-- folder to contain namelists, forcing file lists, and supporting files
        PARENT_DIR/HICAR/output  <-- folder to store model output
        PARENT_DIR/HICAR/restart <-- folder to store model restart files
        PARENT_DIR/HICAR/forcing <-- folder to store forcing files
        PARENT_DIR/HICAR/domains <-- folder to store domain files (static input)
</pre>

This script will automatically populate the input folder with the supporting files needed by HICAR.

## Input data

To setup your own custom runs, you can consult the other documentation here:

<pre>
- domain_generation.md    <-- for creating your own domains to use with HICAR
- forcing_data.md         <-- information on the requirements for forcing data, and how to generate forcing file lists
- namelist_options.md     <-- for creating custom namelists and getting information on namelist options
</pre>

For more information on how to run the model, including a template Slurm script, see [running.md](running.md)
