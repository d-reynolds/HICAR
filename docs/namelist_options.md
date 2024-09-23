## Namelist Options

You can get information about the namelist options by interacting with the executable generated from compiling the model.

Once the model has been succesfully compiled, from the root repository type:

```bash
./bin/HICAR
```

This will list the different user options available:

```bash
 Usage: ./HICAR [-v [variable_name ...|--all]] [--check-nml] [--gen-nml] namelist_file
     -v [variable_name ...|--all]: Print information about the namelist variable(s) variable_name, ... 
                                   --all prints out information for all namelist variables.
     --check-nml:                  Check the namelist file for errors without running the model.
     --gen-nml:                    Generate a namelist file with default values.
     namelist_file:                The name of the namelist file to use.
 
     Example to generate a namelist with default values:  ./HICAR --gen-nml namelist_file.nml
     Example to check namelist:                           ./HICAR --check-nml namelist_file.nml
     Example to run model:                                ./HICAR namelist_file.nml
     Example to learn about a namelist variable:          ./HICAR -v mp
     Example to generate namelist variable documentation: ./HICAR -v --all > namelist_doc.txt
```

In this way, the documentation for the model should stay tied to the version which was compiled. To create a custom namelist, the user is encouraged to follow the steps:

1. Generate default namelist
2. Read over the commented namelist options in the default namelist
3. use `./HICAR -v variable_name` as needed to list more information about namelist options of interest
4. Edit the default namelist, replacing default values with custom values where desired
5. Run `./HICAR --check-nml namelist_file.nml` on your custom namelist
