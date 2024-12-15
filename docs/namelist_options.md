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

## Nested Runs

HICAR also supports one-way nesting of domains. The rules of nesting are that each domain can only have one parent nest, and the forcing data provided to a run must encompass all of the domains. One can also setup multple "nesting chains", where multiple domains are nested within a single parent nest. Nesting is controlled via the `nests` and `parent_nest` options in the `&general` section of the namelist. "Nests" tells HICAR how many domains to run, and `parent_nest` specifies which domain is the parent nest of a given domain. When multiple domains are specified, the user must explicitly provide the path to the initial_conditions_file for each domain as an ordered list of strings. This is done via the `initial_conditions_file` option in the `&domain` section of the namelist. Additionally, the `dx` option in the `&domain` section must be specified for each domain.

For example, to run with 3 nests, the following options could be set:

```
! In the &domain section
&domain

    ! Note that each domain has the domain id corresponding to its position in this list.
    ! For example, '../path/to/file1.nc' will have id=1, '../path/to/file2.nc' will have id=2
    ! And so on...
    initial_conditions_file = '../path/to/file1.nc', '../path/to/file2.nc', '../path/to/file3.nc'

    ! Set the horizontal resolution for each of these domains
    dx = 1000.0, 500.0, 100.0

    !  All other desired domain options...

/

! In the &general section
&general

    ! Run 3 domains
    nests = 3

    ! parent domain id = 0 -- Read in forcing data from files
    !              ⌄
    ! parent domain id = 1 -- Use nest #1 to provide boundary conditions
    !                 ⌄
    ! parent domain id = 2 -- Use nest #2 to provide boundary conditions
    !                    ⌄
    parernt_nest = 0, 1, 2

    ! All other general options...
/
```

The user can also provide different namelist options for each domain by specifying the desired options as an ordered list, similar to the example above with `initial_conditions_file`. If a single option is specified by the user, it is applied to all domains in the run. If no option is specified, then the model-default value is used (see default namelist, mentioned above). Most all namelist options support this convention. Namelist options marked with a "(-)" in their description only support a single value, and cannot differ among the different domains. Namelist options marked with a "(*)" in their description must be explicitly set for all domains.