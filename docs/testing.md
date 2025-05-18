# Testing

## Running Test Cases

To run a standard test case, run

```bash
cd build/
make test_cases
```

The output of the test case will be written to the directory `tests/Test_Cases`

### Test Cases on SLURM

To run the automated test case script on a system using SLURM, a modification to the cmake step is needed. This is because systems running SLURM often restrict direct access to the mpiexec command, and require additional information such as which account to bill node-hours to. To setup the automated test case script to run on a SLURM system, pass the flag `-DSRUN_FLAGS` to the cmake command as such:

```bash
cd build/
cmake ../ -DSRUN_FLAGS='-A 9999'
```

This will setup test cases to run with srun using the user account "9999". Any additional flags which could be passed to srun can also be included here, **with the exception of -N and -n, which are automatically set by the test case script**

## Testing model components (developers)

Some basic integration testing of the model dynamics are implemented, and can be run by calling

```bash
cd build/
make check
```

tests are defined under the `tests/` folder. `test_driver.F90` manages the execution of the different test modules, which are defined as `test_XXXX.F90`. Note that only some of the test cases are integrated with `test_driver.F90`, while most are legacy tests developed to work with ICAR.
