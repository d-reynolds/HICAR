# Running HICAR

## Example

```bash
mpiexec -np 2 ./HICAR/bin/HICAR HICAR_Test_Case.nml
```
In the above example, the number of MPI ranks is set with `-np 2`. The number of ranks must always be greater than 1, as at least 1 processor is needed for I/O. An even number of ranks may lead to inefficient domain decomposition, since an odd number of ranks are used in the domain decomposition (i.e. `-np 6` results in 1 I/O task and 5 compute tasks).

For more information about the namelist, see [namelist_options.md](namelist_options.md)

## Divisioning of I/O processes

HICAR uses asynchronous I/O to overlap read/write tasks with compute tasks. To accomplish this, the model divides the total number of MPI ranks at runtime into two groups: I/O tasks and compute tasks. By default, the model assigns **one I/O task per node**, with the remaining ranks on that node used for computation. This is why a single-node run reserves one rank for I/O and must be launched with at least two ranks (see the example above); a run with no compute ranks left over aborts at startup with an error.

### Controlling the number of I/O tasks per node

The number of I/O tasks assigned per node is controlled by the environment variable **`HICAR_IO_PER_NODE`**. If it is not set, HICAR uses one I/O task per node (the default). To dedicate additional ranks per node to I/O — for example on a run where a single I/O rank per node is a write bottleneck — set it before launching the model:

```bash
export HICAR_IO_PER_NODE=2
mpiexec -np 64 ./bin/HICAR your_namelist.nml
```

The value is validated at startup, and the following rules apply:

- `HICAR_IO_PER_NODE` must **evenly divide the number of MPI ranks per node**. If it does not, HICAR decrements it to the next value that does and prints a warning.
- It is capped at the number of ranks per node; a larger value is reduced to that number (with a warning).
- Every I/O task replaces a compute task, so increasing `HICAR_IO_PER_NODE` reduces the number of compute ranks. At least one compute rank must remain, or the run aborts at startup.

## Running on GPUs

A GPU build (`-DOPENACC=ON`, see [Compiling](compiling.md)) maps **one GPU to
each compute MPI rank**. Launch HICAR with as many compute ranks per node as
there are GPUs on that node (plus the one I/O rank per node that HICAR reserves
automatically). On most systems the MPI launcher or a small binding wrapper is
responsible for giving each rank its own GPU; consult your HPC documentation for
the recommended rank-to-GPU binding. Halo exchanges use NCCL when it is available
in the build.

## Example Slurm Script

If running HICAR on an HPC environment with a Slurm batch scheduler, the following Slurm script can be used as a template. Consult the documentation for your HPC environment should any issues arrise:

```bash
#!/bin/bash -l
#SBATCH --job-name="HICAR"
#SBATCH --time=00:05:00        # : Wall time of run
#SBATCH --output="HICAR.out"   # : File to write standard output to
#SBATCH --error="HICAR.err"    # : File to write standard error to
#SBATCH --mail-type=NONE
#SBATCH --nodes=2
#SBATCH --ntasks-per-core=1    # : Values greater than one turn hyperthreading on
#SBATCH --ntasks-per-node=64   # : Defines the number of MPI ranks per node
#SBATCH --cpus-per-task=1
#SBATCH --partition=PARTITION  # : The desired partition to use
#SBATCH --hint=nomultithread
#SBATCH --account=YOUR_ACCOUNT # : The account from which resources should be deducted
#SBATCH --mem=60GB             # : Amount of RAM per node to request
#SBATCH --begin=now

# These may vary, or just be unnecesarry, depending on your computing environment
export MPICH_OFI_STARTUP_CONNECT=1
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

srun --cpu-bind=verbose,cores ~/HICAR/bin/HICAR HICAR_Test_Case.nml
```
