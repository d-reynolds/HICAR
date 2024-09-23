# Running HICAR

## Example

```bash
mpiexec -np 2 ./HICAR/bin/HICAR HICAR_Test_Case.nml
```
In the above example, the number of MPI ranks is set with `-np 2`. The number of ranks must always be greater than 1, as at least 1 processor is needed for I/O. An even number of ranks may lead to inefficient domain decomposition, since an odd number of ranks are used in the domain decomposition (i.e. `-np 6` results in 1 I/O task and 5 compute tasks).

For more information about the namelist, see [namelist_options.md](namelist_options.md)

## Divisioning of I/O processes

HICAR uses asynchronous I/O to overlap read/write tasks with compute tasks. To accomplish this, the model divides the total number of MPI ranks at runtime into two groups: I/O tasks and compute tasks. Currently, the code will assign one I/O task per node, with the rest of the tasks on a node being compute tasks. This means that local runs will only have one I/O task, and that the user cannot control the division of I/O tasks. This may be added as a feature in the future.

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
