#!/bin/bash

#SBATCH --partition=standard
#SBATCH --job-name=Intrinsic_interface_lammps_job
#SBATCH --time=0:30:0
#SBATCH --account=tc075
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --qos=standard
#SBATCH --export=none

# Set the number of threads to 1
set -e
export OMP_NUM_THREADS=1

#Remove old lammps output
rm -f day2_* slurm*

source /etc/profile
module use /work/y07/shared/cirrus-ex/cirrus-ex-software/spack-cirrus-ex/0.2/cirrus-ex-cse/modules/Core
module load lammps/20250612
LMP_BIN="$(command -v lmp || true)"
if [ -z "$LMP_BIN" ]; then
    LMP_BIN="/mnt/lustre/e1000/home/y07/shared/cirrus-ex/cirrus-ex-software/spack-cirrus-ex/0.2/cirrus-ex-cse/opt/linux-rhel9-zen5/gcc-14.2/lammps-20250612-aqjcjka5m7cg2ck5auwcfqewujlbjpkf/bin/lmp"
fi
if [ ! -x "$LMP_BIN" ]; then
    echo "error: LAMMPS binary not found or not executable: $LMP_BIN" >&2
    exit 1
fi

# Launch the parallel job
srun --export=ALL -n 32 --hint=nomultithread --distribution=block:block "$LMP_BIN" -var rcut 0 -var dump 10 -var T 0.82 -var Lz_t 48 -var Lx_t 24 -var Ly_t 24 -var tslab 12 -var rho_l 0.746 -var rho_v 0.001 -var fmt lammps -var nprod 1000000 < liquid_vapour.in
