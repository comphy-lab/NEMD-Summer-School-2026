#!/usr/bin/env bash
# submit_nequil_matrix.sh -- vanilla conductance equilibration sweep on Cirrus.
#
# Runs independent conductance.in jobs for several nequil values and seeds inside
# one Slurm allocation, renaming the output from each run so the analysis can ask
# whether the heat-bath imbalance decreases monotonically or only on average.
#
# Usage from day1/conductance:
#   ./submit_nequil_matrix.sh
#   ./submit_nequil_matrix.sh "10000 40000 80000 160000" "90187 90287" 40000
set -euo pipefail

if [ ! -f conductance.in ] || [ ! -f ../shared_setup.lmp ]; then
  echo "error: run this from day1/conductance" >&2
  exit 1
fi

NEQUILS="${1:-10000 40000 80000 160000}"
SEEDS="${2:-90187 90287}"
NPROD="${3:-40000}"

# Keep old named matrix results until the new Slurm job starts, then clear only
# this matrix's outputs. Do not delete the base day1_* run outputs here.
JOBID=$(sbatch --parsable <<EOF
#!/bin/bash -l
#SBATCH --job-name=day1-nq
#SBATCH --account=tc075
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --export=none

set -eo pipefail
export OMP_NUM_THREADS=1

source /etc/profile
set -u
module use /work/y07/shared/cirrus-ex/cirrus-ex-software/spack-cirrus-ex/0.2/cirrus-ex-cse/modules/Core
module load lammps/20250612
LMP_BIN="\$(command -v lmp || true)"
if [ -z "\$LMP_BIN" ]; then
    LMP_BIN="/mnt/lustre/e1000/home/y07/shared/cirrus-ex/cirrus-ex-software/spack-cirrus-ex/0.2/cirrus-ex-cse/opt/linux-rhel9-zen5/gcc-14.2/lammps-20250612-aqjcjka5m7cg2ck5auwcfqewujlbjpkf/bin/lmp"
fi
if [ ! -x "\$LMP_BIN" ]; then
    echo "error: LAMMPS binary not found or not executable: \$LMP_BIN" >&2
    exit 1
fi
echo "LAMMPS: \$LMP_BIN"

rm -f day1_nequil_*
printf "# case nequil seed nprod\\n" > day1_nequil_cases.txt

for nq in $NEQUILS; do
  for seed in $SEEDS; do
    case_id="nq\${nq}_seed\${seed}"
    echo "==> \$case_id"
    printf "%s %s %s %s\\n" "\$case_id" "\$nq" "\$seed" "$NPROD" >> day1_nequil_cases.txt

    rm -f day1_Tz.profile day1_heat.dat day1_params.txt day1_heatflux.profile day1_conductance.png log.lammps
    srun --export=ALL --hint=nomultithread --distribution=block:block "\$LMP_BIN" \
      -in conductance.in -var nequil "\$nq" -var nprod "$NPROD" -var seed "\$seed"

    mv day1_Tz.profile "day1_nequil_\${case_id}_Tz.profile"
    mv day1_heat.dat "day1_nequil_\${case_id}_heat.dat"
    mv day1_params.txt "day1_nequil_\${case_id}_params.txt"
    mv log.lammps "day1_nequil_\${case_id}_log.lammps"
  done
done

echo "Done. Analyse with: python analyse_conductance_nequil.py"
EOF
)

echo "Submitted nequil matrix job $JOBID. Analyse after it finishes:"
echo "    python analyse_conductance_nequil.py"
