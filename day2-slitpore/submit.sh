#!/bin/bash
# submit.sh -- submit one Day-2 slit-pore run to Cirrus as a batch job.
#
#   Run it FROM INSIDE the measurement folder (poiseuille/) so the job starts
#   where `include ../shared_setup.lmp` resolves:
#
#     cd poiseuille && ../submit.sh poiseuille.in                 # wide pore (w=10)
#     cd poiseuille && ../submit.sh poiseuille.in -var width 4    # the narrow-pore Push
#
# It writes a short SLURM job and sends it to the queue with sbatch. Check
# progress with `squeue --me`; LAMMPS runs on a compute node and its output
# lands in slurm-<jobid>.out. When the job has finished:
#     python analyse_poiseuille.py
set -e

if [ "$#" -lt 1 ]; then
  echo "usage: cd poiseuille && ../submit.sh poiseuille.in [-var width 4] [-var NAME value ...]"
  exit 1
fi

IN="$1"   # the .in file (any -var arguments follow it)

# Preflight: refuse to submit from the wrong place so a job cannot fail late on a
# compute node (wasting class time) on a missing include.
if [ ! -f "$IN" ]; then
  echo "error: input '$IN' not found here -- cd into the poiseuille/ folder and run ../submit.sh from inside it."
  exit 1
fi
if [ ! -f ../shared_setup.lmp ]; then
  echo "error: ../shared_setup.lmp not found -- run this from INSIDE poiseuille/, not the top level."
  exit 1
fi

# Refuse to submit from /home: Cirrus compute nodes can read /work but not /home, so a job
# launched from a /home folder passes the checks above but fails late on the compute node.
case "$PWD" in
  /work/tc075/tc075/*) ;;
  *)
    echo "error: Cirrus jobs must be submitted from /work/tc075/tc075/USERNAME, not:"
    echo "       $PWD"
    echo "       cd /work/tc075/tc075/$USER and run from the poiseuille/ folder there."
    exit 1
    ;;
esac

# Determine the width for this run (default 10). Clear ONLY this width's outputs: the wide and
# narrow results are compared side by side at the end of the sheet, so a w=4 run must
# not delete the w=10 data; the outputs are width-tagged.
# Stale same-width data would otherwise be read as a silently-wrong result, so the
# current width's outputs are still cleared before submitting.
W=10
prev=""
for a in "$@"; do
  if [ "$prev" = "width" ]; then W="$a"; fi
  prev="$a"
done
rm -f "poiseuille_w${W}.profile" "poiseuille_w${W}.params.txt" "poiseuille_w${W}.png" \
      "poiseuille_w${W}_traj.xyz"
rm -f "log.lammps"

printf -v LMP_ARGS ' %q' "$@"   # input + any -var args, shell-quoted so spaces/globs survive

# The default nequil/nprod production runs need more than Day-1's 20-minute walltime;
# 45 minutes covers them on one core. Single-task keeps run-to-run results reproducible.
JOBID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=day2pore
#SBATCH --account=tc075
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --export=none

# Clean job environment (--export=none) so the Python modules you loaded for the
# analysers cannot leak in and clash with 'module load lammps'.
export OMP_NUM_THREADS=1

module load lammps
srun --hint=nomultithread --distribution=block:block lmp -in$LMP_ARGS
EOF
)
echo "Submitted job $JOBID (width w=${W}). Wait until 'squeue --me' no longer lists it, THEN run:"
echo "    python analyse_poiseuille.py"
echo "(running it before the job finishes reports 'not found'.)"
