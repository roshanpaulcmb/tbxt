#!/bin/bash
#SBATCH --job-name=namd3TbxtRest2
#SBATCH --output=slurm/%x_%j.out
#SBATCH --error=slurm/%x_%j.err
#SBATCH --nodes=1                # 4 replicas
#SBATCH --gres=gpu:4             # 1 GPU per replica
#SBATCH --ntasks-per-node=4      # 4 processes per node
#SBATCH --cpus-per-task=4        # 4 threads per process
#SBATCH --partition=dept_gpu,koes_gpu
##SBATCH --constraint="2080Ti|L40|A4500|A100"
#SBATCH --mail-user=rop174@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00

# -----------------------------
# Restart launcher for the TBXT REST2 run.
# Requires NAMD_HOME in env.
# Submit from the tbxt project root:  sbatch rest2/restart.sh
# -----------------------------
# Env knobs:
#   RUN_NAME  -- must match the initial run. Default: tbxtAf3g177d_run1.
#   NUM_RUNS  -- total run target. Default: 100000 (200 ns).
#                Override at submission: NUM_RUNS=500000 sbatch rest2/restart.sh
RUN_NAME="${RUN_NAME:-tbxtAf3g177d_run1}"
NUM_RUNS="${NUM_RUNS:-100000}"
LAUNCH_CONF="restart.conf"
export RUN_NAME NUM_RUNS

if [ -z "${NAMD_HOME:-}" ]; then
    echo "ERROR: NAMD_HOME is not set. Point it at your netlrts-smp-CUDA NAMD install, e.g.:" >&2
    echo "  export NAMD_HOME=\$HOME/software/NAMD_3.0.2_Linux-x86_64-netlrts-smp-CUDA" >&2
    exit 1
fi
if [ ! -x "$NAMD_HOME/namd3" ]; then
    echo "ERROR: $NAMD_HOME/namd3 not found or not executable" >&2
    exit 1
fi
if [ ! -f "$NAMD_HOME/lib/replica/REST2/rest2_remd.namd" ]; then
    echo "ERROR: REST2 driver missing at $NAMD_HOME/lib/replica/REST2/rest2_remd.namd" >&2
    exit 1
fi
export NAMD_HOME

# -----------------------------
# Working directories
# -----------------------------
SCRDIR=/scr/${SLURM_JOB_ID}
if ! mkdir -p $SCRDIR 2>/dev/null; then
    echo "WARNING: /scr full on $(hostname), falling back to home scratch" >&2
    SCRDIR=$HOME/scratch/${SLURM_JOB_ID}
    mkdir -p $SCRDIR
fi
cd $SCRDIR

rsync -a --exclude='runs/' $SLURM_SUBMIT_DIR/ ./

DEST=$SLURM_SUBMIT_DIR/runs
mkdir -p $DEST
RUNDIR=$SCRDIR/runs/$RUN_NAME

if [ ! -d "$DEST/$RUN_NAME" ]; then
    echo "ERROR: cannot restart -- $DEST/$RUN_NAME not found" >&2
    exit 1
fi
mkdir -p $SCRDIR/runs
cp -r "$DEST/$RUN_NAME" "$SCRDIR/runs/"
LOG_TAG="restart-$(date +%Y%m%d_%H%M%S)"

trap "cp -r $RUNDIR $DEST/ 2>/dev/null" EXIT

cd $SCRDIR/rest2

PPN=$SLURM_CPUS_PER_TASK
$NAMD_HOME/charmrun $NAMD_HOME/namd3 ++local +p $PPN \
    +replicas 4 +devicesperreplica 1 $LAUNCH_CONF \
    +stdout $RUNDIR/%d/$LOG_TAG.log
