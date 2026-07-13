#!/bin/bash
#SBATCH --job-name=namd3TbxtRest2
#SBATCH --output=slurm/%x_%j.out
#SBATCH --error=slurm/%x_%j.err
#SBATCH --nodes=1                # 8 replicas
#SBATCH --gres=gpu:8             # 1 GPU per replica
#SBATCH --ntasks-per-node=8      # 8 processes per node
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
#   RUN_NAME  -- must match the initial run. Default: tbxtAf3g177d_run3.
#   NUM_RUNS  -- total run target. Default: 100000 (200 ns).
#                Override at submission: NUM_RUNS=500000 sbatch rest2/restart.sh
RUN_NAME="${RUN_NAME:-tbxtAf3g177d_run3}"
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
# Run in place on shared storage: every checkpoint the driver writes -- the
# per-replica .coor/.vel/.xsc/.tcl and the global restart .tcl -- lands directly
# under runs/$RUN_NAME as it is written. No node-local /scr staging and no
# copy-back, so a walltime SIGKILL (or any kill) cannot truncate the restart
# files on the way home.
DEST=$SLURM_SUBMIT_DIR/runs
RUNDIR=$DEST/$RUN_NAME

if [ ! -d "$RUNDIR" ]; then
    echo "ERROR: cannot restart -- $RUNDIR not found" >&2
    exit 1
fi
LOG_TAG="restart-$(date +%Y%m%d_%H%M%S)"

cd $SLURM_SUBMIT_DIR/rest2

PPN=$SLURM_CPUS_PER_TASK
$NAMD_HOME/charmrun $NAMD_HOME/namd3 ++local +p $PPN \
    +replicas 8 +devicesperreplica 1 $LAUNCH_CONF \
    +stdout $RUNDIR/%d/$LOG_TAG.log
