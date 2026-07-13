#!/bin/bash
#SBATCH --job-name=namd3TbxtRest2
#SBATCH --output=slurm/%x_%j.out
#SBATCH --error=slurm/%x_%j.err
#SBATCH --nodes=1                # 8 replicas
#SBATCH --gres=gpu:8             # 1 GPU per replica
#SBATCH --ntasks-per-node=8      # 8 process per node
#SBATCH --cpus-per-task=4        # 4 threads mapping to 4 cores per node
#SBATCH --partition=dept_gpu,koes_gpu
##SBATCH --constraint="2080Ti|L40|A4500|A100"
#SBATCH --mail-user=rop174@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00          # adjust based on estimated runtime

# -----------------------------
# Requires NAMD_HOME in env. Recommended one-time setup on each machine:
#   echo 'export NAMD_HOME=$HOME/software/NAMD_3.0.2_Linux-x86_64-netlrts-smp-CUDA' >> ~/.bashrc
# Submit from the tbxt project root:  sbatch rest2/rest2.sh
# -----------------------------
# Other env knobs (used by the python wrapper):
#   RUN_NAME     -- matches $run_name in initTbxtAf3g177d.conf. Use a fresh
#                   value per logical run (e.g. tbxtAf3g177d_$(date +...)).
#   LAUNCH_CONF  -- job0.conf for the initial run, restart.conf to continue.
#                   Default: job0.conf.
#   EQ_PROTEIN   -- (initial run only) subdir under equilibrate/ to source
#                   equinpt.coor + equinpt.xsc from. Default: tbxtAf3g177d.
#   NUM_RUNS     -- (restart only) override the new num_runs target. If unset,
#                   restart.conf extends by the init's num_runs each pass.
RUN_NAME="${RUN_NAME:-tbxtAf3g177d_run3}"
LAUNCH_CONF="${LAUNCH_CONF:-job0.conf}"
EQ_PROTEIN="${EQ_PROTEIN:-tbxtAf3g177d}"
export RUN_NAME

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
export NAMD_HOME   # so job0.conf / restart.conf can source via $env(NAMD_HOME)

# -----------------------------
# Working directories
# -----------------------------
# Run in place on shared storage: every file the driver writes -- per-replica
# checkpoints, the global restart .tcl, DCDs, histories -- lands directly under
# runs/$RUN_NAME as it is written. No node-local /scr staging and no copy-back,
# so a walltime SIGKILL (or any kill) cannot truncate results on the way home.
DEST=$SLURM_SUBMIT_DIR/runs
RUNDIR=$DEST/$RUN_NAME
mkdir -p $DEST

if [ "$LAUNCH_CONF" = "job0.conf" ]; then
    # Stage equinpt.coor + equinpt.xsc from the protein's equilibration dir
    # into rest2/ -- the REST2 driver reads them by hardcoded name from CWD.
    EQ_SRC=$SLURM_SUBMIT_DIR/equilibrate/$EQ_PROTEIN
    if [ ! -f "$EQ_SRC/equinpt.coor" ] || [ ! -f "$EQ_SRC/equinpt.xsc" ]; then
        echo "ERROR: missing equinpt.{coor,xsc} in $EQ_SRC" >&2
        exit 1
    fi
    cp "$EQ_SRC/equinpt.coor" "$EQ_SRC/equinpt.xsc" "$SLURM_SUBMIT_DIR/rest2/"

    # Pre-create empty per-replica output subdirs.
    # Keep this seq in sync with num_replicas in initTbxtAf3g177d.conf.
    for i in $(seq 0 7); do
        mkdir -p $RUNDIR/$i
    done
    LOG_TAG="job0"
else
    # Restart: the prior run tree is already on shared storage; the driver
    # finds the latest checkpoint .tcl + .coor/.vel/.xsc directly under it.
    if [ ! -d "$RUNDIR" ]; then
        echo "ERROR: cannot restart -- $RUNDIR not found" >&2
        exit 1
    fi
    LOG_TAG="restart-$(date +%Y%m%d_%H%M%S)"
fi

# Run NAMD from inside rest2/ so ../structures, ../toppar, ../runs in the
# config files resolve as siblings.
cd $SLURM_SUBMIT_DIR/rest2

PPN=$SLURM_CPUS_PER_TASK
$NAMD_HOME/charmrun $NAMD_HOME/namd3 ++local +p $PPN \
    +replicas 8 +devicesperreplica 1 $LAUNCH_CONF \
    +stdout $RUNDIR/%d/$LOG_TAG.log
