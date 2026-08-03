#!/bin/bash
#SBATCH --job-name=namd3TbxtRest2Mn
#SBATCH --output=slurm/%x_%j.out
#SBATCH --error=slurm/%x_%j.err
#SBATCH --nodes=2                # 8 replicas split 4+4
#SBATCH --gres=gpu:4             # per node; 1 GPU per replica
#SBATCH --ntasks-per-node=4      # 4 replica processes per node
#SBATCH --cpus-per-task=4
#SBATCH --partition=dept_gpu,koes_gpu
# Square brackets force EVERY allocated node to match the SAME option. A bare
# "A|B" lets each node pick independently, which would give a mixed ladder
# (e.g. 4 L40 + 4 2080Ti) -- every exchange then blocks on the slowest replica.
# Only features with >=2 nodes carrying >=4 GPUs can satisfy --nodes=2.
#SBATCH --constraint="[L40|2080Ti|A4500|A4000|Xp]"
#SBATCH --exclude=g021           # g021 GPU4 (pci A5:00.0): 3026 uncorrected SRAM parity ECC
                                 # errors, intermittently refuses CUDA contexts (job 56910571).
#SBATCH --mail-user=rop174@pitt.edu
#SBATCH --mail-type=ALL
# 25000 runs needs ~15.3 h at the ~1640 runs/h measured on 2x4 L40 (job 56913090),
# so 12:00:00 walls out around run ~19600. 24 h covers a full ladder in one shot.
# Shorter jobs backfill sooner -- override on the command line (sbatch --time=...)
# for continuations that only need a few hours.
#SBATCH --time=24:00:00

# -----------------------------------------------------------------------------
# Multi-node REST2. Same science as rest2.sh, but spreads the 8-replica ladder
# over 2 nodes x 4 GPUs instead of demanding 8 GPUs on a single node.
#
# Why: 8-GPU single-node is the rarest request on this cluster (14 eligible
# nodes) and has cost 8-10 day queue waits. 4 GPUs x 2 nodes schedules far more
# easily. REST2 tolerates the split because an exchange moves a scalar energy
# and a parameter ID -- bytes, not coordinates -- so the 10 GbE interconnect
# (there is no InfiniBand on this cluster) is not a bottleneck.
#
# Launch mechanism: charmrun's ++mpiexec delegates process spawning to srun, so
# Slurm starts all 8 processes across both nodes natively. This needs NO ssh
# keys between compute nodes (charmrun's default remote shell is ssh, and
# ~/.ssh/id_rsa is not in authorized_keys here) and NO MPI build of NAMD -- the
# stock netlrts-smp-CUDA binary is fine. netlrts IS a multi-node TCP/IP layer;
# only the ++local flag in rest2.sh made it single-node.
#
# Submit from the tbxt project root:  sbatch rest2/rest2netlrts.sh
# Env knobs are the same as rest2.sh (RUN_NAME / LAUNCH_CONF / EQ_PROTEIN /
# NUM_RUNS) plus MAX_TEMP. RUN_NAME defaults to a DIFFERENT run than rest2.sh so
# the two can be queued simultaneously without both writing into the same tree.
#
# run3 (hotter ladder, 310 -> 376 K, targeting ~20% acceptance):
#   sbatch --export=ALL,RUN_NAME=tbxtAf3g177d_run3,MAX_TEMP=376 rest2/rest2netlrts.sh
# -----------------------------------------------------------------------------
NUM_REPLICAS=8   # keep in sync with num_replicas in initTbxtAf3g177d.conf
RUN_NAME="${RUN_NAME:-tbxtAf3g177d_run4}"
LAUNCH_CONF="${LAUNCH_CONF:-job0.conf}"
EQ_PROTEIN="${EQ_PROTEIN:-tbxtAf3g177d}"
# Top of the ladder. Default 360 reproduces run4 exactly; run3 overrides to 376.
MAX_TEMP="${MAX_TEMP:-360}"
export RUN_NAME MAX_TEMP

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

DEST=$SLURM_SUBMIT_DIR/runs
RUNDIR=$DEST/$RUN_NAME
mkdir -p $DEST

# =============================================================================
# PREFLIGHT -- everything below prints to the slurm .out before NAMD starts, so
# a failed run is diagnosable from the log alone without needing the node back.
# =============================================================================
echo "############ PREFLIGHT $(date -Is) ############"
echo "job            : $SLURM_JOB_ID"
echo "nodes          : $SLURM_JOB_NODELIST"
echo "tasks/node     : $SLURM_NTASKS_PER_NODE   cpus/task: $SLURM_CPUS_PER_TASK"
echo "replicas       : $NUM_REPLICAS"
echo "run name       : $RUN_NAME"
echo "ladder         : 310 -> $MAX_TEMP K (geometric, $NUM_REPLICAS rungs)"
echo "launch conf    : $LAUNCH_CONF"
echo "NAMD_HOME      : $NAMD_HOME"
echo "restart count  : ${SLURM_RESTART_COUNT:-0}"

HOSTS=$(scontrol show hostnames "$SLURM_JOB_NODELIST")
echo "--- node features (must be identical; bracket constraint should guarantee it) ---"
for h in $HOSTS; do
    echo "  $h : $(scontrol show node $h | grep -o 'ActiveFeatures=[^ ]*')"
done
NFEAT=$(for h in $HOSTS; do scontrol show node $h | grep -o 'ActiveFeatures=[^ ]*'; done | sort -u | wc -l)
if [ "$NFEAT" -ne 1 ]; then
    echo "WARNING: allocated nodes are NOT homogeneous -- ladder will run at the speed" >&2
    echo "         of the slowest GPU and exchange stats will be skewed." >&2
fi

# --- GPU health gate -------------------------------------------------------
# Job 56910571 died because one allocated GPU had 3026 uncorrected ECC errors
# and refused a CUDA context. Catching that here costs ~10s and lets us hand the
# allocation back via requeue instead of burning it.
GPUCHK=$RUNDIR.gpucheck.$SLURM_JOB_ID.sh
mkdir -p "$(dirname "$GPUCHK")"
# Quoted heredoc: nothing below is expanded by the submitting shell. Keeping this
# in a file rather than inline avoids nested srun/bash -c/awk quoting hazards.
cat > "$GPUCHK" <<'GPUCHKEOF'
#!/bin/bash
H=$(hostname)
if [ "$1" = "report" ]; then
    echo "  [$H] CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
    nvidia-smi --query-gpu=index,pci.bus_id,name,compute_mode,memory.used,ecc.errors.uncorrected.volatile.total \
               --format=csv,noheader | sed "s/^/  [$H] /"
else
    nvidia-smi --query-gpu=index,ecc.errors.uncorrected.volatile.total \
               --format=csv,noheader,nounits |
    while IFS=, read -r idx ecc; do
        ecc=$(echo "$ecc" | tr -d ' ')
        case "$ecc" in ''|'[N/A]'|'N/A'|0) ;; *) echo "$H:gpu$idx=$ecc" ;; esac
    done
fi
GPUCHKEOF
chmod +x "$GPUCHK"

echo "--- GPU inventory + ECC on every allocated node ---"
srun --ntasks-per-node=1 --ntasks=$SLURM_NNODES "$GPUCHK" report 2>&1

BADGPU=$(srun --ntasks-per-node=1 --ntasks=$SLURM_NNODES "$GPUCHK" check 2>/dev/null)

if [ -n "$BADGPU" ]; then
    echo "ERROR: allocated GPU(s) report uncorrected ECC errors:" >&2
    echo "$BADGPU" | sed 's/^/  /' >&2
    if [ "${SLURM_RESTART_COUNT:-0}" -lt 3 ]; then
        echo "Requeueing (attempt ${SLURM_RESTART_COUNT:-0} of 3) rather than burning the allocation." >&2
        scontrol requeue $SLURM_JOB_ID
        exit 1
    fi
    echo "Requeue limit reached -- running anyway; expect a CUDA init failure." >&2
fi

# --- input staging ---------------------------------------------------------
if [ "$LAUNCH_CONF" = "job0.conf" ]; then
    EQ_SRC=$SLURM_SUBMIT_DIR/equilibrate/$EQ_PROTEIN
    if [ ! -f "$EQ_SRC/equinpt.coor" ] || [ ! -f "$EQ_SRC/equinpt.xsc" ]; then
        echo "ERROR: missing equinpt.{coor,xsc} in $EQ_SRC" >&2
        exit 1
    fi
    cp "$EQ_SRC/equinpt.coor" "$EQ_SRC/equinpt.xsc" "$SLURM_SUBMIT_DIR/rest2/"
    for i in $(seq 0 $((NUM_REPLICAS - 1))); do
        mkdir -p $RUNDIR/$i
    done
    LOG_TAG="job0"
else
    if [ ! -d "$RUNDIR" ]; then
        echo "ERROR: cannot restart -- $RUNDIR not found" >&2
        exit 1
    fi
    LOG_TAG="restart-$(date +%Y%m%d_%H%M%S)"
fi

cd $SLURM_SUBMIT_DIR/rest2

# =============================================================================
# LAUNCH
# =============================================================================
# ++mpiexec        -- spawn via an external launcher instead of charmrun's ssh
# ++remote-shell   -- that launcher is srun, so Slurm places all 8 processes
#                     across both nodes and hands each its share of the GRES
# ++ppn 1          -- 1 PE per process; +p is the TOTAL PE count and must be
#                     divisible by +replicas (this is what broke job 56883716)
# +devicesperreplica 1 -- each replica binds one GPU. All 4 GPUs stay visible to
#                     every task on a node, and NAMD picks by physical rank
#                     within the host, so replicas 4-7 map to local devices 0-3.
echo "############ LAUNCH $(date -Is) ############"
set -x
$NAMD_HOME/charmrun $NAMD_HOME/namd3 ++mpiexec ++remote-shell srun \
    +p $NUM_REPLICAS ++ppn 1 \
    +replicas $NUM_REPLICAS +devicesperreplica 1 $LAUNCH_CONF \
    +stdout $RUNDIR/%d/$LOG_TAG.log
RC=$?
set +x
echo "############ charmrun exited rc=$RC $(date -Is) ############"

# =============================================================================
# POST-MORTEM -- only on failure. Replica logs are the only place NAMD reports
# per-replica GPU binding and CUDA errors; the slurm .err just gets the abort.
# =============================================================================
if [ $RC -ne 0 ]; then
    echo "--- per-replica log tails ---" >&2
    for i in $(seq 0 $((NUM_REPLICAS - 1))); do
        echo "===== replica $i ($RUNDIR/$i/$LOG_TAG.log) =====" >&2
        if [ -f "$RUNDIR/$i/$LOG_TAG.log" ]; then
            grep -E "binding to CUDA device|FATAL|ERROR|Charm\+\+> Running" \
                 "$RUNDIR/$i/$LOG_TAG.log" | head -10 >&2
            echo "  ...last 5 lines:" >&2
            tail -5 "$RUNDIR/$i/$LOG_TAG.log" | sed 's/^/  /' >&2
        else
            echo "  (no log written -- process died before NAMD startup)" >&2
        fi
    done
    echo "--- GPU state after failure ---" >&2
    srun --ntasks-per-node=1 --ntasks=$SLURM_NNODES "$GPUCHK" report >&2 2>&1
fi

rm -f "$GPUCHK"
exit $RC
