#!/usr/bin/env bash
# Launch a Cobaya MCMC run defined by configs/<name>.yaml, under MPI, logging to
# chains/<name>/<name>.log.
#
# Per CLAUDE.md: every run must come from a configs/ YAML (no ad hoc runs), and any
# run expected to take more than 10 minutes must be launched by the user themselves
# (e.g. in tmux), not started automatically by an agent. This script is what gets
# printed as that launch command -- it is not invoked on its own here.
#
# Usage:
#   scripts/run_chain.sh <config-name> <n-mpi-processes> [--resume]
#
# OMP_NUM_THREADS (CAMB's OpenMP thread count per MPI process) defaults to 1 -- see
# the comment below for why. Override by exporting it before calling this script.
#
# Examples:
#   scripts/run_chain.sh lcdm_mnu 4
#   scripts/run_chain.sh lcdm_mnu 4 --resume
#   scripts/run_chain.sh w0wa_mnu 4
#   OMP_NUM_THREADS=2 scripts/run_chain.sh lcdm_mnu 4   # only if cores are free

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <config-name> <n-mpi-processes> [--resume]" >&2
    exit 1
fi

NAME="$1"
NPROC="$2"
shift 2

RESUME_FLAG=()
if [[ "${1:-}" == "--resume" ]]; then
    RESUME_FLAG=(--resume)
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$REPO_ROOT/configs/${NAME}.yaml"
OUT_DIR="$REPO_ROOT/chains/${NAME}"
# Deliberately *not* hardcoded to any machine-specific env path (e.g.
# /home/astro/.conda/envs/nuproj/bin/python) -- this script needs to work on
# whatever machine you actually run long chains on, which may not be the same one
# this project was scaffolded on. Instead, this requires the `nuproj` conda/mamba
# env to already be *activated* in the calling shell (`conda activate nuproj`),
# and just uses `python`/`mpirun` off PATH. Checked explicitly below so a missing
# activation fails loudly here, in the foreground, rather than causing `mpirun`/
# cobaya to silently fail deep inside an MPI subprocess where tmux -d would hide it.
PYTHON="python"
MPIRUN="mpirun"

if ! "$PYTHON" -c "import cobaya" >/dev/null 2>&1; then
    echo "ERROR: 'python -c \"import cobaya\"' failed. Activate the nuproj env" \
         "first, e.g.: conda activate nuproj (or: mamba activate nuproj)" >&2
    exit 1
fi
if ! command -v "$MPIRUN" >/dev/null 2>&1; then
    echo "ERROR: 'mpirun' not found on PATH. It should come bundled with the" \
         "nuproj env's mpi4py conda-forge install -- activate nuproj first, or" \
         "rebuild the env with: mamba env create -f environment.yml" >&2
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "No such config: $CONFIG" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

export PYTHONNOUSERSITE=1

# CAMB's Fortran backend is OpenMP-parallelized and, left unset, OMP_NUM_THREADS
# defaults to using every visible core *per MPI process* -- with NPROC MPI ranks
# each also spawning that many OpenMP threads, you get NPROC x nproc-way
# oversubscription on this machine, which is slower than not threading at all.
# Default here to 1 (pure MPI: one core per chain) since CAMB's own OpenMP scaling
# per evaluation is limited and cobaya's R-1 diagnostic benefits more from extra
# independent chains than from speeding up any single one. Override by exporting
# OMP_NUM_THREADS before calling this script if you have idle cores to spare and
# know what you're doing (e.g. running only one config alone).
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

AVAILABLE_CORES="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)"
REQUESTED_CORES=$((NPROC * OMP_NUM_THREADS))
if (( REQUESTED_CORES > AVAILABLE_CORES )); then
    echo "WARNING: requesting $NPROC MPI processes x $OMP_NUM_THREADS OpenMP" \
         "threads = $REQUESTED_CORES cores, but only $AVAILABLE_CORES are visible" \
         "on this machine. This check only sees this one invocation -- if you're" \
         "launching another config's chain at the same time, account for both" \
         "here manually." >&2
fi

LOG_FILE="$OUT_DIR/${NAME}.log"
echo "Launching '$NAME' with $NPROC MPI processes x OMP_NUM_THREADS=$OMP_NUM_THREADS (resume=${RESUME_FLAG:+yes})"
echo "Cobaya output: $OUT_DIR/${NAME}.<n>.txt etc."
echo "Run log:       $LOG_FILE"

cd "$REPO_ROOT"
# No -x/-genv flags: this is always a local (single-machine) run, and both Open MPI
# and MPICH inherit the launching shell's environment for local ranks by default, so
# OMP_NUM_THREADS/PYTHONNOUSERSITE (already exported above) reach every rank without
# needing implementation-specific propagation flags.
"$MPIRUN" -np "$NPROC" \
    "$PYTHON" -m cobaya run "$CONFIG" "${RESUME_FLAG[@]}" 2>&1 \
    | tee -a "$LOG_FILE"
