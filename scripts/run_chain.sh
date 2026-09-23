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
# Examples:
#   scripts/run_chain.sh lcdm_mnu 4
#   scripts/run_chain.sh lcdm_mnu 4 --resume
#   scripts/run_chain.sh w0wa_mnu 4

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
# Absolute interpreter path -- see CLAUDE.md's environment gotcha: ambient shell
# state on this machine can otherwise silently pick up the wrong Python/env.
PYTHON="/home/astro/.conda/envs/nuproj/bin/python"

if [[ ! -f "$CONFIG" ]]; then
    echo "No such config: $CONFIG" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

export PYTHONNOUSERSITE=1

LOG_FILE="$OUT_DIR/${NAME}.log"
echo "Launching '$NAME' with $NPROC MPI processes (resume=${RESUME_FLAG:+yes})"
echo "Cobaya output: $OUT_DIR/${NAME}.<n>.txt etc."
echo "Run log:       $LOG_FILE"

cd "$REPO_ROOT"
mpirun -np "$NPROC" "$PYTHON" -m cobaya run "$CONFIG" "${RESUME_FLAG[@]}" 2>&1 \
    | tee -a "$LOG_FILE"
