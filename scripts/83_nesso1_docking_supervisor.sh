#!/bin/bash
### Run on nebula, within the nesso conda environment, via nohup.
# Script 80 makes exactly one pass over all structures per invocation, by design, and exits --
# structures that fail (e.g. the dimer OOM'ing even with the trimmed fallback) or were never
# attempted just get skipped, not retried, until something re-invokes the script. This
# supervisor keeps re-invoking it so a long run doesn't silently stall after one pass:
# already-done structures are cheap to skip (Nesso's own resumability), so repeated passes
# mostly cost time only on structures that still need work. Kill this process (not just script
# 80's) to actually stop the run.
#
# Usage (on nebula):
#   conda activate nesso
#   nohup setsid bash scripts/83_nesso1_docking_supervisor.sh >> output/80_nesso1_docking/full_run.log 2>&1 &

cd "$(dirname "$0")/.."

SLEEP_BETWEEN_PASSES_S=60

while true; do
  echo "=== Supervisor: pass started $(date) ==="
  python scripts/80_nesso1_docking.py
  echo "=== Supervisor: pass finished $(date), sleeping ${SLEEP_BETWEEN_PASSES_S}s ==="
  sleep "$SLEEP_BETWEEN_PASSES_S"
done
