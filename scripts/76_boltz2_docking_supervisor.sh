#!/bin/bash
### Run on nebula, within the boltz conda environment, via nohup.
# Script 73 makes exactly one pass over all pockets per invocation, by design, and exits --
# pockets that fail (e.g. a transient MSA-server rate-limit) or were never attempted just get
# skipped, not retried, until something re-invokes the script. This supervisor keeps re-invoking
# it so the multi-day run doesn't silently stall after one pass: already-done pockets are cheap
# to skip (boltz's own resumability + our MSA cache), so repeated passes mostly cost time only on
# pockets that still need work. Kill this process (not just script 73's) to actually stop the run.
#
# Usage (on nebula):
#   conda activate boltz
#   nohup setsid bash scripts/76_boltz2_docking_supervisor.sh >> output/73_boltz2_docking/full_run.log 2>&1 &

cd "$(dirname "$0")/.."

SLEEP_BETWEEN_PASSES_S=300

while true; do
  echo "=== Supervisor: pass started $(date) ==="
  python scripts/73_boltz2_docking.py
  echo "=== Supervisor: pass finished $(date), sleeping ${SLEEP_BETWEEN_PASSES_S}s ==="
  sleep "$SLEEP_BETWEEN_PASSES_S"
done
