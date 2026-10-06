#!/bin/bash
#SBATCH --job-name=p4_qa_watch
#SBATCH --partition=anantuabhg
#SBATCH --account=anantuabhg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH --output=/work/vmo703/igrmonty_logs/p4_qa_watch_%j.out
#SBATCH --error=/work/vmo703/igrmonty_logs/p4_qa_watch_%j.err

# Phase-4 QA watchdog (rescue 2026-09-25): audits every final spectrum as it
# lands in run_2026-09-15 by running the production QA harness, and logs any
# rescue-array task failures. Everything appends to p4_qa_live.log.
# Rescue arrays: 818690 (crit), 818691 (resume, superseded), 818692 (MAD wJET
# x16), 820422 (MAD CRITBETA fixed-bias), 822136 (resume redo x8, 2026-09-29).

RUN_DIR=/work/vmo703/igrmonty_outputs/m87/run_2026-09-15
QA=/work/vmo703/igrmonty_outputs/m87/_qa/p4_production_qa.py
LOG=/work/vmo703/igrmonty_outputs/m87/_qa/p4_qa_live.log
SEEN=/work/vmo703/igrmonty_outputs/m87/_qa/.p4_qa_seen
PY=/work/vmo703/ipole_venv/bin/python
ARRAYS="818690 818692 820422 824109 824666 825220"

touch "$SEEN"
echo "[$(date '+%F %T')] QA watchdog started (job ${SLURM_JOB_ID:-manual}) on ${HOSTNAME}" >> "$LOG"

while true; do
  new=0
  for f in "$RUN_DIR"/spectrum_*.h5; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    case "$base" in *_trial*) continue;; esac
    grep -qxF "$base" "$SEEN" && continue
    echo "[$(date '+%F %T')] NEW FINAL: $base" >> "$LOG"
    echo "$base" >> "$SEEN"
    new=1
  done
  if [ "$new" = "1" ]; then
    sleep 60   # let any concurrent writer finish
    echo "[$(date '+%F %T')] running QA harness over $RUN_DIR ..." >> "$LOG"
    "$PY" "$QA" "$RUN_DIR" >> "$LOG" 2>&1
    "$PY" "$(dirname "$QA")/p4_passfail.py" >> "$LOG" 2>&1
    echo "[$(date '+%F %T')] QA pass done (CSV + pass/fail table refreshed)" >> "$LOG"
  fi
  for j in $ARRAYS; do
    sacct -j "$j" -X -n --format=JobID%20,State%16 2>/dev/null \
      | grep -E 'FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL|CANCELLED' \
      | while read -r jid st _; do
          [ -z "$jid" ] && continue
          grep -qxF "FAIL $jid" "$SEEN" && continue
          echo "[$(date '+%F %T')] TASK TERMINAL-BAD: $jid $st" >> "$LOG"
          echo "FAIL $jid" >> "$SEEN"
        done
  done
  sleep 1200
done
