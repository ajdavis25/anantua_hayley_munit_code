#!/bin/bash
#SBATCH --job-name=ls_freeze_probe
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --output=/work/vmo703/scratch/ls_probe_778297/compute_probe-%j.out

# Diagnose the 778296/778297 freeze: rerun both frozen workloads on a compute1
# node with a liveness sidecar (per-minute CPU ticks, MemAvailable, checkpoint
# mtime, log sizes). Distinguishes "computing but slow" (utime climbs) from
# "blocked" (utime flat) independent of stdio buffering.
set -u
P=/work/vmo703/scratch/ls_probe_778297
R=/work/vmo703/light-scattering
W=$R/monte_carlo_cpp/results/measurement_case_reports/_batched_work
B=frozen_marseille_twilight_20220815_191413z_measurement_strict_subset_fixed__batch_0000_0000
cd "$R"
echo "probe_start node=$(hostname) utc=$(date -u +%FT%TZ)"
free -g | sed 's/^/mem: /'

mkdir -p "$P/cn_out"
cp "$W/${B}_checkpoint.txt" "$P/cn_checkpoint.txt"
sed -e "s|^output_dir=.*|output_dir=$P/cn_out|" "$W/${B}.cfg" > "$P/cn_probe.cfg"
env OMP_NUM_THREADS=1 stdbuf -oL -eL monte_carlo_cpp/build_fix/MeasurementCaseRunner \
    "$P/cn_probe.cfg" --checkpoint-state "$P/cn_checkpoint.txt" --higher-order-block-size 1 \
    > "$P/cn_probe.log" 2>&1 &
M_PID=$!

env OMP_NUM_THREADS=2 stdbuf -oL -eL monte_carlo_cpp/build_cluster_gcc8/ValidationRunner \
    monte_carlo_cpp/config/default_clear_sky.cfg > "$P/cn_default.log" 2>&1 &
V_PID=$!

for i in $(seq 1 345); do
  sleep 60
  line="tick=$i utc=$(date -u +%T)"
  for tag in "M:$M_PID" "V:$V_PID"; do
    pid=${tag#*:}
    if [ -d "/proc/$pid" ]; then
      line="$line ${tag%%:*}=[$(awk '{print $3":"$14}' /proc/$pid/stat 2>/dev/null)]"
    else
      line="$line ${tag%%:*}=[gone]"
    fi
  done
  line="$line memavail_kb=$(awk '/MemAvailable/{print $2}' /proc/meminfo) load=$(cut -d' ' -f1 /proc/loadavg)"
  line="$line ckpt=$(date -u -r "$P/cn_checkpoint.txt" +%T 2>/dev/null) mlog=$(wc -c < "$P/cn_probe.log") vlog=$(wc -c < "$P/cn_default.log")"
  echo "$line"
  if [ ! -d "/proc/$M_PID" ] && [ ! -d "/proc/$V_PID" ]; then break; fi
done
echo "probe_end utc=$(date -u +%FT%TZ)"
wait
