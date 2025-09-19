# low point
set_kv M_unit 1.83E27 "$PAR0"
rm -f "$H50"; /work/vmo703/aricarte/run_ipole.sh -par "$PAR0" >/dev/null
F0=$(readF "$H50"); M0=$(awk '$1=="M_unit"{print $2}' "$PAR0")
echo "low:  M0=$M0  F0=$F0"

# high point (known good for this setup)
set_kv M_unit 6.17572E28 "$PAR0"
rm -f "$H50"; /work/vmo703/aricarte/run_ipole.sh -par "$PAR0" >/dev/null
F1=$(readF "$H50"); M1=$(awk '$1=="M_unit"{print $2}' "$PAR0")
echo "high: M1=$M1  F1=$F1"

# one secant step to T=0.53 (optional polish)
M2=$(python - <<PY
M0=float("$M0"); F0=float("$F0")
M1=float("$M1"); F1=float("$F1"); T=0.53
print("{:.6E}".format(M1 + (T - F1)*(M1 - M0)/(F1 - F0)))
PY
)
set_kv M_unit "$M2" "$PAR0"
rm -f "$H50"; /work/vmo703/aricarte/run_ipole.sh -par "$PAR0" >/dev/null
python /work/vmo703/scratch/paper_metrics.py "$H50"
