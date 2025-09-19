# set or insert: key value
set_kv () {
  k="$1"; v="$2"; f="$3"
  awk -v k="$k" -v v="$v" 'BEGIN{d=0}
    { if(!d && $1==k){print k, v; d=1} else print }
    END{ if(!d) print k, v }' "$f" > "$f.new" && mv "$f.new" "$f"
}

# reset a .par to the clean baseline for Fig 9 RBETA
# args: reset_par <par> <positron_ratio 0|1>
reset_par () {
  par="$1"; pairs="$2"

  # remove confusing/duplicate or unwanted lines
  sed -i -E '/^[[:space:]]*positronRatio[[:space:]]/d' "$par"   # kill camelCase duplicate
  sed -i -E '/^[[:space:]]*rmax_geo[[:space:]]/d' "$par"        # ensure no geo cap

  # enforce the baseline keys/values we want
  set_kv positron_ratio "$pairs" "$par"
  set_kv electronModel RBETA "$par"
  set_kv emission_type 4 "$par"
  set_kv trat_small 1 "$par"
  set_kv trat_large 20 "$par"             # <- back to 20 per PI
  set_kv beta_crit 1 "$par"
  set_kv beta_crit_coefficient 0.50 "$par"
  set_kv sigma_cut 2 "$par"
}

# one run -> read Ftot (Jy)
read_F () {
  python - <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], 'r') as H:
    print(float(H['/Ftot'][()]))
PY
}

# “round to 0.5 Jy” pre-step + N secant steps to final target
# args: secant_flux_tune <par> <h5> <F_target> <N_secant_steps> [F_round=0.5]
secant_flux_tune () {
  par="$1"; h5="$2"; target="$3"; n="$4"; FROUND="${5:-0.5}"

  # make sure we start from a fresh .h5
  rm -f "$h5"
  /work/vmo703/aricarte/run_ipole.sh -par "$par" >/dev/null
  F0=$(read_F "$h5"); M0=$(awk '$1=="M_unit"{print $2}' "$par")

  # pre-step: aim flux at ~0.5 Jy to “round I to 0.5”
  M1=$(python - <<PY
M0=float("$M0"); F0=float("$F0"); Fr=float("$FROUND")
print("{:.6E}".format(M0*Fr/F0))
PY
)
  set_kv M_unit "$M1" "$par"
  rm -f "$h5"
  /work/vmo703/aricarte/run_ipole.sh -par "$par" >/dev/null
  F1=$(read_F "$h5")

  # secant iterations (use (M0,F0) & (M1,F1))
  Mprev="$M0"; Fprev="$F0"; Mcurr="$M1"; Fcurr="$F1"
  for i in $(seq 1 $n); do
    Mnext=$(python - <<PY
M0=float("$Mprev"); F0=float("$Fprev"); M1=float("$Mcurr"); F1=float("$Fcurr"); T=float("$target")
# fallback to linear scaling if degenerate
if abs(F1-F0) < 1e-12:
    print("{:.6E}".format(M1*T/F1))
else:
    print("{:.6E}".format(M1 + (T - F1)*(M1 - M0)/(F1 - F0)))
PY
)
    set_kv M_unit "$Mnext" "$par"
    rm -f "$h5"
    /work/vmo703/aricarte/run_ipole.sh -par "$par" >/dev/null
    Fnext=$(read_F "$h5")
    echo "  step $i: F=${Fnext} Jy   M_unit=${Mnext}"
    Mprev="$Mcurr"; Fprev="$Fcurr"; Mcurr="$Mnext"; Fcurr="$Fnext"
  done

  # final quick report (uses your paper_metrics.py)
  python /work/vmo703/scratch/paper_metrics.py "$h5"
}
