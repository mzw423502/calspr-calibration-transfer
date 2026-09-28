#!/usr/bin/env bash
set -uo pipefail

if [ "$#" -ne 3 ]; then
  echo "usage: $0 RUN_NAME CONFIG SEED" >&2
  exit 64
fi

run_name="$1"
config="$2"
seed="$3"
root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/logs"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="$root/logs/${run_name}_${stamp}.log"
status_file="$root/logs/${run_name}_${stamp}.status"

if [ -e "$log" ] || [ -e "$status_file" ]; then
  echo "refusing to overwrite existing run log" >&2
  exit 73
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED="$seed"
runtime_vendor="${CALSPR_RUNTIME_VENDOR:-}"
if [ -n "$runtime_vendor" ]; then
  export PYTHONPATH="$runtime_vendor${PYTHONPATH:+:$PYTHONPATH}"
fi

if [ "$run_name" = "tests" ]; then
  program="$root/tests/run_tests.py"
  args=()
  input_description="none"
  input_hash="not_applicable"
else
  program="$root/scripts/run_${run_name}.py"
  args=("$root/$config")
  input_description="$root/$config"
  input_hash="$(sha256sum "$root/$config" | awk '{print $1}')"
fi

start_epoch="$(date +%s)"
{
  echo "run_name=$run_name"
  echo "host=$(hostname)"
  echo "user=$(whoami)"
  echo "project_root=$root"
  echo "start_utc=$(date -u --iso-8601=seconds)"
  echo "command=/usr/bin/time -v python3 $program ${args[*]}"
  echo "input=$input_description"
  echo "input_sha256=$input_hash"
  echo "random_seed=$seed"
  echo "runtime_vendor=$runtime_vendor"
  echo "python=$(python3 --version 2>&1)"
  echo "numpy=$(python3 -c 'import numpy; print(numpy.__version__)' 2>&1)"
  if [ "$run_name" = "figure1" ]; then
    echo "matplotlib=$(python3 -c 'import matplotlib; print(matplotlib.__version__)' 2>&1 || echo unavailable)"
  else
    echo "matplotlib=not_requested_for_numerical_run"
  fi
  echo "git=$(git --version 2>&1)"
  echo "threads=OMP:1,OPENBLAS:1,MKL:1,NUMEXPR:1"
  echo "memory_before_begin"
  free -h
  echo "memory_before_end"
} > "$log"

/usr/bin/time -v python3 "$program" "${args[@]}" >> "$log" 2>&1
code=$?
end_epoch="$(date +%s)"
{
  echo "end_utc=$(date -u --iso-8601=seconds)"
  echo "elapsed_wall_seconds=$((end_epoch-start_epoch))"
  echo "exit_status=$code"
  echo "output_files_begin"
  find "$root/outputs" -maxdepth 5 -type f -printf '%p\n' | sort
  echo "output_files_end"
  echo "memory_after_begin"
  free -h
  echo "memory_after_end"
} >> "$log"
printf 'exit_status=%s\nlog=%s\n' "$code" "$log" > "$status_file"
exit "$code"
