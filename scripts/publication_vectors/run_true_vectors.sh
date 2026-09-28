#!/bin/sh
# Run from any working directory. Supply an absent or empty output directory.
set -eu
if test "$#" -ne 1; then
  echo 'Usage: sh scripts/publication_vectors/run_true_vectors.sh /absolute/new/output_dir' >&2
  exit 64
fi
script_dir=$(dirname "$(realpath "$0")")
vector_python=${CALSPR_VECTOR_PYTHON:-python3}
out=$("$vector_python" -c 'import pathlib,sys; root=pathlib.Path(sys.argv[1]).resolve().parents[1]; out=pathlib.Path(sys.argv[2]).expanduser().resolve(); assert out != root and root not in out.parents, "Output must be outside the frozen repository tree"; assert not out.exists() or (out.is_dir() and not any(out.iterdir())), "Output must be absent or empty"; print(out)' "$script_dir" "$1")
log="$out.run.log"
status="$out.exit_status"
if test -e "$log" || test -e "$status"; then
  echo "Refusing to overwrite existing run log or status: $log" >&2
  exit 73
fi
mkdir -p "$(dirname "$out")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
{
  date -u '+START_UTC=%Y-%m-%dT%H:%M:%SZ'
  echo "OUTPUT_DIR=$out"
  echo "LOG_ABSOLUTE_PATH=$log"
  echo 'RANDOM_SEED=NA; deterministic display regeneration; one CPU thread'
  "$vector_python" --version
  echo "COMMAND=$vector_python $script_dir/make_true_vectors.py --output-dir $out"
  set +e
  "$vector_python" "$script_dir/make_true_vectors.py" --output-dir "$out"
  code=$?
  if test "$code" -eq 0; then
    echo "COMMAND=$vector_python $script_dir/audit_true_vectors.py $out"
    "$vector_python" "$script_dir/audit_true_vectors.py" "$out"
    code=$?
  fi
  set -e
  echo "EXIT_STATUS=$code"
  date -u '+END_UTC=%Y-%m-%dT%H:%M:%SZ'
  echo "$code" > "$status"
  exit "$code"
} > "$log" 2>&1
