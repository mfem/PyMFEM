#!/usr/bin/env bash
# Build the MPICH plus GSLIB add-on and validate it with a matching serial wheel.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/wheel_build_common.sh"

usage() {
    echo "Usage: bash $0 NEW_DIR BASE_WHEEL [--plat LINUX_PLATFORM]"
}
[[ ${1:-} != -h && ${1:-} != --help ]] || { usage; exit 0; }
[[ $# -ge 2 ]] || { usage >&2; exit 1; }
work=$1
base=$2
plat=
shift 2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --plat) plat=${2:?Missing platform}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done
[[ -n "$work" && -n "$base" ]] || { usage >&2; exit 1; }
base=$(python -c 'import sys; from pathlib import Path; print(Path(sys.argv[1]).resolve(strict=True))' "$base")
mpi_bin=$(python -c 'import sys; from pathlib import Path; print(Path(sys.executable).parent)')
export MPICC=${MPICC:-$mpi_bin/mpicc}
export MPICXX=${MPICXX:-$mpi_bin/mpicxx}
prepare_build
copy_build_source
cd "$work/source"
set -x

pip freeze
cp pyproject-mpich.toml pyproject.toml
pip wheel . --no-deps --verbose --wheel-dir "$work/raw" \
    -Cwith-parallel=Yes -Cwith-gslib=Yes
cd "$work"
python "$script_dir/select_wheel_payload.py" "$(one_wheel "$work/raw")" \
    --profile mpich --output-dir "$work/selected"
# Use the inventory recorded by this wheel's isolated build environment.
python - "$(one_wheel "$work/selected")" > "$work/mpi-libraries.txt" <<'PY'
import json
import sys
from zipfile import ZipFile

with ZipFile(sys.argv[1]) as wheel:
    config = json.loads(wheel.read('mfem/_par/mpi-runtime.json'))
print('\n'.join(config['libraries']))
PY
excludes=()
while IFS= read -r library; do
    excludes+=(--exclude "$library")
done < "$work/mpi-libraries.txt"
[[ ${#excludes[@]} -gt 0 ]] || { echo 'No supplied MPI libraries found.' >&2; exit 1; }
if [[ "$system" == Linux ]]; then
    auditwheel repair "$(one_wheel "$work/selected")" "${excludes[@]}" \
        --plat "$plat" --wheel-dir "$work/repaired"
else
    delocate-wheel --verbose --require-archs "$(uname -m)" --lib-sdir .dylibs \
        --sanitize-rpaths "${excludes[@]}" --wheel-dir "$work/repaired" \
        "$(one_wheel "$work/selected")"
fi
repaired=$(one_wheel "$work/repaired")
python "$script_dir/select_wheel_payload.py" "$repaired" \
    --profile mpich --normalize-rpaths --output-dir "$work"
final="$work/${repaired##*/}"
python "$script_dir/select_wheel_payload.py" "$final" --profile mpich --check
twine check "$final"
python "$script_dir/validate_mpich_pair.py" "$base" "$final" \
    --work-dir "$work/validation" --examples-dir "$work/source/examples" \
    --mesh "$work/source/data/star.mesh"
set +x
echo "Validated wheel: $final"
