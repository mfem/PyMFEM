#!/usr/bin/env bash
# Build, repair and validate a serial distribution wheel using the active Python.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$script_dir/wheel_build_common.sh"

usage() {
    echo "Usage: bash $0 NEW_DIR [--plat LINUX_PLATFORM]"
}
[[ ${1:-} != -h && ${1:-} != --help ]] || { usage; exit 0; }
[[ $# -ge 1 ]] || { usage >&2; exit 1; }
work=$1
plat=
shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        --plat) plat=${2:?Missing platform}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done
prepare_build
copy_build_source
cd "$work/source"
set -x

pip freeze
pip wheel . --no-deps --verbose --wheel-dir "$work/raw"
cd "$work"
python "$script_dir/select_wheel_payload.py" "$(one_wheel "$work/raw")" \
    --profile serial --output-dir "$work/selected"
if [[ "$system" == Linux ]]; then
    auditwheel repair "$(one_wheel "$work/selected")" \
        --plat "$plat" --wheel-dir "$work/repaired"
else
    delocate-wheel --verbose --require-archs "$(uname -m)" --lib-sdir .dylibs \
        --sanitize-rpaths --wheel-dir "$work/repaired" "$(one_wheel "$work/selected")"
fi
repaired=$(one_wheel "$work/repaired")
python "$script_dir/select_wheel_payload.py" "$repaired" \
    --profile serial --normalize-rpaths --output-dir "$work"
final="$work/${repaired##*/}"
python "$script_dir/select_wheel_payload.py" "$final" --profile serial --check
twine check "$final"

python -m venv "$work/runtime"
# Run installed-wheel checks without build-time Python/library search paths.
unset PYTHONPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH
"$work/runtime/bin/pip" install --only-binary=:all: "$final"
"$work/runtime/bin/pip" check
"$work/runtime/bin/python" -I "$script_dir/check_serial_install.py" \
    "$work/source/examples/ex1.py" "$work/source/data/star.mesh"
set +x
echo "Validated wheel: $final"
