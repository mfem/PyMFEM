# Shared preparation for the Linux and macOS distribution build scripts.
# Sourced after each entry point parses its arguments.
repo=$(cd "$script_dir/../.." && pwd -P)

copy_build_source() {
    # Preserve source edits while forcing native dependencies and wrappers to rebuild.
    python - "$repo" "$work/source" <<'PY'
import shutil
import sys

shutil.copytree(sys.argv[1], sys.argv[2], ignore=shutil.ignore_patterns(
    '.git', 'build', 'dist', 'external', '*.egg-info', '__pycache__',
    '*.pyc', '*.so', '*.so.*', '*.dylib', '*.o', '*.a',
    '*_wrap.cxx', '*_wrap.h', 'setup_local.py', 'mpi-runtime.json'))
PY
}

prepare_build() {
    system=$(uname -s)
    case "$system" in
        Linux)
            if [[ -z "$plat" ]]; then
                # Newer hosts may introduce symbols beyond the release baseline.
                # Auditwheel determines the actual wheel compatibility in auto mode.
                plat=$(python - <<'PY'
import platform

libc, version = platform.libc_ver()
if libc == 'glibc' and tuple(map(int, version.split('.')[:2])) <= (2, 28):
    print('manylinux_2_28_' + platform.machine())
else:
    print('auto')
PY
)
            fi
            ;;
        Darwin)
            export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-14.0}
            export CMAKE_OSX_DEPLOYMENT_TARGET=$MACOSX_DEPLOYMENT_TARGET
            export CMAKE_OSX_ARCHITECTURES=$(uname -m)
            export ARCHFLAGS="-arch $CMAKE_OSX_ARCHITECTURES"
            export _PYTHON_HOST_PLATFORM="macosx-$MACOSX_DEPLOYMENT_TARGET-$CMAKE_OSX_ARCHITECTURES"
            ;;
        *) echo 'These build scripts support Linux and macOS.' >&2; exit 1 ;;
    esac
    [[ -n "$work" ]] || { usage >&2; exit 1; }
    work=$(python -c 'import sys; from pathlib import Path; print(Path(sys.argv[1]).resolve())' "$work")
    [[ "$work/" != "$repo/"* ]] || { echo 'Use a work directory outside the checkout to avoid copying it into itself.' >&2; exit 1; }
    [[ ! -e "$work" ]] || { echo "Use a fresh work directory: $work" >&2; exit 1; }
    mkdir -p -- "$(dirname "$work")"
    mkdir -- "$work"
    exec > >(tee -a "$work/build.log") 2>&1
    trap 'echo "Build failed; see $work/build.log" >&2' ERR
    export PYMFEM_BUILD_NCPU=${PYMFEM_BUILD_NCPU:-2}
    echo "Build transcript: $work/build.log"
    if [[ "$system" == Linux ]]; then
        echo "Linux repair target: $plat (override with --plat PLATFORM)"
    fi
}

one_wheel() {
    local wheels=("$1"/*.whl)
    if [[ ${#wheels[@]} != 1 || ! -f "${wheels[0]}" ]]; then
        echo "Expected exactly one wheel in $1" >&2
        return 1
    fi
    printf '%s\n' "${wheels[0]}"
}
