# Shared preparation for the Linux and macOS distribution build scripts.
# Sourced after each entry point parses its arguments.
repo=$(cd "$script_dir/../.." && pwd -P)

prepare_build() {
    system=$(uname -s)
    case "$system" in
        Linux) plat=${plat:-manylinux_2_28_$(uname -m)} ;;
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
}

one_wheel() {
    local wheels=("$1"/*.whl)
    if [[ ${#wheels[@]} != 1 || ! -f "${wheels[0]}" ]]; then
        echo "Expected exactly one wheel in $1" >&2
        return 1
    fi
    printf '%s\n' "${wheels[0]}"
}
