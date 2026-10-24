#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
CONFIG_FILE="${PYGAMD_INSTALL_CONFIG:-$ROOT/install.txt}"

TARGET_OS=""
GPU_VALUE=""
GPU_MODEL=""
GPU_GENERATION=0
OPEN_MPI="false"
RUN_EXAMINE=""
ENV_NAME=""
CUDA_FLAVOR=""
CUDA_EXTRA=""
NATIVE_GALA=""
NATIVE_MOLGEN=""
STAGED_GALA=""
STAGED_MOLGEN=""
ENV_PREFIX=""
RUNTIME_LD=""
VALIDATION_DIR=""
CLEANED=0
CLEAN_ONLY=0
RECREATE_ENV=0
FORCE_REINSTALL=0

info() {
    printf '[pygamd-install] %s\n' "$*"
}

warn() {
    printf '[pygamd-install] WARNING: %s\n' "$*" >&2
}

die() {
    printf '[pygamd-install] ERROR: %s\n' "$*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage: bash install.sh [OPTION]

Install PYGAMD using the settings in install.txt.

Options:
  --clean-only       Remove generated repository artifacts without installing
  --recreate-env     Remove the configured Conda environment, recreate it, and install
  --force-reinstall  Reinstall PYGAMD and all resolved dependencies in the environment
  -h, --help         Show this help message

With no option, an existing environment is reused or a missing environment is
created, then PYGAMD is installed or upgraded.
EOF
}

parse_args() {
    while (( $# )); do
        case "$1" in
            --clean-only)
                CLEAN_ONLY=1
                ;;
            --recreate-env)
                RECREATE_ENV=1
                ;;
            --force-reinstall)
                FORCE_REINSTALL=1
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                die "unknown option: $1 (use --help for usage)"
                ;;
        esac
        shift
    done

    if (( CLEAN_ONLY + RECREATE_ENV + FORCE_REINSTALL > 1 )); then
        die '--clean-only, --recreate-env, and --force-reinstall are mutually exclusive'
    fi
}

read_global_value() {
    local key="$1"
    awk -v wanted="$key" '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        {
            sub(/\r$/, "")
            line = $0
            sub(/[[:space:]]*#.*/, "", line)
            line = trim(line)
            if (line == "") next
            if (line == "{") { depth++; next }
            if (line == "}") { if (depth > 0) depth--; next }
            if (depth != 0) next
            separator = index(line, "=")
            if (!separator) next
            name = trim(substr(line, 1, separator - 1))
            if (name == wanted) {
                print trim(substr(line, separator + 1))
                exit
            }
        }
    ' "$CONFIG_FILE"
}

read_section_value() {
    local section="$1"
    local key="$2"
    awk -v wanted_section="$section" -v wanted_key="$key" '
        function trim(value) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            return value
        }
        {
            sub(/\r$/, "")
            line = $0
            sub(/[[:space:]]*#.*/, "", line)
            line = trim(line)
            if (line == "") next
            if (line == "{") {
                inside = pending
                pending = 0
                next
            }
            if (line == "}") {
                inside = 0
                next
            }
            separator = index(line, "=")
            if (!separator) next
            name = trim(substr(line, 1, separator - 1))
            value = trim(substr(line, separator + 1))
            if (!inside) {
                pending = (name == wanted_section)
                next
            }
            if (name == wanted_key) {
                print value
                exit
            }
        }
    ' "$CONFIG_FILE"
}

parse_boolean() {
    local value="${1,,}"
    local label="$2"
    case "$value" in
        true|false)
            printf '%s\n' "$value"
            ;;
        *)
            die "$label must be true or false in install.txt"
            ;;
    esac
}

resolve_config() {
    local linux_enabled windows_enabled

    [[ -f "$CONFIG_FILE" ]] || die "configuration file not found: $CONFIG_FILE"

    linux_enabled="$(parse_boolean "$(read_global_value linux)" linux)"
    windows_enabled="$(parse_boolean "$(read_global_value windows)" windows)"
    RUN_EXAMINE="$(parse_boolean "$(read_global_value run-examine)" run-examine)"
    ENV_NAME="$(read_global_value conda-python-name)"
    [[ -n "$ENV_NAME" ]] || ENV_NAME='pygamd-python312'
    [[ "$ENV_NAME" =~ ^[A-Za-z0-9._-]+$ ]] \
        || die "invalid Conda environment name: $ENV_NAME"

    if [[ "$linux_enabled" == "$windows_enabled" ]]; then
        die 'exactly one of linux or windows must be true in install.txt'
    fi
    if [[ "$linux_enabled" == true ]]; then
        TARGET_OS='linux'
        OPEN_MPI="$(parse_boolean "$(read_section_value linux open-mpi)" open-mpi)"
    else
        TARGET_OS='windows'
        OPEN_MPI='false'
    fi

    GPU_VALUE="$(read_section_value "$TARGET_OS" gpu)"
    [[ -n "$GPU_VALUE" ]] || die "the $TARGET_OS section must define gpu=..."
    GPU_MODEL="$(printf '%s\n' "$GPU_VALUE" | sed -nE 's/^[^0-9]*([0-9]{4}).*/\1/p')"
    [[ "$GPU_MODEL" =~ ^[0-9]{4}$ ]] \
        || die "unable to read a four-digit NVIDIA GPU model from: $GPU_VALUE"
    GPU_GENERATION=$((10#$GPU_MODEL / 1000))

    if (( GPU_GENERATION >= 5 )); then
        CUDA_FLAVOR='cuda13'
        CUDA_EXTRA='cuda13'
    else
        CUDA_FLAVOR='cuda12'
        CUDA_EXTRA='cuda12'
    fi

    if [[ "$TARGET_OS" == linux ]]; then
        NATIVE_MOLGEN="$ROOT/install/linux_python312/molgen.so"
        STAGED_MOLGEN="$ROOT/poetry/molgen.so"
        STAGED_GALA="$ROOT/poetry/cu_gala.so"
        if [[ "$CUDA_FLAVOR" == cuda13 ]]; then
            if [[ "$OPEN_MPI" == true ]]; then
                NATIVE_GALA="$ROOT/install/linux_cuda-13_python312_mpi_cant_sm_61/cu_gala.so"
            else
                NATIVE_GALA="$ROOT/install/linux_cuda-13_python312_cant_sm_61/cu_gala.so"
            fi
        elif [[ "$OPEN_MPI" == true ]]; then
            NATIVE_GALA="$ROOT/install/linux_cuda-12-2_python312_mpi_cant_sm_120/cu_gala.so"
        else
            NATIVE_GALA="$ROOT/install/linux_cuda-12-2_python312_cant_sm_120/cu_gala.so"
        fi
    else
        NATIVE_MOLGEN="$ROOT/install/win_python312/molgen.dll"
        STAGED_MOLGEN="$ROOT/poetry/molgen.pyd"
        STAGED_GALA="$ROOT/poetry/cu_gala.pyd"
        if [[ "$CUDA_FLAVOR" == cuda13 ]]; then
            NATIVE_GALA="$ROOT/install/win_cuda-13_python312_cant_sm_61/cu_gala.dll"
            [[ -f "$NATIVE_GALA" ]] \
                || die 'Windows 50-series support is reserved but its cu_gala.dll has not been compiled yet'
        else
            NATIVE_GALA="$ROOT/install/win_cuda-12-2_python312_cant_sm_120/cu_gala.dll"
        fi
    fi

    [[ -f "$NATIVE_MOLGEN" ]] || die "missing native library: $NATIVE_MOLGEN"
    [[ -f "$NATIVE_GALA" ]] || die "missing native library: $NATIVE_GALA"
}

cleanup_repo_artifacts() {
    if (( CLEANED )); then
        return 0
    fi
    CLEANED=1
    find "$ROOT" -type d \( -name '.pytest_cache' -o -name '__pycache__' \) \
        -prune -exec rm -rf -- {} + 2>/dev/null || true
    find "$ROOT" -type f \( -name '*.pyc' -o -name '*.pyo' \) \
        -delete 2>/dev/null || true
    find "$ROOT" -maxdepth 2 -type d \
        \( -name 'build' -o -name 'dist' -o -name '*.egg-info' \) \
        -prune -exec rm -rf -- {} + 2>/dev/null || true
}

cleanup_on_exit() {
    local status=$?
    trap - EXIT
    cleanup_repo_artifacts
    if [[ -n "$VALIDATION_DIR" && -d "$VALIDATION_DIR" ]]; then
        rm -rf -- "$VALIDATION_DIR"
    fi
    exit "$status"
}

env_exists() {
    "$CONDA_EXE" env list 2>/dev/null \
        | awk -v target="$ENV_NAME" '$1 == target { found = 1 } END { exit !found }'
}

run_python() {
    if [[ "$TARGET_OS" == linux && -n "$RUNTIME_LD" ]]; then
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONNOUSERSITE=1 \
        LD_LIBRARY_PATH="$RUNTIME_LD${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
            "$CONDA_EXE" run --no-capture-output -n "$ENV_NAME" python "$@"
    else
        PYTHONDONTWRITEBYTECODE=1 \
        PYTHONNOUSERSITE=1 \
            "$CONDA_EXE" run --no-capture-output -n "$ENV_NAME" python "$@"
    fi
}

stage_native_libraries() {
    rm -f -- \
        "$ROOT/poetry/molgen.so" \
        "$ROOT/poetry/cu_gala.so" \
        "$ROOT/poetry/molgen.pyd" \
        "$ROOT/poetry/cu_gala.pyd"
    cp -- "$NATIVE_MOLGEN" "$STAGED_MOLGEN"
    cp -- "$NATIVE_GALA" "$STAGED_GALA"
    chmod 0644 "$STAGED_MOLGEN" "$STAGED_GALA" 2>/dev/null || true
}

write_linux_runtime_hooks() {
    local activate_dir="$ENV_PREFIX/etc/conda/activate.d"
    local deactivate_dir="$ENV_PREFIX/etc/conda/deactivate.d"
    local activate_hook="$activate_dir/pygamd-runtime.sh"
    local deactivate_hook="$deactivate_dir/pygamd-runtime.sh"
    local runtime_q flavor_q mpi_q

    runtime_q="$(printf '%q' "$RUNTIME_LD")"
    flavor_q="$(printf '%q' "$CUDA_FLAVOR")"
    mpi_q="$(printf '%q' "$OPEN_MPI")"
    mkdir -p "$activate_dir" "$deactivate_dir"
    {
        printf '%s\n' '# Generated by install.sh.'
        printf '%s\n' 'if [[ -z "${_PYGAMD_LD_LIBRARY_PATH_ACTIVE+x}" ]]; then'
        printf '%s\n' '    export _PYGAMD_LD_LIBRARY_PATH_ACTIVE=1'
        printf '%s\n' '    export _PYGAMD_LD_LIBRARY_PATH_WAS_SET="${LD_LIBRARY_PATH+x}"'
        printf '%s\n' '    export _PYGAMD_LD_LIBRARY_PATH_SAVED="${LD_LIBRARY_PATH-}"'
        printf '%s\n' 'fi'
        printf 'export PYGAMD_CUDA_FLAVOR=%s\n' "$flavor_q"
        printf 'export PYGAMD_OPEN_MPI=%s\n' "$mpi_q"
        printf 'export LD_LIBRARY_PATH=%s${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}\n' "$runtime_q"
    } > "$activate_hook"
    {
        printf '%s\n' '# Generated by install.sh.'
        printf '%s\n' 'if [[ "${_PYGAMD_LD_LIBRARY_PATH_ACTIVE:-0}" == 1 ]]; then'
        printf '%s\n' '    if [[ "${_PYGAMD_LD_LIBRARY_PATH_WAS_SET:-0}" == 1 ]]; then'
        printf '%s\n' '        export LD_LIBRARY_PATH="${_PYGAMD_LD_LIBRARY_PATH_SAVED}"'
        printf '%s\n' '    else'
        printf '%s\n' '        unset LD_LIBRARY_PATH'
        printf '%s\n' '    fi'
        printf '%s\n' '    unset _PYGAMD_LD_LIBRARY_PATH_ACTIVE'
        printf '%s\n' '    unset _PYGAMD_LD_LIBRARY_PATH_WAS_SET'
        printf '%s\n' '    unset _PYGAMD_LD_LIBRARY_PATH_SAVED'
        printf '%s\n' 'fi'
        printf '%s\n' 'unset PYGAMD_CUDA_FLAVOR'
        printf '%s\n' 'unset PYGAMD_OPEN_MPI'
    } > "$deactivate_hook"
    chmod 0644 "$activate_hook" "$deactivate_hook"
}

main() {
    local env_output install_extras validation_status=0
    # Prefer compatible binary wheels over newer source distributions.
    # This keeps installations reliable on older Linux systems (for example
    # glibc 2.17) while still selecting the newest compatible wheels on newer
    # Linux and Windows hosts.
    local -a pip_install_args=(--upgrade --no-build-isolation --prefer-binary)

    cd "$ROOT"
    parse_args "$@"
    if (( CLEAN_ONLY )); then
        info 'Removing generated repository artifacts; installation is disabled'
        cleanup_repo_artifacts
        info 'Repository cleanup completed successfully'
        return 0
    fi

    resolve_config
    trap cleanup_on_exit EXIT

    if [[ -n "${CONDA_EXE:-}" ]]; then
        CONDA_EXE="$CONDA_EXE"
    else
        CONDA_EXE="$(command -v conda || true)"
    fi
    [[ -n "$CONDA_EXE" ]] || die 'conda was not found in PATH'

    info "Target: $TARGET_OS"
    info "GPU setting: $GPU_VALUE"
    info "CUDA flavor: $CUDA_FLAVOR"
    [[ "$TARGET_OS" == linux ]] && info "OpenMPI binary: $OPEN_MPI"
    info "Conda environment: $ENV_NAME"
    info "Post-install examination: $RUN_EXAMINE"

    if (( RECREATE_ENV )); then
        [[ "$ENV_NAME" != base ]] || die 'refusing to remove the Conda base environment'
        if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV_NAME" ]]; then
            die "deactivate Conda environment $ENV_NAME before using --recreate-env"
        fi
        if env_exists; then
            info "Removing Conda environment: $ENV_NAME"
            "$CONDA_EXE" env remove -n "$ENV_NAME" -y
        fi
        info "Creating Conda environment: $ENV_NAME (Python 3.12)"
        "$CONDA_EXE" create -n "$ENV_NAME" python=3.12 pip -y
    elif env_exists; then
        info "Reusing existing Conda environment: $ENV_NAME"
    else
        info "Creating Conda environment: $ENV_NAME (Python 3.12)"
        "$CONDA_EXE" create -n "$ENV_NAME" python=3.12 pip -y
    fi

    if ! run_python -c \
        'import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 12) else 1)'; then
        die "environment $ENV_NAME is not Python 3.12"
    fi

    if [[ "$TARGET_OS" == linux ]]; then
        env_output="$(run_python -c 'import sys; print(sys.prefix)')" \
            || die "unable to resolve Conda prefix for $ENV_NAME"
        ENV_PREFIX="$(printf '%s\n' "$env_output" | tail -n 1)"
        [[ -d "$ENV_PREFIX" ]] || die "invalid Conda prefix: $ENV_PREFIX"
        RUNTIME_LD="$ENV_PREFIX/lib"
        write_linux_runtime_hooks
    fi

    info "Staging $(basename -- "$NATIVE_MOLGEN") and $(basename -- "$NATIVE_GALA") into poetry/"
    stage_native_libraries

    if ! run_python -c \
        'import setuptools; raise SystemExit(0 if int(setuptools.__version__.split(".", 1)[0]) >= 77 else 1)'; then
        info 'Updating build tooling'
        run_python -m pip install --upgrade 'setuptools>=77' wheel
    fi

    install_extras="$CUDA_EXTRA,aa2cg,rdf"
    [[ "$RUN_EXAMINE" == true ]] && install_extras+=',test'
    if (( FORCE_REINSTALL )); then
        info 'Force-reinstalling PYGAMD and all resolved dependencies'
        pip_install_args+=(--force-reinstall --upgrade-strategy eager)
    else
        info 'Installing PYGAMD'
    fi
    run_python -m pip install "${pip_install_args[@]}" ".[${install_extras}]"

    if [[ "$RUN_EXAMINE" == false ]]; then
        info 'Installation completed successfully; post-install examination was disabled'
        info "Activate with: conda activate $ENV_NAME"
        return 0
    fi

    VALIDATION_DIR="$(mktemp -d -t pygamd-validation.XXXXXX)"
    info 'Checking installed dependencies and native modules'
    if ! run_python -m pip check; then
        warn 'pip dependency check failed'
        validation_status=1
    fi
    if ! (cd "$VALIDATION_DIR" && run_python -c \
        'from poetry import molgen, cu_gala; print("native poetry modules: ok")'); then
        warn 'native module import check failed'
        validation_status=1
    fi

    info 'Running pytest regression suite'
    if ! PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 run_python -m pytest -q test \
        -p no:cacheprovider --basetemp="$VALIDATION_DIR/pytest"; then
        warn 'pytest regression suite failed'
        validation_status=1
    fi

    info 'Running bundled examples'
    if ! run_python scripts/run_examples.py \
        --gpu 0 \
        --gpu-steps "${PYGAMD_GPU_STEPS:-5}" \
        --output-dir "$VALIDATION_DIR/examples"; then
        warn 'bundled examples failed'
        validation_status=1
    fi

    info 'Removing validation outputs and generated repository artifacts'
    cleanup_repo_artifacts
    rm -rf -- "$VALIDATION_DIR"
    VALIDATION_DIR=""
    (( validation_status == 0 )) || die 'installation completed, but examination failed'

    info 'Installation and examination completed successfully'
    info "Activate with: conda activate $ENV_NAME"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
