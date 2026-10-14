#!/usr/bin/env bash
# Eight-library build environment. P_SPRAY_ROOT keeps existing Make/options compatible.
# Source from Bash. Override discovery with OPENFOAM14_BASHRC if required.
export P_SPRAY_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${OPENFOAM14_BASHRC:-}" ]]; then
    p_spray_bashrc="$OPENFOAM14_BASHRC"
elif [[ "${WM_PROJECT_VERSION:-}" == "14" && -f "${WM_PROJECT_DIR:-}/etc/bashrc" ]]; then
    p_spray_bashrc="$WM_PROJECT_DIR/etc/bashrc"
else
    p_spray_bashrc=""
    for p_spray_candidate in "$HOME/OpenFOAM/OpenFOAM-14/etc/bashrc" /opt/openfoam14/etc/bashrc /usr/lib/openfoam/openfoam14/etc/bashrc; do
        if [[ -f "$p_spray_candidate" ]]; then
            p_spray_bashrc="$p_spray_candidate"
            break
        fi
    done
fi
if [[ ! -f "$p_spray_bashrc" ]]; then
    echo "Set OPENFOAM14_BASHRC to the Foundation OpenFOAM-14 etc/bashrc." >&2
    return 1 2>/dev/null || exit 1
fi
source "$p_spray_bashrc" || { return 1 2>/dev/null || exit 1; }
if [[ "${WM_PROJECT_VERSION:-}" != "14" ]]; then
    echo "This case requires OpenFOAM Foundation 14." >&2
    return 1 2>/dev/null || exit 1
fi
export FOAM_USER_LIBBIN="$P_SPRAY_ROOT/build/lib"
export FOAM_USER_APPBIN="$P_SPRAY_ROOT/build/bin"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:${LD_LIBRARY_PATH:-}"
export PATH="$FOAM_USER_APPBIN:$PATH"
unset p_spray_bashrc p_spray_candidate
