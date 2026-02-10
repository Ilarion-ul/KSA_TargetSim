#!/usr/bin/env bash
# shellcheck disable=SC1090,SC2155
set -euo pipefail

# Helper: source first existing file from arguments.
source_first_found() {
  for p in "$@"; do
    if [[ -f "$p" ]]; then
      # shellcheck source=/dev/null
      source "$p"
      echo "[env] sourced: $p"
      return 0
    fi
  done
  return 1
}

echo "[env] Preparing Geant4/ROOT runtime environment..."

# --- Geant4 block -----------------------------------------------------------
# If user already has Geant4 hints, keep their environment and report it.
if [[ -n "${Geant4_DIR:-}" || -n "${G4INSTALL:-}" ]]; then
  echo "[env] Geant4 appears configured already (Geant4_DIR/G4INSTALL present)."
else
  echo "[env] Geant4 not preconfigured, probing typical setup scripts..."
  source_first_found \
    /opt/geant4/*/bin/geant4.sh \
    "$HOME"/opt/geant4/*/bin/geant4.sh \
    /usr/local/share/Geant4-*/geant4.sh \
    /usr/local/bin/geant4.sh \
    || echo "[env] Geant4 setup script not found automatically."
fi

# --- ROOT block -------------------------------------------------------------
# Prefer existing ROOT environment (root-config in PATH). If absent, probe.
if command -v root-config >/dev/null 2>&1; then
  echo "[env] ROOT already available via root-config."
else
  echo "[env] ROOT not preconfigured, probing thisroot.sh..."
  source_first_found \
    /opt/root/*/bin/thisroot.sh \
    "$HOME"/opt/root/*/bin/thisroot.sh \
    /usr/local/root/bin/thisroot.sh \
    /usr/local/bin/thisroot.sh \
    || echo "[env] ROOT setup script not found automatically."
fi

# --- Report versions --------------------------------------------------------
if command -v geant4-config >/dev/null 2>&1; then
  echo "[env] Geant4 version: $(geant4-config --version)"
else
  # Optional fallback: try to infer from headers if available.
  g4_header=$(find /usr/include /usr/local/include -maxdepth 4 -name G4Version.hh 2>/dev/null | head -n 1 || true)
  if [[ -n "$g4_header" ]]; then
    g4_ver=$(grep -E 'G4VERSION_NUMBER|G4VERSION_TAG' "$g4_header" | head -n 1 || true)
    echo "[env] Geant4 version (header hint): ${g4_ver:-unknown}"
  else
    echo "[env] Geant4 version: unavailable (geant4-config not found)."
  fi
fi

if command -v root-config >/dev/null 2>&1; then
  echo "[env] ROOT version: $(root-config --version)"
else
  echo "[env] ROOT version: unavailable (root-config not found)."
fi

# --- Print Geant4 data variables -------------------------------------------
echo "[env] Geant4 DATA-related environment variables:"
if env | grep -E '^G4.*DATA=' >/dev/null 2>&1; then
  env | grep -E '^G4.*DATA=' | sort
else
  echo "[env] (none found)"
fi

# --- Final hints ------------------------------------------------------------
if ! command -v geant4-config >/dev/null 2>&1; then
  echo "[hint] Source Geant4 manually, e.g. /path/to/geant4.sh"
fi
if ! command -v root-config >/dev/null 2>&1; then
  echo "[hint] ROOT optional: source /path/to/thisroot.sh for ROOT output support"
fi
