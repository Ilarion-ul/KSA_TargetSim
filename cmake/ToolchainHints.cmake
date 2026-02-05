# ToolchainHints.cmake
# -----------------------------------------------------------------------------
# This file is intentionally a human-readable hint sheet, not an active
# toolchain file. Copy/paste examples as needed in your local configure command.

# Typical ways to point CMake to Geant4/ROOT installations:
#   -DGeant4_DIR=/path/to/lib/cmake/Geant4
#   -DROOT_DIR=/path/to/lib/cmake/ROOT
#   -DCMAKE_PREFIX_PATH="/opt/geant4;/opt/root"

# Example configure command (Ninja):
#   cmake -S . -B build -G Ninja \
#     -DGeant4_DIR=/opt/geant4/lib/cmake/Geant4 \
#     -DROOT_DIR=/opt/root/lib/cmake/ROOT \
#     -DKSA_USE_ROOT=ON

# Linux / WSL / cluster note:
# Do not hardcode LD_LIBRARY_PATH edits in CMakeLists. Prefer sourcing runtime
# setup scripts (for example geant4.sh / thisroot.sh) in your shell session.

# Minimal variable examples (disabled by default):
# set(Geant4_DIR "/opt/geant4/lib/cmake/Geant4" CACHE PATH "Path to Geant4 CMake package")
# set(ROOT_DIR "/opt/root/lib/cmake/ROOT" CACHE PATH "Path to ROOT CMake package")
# set(CMAKE_PREFIX_PATH "/opt/geant4;/opt/root" CACHE STRING "Extra package roots")

# Multithreading note:
# Geant4 must be built with G4MULTITHREADED support to use worker threads.
# KSA_TargetSim thread count is configured at runtime via JSON/CLI (-t), where
# nThreads=0 means "auto" according to Geant4 defaults.
