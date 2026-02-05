cmake_minimum_required(VERSION 3.20)

project(KSA_TargetSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(KSA_USE_ROOT "Enable ROOT output support when ROOT is available" ON)
option(KSA_ENABLE_TESTS "Enable CTest-based smoke tests" ON)
option(KSA_BUILD_ANALYSIS "Build analysis helpers/macros" ON)
option(KSA_ENABLE_MT "Enable Geant4 multithreading support when available" ON)

# Keep binaries/libraries in predictable build-local folders.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(Geant4 REQUIRED ui_all vis_all)

# Geant4's CMake module exports include directories/definitions via Geant4_USE_FILE.
# This handles platform-specific compile definitions and environment glue.
include(${Geant4_USE_FILE})

set(KSA_ROOT_ENABLED OFF)
if(KSA_USE_ROOT)
  find_package(ROOT QUIET COMPONENTS Core RIO Hist Tree)
  if(ROOT_FOUND)
    set(KSA_ROOT_ENABLED ON)
    message(STATUS "ROOT found: enabling ROOT-backed outputs")
  else()
    set(KSA_USE_ROOT OFF)
    message(WARNING "KSA_USE_ROOT=ON but ROOT was not found. Falling back to text logs in results/logs.")
  endif()
endif()

add_subdirectory(app)

if(KSA_BUILD_ANALYSIS)
  add_subdirectory(analysis/root)
endif()

if(KSA_ENABLE_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

# Install executable and example runtime inputs.
install(PROGRAMS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ksasim"
        DESTINATION bin
        OPTIONAL)
install(DIRECTORY app/macros/ app/config/
        DESTINATION share/KSA_TargetSim/app)

# Results are written at runtime into results/{logs,root,vis}.
# Geant4 runtime datasets should be prepared by sourcing geant4.sh before running.
