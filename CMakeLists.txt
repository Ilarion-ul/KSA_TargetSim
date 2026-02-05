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
set(KSA_ROOT_VIA_CMAKE OFF)
set(KSA_ROOT_CFLAGS "")
set(KSA_ROOT_LIBFLAGS "")
if(KSA_USE_ROOT)
  # Prefer root-config for distro ROOT packages where ROOTConfig.cmake may be
  # broken (for example cmake_policy PUSH/POP mismatch in packaged files).
  find_program(ROOT_CONFIG_EXECUTABLE root-config)
  if(ROOT_CONFIG_EXECUTABLE)
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
      OUTPUT_VARIABLE KSA_ROOT_CFLAGS
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _ksa_root_cflags_rc)
    execute_process(
      COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
      OUTPUT_VARIABLE KSA_ROOT_LIBFLAGS
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _ksa_root_libs_rc)

    if(_ksa_root_cflags_rc EQUAL 0 AND _ksa_root_libs_rc EQUAL 0)
      set(KSA_ROOT_ENABLED ON)
      message(STATUS "ROOT enabled via root-config fallback")
    else()
      set(KSA_USE_ROOT OFF)
      message(WARNING "root-config found, but cflags/libs could not be queried. Falling back to text logs.")
    endif()
  else()
    # Try modern CMake package only if root-config is unavailable.
    # This avoids hard configure failures from broken distro ROOTConfig.cmake.
    find_package(ROOT QUIET COMPONENTS Core RIO Hist Tree)
    if(ROOT_FOUND)
      set(KSA_ROOT_ENABLED ON)
      set(KSA_ROOT_VIA_CMAKE ON)
      message(STATUS "ROOT found via CMake package: enabling ROOT-backed outputs")
    else()
      set(KSA_USE_ROOT OFF)
      message(WARNING "KSA_USE_ROOT=ON but ROOT was not found. Falling back to text logs in results/logs.")
    endif()
  endif()

  unset(_ksa_root_cflags_rc)
  unset(_ksa_root_libs_rc)
  mark_as_advanced(ROOT_CONFIG_EXECUTABLE)
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
