file(GLOB KSA_APP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(ksasim ${KSA_APP_SOURCES})
target_compile_features(ksasim PRIVATE cxx_std_17)

target_include_directories(ksasim
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Geant4 compile definitions/include setup via legacy use file + libs.
target_link_libraries(ksasim PRIVATE ${Geant4_LIBRARIES})

# Optional nlohmann/json (header-only): prefer package, otherwise external include path.
find_package(nlohmann_json QUIET)
if(nlohmann_json_FOUND)
  target_link_libraries(ksasim PRIVATE nlohmann_json::nlohmann_json)
endif()

if(KSA_ROOT_ENABLED)
  target_compile_definitions(ksasim PRIVATE KSA_USE_ROOT)
  if(KSA_ROOT_VIA_CMAKE)
    target_link_libraries(ksasim PRIVATE ROOT::Core ROOT::RIO ROOT::Hist ROOT::Tree)
  else()
    # root-config fallback path: use flags directly for distro ROOT packages.
    if(KSA_ROOT_CFLAGS)
      separate_arguments(KSA_ROOT_CFLAGS_LIST NATIVE_COMMAND "${KSA_ROOT_CFLAGS}")
      target_compile_options(ksasim PRIVATE ${KSA_ROOT_CFLAGS_LIST})
    endif()
    if(KSA_ROOT_LIBFLAGS)
      separate_arguments(KSA_ROOT_LIBFLAGS_LIST NATIVE_COMMAND "${KSA_ROOT_LIBFLAGS}")
      target_link_options(ksasim PRIVATE ${KSA_ROOT_LIBFLAGS_LIST})
    endif()
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(ksasim PRIVATE -Wall -Wextra -Wpedantic)
endif()

set_target_properties(ksasim PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Convenience runtime copy for macros/config files.
set(KSA_RUNTIME_SHARE "${CMAKE_BINARY_DIR}/share/KSA_TargetSim/app")
add_custom_command(TARGET ksasim POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${KSA_RUNTIME_SHARE}/macros"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${KSA_RUNTIME_SHARE}/config"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/macros" "${KSA_RUNTIME_SHARE}/macros"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/config" "${KSA_RUNTIME_SHARE}/config"
)

install(TARGETS ksasim RUNTIME DESTINATION bin)
