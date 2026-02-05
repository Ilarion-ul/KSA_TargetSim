add_test(NAME ksasim_help
  COMMAND ${CMAKE_BINARY_DIR}/bin/ksasim --help)
set_tests_properties(ksasim_help PROPERTIES TIMEOUT 60)

add_test(NAME ksasim_smoke_batch
  COMMAND ${CMAKE_COMMAND}
    -DKSASIM_BIN=${CMAKE_BINARY_DIR}/bin/ksasim
    -DKSASIM_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_smoke.cmake)
set_tests_properties(ksasim_smoke_batch PROPERTIES TIMEOUT 60)
