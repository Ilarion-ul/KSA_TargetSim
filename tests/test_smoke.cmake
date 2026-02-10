if(NOT DEFINED KSASIM_BIN)
  message(FATAL_ERROR "KSASIM_BIN is not defined")
endif()
if(NOT DEFINED KSASIM_SOURCE_DIR)
  message(FATAL_ERROR "KSASIM_SOURCE_DIR is not defined")
endif()

set(CONFIG "${KSASIM_SOURCE_DIR}/app/config/default_WTa.json")
set(MACRO "${KSASIM_SOURCE_DIR}/app/macros/run_batch.mac")

execute_process(
  COMMAND ${KSASIM_BIN} -c ${CONFIG} -m ${MACRO} -n 50
  WORKING_DIRECTORY ${KSASIM_SOURCE_DIR}
  RESULT_VARIABLE run_code
)
if(NOT run_code EQUAL 0)
  message(FATAL_ERROR "Smoke run failed with exit code ${run_code}")
endif()

set(LOG_FILE "${KSASIM_SOURCE_DIR}/results/logs/run_summary.json")
file(GLOB ROOT_FILES "${KSASIM_SOURCE_DIR}/results/root/*.root")

if(EXISTS "${LOG_FILE}")
  file(SIZE "${LOG_FILE}" LOG_SIZE)
  if(LOG_SIZE GREATER 0)
    message(STATUS "Smoke output log exists and is non-empty: ${LOG_FILE}")
    return()
  endif()
endif()

set(ROOT_OK FALSE)
foreach(rf IN LISTS ROOT_FILES)
  file(SIZE "${rf}" RF_SIZE)
  if(RF_SIZE GREATER 0)
    set(ROOT_OK TRUE)
  endif()
endforeach()

if(NOT ROOT_OK)
  message(FATAL_ERROR "Smoke output missing: expected non-empty results/logs/run_summary.json or results/root/*.root")
endif()
