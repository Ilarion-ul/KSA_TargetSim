# ROOT analysis helpers are distributed as macros; no compilation required.
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/plot_run_summary.C
  ${CMAKE_CURRENT_SOURCE_DIR}/compare_configs.C
  DESTINATION share/KSA_TargetSim/analysis/root)
