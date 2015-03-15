add_executable(do_spe_fit do_spe_fit)  
target_link_libraries(do_spe_fit calin_calib calin_math calin_proto calin_nlopt ${EXTRA_LIBS})
