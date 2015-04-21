add_executable(do_spe_fit do_spe_fit)  
target_link_libraries(do_spe_fit calin_calib calin_math calin_proto calin_nlopt ${EXTRA_LIBS})

add_executable(do_spe_expgauss_fit do_spe_expgauss_fit)  
target_link_libraries(do_spe_expgauss_fit calin_calib calin_math calin_proto calin_nlopt ${EXTRA_LIBS})

install(TARGETS do_spe_fit do_spe_expgauss_fit DESTINATION bin)
