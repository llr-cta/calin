add_library(calin_calib SHARED spe_fit.cpp pmt_model_pg.cpp)
target_link_libraries(calin_calib calin_math) 
target_link_libraries(calin_calib ${EXTRA_LIBS}) 
