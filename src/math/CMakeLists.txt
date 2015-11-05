set(CALIN_TARGET_LIBRARY calin_math)
set(CALIN_SOURCES accumulator.cpp histogram.cpp function.cpp
	pdf_1d.cpp optimizer.cpp nlopt_optimizer.cpp cminpack_optimizer.cpp
	minuit75_optimizer.cpp hessian.cpp brent.cpp
	log_quadratic_spline_pdf_1d.cpp
	special.cpp hex_array.cpp hex_array_vvv.cpp vec3d.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY}
	calin_nlopt calin_cminpack calin_f2c calin_io calin_proto ${EXTRA_LIBS})

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

