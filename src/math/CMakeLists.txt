set(CALIN_TARGET_LIBRARY calin_math)
set(CALIN_SOURCES accumulator.cpp histogram.cpp function.cpp
	pdf_1d.cpp optimizer.cpp nlopt_optimizer.cpp cminpack_optimizer.cpp
	minuit75_optimizer.cpp hessian.cpp brent.cpp)

add_library(${CALIN_TARGET_LIBRARY} SHARED ${CALIN_SOURCES})
target_link_libraries(${CALIN_TARGET_LIBRARY}
	calin_nlopt calin_cminpack calin_f2c calin_io calin_proto)

install(TARGETS ${CALIN_TARGET_LIBRARY} DESTINATION lib)

