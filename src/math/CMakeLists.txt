add_library(calin_math SHARED accumulator.cpp histogram.cpp function.cpp
	pdf_1d.cpp optimizer.cpp nlopt_optimizer.cpp
	minuit75_optimizer.cpp hessian.cpp brent.cpp)
target_link_libraries(calin_math calin_f2c calin_nlopt calin_proto)

install(TARGETS calin_math DESTINATION lib)
