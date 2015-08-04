//-*-mode:swig;-*-

/* 

   calin/math/pdf_1d.i -- Stephen Fegan -- 2015-04-16

*/

%module (package="calin.math") pdf_1d

%{
#include "math/pdf_1d.hpp"
#include "math/log_quadratic_spline_pdf_1d.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"
%include "typemaps.i"

%import "math/function.i"

%include "math/pdf_1d.hpp"
%include "math/log_quadratic_spline_pdf_1d.hpp"
