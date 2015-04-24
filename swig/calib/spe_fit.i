//-*-mode:swig;-*-

/* 

   calin/math/spe_fit.i -- Stephen Fegan -- 2015-04-24

*/

%module (package="calin.calib") spe_fit

%{
#include "calib/spe_fit.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}


%include "package_wide_definitions.i"

%import "math/function.i"
%import "math/histogram.i"

%include "calib/spe_fit.hpp"
