//-*-mode:swig;-*-

/* 

   calin/math/optimizer.i -- Stephen Fegan -- 2015-04-23

*/

%module (package="calin.math") optimizer

%{
#include "math/optimizer.hpp"
#include "math/nlopt_optimizer.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%import "math/function.i"

%include "nlopt/nlopt.hpp"
%include "math/optimizer.hpp"
%include "math/nlopt_optimizer.hpp"
