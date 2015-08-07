//-*-mode:swig;-*-

/* 

   calin/math/special.i -- Stephen Fegan -- 2015-08-06

*/

%module (package="calin.math") special

%{
#include "math/special.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%include "math/special.hpp"
