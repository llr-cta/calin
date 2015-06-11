//-*-mode:swig;-*-

/* 

   calin/math/atmosphere.i -- Stephen Fegan -- 2015-06-11

*/

%module (package="calin.air_shower") atmosphere

%{
#include "air_shower/atmosphere.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%template (VectorAtmSlice) std::vector<calin::air_shower::atmosphere::AtmSlice>;
%template (VectorLayeredAtmosphereLevel) std::vector<calin::air_shower::atmosphere::LayeredAtmosphereLevel>;

%include "air_shower/atmosphere.hpp"
