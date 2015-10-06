//-*-mode:swig;-*-

/* 

   calin/math/atmosphere.i -- Stephen Fegan -- 2015-06-11

*/

%module (package="calin.simulation") atmosphere

%{
#include "simulation/atmosphere.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%template (VectorAtmSlice) std::vector<calin::simulation::atmosphere::AtmSlice>;
%template (VectorLayeredAtmosphereLevel) std::vector<calin::simulation::atmosphere::LayeredAtmosphereLevel>;

%include "simulation/atmosphere.hpp"
