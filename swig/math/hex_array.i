//-*-mode:swig;-*-

/* 

   calin/math/hex_array.i -- Stephen Fegan -- 2015-10-29

*/

%module (package="calin.math") hex_array

%{
#include "math/hex_array.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "package_wide_definitions.i"

%apply unsigned &OUTPUT { unsigned& ringid, unsigned& segid, unsigned& runid };
%apply int &OUTPUT { int& u, int& v };
%apply std::vector<int> &OUTPUT { std::vector<int>& u_neighbors,
       std::vector<int>& v_neighbors,
       std::vector<int>& u, std::vector<int>& v
       };
%apply double &OUTPUT { double& x, double& y };
%apply double &INOUT { double& x_in_dx_out, double& y_in_dy_out };
%apply std::vector<double> &OUTPUT { std::vector<double>& x,
       std::vector<double>& y };

%include "math/hex_array.hpp"
