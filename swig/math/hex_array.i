//-*-mode:swig;-*-

/*

   calin/math/hex_array.i -- Stephen Fegan -- 2015-10-29

   SWIG interface file for calin.math.hex_array

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") hex_array
%feature(autodoc,2);

%{
#include "math/hex_array.hpp"
/* #include "math/hex_array_simd.hpp" */
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply unsigned &OUTPUT { unsigned& ringid, unsigned& segid, unsigned& runid };
%apply int &OUTPUT { int& uout, int& vout };
%apply std::vector<int> &OUTPUT { std::vector<int>& u_neighbors,
       std::vector<int>& v_neighbors,
       std::vector<int>& u, std::vector<int>& v
       };
%apply double &OUTPUT { double& x, double& y };
%apply double &INOUT { double& x_in_dx_out, double& y_in_dy_out };
%apply std::vector<double> &OUTPUT { std::vector<double>& x,
       std::vector<double>& y };

%apply float &OUTPUT { float& x, float& y };
%apply float &INOUT { float& x_in_dx_out, float& y_in_dy_out };

%apply int &INOUT { int& uinout, int& vinout };

%include "math/hex_array.hpp"
/* %include "math/hex_array_simd.hpp" */
