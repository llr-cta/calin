//-*-mode:swig;-*-

/*

   calin/math/spline_interpolation.i -- Stephen Fegan -- 2018-12-18

   SWIG interface file for calin.math.spline_interpolation

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.math") spline_interpolation
%feature(autodoc,2);

%{
#include "math/spline_interpolation.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%apply double &OUTPUT { double& value };
%apply double &OUTPUT { double& first_derivative };

%apply double &OUTPUT { double& value0 };
%apply double &OUTPUT { double& value1 };
%apply double &OUTPUT { double& value2 };
%apply double &OUTPUT { double& value3 };
%apply double &OUTPUT { double& value4 };
%apply double &OUTPUT { double& value5 };
//%apply Eigen::VectorXd &OUTPUT { Eigen::VectorXd& xv, Eigen::VectorXd& yv };

%newobject new_regularized_multi_spline(double dx);
%include "math/spline_interpolation.hpp"

%template(cubic_value) calin::math::spline_interpolation::cubic_value<double>;
%template(cubic_1st_derivative) calin::math::spline_interpolation::cubic_1st_derivative<double>;
%template(cubic_1st_derivative_and_value) calin::math::spline_interpolation::cubic_1st_derivative_and_value<double>;
%template(cubic_2nd_derivative) calin::math::spline_interpolation::cubic_2nd_derivative<double>;
%template(cubic_3rd_derivative) calin::math::spline_interpolation::cubic_3rd_derivative<double>;
%template(cubic_integral) calin::math::spline_interpolation::cubic_integral<double>;

%extend calin::math::spline_interpolation::CubicMultiSpline {
  Eigen::VectorXd xknot() const { return $self->xknot_as_eigen(); }
  Eigen::VectorXd yknot(unsigned ispline) const { return $self->yknot_as_eigen(ispline); }
  Eigen::VectorXd dydxknot(unsigned ispline) const { return $self->dydxknot_as_eigen(ispline); }
};
