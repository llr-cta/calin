/*

   calin/math/data_modeling.i -- Stephen Fegan -- 2017-04-07

   Data modeling functions for various data types

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin.math") data_modeling
%feature(autodoc,2);

%{
#include "math/data_modeling.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

//%include "math/accumulator.i"
%import "math/function.i"
%import "math/histogram.i"

%include "math/m_estimate.hpp"
%include "math/data_modeling.hpp"

%extend calin::math::data_modeling::IID1DDataLikelihoodFunction {
  %template(IID1DDataLikelihoodFunction) IID1DDataLikelihoodFunction<
    calin::math::histogram::BasicHistogram1D<
      calin::math::histogram::DefaultAccumulator>>;
}

%extend calin::math::data_modeling::IID1DDataMEstimateLikelihoodFunction {
  %template(IID1DDataMEstimateLikelihoodFunction) IID1DDataMEstimateLikelihoodFunction<
    calin::math::histogram::BasicHistogram1D<
      calin::math::histogram::DefaultAccumulator>>;
}

%extend calin::math::data_modeling::IID1DDataChi2Function {
  %template(IID1DDataChi2Function) IID1DDataChi2Function<
    calin::math::histogram::BasicHistogram1D<
      calin::math::histogram::DefaultAccumulator>>;
}
