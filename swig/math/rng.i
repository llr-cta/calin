//-*-mode:swig;-*-

/*

   calin/math/rng.i -- Stephen Fegan -- 2015-04-15

   SWIG interface file for calin.math.rng

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

%module (package="calin.math") rng
%feature(autodoc,2);

%{
#include "math/rng.hpp"
#include "math/rng_vcl.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "math/rng.pb.i"

%newobject calin::math::rng::RNG::as_proto() const;
%newobject calin::math::rng::RNGCore::as_proto() const;
%include "math/rng.hpp"

%template(NR3_EmulateSIMD_RNGCore_4) calin::math::rng::NR3_EmulateSIMD_RNGCore<4>;

%include "math/rng_vcl.hpp"
%template (VCLRNG128) calin::math::rng::VCLRNG<calin::util::vcl::VCL128Architecture>;
%template (VCLRNG256) calin::math::rng::VCLRNG<calin::util::vcl::VCL256Architecture>;
%template (VCLRNG512) calin::math::rng::VCLRNG<calin::util::vcl::VCL512Architecture>;

%template (VCLRealRNGFloat128) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL128FloatReal>;
%template (VCLRealRNGFloat256) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL256FloatReal>;
%template (VCLRealRNGFloat512) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL512FloatReal>;

%template (VCLRealRNGDouble128) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL128DoubleReal>;
%template (VCLRealRNGDouble256) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL256DoubleReal>;
%template (VCLRealRNGDouble512) calin::math::rng::VCLRealRNG<calin::util::vcl::VCL512DoubleReal>;
