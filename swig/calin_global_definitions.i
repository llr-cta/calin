//-*-mode:swig;-*-

/*

   calin/calin_global_definitions.i -- Stephen Fegan -- 2015-04-21

   SWIG interface file for calin.package_wide_definitions

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

%module (package="calin") calin_global_definitions

%{
#include <iostream>
#include "calin_global_definitions.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%include "calin_typemaps.i"

%include <std_vector.i>
%include <std_string.i>
%include <stdint.i>

%init %{
  import_array();
%}

%template (VectorDouble)     std::vector<double>;
%template (VectorFloat)      std::vector<float>;
%template (VectorUInt32)     std::vector<uint32_t>;
%template (VectorInt32)      std::vector<int32_t>;
%template (VectorUInt64)     std::vector<uint64_t>;
%template (VectorInt64)      std::vector<int64_t>;
%template (VectorBool)       std::vector<bool>;
%template (VectorString)     std::vector<std::string>;

%import "calin_global_definitions.hpp"
