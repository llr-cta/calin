//-*-mode:swig;-*-

/*

   calin/calin_global_definitions.i -- Stephen Fegan -- 2015-04-21

   SWIG interface file for calin.package_wide_definitions

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

%module (package="calin") calin_global_definitions
%feature(autodoc,2);

%{
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "calin_global_definitions.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%include "calin_typemaps.i"

%include <std_pair.i>
%include <std_vector.i>
%include <std_string.i>
%include <std_map.i>
%include "calin_stdint.i"

%init %{
  import_array();
%}

// alignas not supported by SWIG (yet)
#define alignas(T)

%template (VectorDouble)      std::vector<double>;
%template (VectorFloat)       std::vector<float>;
%template (VectorChar)        std::vector<char>;
%template (VectorUChar)       std::vector<unsigned char>;
%template (VectorShort)       std::vector<short>;
%template (VectorUShort)      std::vector<unsigned short>;
%template (VectorInt)         std::vector<int>;
%template (VectorUInt)        std::vector<unsigned int>;
%template (VectorLong)        std::vector<long>;
%template (VectorULong)       std::vector<unsigned long>;
%template (VectorLongLong)    std::vector<long long>;
%template (VectorULongLong)   std::vector<unsigned long long>;
%template (VectorBool)        std::vector<bool>;
%template (VectorString)      std::vector<std::string>;

%template (PairDoubleDouble) std::pair<double,double>;
%template (VectorPairDoubleDouble) std::vector<std::pair<double,double> >;

%template (PairUnsignedUnsigned) std::pair<unsigned,unsigned>;
%template (VectorPairUnsignedUnsigned) std::vector<std::pair<unsigned,unsigned> >;

%template (MapStringString)      std::map<std::string,std::string>;

%import "calin_global_definitions.hpp"
