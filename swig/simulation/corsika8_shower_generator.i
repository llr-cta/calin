/*

   calin/simulation/geant4_shower_generator.i -- Stephen Fegan -- 201624-08-25

   SWIG interface file for CORSIKA 8 shower generator

   Copyright 2024, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation", threads=1) corsika8_shower_generator
%feature(autodoc,2);

%nothread;

%{
//#include <G4Exception.hh>
#include "simulation/tracker.hpp"
#include "simulation/misc_trackers.hpp"
#include "simulation/straight_track_generator.hpp"
#include "simulation/air_cherenkov_tracker.hpp"
#include "simulation/iact_array_tracker.hpp"
#include "simulation/quadrature_iact_array_integration.hpp"
#include "simulation/vso_quadrature_iact_array_integration.hpp"
#include "simulation/corsika8_shower_generator.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "simulation/atmosphere.i"
%import "simulation/tracker.i"

%import "simulation/corsika8_shower_generator.pb.i"

%thread;
%include "simulation/corsika8_shower_generator.hpp"
%nothread;
