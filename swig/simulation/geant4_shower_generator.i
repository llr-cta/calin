/*

   calin/simulation/geant4_shower_generator.i -- Stephen Fegan -- 2016-11-10

   SWIG interface file for Geant 4 shower generator

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

%module (package="calin.simulation") geant4_shower_generator
%feature(autodoc,2);

%{
//#include <G4Exception.hh>
#include "simulation/tracker.hpp"
#include "simulation/misc_trackers.hpp"
#include "simulation/straight_track_generator.hpp"
#include "simulation/air_cherenkov_tracker.hpp"
#include "simulation/iact_array_tracker.hpp"
#include "simulation/quadrature_iact_array_integration.hpp"
#include "simulation/vso_quadrature_iact_array_integration.hpp"
#include "simulation/geant4_shower_generator.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

/*
%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const G4Exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}
*/

%import "simulation/atmosphere.i"
%import "simulation/tracker.i"

%include "simulation/geant4_shower_generator.hpp"
