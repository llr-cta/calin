//-*-mode:swig;-*-

/*

   calin/simulation/tracker.i -- Stephen Fegan -- 2016-10-06

   SWIG interface file for calin.simulation.tracker

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

%module (package="calin.simulation") tracker
%feature(autodoc,2);

%{
#include "simulation/tracker.hpp"
#include "simulation/misc_trackers.hpp"
#include "simulation/straight_track_generator.hpp"
#include "simulation/bfield_track_generator.hpp"
#include "simulation/air_cherenkov_tracker.hpp"
#include "simulation/iact_array_tracker.hpp"
#include "simulation/quadrature_iact_array_integration.hpp"
#include "simulation/vso_quadrature_iact_array_integration.hpp"
#include "simulation/vcl_ray_processor.hpp"
#define SWIG_FILE_WITH_INIT
  %}

%init %{
  import_array();
%}

%include "calin_typemaps.i"
%import "calin_global_definitions.i"

%import "simulation/atmosphere.i"
%import "simulation/world_magnetic_model.i"
%import "simulation/detector_efficiency.i"
%import "simulation/vs_optics.i"
%import "simulation/tracker.pb.i"
%import "math/histogram.i"
%import "simulation/ray_processor.i"
//%newobject *::dump_as_proto() const;
//%newobject *::create_from_proto;

%include "simulation/tracker.hpp"

%template (VectorSimulationTrack) std::vector<calin::simulation::tracker::Event>;
%template (VectorSimulationEvent) std::vector<calin::simulation::tracker::Track>;

%apply Eigen::Vector3d &OUTPUT { Eigen::Vector3d& x0, Eigen::Vector3d& x1 };
%apply Eigen::MatrixXd &OUTPUT { Eigen::MatrixXd& x0, Eigen::MatrixXd& x1 };

%include "simulation/misc_trackers.hpp"
%include "simulation/straight_track_generator.hpp"
%include "simulation/bfield_track_generator.hpp"
%include "simulation/air_cherenkov_tracker.hpp"
%include "simulation/iact_array_tracker.hpp"
%include "simulation/quadrature_iact_array_integration.hpp"
%include "simulation/vso_quadrature_iact_array_integration.hpp"
