/* 

   calin/air_shower/tracker.hpp -- Stephen Fegan -- 2015-06-23

   Base class for all air shower track visitors

*/

#pragma once

// Units:
// Time, t:      ns
// Height, z:    cm
// Energy, e:    eV

#include"Eigen/Core"

namespace calin { namespace air_shower { namespace tracker {

enum class ParticleType { GAMMA, ELECTRON, POSITRON, MUON_MINUS, MUON_PLUS,
    PROTON_PLUS, PROTON_MINUS, OTHER };

struct Track
{
  ParticleType type;     // Particle type code
  double q;              // Particle charge              [e]
  double e0;             // Particle rest mass           [eV]
  double e;              // Total particle energy        [eV]
  Eigen::Vector3d x;     // Position of start of track   [cm]
  double t;              // Time at start of track       [ns]
  Eigen::Vector3d u;     // Direction of  propogation    [1]
  double ustep;          // Track length                 [cm]
};

class TrackVisitor
{
 public:
  virtual ~TrackVisitor();
  virtual void visitTrack(const Track& track) = 0;
};

} } } // namespace calin::air_shower::tracker
