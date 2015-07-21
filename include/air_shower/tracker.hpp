/* 

   calin/air_shower/tracker.hpp -- Stephen Fegan -- 2015-06-23

   Base class for all air shower track visitors

*/

#pragma once

// Units:
// Time, t:      ns
// Height, z:    cm
// Energy, e:    MeV

#include"Eigen/Core"

namespace calin { namespace air_shower { namespace tracker {

enum class ParticleType { GAMMA, ELECTRON, POSITRON, MUON, ANTI_MUON,
    PROTON, ANTI_PROTON, OTHER };

struct Track
{
  ParticleType type;     // Simplified particle type
  int pdg_type;          // PDG particle type code
  double q;              // PDG particle charge          [e]
  double mass;           // PDG particle rest mass       [MeV]

  Eigen::Vector3d x0;    // Position of start of track   [cm]
  Eigen::Vector3d u0;    // Direction of start of track  [1]
  double e0;             // Total energy at start of trk [MeV]
  double t0;             // Time at start of track       [ns]
  
  Eigen::Vector3d x1;    // Position of end of track     [cm]
  Eigen::Vector3d u1;    // Direction of end of track    [1]
  double e1;             // Total energy at end of track [MeV]
  double t1;             // Time at end of track         [ns]

  double dx;             // Step length                  [cm]

  double weight;         // Track weighting for thinning
};

class TrackVisitor
{
 public:
  virtual ~TrackVisitor();
  virtual void visitTrack(const Track& track, bool& kill_track);
};

} } } // namespace calin::air_shower::tracker
