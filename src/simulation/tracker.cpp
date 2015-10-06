/* 

   calin/simulation/tracker.cpp -- Stephen Fegan -- 2015-07-17

   Base class for all air shower track visitors

*/

#include "simulation/tracker.hpp"

using namespace calin::simulation::tracker;

TrackVisitor::~TrackVisitor()
{
  // nothing to see here
}

void TrackVisitor::visitTrack(const Track& track, bool& kill_track)
{
  // default is to do nothing
}
