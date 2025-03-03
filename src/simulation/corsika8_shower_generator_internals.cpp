/*

   calin/simulation/corsika8_shower_generator_impl.cpp -- Stephen Fegan -- 2024-08-23

   Class to generate extensive air showers using Corsika 8

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

// #include <Eigen/Dense>

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/process/DynamicInteractionProcess.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/IMagneticFieldModel.hpp>
#include <corsika/media/MediumPropertyModel.hpp>
#include <corsika/media/UniformMagneticField.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <corsika/modules/writers/WriterOff.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/Epos.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <util/log.hpp>
#include <simulation/corsika8_shower_generator.hpp>
#include <provenance/chronicle.hpp>

using namespace calin::simulation::corsika8_shower_generator;
using namespace calin::util::log;
using namespace corsika;

namespace {

  class TrackHandoff: public ContinuousProcess<TrackHandoff>, public WriterOff {
  public:
    TrackHandoff(double r_earth): WriterOff("track.parquet"), r_earth_(r_earth) { };
    virtual ~TrackHandoff();

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle> const&, bool const limitFlag);

    template <typename TParticle, typename TTrack>
    LengthType getMaxStepLength(TParticle const&, TTrack const&);

    YAML::Node getConfig() const override;

    void startOfShower(unsigned int const /*showerId*/) override;
    void endOfShower(unsigned int const showerId) override;

    void set_visitor(calin::simulation::tracker::TrackVisitor* visitor) { visitor_ = visitor; }
    void clear_visitor() { visitor_ = nullptr; }
  private:
    double r_earth_ = 0;
    unsigned int shower_id_ = 0;
    calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  }; // class TrackWriterHandoff

  class CORSIKA8ShowerGeneratorImpl: public CORSIKA8ShowerGenerator {
  public:
    using EnvironmentInterface = IMagneticFieldModel<IMediumModel>;
    template <typename T> using MyExtraEnv = UniformMagneticField<T>;
    // using EnvironmentInterface =
    //   IRefractiveIndexModel<IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>>;
    using EnvType = Environment<EnvironmentInterface>;
    using StackType = setup::Stack<EnvType>;

    CALIN_TYPEALIAS(config_type,
                    calin::ix::simulation::corsika8_shower_generator::CORSIKA8ShowerGeneratorConfiguration);

    CORSIKA8ShowerGeneratorImpl(const config_type& config = default_config());
    virtual ~CORSIKA8ShowerGeneratorImpl();

    void generate_showers(calin::simulation::tracker::TrackVisitor* visitor,
                          unsigned num_events,
                          calin::simulation::tracker::ParticleType type,
                          double total_energy,
                          const Eigen::Vector3d& x0 = Eigen::Vector3d(0,0,0),
                          const Eigen::Vector3d& u0 = Eigen::Vector3d(0,0,-1),
                          double weight=1.0) override;

  private:
    config_type config_; 

    EnvType env_;
    CoordinateSystemPtr root_cs_ = env_.getCoordinateSystem();
    Point center_ { root_cs_, 0_m, 0_m, 0_m };
    Point ground_; // set in constructor from value passed in config
    DynamicInteractionProcess<StackType> he_model_;
    std::shared_ptr<corsika::sibyll::Interaction> sibyll_;

    std::unique_ptr<TrackHandoff> track_handoff_;
  };

} // anonymous namespace

TrackHandoff::~TrackHandoff() { 
   // nothing to see here
} 

template <typename TParticle>
ProcessReturn TrackHandoff::doContinuous(Step<TParticle> const& step, bool const limitFlag)
{
  if(!visitor_) return ProcessReturn::Ok;

  const auto& particle_pre { step.getParticlePre() };
  const auto& particle_post { step.getParticlePost() };

  calin::simulation::tracker::Track track;
  track.track_id        = 0; // what to do here
  track.parent_track_id = 0; // what to do here

  track.pdg_type        = particle_pre.getPID();
  track.q               = particle_pre.getChargeNumber();
  track.mass            = particle_pre.getMass()/1_MeV;
  track.type            = calin::simulation::tracker::pdg_type_to_particle_type(track.pdg_type);

  const auto& x0 { particle_pre.getPosition() };
  const auto& u0 { particle_pre.getDirection() };

  track.e0              = particle_pre.getEnergy()/1_MeV;
  track.x0              << x0.getX()/1_cm, x0.getY()/1_cm, (x0.getZ()-r_earth_)/1_cm;
  track.u0              << u0.getX(), u0.getY(), u0.getZ();
  track.t0              = particle_pre->getTime()/1_ns;

  const auto& x1 { particle_post.getPosition() };
  const auto& u1 { particle_post.getDirection() };

  track.e1              = particle_post.getEnergy()/1_MeV;
  track.x1              << x1.getX()/1_cm, x1.getY()/1_cm, (x1.getZ()-r_earth_)/1_cm;
  track.u1              << u1.getX(), u1.getY(), u1.getZ();
  track.t1              = particle_post->getTime()/1_ns;

  track.dx_hat          = track.x1 - track.x0;
  track.dx              = track.dx_hat.norm();
  track.dx_hat          /= track.dx;
  track.de              = track.e1 - track.e0;
  track.dt              = track.t1 - track.t0;
  track.weight          = particle_pre.getWeight();
  
#if 0
  static int nprint=0;
  if(nprint<100 and pdg_info->GetPDGEncoding()==2112) {
    LOG(INFO) << particle_pre.getKineticEnergy()/1_MeV; << ' ' << track.dx;
    ++nprint;
  }
#endif
  
  bool kill_track = false;
  visitor_->visit_track(track, kill_track);
  return kill_track ? ProcessReturn::ParticleAbsorbed : ProcessReturn::Ok;
}

template <typename TParticle, typename TTrack>
LengthType TrackHandoff::getMaxStepLength(TParticle const&, TTrack const&)
{
  return 1_m * std::numeric_limits<double>::infinity();
}

YAML::Node TrackHandoff::getConfig() const
{
  return YAML::Node{};
}

void TrackHandoff::startOfShower(unsigned int const showerId)
{
  shower_id_ = showerId;
}

void TrackHandoff::endOfShower(unsigned int const showerId)
{
  shower_id_ = 0;
}

CORSIKA8ShowerGeneratorImpl::
CORSIKA8ShowerGeneratorImpl(const CORSIKA8ShowerGeneratorImpl::config_type& config):
  CORSIKA8ShowerGenerator(), config_(config),  
  ground_ {root_cs_, 0_m, 0_m, (config_.earth_radius() + config_.zground())*1_cm}
{

  // ==========================================================================  
  // Initialize random number sequence(s)
  // ==========================================================================  

  uint32_t seed = config_.seed();
  while(seed == 0)seed = calin::math::rng::RNG::uint32_from_random_device();
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  // RNGManager<>::getInstance().registerRandomStream("fluka");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("thinning");
  // RNGManager<>::getInstance().registerRandomStream("primary_particle");
  RNGManager<>::getInstance().setSeed(seed);
  calin::provenance::chronicle::register_external_rng_open(seed, "CORSIKA8::RNGManager",
    __PRETTY_FUNCTION__);

  // ==========================================================================
  // SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM
  // ==========================================================================
  
  if(config_.atmospheric_model() < calin::ix::simulation::corsika8_shower_generator::ATM_CORSIKA_MAX) {
    // Build a standard CORSIKA atmosphere into env_
    create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env_, static_cast<AtmosphereId>(config_.atmospheric_model()), center_, 
      MagneticFieldVector{root_cs_, 
        config_.uniform_magnetic_field().x()*1_nT, 
        config_.uniform_magnetic_field().y()*1_nT, 
        config_.uniform_magnetic_field().z()*1_nT });
  } else {
    throw std::runtime_error("Custom atmospheric model not supported yet");
  }

  // ==========================================================================
  // SETUP HADRONIC INTERACTIONS
  // ==========================================================================

  const auto all_elements = corsika::get_all_elements_in_universe(env_);
  // have SIBYLL always for PROPOSAL photo-hadronic interactions
  sibyll_ = std::make_shared<corsika::sibyll::Interaction>(all_elements, corsika::setup::C7trackedParticles);

  switch(config_.he_hadronic_model()) {
  case calin::ix::simulation::corsika8_shower_generator::SIBYLL:
  default:
    he_model_ = DynamicInteractionProcess<StackType>{sibyll_};
    break;
  case calin::ix::simulation::corsika8_shower_generator::QGSJet:
    he_model_ = DynamicInteractionProcess<StackType>{
      std::make_shared<corsika::qgsjetII::Interaction>()};
    break;
  case calin::ix::simulation::corsika8_shower_generator::EPOS_LHC:
    he_model_ = DynamicInteractionProcess<StackType>{
      std::make_shared<corsika::epos::Interaction>(corsika::setup::C7trackedParticles)};
    break;    
  };

  track_handoff_ = std::make_unique<TrackHandoff>(config_.earth_radius());
}

CORSIKA8ShowerGeneratorImpl::~CORSIKA8ShowerGeneratorImpl()
{
  // nothing to see here
}

void CORSIKA8ShowerGeneratorImpl::
generate_showers(calin::simulation::tracker::TrackVisitor* visitor,
                 unsigned num_events,
                 calin::simulation::tracker::ParticleType type,
                 double total_energy,
                 const Eigen::Vector3d& x0,
                 const Eigen::Vector3d& u0,
                 double weight)
{
  track_handoff_->set_visitor(visitor);
  
  track_handoff_->clear_visitor();
}
