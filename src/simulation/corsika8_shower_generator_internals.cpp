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
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/geometry/CoordinateSystem.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/process/DynamicInteractionProcess.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>

#include <corsika/media/Environment.hpp>
#include <corsika/media/IMagneticFieldModel.hpp>
#include <corsika/media/MediumPropertyModel.hpp>
#include <corsika/media/UniformMagneticField.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/Epos.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
// #include <corsika/modules/TAUOLA.hpp>
#include <corsika/modules/UrQMD.hpp>
#include <corsika/modules/writers/WriterOff.hpp>

#include <corsika/setup/SetupC7trackedParticles.hpp>
#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <Eigen/Core>

#include <util/log.hpp>
#include <math/special.hpp>
#include <math/geometry.hpp>
#include <simulation/corsika8_shower_generator.hpp>
#include <provenance/chronicle.hpp>

using namespace calin::simulation::corsika8_shower_generator;
using namespace calin::util::log;
using namespace corsika;

using calin::math::geometry::box_has_future_intersection;
using calin::math::special::SQR;
using calin::util::log::LOG;

namespace {

  // Class to handle the handoff of a track to a visitor and also a detector boundary box.
  // We don't use the C8 ObservationPlane as it requires inclusion of Apache Parquet which
  // is not part of the calin singularity container.
  class TrackHandoff: public ContinuousProcess<TrackHandoff> {
  public:
    TrackHandoff(double r_earth, double z_bottom, double z_top, double xy_side): 
      r_earth_(r_earth), z_bottom_(z_bottom), z_top_(z_top), xy_side_(xy_side),
      max_step_length_(1_cm * std::sqrt(2.0*SQR(xy_side_) + SQR(z_top_-z_bottom_))),
      min_corner_(-0.5*xy_side_, -0.5*xy_side_, r_earth_+z_bottom_),
      max_corner_( 0.5*xy_side_,  0.5*xy_side_, r_earth_+z_top_) { };
    virtual ~TrackHandoff();

    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle> const&, bool const limitFlag);

    template <typename TParticle, typename TTrack>
    LengthType getMaxStepLength(TParticle const&, TTrack const&);

    void set_visitor(calin::simulation::tracker::TrackVisitor* visitor) { visitor_ = visitor; }
    void clear_visitor() { visitor_ = nullptr; }
  private:
    double r_earth_  = 6371.0 * 1e5;  // Earth radius in cm
    double z_bottom_ = 0.0;
    double z_top_    = 125.0  * 1e5;  // 125km "top of atmosphere" in cm
    double xy_side_  = 1000.0 * 1e5;  // 1000km side of detector in cm
    LengthType max_step_length_ = 1_km;
    Eigen::Vector3d min_corner_;
    Eigen::Vector3d max_corner_;  
    calin::simulation::tracker::TrackVisitor* visitor_ = nullptr;
  }; // class TrackHandoff

} // anonymous namespace

TrackHandoff::~TrackHandoff() { 
  // nothing to see here
} 

template <typename TParticle>
ProcessReturn TrackHandoff::doContinuous(Step<TParticle> const& step, bool const limitFlag)
{
  calin::simulation::tracker::Track track;

  const auto& particle_post { step.getParticlePost() };

  const auto& x1 { particle_post.getPosition() };
  const auto& u1 { particle_post.getDirection() };

  track.x1              << x1.getX()/1_cm, x1.getY()/1_cm, (x1.getZ()-r_earth_)/1_cm;
  track.u1              << u1.getX(), u1.getY(), u1.getZ();

  // If the particle is not heading for an interaction with the detector box (i.e.
  // if is outside the box and heading away from it) then we absorb the particle.
  if(not box_has_future_intersection(min_corner_, max_corner_, track.x0, track.u0)) {
    return ProcessReturn::ParticleAbsorbed;
  }
  
  if(!visitor_) return ProcessReturn::Ok;

  const auto& particle_pre { step.getParticlePre() };

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

  track.e1              = particle_post.getEnergy()/1_MeV;
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
  return max_step_length_;
}
 
namespace {

  // Workaround for UrQMD not having InteractionProcess as a base class
  class MyUrQMD: public corsika::urqmd::UrQMD, public InteractionProcess<MyUrQMD> {};

  class CORSIKA8ShowerGeneratorImpl: public CORSIKA8ShowerGenerator {
  public:
    using EnvironmentInterface = IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>;
    template <typename T> using MyExtraEnv = MediumPropertyModel<UniformMagneticField<T>>;
    // using EnvironmentInterface =
    //   IRefractiveIndexModel<IMagneticFieldModel<IMediumModel>>>;
    using EnvType = Environment<EnvironmentInterface>;
    using StackType = setup::Stack<EnvType>;
    using TrackingType = setup::Tracking;
    using Particle = StackType::particle_type;

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

    // struct IsTauSwitch {
    //   bool operator()(const Particle& p) const {
    //     return (p.getPID() == Code::TauMinus || p.getPID() == Code::TauPlus);
    //   }
    // };

    struct EMHadronSwitch {
      EMHadronSwitch() = default;
      bool operator()(const Particle& p) const { return is_hadron(p.getPID()); }
    };
  
    struct EnergySwitch {
      HEPEnergyType cutE_;
      EnergySwitch(HEPEnergyType cutE)
          : cutE_(cutE) {}
      bool operator()(const Particle& p) const { return (p.getKineticEnergy() < cutE_); }
    };

    EnvType env_;
    CoordinateSystemPtr root_cs_ = env_.getCoordinateSystem();
    Point center_ { root_cs_, 0_m, 0_m, 0_m };
    Point ground_; // set in constructor from value passed in config

    // Hadronic interaction
    using HEModelType = DynamicInteractionProcess<StackType>;
    using SibyllType = corsika::sibyll::Interaction;
    using QGSJetType = corsika::qgsjetII::Interaction;
    using EPOSType = corsika::epos::Interaction;
    std::shared_ptr<SibyllType> sibyll_;
    std::shared_ptr<HEModelType> he_model_;
    
    // Decay
    using PythiaType = corsika::pythia8::Decay;
    // using TauolaType = corsika::tauola::Decay;
    // using DecayType = SwitchProcessSequence<IsTauSwitch, corsika::tauola::Decay, corsika::pythia8::Decay>;
    using DecaySequenceType = PythiaType;
    std::shared_ptr<PythiaType> decay_pythia_;
    // std::shared_ptr<TauolaType> decay_tauola_;
    std::shared_ptr<DecaySequenceType> decay_sequence_;

    // Photo hadronic interactions
    using SophiaType = corsika::sophia::InteractionModel;
    using ParticleCutType = ParticleCut<>;
    using EMCascadeType = corsika::proposal::Interaction<SophiaType, corsika::sibyll::HadronInteractionModel> ;
    std::shared_ptr<SophiaType> sophia_;
    std::shared_ptr<ParticleCutType> particle_cut_;
    std::shared_ptr<EMCascadeType> em_cascade_;

    // Continuous losses
    using BetheBlochLossType = BetheBlochPDG<>;
    using ProposalLossType = corsika::proposal::ContinuousProcess<>;
    using ContinuousLossSequenceType = SwitchProcessSequence<EMHadronSwitch, BetheBlochLossType&, ProposalLossType& >;
    std::shared_ptr<ProposalLossType > em_continuous_proposal_;
    std::shared_ptr<BetheBlochLossType> em_continuous_bethe_;
    std::shared_ptr<ContinuousLossSequenceType> em_continuous_;

    // Low energy interactions
    using UrQMDType = MyUrQMD;
    std::shared_ptr<UrQMDType> le_int_model_;

    // Hadron sequence
    using HadronSequenceType = SwitchProcessSequence<EnergySwitch, MyUrQMD&, DynamicInteractionProcess<StackType>&>;
    std::shared_ptr<HadronSequenceType> hadron_sequence_;

    // Track handoff
    using TrackHandoffType = TrackHandoff;
    std::shared_ptr<TrackHandoffType> track_handoff_;

    // Final process sequence
    using FinalProcessSequenceType = 
      decltype(make_sequence(*hadron_sequence_, *decay_sequence_, *em_cascade_, *em_continuous_, 
                             *track_handoff_, *particle_cut_));
    std::shared_ptr<FinalProcessSequenceType> process_sequence_;
  };

} // anonymous namespace

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
      Medium::AirDry1Atm,
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

  LOG(INFO) << "Setting up SIBYLL";
  const auto all_elements = corsika::get_all_elements_in_universe(env_);
  // have SIBYLL always for PROPOSAL photo-hadronic interactions
  sibyll_ = std::make_shared<SibyllType>(all_elements, corsika::setup::C7trackedParticles);

  switch(config_.he_hadronic_model()) {
  case calin::ix::simulation::corsika8_shower_generator::SIBYLL:
  default:
    he_model_ = std::make_shared<HEModelType>(sibyll_);
    break;
  case calin::ix::simulation::corsika8_shower_generator::QGSJet:
    LOG(INFO) << "Setting up QGSJet";
    he_model_ = std::make_shared<HEModelType>(std::make_shared<QGSJetType>());
    break;
  case calin::ix::simulation::corsika8_shower_generator::EPOS_LHC:
    LOG(INFO) << "Setting up EPOS";
    he_model_ = std::make_shared<HEModelType>(
      std::make_shared<EPOSType>(corsika::setup::C7trackedParticles));
    break;    
  };

  // ==========================================================================
  // SETUP DECAY MODELS
  // ==========================================================================

  LOG(INFO) << "Setting up Pythia";
  decay_pythia_ = std::make_shared<PythiaType>();
  // decay_tauola_ = std::make_shared<TauolaType>(corsika::tauola::Helicity::LeftHanded);
  // decay_sequence_ = std::make_shared<DecaySequenceType>(IsTauSwitch(), decay_tauola_, decay_pythia_);
  decay_sequence_ = decay_pythia_;

  // ==========================================================================
  // Hadronic photon interactions in resonance region
  // ==========================================================================

  LOG(INFO) << "Setting up Sophia";
  sophia_ = std::make_shared<SophiaType>();
  
  // ==========================================================================
  // PARTICLE CUTS
  // ==========================================================================

  HEPEnergyType emcut  = config_.electrion_photon_cut() * 1_MeV;
  HEPEnergyType hadcut = config_.hadronic_cut() * 1_MeV;
  HEPEnergyType mucut  = config_.muon_cut() * 1_MeV;
  HEPEnergyType taucut = config_.tau_cut() * 1_MeV;
  LOG(INFO) << "Setting particle cuts (MeV): em=" << emcut/1_MeV << " had=" << hadcut/1_MeV << " mu=" << mucut/1_MeV << " tau=" << taucut/1_MeV;
  particle_cut_ = std::make_shared<ParticleCutType>(emcut, emcut, hadcut, mucut, taucut, true);

  // ==========================================================================
  // EM CASCADE
  // ==========================================================================

  LOG(INFO) << "Setting EM production thresholds";
  // tell proposal that we are interested in all energy losses above the particle cut
  auto prod_threshold = std::min({emcut, hadcut, mucut, taucut});
  set_energy_production_threshold(Code::Electron, prod_threshold);
  set_energy_production_threshold(Code::Positron, prod_threshold);
  set_energy_production_threshold(Code::Photon, prod_threshold);
  set_energy_production_threshold(Code::MuMinus, prod_threshold);
  set_energy_production_threshold(Code::MuPlus, prod_threshold);
  set_energy_production_threshold(Code::TauMinus, prod_threshold);
  set_energy_production_threshold(Code::TauPlus, prod_threshold);

  // energy threshold for high energy hadronic model. Affects LE/HE switch for
  // hadron interactions and the hadronic photon model in proposal
  HEPEnergyType const he_hadron_model_threshold = config_.he_hadronic_transition_energy() * 1_MeV;

  LOG(INFO) << "Setting up EM hadronic cascade";
  em_cascade_ = std::make_shared<EMCascadeType>(
    env_, *sophia_, sibyll_->getHadronInteractionModel(), he_hadron_model_threshold);

  // ==========================================================================
  // Continuous losses
  // ==========================================================================

  LOG(INFO) << "Setting up Bethe-Bloch continuous losses";
  em_continuous_bethe_ = std::make_shared<BetheBlochLossType>();
  LOG(INFO) << "Setting up Proposal continuous losses";
  em_continuous_proposal_ = std::make_shared<ProposalLossType>(env_);
  LOG(INFO) << "Setting up switched continuous cascade";
  em_continuous_ = std::make_shared<ContinuousLossSequenceType>(
      EMHadronSwitch(), *em_continuous_bethe_, *em_continuous_proposal_);

  // ==========================================================================
  // LOW ENERGY INTERACTIONS
  // ==========================================================================

  LOG(INFO) << "Setting up UrQMD";
  le_int_model_ = std::make_shared<UrQMDType>();
  
  // ==========================================================================
  // HADRON SEQUENCE
  // ==========================================================================

  LOG(INFO) << "Setting up hadron sequence";
  hadron_sequence_ = std::make_shared<HadronSequenceType>(
      EnergySwitch(he_hadron_model_threshold), *le_int_model_, *he_model_);

  // ==========================================================================
  // TRACK HANDOFF
  // ==========================================================================

  LOG(INFO) << "Setting up track handoff";
  track_handoff_ = std::make_shared<TrackHandoff>(config_.earth_radius(),
    config_.zground(), config_.ztop(), config_.detector_box_side());

  // ==========================================================================
  // FINAL PROCESS SEQUENCE
  // ==========================================================================

  LOG(INFO) << "Setting up process sequence";
  process_sequence_ = std::make_shared<FinalProcessSequenceType>(
    make_sequence(*hadron_sequence_, *decay_sequence_, *em_cascade_, *em_continuous_, 
      *track_handoff_, *particle_cut_));
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

CORSIKA8ShowerGenerator* CORSIKA8ShowerGenerator::new_instance(const config_type& config)
{
  return new CORSIKA8ShowerGeneratorImpl(config);
}
