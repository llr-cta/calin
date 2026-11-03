# calin/python/simulation/iact_factory.py -- Stephen Fegan -- 2026-09-29
#
# Unified factory functions for setting up IACT arrays, optics, atmosphere,
# magnetic field, detector efficiencies, and Geant4 generators across simulation scripts.
#
# Copyright 2026, Stephen Fegan <sfegan@llr.in2p3.fr>
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License version 2 or later, as published by
# the Free Software Foundation.
#
# "calin" is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

import numpy

import calin.math.geometry
import calin.simulation.atmosphere
import calin.simulation.detector_efficiency
import calin.simulation.geant4_shower_generator
import calin.simulation.ray_processor
import calin.simulation.tracker
import calin.simulation.vcl_iact
import calin.simulation.vs_cta
import calin.simulation.vs_optics
import calin.simulation.vs_panoseti
import calin.simulation.world_magnetic_model
import calin.ix.simulation.vcl_iact
import calin.ix.simulation.simulated_event

class SiteEnvironment:
    """Container holding site-level environment, atmosphere, and detector responses."""
    def __init__(self, site, zobs, atm, atm_abs, det_eff, cone_eff, pe_gen,
                 lens_refractive_index_spline, wmm, bfield, array_origin, array_config_fn):
        self.site = site
        self.zobs = zobs
        self.atm = atm
        self.atm_abs = atm_abs
        self.det_eff = det_eff
        self.cone_eff = cone_eff
        self.pe_gen = pe_gen
        self.lens_refractive_index_spline = lens_refractive_index_spline
        self.wmm = wmm
        self.bfield = bfield
        self.array_origin = array_origin
        self.array_config_fn = array_config_fn

def load_site_environment(site = 'ctan', enable_pe_spectrum = False,
                          no_bfield = False, quiet = True, avx = 512):
    """
    Load site atmosphere, atmospheric absorption, observation level, detector efficiencies,
    and geomagnetic field.
    """
    if site == 'dark100':
        zobs = calin.simulation.vs_panoseti.palomar_observation_level()
        atm = calin.simulation.vs_panoseti.palomar_atmosphere(quiet=quiet)
        atm_abs = calin.simulation.vs_panoseti.palomar_atmospheric_absorption(quiet=quiet)
        dark100_datapack = calin.simulation.vs_panoseti.read_optical_model_datapack()
        array_config_fn = calin.simulation.vs_panoseti.dark100_palomar_config
        det_eff = calin.simulation.vs_panoseti.detection_efficiency_from_datapack(dark100_datapack, quiet=quiet)
        cone_eff = None
        pe_gen = None
        lens_refractive_index_spline = calin.simulation.vs_panoseti.lens_refractive_index_spline_from_datapack(dark100_datapack)
        dummy_array_origin = array_config_fn(elevation=0).array_origin()
    elif site == 'ctan':
        zobs = calin.simulation.vs_cta.ctan_observation_level()
        atm = calin.simulation.vs_cta.ctan_atmosphere(quiet=quiet)
        atm_abs = calin.simulation.vs_cta.ctan_atmospheric_absorption(quiet=quiet)
        array_config_fn = calin.simulation.vs_cta.mstn1_config
        det_eff = calin.simulation.vs_cta.mstn_detection_efficiency(quiet=quiet)
        cone_eff = calin.simulation.vs_cta.mstn_cone_efficiency(quiet=quiet)
        pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(
            avx=avx, quiet=quiet) if enable_pe_spectrum else None
        lens_refractive_index_spline = None
        dummy_array_origin = array_config_fn(elevation=0).array_origin()
    elif site == 'ctas':
        zobs = calin.simulation.vs_cta.ctas_observation_level()
        atm = calin.simulation.vs_cta.ctas_atmosphere(quiet=quiet)
        atm_abs = calin.simulation.vs_cta.ctas_atmospheric_absorption(quiet=quiet)
        array_config_fn = calin.simulation.vs_cta.msts1_config
        det_eff = calin.simulation.vs_cta.mstn_detection_efficiency(quiet=quiet)
        cone_eff = calin.simulation.vs_cta.mstn_cone_efficiency(quiet=quiet)
        pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(
            avx=avx, quiet=quiet) if enable_pe_spectrum else None
        lens_refractive_index_spline = None
        dummy_array_origin = array_config_fn(elevation=0).array_origin()
    else:
        raise ValueError(f'Unknown site: {site}')

    if no_bfield:
        wmm = None
        bfield = None
    else:
        wmm = calin.simulation.world_magnetic_model.WMM()
        bfield = wmm.field_vs_elevation(dummy_array_origin.latitude(), dummy_array_origin.longitude())

    return SiteEnvironment(
        site=site,
        zobs=zobs,
        atm=atm,
        atm_abs=atm_abs,
        det_eff=det_eff,
        cone_eff=cone_eff,
        pe_gen=pe_gen,
        lens_refractive_index_spline=lens_refractive_index_spline,
        wmm=wmm,
        bfield=bfield,
        array_origin=dummy_array_origin,
        array_config_fn=array_config_fn
    )

def setup_telescope_array(site, el_deg = 60.0, zero_scope_0 = True):
    """
    Construct telescope array parameters for the given elevation and resolve
    channel and telescope counts.
    """
    if site == 'dark100':
        dark100 = calin.simulation.vs_panoseti.dark100_palomar_config(elevation=el_deg)
        # Note: keep true positions for PANOSETI telescopes
        nscope = dark100.scope_positions_size()
        nchan = dark100.num_pixels_per_axis() ** 2
        return dark100, nscope, nchan, 'PANOSETI/Dark100'
    elif site in ('ctan', 'ctas'):
        mst_config = calin.simulation.vs_cta.mstn1_config if site == 'ctan' else calin.simulation.vs_cta.msts1_config
        mst = mst_config(elevation=el_deg)
        if zero_scope_0 and mst.prescribed_array_layout().scope_positions_size() > 0:
            mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_x(0)
            mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_y(0)
        nscope = mst.prescribed_array_layout().scope_positions_size()
        array = calin.simulation.vs_optics.VSOArray()
        array.generateFromArrayParameters(mst)
        scope = array.telescope(0)
        telescope_layout = scope.convert_to_telescope_layout()
        nchan = telescope_layout.camera().channel_size()
        return mst, nscope, nchan, 'MST/NC'
    else:
        raise ValueError(f'Unknown site: {site}')

def create_iact_array(atm, atm_abs, avx = 512, no_refraction = False):
    """
    Select simulation class based on AVX architecture and instantiate VCLIACTArray.
    """
    if avx == 128:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray128
    elif avx == 256:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray256
    else:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray512

    iact_cfg = iact_class.default_config()
    if no_refraction:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_NO_RAYS)
    else:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_ONLY_CLOSE_RAYS)

    iact = iact_class(atm, atm_abs, iact_cfg)
    return iact, iact_cfg, iact_class

def attach_iact_propagators(iact, site, array_params, bmax_polynomial,
                            reuse = 1, nscope = 1, nchan = 0,
                            det_eff = None, cone_eff = None, pe_gen = None,
                            tts = 0.0, lens_spline = None,
                            detector_type_name = 'MST/NC',
                            pe_processor_factory = None):
    """
    Attach propagator sets and ray propagators (Davies-Cotton or PANOSETI).

    ``bmax_polynomial`` is in metres with coefficients in ascending order,
    matching the command-line polynomial convention used by the simulation scripts.
    """
    if pe_processor_factory is None:
        pe_processor_factory = lambda s, c: calin.simulation.ray_processor.SimpleListPEProcessor(s, c)

    # Inputs are in metres, with polynomial coefficients in ascending order.
    # IACT stores coefficients in numpy.polyval order and distances in cm.
    bmax_coefficients_m = numpy.atleast_1d(numpy.asarray(bmax_polynomial, dtype=float))
    if bmax_coefficients_m.size == 0:
        bmax_poly_cm = numpy.asarray([0.0])
    else:
        bmax_poly_cm = numpy.flipud(bmax_coefficients_m) * 100.0

    all_pe_processor = []
    all_prop = []
    n_reuse = max(1, reuse)

    for i in range(n_reuse):
        iact.add_propagator_set(numpy.flipud(bmax_poly_cm), f"Super array {i}")
        pe_processor = pe_processor_factory(nscope, nchan)
        if site == 'dark100':
            prop = iact.add_panoseti_propagator(
                array_params, lens_spline, pe_processor, det_eff,
                detector_type_name, pe_gen, tts)
        else:
            prop = iact.add_davies_cotton_propagator(
                array_params, pe_processor, det_eff, cone_eff, pe_gen, tts, detector_type_name)
        all_pe_processor.append(pe_processor)
        all_prop.append(prop)

    return all_pe_processor, all_prop

def set_iact_pointing(iact, el_deg, az_deg, theta_deg = 0.0, phi_deg = 0.0,
                      apply_viewcone_cut = True):
    """
    Set pointing direction for all telescopes, apply optional viewcone cut,
    and compute the primary viewcone direction vector.
    """
    iact.point_all_telescopes_az_el_deg(az_deg, el_deg)
    el = el_deg * numpy.pi / 180.0
    az = az_deg * numpy.pi / 180.0
    pt_dir = numpy.asarray([
        numpy.cos(el) * numpy.sin(az),
        numpy.cos(el) * numpy.cos(az),
        numpy.sin(el)
    ])
    theta = theta_deg * numpy.pi / 180.0
    phi   = phi_deg   * numpy.pi / 180.0
    vc_dir = numpy.asarray([
        numpy.sin(theta) * numpy.cos(phi),
        numpy.sin(theta) * numpy.sin(phi),
        numpy.cos(theta)
    ])
    vc_dir = calin.math.geometry.rotate_vec_z_to_u_Rzy(vc_dir, -pt_dir)

    if apply_viewcone_cut:
        iact.set_viewcone_from_telescope_fields_of_view()

    return pt_dir, vc_dir

def create_geant4_generator(atm, bfield = None, multiple_scattering = 'normal', primary = 'gamma'):
    """
    Configure and instantiate the Geant4 shower generator with requested MSC parameters.
    """
    geant4_cfg = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator.customized_config(
        1000, 0, atm.top_of_atmosphere(),
        calin.simulation.geant4_shower_generator.VerbosityLevel_SUPRESSED_STDOUT)

    if multiple_scattering == 'minimal':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit Minimal')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad Minimal')
    elif multiple_scattering == 'simple':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseSafety')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseSafety')
    elif multiple_scattering == 'normal':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
    elif multiple_scattering == 'better':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactor 0.01')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactorMuHad 0.01')
    elif multiple_scattering == 'insane':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactor 0.001')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactorMuHad 0.001')

    if primary == 'iron':
        geant4_cfg.set_enable_ions(True)

    generator = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator(atm, geant4_cfg, bfield)
    return generator, geant4_cfg

def get_tracker_particle_type(primary):
    """Map primary particle name string to calin tracker ParticleType enum."""
    particle_map = {
        'gamma':    calin.simulation.tracker.ParticleType_GAMMA,
        'muon':     calin.simulation.tracker.ParticleType_MUON,
        'electron': calin.simulation.tracker.ParticleType_ELECTRON,
        'proton':   calin.simulation.tracker.ParticleType_PROTON,
        'helium':   calin.simulation.tracker.ParticleType_HELIUM,
        'iron':     calin.simulation.tracker.ParticleType_IRON,
    }
    if primary not in particle_map:
        raise ValueError(f'Unknown primary particle type: {primary}')
    return particle_map[primary]

def get_simulated_event_particle_type(primary):
    """Map primary particle name string to simulated_event protobuf enum."""
    particle_map = {
        'gamma':    calin.ix.simulation.simulated_event.GAMMA,
        'muon':     calin.ix.simulation.simulated_event.MUON,
        'electron': calin.ix.simulation.simulated_event.ELECTRON,
        'proton':   calin.ix.simulation.simulated_event.PROTON,
        'helium':   calin.ix.simulation.simulated_event.HELIUM,
        'iron':     calin.ix.simulation.simulated_event.IRON,
    }
    if primary not in particle_map:
        raise ValueError(f'Unknown primary particle type: {primary}')
    return particle_map[primary]
