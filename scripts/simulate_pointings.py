# calin/scripts/sim_pointings.py -- Stephen Fegan - 2026-03-16
#
# Simulate events from given spectrum with GEANT4 accounting for variable
# scattering radius and viewcone. Write events to HDF5 files, one or more per
# pointing, iterating over a grid of pointings supplied by the user-provided
# pointing grid functions "num_grid_pointings" and "get_grid_pointing".
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

import os
import argparse
import signal
import concurrent.futures
import datetime
import time
import numpy
import scipy.interpolate
import calin.math.geometry
import calin.math.healpix_array
import calin.ix.simulation.vcl_iact
import calin.simulation.tracker
import calin.simulation.vs_cta
import calin.simulation.vs_optics
import calin.simulation.ray_processor
import calin.simulation.world_magnetic_model
import calin.simulation.geant4_shower_generator
import calin.simulation.vcl_iact
import calin.ix.simulation.simulated_event
import calin.iact_data.instrument_layout
import calin.iact_data.nectarcam_layout
import calin.provenance.anthology

# =============================================================================
# POINTING GRID INTERFACE
#
# Replace these two functions with your healpix (or other) grid implementation.
# - num_grid_pointings(nside) : returns the total number of pointings for the given
#                          healpix nside parameter, ignoring any zn cut
# - get_grid_pointing(nside,i): returns (index, zn_deg, az_deg) for the i-th pixel
#                          in the grid; index == i for a healpix grid
#
# num_grid_pointings() is used only to determine the zero-padding width for output
# filenames, so the mapping index<->pointing is stable regardless of grid_znmax.
# =============================================================================

def num_grid_pointings(nside):
    """Return the total number of pointings for the given healpix nside."""
    return calin.math.healpix_array.npixel(nside)

def get_grid_pointing(nside, i):
    """Return (index, zn_deg, az_deg) for the i-th pixel of the healpix grid."""
    zn_rad, az_rad = calin.math.healpix_array.pixid_to_ang(nside, i)
    return (i, numpy.rad2deg(zn_rad), numpy.rad2deg(az_rad))

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def comma_separated_floats(value):
    try:
        return [float(x) for x in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a comma-separated list of floats (e.g. 1.0,2.3,4.2)")

def zn_polynomial_list(value):
    """
    Parse a zenith-angle keyed polynomial entry of the form 'zn:P0,P1,P2,...'
    or a bare polynomial 'P0,P1,P2,...' (which is treated as zn=0).

    Returns (zn_deg, [P0, P1, ...]).

    Multiple --bmax_polynomial / --viewcone_polynomial options may be given on
    the command line; argparse collects them into a list via action='append'.
    At runtime the entry whose zenith angle is closest to the pointing zenith
    angle is selected (nearest-neighbour lookup).

    Examples
    --------
      --bmax_polynomial 1000.0
      --bmax_polynomial 0:1000.0
      --bmax_polynomial 20:1200.0,50.0
      --bmax_polynomial 0:800.0  --bmax_polynomial 20:1200.0,50.0  --bmax_polynomial 40:1500.0,80.0,3.0
    """
    if ':' in value:
        zn_str, poly_str = value.split(':', 1)
        try:
            zn = float(zn_str)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Zenith angle before ':' must be a float, got '{zn_str}'")
    else:
        zn = 0.0
        poly_str = value
    try:
        poly = [float(x) for x in poly_str.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Polynomial coefficients must be comma-separated floats, got '{poly_str}'")
    return (zn, poly)

def zn_interpolation_metric(zn_deg):
    """
    Map a zenith angle in degrees to the scalar used as the interpolation axis
    when blending polynomial entries from the zn-keyed polynomial lists.

    The default is cos(zn), which gives uniform weighting in solid angle and
    matches the natural dependence of shower physics on atmospheric depth
    (~1/cos(zn)).  Replace this function to change the interpolation metric
    without touching any other code.
    """
    return numpy.cos(numpy.deg2rad(zn_deg))

def select_polynomial_for_zn(zn_poly_list, zn_deg):
    """
    Given a list of (zn_deg, coefficients) pairs (as parsed by zn_polynomial_list),
    return a coefficient list interpolated linearly in zn_interpolation_metric()
    between the two bracketing entries.

    Edge behaviour:
      - If zn_deg is below the lowest entry's zenith angle, the lowest entry
        is returned unchanged (no extrapolation).
      - If zn_deg is above the highest entry's zenith angle, the highest entry
        is returned unchanged (no extrapolation).
      - If only one entry is present, it is returned as-is.

    The interpolated polynomial has the same degree as the longer of the two
    bracketing polynomials (the shorter one is zero-padded on the high-degree
    end before blending).
    """
    if len(zn_poly_list) == 1:
        return zn_poly_list[0][1]

    # Build (metric_value, coefficients) pairs and sort ascending by metric.
    # cos(zn) decreases as zn increases, so ascending metric = ascending zn.
    metric = [(zn_interpolation_metric(zn), poly) for zn, poly in zn_poly_list]
    metric.sort(key=lambda t: t[0])

    t_target = zn_interpolation_metric(zn_deg)

    # Clamp to range -- no extrapolation beyond the supplied entries
    if t_target <= metric[0][0]:
        return metric[0][1]
    if t_target >= metric[-1][0]:
        return metric[-1][1]

    # Find the bracketing pair and interpolate
    for i in range(len(metric) - 1):
        t0, p0 = metric[i]
        t1, p1 = metric[i + 1]
        if t0 <= t_target <= t1:
            alpha = (t_target - t0) / (t1 - t0)   # 0 at t0, 1 at t1
            # Zero-pad the shorter polynomial to equal length before blending
            n = max(len(p0), len(p1))
            a0 = numpy.array(p0 + [0.0] * (n - len(p0)))
            a1 = numpy.array(p1 + [0.0] * (n - len(p1)))
            return list((1.0 - alpha) * a0 + alpha * a1)

    return metric[-1][1]   # unreachable, but satisfies linters

def handle_sigint_quietly(signum, frame):
    global stop_requested
    stop_requested = True

def handle_sigint(signum, frame):
    global stop_requested
    if not stop_requested:
        print("\nStop requested - terminating after current pointings complete")
    stop_requested = True

def format_energy(etev):
    if etev < 0.1:
        return f'{etev*1000:,.2f} GeV'
    elif etev < 1.0:
        return f'{etev*1000:,.1f} GeV'
    elif etev < 10:
        return f'{etev:,.3f} TeV'
    elif etev < 100:
        return f'{etev:,.2f} TeV'
    else:
        return f'{etev:,.1f} TeV'

def format_duration(seconds: int) -> str:
    seconds = int(seconds)
    minutes = seconds // 60
    hours = minutes // 60
    days = hours // 24
    if seconds < 3600:
        return f"{minutes}m{seconds % 60:02d}"
    elif seconds < 86400:
        return f"{hours}h{minutes % 60:02d}"
    else:
        return f"{days}d{hours % 24:02d}h{minutes % 60:02d}"

def pointing_filename(template, pointing_index, seq_num, pointing_index_width):
    """
    Format an output filename. Supported tokens:
      {id}       - random 128-bit hex identifier
      {seq}      - sequential file number within this pointing (0-based)
      {pointing} - zero-padded pointing index
    """
    rand_id = f'{numpy.random.bit_generator.randbits(128):032X}'
    pointing_str = str(pointing_index).zfill(pointing_index_width)
    return template.format(id=rand_id, seq=seq_num, pointing=pointing_str)

# =============================================================================
# WORKER INITIALISATION
#
# Sets up everything that is constant across all pointings: atmosphere, optics,
# detector efficiency, spectrum transform, Geant4 generator, magnetic field.
# The IACT array itself is NOT constructed here because it encodes the pointing
# direction; it is (re-)constructed inside gen_batch() for each pointing.
# =============================================================================

def init(args):
    global saved_args
    saved_args = args

    numpy.random.seed()

    global instance_id
    global instance_event_id
    instance_id = numpy.random.bit_generator.randbits(64)
    instance_event_id = 0

    # Select simulation classes based on AVX size requested
    global iact_class
    if args.avx == 128:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray128
    elif args.avx == 256:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray256
    else:
        iact_class = calin.simulation.vcl_iact.VCLIACTArray512

    # Load site-specific atmosphere, observation level and array layout
    global zobs
    global atm
    global atm_abs
    global mst_config
    if args.site == 'ctan':
        zobs = calin.simulation.vs_cta.ctan_observation_level()
        atm = calin.simulation.vs_cta.ctan_atmosphere(quiet=True)
        atm_abs = calin.simulation.vs_cta.ctan_atmospheric_absorption(quiet=True)
        mst_config = calin.simulation.vs_cta.mstn1_config
    else:
        zobs = calin.simulation.vs_cta.ctas_observation_level()
        atm = calin.simulation.vs_cta.ctas_atmosphere(quiet=True)
        atm_abs = calin.simulation.vs_cta.ctas_atmospheric_absorption(quiet=True)
        mst_config = calin.simulation.vs_cta.msts1_config
    mst = mst_config(elevation=0)  # dummy config to get array origin for bfield

    # IACT config object (reused for every pointing, only instantiation deferred)
    global iact_cfg
    iact_cfg = iact_class.default_config()
    if args.no_refraction:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_NO_RAYS)
    else:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_ONLY_CLOSE_RAYS)

    # Load detector and efficiency models and the SPE generator
    global det_eff
    global cone_eff
    global pe_gen
    global store_pe_weights
    det_eff = calin.simulation.vs_cta.mstn_detection_efficiency(quiet=True)
    cone_eff = calin.simulation.vs_cta.mstn_cone_efficiency(quiet=True)
    pe_gen = None
    store_pe_weights = False
    if args.enable_pe_spectrum:
        pe_gen = calin.simulation.vs_cta.mstn_spe_amplitude_generator(quiet=True)
        store_pe_weights = True

    global store_times_as_integer
    store_times_as_integer = args.store_times_as_integer

    # bmax and viewcone polynomials are resolved per-pointing in setup_pointing();
    # nothing to compute here.

    # Magnetic field (if not disabled) — depends only on site, not pointing
    global wmm
    global bfield
    if args.no_bfield:
        bfield = None
    else:
        wmm = calin.simulation.world_magnetic_model.WMM()
        bfield = wmm.field_vs_elevation(mst.array_origin().latitude(), mst.array_origin().longitude())

    # Configure Geant4 shower generator
    geant4_cfg = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator.customized_config(
        1000, 0, atm.top_of_atmosphere(),
        calin.simulation.geant4_shower_generator.VerbosityLevel_SUPRESSED_STDOUT)
    if args.multiple_scattering == 'minimal':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit Minimal')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad Minimal')
    elif args.multiple_scattering == 'simple':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseSafety')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseSafety')
    elif args.multiple_scattering == 'normal':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
    elif args.multiple_scattering == 'better':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactor 0.01')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactorMuHad 0.01')
    elif args.multiple_scattering == 'insane':
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimit UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/StepLimitMuHad UseDistanceToBoundary')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactor 0.001')
        geant4_cfg.add_pre_init_commands('/process/msc/RangeFactorMuHad 0.001')
    if args.primary == 'iron':
        geant4_cfg.set_enable_ions(True)

    global generator
    global saved_geant4_cfg
    saved_geant4_cfg = geant4_cfg
    generator = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator(atm, geant4_cfg, bfield)

    # Particle type
    global particle_type
    if args.primary == 'gamma':
        particle_type = calin.simulation.tracker.ParticleType_GAMMA
    elif args.primary == 'muon':
        particle_type = calin.simulation.tracker.ParticleType_MUON
    elif args.primary == 'electron':
        particle_type = calin.simulation.tracker.ParticleType_ELECTRON
    elif args.primary == 'proton':
        particle_type = calin.simulation.tracker.ParticleType_PROTON
    elif args.primary == 'helium':
        particle_type = calin.simulation.tracker.ParticleType_HELIUM
    elif args.primary == 'iron':
        particle_type = calin.simulation.tracker.ParticleType_IRON
    else:
        raise ValueError(f'Unknown primary particle type: {args.primary}')

    # The spectral transform depends on the per-pointing bmax and viewcone
    # polynomials, so it is computed inside setup_pointing() for each pointing.

    signal.signal(signal.SIGINT, handle_sigint_quietly)

# =============================================================================
# POINTING SETUP  (called once per pointing, inside each worker)
#
# Constructs a fresh IACT instance pointed in the given direction, resolves the
# bmax and viewcone polynomials for this zenith angle, builds the spectral
# transform, and returns everything needed by gen_event() / gen_batch().
# =============================================================================

def make_bmax_polynomial(zn_deg):
    # Return bmax polynomial (numpy.polyval order, internal units cm) for zn_deg.
    raw = select_polynomial_for_zn(saved_args.bmax_polynomial, zn_deg)
    if len(raw) == 0:
        return numpy.asarray([0.0])
    return numpy.flipud(raw) * 100.0   # m -> cm

def make_viewcone_polynomial(zn_deg):
    # Return viewcone polynomial (numpy.polyval order, radians) for zn_deg.
    raw = select_polynomial_for_zn(saved_args.viewcone_polynomial, zn_deg)
    if len(raw) == 0:
        return numpy.asarray([0.0])
    return numpy.flipud(raw) * numpy.pi / 180.0   # deg -> rad

def make_spectral_transform(bmax_polynomial, viewcone_polynomial):
    # Build CDF-inversion spectral transform for the given bmax and viewcone
    # polynomials. Returns callable: uniform deviate -> log10(E/1TeV).
    args = saved_args
    if args.emax > args.emin:
        if len(args.spectral_polynomial) == 0:
            sp = numpy.asarray([1.0, 0.0])
        else:
            sp = numpy.flipud(args.spectral_polynomial).copy()
            sp[-1] += 1.0
            sp = numpy.append(sp, [0.0])
        xmin = numpy.log10(args.emin)
        xmax = numpy.log10(args.emax)
        x = numpy.linspace(xmin, xmax, 10000)
        y = 10 ** (numpy.polyval(sp, x))
        if len(bmax_polynomial) > 1:
            y *= numpy.polyval(bmax_polynomial, x) ** 2
        if len(viewcone_polynomial) > 1:
            y *= 1 - numpy.cos(numpy.polyval(viewcone_polynomial, x))
        p = numpy.append([0.0], numpy.cumsum(y[1:] + y[:-1]) / 2)
        p /= p[-1]
        if 0 < args.spectral_constant <= 1.0:
            p = (p * (1 - args.spectral_constant)
                 + numpy.linspace(0.0, 1.0, len(p)) * args.spectral_constant)
        return scipy.interpolate.interp1d(
            p, x, kind='linear', bounds_error=False, fill_value=(xmin, xmax))
    else:
        lx = numpy.log10(args.emin)
        return lambda r: lx

def setup_pointing(el_deg, az_deg):
    # Generate array parameters
    mst = mst_config(elevation = el_deg)
    mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_x(0)
    mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_y(0)
    nscope = mst.prescribed_array_layout().scope_positions_size()

    # Instantiate array
    global array
    array = calin.simulation.vs_optics.VSOArray()
    array.generateFromArrayParameters(mst)
    scope = array.telescope(0)
    telescope_layout = scope.convert_to_telescope_layout()
    nchan = telescope_layout.camera().channel_size()

    # Build and return (iact, vc_dir, bmax_poly, viewcone_poly,
    # spectral_tf, sim_config) for the pointing (el_deg, az_deg).
    zn_deg = 90.0 - el_deg
    detector_type_name = 'MST/NC'

    # Resolve polynomials for this zenith angle
    bmax_poly     = make_bmax_polynomial(zn_deg)
    viewcone_poly = make_viewcone_polynomial(zn_deg)
    spectral_tf   = make_spectral_transform(bmax_poly, viewcone_poly)

    # Fresh IACT instance for this pointing
    iact = iact_class(atm, atm_abs, iact_cfg)

    global all_pe_processor
    global all_prop
    all_pe_processor = []
    all_prop = []
    for i in range(numpy.max([1, saved_args.reuse])):
        iact.add_propagator_set(numpy.flipud(bmax_poly), f"Super array {i}")
        pe_processor = calin.simulation.ray_processor.SimpleListPEProcessor(nscope, nchan)
        prop = iact.add_davies_cotton_propagator(
            mst, pe_processor, det_eff, cone_eff, pe_gen, saved_args.tts, detector_type_name)
        all_pe_processor.append(pe_processor)
        all_prop.append(prop)


    # Pointing direction and viewcone axis
    iact.point_all_telescopes_az_el_deg(az_deg, el_deg)
    el = el_deg * numpy.pi / 180.0
    az = az_deg * numpy.pi / 180.0
    pt_dir = numpy.asarray([
        numpy.cos(el) * numpy.sin(az),
        numpy.cos(el) * numpy.cos(az),
        numpy.sin(el)
    ])
    theta = saved_args.theta * numpy.pi / 180
    phi   = saved_args.phi   * numpy.pi / 180
    vc_dir = numpy.asarray([
        numpy.sin(theta) * numpy.cos(phi),
        numpy.sin(theta) * numpy.sin(phi),
        numpy.cos(theta)
    ])
    vc_dir = calin.math.geometry.rotate_vec_z_to_u_Rzy(vc_dir, -pt_dir)

    if not saved_args.no_viewcone_cut:
        iact.set_viewcone_from_telescope_fields_of_view()

    # Simulation configuration protobuf for this pointing
    sim_config = calin.ix.simulation.simulated_event.SimulationConfiguration()
    particle_map = {
        'gamma':    calin.ix.simulation.simulated_event.GAMMA,
        'muon':     calin.ix.simulation.simulated_event.MUON,
        'electron': calin.ix.simulation.simulated_event.ELECTRON,
        'proton':   calin.ix.simulation.simulated_event.PROTON,
        'helium':   calin.ix.simulation.simulated_event.HELIUM,
        'iron':     calin.ix.simulation.simulated_event.IRON,
    }
    sim_config.set_particle_type(particle_map[saved_args.primary])
    sim_config.set_energy_lo(saved_args.emin)
    if saved_args.emax > saved_args.emin:
        sim_config.set_energy_hi(saved_args.emax)
        sim_config.set_energy_spectrum_polynomial(saved_args.spectral_polynomial)
    else:
        sim_config.set_energy_hi(saved_args.emin)
    sim_config.set_elevation(el_deg)
    sim_config.set_azimuth(az_deg)
    sim_config.set_theta(saved_args.theta)
    sim_config.set_phi(saved_args.phi)
    sim_config.mutable_primary_axis().set_x(vc_dir[0])
    sim_config.mutable_primary_axis().set_y(vc_dir[1])
    sim_config.mutable_primary_axis().set_z(vc_dir[2])
    sim_config.set_viewcone_halfangle_polynomial(
        numpy.flipud(viewcone_poly) * 180.0 / numpy.pi)
    sim_config.set_scattering_radius_polynomial(
        numpy.flipud(bmax_poly) * 0.01)
    sim_config.set_banner(get_banner(el_deg=el_deg, az_deg=az_deg,
                                     bmax_poly=bmax_poly,
                                     viewcone_poly=viewcone_poly,
                                     spectral_tf=spectral_tf))

    atm_abs_zmin = numpy.min(atm_abs.levels_cm())
    atm_abs_zmax = numpy.max(atm_abs.levels_cm())
    for level in atm.get_levels():
        proto_level = sim_config.add_atmospheric_level()
        proto_level.set_altitude(level.z)
        proto_level.set_thickness(level.t)
        proto_level.set_density(level.rho)
        proto_level.set_n_minus_one(level.nmo)
        if bfield is not None:
            b = bfield.field_nT(level.z)
            proto_level.mutable_bfield().set_x(b[0])
            proto_level.mutable_bfield().set_y(b[1])
            proto_level.mutable_bfield().set_z(b[2])
        if atm_abs_zmin < level.z < atm_abs_zmax:
            proto_level.set_optical_depth_1d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 1.5))
            proto_level.set_optical_depth_2d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.0))
            proto_level.set_optical_depth_2d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.5))
            proto_level.set_optical_depth_3d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.0))
            proto_level.set_optical_depth_3d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.5))
            proto_level.set_optical_depth_4d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 4.0))

    sim_config.mutable_geant4_shower_generator_config().CopyFrom(saved_geant4_cfg)
    sim_config.mutable_iact_array_config().CopyFrom(iact_cfg)

    for ipropagatorset in range(iact.num_propagator_sets()):
        array_config = sim_config.add_detector_array_config()
        array_config.set_name(iact.propagator_set_name(ipropagatorset))
        array_config.set_scattering_radius_polynomial(
            iact.scattering_radius_polynomial(ipropagatorset))
        for ipropagator in range(iact.propagator_set_size(ipropagatorset)):
            propagator = iact.propagator_set_dc_element(ipropagatorset, ipropagator)
            dc_array = None
            if propagator:
                dc_array = propagator.array()
            else:
                propagator = iact.propagator_set_element(ipropagatorset, ipropagator)
            detector_spheres = propagator.detector_spheres()
            group_config = array_config.add_detector_group_config()
            group_config.set_type_name(
                iact.propagator_set_element_name(ipropagatorset, ipropagator))
            group_config.set_ndetector(len(detector_spheres))
            group_config.set_propagator_banner(propagator.banner())
            if dc_array:
                dc_array.dump_as_proto(group_config.mutable_dc_array_config())

    detector_type_config = sim_config.mutable_detector_type_config(detector_type_name)
    detector_type_config.set_type_name(detector_type_name)
    detector_type_config.set_detector_efficiency_banner(det_eff.banner())
    detector_type_config.set_detector_efficiency_energy(numpy.asarray(det_eff.all_xi()))
    detector_type_config.set_detector_efficiency_efficiency(numpy.asarray(det_eff.all_yi()))
    detector_type_config.set_angular_response_banner(cone_eff.banner())
    detector_type_config.set_angular_response_costheta(numpy.asarray(det_eff.all_xi()))
    detector_type_config.set_angular_response_efficiency(numpy.asarray(det_eff.all_yi()))
    if pe_gen:
        detector_type_config.set_pe_spectrum_banner(pe_gen.banner())
        detector_type_config.set_pe_spectrum_banner_q(pe_gen.raw_q())
        detector_type_config.set_pe_spectrum_banner_dp_dq(pe_gen.raw_dp_dq())
    detector_type_config.set_pe_time_spread(saved_args.tts)
    detector_type_config.mutable_dc_array_parameters().CopyFrom(mst)

    args_dict = vars(saved_args)
    for arg in args_dict:
        sim_config.set_command_line_args(arg, str(args_dict[arg]))

    return iact, vc_dir, bmax_poly, viewcone_poly, spectral_tf, sim_config

# =============================================================================
# EVENT AND BATCH GENERATION
# =============================================================================

def gen_event(iact, vc_dir, viewcone_polynomial, spectral_transform):
    # Generate a single shower event and return a SimulatedEvent proto.
    x = spectral_transform(numpy.random.uniform())   # x = log10(E/1TeV)
    e = 10 ** (x + 6)                               # convert TeV to MeV

    vc_halfangle = numpy.polyval(viewcone_polynomial, x)
    costheta = 1.0 - (1.0 - numpy.cos(vc_halfangle)) * numpy.random.uniform()
    theta = numpy.arccos(costheta)
    phi   = numpy.random.uniform() * 2 * numpy.pi
    u = numpy.asarray([
        numpy.sin(theta) * numpy.cos(phi),
        numpy.sin(theta) * numpy.sin(phi),
        numpy.cos(theta)
    ])
    u = calin.math.geometry.rotate_vec_z_to_u_Rzy(u, vc_dir)

    ct0 = 1.0 / u[2] * (atm.top_of_atmosphere() - atm.zobs(0))
    x0  = numpy.asarray([0, 0, atm.zobs(0)]) + ct0 * u

    generator.generate_showers(iact, 1, particle_type, e, x0, u, ct0)

    sim_event = calin.ix.simulation.simulated_event.SimulatedEvent()
    iact.save_to_simulated_event(sim_event, store_pe_weights, store_times_as_integer)

    sim_event.mutable_viewcone_axis().set_x(vc_dir[0])
    sim_event.mutable_viewcone_axis().set_y(vc_dir[1])
    sim_event.mutable_viewcone_axis().set_z(vc_dir[2])
    sim_event.set_viewcone_costheta(costheta)

    return sim_event


def gen_batch(filename, pointing_index, el_deg, az_deg, nevents):
    # Worker entry point. Generates nevents events for the pointing
    # (el_deg, az_deg) and writes them to filename.
    # Returns (filename, pointing_index, el_deg, az_deg,
    #          nevents, ntracks, nsteps, nrays, nbytes) or None if stopping.
    global instance_event_id

    if stop_requested:
        return None

    # Build a fresh IACT and resolve all pointing-dependent objects
    iact, vc_dir, bmax_poly, viewcone_poly, spectral_tf, sim_config = setup_pointing(el_deg, az_deg)
 
    # Write simulation configuration
    h5 = type(sim_config).NewHDFStreamWriter(filename, 'config', truncate=True)
    h5.write(sim_config)
    del h5

    h5 = calin.ix.simulation.simulated_event.SimulatedEvent.NewHDFStreamWriter(
        filename, 'events', truncate=False)

    num_events = 0
    num_tracks = 0
    num_steps  = 0
    num_rays   = 0

    for _ in range(nevents):
        if stop_requested:
            break
        sim_event = gen_event(iact, vc_dir, viewcone_poly, spectral_tf)
        sim_event.set_instance_id(instance_id)
        sim_event.set_event_id(instance_event_id)
        instance_event_id += 1
        num_events += 1
        h5.write(sim_event)
        num_tracks += iact.num_tracks()
        num_steps  += iact.num_steps()
        num_rays   += iact.num_rays()

    del h5

    # Write provenance anthology
    provenance_anthology = calin.provenance.anthology.get_current_anthology()
    h5 = type(provenance_anthology).NewHDFStreamWriter(filename, 'provenance', truncate=False)
    h5.write(provenance_anthology)
    del h5

    num_bytes = os.path.getsize(filename)
    return filename, pointing_index, el_deg, az_deg, num_events, num_tracks, num_steps, num_rays, num_bytes

# =============================================================================
# BANNER AND PROGRESS OUTPUT
# =============================================================================

def get_banner(sleep_time=0, el_deg=None, az_deg=None,
               bmax_poly=None, viewcone_poly=None, spectral_tf=None):
    # Use nominal zn=0 polynomials for the startup banner (before any pointing
    # is known); each per-pointing banner will carry the resolved values.
    if el_deg is None:
        zn_deg = 90.0
    else:
        zn_deg = 90.0 - el_deg

    if bmax_poly is None:
        bmax_poly = make_bmax_polynomial(zn_deg)
    if viewcone_poly is None:
        viewcone_poly = make_viewcone_polynomial(zn_deg)
    if spectral_tf is None:
        spectral_tf = make_spectral_transform(bmax_poly, viewcone_poly)

    banner = "Command line arguments :\n"
    args_dict = vars(saved_args)
    for arg in args_dict:
        banner += f'- {arg}: {args_dict[arg]}\n'

    banner += f'\nGrid parameters :\n'
    banner += f'- HealPix nside: {saved_args.grid_nside}\n'
    banner += f'- Number of cells: {num_grid_pointings(saved_args.grid_nside):,d}\n'
    banner += f'- Zenith cut: {saved_args.grid_znmax} deg\n'
    banner += f'- Cell area: {calin.math.healpix_array.cell_area(saved_args.grid_nside) * (180/numpy.pi)**2:.2f} deg^2\n'
    banner += f'- Cell dimension: {calin.math.healpix_array.cell_dimension(saved_args.grid_nside) * (180/numpy.pi):.2f} deg\n'
    banner += f'- Number of pointings: {total_pointings:,d}\n'

    if el_deg is not None:
        banner += f'\nPointing : El={el_deg:.2f} deg (Zn={zn_deg:.2f} deg), Az={az_deg:.2f} deg\n'

    banner += '\nPrimary viewcone :'
    if len(viewcone_poly) == 0 or (len(viewcone_poly) == 1 and viewcone_poly[0] == 0):
        banner += " DISABLED\n"
    elif len(viewcone_poly) == 1:
        banner += f' FIXED - {viewcone_poly[0]/numpy.pi*180.0:.1f} degrees\n'
    else:
        banner += " VARIABLE\n"
        banner += "  10GeV, 100GeV, 1TeV, 10TeV, 100TeV :"
        for i, x in enumerate([-2, -1, 0, 1, 2]):
            banner += f'{", " if i > 0 else " "}{numpy.polyval(viewcone_poly, x)/numpy.pi*180.0:.1f}'
        banner += " degrees\n"

    banner += '\nShower energy :'
    if saved_args.emin >= saved_args.emax:
        banner += f' MONOCHROMATIC at {format_energy(saved_args.emin)}'
    else:
        banner += ' SPECTRUM\n'
        x0 = spectral_tf(0)
        p0 = 0
        for p1 in 1 - 0.5 ** numpy.arange(1, 10):
            x1 = spectral_tf(p1)
            banner += f'  {format_energy(10**x0)} to {format_energy(10**x1)} : {(p1-p0)*100:.4f} %\n'
            x0 = x1
            p0 = p1
        p1 = 1.0
        x1 = spectral_tf(p1)
        banner += f'  {format_energy(10**x0)} to {format_energy(10**x1)} : {(p1-p0)*100:.4f} %'

    if sleep_time > 0:
        time.sleep(sleep_time)
    return banner


def print_line(filename, el_deg, az_deg,
               num_files_generated,
               num_events_total, num_rays_total, num_steps_total,
               num_tracks_total, num_bytes_total,
               begin_utc):
    runtime = (datetime.datetime.now(datetime.timezone.utc) - begin_utc).total_seconds()

    line = f'{filename} [El={el_deg:04.1f} Az={az_deg:05.1f}]:'
    
    # line += f' {num_files_generated:,d} / {total_files_to_generate:,d}'

    fraction = num_events_total / total_events_to_generate
    line += f' ; {num_events_total:,d} / {total_events_to_generate:,d} = {fraction*100.0:.2f}%'

    line += f' ; {format_duration(runtime)} / {format_duration(runtime/fraction)}'
    line += f' ({num_events_total/runtime:,.1f} Hz)'
    line += f' ; {num_rays_total:,d} rays'
    if num_tracks_total > 0:
        line += f' ; {num_rays_total/num_steps_total:.2f} {num_steps_total/num_tracks_total:.2f}'
    line += f' ; {num_bytes_total*1e-9:,.3f} / {num_bytes_total/fraction*1e-9:,.3f} GB'
    if stop_requested:
        line += ' (STOPPING)'
    print(line)

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Nectarcam multi-pointing shower simulation over a healpix pointing grid')

    # --- Output ---
    parser.add_argument('-o', '--output', type=str, default='events_{pointing}_{seq}.h5',
        help='Output filename template. Tokens: {id} random hex id, {seq} sequential file '
             'number within the pointing, {pointing} zero-padded pointing index '
             '(default: events_{pointing}_{seq}.h5)')
    parser.add_argument('--store_times_as_integer', action='store_true',
        help='Store times as 16-bit integers in output file, where possible, with 1-ps '
             'resolution (default: disabled)')

    # --- Pointing grid ---
    parser.add_argument('--grid_nside', type=int, default=8,
        help='Healpix nside parameter controlling the pointing grid spacing (default: 8)')
    parser.add_argument('--grid_znmax', type=float, default=65.0,
        help='Maximum zenith angle in degrees; pointings with zn > this value are skipped. '
             'The index-to-pointing mapping is unaffected (default: 65.0)')
    parser.add_argument('--nfiles_per_pointing', type=int, default=1,
        help='Number of HDF5 output files to generate per pointing (default: 1)')
    parser.add_argument('--block_size', type=int, default=1000,
        help='Number of shower events per output file (default: 1000)')

    # --- Site ---
    parser.add_argument('--site', type=str, default='ctan', choices=['ctan', 'ctas'],
        help='Site to simulate (default: ctan)')

    # --- Shower reuse ---
    parser.add_argument('--reuse', type=int, default=10,
        help='Number of times to reuse each shower (default: 10)')

    # --- Impact parameter (zn-keyed polynomial) ---
    parser.add_argument('-b', '--bmax_polynomial', type=zn_polynomial_list,
        action='append', default=None, metavar='[ZN:]P0[,P1,...]',
        help='Polynomial defining maximum shower impact parameter in meters as a function '
             'of energy in TeV. x=log10(E/1TeV), b(x)=P0+P1*x+P2*x^2+... '
             'May optionally be prefixed with a zenith angle in degrees and a colon '
             '(e.g. "20:1000.0,50.0") to apply that polynomial at that zenith angle. '
             'Multiple -b options may be given for different zenith angles; the nearest '
             'zenith angle entry is used for each pointing. (default: 1000.0)')

    # --- Viewcone (zn-keyed polynomial) ---
    parser.add_argument('--viewcone_polynomial', type=zn_polynomial_list,
        action='append', default=None, metavar='[ZN:]P0[,P1,...]',
        help='Polynomial defining viewcone half-angle in degrees as a function of energy '
             'in TeV. x=log10(E/1TeV), theta(x)=P0+P1*x+... '
             'May optionally be prefixed with a zenith angle in degrees and a colon '
             '(e.g. "20:3.0,0.5"). Multiple options may be given for different zenith '
             'angles. (default: 0.0)')
    parser.add_argument('--theta', type=float, default=0.0,
        help='Fixed offset between primary direction and telescope pointing direction in degrees')
    parser.add_argument('--phi', type=float, default=0.0,
        help='Fixed polar angle of primary direction around telescope pointing direction in degrees')
    parser.add_argument('--no_viewcone_cut', action='store_true',
        help='Do not apply viewcone cut to ray generation (default: cut enabled)')

    # --- Particle and spectrum ---
    parser.add_argument('-p', '--primary', type=str, default='gamma',
        choices=['gamma', 'muon', 'electron', 'proton', 'helium', 'iron'],
        help='Primary particle type (default: gamma)')
    parser.add_argument('--emin', type=float, default=0.01,
        help='Lower bound of energy spectrum in TeV (default: 0.01)')
    parser.add_argument('--emax', type=float, default=100.0,
        help='Upper bound of energy spectrum in TeV (default: 100.0)')
    parser.add_argument('--spectral_polynomial', type=comma_separated_floats, default=[-2.7],
        help='Spectral shape polynomial in log10(E/1TeV). dN/dE(x)=10^(P0*x+P1*x^2+...) '
             '(default: -2.7)')
    parser.add_argument('--spectral_constant', type=float, default=0.0,
        help='Fraction of simulated events using a flat spectrum in dN/dlog(E) (default: 0.0)')

    # --- PE spectrum and timing ---
    parser.add_argument('--tts', type=float, default=0.00,
        help='Transit time spread RMS in ns (default: 0.0)')
    parser.add_argument('--enable_pe_spectrum', action='store_true',
        help='Enable PE charge spectrum and store PE charges in output file (default: disabled)')

    # --- Physics options ---
    parser.add_argument('--no_bfield', action='store_true',
        help='Disable Earth magnetic field (default: enabled)')
    parser.add_argument('--no_refraction', action='store_true',
        help='Disable atmospheric refraction of rays (default: enabled)')
    parser.add_argument('--multiple_scattering', type=str, default='normal',
        choices=['minimal', 'simple', 'normal', 'better', 'insane'],
        help='Multiple scattering model (default: normal)')

    # --- Parallelism ---
    parser.add_argument('--nthread', type=int, default=0,
        help='Number of worker processes (default: 0 = number of CPUs available)')
    parser.add_argument('--avx', type=int, default=512, choices=[128, 256, 512],
        help='AVX vector size in bits (default: 512)')

    global args
    args = parser.parse_args()

    # Apply defaults for polynomial lists (argparse append starts from None)
    if args.bmax_polynomial is None:
        args.bmax_polynomial = [(0.0, [1000.0])]
    if args.viewcone_polynomial is None:
        args.viewcone_polynomial = [(0.0, [0.0])]

    # Sort both polynomial lists by zenith angle for clarity in banners/logs
    args.bmax_polynomial     = sorted(args.bmax_polynomial,     key=lambda t: t[0])
    args.viewcone_polynomial = sorted(args.viewcone_polynomial, key=lambda t: t[0])

    # Determine total pointings (ignoring znmax) for zero-padding width
    num_grid_cells = num_grid_pointings(args.grid_nside)
    pointing_index_width = len(str(num_grid_cells - 1)) if num_grid_cells > 0 else 1

    args.nfiles_per_pointing = max(1, args.nfiles_per_pointing)

    global total_pointings
    global total_files_to_generate
    global total_events_to_generate

    # Count how many pointings survive the znmax cut (for progress reporting)
    total_grid_pointings = num_grid_pointings(args.grid_nside)
    total_pointings = sum(
        1 for i in range(total_grid_pointings)
        if get_grid_pointing(args.grid_nside, i)[1] <= args.grid_znmax)

    total_files_to_generate = total_pointings * args.nfiles_per_pointing
    total_events_to_generate = total_files_to_generate * args.block_size

    max_workers = args.nthread or os.cpu_count() or 1

    global begin_utc
    global num_files_generated
    global num_events_total
    global num_tracks_total
    global num_steps_total
    global num_rays_total
    global num_bytes_total
    global stop_requested

    begin_utc        = datetime.datetime.now(datetime.timezone.utc)
    num_files_generated = 0
    num_events_total = 0
    num_tracks_total = 0
    num_steps_total  = 0
    num_rays_total   = 0
    num_bytes_total  = 0
    stop_requested   = False


    def process_results(result):
        if result is None:
            return
        global num_files_generated, num_events_total, num_tracks_total
        global num_steps_total, num_rays_total, num_bytes_total

        filename, p_idx, el_deg, az_deg, nevents, ntracks, nsteps, nrays, nbytes = result
        num_files_generated += 1
        num_events_total    += nevents
        num_tracks_total    += ntracks
        num_steps_total     += nsteps
        num_rays_total      += nrays
        num_bytes_total     += nbytes
        print_line(filename, el_deg, az_deg,
                   num_files_generated,
                   num_events_total, num_rays_total, num_steps_total, num_tracks_total,
                   num_bytes_total, begin_utc)

    def pointing_jobs():
        # Yield one job per (pointing x file), skipping pointings above grid_znmax.
        for i in range(total_grid_pointings):
            p_idx, zn_deg, az_deg = get_grid_pointing(args.grid_nside, i)
            if zn_deg > args.grid_znmax:
                continue
            el_deg = 90.0 - zn_deg
            for seq in range(args.nfiles_per_pointing):
                fname = pointing_filename(
                    args.output, p_idx, seq, pointing_index_width)
                yield fname, p_idx, el_deg, az_deg, args.block_size

    if max_workers == 1:
        init(args)
        print(get_banner())
        print()
        signal.signal(signal.SIGINT, handle_sigint)

        for fname, p_idx, el_deg, az_deg, nevents in pointing_jobs():
            if stop_requested:
                break
            result = gen_batch(fname, p_idx, el_deg, az_deg, nevents)
            process_results(result)

    else:
        in_flight_cap = 4 * max_workers
        jobs = pointing_jobs()

        print(f'Launching {max_workers} workers in process pool')
        with concurrent.futures.ProcessPoolExecutor(
                initializer=init, initargs=(args,),
                max_workers=max_workers) as executor:

            futures = {}

            # Wait for all workers to initialise before printing the banner
            warmup = set()
            for _ in range(max_workers):
                warmup.add(executor.submit(get_banner, 5))
            banner = None
            while warmup:
                done, _ = concurrent.futures.wait(
                    warmup, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    warmup.discard(fut)
                    banner = fut.result()
            if banner is not None:
                print(banner)
                print()

            signal.signal(signal.SIGINT, handle_sigint)

            # Fill initial queue up to cap
            for job in jobs:
                if len(futures) >= in_flight_cap or stop_requested:
                    break
                fname, p_idx, el_deg, az_deg, nevents = job
                fut = executor.submit(gen_batch, fname, p_idx, el_deg, az_deg, nevents)
                futures[fut] = p_idx

            # Drain completions and refill
            while futures:
                done, _ = concurrent.futures.wait(
                    futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    futures.pop(fut)
                    process_results(fut.result())
                    if not stop_requested:
                        try:
                            fname, p_idx, el_deg, az_deg, nevents = next(jobs)
                            new_fut = executor.submit(
                                gen_batch, fname, p_idx, el_deg, az_deg, nevents)
                            futures[new_fut] = p_idx
                        except StopIteration:
                            pass
