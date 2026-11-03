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
import sys
import re
import glob
import argparse
import signal
import collections
import itertools
import traceback
import multiprocessing
import multiprocessing.connection
import datetime
import time
import numpy
import scipy.interpolate
import calin.math.geometry
import calin.math.healpix_array
import calin.simulation.iact_factory
import calin.ix.simulation.simulated_event
import calin.provenance.anthology

# =============================================================================
# PROCESS-WIDE STATE
#
# stop_requested : set by the SIGINT handler in the main process (and polled by
#                  the in-process generation loop when running with one thread)
# stop_flag      : multiprocessing.RawValue shared between the main process and
#                  the workers; main sets it to 1 on the first Ctrl-C and workers
#                  poll it between events (a plain shared byte, so it is safe to
#                  set from a signal handler)
# heartbeat_fn   : set in worker processes, called from gen_batch() so the
#                  supervisor can detect hung workers
# =============================================================================

stop_requested = False
sigint_count = 0
stop_flag = None
heartbeat_fn = None

HEARTBEAT_INTERVAL = 5.0            # seconds between heartbeats from a busy worker
MAX_CONSECUTIVE_INIT_FAILURES = 10  # abort if workers repeatedly die during init

# Command line arguments that may legitimately differ between the original run
# and a --resume run. All other arguments must match those stored in the
# configuration of each existing output file.
RESUME_IGNORED_ARGS = frozenset([
    'output', 'block_size', 'nfiles_per_pointing', 'grid_znmax', 'pointing',
    'nthread', 'avx', 'resume', 'prune_dead_files', 'ignore_config_mismatch',
    'max_attempts', 'job_timeout', 'recycle_after', 'dry_run'])

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

def handle_sigint(signum, frame):
    # First Ctrl-C : stop gracefully, workers finish their current event and
    #                close their files (events written so far are kept).
    # Second Ctrl-C: raise KeyboardInterrupt in the main process, which kills the
    #                workers immediately (files being written are left incomplete).
    global stop_requested, sigint_count
    sigint_count += 1
    if sigint_count == 1:
        print("\nStop requested - terminating after current files complete "
              "(Ctrl-C again to terminate immediately)", flush=True)
    elif sigint_count == 2:
        raise KeyboardInterrupt
    stop_requested = True
    if stop_flag is not None:
        stop_flag.value = 1

def stopping():
    return stop_requested or (stop_flag is not None and stop_flag.value != 0)

def beat():
    if heartbeat_fn is not None:
        heartbeat_fn()

def log(msg):
    print(f'[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}', flush=True)

def describe_exit(code):
    if code is None:
        return 'unknown exit status'
    if code < 0:
        try:
            return f'killed by signal {-code} ({signal.Signals(-code).name})'
        except ValueError:
            return f'killed by signal {-code}'
    return f'exit status {code}'

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

    # Load site environment and detector response models once per worker.
    global site_env
    global store_pe_weights
    site_env = calin.simulation.iact_factory.load_site_environment(
        args.site, enable_pe_spectrum=args.enable_pe_spectrum,
        no_bfield=args.no_bfield, quiet=True, avx=args.avx)
    store_pe_weights = site_env.pe_gen is not None

    global store_times_as_integer
    store_times_as_integer = args.store_times_as_integer

    # bmax and viewcone polynomials are resolved per-pointing in setup_pointing();
    # nothing to compute here.

    # The site environment also resolves the magnetic field at the array origin.
    # Keep site_env alive for the C++ generator and propagators that use its models.
    global generator
    global saved_geant4_cfg
    generator, saved_geant4_cfg = calin.simulation.iact_factory.create_geant4_generator(
        site_env.atm, bfield=site_env.bfield,
        multiple_scattering=args.multiple_scattering,
        primary=args.primary)

    # Particle type
    global particle_type
    particle_type = calin.simulation.iact_factory.get_tracker_particle_type(args.primary)

    # The spectral transform depends on the per-pointing bmax and viewcone
    # polynomials, so it is computed inside setup_pointing() for each pointing.

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
    global iact_cfg
    global all_pe_processor
    global all_prop
    array_params, nscope, nchan, detector_type_name = \
        calin.simulation.iact_factory.setup_telescope_array(
            saved_args.site, el_deg=el_deg)

    # Build and return (iact, vc_dir, bmax_poly, viewcone_poly,
    # spectral_tf, sim_config) for the pointing (el_deg, az_deg).
    zn_deg = 90.0 - el_deg

    # Resolve polynomials for this zenith angle
    bmax_poly     = make_bmax_polynomial(zn_deg)
    viewcone_poly = make_viewcone_polynomial(zn_deg)
    spectral_tf   = make_spectral_transform(bmax_poly, viewcone_poly)

    # Build a fresh IACT instance for this pointing and attach its propagators.
    iact, iact_cfg, _ = calin.simulation.iact_factory.create_iact_array(
        site_env.atm, site_env.atm_abs, avx=saved_args.avx,
        no_refraction=saved_args.no_refraction)
    all_pe_processor, all_prop = calin.simulation.iact_factory.attach_iact_propagators(
        iact, saved_args.site, array_params,
        bmax_polynomial=numpy.flipud(bmax_poly) * 0.01,
        reuse=saved_args.reuse, nscope=nscope, nchan=nchan,
        det_eff=site_env.det_eff, cone_eff=site_env.cone_eff,
        pe_gen=site_env.pe_gen,
        tts=saved_args.tts,
        lens_spline=site_env.lens_refractive_index_spline,
        detector_type_name=detector_type_name)

    pt_dir, vc_dir = calin.simulation.iact_factory.set_iact_pointing(
        iact, el_deg, az_deg, saved_args.theta, saved_args.phi,
        apply_viewcone_cut=not saved_args.no_viewcone_cut)

    # Simulation configuration protobuf for this pointing
    sim_config = calin.ix.simulation.simulated_event.SimulationConfiguration()
    sim_config.set_particle_type(
        calin.simulation.iact_factory.get_simulated_event_particle_type(
            saved_args.primary))
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

    atm_abs_zmin = numpy.min(site_env.atm_abs.levels_cm())
    atm_abs_zmax = numpy.max(site_env.atm_abs.levels_cm())
    for level in site_env.atm.get_levels():
        proto_level = sim_config.add_atmospheric_level()
        proto_level.set_altitude(level.z)
        proto_level.set_thickness(level.t)
        proto_level.set_density(level.rho)
        proto_level.set_n_minus_one(level.nmo)
        if site_env.bfield is not None:
            b = site_env.bfield.field_nT(level.z)
            proto_level.mutable_bfield().set_x(b[0])
            proto_level.mutable_bfield().set_y(b[1])
            proto_level.mutable_bfield().set_z(b[2])
        if atm_abs_zmin < level.z < atm_abs_zmax:
            proto_level.set_optical_depth_1d5ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 1.5))
            proto_level.set_optical_depth_2d0ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.0))
            proto_level.set_optical_depth_2d5ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.5))
            proto_level.set_optical_depth_3d0ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.0))
            proto_level.set_optical_depth_3d5ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.5))
            proto_level.set_optical_depth_4d0ev(site_env.atm_abs.optical_depth_for_altitude_and_energy(level.z, 4.0))

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
    detector_type_config.set_detector_efficiency_banner(site_env.det_eff.banner())
    detector_type_config.set_detector_efficiency_energy(numpy.asarray(site_env.det_eff.all_xi()))
    detector_type_config.set_detector_efficiency_efficiency(numpy.asarray(site_env.det_eff.all_yi()))
    if site_env.cone_eff:
        detector_type_config.set_angular_response_banner(site_env.cone_eff.banner())
        detector_type_config.set_angular_response_costheta(numpy.asarray(site_env.cone_eff.all_xi()))
        detector_type_config.set_angular_response_efficiency(numpy.asarray(site_env.cone_eff.all_yi()))
    if site_env.pe_gen:
        detector_type_config.set_pe_spectrum_banner(site_env.pe_gen.banner())
        detector_type_config.set_pe_spectrum_banner_q(site_env.pe_gen.raw_q())
        detector_type_config.set_pe_spectrum_banner_dp_dq(site_env.pe_gen.raw_dp_dq())
    detector_type_config.set_pe_time_spread(saved_args.tts)
    if saved_args.site != 'dark100':
        detector_type_config.mutable_dc_array_parameters().CopyFrom(array_params)

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

    ct0 = 1.0 / u[2] * (site_env.atm.top_of_atmosphere() - site_env.atm.zobs(0))
    x0  = numpy.asarray([0, 0, site_env.atm.zobs(0)]) + ct0 * u

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
    #          nevents, ntracks, nsteps, nrays, nbytes) or None if stopping
    # before any event was generated (no file is left behind in that case).
    global instance_event_id

    if stopping():
        return None
    beat()

    # Build a fresh IACT and resolve all pointing-dependent objects
    iact, vc_dir, bmax_poly, viewcone_poly, spectral_tf, sim_config = setup_pointing(el_deg, az_deg)
    beat()
    if stopping():
        return None

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
        if stopping():
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
        beat()

    del h5

    if num_events == 0:
        os.remove(filename)
        return None

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
        if saved_args.pointing is not None:
            _, zn_deg, az_deg = get_grid_pointing(saved_args.grid_nside, saved_args.pointing)
            el_deg = 90.0 - zn_deg
        else:
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
    if saved_args.pointing is not None:
        banner += f'- Single pointing mode: Index {saved_args.pointing}\n'
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

    fraction = max(num_events_total / total_events_to_generate, 1e-12)
    line += f' ; {num_events_total:,d} / {total_events_to_generate:,d} = {fraction*100.0:.2f}%'

    line += f' ; {format_duration(runtime)} / {format_duration(runtime/fraction)}'
    line += f' ({num_events_total/runtime:,.1f} Hz)'
    line += f' ; {num_rays_total:,d} rays'
    if num_tracks_total > 0 and num_steps_total > 0:
        line += f' ; {num_rays_total/num_steps_total:.2f} {num_steps_total/num_tracks_total:.2f}'
    line += f' ; {num_bytes_total*1e-9:,.3f} / {num_bytes_total/fraction*1e-9:,.3f} GB'
    if stop_requested:
        line += ' (STOPPING)'
    print(line)

# =============================================================================
# WORKER PROCESS AND SUPERVISOR
#
# Workers are plain multiprocessing.Process objects, each with its own pipe to
# the main process. Unlike concurrent.futures.ProcessPoolExecutor, the death of
# one worker (SIGSEGV/SIGABRT in the C++ code, OOM kill, ...) does not affect the
# others: the Supervisor notices the dead process, requeues the job it was
# running (same output filename, which is simply overwritten) and starts a
# replacement worker. A job that fails max_attempts times is skipped and
# reported at the end.
#
# Messages worker -> main : ('ready',) ('hb',) ('done', result) ('error', tb)
#                           ('init_error', tb)
# Messages main -> worker : job tuple (fname, pointing, el_deg, az_deg, nevents)
#                           or None (exit)
# =============================================================================

def worker_main(conn, args, shared_stop_flag, parent_pid, n_total_pointings):
    global stop_flag, heartbeat_fn, total_pointings, stop_requested
    signal.signal(signal.SIGINT, signal.SIG_IGN)    # the main process owns Ctrl-C
    stop_flag = shared_stop_flag
    total_pointings = n_total_pointings             # used by get_banner()
    stop_requested = False
    last_beat = [0.0]

    def beat_():
        if os.getppid() != parent_pid:              # main process is gone
            os._exit(1)
        now = time.monotonic()
        if now - last_beat[0] >= HEARTBEAT_INTERVAL:
            last_beat[0] = now
            conn.send(('hb',))
    heartbeat_fn = beat_

    try:
        init(args)
    except BaseException:
        conn.send(('init_error', traceback.format_exc()))
        return
    conn.send(('ready',))

    while True:
        if not conn.poll(5.0):
            if os.getppid() != parent_pid:
                return
            continue
        try:
            job = conn.recv()
        except EOFError:
            return
        if job is None:
            return
        last_beat[0] = 0.0
        try:
            result = gen_batch(*job)
        except Exception:
            # State of the C++ objects is unknown after an exception, so report
            # and exit; the supervisor starts a replacement.
            conn.send(('error', traceback.format_exc()))
            return
        conn.send(('done', result))


class Supervisor:
    def __init__(self, nworkers, jobs, njobs, args, on_result, on_started,
                 max_attempts=3, job_timeout=0.0, recycle_after=0):
        self.ctx = multiprocessing.get_context('fork')
        self.nworkers = nworkers
        self.jobs = iter(jobs)
        self.outstanding = njobs            # jobs not yet done or finally failed
        self.args = args
        self.on_result = on_result
        self.on_started = on_started
        self.max_attempts = max_attempts
        self.job_timeout = job_timeout      # seconds without a message; 0 = off
        self.recycle_after = recycle_after  # jobs per worker; 0 = never
        self.workers = {}
        self.next_wid = 0
        self.lookahead = collections.deque()
        self.retry = collections.deque()
        self.attempts = collections.Counter()
        self.failed = []                    # (job, reason)
        self.initial_workers = min(nworkers, njobs)
        self.nready = 0
        self.started = False
        self.init_failures = 0

    def run(self):
        try:
            self._run()
        finally:
            self._cleanup()
        return self.failed

    # --- main loop -----------------------------------------------------------

    def _run(self):
        while True:
            self._rebalance()
            if not self.workers:
                break
            if not self.started and self.nready >= self.initial_workers:
                self.started = True
                self.on_started()
            self._dispatch()
            waitables = [w['conn'] for w in self.workers.values()]
            waitables += [w['proc'].sentinel for w in self.workers.values()]
            multiprocessing.connection.wait(waitables, timeout=1.0)
            for wid, w in list(self.workers.items()):
                self._drain(w)
                if w['proc'].exitcode is not None:
                    self._drain(w)
                    self._reap(wid)
            self._check_timeouts()

    def _rebalance(self):
        # Keep min(nworkers, outstanding jobs) workers alive (not counting those
        # that have been told to exit).
        if stop_requested:
            return
        want = min(self.nworkers, self.outstanding)
        have = sum(1 for w in self.workers.values() if w['state'] != 'retiring')
        for _ in range(want - have):
            self._spawn()

    def _spawn(self):
        wid = self.next_wid
        self.next_wid += 1
        parent_conn, child_conn = self.ctx.Pipe()
        proc = self.ctx.Process(
            target=worker_main, name=f'worker-{wid}', daemon=True,
            args=(child_conn, self.args, stop_flag, os.getpid(), total_pointings))
        proc.start()
        child_conn.close()
        self.workers[wid] = dict(proc=proc, conn=parent_conn, job=None,
                                 state='starting', last_msg=time.monotonic(),
                                 njobs=0, kill_reason=None)

    def _has_work(self):
        if self.retry or self.lookahead:
            return True
        try:
            self.lookahead.append(next(self.jobs))
        except StopIteration:
            return False
        return True

    def _take_job(self):
        return self.retry.popleft() if self.retry else self.lookahead.popleft()

    def _dispatch(self):
        for w in self.workers.values():
            if w['state'] == 'starting' and stop_requested:
                w['state'] = 'retiring'
                w['proc'].terminate()
            elif w['state'] == 'idle':
                if stop_requested or self.outstanding <= 0:
                    self._retire(w)
                elif self.started and self._has_work():
                    self._assign(w, self._take_job())

    def _retire(self, w):
        w['state'] = 'retiring'
        try:
            w['conn'].send(None)
        except (OSError, ValueError):
            pass

    def _assign(self, w, job):
        w['job'] = job
        w['state'] = 'busy'
        w['last_msg'] = time.monotonic()
        try:
            w['conn'].send(job)
        except (OSError, ValueError):
            # Worker already dead; the job was not started, so no attempt is used.
            w['job'] = None
            w['state'] = 'retiring'
            self.retry.appendleft(job)

    # --- messages from workers ---------------------------------------------------

    def _drain(self, w):
        while True:
            try:
                if not w['conn'].poll():
                    return
                msg = w['conn'].recv()
            except (EOFError, OSError):
                return
            self._handle(w, msg)

    def _handle(self, w, msg):
        w['last_msg'] = time.monotonic()
        kind = msg[0]
        if kind == 'ready':
            w['state'] = 'idle'
            self.nready += 1
            self.init_failures = 0
        elif kind == 'done':
            w['job'] = None
            w['state'] = 'idle'
            w['njobs'] += 1
            self.outstanding -= 1
            self.on_result(msg[1])
            if self.recycle_after and w['njobs'] >= self.recycle_after:
                self._retire(w)
        elif kind == 'error':
            job = w['job']
            w['job'] = None
            w['state'] = 'retiring'         # worker exits after reporting
            log(f'Exception in worker pid {w["proc"].pid}:\n{msg[1]}')
            self._job_failed(job, 'exception in worker')
        elif kind == 'init_error':
            log(f'Worker pid {w["proc"].pid} failed during initialisation:\n{msg[1]}')

    def _job_failed(self, job, reason):
        fname, p_idx, el_deg, az_deg, nevents = job
        self.attempts[fname] += 1
        n = self.attempts[fname]
        log(f'FAILED {fname} (pointing {p_idx}, El={el_deg:.2f} Az={az_deg:.2f}, '
            f'{nevents:,d} events): {reason}; attempt {n}/{self.max_attempts}')
        if stop_requested:
            self.failed.append((job, reason + ' (not retried: stop requested)'))
        elif n < self.max_attempts:
            self.retry.append(job)
        else:
            self.failed.append((job, reason))
            self.outstanding -= 1

    # --- worker exit -----------------------------------------------------------------

    def _reap(self, wid):
        w = self.workers.pop(wid)
        w['conn'].close()
        w['proc'].join()
        code = w['proc'].exitcode
        pid = w['proc'].pid
        if w['job'] is not None:
            reason = w['kill_reason'] or describe_exit(code)
            self._job_failed(w['job'], f'worker pid {pid} died: {reason}')
        elif w['state'] == 'starting':
            self.init_failures += 1
            log(f'Worker pid {pid} died during initialisation: {describe_exit(code)}')
            if self.init_failures >= MAX_CONSECUTIVE_INIT_FAILURES:
                raise RuntimeError(
                    f'{self.init_failures} consecutive worker initialisation failures')
        elif w['state'] == 'idle':
            log(f'Idle worker pid {pid} died: {describe_exit(code)}')

    def _check_timeouts(self):
        if not self.job_timeout:
            return
        now = time.monotonic()
        for w in self.workers.values():
            if (w['job'] is not None and w['kill_reason'] is None
                    and now - w['last_msg'] > self.job_timeout):
                w['kill_reason'] = (f'no heartbeat for {self.job_timeout/60:.1f} min '
                                    f'(killed by supervisor)')
                w['proc'].kill()

    def _cleanup(self):
        for w in self.workers.values():
            if w['proc'].is_alive():
                w['proc'].terminate()
        for w in self.workers.values():
            w['proc'].join(5.0)
            if w['proc'].is_alive():
                w['proc'].kill()
                w['proc'].join()
            w['conn'].close()
        self.workers.clear()

# =============================================================================
# RESUME
#
# With --resume, existing output files matching the --output template are
# inspected with the HDF reader. Each file that opens contributes the number of
# rows in its 'events' dataset to its pointing (taken from the {pointing} token
# in the filename); files that cannot be read (events, provenance or config
# missing/unreadable, e.g. because the worker died while writing) are "dead" and
# count as zero events. Pointings that already have enough events are skipped,
# others are topped up with the missing events.
# =============================================================================

def template_to_patterns(template):
    # Returns (glob pattern, compiled regex whose group 1 is the pointing index)
    glob_pat = ''
    rx = ''
    for part in re.split(r'(\{id\}|\{seq\}|\{pointing\})', template):
        if part == '{id}':
            glob_pat += '*'
            rx += '[0-9A-Fa-f]{32}'
        elif part == '{seq}':
            glob_pat += '*'
            rx += r'\d+'
        elif part == '{pointing}':
            glob_pat += '*'
            rx += r'(\d+)'
        else:
            glob_pat += glob.escape(part)
            rx += re.escape(part)
    return glob_pat, re.compile(rx)

def inspect_file(fname, anthology_type):
    # Returns (number of events, stored command line args as dict of str) or
    # raises if any of the events, provenance or config groups cannot be read.
    event_type = calin.ix.simulation.simulated_event.SimulatedEvent
    h5 = event_type.NewHDFStreamReader(fname, 'events')
    nevents = h5.nrow()
    del h5

    h5 = anthology_type.NewHDFStreamReader(fname, 'provenance')
    if h5.nrow() < 1:
        raise RuntimeError('provenance group is empty')
    del h5

    config = calin.ix.simulation.simulated_event.SimulationConfiguration()
    h5 = type(config).NewHDFStreamReader(fname, 'config')
    if h5.nrow() < 1:
        raise RuntimeError('config group is empty')
    h5.read(0, config)
    del h5
    stored_args = {k: config.command_line_args(k) for k in config.command_line_args_keys()}
    return nevents, stored_args

def scan_existing(args):
    # Returns (events per pointing, dead files, config mismatches, number of good files)
    #   dead files       : list of (filename, pointing, error)
    #   config mismatches: dict (arg, stored value, current value) -> [nfiles, example file]
    glob_pat, rx = template_to_patterns(args.output)
    found = []
    for fname in sorted(glob.glob(glob_pat)):
        m = rx.fullmatch(fname)
        if m is None:
            continue
        p_idx = int(m.group(1))
        if args.pointing is not None and p_idx != args.pointing:
            continue
        found.append((fname, p_idx))

    print(f'Resume: scanning {len(found):,d} existing files')
    anthology_type = type(calin.provenance.anthology.get_current_anthology())
    current_args = {k: str(v) for k, v in vars(args).items() if k not in RESUME_IGNORED_ARGS}

    counts = collections.Counter()
    dead = []
    mismatches = {}
    ngood = 0
    t_last = time.monotonic()
    for i, (fname, p_idx) in enumerate(found):
        try:
            nevents, stored_args = inspect_file(fname, anthology_type)
        except Exception as e:
            dead.append((fname, p_idx, f'{type(e).__name__}: {e}'))
            continue
        ngood += 1
        counts[p_idx] += nevents
        for k, v in current_args.items():
            if stored_args.get(k) != v:
                entry = mismatches.setdefault((k, stored_args.get(k), v), [0, fname])
                entry[0] += 1
        if time.monotonic() - t_last > 10.0:
            t_last = time.monotonic()
            print(f'  scanned {i+1:,d} / {len(found):,d} files', flush=True)
    return counts, dead, mismatches, ngood

def prune_dead_files(dead):
    # Delete the dead files, but only after the user types "yes".
    print(f'{len(dead):,d} dead files will be removed, e.g.:')
    for fname, p_idx, err in dead[:5]:
        print(f'  {fname} ({err})')
    try:
        answer = input(f'Type "yes" to delete {len(dead):,d} files: ')
    except EOFError:
        answer = ''
    if answer.strip() != 'yes':
        print('Not deleting dead files; new files will use different names')
        return
    nremoved = 0
    for fname, p_idx, err in dead:
        try:
            os.remove(fname)
            nremoved += 1
        except OSError as e:
            print(f'  could not remove {fname}: {e}')
    print(f'Removed {nremoved:,d} dead files')

# =============================================================================
# MAIN
# =============================================================================

def build_parser():
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
    parser.add_argument('--pointing', type=int, default=None,
        help='Single pointing index to simulate. When specified, multi-process threading '
             'is automatically disabled and only this pointing is processed (default: all valid grid pointings)')
    parser.add_argument('--nfiles_per_pointing', type=int, default=1,
        help='Number of HDF5 output files to generate per pointing (default: 1)')
    parser.add_argument('--block_size', type=int, default=1000,
        help='Number of shower events per output file (default: 1000)')

    # --- Site ---
    parser.add_argument('--site', type=str, default='ctan', choices=['ctan', 'ctas', 'dark100'],
        help='Site / instrument to simulate: ctan, ctas, or dark100 (default: ctan)')

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

    # --- Fault tolerance and resume ---
    parser.add_argument('--resume', action='store_true',
        help='Inspect existing output files matching the --output template (which must '
             'contain {pointing}) and only generate the events still missing for each '
             'pointing. Files that cannot be read count as zero events. Do not resume '
             'while another instance is writing to the same files (default: disabled)')
    parser.add_argument('--prune_dead_files', action='store_true',
        help='With --resume, offer to delete unreadable files found by the scan; '
             'nothing is deleted unless "yes" is typed at the prompt (default: disabled)')
    parser.add_argument('--ignore_config_mismatch', action='store_true',
        help='With --resume, only warn (instead of aborting) if the command line '
             'arguments stored in existing files differ from the current ones, '
             'excluding output, block_size, nfiles_per_pointing, grid_znmax, pointing, '
             'nthread, avx and the resume/fault-tolerance options (default: abort)')
    parser.add_argument('--max_attempts', type=int, default=3,
        help='Maximum number of attempts for a file whose worker crashes or raises; '
             'after this the file is skipped and reported at the end (default: 3)')
    parser.add_argument('--job_timeout', type=float, default=0.0,
        help='Kill a worker, and count a failed attempt, if it sends no heartbeat for '
             'this many minutes. Heartbeats are sent between events, so this must exceed '
             'the longest single-event time (default: 0 = disabled)')
    parser.add_argument('--recycle_after', type=int, default=0,
        help='Replace each worker process after it has completed this many files, '
             'as protection against memory leaks (default: 0 = never)')

    # --- Dry run ---
    parser.add_argument('--dry_run', action='store_true',
        help='Do everything except launch the simulation: parse the arguments, scan existing '
             'files and prune dead files if requested (--resume, --prune_dead_files), then '
             'print the banner and the work to be done and exit (default: disabled)')

    return parser

def main(argv=None):
    global saved_args, stop_flag, stop_requested, sigint_count
    global total_pointings, total_files_to_generate, total_events_to_generate
    global begin_utc, num_files_generated, num_events_total, num_tracks_total
    global num_steps_total, num_rays_total, num_bytes_total

    parser = build_parser()
    args = parser.parse_args(argv)

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

    if args.prune_dead_files and not args.resume:
        parser.error('--prune_dead_files requires --resume')
    if args.resume and '{pointing}' not in args.output:
        parser.error('--resume requires the {pointing} token in the --output template')
    if args.max_attempts < 1:
        parser.error('--max_attempts must be at least 1')

    # get_banner() and friends read saved_args; set it here so that the main
    # process can print the banner without initialising the simulation.
    saved_args = args

    total_grid_pointings = num_grid_pointings(args.grid_nside)

    if args.pointing is not None:
        if args.pointing < 0 or args.pointing >= total_grid_pointings:
            parser.error(f'--pointing {args.pointing} is out of range [0, {total_grid_pointings - 1}] for grid_nside={args.grid_nside}')
        p_idx, p_zn, p_az = get_grid_pointing(args.grid_nside, args.pointing)
        if p_zn > args.grid_znmax:
            parser.error(f'--pointing {args.pointing} has zenith angle {p_zn:.2f} deg which exceeds --grid_znmax {args.grid_znmax:.2f} deg')
        valid_pointings = [(p_idx, p_zn, p_az)]
        max_workers = 1
    else:
        # Pointings that survive the znmax cut
        valid_pointings = []
        for i in range(total_grid_pointings):
            p_idx, p_zn, p_az = get_grid_pointing(args.grid_nside, i)
            if p_zn <= args.grid_znmax:
                valid_pointings.append((p_idx, p_zn, p_az))
        max_workers = args.nthread or os.cpu_count() or 1
    total_pointings = len(valid_pointings)

    # Events needed per pointing, minus those already in readable existing files
    events_per_pointing = args.nfiles_per_pointing * args.block_size
    existing = collections.Counter()
    if args.resume:
        existing, dead, mismatches, ngood = scan_existing(args)
        print(f'Resume: {ngood:,d} readable files with {sum(existing.values()):,d} events, '
              f'{len(dead):,d} dead files')
        if mismatches:
            print('Resume: command line arguments differ from those stored in existing files:')
            for (k, stored, current), (nfiles, example) in list(mismatches.items())[:30]:
                print(f'  {k}: stored={stored!r} current={current!r} '
                      f'({nfiles:,d} files, e.g. {example})')
            if not args.ignore_config_mismatch:
                print('Aborting; use --ignore_config_mismatch to continue anyway')
                return 2
            print('Continuing because of --ignore_config_mismatch')
        if dead:
            if args.prune_dead_files:
                prune_dead_files(dead)
            else:
                print('Dead files are ignored and left in place (see --prune_dead_files)')

    plan = [(p_idx, zn, az, events_per_pointing - existing[p_idx])
            for p_idx, zn, az in valid_pointings if existing[p_idx] < events_per_pointing]
    num_jobs = sum(-(-remaining // args.block_size) for _, _, _, remaining in plan)

    total_files_to_generate = num_jobs
    total_events_to_generate = sum(remaining for _, _, _, remaining in plan)

    if args.resume:
        ncomplete = sum(1 for p_idx, _, _ in valid_pointings if existing[p_idx] >= events_per_pointing)
        nuntouched = sum(1 for p_idx, _, _ in valid_pointings if existing[p_idx] == 0)
        npartial = total_pointings - ncomplete - nuntouched
        print(f'Resume: {ncomplete:,d} pointings complete, {npartial:,d} partial, '
              f'{nuntouched:,d} untouched; {total_events_to_generate:,d} events in '
              f'{num_jobs:,d} files to generate')
    if num_jobs == 0:
        print('Nothing to generate')
        return 0

    begin_utc        = datetime.datetime.now(datetime.timezone.utc)
    num_files_generated = 0
    num_events_total = 0
    num_tracks_total = 0
    num_steps_total  = 0
    num_rays_total   = 0
    num_bytes_total  = 0
    stop_requested   = False
    sigint_count     = 0

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
        # Yield one job per file still to be written: (fname, pointing, el, az, nevents)
        for p_idx, zn_deg, az_deg, remaining in plan:
            el_deg = 90.0 - zn_deg
            seq = 0
            while remaining > 0:
                nevents = min(args.block_size, remaining)
                fname = pointing_filename(args.output, p_idx, seq, pointing_index_width)
                if args.resume:
                    # Never overwrite an existing file (dead files that were not pruned)
                    while os.path.exists(fname):
                        seq += 1
                        fname = pointing_filename(args.output, p_idx, seq, pointing_index_width)
                seq += 1
                remaining -= nevents
                yield fname, p_idx, el_deg, az_deg, nevents

    if args.dry_run:
        print(get_banner())
        print()
        mode = 'in-process' if max_workers == 1 else f'{min(max_workers, num_jobs):,d} worker processes'
        print(f'Dry run: would generate {total_events_to_generate:,d} events in {num_jobs:,d} files '
              f'for {len(plan):,d} pointings, {mode}')
        print('Example files (names containing {id} are random):')
        for fname, p_idx, el_deg, az_deg, nevents in itertools.islice(pointing_jobs(), 5):
            print(f'  {fname} [pointing {p_idx}, El={el_deg:.2f} Az={az_deg:.2f}, {nevents:,d} events]')
        noverwrite = sum(1 for job in pointing_jobs() if os.path.exists(job[0]))
        if noverwrite:
            print(f'WARNING: {noverwrite:,d} of these files already exist and would be overwritten')
        print('Nothing launched')
        return 0

    failed = []
    try:
        if max_workers == 1:
            # In-process: no protection against crashes, but easy to debug
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
            nworkers = min(max_workers, num_jobs)
            print(f'Launching {nworkers} workers')
            stop_flag = multiprocessing.get_context('fork').RawValue('b', 0)
            signal.signal(signal.SIGINT, handle_sigint)

            def on_started():
                print(get_banner())
                print()

            supervisor = Supervisor(
                nworkers, pointing_jobs(), num_jobs, args, process_results, on_started,
                max_attempts=args.max_attempts, job_timeout=args.job_timeout*60.0,
                recycle_after=args.recycle_after)
            failed = supervisor.run()
    except KeyboardInterrupt:
        print('\nTerminated immediately; files being written were left incomplete '
              '(use --resume --prune_dead_files to clean up and continue)', flush=True)
        return 130

    rc = 0
    if failed:
        rc = 1
        print(f'\n{len(failed):,d} files FAILED and were not generated:')
        for (fname, p_idx, el_deg, az_deg, nevents), reason in failed:
            print(f'  {fname} (pointing {p_idx}, El={el_deg:.2f} Az={az_deg:.2f}, '
                  f'{nevents:,d} events): {reason}')
        print('Any partial files are unreadable and will be regenerated by --resume')
    if stop_requested:
        print(f'\nStopped early with {num_events_total:,d} / {total_events_to_generate:,d} '
              f'events generated; rerun with --resume to finish')
    return rc

if __name__ == '__main__':
    sys.exit(main())
