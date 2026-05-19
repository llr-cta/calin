# Example usage: python3.12 shower_trigger_threshold.py -o test.pickle -n 6000000 -b 6500 -e 0.0175 --nsb 0.0 -p proton --viewcone 13.0 --write_batch=100000 --omit_untriggered --reuse 200 --enable_viewcone_cut

# calin/scripts/simulate_events.py -- Stephen Fegan - 2026-01-15
#
# Simulate events from given spectrum with GEANT4 accountig for variable scattering 
# radius and viewcone. Write events to an HDF5 file
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
import concurrent
import datetime
import time
import numpy
import scipy.interpolate
import calin.math.geometry
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

def comma_separated_floats(value):
    try:
        return [float(x) for x in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a comma-separated list of floats (e.g. 1.0,2.3,4.2)")

def handle_sigint_quietly(signum, frame):
    global stop_requested
    stop_requested = True

def init(args):
    global saved_args
    saved_args = args

    numpy.random.seed()

    global instance_id
    global instance_event_id
    instance_id = numpy.random.bit_generator.randbits(64)
    instance_event_id = 0

    global block_size
    block_size = args.block_size

    # Select simulation classes based on AVX size requested
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
    global mst
    global nscope
    if args.site == 'ctan':
        zobs = calin.simulation.vs_cta.ctan_observation_level()
        atm = calin.simulation.vs_cta.ctan_atmosphere(quiet=True)
        atm_abs = calin.simulation.vs_cta.ctan_atmospheric_absorption(quiet=True)
        mst = calin.simulation.vs_cta.mstn1_config()
    else:
        zobs = calin.simulation.vs_cta.ctas_observation_level()
        atm = calin.simulation.vs_cta.ctas_atmosphere(quiet=True)
        atm_abs = calin.simulation.vs_cta.ctas_atmospheric_absorption(quiet=True)
        mst = calin.simulation.vs_cta.msts1_config()
    mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_x(0)
    mst.mutable_prescribed_array_layout().mutable_scope_positions(0).set_y(0)
    nscope = mst.prescribed_array_layout().scope_positions_size()

    # Instantiate array
    global telescope_layout
    global nchan
    array = calin.simulation.vs_optics.VSOArray()
    array.generateFromArrayParameters(mst)
    scope = array.telescope(0)
    telescope_layout = scope.convert_to_telescope_layout()
    nchan = telescope_layout.camera().channel_size()

    # Configure IACT array
    global iact
    iact_cfg = iact_class.default_config()
    if args.no_refraction:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_NO_RAYS)
    else:
        iact_cfg.set_refraction_mode(calin.ix.simulation.vcl_iact.REFRACT_ONLY_CLOSE_RAYS)
    iact = iact_class(atm, atm_abs, iact_cfg)

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

    # Add telescope arrays
    global all_pe_processor
    global all_prop
    global bmax_polynomial
    all_pe_processor = []
    all_prop = []
    if(len(args.bmax_polynomial) == 0):
        bmax_polynomial = numpy.asarray([0.0])
    else:
        bmax_polynomial = numpy.flipud(args.bmax_polynomial) * 100.0
    detector_type_name = 'MST/NC'
    for i in range(numpy.max([1,  args.reuse])):
        iact.add_propagator_set(numpy.flipud(bmax_polynomial), f"Super array {i}")
        pe_processor = calin.simulation.ray_processor.SimpleListPEProcessor(nscope,nchan)
        prop = iact.add_davies_cotton_propagator(mst, pe_processor, det_eff, cone_eff, pe_gen, args.tts, detector_type_name)
        all_pe_processor.append(pe_processor)
        all_prop.append(prop)

    # Set telescope pointing direction and viewcone parameters
    global pt_dir
    global vc_dir
    global viewcone_polynomial
    iact.point_all_telescopes_az_el_deg(args.az, args.el)
    el = args.el * numpy.pi/180.0
    az = args.az * numpy.pi/180.0
    pt_dir  = numpy.asarray([numpy.cos(el)*numpy.sin(az), numpy.cos(el)*numpy.cos(az), numpy.sin(el)])
    theta = args.theta * numpy.pi/180
    phi = args.phi * numpy.pi/180
    vc_dir= numpy.asarray([numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)])
    vc_dir = calin.math.geometry.rotate_vec_z_to_u_Rzy(vc_dir, -pt_dir)
    if(len(args.viewcone_polynomial)==0):
        viewcone_polynomial = numpy.asarray([0.0])
    else:
        viewcone_polynomial = numpy.flipud(args.viewcone_polynomial) * numpy.pi/180.0
    if not args.no_viewcone_cut:
        iact.set_viewcone_from_telescope_fields_of_view()

    # Magnetic field (if not disabled)
    global bfield
    if args.no_bfield:
        bfield = None
    else:
        wmm = calin.simulation.world_magnetic_model.WMM()
        bfield = wmm.field_vs_elevation(mst.array_origin().latitude(), mst.array_origin().longitude())

    # Configure Geant4 shower generator
    geant4_cfg = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator.customized_config(
        1000, 0, atm.top_of_atmosphere(), calin.simulation.geant4_shower_generator.VerbosityLevel_SUPRESSED_STDOUT)
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
        cfg.set_enable_ions(True)

    # Instantiate Geant4 shower generator
    global generator
    generator = calin.simulation.geant4_shower_generator.Geant4ShowerGenerator(atm, geant4_cfg, bfield);

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

    # Spectrum RNG transformation - convert uniform deviate to log10(energy)
    global spectral_transform
    if args.emax > args.emin:
        if(len(args.spectral_polynomial)==0):
            spectral_polynomial = numpy.asarray([1.0, 0.0])
        else:
            spectral_polynomial = numpy.flipud(args.spectral_polynomial)
            spectral_polynomial[-1] += 1.0
            spectral_polynomial = numpy.append(spectral_polynomial, [0.0])
        xmin = numpy.log10(args.emin)
        xmax = numpy.log10(args.emax)
        x = numpy.linspace(xmin, xmax, 10000) 
        y = 10**(numpy.polyval(spectral_polynomial, x))
        if len(bmax_polynomial) > 1:
            # If the impact parameter is not constant then must account for area
            y *= numpy.polyval(bmax_polynomial, x)**2
        if len(viewcone_polynomial) > 1:
            # If the viewcone is not constant then must account for solid angle
            y *= 1 - numpy.cos(numpy.polyval(viewcone_polynomial, x))
        p = numpy.append([0.0], numpy.cumsum(y[1:]+y[:-1])/2)
        p /= p[-1]
        if(0 < args.spectral_constant <= 1.0):
            p = p * (1-args.spectral_constant) + numpy.linspace(0.0,1.0,len(p))*args.spectral_constant
        spectral_transform = scipy.interpolate.interp1d(p, x, kind='linear', bounds_error=False, 
                                                        fill_value=(xmin,xmax))            
    else:
        x = numpy.log10(args.emin)
        spectral_transform = lambda r: x

    # Simulation configuration
    global sim_config
    sim_config = calin.ix.simulation.simulated_event.SimulationConfiguration()
    if args.primary == 'gamma':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.GAMMA)
    elif args.primary == 'muon':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.MUON)
    elif args.primary == 'electron':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.ELECTRON)
    elif args.primary == 'proton':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.PROTON)
    elif args.primary == 'helium':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.HELIUM)
    elif args.primary == 'iron':
        sim_config.set_particle_type(calin.ix.simulation.simulated_event.IRON)
    sim_config.set_energy_lo(args.emin)
    if(args.emax > args.emin):
        sim_config.set_energy_hi(args.emax)
        sim_config.set_energy_spectrum_polynomial(args.spectral_polynomial)
    else:
        sim_config.set_energy_hi(args.emin)
    sim_config.set_elevation(args.el)
    sim_config.set_azimuth(args.az)
    sim_config.set_theta(args.theta)
    sim_config.set_phi(args.phi)
    sim_config.mutable_primary_axis().set_x(vc_dir[0])
    sim_config.mutable_primary_axis().set_y(vc_dir[1])
    sim_config.mutable_primary_axis().set_z(vc_dir[2])
    sim_config.set_viewcone_halfangle_polynomial(numpy.flipud(viewcone_polynomial) * 180.0/numpy.pi)
    sim_config.set_scattering_radius_polynomial(numpy.flipud(bmax_polynomial)*0.01)
    sim_config.set_banner(get_banner())
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
        if(atm_abs_zmin < level.z < atm_abs_zmax):
            proto_level.set_optical_depth_1d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 1.5))
            proto_level.set_optical_depth_2d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.0))
            proto_level.set_optical_depth_2d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 2.5))
            proto_level.set_optical_depth_3d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.0))
            proto_level.set_optical_depth_3d5ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 3.5))
            proto_level.set_optical_depth_4d0ev(atm_abs.optical_depth_for_altitude_and_energy(level.z, 4.0))
    sim_config.mutable_geant4_shower_generator_config().CopyFrom(geant4_cfg)
    sim_config.mutable_iact_array_config().CopyFrom(iact_cfg)
    for ipropagatorset in range(iact.num_propagator_sets()):
        array_config = sim_config.add_detector_array_config()
        array_config.set_name(iact.propagator_set_name(ipropagatorset))
        array_config.set_scattering_radius_polynomial(iact.scattering_radius_polynomial(ipropagatorset))
        for ipropagator in range(iact.propagator_set_size(ipropagatorset)):
            propagator = iact.propagator_set_dc_element(ipropagatorset,ipropagator)
            dc_array = None
            if propagator:
                dc_array = propagator.array()
            else:
                propagator = iact.propagator_set_element(ipropagatorset,ipropagator)
            detector_spheres = propagator.detector_spheres()
            group_config = array_config.add_detector_group_config()
            group_config.set_type_name(iact.propagator_set_element_name(ipropagatorset,ipropagator))
            group_config.set_ndetector(len(detector_spheres))
            group_config.set_propagator_banner(propagator.banner())
            if(dc_array):
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
    detector_type_config.set_pe_time_spread(args.tts)
    detector_type_config.mutable_dc_array_parameters().CopyFrom(mst)

    args_dict = vars(args)
    for arg in args_dict:
        val = str(args_dict[arg])
        sim_config.set_command_line_args(arg, val)

    signal.signal(signal.SIGINT, handle_sigint_quietly)

def gen_event(args):
    x = spectral_transform(numpy.random.uniform()) # x=log10(E/1TeV)
    e = 10**(x+6) # Convert TeV to MeV

    vc_halfangle = numpy.polyval(viewcone_polynomial, x)
    costheta = 1.0 - (1.0 - numpy.cos(vc_halfangle))*numpy.random.uniform()
    theta = numpy.arccos(costheta)
    phi = numpy.random.uniform() * 2*numpy.pi
    u = numpy.asarray([numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)])
    u = calin.math.geometry.rotate_vec_z_to_u_Rzy(u, vc_dir)

    ct0 = 1.0/u[2]*(atm.top_of_atmosphere() - atm.zobs(0))
    x0 = numpy.asarray([0,0,atm.zobs(0)]) + ct0*u

    generator.generate_showers(iact, 1, particle_type, e, x0, u, ct0)

    sim_event = calin.ix.simulation.simulated_event.SimulatedEvent()
    iact.save_to_simulated_event(sim_event, store_pe_weights, store_times_as_integer)

    sim_event.mutable_viewcone_axis().set_x(vc_dir[0])
    sim_event.mutable_viewcone_axis().set_y(vc_dir[1])
    sim_event.mutable_viewcone_axis().set_z(vc_dir[2])
    sim_event.set_viewcone_costheta(costheta)

    return sim_event

def gen_batch(filename):
    global instance_event_id
    if stop_requested:
        return None
    # Write simulation configuration
    h5 = type(sim_config).NewHDFStreamWriter(filename, 'config', truncate=True)
    h5.write(sim_config)
    del h5
    h5 = calin.ix.simulation.simulated_event.SimulatedEvent.NewHDFStreamWriter(filename, 'events', truncate=False)
    num_events = 0
    num_tracks = 0
    num_steps = 0
    num_rays = 0
    
    for i in range(block_size):
        if stop_requested:
            break
        sim_event = gen_event(args)
        sim_event.set_instance_id(instance_id)
        sim_event.set_event_id(instance_event_id)
        instance_event_id += 1
        num_events += 1
        h5.write(sim_event)
        num_tracks += iact.num_tracks()
        num_steps += iact.num_steps()
        num_rays += iact.num_rays()
    del h5
    # Write provenance anthology
    provenance_antology = calin.provenance.anthology.get_current_anthology()    
    h5 = type(provenance_antology).NewHDFStreamWriter(filename, 'provenance', truncate=False)
    h5.write(provenance_antology)
    del h5
    num_bytes = os.path.getsize(filename)
    return filename, num_events, num_tracks, num_steps, num_rays, num_bytes

def format_energy(etev):
    if(etev < 0.1):
        return f'{etev*1000:,.2f} GeV'    
    elif(etev < 1.0):
        return f'{etev*1000:,.1f} GeV'    
    elif(etev<10):
        return f'{etev:,.3f} TeV'
    elif(etev<100):
        return f'{etev:,.2f} TeV'
    else:
        return f'{etev:,.1f} TeV'

def get_banner(sleep_time = 0):
    banner = "Command line arguments :\n"
    args_dict = vars(saved_args)
    for arg in args_dict:
        banner += f'- {arg}: {args_dict[arg]}\n'
    banner += '\nIACT configuration :\n'
    banner += iact.banner() + '\n'
    banner += '\nPrimary viewcone :'
    if len(viewcone_polynomial)==0 or (len(viewcone_polynomial)==1 and viewcone_polynomial[0]==0):
        banner += " DISABLED\n"
    elif len(viewcone_polynomial)==1:
        banner += f' FIXED - {viewcone_polynomial/numpy.pi*180.0:.1f} degrees\n'
    else:
        banner += " VARIABLE\n"
        banner += "  10GeV, 100GeV, 1TeV, 10TeV, 100TeV :"
        for i,x in enumerate([-2,-1,0,1,2]):
            banner += f'{", " if i>0 else " "}{numpy.polyval(viewcone_polynomial,x)/numpy.pi*180.0:.1f}'
        banner += " degrees\n"

    banner += '\nShower energy :'
    if saved_args.emin >= saved_args.emax:
        banner += f' MONOCHROMATIC at {format_energy(saved_args.emin)}'
    else:
        banner += ' SPECTRUM\n'
        x0 = spectral_transform(0)
        p0 = 0
        for p1 in 1-0.5**numpy.arange(1,10):
            x1 = spectral_transform(p1)
            banner += f'  {format_energy(10**x0)} to {format_energy(10**x1)} : {(p1-p0)*100:.4f} %\n'
            x0 = x1
            p0 = p1
    p1 = 1.0
    x1 = spectral_transform(p1)
    banner += f'  {format_energy(10**x0)} to {format_energy(10**x1)} : {(p1-p0)*100:.4f} %'

    if sleep_time>0:
        time.sleep(sleep_time)
    return banner

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

def print_line(filename):    
    runtime = (datetime.datetime.now(datetime.timezone.utc) - begin_utc).total_seconds()
    line = f'{filename}: {num_events:,d}'
    fraction = num_events/(args.n*args.block_size) if args.n > 0 else 0.0
    if args.n > 0:
        line += f' / {args.n*args.block_size:,d} = {fraction*100.0:.2f}%'
    line += f' ; {format_duration(runtime)}'
    if args.n > 0:
        totaltime = runtime/fraction
        line += f' / {format_duration(totaltime)}'
    line += f' ({num_events/runtime:,.1f} Hz)'
    line += f' ; {num_rays:,d} rays ; {num_rays/num_steps:.2f} {num_steps/num_tracks:.2f}';
    line += f' ; {num_bytes*1e-9:,.3f}'
    if args.n > 0:
        totalsize = num_bytes/fraction
        line += f' / {totalsize*1e-9:,.3f}'
    line += ' GB'
    if stop_requested:
        line += ' (STOPPING)'
    print(line)
    
    if num_batch % 50 == 0 and num_batch != args.n and not stop_requested:
        if args.emin < args.emax:
            spectline = f'Spectrum: {format_energy(args.emin)} - {format_energy(args.emax)}, Index: {",".join([str(x) for x in args.spectral_polynomial or [0.0]])}'
        else:
            spectline = f'Monochromatic: {format_energy(args.emin)}'

        if len(numpy.flipud(args.bmax_polynomial))>1 and args.emin < args.emax:
            bmaxline = f'{numpy.polyval(numpy.flipud(args.bmax_polynomial), numpy.log10(args.emin)):,.1f} - ' \
                + f'{numpy.polyval(numpy.flipud(args.bmax_polynomial), numpy.log10(args.emax)):,.1f} m'
        elif len(args.bmax_polynomial)>0:
            bmaxline = f'{numpy.polyval(numpy.flipud(args.bmax_polynomial), numpy.log10(args.emin)):,.1f} m'
        else:
            bmaxline = '0.0 m'

        if len(numpy.flipud(args.viewcone_polynomial))>1 and args.emin < args.emax:
            vcline = f'{numpy.polyval(numpy.flipud(args.viewcone_polynomial), numpy.log10(args.emin)):,.1f} - ' \
                + f'{numpy.polyval(numpy.flipud(args.viewcone_polynomial), numpy.log10(args.emax)):,.1f} deg'
        elif len(args.viewcone_polynomial)>0:
            vcline = f'{numpy.polyval(numpy.flipud(args.viewcone_polynomial), numpy.log10(args.emin)):,.1f} deg'
        else:
            vcline = '0.0 deg'

        filesline = f'{num_batch}:,d'
        if args.n > 0:
            filesline += f' / {args.n;,d}'

        print(f'\n===== Particle: {args.primary} ; Site: {args.site} ; El: {args.el:.1f}, Az: {args.az:.1f} ; {spectline} ; Bmax: {bmaxline} ; Viewcone: {vcline} ; Files: {filesline} =====\n')
    
def process_results(results):
    if results is None:
        return
    
    global num_batch
    global num_events
    global num_tracks
    global num_steps
    global num_rays
    global num_bytes

    filename, nevents, ntracks, nsteps, nrays, nbytes = results
    num_batch  += 1
    num_events += nevents
    num_tracks += ntracks
    num_steps  += nsteps
    num_rays   += nrays
    num_bytes  += nbytes
    print_line(filename)

def unique_filename(template, seq_num):
    id = f'{numpy.random.bit_generator.randbits(128):032X}'
    return template.format(id=id, seq=seq_num)

def handle_sigint(signum, frame):
    global stop_requested
    if not stop_requested:
        print("\nStop requested - terminating current run")
    stop_requested = True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Nectarcam single-telescope shower simulation')
    
    parser.add_argument('-n', type=int, default=0,
        help='Specify the number of file blocks to simulate, or zero to simulate indefinitely until killed')
    parser.add_argument('--block_size', type=int, default=1000,
        help='Specify the number of shower events per file block')
    parser.add_argument('--reuse', type=int, default=10,
        help='Specify the number of times to reuse each shower')

    parser.add_argument('-o', '--output', type=str, default='events_{id}.h5',
        help='Output filename. May contain the token "{id}", which will be replaced by a random identifier, or "{seq}" giving a sequential number (default: events_{id}.h5)')
    parser.add_argument('--store_times_as_integer', action='store_true',
        help='Store times as 16-bit integers in output file, where possible, with 1-ps resulution (default: disabled)')


    parser.add_argument('--site', type=str, default='ctan', choices=['ctan','ctas'],
        help='Site to simulate (default: ctan)')

    # Shower impact parameter parameters
    parser.add_argument('-b', '--bmax_polynomial', type=comma_separated_floats, default=[1000.0],
        help='Specify a polynomial defining maximum shower impact parameter in meters as a function of energy in TeV. If x=log10(E/1TeV) then b(x)=P0+P1*x+P2*x^2+... The polynomial should be specified as a comma-separated list of coefficients, e.g. P0,P1,P2... (default: 1000.0)')

    # Pointing and viewcone parameters
    parser.add_argument('--viewcone_polynomial', type=comma_separated_floats, default=[0.0],
        help='Specify a polynomial defining viewcone half-angle in degrees as a function of energy in TeV. If x=log10(E/1TeV) then theta(x)=P0+P1*x+P2*x^2+... The polynomial should be specified as a comma-separated list of coefficients, e.g. P0,P1,P2... (default: 0.0)')
    parser.add_argument('--az', type=float, default=0.0,
        help='Specify the telescope azimuth angle in degrees')
    parser.add_argument('--el', type=float, default=70.0,
        help='Specify the telescope elevation angle in degrees')
    parser.add_argument('--theta', type=float, default=0.0,
        help='Specify the fixed offset between the primary direction and the telescope pointing direction in degrees')
    parser.add_argument('--phi', type=float, default=0.0,
        help='Specify the fixed polar angle of the primary direction around the telescope pointing direction in degrees')
    parser.add_argument('--no_viewcone_cut', action='store_true', 
        help='Do not generate photons on tracks that are outside the detector pointing viewcone (default: enabled)')

    # Particle type and energy
    parser.add_argument('-p', '--primary', type=str, default='gamma', 
        choices=['gamma','muon','electron','proton','helium','iron'],
        help='Specify the primary particle type')
    parser.add_argument('--emin', type=float, default=0.01,
        help='Specify the lower bound on the energy spectrum in TeV')
    parser.add_argument('--emax', type=float, default=100.0,
        help='Specify the upper bound on the energy spectrum in TeV')
    parser.add_argument('--spectral_polynomial', type=comma_separated_floats, default=[-2.7],
        help='Specify the spectral shape as a polynomial in log10(E/1TeV). If x=log10(E/1TeV) then dN/dE(x)=10^(P0*x+P1*x^2+...). The polynomial should be specified as a comma-separated list of coefficients, e.g. P0,P1,P2... (default: -2.7)')
    parser.add_argument('--spectral_constant', type=float, default=0.0,
        help='Fraction of simulated events to produce using flat spectrum in dN/dlog(E).')

    # PE spectrum and transit time spread
    parser.add_argument('--tts', type=float, default=0.00,
        help='Specify the transit time spread RMS in ns.')
    parser.add_argument('--enable_pe_spectrum', action='store_true',
        help='Enable PE charge spectrum and store PE charges in output file (default: disabled)')

    # Low-level simulation options
    parser.add_argument('--no_bfield', action='store_true', 
        help='Disable the magnetic field (default: enabled)')
    parser.add_argument('--no_refraction', action='store_true', 
        help='Disable refraction of rays in the atmosphere (default: enabled)')
    parser.add_argument('--multiple_scattering', type=str, default='normal',
        choices=['minimal','simple','normal','better','insane'],
        help='Specify the multiple scattering model (default: normal)')

    parser.add_argument('--nthread', type=int, default=0,
        help='Number of threads to use (default: 0 = number of CPUs available)')
    parser.add_argument('--avx', type=int, default=512, choices=[128,256,512],
        help='Set the AVX vector size in bits (default: 512)')

    global args
    args = parser.parse_args()

    max_workers = args.nthread or os.cpu_count() or 1

    global begin_utc
    global num_queued
    global num_batch
    global num_events
    global num_tracks
    global num_steps
    global num_rays
    global num_bytes
    global stop_requested

    begin_utc = datetime.datetime.now(datetime.timezone.utc)
    num_queued = 0
    num_batch = 0
    num_events = 0
    num_rays = 0
    num_steps = 0
    num_tracks = 0
    num_bytes = 0
    stop_requested = False

    # Run the simulations in this thread
    if max_workers == 1:
        init(args)
        print(get_banner())
        print()
        signal.signal(signal.SIGINT, handle_sigint)
        while (args.n==0 or num_batch<args.n) and not stop_requested:
            results = gen_batch(unique_filename(args.output, num_batch))
            process_results(results)
    else:
        # Use a process pool for parallelism
        batch_size = 4 * max_workers
        print(f'Launching {max_workers} workers in process pool')
        with concurrent.futures.ProcessPoolExecutor(initializer=init, initargs=(args,), max_workers=max_workers) as executor:
            futures = set()

            # This crazyness is here to launch all workers and wait for them to print
            # there initialization text before we print the banner for the simulation
            for i in range(max_workers):
                futures.add(executor.submit(get_banner,5)) # Sleep 5 seconds to allow all threads to start
            while futures:
                done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    futures.remove(fut)
                    banner = fut.result()
            if banner is not None:
                print(banner)
                print()
            
            signal.signal(signal.SIGINT, handle_sigint)

            remaining = args.n

            # Submit initial batch
            first_batch = min(batch_size, remaining or batch_size)
            for i in range(first_batch):
                futures.add(executor.submit(gen_batch, unique_filename(args.output, num_queued)))
                num_queued += 1
            if remaining:
                remaining -= first_batch

            # Keep a bounded number of in-flight futures; as one completes, submit another
            while futures:
                done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    futures.remove(fut)
                    num_queued -= 1
                    process_results(fut.result())
                    if (args.n == 0 or remaining > 0) and not stop_requested:
                        futures.add(executor.submit(gen_batch, unique_filename(args.output, num_batch+num_queued)))
                        num_queued += 1
                        if(remaining):
                            remaining -= 1
