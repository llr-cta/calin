# Example usage: python3.12 shower_trigger_threshold.py -o test.pickle -n 6000000 -b 6500 -e 0.0175 --nsb 0.0 -p proton --viewcone 13.0 --write_batch=100000 --omit_untriggered --reuse 200 --enable_viewcone_cut

# calin/scripts/shower_trigger_threshold.py -- Stephen Fegan - 2026-01-15
#
# Simulate mono-energetic events with GEANT4 and calculate trigger threshold.
# Write results to a pickle file.
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
import concurrent
import pickle
import datetime
import platform
import numpy
import calin.math.geometry
import calin.simulation.vs_cta
import calin.simulation.vs_panoseti
import calin.simulation.ray_processor
import calin.simulation.iact_factory
import calin.iact_data.instrument_layout
import calin.iact_data.nectarcam_layout

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Nectarcam single-telescope shower trigger threshold calculation')
parser.add_argument('-n', type=int, default=1000,
                   help='Specify the number of showers to simulate')
parser.add_argument('--reuse', type=int, default=10,
                   help='Specify the number of times to reuse each shower')
parser.add_argument('--write_batch', type=int, default=100,
                   help='Specify the number of events to write per batch to output file')

parser.add_argument('-o', '--output', type=str, default='tt.pickle',
                    help='Write trigger thresholds to this file')
parser.add_argument('--omit_untriggered', action='store_true',
                    help='Reduce file size by omitting untriggered events')

parser.add_argument('--site', type=str, default='ctan', choices=['ctan','ctas','dark100'],
                    help='Site to simulate (default: ctan)')
parser.add_argument('-b', '--bmax', type=float, default=1000.0,
                   help='Specify the maximum shower impact parameter in meters')
parser.add_argument('--tts', type=float, default=0.75,
                   help='Specify the transit time spread RMS in ns')
parser.add_argument('--az', type=float, default=0.0,
                   help='Specify the telescope azimuth angle in degrees')
parser.add_argument('--el', type=float, default=70.0,
                   help='Specify the telescope elevation angle in degrees')
parser.add_argument('-e', '--energy', type=float, default=0.3,
                   help='Specify the telescope energy in TeV')
parser.add_argument('-p', '--primary', type=str, default='gamma', 
                    choices=['gamma','muon','electron','proton','helium','iron'],
                    help='Specify the primary particle type')
parser.add_argument('--viewcone', type=float, default=0.0,
                   help='Specify the random sampling offset between the primary direction and the telescope pointing direction in degrees')
parser.add_argument('--theta', type=float, default=0.0,
                   help='Specify the fixed offset between the primary direction and the telescope pointing direction in degrees')
parser.add_argument('--phi', type=float, default=0.0,
                   help='Specify the fixed polar angle of the primary direction around the telescope pointing direction in degrees')
parser.add_argument('--no_bfield', action='store_true', 
                   help='Disable the magnetic field (default: enabled)')
parser.add_argument('--no_refraction', action='store_true', 
                   help='Disable refraction of rays in the atmosphere (default: enabled)')
parser.add_argument('--enable_viewcone_cut', action='store_true', 
                   help='Do not generate photons on tracks that are outside the viewcone (default: disabled)')

parser.add_argument('--multiple_scattering', type=str, 
                    choices=['minimal','simple','normal','better','insane'],
                    default='normal',
                    help='Specify the multiple scattering model (default: normal)')
parser.add_argument('--nsb', type=float, default=0.30,
                   help='Specify the NSB rate in GHz')
parser.add_argument('--noise', action='store_true', help='Add electronics noise')
parser.add_argument('--coupling', choices=['ac', 'dc'], default='ac',
                    help='Camera coupling mode (default: ac)')
parser.add_argument('--no_after_pulsing', action='store_true', help='Disable after-pulsing in SPE spectrum (default: enabled)')
parser.add_argument('-t', '--trigger', type=str, default='3nn', choices=['3nn','4nn','m2','m3','m4','multiplicity'],
                    help='Trigger algorithm to use (default: 3nn)')
parser.add_argument('-m', '--multiplicity', type=int, default=3,
                    help='Channel multiplicity if "multiplicity" algorithm is selected')
parser.add_argument('-c', '--coincidence', type=int, default=24,
                    help='Set the L1 trigger coincidence time in samples')

parser.add_argument('--threshold_min', type=float, default=40.0,
                    help='Minimum threshold for search (default: 40.0)')

parser.add_argument('--nthread', type=int, default=0,
                    help='Number of threads to use (default: 0 = number of CPUs available)')
parser.add_argument('--avx', type=int, default=512, choices=[128,256,512],
                    help='Set the AVX vector size in bits (default: 512)')

args = parser.parse_args()

tcoincidence = args.coincidence
iact = None
has_one_event = False

begin_utc = datetime.datetime.now(datetime.timezone.utc)

# Prepare a JSON header describing the run
config = vars(args).copy()
config['_begin_utc'] = begin_utc.isoformat()
config['_host'] = platform.node()
config['_num_events'] = 0
config['_num_tracks'] = 0
config['_num_steps'] = 0
config['_num_rays'] = 0

def init():
    numpy.random.seed()

    # Select simulation classes based on AVX size requested
    if args.avx == 128:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat128
    elif args.avx == 256:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat256
    else:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat512

    # Load site environment (atmosphere, absorption, efficiencies, B-field)
    global site_env
    site_env = calin.simulation.iact_factory.load_site_environment(
        args.site, enable_pe_spectrum=True, no_bfield=args.no_bfield,
        quiet=True, avx=args.avx)

    # Set up telescope array (site-specific; PANOSETI keeps real positions)
    array_params, nscope, nchan_from_array, detector_type_name = \
        calin.simulation.iact_factory.setup_telescope_array(args.site, el_deg=args.el)

    # Use the site camera channel count, and retain every telescope in the array.
    global nchan
    global nscope
    global scam
    if args.site in ('ctan', 'ctas'):
        ncam = calin.iact_data.nectarcam_layout.nectarcam_layout()
        scam = calin.iact_data.instrument_layout.reorder_camera_channels(
            ncam, ncam.pixel_spiral_channel_index())
        nchan = scam.channel_size()
    else:
        scam = None
        nchan = nchan_from_array

    # Configure IACT array
    global iact
    iact, _, _ = calin.simulation.iact_factory.create_iact_array(
        site_env.atm, site_env.atm_abs, avx=args.avx,
        no_refraction=args.no_refraction)

    # Load the selected site/shaper impulse response for waveform simulation.
    global isample0
    global dtsample
    global nsample
    if args.site == 'dark100':
        pulse = calin.simulation.vs_panoseti.panoseti_fast_impulse_response()
    else:
        pulse = calin.simulation.vs_cta.mstn_impulse_response()
    hg = pulse['hg']
    dtsample = pulse['dt']
    isample0 = len(hg)
    nsample = 1 << ((len(hg)-1).bit_length() + 1)  # Next power of two greater than 2*len(hg)

    # Each propagator reuse needs its own PE list and waveform state.
    global all_electronics_sim
    all_electronics_sim = []

    def make_waveform_pe_processor(ns, nc):
        electronics_sim = electronics_sim_class(ns, nc, nsample, dtsample, isample0)
        all_electronics_sim.append(electronics_sim)
        return electronics_sim

    # Attach propagators
    global all_pe_processor
    global all_prop
    all_pe_processor, all_prop = calin.simulation.iact_factory.attach_iact_propagators(
        iact, args.site, array_params,
        bmax_polynomial=args.bmax,   # scalar in metres
        reuse=args.reuse, nscope=nscope, nchan=nchan,
        det_eff=site_env.det_eff, cone_eff=site_env.cone_eff,
        pe_gen=site_env.pe_gen,
        tts=args.tts,
        lens_spline=site_env.lens_refractive_index_spline,
        detector_type_name=detector_type_name,
        pe_processor_factory=make_waveform_pe_processor)

    # Set telescope pointing direction and optional viewcone cut.
    global pt_dir
    pt_dir, _ = calin.simulation.iact_factory.set_iact_pointing(
        iact, args.el, args.az, apply_viewcone_cut=args.enable_viewcone_cut)

    # Configure and instantiate Geant4 shower generator
    global generator
    generator, _ = calin.simulation.iact_factory.create_geant4_generator(
        site_env.atm, bfield=site_env.bfield,
        multiple_scattering=args.multiple_scattering, primary=args.primary)
    generator.set_minimum_energy_cut(20)  # 20 MeV cut on KE (e-,p+,n,ions) or Etot

    # Select trigger algorithm once, then configure each independent reuse slot.
    global trigger_method
    if args.trigger in ('multiplicity', 'm2', 'm3', 'm4'):
        trigger_method = ('trigger_panoseti_multiplicity_cr'
            if args.site == 'dark100' else 'trigger_multiplicity_cr')
    elif args.trigger == '3nn':
        trigger_method = 'trigger_3nn_cr'
    elif args.trigger == '4nn':
        trigger_method = 'trigger_4nn_cr'
    else:
        raise ValueError(f'Unknown trigger algorithm: {args.trigger}')

    neighbors = None
    if args.site in ('ctan', 'ctas'):
        neighbors = scam
    elif args.trigger in ('3nn', '4nn'):
        # Dark100 pixels form a square grid; use four-connected neighbors.
        npixel_axis = int(round(numpy.sqrt(nchan)))
        if npixel_axis * npixel_axis != nchan:
            raise ValueError('Dark100 camera channel count is not a square grid')
        neighbors = numpy.full((4, nchan), -1, dtype=numpy.int32)
        for iy in range(npixel_axis):
            for ix in range(npixel_axis):
                ichan = iy * npixel_axis + ix
                adjacent = []
                if ix > 0:
                    adjacent.append(ichan - 1)
                if ix + 1 < npixel_axis:
                    adjacent.append(ichan + 1)
                if iy > 0:
                    adjacent.append(ichan - npixel_axis)
                if iy + 1 < npixel_axis:
                    adjacent.append(ichan + npixel_axis)
                neighbors[:len(adjacent), ichan] = adjacent
    for electronics_sim in all_electronics_sim:
        electronics_sim.register_impulse_response(hg, 'DC')
        electronics_sim.add_camera_response(
            numpy.asarray([0]), args.coupling == 'ac')
        if neighbors is not None:
            electronics_sim.set_cr_neighbors(0, neighbors)

        if trigger_method in ('trigger_multiplicity_cr', 'trigger_panoseti_multiplicity_cr'):
            multiplicity = {
                'multiplicity': args.multiplicity,
                'm2': 2,
                'm3': 3,
                'm4': 4,
            }[args.trigger]
            electronics_sim.set_cr_multiplicity(0, multiplicity)

        if args.no_after_pulsing:
            pe_gen_nsb = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, quiet=True)
        else:
            pe_gen_nsb = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, afterpulsing=True, quiet=True)
        pe_gen_nsb.this.disown()  # Let this electronics processor own it

        if args.nsb > 0:
            nsb = numpy.zeros(nchan) + args.nsb
            electronics_sim.set_cr_nsb_rate(0, nsb, pe_gen_nsb, True)

        if args.noise:
            def noise(f):
                fhi = 330.0
                return (1+0.2*(f/275.0)**2)*(numpy.tanh((fhi-f)/20.0)+1)*0.9 + 0.6
            freq = electronics_sim.spectral_frequencies_ghz(False)
            noise_spectrum = noise(freq*1000)
            noise_spectrum[0] = 0
            noise_spectrum *= numpy.sqrt(10.5**2/16/electronics_sim.noise_spectrum_var(noise_spectrum))
            electronics_sim.set_cr_noise_spectrum(0, noise_spectrum)

def gen_event():
    e = args.energy * 1e6 # Convert TeV to MeV

    costheta = 1.0 - (1.0 - numpy.cos(args.viewcone * numpy.pi/180))*numpy.random.uniform()
    theta = numpy.arccos(costheta)
    phi = numpy.random.uniform() * 2*numpy.pi
    u = numpy.asarray([numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)])

    theta = args.theta * numpy.pi/180
    phi = args.phi * numpy.pi/180
    v = numpy.asarray([numpy.sin(theta)*numpy.cos(phi), numpy.sin(theta)*numpy.sin(phi), numpy.cos(theta)])
    u = calin.math.geometry.rotate_vec_z_to_u_Rzy(u, v)
    u = calin.math.geometry.rotate_vec_z_to_u_Rzy(u, -pt_dir)

    x0 = numpy.asarray([0,0,site_env.atm.zobs(0)]) + u/u[2]*(
        site_env.atm.top_of_atmosphere() - site_env.atm.zobs(0))
    pt = calin.simulation.iact_factory.get_tracker_particle_type(args.primary)

    generator.generate_showers(iact, 1, pt, e, x0, u)

    return e,pt,u,x0,costheta

def find_threshold(iarray, iscope):
    # Clear previous waveforms, transfer the PEs, add NSB, convolve impulse response
    electronics_sim = all_electronics_sim[iarray]
    electronics_sim.clear_waveforms()
    electronics_sim.transfer_scope_pes_to_waveform(all_pe_processor[iarray], iscope)
    if args.nsb>0:
        electronics_sim.add_nsb_noise_to_waveform_cr(0)
    electronics_sim.convolve_impulse_response_fftw_codelet_cr(0)

    # Now perform threshold search
    trigger_fn = getattr(electronics_sim, trigger_method)
    
    # Start at threshold_min
    threshold = args.threshold_min
    electronics_sim.set_cr_threshold(0, numpy.zeros(nchan) + threshold, tcoincidence)
    itrig = trigger_fn(0, isample0)
    if itrig == -1:
        return -1  # No trigger even at min threshold
    
    # Double until it fails, updating lower bound
    lower = threshold
    upper = threshold
    while True:
        upper *= 2
        electronics_sim.set_cr_threshold(0, numpy.zeros(nchan) + upper, tcoincidence)
        itrig = trigger_fn(0, isample0)
        if itrig == -1:
            break
        lower = upper
    while (upper - lower) / lower > 0.01:
        mid = (lower + upper) / 2
        electronics_sim.set_cr_threshold(0, numpy.zeros(nchan) + mid, tcoincidence)
        itrig = trigger_fn(0, isample0)
        if itrig == -1:
            upper = mid
        else:
            lower = mid
    
    return lower

def one_event():
    global has_one_event
    event_results = []
    try:
        e,pt,u,x0,costheta = gen_event()
        for iarray in range(args.reuse):
            # A common threshold triggers the array if any telescope triggers.
            # The OR trigger boundary is the largest per-telescope limit.
            scope_thresholds = [find_threshold(iarray, iscope)
                                for iscope in range(nscope)]
            threshold = max(scope_thresholds)
            event_results.append(dict(
                iarray         = iarray,
                e              = e,
                pt             = int(pt),
                u0             = u.tolist(),
                x0             = x0.tolist(),
                costheta       = costheta,
                b              = iact.scattered_distance(iarray),
                offset         = iact.scattered_offset(iarray)[0:2].tolist(),
                threshold      = threshold))
    except Exception as ex:
        print(f'Error simulating event: {ex}')
        raise
    if not has_one_event:
        event_results[0]['_banner'] = iact.banner()
        has_one_event = True
    event_results[0]['_num_tracks'] = iact.num_tracks()
    event_results[0]['_num_steps'] = iact.num_steps()
    event_results[0]['_num_rays'] = iact.num_rays()
    return event_results

def save_results(results, num_events, filehandle):
    now_utc = datetime.datetime.now(datetime.timezone.utc)
    config['_num_events'] = num_events
    config['_end_utc'] = now_utc.isoformat()
    config['_run_time'] = (now_utc - begin_utc).total_seconds()
    output = dict(
        config = config,
        results = results)
    pickle.dump(output, filehandle)
    filehandle. flush()

num_events = 0
num_rays = 0
num_steps = 0
num_tracks = 0
events_written = 0
batch_start = 0
all_results = []

def print_line():    
    print(f'{args.output}: {len(all_results)} ;',
            f'{num_events:,d} / {args.n*args.reuse:,d} =',
            f'{num_events/(args.n*args.reuse)*100:.2f} % ;',
            f'{config["_run_time"]/3600:.2f} /',
            f'{args.n*args.reuse/num_events*config["_run_time"]/3600:.2f} hr ;',
            f'{num_events/config["_run_time"]:,.2f} Hz ;',
            f'{num_rays:,d} rays ;',
            f'{num_rays/max(num_steps,1):.2f}',
            f'{num_steps/max(num_tracks,1):.2f}')
    
def process_results(results):
    global num_events
    global events_written
    global batch_start
    global all_results
    global config
    global num_rays
    global num_steps
    global num_tracks

    for r in results:
        if '_banner' in r:
            config['_banner'] = r['_banner']
            del r['_banner']
        if '_num_tracks' in r:
            config['_num_tracks'] += r['_num_tracks']
            config['_num_steps'] += r['_num_steps']
            config['_num_rays'] += r['_num_rays']
            num_tracks += r['_num_tracks']
            num_steps += r['_num_steps']
            num_rays += r['_num_rays']
            del r['_num_tracks']
            del r['_num_steps']
            del r['_num_rays']
        num_events += 1
        if args.omit_untriggered and r['threshold'] < 0:
            continue
        all_results.append(r)

    if (num_events-events_written)>=args.write_batch:
        save_results(all_results[batch_start:], num_events-events_written, f)
        print_line()
        events_written = num_events
        batch_start = len(all_results)
        config['_num_tracks'] = 0
        config['_num_steps'] = 0
        config['_num_rays'] = 0

max_workers = args.nthread or os.cpu_count() or 1
with open(args.output, 'wb') as f:
    # Run the simulations in this thread
    if max_workers == 1:
        init()
        print(iact.banner())
        for _ in range(args.n):
            process_results(one_event())
    else:
        # Use a process pool for parallelism
        batch_size = max_workers * 10
        with concurrent.futures.ProcessPoolExecutor(initializer=init, max_workers=max_workers) as executor:
            remaining = args.n
            futures = set()

            # Submit initial batch
            first_batch = min(batch_size, remaining)
            for _ in range(first_batch):
                futures.add(executor.submit(one_event))
            remaining -= first_batch

            # Keep a bounded number of in-flight futures; as one completes, submit another
            while futures:
                done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    futures.remove(fut)
                    process_results(fut.result())
                    if remaining > 0:
                        futures.add(executor.submit(one_event))
                        remaining -= 1

    if num_events>events_written:
        save_results(all_results[batch_start:], num_events-events_written, f)
        print_line()
        events_written = num_events
        batch_start = len(all_results)

th = args.threshold_min
pc = 98
bsim = config['bmax']
filtered_results = all_results
while True:
    filtered_results = [r for r in filtered_results if r['threshold'] >= th]
    ntrig = len(filtered_results)
    if(ntrig < 10 or th > 10*config['threshold_min']):
        break
    if(ntrig < 0.9*num_events):
        bmax = numpy.percentile([r['b'] for r in filtered_results],pc)
        thmax = numpy.percentile([numpy.arccos(r['costheta'])/numpy.pi*180.0 for r in filtered_results],pc)        
        print(f'{th:6.1f} {th/20:6.2f} | {ntrig:6d} {ntrig/num_events:5.3f} | {ntrig/num_events*numpy.pi*bsim**2:.3e} {bmax*1e-2:6.1f} | {thmax:6.4f}')
    th = th + 10.0
