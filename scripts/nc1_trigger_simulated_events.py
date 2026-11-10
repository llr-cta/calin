# Example usage: python3.12 shower_trigger_threshold.py -o test.pickle -n 6000000 -b 6500 -e 0.0175 --nsb 0.0 -p proton --viewcone 13.0 --write_batch=100000 --omit_untriggered --reuse 200 --enable_viewcone_cut

# calin/scripts/nc_trigger_simulated_events.py -- Stephen Fegan - 2026-02-19
#
# Read simulated single-telescope NectarCam events from HDF5 files, add noise and 
# electronice response, apply trigger logic and write passing events to a pickle file
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

import glob
import os
import argparse
import pickle
import concurrent
import datetime
import numpy
import calin.simulation.vs_cta
import calin.simulation.ray_processor
import calin.ix.simulation.simulated_event
import calin.iact_data.instrument_layout
import calin.iact_data.nectarcam_layout

def comma_separated_floats(value):
    try:
        return [float(x) for x in value.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a comma-separated list of floats (e.g. 1.0,2.3,4.2)")

def handle_sigint_quietly(signum, frame):
    global stop_requested
    stop_requested = True

def init(args):
    # Select simulation classes based on AVX size requested
    if args.avx == 128:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat128
    elif args.avx == 256:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat256
    else:
        electronics_sim_class = calin.simulation.ray_processor.VCLWaveformPEProcessorFloat512

    # Load camera geometry and sort into spiral order
    global nchan
    ncam = calin.iact_data.nectarcam_layout.nectarcam_layout()
    scam = calin.iact_data.instrument_layout.reorder_camera_channels(ncam, ncam.pixel_spiral_channel_index())
    nchan = scam.channel_size()

    # Load pulse shape and determine number of samples and sample time
    global isample0
    global nsample
    pulse = calin.simulation.vs_cta.mstn_impulse_response()
    hg = pulse['hg']
    dtsample = pulse['dt']
    isample0 = len(hg)
    nsample = 1 << (len(hg)-1).bit_length() + 1  # Next power of two greater than 2*len(hg)

    # Generate electronics simulation
    global electronics_sim
    electronics_sim = electronics_sim_class(1,nchan,nsample,dtsample,isample0)

    # Register pulse shape and configure electronics simulation
    electronics_sim.register_impulse_response(hg, 'DC')
    
    # Add camera response
    electronics_sim.add_camera_response(
        numpy.asarray([0]), args.coupling == 'ac')
    
    # Configure the camera neighbbor map
    electronics_sim.set_cr_neighbors(0, scam)
    
    # Load PE spectrum for prompt photons
    global pe_gen
    pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, quiet=True)

    # Define NSB PE generator and rate and register them
    if args.nsb>0:
        pe_gen_nsb = None
        if args.no_after_pulsing:
            pe_gen_nsb = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, quiet=True)
        else:
            pe_gen_nsb = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, afterpulsing=True, quiet=True)
        pe_gen_nsb.this.disown() # Let electronics_sim own it

        nsb = numpy.zeros(scam.channel_size()) + args.nsb
        electronics_sim.set_cr_nsb_rate(0, nsb, pe_gen_nsb, True)

    # Define noise spectrum
    if(args.noise):
        def noise(f):
            fhi = 330.0
            flo = 0.0
            return (1+0.2*(f/275.0)**2)*(numpy.tanh((fhi-f)/20.0)+1)*0.9 + 0.6  # *(numpy.tanh((f-flo)/1.0)+1)/16
        freq = electronics_sim.spectral_frequencies_ghz(False)
        noise_spectrum = noise(freq*1000)
        noise_spectrum[0] = 0
        noise_spectrum *= numpy.sqrt(10.5**2/16/electronics_sim.noise_spectrum_var(noise_spectrum))
        electronics_sim.set_cr_noise_spectrum(0, noise_spectrum)

    # Select trigger algorithm
    global trigger_method
    if args.trigger == 'multiplicity':
        trigger_method = 'trigger_multiplicity_cr'
        electronics_sim.set_cr_multiplicity(0, args.multiplicity)
    elif args.trigger == 'm3':
        trigger_method = 'trigger_multiplicity_cr'
        electronics_sim.set_cr_multiplicity(0, 3)
    elif args.trigger == 'm4':
        trigger_method = 'trigger_multiplicity_cr'
        electronics_sim.set_cr_multiplicity(0, 4)
    elif args.trigger == '3nn':
        trigger_method = 'trigger_3nn_cr'
    elif args.trigger == '4nn':
        trigger_method = 'trigger_4nn_cr'
    else:
        raise ValueError(f'Unknown trigger algorithm: {args.trigger}')

    # Set the threshold and coincidence time for the trigger
    electronics_sim.set_cr_threshold(0, numpy.zeros(nchan) + args.threshold, args.coincidence)

    # Transit time spread
    global tts
    tts = args.tts

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
    fraction = num_files_processed/num_files
    line = f'{filename}: {num_files_processed:,d} / {num_files:,d} ({fraction*100:.1f}%)'
    runtime = (datetime.datetime.now(datetime.timezone.utc) - begin_utc).total_seconds()
    totaltime = runtime/fraction
    line += f' ; {format_duration(runtime)} / {format_duration(totaltime)}'
    line += f' ; {num_events_triggered:,d} / {num_events_processed:,d} ({num_events_processed/runtime:,.1f} Hz)'
    print(line)
    
def scan_file(file):
    results = []
    sim_event = calin.ix.simulation.simulated_event.SimulatedEvent()
    reader = type(sim_event).NewHDFStreamReader(file, 'events')
    nevent = reader.nrow()
    nprocessed = 0
    trigger_fn = getattr(electronics_sim, trigger_method)
    for ievent in range(nevent):
        reader.read(ievent, sim_event)
        narray = sim_event.array_size()
        nprocessed += narray
        for iarray in range(narray):
            array = sim_event.array(iarray)
            group = array.group(0)
            if group.detector_size()==0:
                continue
            detector = group.detector(0)
            if detector.pixel_size()<3:
                continue
            electronics_sim.load_from_simulated_event_with_pmt_noise(group, pe_gen, tts)
            electronics_sim.clear_waveforms()
            electronics_sim.transfer_scope_pes_to_waveform(0)
            if args.nsb>0:
                electronics_sim.add_nsb_noise_to_waveform_cr(0)
            electronics_sim.convolve_impulse_response_fftw_codelet_cr(0)

            itrig = trigger_fn(0, isample0)
            if itrig >= 0:
                npe = sorted([x.integer_time_size() for x in detector.pixel()])
                results.append([ievent, iarray,                                  # 0, 1
                                sim_event.energy(),                              # 2
                                sim_event.u0().x(), sim_event.u0().y(), sim_event.u0().z(), # 3, 4, 5
                                sim_event.viewcone_costheta(),                   # 6
                                array.scattered_offset().x(),array.scattered_offset().y(), # 7, 8
                                array.scattered_distance(),                      # 9
                                npe[-3], npe[-2], npe[-1],                       # 10, 11, 12
                                detector.pixel_size(), sum(npe), itrig,          # 13, 14, 15
                                detector.reference_time(), detector.time_max()]) # 16, 17
    return file, nprocessed, results

def save_results(file, nevent, results, filehandle):
    output = dict(
        file = file,
        nevent = nevent,
        results = results)
    pickle.dump(output, filehandle)
    filehandle.flush()

def process_results(results, filehandle):
    global num_files_processed
    global num_events_processed
    global num_events_triggered

    file, nevent, results = results
    num_events_processed += nevent
    num_events_triggered += len(results)
    num_files_processed += 1
    print_line(file)

    save_results(file, nevent, results, filehandle)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Nectarcam single-telescope shower simulation')
    
    # Simulation parameters
    parser.add_argument('-o', '--output', type=str, default='triggered_events.pickle',
        help='Output filename.')
    parser.add_argument('-g', '--glob', type=str, default='*.h5',
        help='Input file glob pattern.')

    # Noise and trigger parameters
    parser.add_argument('--nsb', type=float, default=0.30,
                    help='Specify the NSB rate in GHz')
    parser.add_argument('--noise', action='store_true', help='Add electronics noise')
    parser.add_argument('--coupling', choices=['ac', 'dc'], default='ac',
                        help='Camera coupling mode (default: ac)')
    parser.add_argument('--no_after_pulsing', action='store_true', help='Disable after-pulsing in SPE spectrum (default: enabled)')
    parser.add_argument('-t', '--trigger', type=str, default='3nn', choices=['3nn','4nn','m3','m4','multiplicity'],
                        help='Trigger algorithm to use (default: 3nn)')
    parser.add_argument('-m', '--multiplicity', type=int, default=3,
                        help='Channel multiplicity if "multiplicity" algorithm is selected')
    parser.add_argument('-c', '--coincidence', type=int, default=24,
                        help='Set the L1 trigger coincidence time in samples')
    parser.add_argument('--threshold', type=float, default=140.0,
                        help='Threshold for trigger (default: 140.0)')

    # PE spectrum and transit time spread
    parser.add_argument('--tts', type=float, default=0.00,
        help='Specify the transit time spread RMS in ns.')

    # Event processing options
    parser.add_argument('--nthread', type=int, default=0,
        help='Number of threads to use (default: 0 = number of CPUs available)')
    parser.add_argument('--avx', type=int, default=512, choices=[128,256,512],
        help='Set the AVX vector size in bits (default: 512)')

    global args
    args = parser.parse_args()

    files = sorted(glob.glob(args.glob))
    if not files:
        print(f'No input files found matching glob pattern: {args.glob}')
        exit(1)
    print(f'Found {len(files)} input files matching glob pattern: {args.glob}')

    global num_files
    global num_files_processed
    global num_events_processed
    global num_events_triggered
    num_files = len(files)
    num_files_processed = 0
    num_events_processed = 0
    num_events_triggered = 0

    global begin_utc
    begin_utc = datetime.datetime.now(datetime.timezone.utc)

    with open(args.output, 'wb') as f:
        max_workers = args.nthread or os.cpu_count() or 1
    
        if max_workers == 1:
            # Run the simulations in this thread
            init(args)
            while len(files) > 0:
                file = files.pop(0)
                results = scan_file(file)
                process_results(results, filehandle=f)
        else:
            # Use a process pool for parallelism
            batch_size = 4 * max_workers
            print(f'Launching {max_workers} workers in process pool')
            with concurrent.futures.ProcessPoolExecutor(initializer=init, initargs=(args,), max_workers=max_workers) as executor:
                futures = set()

                # Submit initial batch
                first_batch = min(batch_size, len(files))
                for i in range(first_batch):
                    file = files.pop(0)
                    futures.add(executor.submit(scan_file, file))

                # Keep a bounded number of in-flight futures; as one completes, submit another
                while futures:
                    done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                    for fut in done:
                        futures.remove(fut)
                        process_results(fut.result(), filehandle=f)
                        if len(files) > 0:
                            file = files.pop(0)
                            futures.add(executor.submit(scan_file, file))
