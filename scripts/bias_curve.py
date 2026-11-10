#!/usr/bin/env python3

# calin/scripts/bias_curve.py -- Stephen Fegan - 2025-11-24
#
# Caclculate trigger delta-T times for given threhold / NSB rate / algorithm
#
# Copyright 2025, Stephen Fegan <sfegan@llr.in2p3.fr>
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
import json
import datetime
import platform
import numpy
import calin.simulation.detector_efficiency
import calin.simulation.vs_cta
import calin.simulation.ray_processor
import calin.iact_data.instrument_layout
import calin.iact_data.nectarcam_layout
import calin.simulation.vs_panoseti

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Bias curve generator')
parser.add_argument('-n', type=int, default=1000,
                   help='Specify the number of triggers to simulate')
parser.add_argument('--nsb', type=float, default=0.30,
                   help='Specify the NSB rate in GHz')
parser.add_argument('--star_rate', type=float, nargs='*', default=[],
                   help='Specify per-channel additional PE rates from star (one float per channel)')
parser.add_argument('-t', '--threshold', type=float, default=120.0,
                   help='Specify the threshold in DC')
parser.add_argument('--noise', action='store_true', help='Add electronics noise')
parser.add_argument('--coupling', choices=['ac', 'dc'], default='ac',
                    help='Camera coupling mode (default: ac)')
parser.add_argument('--avx', type=int, default=512, choices=[128, 256, 512],
                    help='SIMD vector width in bits (default: 512)')
parser.add_argument('--no_after_pulsing', action='store_true', help='Disable after-pulsing in SPE spectrum (default: enabled)')
parser.add_argument('--crosstalk-mean', type=float, default=0.05,
                    help='PANOSETI SiPM mean crosstalk per triggered cell (default: 0.05)')
parser.add_argument('--spe-resolution', type=float, default=0.1,
                    help='PANOSETI SiPM single-cell charge resolution (sigma/mean; default: 0.1)')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='Write trigger times to this file')
parser.add_argument('--camera', type=str, default='nectarcam', choices=['nectarcam','panoseti','quabo'],
                    help='Camera configuration to simulate (default: nectarcam; quabo simulates one 16x16 PANOSETI Quabo)')
parser.add_argument('-a', '--algorithm', type=str, default='3nn', choices=['3nn','4nn','m1','m2','m3','m4','multiplicity'],
                    help='Trigger algorithm to use (default: 3nn; PANOSETI and quabo support m1, m2, m3, m4, or multiplicity)')
parser.add_argument('-m', '--multiplicity', type=int, default=3,
                    help='Channel multiplicity if "multiplicity" algorithm is selected')
parser.add_argument('-c', '--coincidence', type=int, default=24,
                    help='Set the trigger coincidence time in samples')
args = parser.parse_args()
if args.camera in ('panoseti', 'quabo') and args.algorithm not in ('multiplicity', 'm1', 'm2', 'm3', 'm4'):
    parser.error('PANOSETI and quabo support only multiplicity, m1, m2, m3, or m4 trigger algorithms')

nsb_rate = args.nsb
threshold = args.threshold
isample0 = 480
tcoincidence = args.coincidence

def load_impulse_response():
    if args.camera in ('panoseti', 'quabo'):
        return calin.simulation.vs_panoseti.panoseti_fast_impulse_response()
    return calin.simulation.vs_cta.mstn_impulse_response()

def init(pulse):
    global pe_list_processor, trigger_method, trigger_threshold, trigger_multiplicity, isample0

    hg = pulse['hg']
    dtsample = pulse['dt']
    isample0 = len(hg)

    if args.camera == 'panoseti':
        nchan = 1024
        scam = None
    elif args.camera == 'quabo':
        nchan = 256
        scam = None
    else:
        ncam = calin.iact_data.nectarcam_layout.nectarcam_layout()
        scam = calin.iact_data.instrument_layout.reorder_camera_channels(
            ncam, ncam.pixel_spiral_channel_index())
        nchan = scam.channel_size()

    # Instantiate PE list processor
    pe_processor_classes = {
        128: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat128,
        256: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat256,
        512: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat512,
    }
    pe_list_processor = pe_processor_classes[args.avx](
        1, nchan, 1024, dtsample, isample0)
    pe_list_processor.register_impulse_response(hg, 'DC')

    # Register camera response
    pe_list_processor.add_camera_response(
        numpy.asarray([0]), args.coupling == 'ac')

    # Configure the NectarCam neighbors matrix when simulating NectarCam.
    if scam is not None:
        pe_list_processor.set_cr_neighbors(0, scam)

    # Select trigger algorithm
    trigger_threshold = numpy.zeros(nchan) + threshold
    if args.algorithm == 'multiplicity':
        trigger_method = 'trigger_multiplicity_cr'
        pe_list_processor.set_cr_multiplicity(0, args.multiplicity)
    elif args.algorithm == 'm1':
        trigger_method = 'trigger_multiplicity_cr'
        pe_list_processor.set_cr_multiplicity(0, 1)
    elif args.algorithm == 'm2':
        trigger_method = 'trigger_multiplicity_cr'
        pe_list_processor.set_cr_multiplicity(0, 2)
    elif args.algorithm == 'm3':
        trigger_method = 'trigger_multiplicity_cr'
        pe_list_processor.set_cr_multiplicity(0, 3)
    elif args.algorithm == 'm4':
        trigger_method = 'trigger_multiplicity_cr'
        pe_list_processor.set_cr_multiplicity(0, 4)
    elif args.algorithm == '3nn':
        trigger_method = 'trigger_3nn_cr'
    elif args.algorithm == '4nn':
        trigger_method = 'trigger_4nn_cr'
    else:
        raise ValueError(f'Unknown trigger algorithm: {args.algorithm}')

    if args.camera == 'panoseti':
        if trigger_method != 'trigger_multiplicity_cr':
            raise ValueError(f'PANOSETI supports only multiplicity trigger algorithm')
        trigger_method = 'trigger_panoseti_multiplicity_cr'

    # Instantiate PE generator
    pe_gen = None
    if args.camera in ('panoseti', 'quabo'):
        pe_gen = calin.simulation.vs_panoseti.panoseti_pe_amplitude_generator(
            crosstalk_mean=args.crosstalk_mean,
            spe_resolution=args.spe_resolution,
            avx=args.avx)
    else:
        if args.no_after_pulsing:
            pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, quiet=True)
        else:
            pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(avx=args.avx, afterpulsing=True, quiet=True)
    # The waveform processor uses an architecture-specific SIMD PE generator.
    if pe_gen is not None:
        pe_gen.this.disown() # Let pe_list_processor own it

    # Define NSB
    nsb = numpy.zeros(nchan) + nsb_rate
    # Add per-channel star rates
    for i, star_pe_rate in enumerate(args.star_rate):
        if i < len(nsb):
            nsb[i] += star_pe_rate
    pe_list_processor.set_cr_nsb_rate(0, nsb, pe_gen, pe_gen is not None)

    # Define noise spectrum
    if(args.noise):
        def noise(f):
            fhi = 330.0
            flo = 0.0
            return (1+0.2*(f/275.0)**2)*(numpy.tanh((fhi-f)/20.0)+1)*0.9 + 0.6  # *(numpy.tanh((f-flo)/1.0)+1)/16
        freq = pe_list_processor.spectral_frequencies_ghz(False)
        noise_spectrum = noise(freq*1000)
        noise_spectrum[0] = 0
        noise_spectrum *= numpy.sqrt(10.5**2/16/pe_list_processor.noise_spectrum_var(noise_spectrum))
        pe_list_processor.set_cr_noise_spectrum(0, noise_spectrum)

    # Define threshold
    pe_list_processor.set_cr_threshold(0, trigger_threshold, tcoincidence)

# Load the pulse once in the parent process. Worker initializers receive the
# loaded response so they can set up their own processors without reloading it.
pulse_response = load_impulse_response()
sampling_interval_ns = pulse_response['dt']
sampling_rate_ghz = 1.0 / sampling_interval_ns

# Prepare a JSON header describing the run
config = vars(args).copy()
config['sampling_interval_ns'] = sampling_interval_ns
config['sampling_rate_ghz'] = sampling_rate_ghz
config['_generated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
config['_host'] = platform.node()
header = 'Configuration (JSON):\n' + json.dumps(config, indent=2)

def trigger():
    trigger_fn = getattr(pe_list_processor, trigger_method)
    return trigger_fn(0, isample0)

def one_trigger():
    pe_list_processor.clear_waveforms()
    pe_list_processor.add_nsb_noise_to_waveform_cr(0)
    pe_list_processor.convolve_impulse_response_fftw_codelet_cr(0)
    ttrig = 0
    itrig = trigger()
    while itrig == -1:
        ttrig += 1024-isample0
        pe_list_processor.shift_and_clear_waveforms(isample0)
        pe_list_processor.add_nsb_noise_to_waveform_cr(0, isample0)
        pe_list_processor.convolve_impulse_response_fftw_codelet_cr(0)
        itrig = trigger()
    ttrig += itrig-isample0
    # vwf = pe_list_processor.v_waveform();
    # print(numpy.mean(vwf,axis=0)[:5],numpy.var(vwf,axis=0)[:5])
    return ttrig

trigger_times = []
max_workers = os.cpu_count() or 4
batch_size = max_workers * 10
with concurrent.futures.ProcessPoolExecutor(initializer=init,
        initargs=(pulse_response,), max_workers=max_workers) as executor:
    remaining = args.n
    futures = set()

    # Submit initial batch
    first_batch = min(batch_size, remaining)
    for _ in range(first_batch):
        futures.add(executor.submit(one_trigger))
    remaining -= first_batch

    # Keep a bounded number of in-flight futures; as one completes, submit another
    while futures:
        done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
        for fut in done:
            futures.remove(fut)
            trigger_times.append(fut.result())
            elapsed_ms = numpy.sum(trigger_times) * sampling_interval_ns / 1e6
            rate_khz = 1e6 / (numpy.mean(trigger_times) * sampling_interval_ns)
            print(f'{len(trigger_times)} : {elapsed_ms:.2f} ms {rate_khz:.3f} kHz')

            if args.output and len(trigger_times)>0 and (len(trigger_times)%100)==0:
                numpy.savetxt(args.output, trigger_times, header=header, fmt='%d')
                print(f'Wrote {len(trigger_times)} trigger times to {args.output}')

            if remaining > 0:
                futures.add(executor.submit(one_trigger))
                remaining -= 1

if args.output and (len(trigger_times)==0 or (len(trigger_times)%100)!=0):
    numpy.savetxt(args.output, trigger_times, header=header, fmt='%d')
    print(f'Wrote {len(trigger_times)} trigger times to {args.output}')
