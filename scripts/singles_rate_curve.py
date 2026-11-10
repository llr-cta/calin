#!/usr/bin/env python3

# calin/scripts/singles_curve.py -- Stephen Fegan - 2025-11-24
#
# Caclculate trogger delta-T times for given threhold / NSB rate.
#
# Note this accumualtes stats across multiple channele in parallel
# to take advantage of vectorization. The algorithm is therefore
# somewhat complicated by the need to track trigger times per channel.
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

import argparse
import concurrent.futures
import json
import datetime
import platform
import numpy
import calin.simulation.detector_efficiency
import calin.simulation.vs_cta
import calin.simulation.vs_panoseti
import calin.simulation.ray_processor
import calin.iact_data.instrument_layout
import calin.iact_data.nectarcam_layout

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Bias curve generator')
parser.add_argument('-n', type=int, default=1000,
                   help='Specify the number of triggers to simulate per channel')
parser.add_argument('--nthread', type=int, default=1,
                   help='Number of worker processes (1 runs without a process pool)')
parser.add_argument('--nchannels', type=int, default=16,
                    help='Number of independent channels to simulate (default: 16; use a SIMD-sized batch for PANOSETI)')
parser.add_argument('--nsb', type=float, default=0.30,
                   help='Specify the NSB rate in GHz')
parser.add_argument('-t', '--threshold', type=float, default=120.0,
                   help='Specify the threshold in DC')
parser.add_argument('--noise', action='store_true', help='Add electronics noise')
parser.add_argument('--coupling', choices=['ac', 'dc'], default='ac',
                    help='Camera coupling mode (default: ac)')
parser.add_argument('--camera', type=str, default='nectarcam', choices=['nectarcam', 'panoseti'],
                    help='Camera configuration to simulate (default: nectarcam)')
parser.add_argument('--avx', type=int, default=512, choices=[128, 256, 512],
                    help='SIMD vector width in bits (default: 512)')
parser.add_argument('--crosstalk-mean', type=float, default=0.05,
                    help='PANOSETI SiPM mean crosstalk per triggered cell (default: 0.05)')
parser.add_argument('--spe-resolution', type=float, default=0.1,
                    help='PANOSETI SiPM single-cell charge resolution (sigma/mean; default: 0.1)')
parser.add_argument('-o', '--output', type=str, default=None,
                   help='Write trigger times to this file')
args = parser.parse_args()
if args.n < 0:
    parser.error('-n must be nonnegative')
if args.nthread < 1:
    parser.error('--nthread must be at least 1')

nsb_rate = args.nsb
threshold = args.threshold
isample0 = 480

def load_impulse_response():
    if args.camera == 'panoseti':
        return calin.simulation.vs_panoseti.panoseti_fast_impulse_response()
    return calin.simulation.vs_cta.mstn_impulse_response()

def init(pulse):
    global pe_list_processor, isample0

    hg = pulse['hg']
    dtsample = pulse['dt']
    isample0 = len(hg)

    pe_processor_classes = {
        128: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat128,
        256: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat256,
        512: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat512,
    }
    pe_list_processor = pe_processor_classes[args.avx](
        1, args.nchannels, 1024, dtsample, isample0)
    pe_list_processor.register_impulse_response(hg, 'DC')

    # Register camera response
    pe_list_processor.add_camera_response(
        numpy.asarray([0]), args.coupling == 'ac')

    # Instantiate PE generator
    if args.camera == 'panoseti':
        pe_gen = calin.simulation.vs_panoseti.panoseti_pe_amplitude_generator(
            crosstalk_mean=args.crosstalk_mean,
            spe_resolution=args.spe_resolution,
            avx=args.avx)
    else:
        pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(
            avx=args.avx, afterpulsing=True, quiet=True)
    pe_gen.this.disown() # Let pe_list_processor own it

    # Define NSB
    nsb = numpy.zeros(args.nchannels) + nsb_rate
    pe_list_processor.set_cr_nsb_rate(0, nsb, pe_gen, True)

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
    pe_list_processor.set_cr_threshold(0, numpy.zeros(args.nchannels) + threshold, 24)

def gen_triggers(ntrigger):
    all_ttrig = []
    ntrig = numpy.zeros(args.nchannels, dtype=int)
    ttrig = numpy.zeros(args.nchannels, dtype=int)

    pe_list_processor.clear_waveforms()
    pe_list_processor.add_nsb_noise_to_waveform_cr(0)
    pe_list_processor.convolve_impulse_response_fftw_codelet_cr(0)

    while(numpy.count_nonzero(ntrig<ntrigger) > 0):
        _, itrig = pe_list_processor.trigger_singles_cr(0, isample0)
    
        trigmask = (ttrig>=0) * (itrig>=isample0)
        ttrig[(ttrig>=0) * (itrig<isample0)] += 1024-isample0
        ttrig[ttrig<0] = 0
        ttrig[ntrig>=ntrigger] = -1
        for ichan in numpy.where(trigmask)[0]:
            all_ttrig.append(ttrig[ichan] + itrig[ichan] - isample0)
            ttrig[ichan] = -1
            ntrig[ichan] += 1

        pe_list_processor.shift_and_clear_waveforms(isample0)
        pe_list_processor.add_nsb_noise_to_waveform_cr(0, isample0)
        pe_list_processor.convolve_impulse_response_fftw_codelet_cr(0)
    return all_ttrig

impulse_response = load_impulse_response()
sampling_interval_ns = impulse_response['dt']
sampling_rate_ghz = 1.0 / sampling_interval_ns

# Prepare a JSON header describing the run
config = vars(args).copy()
config['sampling_interval_ns'] = sampling_interval_ns
config['sampling_rate_ghz'] = sampling_rate_ghz
config['_generated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
config['_host'] = platform.node()
header = 'Configuration (JSON):\n' + json.dumps(config, indent=2)

if args.nthread == 1 or args.n == 0:
    init(impulse_response)
    trigger_times = gen_triggers(args.n)
else:
    nworker = min(args.nthread, args.n)
    events_per_worker, extra_events = divmod(args.n, nworker)
    worker_event_counts = [events_per_worker + (iworker < extra_events)
        for iworker in range(nworker)]
    trigger_times = []
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=nworker, initializer=init, initargs=(impulse_response,)) as executor:
        for worker_trigger_times in executor.map(gen_triggers, worker_event_counts):
            trigger_times.extend(worker_trigger_times)

elapsed_ms = numpy.sum(trigger_times) * sampling_interval_ns / 1e6
rate_khz = 1e6 / (numpy.mean(trigger_times) * sampling_interval_ns)
print(f'{len(trigger_times)} : {elapsed_ms:.2f} ms {rate_khz:.3f} kHz')

if args.output:
    numpy.savetxt(args.output, trigger_times, header=header, fmt='%d')
    print(f'Wrote {len(trigger_times)} trigger times to {args.output}')
