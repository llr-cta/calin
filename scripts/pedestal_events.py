#!/usr/bin/env python3

# calin/scripts/pedestal_events.py
#
# Generate pedestal samples from NSB-only waveform simulations.
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 or later.

import argparse
import concurrent.futures
import datetime
import json
import platform
import numpy
import calin.simulation.vs_cta
import calin.simulation.vs_panoseti
import calin.simulation.ray_processor


parser = argparse.ArgumentParser(description='Generate pedestal events')
parser.add_argument('-n', type=int, default=1000,
                    help='Number of events to simulate')
parser.add_argument('--nthread', type=int, default=1,
                    help='Number of worker processes (1 runs without a process pool)')
parser.add_argument('--nchannels', type=int, default=16,
                    help='Number of independent channels to simulate')
parser.add_argument('--nsb', type=float, default=0.30,
                    help='NSB rate in GHz')
parser.add_argument('--noise', action='store_true',
                    help='Add electronics noise')
parser.add_argument('--coupling', choices=['ac', 'dc'], default='ac',
                    help='Camera coupling mode (default: ac)')
parser.add_argument('--camera', choices=['nectarcam', 'panoseti'], default='nectarcam',
                    help='Camera response to simulate (default: nectarcam)')
parser.add_argument('--avx', type=int, choices=[128, 256, 512], default=512,
                    help='SIMD vector width in bits (default: 512)')
parser.add_argument('--crosstalk-mean', type=float, default=0.05,
                    help='PANOSETI SiPM mean crosstalk per triggered cell (default: 0.05)')
parser.add_argument('--spe-resolution', type=float, default=0.1,
                    help='PANOSETI single-cell charge resolution, sigma/mean (default: 0.1)')
parser.add_argument('--sample', type=int, default=512,
                    help='First waveform sample in the pedestal sum (default: 512)')
parser.add_argument('--nsamples', type=int, default=1,
                    help='Number of samples to sum for each pedestal (default: 1)')
parser.add_argument('--stride', type=int, default=1,
                    help='Stride between samples in the pedestal sum (default: 1)')
parser.add_argument('-o', '--output', type=str, default='pedestal_events.dat',
                    help='Output table path (default: pedestal_events.dat)')
args = parser.parse_args()

if args.n < 0:
    parser.error('-n must be nonnegative')
if args.nthread < 1:
    parser.error('--nthread must be at least 1')
if args.nchannels < 1:
    parser.error('--nchannels must be at least 1')
if args.nsamples < 1:
    parser.error('--nsamples must be at least 1')
if args.stride < 1:
    parser.error('--stride must be at least 1')
if args.sample < 0:
    parser.error('--sample must be nonnegative')

sample_indices = args.sample + numpy.arange(args.nsamples) * args.stride
if sample_indices[-1] >= 1024:
    parser.error('Pedestal sample indices must be less than 1024')


def load_impulse_response():
    if args.camera == 'panoseti':
        return calin.simulation.vs_panoseti.panoseti_slow_impulse_response()
    return calin.simulation.vs_cta.mstn_impulse_response()


def init(pulse):
    global pe_list_processor, isample_advance

    impulse_response = pulse['hg']
    sample_interval_ns = pulse['dt']
    isample_advance = len(impulse_response)

    pe_processor_classes = {
        128: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat128,
        256: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat256,
        512: calin.simulation.ray_processor.VCLWaveformPEProcessorFloat512,
    }
    pe_list_processor = pe_processor_classes[args.avx](
        1, args.nchannels, 1024, sample_interval_ns, isample_advance)
    pe_list_processor.register_impulse_response(impulse_response, 'DC')
    pe_list_processor.add_camera_response(
        numpy.asarray([0]), args.coupling == 'ac')

    if args.camera == 'panoseti':
        pe_gen = calin.simulation.vs_panoseti.panoseti_pe_amplitude_generator(
            crosstalk_mean=args.crosstalk_mean,
            spe_resolution=args.spe_resolution,
            avx=args.avx)
    else:
        pe_gen = calin.simulation.vs_cta.vcl_mstn_spe_amplitude_generator(
            avx=args.avx, afterpulsing=True, quiet=True)
    pe_gen.this.disown()

    nsb = numpy.full(args.nchannels, args.nsb)
    pe_list_processor.set_cr_nsb_rate(0, nsb, pe_gen, True)

    if args.noise:
        def noise(f):
            fhi = 330.0
            return (1 + 0.2 * (f / 275.0)**2) * (numpy.tanh((fhi - f) / 20.0) + 1) * 0.9 + 0.6
        freq = pe_list_processor.spectral_frequencies_ghz(False)
        noise_spectrum = noise(freq * 1000)
        noise_spectrum[0] = 0
        noise_spectrum *= numpy.sqrt(
            10.5**2 / 16 / pe_list_processor.noise_spectrum_var(noise_spectrum))
        pe_list_processor.set_cr_noise_spectrum(0, noise_spectrum)


def generate_pedestal_events(nevent):
    pedestal_events = numpy.empty((nevent, args.nchannels), dtype=float)
    for ievent in range(nevent):
        pe_list_processor.clear_waveforms()
        pe_list_processor.add_nsb_noise_to_waveform_cr(0)
        pe_list_processor.convolve_impulse_response_fftw_codelet_cr(0)
        waveform = numpy.asarray(pe_list_processor.v_waveform())
        pedestal_events[ievent, :] = waveform[sample_indices, :].sum(axis=0)
    return pedestal_events


def main():
    pulse = load_impulse_response()
    sampling_interval_ns = pulse['dt']
    config = vars(args).copy()
    config['sampling_interval_ns'] = sampling_interval_ns
    config['sampling_rate_ghz'] = 1.0 / sampling_interval_ns
    config['pedestal_sample_indices'] = sample_indices.tolist()
    config['data_shape'] = [args.n, args.nchannels]
    config['_generated_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    config['_host'] = platform.node()
    header = 'Configuration (JSON):\n' + json.dumps(config, indent=2)

    if args.nthread == 1 or args.n == 0:
        init(pulse)
        pedestal_events = generate_pedestal_events(args.n)
    else:
        nworker = min(args.nthread, args.n)
        events_per_worker, extra_events = divmod(args.n, nworker)
        worker_event_counts = [events_per_worker + (iworker < extra_events)
                               for iworker in range(nworker)]
        event_blocks = []
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=nworker, initializer=init, initargs=(pulse,)) as executor:
            for worker_events in executor.map(generate_pedestal_events, worker_event_counts):
                event_blocks.append(worker_events)
        pedestal_events = numpy.concatenate(event_blocks, axis=0)

    numpy.savetxt(args.output, pedestal_events, header=header, fmt='%.9g')
    print(f'Wrote {args.n} events x {args.nchannels} channels to {args.output}')


if __name__ == '__main__':
    main()
