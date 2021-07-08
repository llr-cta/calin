# calin/python/diagnostics/stage1_analysis.py -- Stephen Fegan -- 2021-06-30
#
# Analysis of stage-one data
#
# Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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
import calin.diagnostics.stage1
import calin.iact_data.instrument_layout

def summarize_camera_clock_regressions(stage1):
    cam_freq_offset_ppm = numpy.zeros(stage1.const_clock_regression().camera_clock_size())
    cam_freq_spread_ppm = numpy.zeros_like(cam_freq_offset_ppm)
    cam_time_offset_ns = numpy.zeros_like(cam_freq_offset_ppm)
    cam_d2_per_event = numpy.zeros_like(cam_freq_offset_ppm)

    print("summarize_camera_clock_regressions: don't forget to remove temporary code")
    cl = stage1.const_run_config().const_camera_layout()
    if(cl.camera_clock_name_size() != cl.camera_clock_frequency_size()):
        # Temporary : remove when all data has been reprocessed with new version
        cl = calin.iact_data.instrument_layout.camera_layout(cl.camera_type())

    master_clock_id = stage1.const_clock_regression().master_clock_id();
    if(master_clock_id == 0):
        # Temporary : remove when all data has been reprocessed with new version
        master_clock_id = stage1.const_config().const_clock_regression().master_clock_id()
    master_clock_freq = cl.camera_clock_frequency(master_clock_id)

    for iclock in range(len(cam_freq_offset_ppm)):
        clock = stage1.const_clock_regression().const_camera_clock(iclock)

        clock_id = clock.clock_id()
        if(clock_id == 0):
            # Temporary : remove when all data has been reprocessed with new version
            clock_id = stage1.const_config().const_clock_regression().default_nectarcam_camera_clocks(iclock).clock_id()
        clock_freq = cl.camera_clock_frequency(clock_id)
        clock_nominal_a = clock_freq/master_clock_freq
        print(iclock, clock_id, clock_freq, master_clock_freq, clock_nominal_a)

        all_t = sorted(clock.bins_keys())
        all_x0 = numpy.zeros(len(all_t), dtype=numpy.int64)
        all_y0 = numpy.zeros(len(all_t), dtype=numpy.int64)
        all_a = numpy.zeros(len(all_t))
        all_b = numpy.zeros(len(all_t))
        all_d2 = numpy.zeros(len(all_t))
        all_n = numpy.zeros(len(all_t))

        for ikey, key in enumerate(all_t):
            bin = clock.const_bins(key)
            all_x0[ikey] = bin.x0()
            all_y0[ikey] = bin.y0()
            all_a[ikey] = bin.a()
            all_b[ikey] = bin.b()
            all_d2[ikey] = bin.d2()
            all_n[ikey] = bin.num_entries()

        # all_y_at_0 = (numpy.remainder(all_x0,1000000000) - all_y0) + all_y0*(all_a-1) - all_b/all_a
        all_y_at_0 = numpy.remainder(all_x0,1000000000) - all_y0/all_a - all_b/all_a

        mask_n = all_n>5

        cam_freq_offset_ppm[iclock] = 1e6*(numpy.mean(all_a[mask_n])/clock_nominal_a-1) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_freq_spread_ppm[iclock] = 1e6*(numpy.max(all_a[mask_n])-numpy.min(all_a[mask_n]))/clock_nominal_a if numpy.count_nonzero(mask_n) else numpy.nan
        cam_time_offset_ns[iclock] = 1e9/clock_freq*numpy.mean(all_y_at_0[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_d2_per_event[iclock] = (1e9/clock_freq)**2*numpy.mean(all_d2[mask_n]/all_n[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan

    return cam_freq_offset_ppm, cam_freq_spread_ppm, cam_time_offset_ns, cam_d2_per_event

def summarize_module_clock_regression(stage1, iclock=0):
    mod_freq_offset_ppm = numpy.zeros(stage1.run_config().configured_module_id_size())
    mod_freq_spread_ppm = numpy.zeros_like(mod_freq_offset_ppm)
    mod_time_offset_ns = numpy.zeros_like(mod_freq_offset_ppm)
    mod_d2_per_event = numpy.zeros_like(mod_freq_offset_ppm)
    mod_problem_bins = numpy.zeros_like(mod_freq_offset_ppm)

    for imod in range(stage1.run_config().configured_module_id_size()):
        clock = stage1.const_clock_regression().const_module_clock(iclock).const_modules(imod)

        all_t = sorted(clock.bins_keys())
        all_x0 = numpy.zeros(len(all_t), dtype=numpy.int64)
        all_y0 = numpy.zeros(len(all_t), dtype=numpy.int64)
        all_a = numpy.zeros(len(all_t))
        all_b = numpy.zeros(len(all_t))
        all_d2 = numpy.zeros(len(all_t))
        all_n = numpy.zeros(len(all_t))

        for ikey, key in enumerate(all_t):
            bin = clock.const_bins(key)
            all_x0[ikey] = bin.x0()
            all_y0[ikey] = bin.y0()
            all_a[ikey] = bin.a()
            all_b[ikey] = bin.b()
            all_d2[ikey] = bin.d2()
            all_n[ikey] = bin.num_entries()

        # all_y_at_0 = (numpy.remainder(all_x0,1000000000) - all_y0) + all_y0*(all_a-1) - all_b/all_a
        all_y_at_0 = numpy.remainder(all_x0,1000000000) - all_y0/all_a - all_b/all_a

        mask_oob_intercept = numpy.bitwise_or(all_y_at_0>250, all_y_at_0<40)
        mask_oob_residual = numpy.sqrt(all_d2) > 100*numpy.sqrt(all_n)
        mask_n = all_n>5

        mask = numpy.bitwise_and(numpy.bitwise_or(mask_oob_intercept, mask_oob_residual), mask_n)

        mod_freq_offset_ppm[imod] = 1e6*(numpy.mean(all_a[mask_n])-1) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_freq_spread_ppm[imod] = 1e6*(numpy.max(all_a[mask_n])-numpy.min(all_a[mask_n])) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_time_offset_ns[imod] = numpy.mean(all_y_at_0[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_d2_per_event[imod] = numpy.mean(all_d2[mask_n]/all_n[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_problem_bins[imod] = numpy.count_nonzero(mask)

    return mod_freq_offset_ppm, mod_freq_spread_ppm, mod_time_offset_ns, mod_d2_per_event, mod_problem_bins
