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
    cam_time_spread_ns = numpy.zeros_like(cam_freq_offset_ppm)
    cam_d2_per_event = numpy.zeros_like(cam_freq_offset_ppm)

    print("summarize_camera_clock_regressions: don't forget to remove temporary code")
    cl = stage1.const_run_config().const_camera_layout()
    if(cl.camera_clock_name_size() != cl.camera_clock_frequency_size()):
        # Temporary : remove when all data has been reprocessed with new version
        cl = calin.iact_data.instrument_layout.camera_layout(cl.camera_type())

    principal_clock_id = stage1.const_clock_regression().principal_clock_id();
    if(principal_clock_id == 0):
        # Temporary : remove when all data has been reprocessed with new version
        principal_clock_id = stage1.const_config().const_clock_regression().principal_clock_id()
    principal_clock_freq = cl.camera_clock_frequency(principal_clock_id)

    for iclock in range(len(cam_freq_offset_ppm)):
        clock = stage1.const_clock_regression().const_camera_clock(iclock)

        clock_id = clock.clock_id()
        if(clock_id == 0):
            # Temporary : remove when all data has been reprocessed with new version
            clock_id = stage1.const_config().const_clock_regression().default_nectarcam_camera_clocks(iclock).clock_id()
        clock_freq = cl.camera_clock_frequency(clock_id)
        clock_nominal_a = clock_freq/principal_clock_freq

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

        mask_n = all_n>5

        # all_x_at_0 = (numpy.remainder(all_x0,1000000000) - all_y0) + all_y0*(1-1/all_a) - all_b/all_a
        all_x_at_0 = numpy.remainder(all_x0[mask_n],1000000000) - all_y0[mask_n]/all_a[mask_n] - all_b[mask_n]/all_a[mask_n]

        cam_freq_offset_ppm[iclock] = 1e6*(numpy.mean(all_a[mask_n])/clock_nominal_a-1) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_freq_spread_ppm[iclock] = 1e6*(numpy.max(all_a[mask_n])-numpy.min(all_a[mask_n]))/clock_nominal_a if numpy.count_nonzero(mask_n) else numpy.nan
        cam_time_offset_ns[iclock] = numpy.mean(all_x_at_0) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_time_spread_ns[iclock] = (numpy.max(all_x_at_0)-numpy.min(all_x_at_0)) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_d2_per_event[iclock] = (1e9/clock_freq)**2*numpy.mean(all_d2[mask_n]/all_n[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan

    return cam_freq_offset_ppm, cam_freq_spread_ppm, cam_time_offset_ns, cam_time_spread_ns, cam_d2_per_event

def summarize_module_clock_regression(stage1, iclock=0):
    mod_freq_offset_ppm = numpy.zeros(stage1.run_config().configured_module_id_size())
    mod_freq_spread_ppm = numpy.zeros_like(mod_freq_offset_ppm)
    mod_time_offset_ns = numpy.zeros_like(mod_freq_offset_ppm)
    mod_time_spread_ns = numpy.zeros_like(mod_freq_offset_ppm)
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

        mask_n = all_n>5

        # all_x_at_0 = (numpy.remainder(all_x0,1000000000) - all_y0) + all_y0*(1-1/all_a) - all_b/all_a
        all_x_at_0 = numpy.remainder(all_x0[mask_n],1000000000) - all_y0[mask_n]/all_a[mask_n] - all_b[mask_n]/all_a[mask_n]

        mask_oob_intercept = numpy.bitwise_or(all_x_at_0>250, all_x_at_0<40)
        mask_oob_residual = numpy.sqrt(all_d2[mask_n]) > 100*numpy.sqrt(all_n[mask_n])
        mask_oob = numpy.bitwise_or(mask_oob_intercept, mask_oob_residual)

        mod_freq_offset_ppm[imod] = 1e6*(numpy.mean(all_a[mask_n])-1) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_freq_spread_ppm[imod] = 1e6*(numpy.max(all_a[mask_n])-numpy.min(all_a[mask_n])) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_time_offset_ns[imod] = numpy.mean(all_x_at_0) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_time_spread_ns[imod] = (numpy.max(all_x_at_0)-numpy.min(all_x_at_0)) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_d2_per_event[imod] = numpy.mean(all_d2[mask_n]/all_n[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan
        mod_problem_bins[imod] = numpy.count_nonzero(mask_oob)

    return mod_freq_offset_ppm, mod_freq_spread_ppm, mod_time_offset_ns, mod_time_spread_ns, mod_d2_per_event, mod_problem_bins

def estimate_run_pedestal(stage1, low_gain=False):
    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()

    charge_stats = stage1.const_charge_stats().const_low_gain() if low_gain \
        else stage1.const_charge_stats().const_high_gain()

    if(charge_stats.ped_trigger_event_count()):
        nsamp = rc.num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_mean()/nsamp
    else:
        nsamp = stage1.config().low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.config().has_low_gain_opt_sum()) \
            else stage1.config().high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_mean()/nsamp

    values[nevent==0] = numpy.nan
    return values
