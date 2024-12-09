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
import scipy.special
import scipy.stats
import scipy.optimize
import collections.abc

import calin.diagnostics.stage1
import calin.iact_data.instrument_layout
import calin.ix.math.histogram

def get_camera_clock_regression_data(stage1, iclock=0):
    cl = stage1.const_run_config().const_camera_layout()
    clock = stage1.const_clock_regression().const_camera_clock(iclock)

    clock_id = clock.clock_id()
    clock_freq = cl.camera_clock_frequency(clock_id)

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

    return clock_id, clock_freq, all_t, all_x0, all_y0, all_a, all_b, all_d2, all_n

def summarize_camera_clock_regressions(stage1):
    cam_freq_offset_ppm = numpy.zeros(stage1.const_clock_regression().camera_clock_size())
    cam_freq_spread_ppm = numpy.zeros_like(cam_freq_offset_ppm)
    cam_time_offset_ns = numpy.zeros_like(cam_freq_offset_ppm)
    cam_time_spread_ns = numpy.zeros_like(cam_freq_offset_ppm)
    cam_d2_per_event = numpy.zeros_like(cam_freq_offset_ppm)

    cl = stage1.const_run_config().const_camera_layout()

    principal_clock_id = stage1.const_clock_regression().principal_clock_id();
    principal_clock_freq = cl.camera_clock_frequency(principal_clock_id)

    for iclock in range(len(cam_freq_offset_ppm)):
        _, clock_freq, _, all_x0, all_y0, all_a, all_b, all_d2, all_n = get_camera_clock_regression_data(stage1, iclock)
        clock_nominal_a = clock_freq/principal_clock_freq
        mask_n = all_n>5

        # all_x_at_0 = (numpy.remainder(all_x0,1000000000) - all_y0) + all_y0*(1-1/all_a) - all_b/all_a
        all_x_at_0 = numpy.remainder(all_x0[mask_n],1000000000) - all_y0[mask_n]/all_a[mask_n] - all_b[mask_n]/all_a[mask_n]

        cam_freq_offset_ppm[iclock] = 1e6*(numpy.mean(all_a[mask_n])/clock_nominal_a-1) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_freq_spread_ppm[iclock] = 1e6*(numpy.max(all_a[mask_n])-numpy.min(all_a[mask_n]))/clock_nominal_a if numpy.count_nonzero(mask_n) else numpy.nan
        cam_time_offset_ns[iclock] = numpy.mean(all_x_at_0) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_time_spread_ns[iclock] = (numpy.max(all_x_at_0)-numpy.min(all_x_at_0)) if numpy.count_nonzero(mask_n) else numpy.nan
        cam_d2_per_event[iclock] = (1e9/clock_freq)**2*numpy.mean(all_d2[mask_n]/all_n[mask_n]) if numpy.count_nonzero(mask_n) else numpy.nan

    return cam_freq_offset_ppm, cam_freq_spread_ppm, cam_time_offset_ns, cam_time_spread_ns, cam_d2_per_event

def get_module_clock_regression_data(stage1, imod, iclock=0):
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

    return all_t, all_x0, all_y0, all_a, all_b, all_d2, all_n

def summarize_module_clock_regression(stage1, iclock=0):
    mod_freq_offset_ppm = numpy.zeros(stage1.run_config().configured_module_id_size())
    mod_freq_spread_ppm = numpy.zeros_like(mod_freq_offset_ppm)
    mod_time_offset_ns = numpy.zeros_like(mod_freq_offset_ppm)
    mod_time_spread_ns = numpy.zeros_like(mod_freq_offset_ppm)
    mod_d2_per_event = numpy.zeros_like(mod_freq_offset_ppm)
    mod_problem_bins = numpy.zeros_like(mod_freq_offset_ppm)

    for imod in range(stage1.run_config().configured_module_id_size()):
        all_t, all_x0, all_y0, all_a, all_b, all_d2, all_n = get_module_clock_regression_data(stage1, imod, iclock)
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
    charge_stats = stage1.const_charge_stats().const_low_gain() if low_gain \
        else stage1.const_charge_stats().const_high_gain()
    if(numpy.max(charge_stats.ped_trigger_event_count()) > 1000):
        nsamp = stage1.const_run_config().num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_mean()/nsamp
    elif True:
        nevent = numpy.zeros_like(charge_stats.all_trigger_event_count())
        values = numpy.zeros_like(charge_stats.all_trigger_event_count(), dtype=float)
        all_mwf = []
        if stage1.has_mean_wf_pedestal(): all_mwf.append(stage1.const_mean_wf_pedestal())
        if stage1.has_mean_wf_physics(): all_mwf.append(stage1.const_mean_wf_physics())
        if stage1.has_mean_wf_external_flasher(): all_mwf.append(stage1.const_mean_wf_external_flasher())
        if stage1.has_mean_wf_internal_flasher(): all_mwf.append(stage1.const_mean_wf_internal_flasher())
        for mwf in all_mwf:
            if low_gain:
                all_chan = [ mwf.channel_low_gain(ichan) for ichan in range(mwf.channel_low_gain_size()) ]
            else:
                all_chan = [ mwf.channel_high_gain(ichan) for ichan in range(mwf.channel_high_gain_size()) ]
            if len(all_chan) == len(nevent):
                nevent += numpy.asarray([ c.num_entries() for c in all_chan ])
                values += numpy.asarray([ c.num_entries()*(c.mean_waveform(4) if c.num_entries() else 0) for c in all_chan ])
        values[nevent>0] /= nevent[nevent>0]
    else:
        nsamp = stage1.const_config().const_low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.const_config().has_low_gain_opt_sum()) \
            else stage1.const_config().const_high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_mean()/nsamp

    values[nevent==0] = numpy.nan
    return values

def median_measured_voltage(stage1):
    if(not stage1.has_nectarcam() or not stage1.const_nectarcam().has_ancillary_data()):
        return None
    nca_data = stage1.const_nectarcam().const_ancillary_data()
    k = nca_data.hvpa_voltage_keys();
    if(len(k) == 0):
        return None
    v = numpy.zeros(len(k))
    for i,k in enumerate(k):
        all_m = nca_data.const_hvpa_voltage(k)
        vm = []
        for j in range(all_m.measurement_size()):
            vm.append(all_m.const_measurement(j).voltage())
        v[i] = numpy.median(vm)
    return numpy.median(v)

def median_voltage(stage1):
    if(stage1.const_run_config().has_nectarcam() and stage1.const_run_config().const_nectarcam().module_size() > 0):
        nc = stage1.const_run_config().const_nectarcam()
        return numpy.median([nc.const_module(im).hvpa_voltage() for im in range(nc.module_size())])
    elif(stage1.const_nectarcam().const_ancillary_data().num_hvpa_voltage_measurements()>0):
        return median_measured_voltage(stage1)
    else:
        return -1

def median_feb_temp(stage1, temperature_set=1):
    if(not stage1.has_nectarcam() or not stage1.const_nectarcam().has_ancillary_data()):
        return None
    nca_data = stage1.const_nectarcam().const_ancillary_data()
    k = nca_data.feb_temperature_keys();
    if(len(k) == 0):
        return None
    t = numpy.zeros(len(k))
    for i,k in enumerate(k):
        all_m = nca_data.const_feb_temperature(k)
        tm = []
        for j in range(all_m.measurement_size()):
            if(temperature_set == 1):
                tm.append(all_m.const_measurement(j).tfeb1())
            else:
                tm.append(all_m.const_measurement(j).tfeb2())
        t[i] = numpy.median(tm)
    return numpy.median(t)

def spread_feb_temp(stage1, temperature_set=1):
    if(not stage1.has_nectarcam() or not stage1.const_nectarcam().has_ancillary_data()):
        return None
    nca_data = stage1.const_nectarcam().const_ancillary_data()
    k = nca_data.feb_temperature_keys();
    if(len(k) == 0):
        return None
    t = numpy.zeros(len(k))
    for i,k in enumerate(k):
        all_m = nca_data.const_feb_temperature(k)
        tm = []
        for j in range(all_m.measurement_size()):
            if(temperature_set == 1):
                tm.append(all_m.const_measurement(j).tfeb1())
            else:
                tm.append(all_m.const_measurement(j).tfeb2())
        t[i] = numpy.median(tm)
    return numpy.max(t)-numpy.min(t)

def run_duration(stage1):
    run_info = stage1.const_run_info()
    camera_layout = stage1.const_run_config().const_camera_layout()
    camera_clocks = [ 'UCTS timestamp',
                      'UCTS combined 10MHz and pps counter', 
                      'TIB combined 10MHz and pps counter', 
                      'FEB local pps and 2ns TDC counter sum',
                      'SWAT timestamp',
                      'EVB timestamp' ]
    camera_clocks_run_duration = []
    for cn in camera_clocks:
        try:
            ic = camera_layout.camera_clock_name().index(cn)
            if(run_info.camera_clock_presence(ic) > run_info.num_events_found()//2 and 
               run_info.camera_clock_min_time(ic) > 0):
                camera_clocks_run_duration.append(
                    (run_info.camera_clock_max_time(ic)-run_info.camera_clock_min_time(ic))/
                        camera_layout.camera_clock_frequency(ic))
        except:
            pass
    if(run_info.min_event_time()>0):
        camera_clocks_run_duration.append((run_info.max_event_time()-run_info.min_event_time())*1e-9)
    camera_clocks_run_duration = numpy.asarray(camera_clocks_run_duration)
    camera_clocks_run_duration = camera_clocks_run_duration[camera_clocks_run_duration>0]
    if(len(camera_clocks_run_duration) == 0):
        return None
    median_run_duration = numpy.median(camera_clocks_run_duration)
    for the_run_duration in camera_clocks_run_duration:
        if(numpy.abs(the_run_duration - median_run_duration) < 1):
            return the_run_duration
    return None

def num_wf(stage1):
    nwfp = calin.ix.math.histogram.Histogram1DData()
    if(stage1.has_wf_hists_pedestal()):
        wfh = stage1.const_wf_hists_pedestal()
        if(wfh.has_dual_gain_camera() and wfh.const_dual_gain_camera().has_nchan_present()):
            nwfp.IntegrateFrom(wfh.const_dual_gain_camera().const_nchan_present())
    if(stage1.has_wf_hists_physics()):
        wfh = stage1.const_wf_hists_physics()
        if(wfh.has_dual_gain_camera() and wfh.const_dual_gain_camera().has_nchan_present()):
            nwfp.IntegrateFrom(wfh.const_dual_gain_camera().const_nchan_present())
    if(stage1.has_wf_hists_external_flasher()):
        wfh = stage1.const_wf_hists_external_flasher()
        if(wfh.has_dual_gain_camera() and wfh.const_dual_gain_camera().has_nchan_present()):
            nwfp.IntegrateFrom(wfh.const_dual_gain_camera().const_nchan_present())
    if(stage1.has_wf_hists_internal_flasher()):
        wfh = stage1.const_wf_hists_internal_flasher()
        if(wfh.has_dual_gain_camera() and wfh.const_dual_gain_camera().has_nchan_present()):
            nwfp.IntegrateFrom(wfh.const_dual_gain_camera().const_nchan_present())
    if(nwfp.sum_w() == 0 or stage1.const_run_config().configured_channel_id_size() == 0):
        return 0
    return nwfp.sum_wx()/nwfp.sum_w()/stage1.const_run_config().configured_channel_id_size()

def analyze_charge_hists(all_hist, ped=None, pedvar0=None, evf=1.2, flasher_resolution=0,
        dataset = 'opt_win', pxl=0.02, pxr = 0.98):
    all_xl = []
    all_xc = []
    all_xr = []

    xscale = (scipy.special.erfinv(pxr*2-1) - scipy.special.erfinv(pxl*2-1))*numpy.sqrt(2)

    for i in range(len(all_hist)):
        if(dataset == 'opt_win' and all_hist[i].has_opt_win_qsum()):
            h = calin.math.histogram.densify(all_hist[i].const_opt_win_qsum())
        elif(dataset == 'sig_win' and all_hist[i].has_sig_win_qsum()):
            h = calin.math.histogram.densify(all_hist[i].const_sig_win_qsum())
        elif(dataset == 'ped_win' and all_hist[i].has_ped_win_qsum()):
            h = calin.math.histogram.densify(all_hist[i].const_ped_win_qsum())
        else:
            raise RuntimeError('Dataset not found : '+dataset)

        x = h.xval0() + h.dxval()*numpy.arange(h.bins_size()+1)
        y = numpy.append(0, numpy.cumsum(h.bins()))
        xl,xc,xr = numpy.interp(y[-1]*numpy.asarray([pxl,0.5,pxr]),y,x)

        all_xl.append(xl)
        all_xc.append(xc)
        all_xr.append(xr)

    all_xl = numpy.asarray(all_xl)
    all_xc = numpy.asarray(all_xc)
    all_xr = numpy.asarray(all_xr)
    all_rms = (all_xr-all_xl)/xscale

    all_gain = None
    all_intensity = None
    if ped is not None:
        all_gain = []
        all_intensity = []

        all_var = all_rms**2
        all_mean = all_xc - ped

        if pedvar0 is not None:
            all_var -= pedvar0
        if flasher_resolution > 0:
            all_var -= (all_mean*flasher_resolution)**2

        all_gain = all_var/all_mean/evf
        all_intensity = all_mean/all_gain

    return all_xl, all_xc, all_xr, all_rms, all_gain, all_intensity

def dhgauss_cdf(x, loc, lscale, rscale):
    xx = x if isinstance(x, collections.abc.Iterable) else numpy.asarray([x])
    lnorm = 2/(1+rscale/lscale)
    rnorm = rscale/lscale*lnorm
    y = numpy.zeros_like(xx,dtype=float)
    y[xx<loc] = scipy.stats.norm.cdf(xx[xx<loc],loc=loc,scale=lscale)*lnorm
    y[xx>=loc] = (scipy.stats.norm.cdf(xx[xx>=loc],loc=loc,scale=rscale)-0.5)*rnorm + 0.5*lnorm
    return y if isinstance(x, collections.abc.Iterable) else y[0]

def dhgauss_sf(x, loc, lscale, rscale):
    xx = x if isinstance(x, collections.abc.Iterable) else numpy.asarray([x])
    lnorm = 2/(1+rscale/lscale)
    rnorm = rscale/lscale*lnorm
    y = numpy.zeros_like(xx,dtype=float)
    y[xx<loc] = (scipy.stats.norm.sf(xx[xx<loc],loc=loc,scale=lscale)-0.5)*lnorm + 0.5*rnorm
    y[xx>=loc] = scipy.stats.norm.sf(xx[xx>=loc],loc=loc,scale=rscale)*rnorm
    return y if isinstance(x, collections.abc.Iterable) else y[0]

def dhgauss_percentile(frac, loc, lscale, rscale):
    lnorm = 2/(1+rscale/lscale)
    rnorm = rscale/lscale*lnorm
    if(2*frac < lnorm):
        return scipy.stats.norm.isf(1 - frac/lnorm,loc=loc,scale=lscale)
    else:
        return scipy.stats.norm.isf(1 - (frac - 0.5*lnorm)/rnorm - 0.5,loc=loc,scale=rscale)
        
def dhgauss_median(loc, lscale, rscale):
    return dhgauss_percentile(0.5, loc, lscale, rscale)

def analyze_trigger_thresholds(hset, hclr, ped, nmin=100, do_mle=False, mle_pmin=1e-200):
    xsetl = int(hset.xval0()/hset.dxval())
    xclrl = int(hclr.xval0()/hclr.dxval())
    xsetr = xsetl + hset.bins_size()
    xclrr = xclrl + hclr.bins_size()

    xl = min(xsetl, xclrl)
    xr = max(xsetr, xclrr)

    x = numpy.arange(xl,xr)
    yset = numpy.zeros_like(x)
    yclr = numpy.zeros_like(x)
    yset[xsetl-xl:xsetr-xl] = hset.bins()
    yclr[xclrl-xl:xclrr-xl] = hclr.bins()
    ytot = yset+yclr
    
    m = numpy.bitwise_and(numpy.bitwise_and(x>=xsetl-5, x<=xclrr+5), ytot>0)
    xp = x[m]*hset.dxval() - ped
    p = yset[m]/ytot[m]

    nset = sum(yset[x<=xclrr])
    nclr = sum(yclr[x>=xsetl])

    xmed = numpy.nan
    xiqr = numpy.nan
    res_sq = numpy.nan
    xfit = None
    if(nset>=nmin and nclr>=nmin):
        fitbounds=(0,numpy.inf)

        x0_50 = xp[numpy.argmin((p-0.5)**2)]
        x0_iqr = min(numpy.abs(xp[numpy.argmin((p-0.75)**2)]-xp[numpy.argmin((p-0.25)**2)]),4)
        
        x0 = (x0_50, 0.5*x0_iqr, 0.5*x0_iqr)
        def residual(x):
            return p - dhgauss_cdf(xp, *x)
        optres = scipy.optimize.least_squares(residual, x0, bounds=fitbounds)

        if(do_mle):
            x1 = optres.x
            def logL(x):
                b = numpy.maximum(dhgauss_cdf(xp, *x), mle_pmin)
                omb = numpy.maximum(dhgauss_sf(xp, *x), mle_pmin)
                lp = yset[m]*numpy.log(b) + yclr[m]*numpy.log(omb)
                return -sum(lp)
            optres = scipy.optimize.minimize(logL, x1, tol=1e-3, bounds=(fitbounds,)*3)

            if(optres.success==False):
                optres = scipy.optimize.minimize(logL, x1, tol=1e-3, method='Nelder-Mead', bounds=(fitbounds,)*3)

        if(optres.success==True):
            xmed = dhgauss_median(*optres.x)
            xiqr = dhgauss_percentile(0.75,*optres.x) - dhgauss_percentile(0.25,*optres.x)
            res_sq = sum(residual(optres.x)**2)
            xfit = optres.x

    return xp, p, nset, nclr, xfit, xmed, xiqr, res_sq