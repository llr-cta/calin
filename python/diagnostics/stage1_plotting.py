# calin/python/diagnostics/stage1_plotting.py -- Stephen Fegan -- 2021-06-23
#
# Rendering of custom stage-one diagnostics plots
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

import matplotlib.pyplot
import matplotlib.collections
import matplotlib.patches
import matplotlib.colors
import numpy
import scipy.special

import calin.plotting
import calin.diagnostics.stage1
import calin.diagnostics.stage1_analysis
import calin.iact_data.instrument_layout
import calin.math.fftw_util
import calin.util.string

def trigger_type_title(trigger_type):
    if(trigger_type == 'physics'):
        return 'Phys'
    elif(trigger_type == 'pedestal'):
        return 'Ped'
    elif(trigger_type == 'external_flasher'):
        return 'ExtFlash'
    elif(trigger_type == 'internal_flasher'):
        return 'IntFlash'
    elif(trigger_type == 'all'):
        return 'All'
    return ''

def trigger_type_and_gain_title(trigger_type, low_gain = False):
    gain_title = ', LG' if low_gain else ', HG'
    if(trigger_type == 'physics'):
        return 'Phys' + gain_title
    elif(trigger_type == 'pedestal'):
        return 'Ped' + gain_title
    elif(trigger_type == 'external_flasher'):
        return 'ExtFlash' + gain_title
    elif(trigger_type == 'internal_flasher'):
        return 'IntFlash' + gain_title
    elif(trigger_type == 'all'):
        return 'All' + gain_title
    return ''

def dc_units(nsamp=0, low_gain = False):
    return 'DC$_{%d}$'%(nsamp) if nsamp>1 else 'DC'

def draw_channel_event_fraction(stage1, channel_count, cb_label=None, log_scale=True,
        lin_scale_suppress_zero=False, cmap = 'CMRmap_r', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75, stat_format='%.3e'):
    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    ri = stage1.const_run_info()

    max_xy = max(numpy.max(numpy.abs(cl.outline_polygon_vertex_x())),
                 numpy.max(numpy.abs(cl.outline_polygon_vertex_y())))

    evts_ondisk    = int(ri.num_events_found())

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    pc = calin.plotting.plot_camera_image(channel_count/evts_ondisk, cl,
        configured_channels=rc.configured_channel_id(), cmap=cmap,
        draw_outline=draw_outline, pix_lw=pix_lw,
        outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True,
        stats_format=stat_format, stats_fontsize=stat_label_fontsize)

    if(log_scale):
        vmin = 1/evts_ondisk**1.1
        vmin_color_change = vmin**0.45

        pc.set_norm(matplotlib.colors.LogNorm(vmin=vmin,vmax=1.0))
        axis.axis(numpy.asarray([-1,1,-1,1])*1.05*max_xy)

        cb = axis.get_figure().colorbar(pc, ax=axis, label=cb_label)
        cb.ax.plot([-max_xy, max_xy], [1.0/evts_ondisk, 1.0/evts_ondisk], 'g-', lw=0.75)
    else:
        if(not lin_scale_suppress_zero):
            pc.set_clim(0,pc.get_clim()[1])
        cb = axis.get_figure().colorbar(pc, ax=axis, label=cb_label)
        cb.formatter.set_powerlimits((-2, 2))
        cb.formatter.set_useMathText(True)
        cb.update_ticks()

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

def draw_log_delta_t_histogram(stage1, event_set = 'all', axis = None):
    if(axis is None):
        axis = matplotlib.pyplot.gca()

    if(type(event_set) is str):
        event_set = [event_set]

    ri = stage1.const_run_info()
    for set in event_set:
        if(set == 'physics'):
            dt_h = ri.const_log10_delta_t_histogram_trigger_physics()
            dt2_h = ri.const_log10_delta2_t_histogram_trigger_physics()
        elif(set == 'consecutive'):
            dt_h = ri.const_log10_delta_t_histogram()
            dt2_h = ri.const_log10_delta2_t_histogram()
        else:
            dt_h = ri.const_log10_delta_t_histogram_all_recorded()
            dt2_h = None

        calin.plotting.plot_histogram(dt_h,axis=axis,xoffset=6,xscale_as_log10=True,
            normalise=True,density=True,label='Delta-T : $t_{i+1}-t_i$')
        if(dt2_h is not None and dt2_h.sum_w()>0):
            calin.plotting.plot_histogram(dt2_h,axis=axis,xoffset=6,xscale_as_log10=True,
                normalise=True,density=True,label='Delta-2T : $t_{i+2}-t_i$')

    axis.set_yscale('log')
    axis.set_xlabel('Time difference [us]')
    axis.set_ylabel('Density [1]')
    axis.legend(loc='lower center')
    axis.grid()

    if(dt_h.sum_w()>0): # protect against calling weight(0) on empty histogram
        dt_sh = calin.math.histogram.SimpleHist(dt_h)
        label_x = 10**(dt_sh.min_xval()+6)
        label_y = dt_sh.weight(0)/dt_sh.sum_w()/dt_sh.dxval()*0.4
        axis.text(label_x, label_y, '%.3fus'%(10**(dt_sh.min_xval()+6)), ha='left', va='bottom')
        if(label_y < axis.get_ylim()[0]):
            axis.set_ylim(label_y, axis.get_ylim()[1])

    if(dt2_h is not None and dt2_h.sum_w()>0): # protect against calling weight(0) on empty histogram
        dt2_sh = calin.math.histogram.SimpleHist(dt2_h)
        label_x = 10**(dt2_sh.min_xval()+6)
        label_y = dt2_sh.weight(0)/dt2_sh.sum_w()/dt2_sh.dxval()*0.5
        axis.text(label_x, label_y, '%.2fus'%(10**(dt2_sh.min_xval()+6)), ha='left', va='bottom')
        if(label_y < axis.get_ylim()[0]):
            axis.set_ylim(label_y, axis.get_ylim()[1])

    def to_khz(x):
        x = numpy.array(x).astype(float)
        near_zero = numpy.isclose(x, 0)
        x[near_zero] = numpy.inf
        x[~near_zero] = 1e3 / x[~near_zero]
        return x

    def from_khz(x):
        x = numpy.array(x).astype(float)
        near_zero = numpy.isclose(x, 0)
        x[near_zero] = numpy.inf
        x[~near_zero] = 1e-3 / x[~near_zero]
        return x

    secax = axis.secondary_xaxis('top', functions=(to_khz, from_khz))
    secax.set_xlabel('Frequency [kHz]')

def draw_pedestal_value(stage1, all_events_ped_win=False, low_gain=False,
        cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    charge_stats = stage1.const_charge_stats().const_low_gain() if low_gain \
        else stage1.const_charge_stats().const_high_gain()

    if(all_events_ped_win):
        nsamp = stage1.const_config().const_low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.config().has_low_gain_opt_sum()) \
            else stage1.const_config().const_high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_mean()/nsamp
    else:
        nsamp = rc.num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_mean()/nsamp

    DC = dc_units(1, low_gain)

    pc = calin.plotting.plot_camera_image(
        values, cl, channel_mask=nevent>0, cmap=cmap,
        configured_channels=rc.configured_channel_id(),
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True, stats_format='%.2f '+DC)

    cb = axis.get_figure().colorbar(pc, ax=axis, label='Pedestal mean [%s]'%DC)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

def draw_pedestal_rms(stage1, all_events_ped_win=False, low_gain=False,
        cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    charge_stats = stage1.const_charge_stats().const_low_gain() if low_gain \
        else stage1.const_charge_stats().const_high_gain()

    if(all_events_ped_win):
        nsamp = stage1.const_config().const_low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.config().has_low_gain_opt_sum()) \
            else stage1.const_config().const_high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_var()
    else:
        nsamp = rc.num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_var()

    data = numpy.sqrt(values)
    mask = nevent>0

    DC = dc_units(nsamp, low_gain)

    pc = calin.plotting.plot_camera_image(
        data, cl, channel_mask=mask, cmap=cmap,
        configured_channels=rc.configured_channel_id(),
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True, stats_format=f'%.2f '+DC,
        draw_top12_val=True)

    cb = calin.plotting.add_colorbar_and_clipping(axis, pc, data, mask=mask, percentile=99.5,
            camera_layout=cl, configured_channels=rc.configured_channel_id(),
            cb_label='Pedestal %d-sample RMS [%s]'%(nsamp,DC))

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc, cb

def draw_pedestal_trend(stage1, pedvar=False, all_events_ped_win=False, low_gain=False, axis = None):
    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    charge_stats = stage1.const_charge_stats().const_low_gain() if low_gain \
        else stage1.const_charge_stats().const_high_gain()

    if(all_events_ped_win):
        nsamp = stage1.const_config().const_low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.const_config().has_low_gain_opt_sum()) \
            else stage1.const_config().const_high_gain_opt_sum().integration_n()
        nchan = charge_stats.all_trigger_ped_win_count_vs_time_size();
    else:
        nsamp = rc.num_samples()
        nchan = charge_stats.ped_trigger_full_wf_count_vs_time_size();

    for ichan in range(nchan):
        if(all_events_ped_win):
            nevent = charge_stats.const_all_trigger_ped_win_count_vs_time(ichan)
            values = charge_stats.const_all_trigger_ped_win_var_vs_time(ichan) if pedvar \
                else charge_stats.const_all_trigger_ped_win_mean_vs_time(ichan)
        else:
            nevent = charge_stats.const_ped_trigger_full_wf_count_vs_time(ichan)
            values = charge_stats.const_ped_trigger_full_wf_var_vs_time(ichan) if pedvar \
                else charge_stats.const_ped_trigger_full_wf_mean_vs_time(ichan)

        nevent = calin.math.histogram.densify(nevent)
        values = calin.math.histogram.densify(values)

        v = numpy.sqrt(values.bins()) if pedvar else values.bins()/nsamp
        t = values.xval0() + values.dxval()*(numpy.arange(values.bins_size()) + 0.5)
        i0 = int(round((values.xval0()-nevent.xval0())/values.dxval()))
        m = nevent.bins()[i0:i0+values.bins_size()] > 5

        if(numpy.count_nonzero(m) > 0):
            axis.plot(t[m], v[m], 'k', alpha=0.1)

    if(all_events_ped_win):
        nevent = charge_stats.const_camera_all_trigger_ped_win_count_vs_time()
        values = charge_stats.const_camera_all_trigger_ped_win_var_vs_time() if pedvar \
            else charge_stats.const_camera_all_trigger_ped_win_mean_vs_time()
    else:
        nevent = charge_stats.const_camera_ped_trigger_full_wf_count_vs_time()
        values = charge_stats.const_camera_ped_trigger_full_wf_var_vs_time() if pedvar \
            else charge_stats.const_camera_ped_trigger_full_wf_mean_vs_time()

    nevent = calin.math.histogram.densify(nevent)
    values = calin.math.histogram.densify(values)

    v = numpy.sqrt(values.bins()/nchan) if pedvar else values.bins()/nsamp/nchan
    t = values.xval0() + values.dxval()*(numpy.arange(values.bins_size()) + 0.5)
    i0 = int(round((values.xval0()-nevent.xval0())/values.dxval()))
    m = nevent.bins()[i0:i0+values.bins_size()] > 5
    
    axis.plot(t[m], v[m], 'r')
    axis.set_xlabel('Elapsed time [s]')
    if(pedvar):
        DC = dc_units(nsamp, low_gain)
        axis.set_ylabel('Pedestal %d-sample RMS [%s]'%(nsamp,DC))
    else:
        DC = dc_units(1, low_gain)
        axis.set_ylabel('Pedestal mean [%s]'%DC)

def draw_pedestal_plots(stage1, figure_factory = calin.plotting.PyPlotFigureFactory(),
        cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    fig_dict = dict()

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_high_gain()
            and max(stage1.const_charge_stats().const_high_gain().ped_trigger_event_count())>0):
        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_value(stage1,all_events_ped_win=False,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal mean (ped events), run : %d'%stage1.run_number())
        fig_dict['pedestal_mean_hg_ped_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_rms(stage1,all_events_ped_win=False,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal rms (ped events), run : %d'%stage1.run_number())
        fig_dict['pedestal_rms_hg_ped_evt'] = [ fig, axis ]

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_low_gain()
            and max(stage1.const_charge_stats().const_low_gain().ped_trigger_event_count())>0):
        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_value(stage1,all_events_ped_win=False,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal mean (ped events), run : %d'%stage1.run_number())
        fig_dict['pedestal_mean_lg_ped_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_rms(stage1,all_events_ped_win=False,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal rms (ped events), run : %d'%stage1.run_number())
        fig_dict['pedestal_rms_lg_ped_evt'] = [ fig, axis ]

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_high_gain()
            and max(stage1.const_charge_stats().const_high_gain().all_trigger_event_count())>0):
        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_value(stage1,all_events_ped_win=True,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal mean (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_mean_hg_all_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_rms(stage1,all_events_ped_win=True,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal rms (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_rms_hg_all_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_trend(stage1, pedvar=False, all_events_ped_win=True, low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal mean trend (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_trend_mean_hg_all_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_trend(stage1, pedvar=True, all_events_ped_win=True, low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal rms trend (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_trend_rms_hg_all_evt'] = [ fig, axis ]

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_low_gain()
            and max(stage1.const_charge_stats().const_low_gain().all_trigger_event_count())>0):
        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_value(stage1,all_events_ped_win=True,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal mean (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_mean_lg_all_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_rms(stage1,all_events_ped_win=True,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal rms (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_rms_lg_all_evt'] = [ fig, axis ]

        fig, ax = figure_factory.new_camera_figure()
        draw_pedestal_trend(stage1, pedvar=False, all_events_ped_win=True, low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal mean trend (all events), run : %d'%stage1.run_number())
        fig_dict['pedestal_trend_mean_lg_all_evt'] = [ fig, axis ]

    return fig_dict

def draw_missing_components_fraction(stage1, cmap = 'CMRmap_r',
        figure_factory = calin.plotting.PyPlotFigureFactory(),
        draw_outline=True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, aux_label_fontsize=5.5, stat_label_fontsize=4.75):
    SQRT3_2 = numpy.sqrt(3)/2
    RAD_30 = 30.0/180.0*numpy.pi

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    ri = stage1.const_run_info()

    fig_dict = dict()

    if(ri.num_events_found() == 0):
        return fig_dict

    fig, axis = figure_factory.new_camera_figure()
    fig_dict['missing_components'] = [ fig, axis ]

    max_xy = max(numpy.max(numpy.abs(cl.outline_polygon_vertex_x())),
                 numpy.max(numpy.abs(cl.outline_polygon_vertex_y())))

    evts_ondisk    = int(ri.num_events_found())
    evts_maxid     = int(ri.event_numbers_found().end_index()[-1])
    evts_duplicate = int(sum(ri.duplicate_event_numbers().count()))

    num_event_missing = evts_maxid-evts_ondisk-1+evts_duplicate
    frac_event_missing = float(num_event_missing)/float(evts_maxid)

    additional_polygons = [ ]
    additional_polygon_data = [ ri.num_events_missing_tib()/evts_ondisk,
                                ri.num_events_missing_cdts()/evts_ondisk,
                                ri.num_events_missing_tib_and_cdts()/evts_ondisk,
                                frac_event_missing]

    polygon_name = [ 'tib', 'cdts', 'both', 'lost\nevent' ]

    rad = 10/120*max_xy
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad,max_xy-rad*SQRT3_2),
                                          numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad,max_xy-rad*SQRT3_2-2*rad),
                                          numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad+2*rad*SQRT3_2,max_xy-rad*SQRT3_2-rad),
                                          numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.Circle(xy=(max_xy-1.5*rad,-max_xy+1.5*rad), radius=1.5*rad))

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    modmissing = numpy.asarray([max(float(evts_ondisk-ri.module(i).num_events_present()),1e-10)/evts_ondisk \
        for i in range(ri.module_size())],dtype=float)

    pc = calin.plotting.plot_camera_module_image(modmissing, cl,
        configured_modules=rc.configured_module_id(), cmap=cmap,
        additional_polygons=additional_polygons,
        additional_polygon_data=additional_polygon_data,
        draw_outline=draw_outline, mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_modules=True)

    vmin = 1/evts_ondisk**1.1
    vmin_color_change = vmin**0.45
    def label_color(value):
        return 'k' if value < vmin_color_change else 'w'

    if(mod_label_fontsize):
        for imod, cmid in enumerate(rc.configured_module_id()):
            m = rc.camera_layout().module(int(cmid))
            axis.text(m.x(), m.y(), '%d'%(m.module_index()), ha='center', va='center',
                fontsize=mod_label_fontsize, color=label_color(modmissing[imod]))

    for ip, pp in enumerate(additional_polygons):
        if(draw_outline):
            if(isinstance(pp, matplotlib.patches.RegularPolygon)):
                xy = pp.xy
                p = matplotlib.patches.RegularPolygon(xy=xy, numVertices=pp.numvertices, 
                                                      radius=pp.radius, orientation=pp.orientation)
            else:
                xy = pp.center
                p = matplotlib.patches.Circle(xy=xy, radius=pp.radius)
            p.set_linewidth(outline_lw)
            p.set_edgecolor(outline_color)
            p.set_fill(False)
            axis.add_patch(p)
        if(aux_label_fontsize):
            axis.text(xy[0], xy[1], polygon_name[ip], ha='center', va='center',
            fontsize=aux_label_fontsize, color=label_color(additional_polygon_data[ip]))

    if(stat_label_fontsize):
        axis.text(-max_xy,-max_xy,
            f'Event counts\nOn disk:   {evts_ondisk:,}\nTriggered: {evts_maxid-1:,}\nDuplicate: {evts_duplicate:,}',
            ha='left', va='bottom', fontfamily='monospace',
            fontsize=stat_label_fontsize);

    pc.set_norm(matplotlib.colors.LogNorm(vmin=vmin,vmax=1.0))
    axis.axis(numpy.asarray([-1,1,-1,1])*1.05*max_xy)
    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)
    axis.set_title('Missing components, run : %d'%stage1.run_number())

    cb = fig.colorbar(pc, ax=axis, label='Component missing fraction')
    cb.ax.plot([-max_xy, max_xy], [1.0/evts_ondisk, 1.0/evts_ondisk], 'g-', lw=0.75)

    return fig_dict

def draw_module_dataorder(stage1, cmap = 'inferno', axis=None,
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4):

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    data = numpy.arange(stage1.run_config().configured_module_id_size())

    pc = calin.plotting.plot_camera_module_image(data, stage1.run_config().camera_layout(),
                    configured_modules=stage1.run_config().configured_module_id(),
                    axis=axis, cmap=cmap, draw_outline=True, draw_stats=False,
                    mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                    hatch_missing_modules=True)
    cb = axis.get_figure().colorbar(pc, label='Module data order')

    if(mod_label_fontsize is not None and mod_label_fontsize>0):
        calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                configured_modules=stage1.run_config().configured_module_id(),
                                pc=pc, module_values = data, fontsize=mod_label_fontsize)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

def draw_channel_dataorder(stage1, cmap = 'inferno', axis=None,
        draw_outline = True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888'):

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    data = numpy.arange(stage1.run_config().configured_channel_id_size())

    pc = calin.plotting.plot_camera_image(data, stage1.run_config().camera_layout(),
                    configured_channels=stage1.run_config().configured_channel_id(),
                    axis=axis, cmap=cmap, draw_outline=True, draw_stats=False,
                    pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                    hatch_missing_channels=True)
    cb = axis.get_figure().colorbar(pc, label='Channel data order')

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

def draw_nectarcam_feb_temperatures(stage1, temperature_set=1, 
        figure_factory = calin.plotting.PyPlotFigureFactory(),
        tcmap = 'coolwarm', dtcmap = 'inferno', draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):

    # We are assured on entry that there is some FEB data temperature in the NectarCAM ancillary data
    ncad = stage1.const_nectarcam().const_ancillary_data()

    fig_trend, axis_trend = figure_factory.new_histogram_figure()

    def minmax(x):
        return numpy.max(x)-numpy.min(x)

    tfeb = []
    dtfeb = []
    mask = []

    t0mt = stage1.const_run_config().const_run_start_time().time_ns()//1000000000

    method = 'tfeb1' if temperature_set==1 else 'tfeb2'
    for modid in stage1.const_run_config().configured_module_id() :
        if(ncad.feb_temperature_has_key(int(modid))):
            measurement_set = ncad.const_feb_temperature(int(modid))
            num_measurments = measurement_set.measurement_size()
            mt = numpy.asarray([getattr(measurement_set.measurement(i), method)() for i in range(num_measurments)])
            tmt = numpy.asarray([measurement_set.measurement(i).time() for i in range(num_measurments)])
            mmask = numpy.bitwise_and(mt>0, mt<100.0)
            if(numpy.count_nonzero(mmask)):
                axis_trend.plot(tmt[mmask] - t0mt, mt[mmask],'k',alpha=0.1)
                tfeb.append(numpy.mean(mt[mmask]))
                dtfeb.append(minmax(mt[mmask]))
                mask.append(True)
            else:
                tfeb.append(numpy.nan)
                dtfeb.append(numpy.nan)
                mask.append(False)
        else:
            tfeb.append(numpy.nan)
            dtfeb.append(numpy.nan)
            mask.append(False)
    axis_trend.set_xlabel('Elapsed run time [s]')
    axis_trend.set_ylabel('FEB temperature %d [\u2103]'%temperature_set)
    axis_trend.set_title('Temperature trend (FEB %d), run : %d'%(temperature_set, stage1.run_number()))

    tecc = None
    dtecc = None

    if(ncad.has_ecc_measurements() and \
            ncad.const_ecc_measurements().measurement_size()>0):
        measurement_set = ncad.const_ecc_measurements()
        num_measurments = measurement_set.measurement_size()

        tecc = []
        dtecc = []

        for itemp in range(16):
            method = f'temp_{itemp+1:02}'
            et = numpy.asarray([getattr(measurement_set.measurement(i), method)() for i in range(num_measurments)])
            emask = numpy.bitwise_and(et!=0, numpy.bitwise_and(et>-30.0, et<50.0))
            if(numpy.count_nonzero(mmask)):
                tecc.append(' %4.1f'%numpy.mean(et[emask]))
                dtecc.append(' %4.1f'%minmax(et[emask]))
            else:
                tecc.append(' ====')
                dtecc.append(' ====')

    camera_layout = stage1.const_run_config().const_camera_layout()
    configured_modules = stage1.const_run_config().configured_module_id()

    def do_draw_tfeb_camera_plot(feb_dataset, ecc_dataset, label, cmap):
        def measurement_val(m, fn):
            m = numpy.asarray(m)
            return ' ====' if(numpy.all(m==0) or numpy.any(m<-30) or numpy.any(m>50)) else ' %4.1f'%fn(m)

        fig, axis = figure_factory.new_camera_figure()
        pc = calin.plotting.plot_camera_module_image(feb_dataset, camera_layout,
                    configured_modules=configured_modules,
                    module_mask=mask, axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                    mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                    hatch_missing_modules=True, stats_format=u'%4.2f\u00b0C',
                    stats_fontsize=stat_label_fontsize)
        fig.colorbar(pc, label=u'FEB temperature %d %s[\u2103]'%(temperature_set,label))
        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, camera_layout,
                                    configured_modules=configured_modules,
                                    pc=pc, module_values = feb_dataset, fontsize=mod_label_fontsize)
        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

        if(ecc_dataset is not None):
            text_ecc = 'ECC'
            text_ecc += ecc_dataset[0] + ecc_dataset[1] + ecc_dataset[2] + ecc_dataset[3]
            text_ecc += u'\n\u00b0C'
            text_ecc += ecc_dataset[4] + ecc_dataset[5] + ecc_dataset[6] + ecc_dataset[7]
            text_ecc += '\n'
            text_ecc += ecc_dataset[8] + ecc_dataset[9] + ecc_dataset[10] + ecc_dataset[11]

            text_ecc4 = ecc_dataset[12] + ecc_dataset[13] + ecc_dataset[14] + ecc_dataset[15]
            if(text_ecc4 != '==== ==== ==== ===='):
                text_ecc += '\n' + text_ecc4

            max_xy = max(numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_x())),
                            numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_y())))
            axis.text(max_xy,max_xy,text_ecc,ha='right',va='top',fontfamily='monospace',fontsize=stat_label_fontsize)

        return fig, axis, pc
    
    fig_tfeb, axis_tfeb, _ = do_draw_tfeb_camera_plot(tfeb,tecc,'',tcmap)
    axis_tfeb.set_title('Temperature (FEB %d), run : %d'%(temperature_set, stage1.run_number()))

    fig_dtfeb, axis_dtfeb, _ = do_draw_tfeb_camera_plot(dtfeb,dtecc,'range (max-min) ',dtcmap)
    axis_dtfeb.set_title('Temperature spread (FEB %d), run : %d'%(temperature_set, stage1.run_number()))

    fig_dict = dict()
    fig_dict['temperature_feb%d'%temperature_set] = [ fig_tfeb, axis_tfeb ]
    fig_dict['temperature_spread_feb%d'%temperature_set] = [ fig_dtfeb, axis_dtfeb ]
    fig_dict['temperature_trend_feb%d'%temperature_set] = [ fig_trend, axis_trend ]
    return fig_dict

def draw_all_clock_regression(stage1, figure_factory = calin.plotting.PyPlotFigureFactory(),
        clockid=0, cmap = 'inferno',
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    max_xy = max(numpy.max(numpy.abs(cl.outline_polygon_vertex_x())),
                 numpy.max(numpy.abs(cl.outline_polygon_vertex_y())))

    freq_offset_ppm, freq_spread_ppm, time_offset_ns, time_spread_ns, d2_per_event, n_problem_bins = \
        calin.diagnostics.stage1_analysis.summarize_module_clock_regression(stage1, clockid)
    mask = numpy.isfinite(freq_offset_ppm)

    cam_freq_offset_ppm, cam_freq_spread_ppm, cam_time_offset_ns, cam_time_spread_ns, cam_d2_per_event = \
        calin.diagnostics.stage1_analysis.summarize_camera_clock_regressions(stage1)

    def draw_it(axis, mod_data, cam_data, stats_format, cb_label):
        pc = calin.plotting.plot_camera_module_image(mod_data, stage1.run_config().camera_layout(),
                        configured_modules=stage1.run_config().configured_module_id(),
                        module_mask=mask, axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                        mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_modules=True, stats_format=stats_format,
                        stats_fontsize=stat_label_fontsize, draw_top12_val=True)
        # cb = axis.get_figure().colorbar(pc, label=cb_label)
        cb = calin.plotting.add_colorbar_and_clipping(axis, pc, mod_data, mask=mask, percentile=99.5,
                camera_layout=cl, configured_modules=rc.configured_module_id(),
                cb_label=cb_label)

        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                    configured_modules=stage1.run_config().configured_module_id(),
                                    pc=pc, module_values=mod_data,
                                    fontsize=mod_label_fontsize)

        if(stat_label_fontsize is not None and stat_label_fontsize>0):
            val_text = [stats_format%val if not numpy.isnan(val) else 'missing' for val in cam_data]
            fmt_string = 'UCTS 10MHz : %s\nTIB 10MHz : %s\nFEB sum : %s'
            axis.text(max_xy,max_xy,fmt_string%(val_text[0], val_text[1], val_text[2]),
                fontsize=stat_label_fontsize, fontfamily='monospace',
                ha='right', va='top')

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

        return pc

    fig_dict = dict()

    if(numpy.count_nonzero(~numpy.isnan(freq_offset_ppm))):
        fig_freq, axis_freq = figure_factory.new_camera_figure()
        draw_it(axis_freq, freq_offset_ppm, cam_freq_offset_ppm, '%.2f ppm',
                'Mean oscillator frequency offset [ppm]')
        axis_freq.set_title('Clock frequency error, run : %d'%stage1.run_number())
        fig_dict['clock_frequency_error'] = [ fig_freq, axis_freq ]

    if(numpy.count_nonzero(~numpy.isnan(time_offset_ns))):
        fig_t0, axis_t0 = figure_factory.new_camera_figure()
        draw_it(axis_t0, time_offset_ns, cam_time_offset_ns, '%.2f ns',
                'Mean UCTS time at counter reset [ns]')
        axis_t0.set_title('Clock offset from UCTS, run : %d'%stage1.run_number())
        fig_dict['clock_offset'] = [ fig_t0, axis_t0 ]

    if(numpy.count_nonzero(~numpy.isnan(d2_per_event))):
        fig_chi2, axis_chi2 = figure_factory.new_camera_figure()
        draw_it(axis_chi2, numpy.sqrt(d2_per_event), numpy.sqrt(cam_d2_per_event), '%.2f ns',
                'Linear-fit RMS residual per event [ns]')
        axis_chi2.set_title('Clock vs UCTS fit RMS residual, run : %d'%stage1.run_number())
        fig_dict['clock_residual'] = [ fig_chi2, axis_chi2 ]

    if(numpy.count_nonzero(~numpy.isnan(freq_spread_ppm))):
        fig_freq_spread, axis_freq_spread = figure_factory.new_camera_figure()
        draw_it(axis_freq_spread, freq_spread_ppm, cam_freq_spread_ppm, '%.3f ppm',
                'Drift in oscillator frequency [ppm]')
        axis_freq_spread.set_title('Clock frequency spread, run : %d'%stage1.run_number())
        fig_dict['clock_frequency_spread'] = [ fig_freq_spread, axis_freq_spread ]

    if(numpy.count_nonzero(~numpy.isnan(time_spread_ns))):
        fig_t0_spread, axis_t0_spread = figure_factory.new_camera_figure()
        draw_it(axis_t0_spread, time_spread_ns, cam_time_spread_ns, '%.3f ns',
                'Drift in UCTS time at counter reset [ns]')
        axis_t0_spread.set_title('Clock offset from UCTS spread, run : %d'%stage1.run_number())
        fig_dict['clock_offset_spread'] = [ fig_t0_spread, axis_t0_spread ]

    if(numpy.count_nonzero(~numpy.isnan(freq_spread_ppm))):
        t0 = stage1.const_run_config().const_run_start_time().time_ns()*1e-9
        fig_freq_trend, axis_freq_trend = figure_factory.new_histogram_figure()
        for imod in range(stage1.run_config().configured_module_id_size()):
            _, all_x0, _, all_a, _, _, all_n = \
                calin.diagnostics.stage1_analysis.get_module_clock_regression_data(stage1, imod, clockid)
            mask_n = all_n>5
            if(numpy.count_nonzero(mask_n)):
                axis_freq_trend.plot(all_x0[mask_n]*1e-9 - t0, 1e6*(all_a[mask_n]-1), 'k', alpha=0.1)

        _, clock_freq, _, all_x0, _, all_a, _, _, all_n = \
                calin.diagnostics.stage1_analysis.get_camera_clock_regression_data(stage1, iclock=2)
        mask_n = all_n>5
        if(numpy.count_nonzero(mask_n)):
            principal_clock_id = stage1.const_clock_regression().principal_clock_id();
            cl = stage1.const_run_config().const_camera_layout()
            principal_clock_freq = cl.camera_clock_frequency(principal_clock_id)
            clock_nominal_a = clock_freq/principal_clock_freq
            axis_freq_trend.plot(all_x0[mask_n]*1e-9 - t0, 1e6*((all_a[mask_n])/clock_nominal_a-1), 'r')
        axis_freq_trend.set_xlabel('Elapsed run time [s]')
        axis_freq_trend.set_ylabel('Clock frequency error [ppm]')
        axis_freq_trend.set_title('Clock frequency trend, run : %d'%stage1.run_number())
        fig_dict['clock_frequency_trend'] = [ fig_freq_trend, axis_freq_trend ]

    return fig_dict

def draw_nectarcam_fpm_measurements(stage1,
        axis_voltage=None, axis_cw_current=None, axis_board_current=None,
        axis_voltage_spread=None, axis_cw_current_spread=None, axis_board_current_spread = None,
        cmap = 'inferno',
        draw_outline = True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):

    v_mean = []
    v_spread = []
    v_mask = []

    iboard_mean = []
    iboard_spread = []
    icw_mean = []
    icw_spread = []
    i_mask = []
    for chanid in stage1.run_config().configured_channel_id() :
        if(stage1.nectarcam().ancillary_data().hvpa_voltage_has_key(int(chanid))):
            measurement_set = stage1.nectarcam().ancillary_data().hvpa_voltage(int(chanid))
            measurements = [measurement_set.measurement(i).voltage() for i in range(measurement_set.measurement_size())]
            v_mean.append(numpy.mean(measurements))
            v_spread.append(numpy.max(measurements)-numpy.min(measurements))
            v_mask.append(True)
        else:
            v_mean.append(numpy.nan)
            v_spread.append(numpy.nan)
            v_mask.append(False)

        if(stage1.nectarcam().ancillary_data().hvpa_current_has_key(int(chanid))):
            measurement_set = stage1.nectarcam().ancillary_data().hvpa_current(int(chanid))
            measurements = [measurement_set.measurement(i).current() for i in range(measurement_set.measurement_size())]
            icw_mean.append(numpy.mean(measurements))
            icw_spread.append(numpy.max(measurements)-numpy.min(measurements))
            measurements = [measurement_set.measurement(i).load_current() for i in range(measurement_set.measurement_size())]
            iboard_mean.append(numpy.mean(measurements))
            iboard_spread.append(numpy.max(measurements)-numpy.min(measurements))
            i_mask.append(True)
        else:
            icw_mean.append(numpy.nan)
            icw_spread.append(numpy.nan)
            iboard_mean.append(numpy.nan)
            iboard_spread.append(numpy.nan)
            i_mask.append(False)

    def draw_it(axis, chan_data, chan_mask, stats_format, cb_label):
        pc = calin.plotting.plot_camera_image(chan_data, stage1.run_config().camera_layout(),
                        configured_channels=stage1.run_config().configured_channel_id(),
                        channel_mask=chan_mask, axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                        pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_channels=True, stats_format=stats_format,
                        stats_fontsize=stat_label_fontsize)
        cb = axis.get_figure().colorbar(pc, label=cb_label)

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

        return pc

    all_pc = []
    if(axis_voltage is not None):
        all_pc.append(draw_it(axis_voltage, v_mean, v_mask,
            '%.1f V', 'Mean measured pixel voltage [V]'))

    if(axis_voltage_spread is not None):
        all_pc.append(draw_it(axis_voltage_spread, v_spread, v_mask,
            '%.1f V', 'Measured pixel voltage spread (max-min) [V]'))

    if(axis_cw_current is not None):
        all_pc.append(draw_it(axis_cw_current, icw_mean, i_mask,
            '%.1f uA', 'Mean measured CW current [uA]'))

    if(axis_cw_current_spread is not None):
        all_pc.append(draw_it(axis_cw_current_spread, icw_spread, i_mask,
            '%.1f uA', 'Measured CW current spread (max-min) [uA]'))

    if(axis_board_current is not None):
        all_pc.append(draw_it(axis_board_current, iboard_mean, i_mask,
            '%.1f uA', 'Mean measured HVPA current[uA]'))

    if(axis_board_current_spread is not None):
        all_pc.append(draw_it(axis_board_current_spread, iboard_spread, i_mask,
            '%.1f uA', 'Measured HVPA current spread (max-min) [uA]'))

    return all_pc

def draw_high_gain_channel_event_fraction(stage1, cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().const_high_gain().all_trigger_event_count()
    pc = draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)
    pc.set_clim([0,1])
    return pc

def draw_low_gain_channel_event_fraction(stage1, cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().const_low_gain().all_trigger_event_count()
    pc = draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)
    pc.set_clim([0,1])
    return pc

def draw_trigger_event_fraction(stage1, cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().channel_triggered_count()
    return draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)

def draw_muon_candidate_trigger_event_fraction(stage1, cmap = 'inferno', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().muon_candidate_channel_triggered_count()
    return draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)

def draw_nn_failed_phy_trigger_event_fraction(stage1, cmap = 'CMRmap_r', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().phy_trigger_few_neighbor_channel_triggered_count()
    return draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)

def draw_num_channel_triggered_hist(stage1, axis = None, phys_trigger = False, muon_candidate = False, 
                                    zoom = False):
    axis = axis if axis is not None else matplotlib.pyplot.gca()
    mult_hist = stage1.const_charge_stats().const_num_channel_triggered_hist()
    nn_hist = stage1.const_charge_stats().const_num_contiguous_channel_triggered_hist()
    if(phys_trigger):
        mult_hist = stage1.const_charge_stats().const_phy_trigger_num_channel_triggered_hist()
        nn_hist = stage1.const_charge_stats().const_phy_trigger_num_contiguous_channel_triggered_hist()
    elif(muon_candidate):
        mult_hist = stage1.const_charge_stats().const_muon_candidate_num_channel_triggered_hist()
        nn_hist = None
    calin.plotting.plot_histogram(mult_hist, hatch='/', color='C0', lw=2, label='Multiplicity',
         histtype='step', axis=axis, normalise=True)
    if(nn_hist is not None):
        calin.plotting.plot_histogram(nn_hist, hatch='\\\\', color='C1', lw=1, label='Neighbours',
            histtype='step', axis=axis, normalise=True)
    if(zoom):
        axis.set_xlim(-0.5,9.5)
        axis.set_xticks(numpy.arange(0,10))
    else:
        axis.set_xlim(0.25,stage1.run_config().configured_channel_index_size())
        axis.set_xscale('log')
    axis.set_yscale('log')
    axis.legend()
    axis.set_xlabel('Number of channels')
    axis.set_ylabel('Probability')

def draw_mean_wf(stage1, dataset='pedestal', low_gain = False, pedestals = None,
        subtract_pedestal = False, figure_factory = calin.plotting.PyPlotFigureFactory()):

    mwf = None
    if(dataset == 'pedestal'):
        mwf = stage1.const_mean_wf_pedestal() if stage1.has_mean_wf_pedestal() else None
    elif(dataset == 'physics'):
        mwf = stage1.const_mean_wf_physics() if stage1.has_mean_wf_physics() else None
    elif(dataset == 'external_flasher'):
        mwf = stage1.const_mean_wf_external_flasher() if stage1.has_mean_wf_external_flasher() else None
    elif(dataset == 'internal_flasher'):
        mwf = stage1.const_mean_wf_internal_flasher() if stage1.has_mean_wf_internal_flasher() else None
    else:
        raise RuntimeError('Unknown data set : '+dataset)

    if mwf is None:
        return None

    if(low_gain):
        nchan = mwf.channel_low_gain_size()
        chan_nentries = [ mwf.channel_low_gain(i).num_entries() for i in range(nchan) ]
        chan_traces = [ mwf.channel_low_gain(i).mean_waveform() for i in range(nchan) ]
        cam_nentries = mwf.const_camera_low_gain().num_entries() if mwf.has_camera_low_gain() else 0
        cam_trace = mwf.const_camera_low_gain().mean_waveform() if mwf.has_camera_high_gain() else None
    else:
        nchan = mwf.channel_high_gain_size()
        chan_nentries = [ mwf.channel_high_gain(i).num_entries() for i in range(nchan) ]
        chan_traces = [ mwf.channel_high_gain(i).mean_waveform() for i in range(nchan) ]
        cam_nentries = mwf.const_camera_high_gain().num_entries() if mwf.has_camera_high_gain() else 0
        cam_trace = mwf.const_camera_high_gain().mean_waveform() if mwf.has_camera_high_gain() else None

    if(numpy.max(chan_nentries) == 0):
        return None

    DC = dc_units(1, low_gain)

    fig_mwf, axis_mwf = figure_factory.new_histogram_figure()
    has_label = False
    for ichan in range(nchan):
        if(chan_nentries[ichan]):
            wf = chan_traces[ichan]
            if(subtract_pedestal):
                if(pedestals is None):
                    wf -= numpy.mean(wf)
                else:
                    wf -= pedestals[ichan]
            if(has_label):
                axis_mwf.plot(wf,'k',alpha=0.1)
            else:
                axis_mwf.plot(wf,'k',alpha=0.1,label='Channels')
                has_label = True
    if(cam_nentries):
        wf = cam_trace
        if(subtract_pedestal):
            if(pedestals is None):
                wf -= numpy.mean(wf)
            else:
                wf -= numpy.mean(pedestals)
        axis_mwf.plot(wf,'C1',label='Camera average')
    axis_mwf.legend()
    axis_mwf.grid()
    axis_mwf.set_xlabel('Waveform sample number')
    axis_mwf.set_ylabel('Mean waveform amplitude [%s]'%DC)
    axis_mwf.set_title('Mean waveform %s(%s), run : %d'%('offset ' if subtract_pedestal else '',trigger_type_and_gain_title(dataset, low_gain), stage1.run_number()))

    fig_dict = dict()
    fig_dict['waveform_mean_'+dataset+('_lg' if low_gain else '_hg')+('_offset' if subtract_pedestal else '')] = [ fig_mwf, axis_mwf ]
    return fig_dict

def draw_mean_wf_deviation_from_camera_mean(stage1, dataset='pedestal',
        pedestals = None, low_gain=False, axis = None, cmap = 'inferno',
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75, stat_format=None):

    axis = axis if axis is not None else matplotlib.pyplot.gca()

    if(dataset == 'pedestal'):
        mwf = stage1.const_mean_wf_pedestal()
    elif(dataset == 'physics'):
        mwf = stage1.const_mean_wf_physics()
    elif(dataset == 'external_flasher'):
        mwf = stage1.const_mean_wf_external_flasher()
    elif(dataset == 'internal_flasher'):
        mwf = stage1.const_mean_wf_internal_flasher()
    else:
        raise RuntimeError('Unknown data set : '+dataset)

    if(low_gain):
        cwf = mwf.const_camera_low_gain().mean_waveform()
    else:
        cwf = mwf.const_camera_high_gain().mean_waveform()
    if(pedestals is None):
        cwf -= numpy.mean(cwf)
    else:
        cwf -= numpy.mean(pedestals)

    DC = dc_units(1, low_gain)
    if(stat_format is None):
        stat_format = '%.2f ' + DC

    mask = numpy.zeros(mwf.channel_high_gain_size(),dtype=bool)
    chi2 = numpy.zeros(mwf.channel_high_gain_size())
    for ichan in range(mwf.channel_high_gain_size()):
        if(low_gain):
            chan = mwf.const_channel_low_gain(ichan)
        else:
            chan = mwf.const_channel_high_gain(ichan)
        if(chan.num_entries()):
            wf = chan.mean_waveform()
            if(pedestals is None):
                wf -= numpy.mean(wf)
            else:
                wf -= pedestals[ichan]
            mask[ichan] = True
            chi2[ichan] = sum((wf-cwf)**2)
        else:
            chi2[ichan] = numpy.nan
    data = numpy.sqrt(chi2)

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()

    pc = calin.plotting.plot_camera_image(data, cl,
        configured_channels=rc.configured_channel_id(), channel_mask = mask,
        cmap=cmap, draw_outline=draw_outline, pix_lw=pix_lw,
        outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True,
        stats_format=stat_format, stats_fontsize=stat_label_fontsize)

    cb = calin.plotting.add_colorbar_and_clipping(axis, pc, data, mask=mask, percentile=99.5,
            camera_layout=cl, configured_channels=rc.configured_channel_id(),
            cb_label='Waveform RMS from camera mean [%s]'%DC)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

def draw_elapsed_time_hist(stage1, axis = None):
    axis = axis if axis is not None else matplotlib.pyplot.gca()
    ri = stage1.const_run_info()
    so = []
    so.append(calin.plotting.plot_histogram(ri.elapsed_time_histogram(),
        # xright=ri.elapsed_time_histogram().xval_max(),
        color='C0', density=True, lw=2, label='All events', axis=axis))
    if(ri.elapsed_time_histogram_trigger_physics().sum_w()):
        so.append(calin.plotting.plot_histogram(ri.elapsed_time_histogram_trigger_physics(),
            # xright=ri.elapsed_time_histogram().xval_max(),
            color='C1', density=True, label='Physics', axis=axis))
    if(ri.elapsed_time_histogram_trigger_pedestal().sum_w()):
        so.append(calin.plotting.plot_histogram(ri.elapsed_time_histogram_trigger_pedestal(),
            # xright=ri.elapsed_time_histogram().xval_max(),
            color='C2', density=True, label='Pedestal', axis=axis))
    if(ri.elapsed_time_histogram_trigger_external_flasher().sum_w()):
        so.append(calin.plotting.plot_histogram(ri.elapsed_time_histogram_trigger_external_flasher(),
            # xright=ri.elapsed_time_histogram().xval_max(),
            color='C3', density=True, label='External flasher', axis=axis))
    if(ri.elapsed_time_histogram_trigger_internal_flasher().sum_w()):
        so.append(calin.plotting.plot_histogram(ri.elapsed_time_histogram_trigger_internal_flasher(),
            # xright=ri.elapsed_time_histogram().xval_max(),
            color='C4', density=True, label='Internal flasher', axis=axis))

    fraclost = ri.elapsed_time_histogram().overflow_hi()/(ri.elapsed_time_histogram().sum_w()+ri.elapsed_time_histogram().overflow_hi()+ri.elapsed_time_histogram().overflow_lo())
    if(fraclost > 0.01):
        axis.text(0.02, 0.98, 'Warning: %.1f%% of event times exceed histogram limit'%(fraclost*100), 
                  color='k', fontsize=10, ha='left', va='top', transform=axis.transAxes)

    axis.set_ylim([0,axis.get_ylim()[1]*1.15])
    axis.set_xlabel('Elapsed time [s]')
    axis.set_ylabel('Event rate on disk [Hz]')
    axis.legend(loc=4)
    axis.grid()

    return so

def draw_event_number_histogram(stage1, axis = None):
    axis = axis if axis is not None else matplotlib.pyplot.gca()
    ri = stage1.const_run_info()

    so = calin.plotting.plot_histogram(ri.const_event_number_histogram(),lw=2,density=True,
        xleft=1,xright=ri.const_event_number_histogram().xval_max(), axis=axis)

    axis.set_ylim([0,axis.get_ylim()[1]*1.15])
    axis.set_xlabel('Event number')
    axis.set_ylabel('Fraction of events on disk')
    axis.grid()

    return so

def draw_charge_spectrum(stage1, dataset = 'external_flasher', low_gain = False,
        ped = None, pedvarbase = None, evf = 1.2, draw_median=False, draw_scale=False, draw_gain=False,
        figure_factory = calin.plotting.PyPlotFigureFactory(),
        cmap = 'inferno', stat_label_fontsize=4.75, pix_lw = 0, outline_lw = 0.5,
        outline_color = '#888888'):

    if(dataset == 'physics'):
        wfh = stage1.const_wf_hists_physics() if stage1.has_wf_hists_physics() else None
    elif(dataset == 'pedestal'):
        wfh = stage1.const_wf_hists_pedestal() if stage1.has_wf_hists_pedestal() else None
    elif(dataset == 'external_flasher'):
        wfh = stage1.const_wf_hists_external_flasher() if stage1.has_wf_hists_external_flasher() else None
    elif(dataset == 'internal_flasher'):
        wfh = stage1.const_wf_hists_internal_flasher() if stage1.has_wf_hists_internal_flasher() else None
    else:
        raise RuntimeError('Unknown dataset type :',dataset)

    if wfh is None:
        return None

    all_hist = []
    if(low_gain):
        csg = stage1.const_charge_stats().const_low_gain()
        for i in range(wfh.low_gain_channel_size()):
            all_hist.append(wfh.const_low_gain_channel(i))
        cam_hist = wfh.const_low_gain_camera()
        fig_name = dataset+'_lg'
    else:
        csg = stage1.const_charge_stats().const_high_gain()
        for i in range(wfh.high_gain_channel_size()):
            all_hist.append(wfh.const_high_gain_channel(i))
        cam_hist = wfh.const_high_gain_camera()
        fig_name = dataset+'_hg'

    if(len(all_hist)==0 or numpy.max([h.const_opt_win_qsum().sum_w() for h in all_hist])==0):
        return None

    nsamp = stage1.const_config().const_low_gain_opt_sum().integration_n() if \
            (low_gain and stage1.config().has_low_gain_opt_sum()) \
        else stage1.const_config().const_high_gain_opt_sum().integration_n()
    if ped is None:
        ped = calin.diagnostics.stage1_analysis.estimate_run_pedestal(stage1, low_gain)
        ped *= nsamp

    fig_dict = dict()

    fig_hist, axis_hist = figure_factory.new_histogram_figure()
    fig_dict['charge_spectrum_'+fig_name] = [ fig_hist, axis_hist ]

    p_x_lr = 0.005
    all_xl = []
    all_xr = []
    all_xc = []
    all_hist_x = []
    all_hist_y = []

    xl_min = numpy.inf
    xr_max = -numpy.inf
    x_min = numpy.inf
    x_max = -numpy.inf
    has_label = False

    nevent = 0
    nchanevent = 0
    chan_hist_x = []
    chan_hist_y = []

    for i in range(len(all_hist)):
        offset = ped[i]
        if(all_hist[i].has_opt_win_qsum()):
            h = calin.math.histogram.densify(all_hist[i].const_opt_win_qsum())
        else:
            h = calin.math.histogram.densify(all_hist[i].const_sig_win_qsum())
        x = h.xval0() + h.dxval()*numpy.arange(h.bins_size()+1) - offset
        y = numpy.append(0, numpy.cumsum(h.bins()))
        xl,xc,xr = numpy.interp(y[-1]*numpy.asarray([p_x_lr,0.5,1-p_x_lr]),y,x)
        xl_min = min(xl_min, max(1.5*(xl-xc) + xc, numpy.min(x)))
        xr_max = max(xr_max, min(1.5*(xr-xc) + xc, numpy.max(x)))
        x_min = min(x_min, numpy.min(x))
        x_max = max(x_max, numpy.max(x))
        all_xl.append(xl)
        all_xc.append(xc)
        all_xr.append(xr)
        all_hist_x.append(x)
        all_hist_y.append(y)
        nevent = max(nevent, h.sum_w())
        nchanevent += h.sum_w()
        args = dict()
        if(not has_label):
            args['label'] = 'Channels'
            has_label = True
        calin.plotting.plot_histogram(h, xoffset=-offset, histtype='floating',
            color='k', alpha=0.1, density=True, normalise=True, axis=axis_hist, **args)

    h = None
    if(cam_hist.has_opt_win_qsum()):
        h = calin.math.histogram.densify(cam_hist.const_opt_win_qsum())
    elif(cam_hist.has_sig_win_qsum()):
        h = calin.math.histogram.densify(cam_hist.const_sig_win_qsum())
    if(h is not None):
        scale = 1/len(all_hist)
        offset = sum(ped) * scale
        h = calin.math.histogram.rebin(h, int(numpy.sqrt(len(all_hist))))
        calin.plotting.plot_histogram(h, xscale=scale, xoffset=-offset, histtype='floating',
            color='C1', density=True, normalise=True, label='Camera average', axis=axis_hist)
        x = (h.xval0() + h.dxval()*numpy.arange(h.bins_size()+1)) * scale - offset
        y = numpy.append(0, numpy.cumsum(h.bins()))
        xl,xc,xr = numpy.interp(y[-1]*numpy.asarray([p_x_lr,0.5,1-p_x_lr]),y,x)
        xl_min = min(xl_min, max(1.5*(xl-xc) + xc, numpy.min(x)))
        xr_max = max(xr_max, min(1.5*(xr-xc) + xc, numpy.max(x)))
        cam_hist_x = x
        cam_hist_y = y

    xlim_l = xl_min - 0.025*(xr_max-xl_min)
    xlim_r = xr_max + 0.025*(xr_max-xl_min)

    nchanevent_shown = 0
    for i in range(len(all_hist)):
        x = all_hist_x[i]
        y = all_hist_y[i]
        yl,yr = numpy.interp([xlim_l, xlim_r],x,y,left=0,right=y[-1])
        nchanevent_shown += yr-yl

    DC = dc_units(nsamp, low_gain)

    t = 'Number of events : %s'%calin.util.string.int64_to_string_with_commas(int(nevent))
    if(nchanevent_shown < 0.9999*nchanevent):
        t += ' (%.2f%% shown)'%(100*nchanevent_shown/nchanevent)
    elif(nchanevent_shown < nchanevent):
        t += ' (>99.99%% shown)'%(100*nchanevent_shown/nchanevent)
    t += '\nChannel charge range : %s to %s %s'%(
        calin.util.string.double_to_string_with_commas(x_min,1),
        calin.util.string.double_to_string_with_commas(x_max,1), DC)
    if(h is not None):
        t += '\nCamera mean charge range : %s to %s %s'%(
            calin.util.string.double_to_string_with_commas(numpy.min(cam_hist_x),1),
            calin.util.string.double_to_string_with_commas(numpy.max(cam_hist_x),1), DC)

    axis_hist.set_yscale('log')
    axis_hist.set_xlabel(f'Summed waveform [{DC}]')
    axis_hist.set_ylabel(f'Density [1/{DC}]')
    axis_hist.legend(loc=1)
    axis_hist.set_title('Charge spectrum (' + trigger_type_and_gain_title(dataset, low_gain) + '), run: %d'%stage1.run_number())

    ylo,yhi = axis_hist.get_ylim()
    axis_hist.set_ylim(ylo, yhi*3)
    axis_hist.set_xlim(xlim_l, xlim_r)
    axis_hist.text(0.02, 0.98, t, fontsize=6, ha='left', va='top', transform=axis_hist.transAxes)

    if(draw_median):
        fig_median, axis_median = figure_factory.new_camera_figure()
        fig_dict['charge_median_'+fig_name] = [ fig_median, axis_median ]

        all_median = numpy.asarray(all_xc)
        ref_value = numpy.nanmedian(all_median)
        data = all_median/ref_value
        mask = numpy.ones_like(all_median, dtype=bool)
        rc = stage1.const_run_config()
        cl = rc.const_camera_layout()

        pc = calin.plotting.plot_camera_image(data, cl, channel_mask=mask,
                        configured_channels=rc.configured_channel_id(),
                        axis=axis_median, cmap=cmap, draw_outline=True, draw_stats=True,
                        pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_channels=True, stats_format='%.3f',
                        stats_fontsize=stat_label_fontsize, draw_top12_val=True)

        cb = calin.plotting.add_colorbar_and_clipping(axis_median, pc, data, mask=mask, percentile=99.5,
                camera_layout=cl, configured_channels=rc.configured_channel_id(),
                cb_label='Median charge relative to reference')

        max_xy = cl.camera_boundary_maxabs_xy()
        axis_median.text(max_xy,max_xy,'Reference : %.1f %s'%(ref_value, DC),
                fontsize=stat_label_fontsize, fontfamily='monospace',
                ha='right', va='top')

        # cb = axis_median.get_figure().colorbar(pc, ax=axis_median, label='Median relative signal')

        axis_median.get_xaxis().set_visible(False)
        axis_median.get_yaxis().set_visible(False)
        axis_median.set_title('Charge spectrum median (' + trigger_type_and_gain_title(dataset, low_gain) + '), run: %d'%stage1.run_number())

    if(draw_scale):
        fig_scale, axis_scale = figure_factory.new_camera_figure()
        fig_dict['charge_scale_'+fig_name] = [ fig_scale, axis_scale ]

        all_rms = (numpy.asarray(all_xr)-numpy.asarray(all_xl))/(-2*scipy.special.erfinv(0.005*2-1)*numpy.sqrt(2))
        data = all_rms
        mask = numpy.ones_like(all_rms, dtype=bool)
        rc = stage1.const_run_config()
        cl = rc.const_camera_layout()

        pc = calin.plotting.plot_camera_image(data, cl, channel_mask=mask,
                        configured_channels=rc.configured_channel_id(),
                        axis=axis_scale, cmap=cmap, draw_outline=True, draw_stats=True,
                        pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_channels=True, stats_format='%.3f',
                        stats_fontsize=stat_label_fontsize, draw_top12_val=True)

        cb = calin.plotting.add_colorbar_and_clipping(axis_scale, pc, data, mask=mask, percentile=99.5,
                camera_layout=cl, configured_channels=rc.configured_channel_id(),
                cb_label='Charge spectrum scale [%s]'%DC)

        # cb = axis_scale.get_figure().colorbar(pc, ax=axis_scale, label='Median relative signal')

        axis_scale.get_xaxis().set_visible(False)
        axis_scale.get_yaxis().set_visible(False)
        axis_scale.set_title('Charge spectrum scale (' + trigger_type_and_gain_title(dataset, low_gain) + '), run: %d'%stage1.run_number())

    if(draw_gain):
        fig_gain, axis_gain = figure_factory.new_camera_figure()
        fig_dict['charge_gain_'+fig_name] = [ fig_gain, axis_gain ]

        rc = stage1.const_run_config()
        cl = rc.const_camera_layout()
        all_rms = (numpy.asarray(all_xr)-numpy.asarray(all_xl))/(-2*scipy.special.erfinv(0.005*2-1)*numpy.sqrt(2))
        mask = numpy.ones_like(all_rms, dtype=bool)
        all_var = all_rms**2
        if pedvarbase is not None:
            all_var -= pedvarbase
        all_gain = all_var/numpy.asarray(all_xc)/evf

        pc = calin.plotting.plot_camera_image(all_gain, cl, channel_mask=mask,
                        configured_channels=rc.configured_channel_id(),
                        axis=axis_gain, cmap=cmap, draw_outline=True, draw_stats=True,
                        pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_channels=True, stats_format='%.2f '+DC,
                        stats_fontsize=stat_label_fontsize, draw_top12_val=True)

        cb = calin.plotting.add_colorbar_and_clipping(axis_gain, pc, all_gain, mask=mask, percentile=99.5,
                camera_layout=cl, configured_channels=rc.configured_channel_id(),
                cb_label='Gain estimate [%s/PE]'%DC)

        # cb = axis_gain.get_figure().colorbar(pc, ax=axis_gain, label='Gain estimate [%s]'%DC)

        max_xy = cl.camera_boundary_maxabs_xy()
        axis_gain.text(max_xy,max_xy,'EVF=$1+\\beta^2$ : %.3f'%numpy.median(evf),
                fontsize=stat_label_fontsize, fontfamily='monospace',
                ha='right', va='top')

        axis_gain.get_xaxis().set_visible(False)
        axis_gain.get_yaxis().set_visible(False)
        axis_gain.set_title('Photostatistics gain (' + trigger_type_and_gain_title(dataset, low_gain) + '), run: %d'%stage1.run_number())

        fig_intensity, axis_intensity = figure_factory.new_camera_figure()
        fig_dict['charge_intensity_'+fig_name] = [ fig_intensity, axis_intensity ]

        all_intensity = numpy.asarray(all_xc)/all_gain

        pc = calin.plotting.plot_camera_image(all_intensity, cl, channel_mask=mask,
                        configured_channels=rc.configured_channel_id(),
                        axis=axis_intensity, cmap=cmap, draw_outline=True, draw_stats=True,
                        pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_channels=True, stats_format='%.3f PE',
                        stats_fontsize=stat_label_fontsize, draw_top12_val=True)

        cb = calin.plotting.add_colorbar_and_clipping(axis_intensity, pc, all_intensity, mask=mask, percentile=99.5,
                camera_layout=cl, configured_channels=rc.configured_channel_id(),
                cb_label='Approximate intensity [PE/flash]')

        # cb = axis_intensity.get_figure().colorbar(pc, ax=axis_intensity, label='Intensity estimate [PE]')

        axis_intensity.get_xaxis().set_visible(False)
        axis_intensity.get_yaxis().set_visible(False)
        axis_intensity.set_title('Estimated flash intensity (' + trigger_type_and_gain_title(dataset, low_gain) + '), run: %d'%stage1.run_number())

    return fig_dict

def draw_high_gain_low_gain(stage1, dataset='max_sample', subtract_pedestal=False, draw_P1 = False,
    figure_factory = calin.plotting.PyPlotFigureFactory()):

    if(not stage1.has_charge_stats() or not stage1.const_charge_stats().has_dual_gain() or \
            stage1.const_charge_stats().const_dual_gain().all_max_sample_count_size() == 0):
        return None, None, None

    dg = stage1.const_charge_stats().const_dual_gain()

    all_h_c = []
    all_h_m = []
    all_h_v = []
    if(dataset == 'max_sample'):
        hg_scale = 1.0
        lg_scale = 1.0
        dataset_label = 'maximum sample'
        figure_name = 'high_low_max_sample'
        figure_title = 'High-gain vs low-gain max sample'
        xcut_base = 30
        ycut_base = 4000
        for ichan in range(dg.all_opt_sum_count_size()):
            all_h_c.append(dg.all_max_sample_count(ichan))
            all_h_m.append(dg.all_max_sample_mean(ichan))
            all_h_v.append(dg.all_max_sample_var(ichan))
        DChg = dc_units(1)
        DClg = dc_units(1, True)
    elif(dataset == 'opt_sum'):
        conf_hg = stage1.const_config().const_high_gain_opt_sum()
        conf_lg = stage1.const_config().const_low_gain_opt_sum() if stage1.const_config().has_low_gain_opt_sum() else stage1.const_config().const_high_gain_opt_sum()
        hg_scale = 1.0/float(conf_hg.integration_n())
        lg_scale = 1.0/float(conf_lg.integration_n())
        dataset_label = '%d-sample average'%conf_hg.integration_n()
        figure_name = 'high_low_window_sum'
        figure_title = 'High-gain vs low-gain window sum'
        xcut_base = 8
        ycut_base = 4000
        for ichan in range(dg.all_opt_sum_count_size()):
            all_h_c.append(dg.all_opt_sum_count(ichan))
            all_h_m.append(dg.all_opt_sum_mean(ichan))
            all_h_v.append(dg.all_opt_sum_var(ichan))
        DChg = dc_units(conf_hg.integration_n())
        DClg = dc_units(conf_lg.integration_n(), True)
    else:
        raise RuntimeError('Unknown dataset : '+dataset)

    if(numpy.max([h.sum_w() for h in all_h_c]) == 0):
        return None, None, None

    fig_dict = dict()

    fig_hist, axis_hist = figure_factory.new_histogram_figure()
    fig_dict[figure_name] = [ fig_hist, axis_hist ]

    all_hg_ped = calin.diagnostics.stage1_analysis.estimate_run_pedestal(stage1, low_gain=False)
    all_lg_ped = calin.diagnostics.stage1_analysis.estimate_run_pedestal(stage1, low_gain=True)

    if subtract_pedestal and \
            numpy.all(numpy.bitwise_and(numpy.isnan(all_hg_ped),numpy.isnan(all_lg_ped))):
        return None, None, None

    all_P0 = []
    all_P1 = []
    fit_mask = []

    minx = numpy.inf
    maxx = -numpy.inf
    for h_c, h_m, h_v, hg_ped, lg_ped in zip(all_h_c, all_h_m, all_h_v, all_hg_ped, all_lg_ped):
        h_c = calin.math.histogram.densify(h_c)
        h_m = calin.math.histogram.densify(h_m)
        h_v = calin.math.histogram.densify(h_v)

        m = h_c.bins()>0
        if numpy.count_nonzero(m) == 0:
            all_P0.append(numpy.nan)
            all_P1.append(numpy.nan)
            fit_mask.append(False)
            continue

        x = h_c.xval0() + (numpy.arange(h_c.bins_size())+0.5)*h_c.dxval()

        vindex = int((h_m.xval0()-h_c.xval0())/h_v.dxval())
        y = x*0
        y[vindex:vindex+h_m.bins_size()] = h_m.bins()

        vindex = int((h_v.xval0()-h_c.xval0())/h_v.dxval())
        dy = y*0
        if(h_v.bins_size() > 0):
            dy[vindex:vindex+h_v.bins_size()] = numpy.sqrt(numpy.max(h_v.bins(),0))
        dy = dy/(numpy.sqrt(h_c.bins())+1e-9)

        xcut = xcut_base
        ycut = ycut_base
        if(subtract_pedestal):
            if(numpy.isnan(lg_ped) or numpy.isnan(hg_ped)):
                all_P0.append(numpy.nan)
                all_P1.append(numpy.nan)
                fit_mask.append(False)
                continue

            x = lg_scale * x - lg_ped
            y = hg_scale * y - hg_ped
            ycut -= hg_ped
        else:
            x = lg_scale * x
            y = hg_scale * y
            xcut += lg_ped
        dy = hg_scale * dy

        mm = numpy.bitwise_and(numpy.bitwise_and(h_c.bins()>10, dy>0),
            numpy.bitwise_and(x>xcut , y<ycut))

        if(numpy.count_nonzero(mm) > 5):
            try:
                P = numpy.polyfit(x[mm], y[mm], 1, w=dy[mm])
                all_P0.append(P[0])
                all_P1.append(P[1])
                fit_mask.append(True)
            except:
                all_P0.append(numpy.nan)
                all_P1.append(numpy.nan)
                fit_mask.append(False)
        else:
            all_P0.append(numpy.nan)
            all_P1.append(numpy.nan)
            fit_mask.append(False)

        minx = min(minx, numpy.min(x[m]))
        maxx = max(maxx, numpy.max(x[m]))

        if(len(x) > 0 and numpy.count_nonzero(m) > 0):
            axis_hist.plot(x[m], y[m], 'k.', markersize=2, alpha=0.2)


    axis_hist.set_xlim(minx-10, min(maxx+10,700))
    axis_hist.set_xlabel('Low-gain %s [%s]'%(dataset_label,DClg))
    axis_hist.set_ylabel('Mean high-gain %s [%s]'%(dataset_label,DChg))
    axis_hist.grid()
    axis_hist.set_title(figure_title+', run : %d'%stage1.run_number())

    mask = ~numpy.isnan(all_P0)
    if(numpy.count_nonzero(mask)):
        fig_P0, axis_P0 = figure_factory.new_camera_figure()
        fig_dict[figure_name+'_coefficient'] = [ fig_P0, axis_P0 ]

        cl = stage1.const_run_config().const_camera_layout()
        ccid = stage1.const_run_config().configured_channel_id()
        pc = calin.plotting.plot_camera_image(all_P0, cl, channel_mask=mask,
            configured_channels=ccid, cmap='inferno', axis=axis_P0,
            draw_outline=True, hatch_missing_channels=True)
        cb = axis_P0.get_figure().colorbar(pc, ax=axis_P0, label='High/Low ratio')
        # calin.plotting.add_colorbar_and_clipping(gca(), pc, asarray(all_P0), camera_layout=cl,
        #             configured_channels=s1.run_config().configured_channel_id(),
        #             percentile=99.5, percentile_factor=2.0, cb_label='High/Low gain')
        calin.plotting.add_stats(axis_P0, cl.camera_boundary_maxabs_xy(), all_P0,
             ccid, mask=mask, draw_top12_val=True)

        axis_P0.get_xaxis().set_visible(False)
        axis_P0.get_yaxis().set_visible(False)
        axis_P0.set_title(figure_title+' coefficient, run : %d'%stage1.run_number())

    if(draw_P1):
        mask = ~numpy.isnan(all_P1)
        if(numpy.count_nonzero(mask)):
            fig_P1, axis_P1 = figure_factory.new_camera_figure()
            fig_dict[figure_name+'_intercept'] = [ fig_P1, axis_P1 ]

            cl = stage1.const_run_config().const_camera_layout()
            ccid = stage1.const_run_config().configured_channel_id()
            pc = calin.plotting.plot_camera_image(all_P1, cl, channel_mask=mask,
                configured_channels=ccid, cmap='inferno', axis=axis_P1,
                draw_outline=True, hatch_missing_channels=True)
            cb = axis_P1.get_figure().colorbar(pc, ax=axis_P1, label='High/Low intercept')
            # calin.plotting.add_colorbar_and_clipping(gca(), pc, asarray(all_P0), camera_layout=cl,
            #             configured_channels=s1.run_config().configured_channel_id(),
            #             percentile=99.5, percentile_factor=2.0, cb_label='High/Low gain')
            calin.plotting.add_stats(axis_P1, cl.camera_boundary_maxabs_xy(), all_P1,
                 ccid, mask=mask, draw_top12_val=True)

            axis_P1.get_xaxis().set_visible(False)
            axis_P1.get_yaxis().set_visible(False)
            axis_P0.set_title(figure_title+' intercept, run : %d'%stage1.run_number())

    return fig_dict, all_P0, all_P1

def draw_psd(stage1, dataset='all', low_gain=False, draw_camera_plots = True,
        min_peak_freq=200, figure_factory = calin.plotting.PyPlotFigureFactory(),
        cmap = 'inferno', stat_label_fontsize=4.75, pix_lw = 0, outline_lw = 0.5,
        outline_color = '#888888'):

    rc = stage1.const_run_config().Clone()
    cl = rc.const_camera_layout()
    ccid = rc.configured_channel_id()

    dataset_psd = None
    if(dataset == 'physics'):
        dataset_psd = stage1.const_psd_wf_physics().Clone() if stage1.has_psd_wf_physics() else None
    elif(dataset == 'pedestal'):
        dataset_psd = stage1.const_psd_wf_pedestal().Clone() if stage1.has_psd_wf_pedestal() else None
    elif(dataset == 'external_flasher'):
        dataset_psd = stage1.const_psd_wf_external_flasher().Clone() if stage1.has_psd_wf_external_flasher() else None
    elif(dataset == 'internal_flasher'):
        dataset_psd = stage1.const_psd_wf_internal_flasher().Clone() if stage1.has_psd_wf_internal_flasher() else None
    elif(dataset == 'all'):
        all_dataset_psd = []
        if stage1.has_psd_wf_physics(): all_dataset_psd.append(stage1.const_psd_wf_physics())
        if stage1.has_psd_wf_pedestal(): all_dataset_psd.append(stage1.const_psd_wf_pedestal())
        if stage1.has_psd_wf_external_flasher(): all_dataset_psd.append(stage1.const_psd_wf_external_flasher())
        if stage1.has_psd_wf_internal_flasher(): all_dataset_psd.append(stage1.const_psd_wf_internal_flasher())
        for ipsd in all_dataset_psd:
            if dataset_psd is None:
                dataset_psd = ipsd.Clone()
            else:
                dataset_psd.IntegrateFrom(ipsd)
    else:
        raise RuntimeError('Unknown dataset type :',dataset)

    if dataset_psd is None:
        return None

    nchan = dataset_psd.low_gain_size() if low_gain else dataset_psd.high_gain_size()
    if nchan==0:
        return None

    max_num_entries = numpy.max([dataset_psd.low_gain(i).num_entries() if low_gain \
        else dataset_psd.high_gain(i).num_entries() for i in range(nchan)])
    if(max_num_entries < 10):
        return None

    nsamp = rc.num_samples()
    nfreq = calin.math.fftw_util.hcvec_num_real(nsamp)
    freq = (rc.nominal_sampling_frequency() if rc.nominal_sampling_frequency() else 1000) \
            * numpy.arange(nfreq)/rc.num_samples()

    fig_psd = None
    axis_psd = None

    sum_psd_sum = numpy.zeros_like(freq)
    sum_psd_count = 0

    all_lf_var = numpy.zeros(nchan)
    all_hf_var = numpy.zeros(nchan)
    all_pk_var = numpy.zeros(nchan)

    for ichan in range(nchan):
        chan_psd = dataset_psd.const_low_gain(ichan) if low_gain else dataset_psd.const_high_gain(ichan)
        if(chan_psd.num_entries()):
            plot_opt_dict = dict()
            if(fig_psd is None):
                fig_psd, axis_psd = figure_factory.new_histogram_figure()
                plot_opt_dict['label'] = 'Channels'
            psd = chan_psd.psd_sum()/nsamp**2/chan_psd.num_entries()
            sum_psd_sum += psd
            sum_psd_count += 1
            axis_psd.plot(freq[1:], psd[1:],'k',alpha=0.2,**plot_opt_dict)

            all_lf_var[ichan] = psd[0] - (chan_psd.dc_sum()/nsamp/chan_psd.num_entries())**2
            all_hf_var[ichan] = numpy.sum(psd[1:])
            all_pk_var[ichan] = numpy.max(psd[freq>=min_peak_freq])
        else:
            all_lf_var[ichan] = numpy.nan
            all_hf_var[ichan] = numpy.nan
            all_pk_var[ichan] = numpy.nan

    if(fig_psd is None):
        return None

    DC = dc_units(1, low_gain)
    axis_psd.plot(freq[1:], sum_psd_sum[1:]/sum_psd_count,'C1', label='Camera average')
    axis_psd.grid()
    axis_psd.set_xlabel('Frequency [MHz]')
    axis_psd.set_ylabel('Power [%s$^2$]'%DC)
    axis_psd.legend()
    axis_psd.set_title('Mean power spectrum (%s), run : %d'%(trigger_type_and_gain_title(dataset, low_gain), stage1.run_number()))

    fig_dict = dict()
    fig_dict['psd_'+dataset+('_lg' if low_gain else '_hg')] = [ fig_psd, axis_psd ]

    all_lf_var = numpy.asarray(all_lf_var)
    all_hf_var = numpy.asarray(all_hf_var)
    all_pk_var = numpy.asarray(all_pk_var)

    mask = ~numpy.isnan(all_lf_var)
    if(draw_camera_plots and numpy.count_nonzero(mask)):
        fig_lf, axis_lf = figure_factory.new_camera_figure()
        fig_dict['psd_lf_rms_'+dataset+('_lg' if low_gain else '_hg')] = [ fig_lf, axis_lf ]
        data = numpy.sqrt(all_lf_var)
        pc = calin.plotting.plot_camera_image(data, cl, channel_mask=mask,
            configured_channels=ccid, cmap='inferno', axis=axis_lf,
            draw_outline=True, hatch_missing_channels=True)
        calin.plotting.add_colorbar_and_clipping(axis_lf, pc, data, camera_layout=cl,
                    configured_channels=ccid,
                    percentile=99.5, percentile_factor=2.0, cb_label='RMS [%s]'%DC)
        calin.plotting.add_stats(axis_lf, cl.camera_boundary_maxabs_xy(), data,
             ccid, mask=mask, draw_top12_val=True)
        axis_lf.get_xaxis().set_visible(False)
        axis_lf.get_yaxis().set_visible(False)
        axis_lf.set_title('Low-frequency RMS (%s), run : %d'%(trigger_type_and_gain_title(dataset, low_gain), stage1.run_number()))


        fig_hf, axis_hf = figure_factory.new_camera_figure()
        fig_dict['psd_hf_rms_'+dataset+('_lg' if low_gain else '_hg')] = [ fig_hf, axis_hf ]
        data = numpy.sqrt(all_hf_var)
        pc = calin.plotting.plot_camera_image(data, cl, channel_mask=mask,
            configured_channels=ccid, cmap='inferno', axis=axis_hf,
            draw_outline=True, hatch_missing_channels=True)
        calin.plotting.add_colorbar_and_clipping(axis_hf, pc, data, camera_layout=cl,
                    configured_channels=ccid,
                    percentile=99.5, percentile_factor=2.0, cb_label='RMS [%s]'%DC)
        calin.plotting.add_stats(axis_hf, cl.camera_boundary_maxabs_xy(), data,
             ccid, mask=mask, draw_top12_val=True)
        axis_hf.get_xaxis().set_visible(False)
        axis_hf.get_yaxis().set_visible(False)
        axis_hf.set_title('High-frequency RMS (%s), run : %d'%(trigger_type_and_gain_title(dataset, low_gain), stage1.run_number()))


        fig_pk, axis_pk = figure_factory.new_camera_figure()
        fig_dict['psd_peak_'+dataset+('_lg' if low_gain else '_hg')] = [ fig_pk, axis_pk ]
        data = all_pk_var
        pc = calin.plotting.plot_camera_image(data, cl, channel_mask=mask,
            configured_channels=ccid, cmap='inferno', axis=axis_pk,
            draw_outline=True, hatch_missing_channels=True)
        calin.plotting.add_colorbar_and_clipping(axis_pk, pc, data, camera_layout=cl,
                    configured_channels=ccid,
                    percentile=99.5, percentile_factor=2.0, cb_label='Power [%s$^2$]'%DC)
        calin.plotting.add_stats(axis_pk, cl.camera_boundary_maxabs_xy(), data,
             ccid, mask=mask, draw_top12_val=True)
        axis_pk.get_xaxis().set_visible(False)
        axis_pk.get_yaxis().set_visible(False)
        axis_pk.set_title(u'Peak power \u2265%gMHz (%s), run : %d'%(min_peak_freq,
            trigger_type_and_gain_title(dataset, low_gain), stage1.run_number()))

    return fig_dict

def draw_trigger_threshold(stage1, draw_camera_plots = True, 
        waveform_sum = False, ped = None, 
        do_mle = True, nmin = 100, mle_pmin = 1e-200,
        figure_factory = calin.plotting.PyPlotFigureFactory(),
        cmap = 'inferno', stat_label_fontsize=4.75, pix_lw = 0, outline_lw = 0.5,
        outline_color = '#888888'):

    ch_set = stage1.const_wf_hists_l0_trigger_bit_set()
    ch_clr = stage1.const_wf_hists_l0_trigger_bit_clear()

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    ri = stage1.const_run_info()

    if ped is None:
        ped = calin.diagnostics.stage1_analysis.estimate_run_pedestal(stage1)

    if(waveform_sum):
        nsamp = stage1.const_config().const_high_gain_opt_sum().integration_n()
        dataset = '%d-sample sum'%nsamp
        filename = 'sum'
        ped *= nsamp
        aligned_range = 500
        DC = dc_units(nsamp)
    else:
        dataset = 'waveform amplitude'
        filename = 'amplitude'
        aligned_range = 25
        DC = dc_units(1)

    allfit = []
    nchan = ch_set.high_gain_channel_size()
    for ichan in range(nchan):
        if(waveform_sum):
            hset = calin.math.histogram.densify(ch_set.const_high_gain_channel(ichan).opt_win_qsum())
            hclr = calin.math.histogram.densify(ch_clr.const_high_gain_channel(ichan).opt_win_qsum())
        else:
            hset = calin.math.histogram.densify(ch_set.const_high_gain_channel(ichan).full_wf_max())
            hclr = calin.math.histogram.densify(ch_clr.const_high_gain_channel(ichan).full_wf_max())
        allfit.append(calin.diagnostics.stage1_analysis.analyze_trigger_thresholds(
            hset,hclr,ped[ichan],do_mle=do_mle,nmin=nmin,mle_pmin=mle_pmin))

    cal_med = [fitres[5] for fitres in allfit]
    cal_iqr = [fitres[6] for fitres in allfit]

    if(numpy.count_nonzero(~numpy.isnan(cal_med)) < nchan//10):
        return None

    fig_dict = dict()

    fig_frac, axis_frac = figure_factory.new_camera_figure()
    for fitres in allfit:
        axis_frac.plot(fitres[0],fitres[1],'k',alpha=0.1)
    axis_frac.set_xlabel('%s [%s]'%(dataset.capitalize(),DC))
    axis_frac.set_ylabel('Trigger efficiency')
    axis_frac.grid()
    axis_frac.set_title('Trigger efficiency vs %s, run : %d'%(dataset, stage1.run_number()))

    fig_dict['trigger_l0_efficiency_'+filename] = [ fig_frac, axis_frac ]

    fig_frac_aligned, axis_frac_aligned = figure_factory.new_camera_figure()
    for fitres in allfit:
        if not numpy.isnan(fitres[5]):
            axis_frac_aligned.plot(fitres[0]-fitres[5],fitres[1],'k',alpha=0.1)
    axis_frac_aligned.set_xlabel('%s offset [%s]'%(dataset.capitalize(),DC))
    axis_frac_aligned.set_ylabel('Trigger fraction')
    axis_frac_aligned.set_xlim(-aligned_range,aligned_range)
    axis_frac_aligned.grid()
    axis_frac_aligned.set_title('Trigger efficiency vs %s offset, run : %d'%(dataset, stage1.run_number()))

    fig_dict['trigger_l0_efficiency_aligned_'+filename] = [ fig_frac_aligned, axis_frac_aligned ]

    if(draw_camera_plots):
        fig_threshold, axis_threshold = figure_factory.new_camera_figure()
        pc = calin.plotting.plot_camera_image(cal_med, cl,
                            channel_mask=~numpy.isnan(cal_med), axis=axis_threshold,
                            cmap=cmap, draw_outline=True, draw_stats=True, hatch_missing_channels=True,
                            outline_lw=outline_lw, outline_color=outline_color,
                            stats_fontsize=stat_label_fontsize, draw_top12_val=False,
                            stats_format='%.1f ' + DC,
                            draw_top12=True)
        fig_threshold.colorbar(pc, label='Trigger threshold as %s [%s]'%(dataset,DC))
        axis_threshold.get_xaxis().set_visible(False)
        axis_threshold.get_yaxis().set_visible(False)
        axis_threshold.set_title('Trigger threshold as %s, run : %d'%(dataset, stage1.run_number()))
       
        fig_dict['trigger_l0_threshold_'+filename] = [ fig_threshold, axis_threshold ]

        fig_switchon, axis_switchon = figure_factory.new_camera_figure()
        cal_iqr = [fitres[6] for fitres in allfit]
        pc = calin.plotting.plot_camera_image(cal_iqr, cl,
                            channel_mask=~numpy.isnan(cal_med), axis=axis_switchon,
                            cmap=cmap, draw_outline=True, draw_stats=True, hatch_missing_channels=True,
                            outline_lw=outline_lw, outline_color=outline_color,
                            stats_fontsize=stat_label_fontsize, draw_top12_val=False,
                            stats_format='%.1f ' + DC,
                            draw_top12=True)
        fig_switchon.colorbar(pc, label='Trigger switch-on IQR as %s [%s]'%(dataset,DC))
        axis_switchon.get_xaxis().set_visible(False)
        axis_switchon.get_yaxis().set_visible(False)
        axis_switchon.set_title('Trigger switch-on IQR as %s, run : %d'%(dataset, stage1.run_number()))

        fig_dict['trigger_l0_switchon_'+filename] = [ fig_switchon, axis_switchon ]

    return fig_dict
