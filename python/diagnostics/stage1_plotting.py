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

import calin.plotting
import calin.diagnostics.stage1
import calin.diagnostics.stage1_analysis
import calin.iact_data.instrument_layout

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
        nsamp = stage1.config().low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.config().has_low_gain_opt_sum()) \
            else stage1.config().high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_mean()/nsamp
    else:
        nsamp = rc.num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_mean()/nsamp

    pc = calin.plotting.plot_camera_image(
        values, cl, channel_mask=nevent>0, cmap=cmap,
        configured_channels=rc.configured_channel_id(),
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True, stats_format='%.2f DC')

    cb = axis.get_figure().colorbar(pc, ax=axis, label='Pedestal mean [DC]')

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
        nsamp = stage1.config().low_gain_opt_sum().integration_n() if \
                (low_gain and stage1.config().has_low_gain_opt_sum()) \
            else stage1.config().high_gain_opt_sum().integration_n()
        nevent = charge_stats.all_trigger_event_count()
        values = charge_stats.all_trigger_ped_win_var()
    else:
        nsamp = rc.num_samples()
        nevent = charge_stats.ped_trigger_event_count()
        values = charge_stats.ped_trigger_full_wf_var()

    data = numpy.sqrt(values)
    mask = nevent>0

    pc = calin.plotting.plot_camera_image(
        data, cl, channel_mask=mask, cmap=cmap,
        configured_channels=rc.configured_channel_id(),
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True, stats_format='%.2f DC')

    cb = calin.plotting.add_colorbar_and_clipping(axis, pc, data, mask=mask, percentile=99.5,
            camera_layout=cl, configured_channels=rc.configured_channel_id(),
            cb_label='Pedestal %d-sample RMS [DC]'%nsamp)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc, cb

def draw_missing_components_fraction(stage1, cmap = 'CMRmap_r', axis = None,
        draw_outline=True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, aux_label_fontsize=5.5, stat_label_fontsize=4.75):
    SQRT3_2 = numpy.sqrt(3)/2
    RAD_30 = 30.0/180.0*numpy.pi

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    ri = stage1.const_run_info()

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
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad,max_xy-rad*SQRT3_2),numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad,max_xy-rad*SQRT3_2-2*rad),numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(-max_xy+rad+2*rad*SQRT3_2,max_xy-rad*SQRT3_2-rad),numVertices=6, radius=rad, orientation=RAD_30))
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
            if(type(pp) is matplotlib.patches.RegularPolygon):
                xy = pp.xy
                p = matplotlib.patches.RegularPolygon(xy, pp.numvertices, pp.radius, pp.orientation)
            else:
                xy = pp.center
                p = matplotlib.patches.Circle(xy, pp.radius)
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

    cb = axis.get_figure().colorbar(pc, ax=axis, label='Component missing fraction')
    cb.ax.plot([-max_xy, max_xy], [1.0/evts_ondisk, 1.0/evts_ondisk], 'g-', lw=0.75)

    return pc

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

def draw_nectarcam_feb_temperatures(stage1, temperature_set=1, cmap = 'coolwarm', axis=None,
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):
    tfeb = []
    mask = []
    for modid in stage1.run_config().configured_module_id() :
        if(stage1.nectarcam().ancillary_data().feb_temperature_has_key(int(modid))):
            measurement_set = stage1.nectarcam().ancillary_data().feb_temperature(int(modid))
            mt = numpy.asarray([measurement_set.measurement(i).tfeb1() for i in range(measurement_set.measurement_size())]) if temperature_set==1 \
                else numpy.asarray([measurement_set.measurement(i).tfeb2() for i in range(measurement_set.measurement_size())])
            mmask = mt != 0
            if(numpy.count_nonzero(mmask)):
                tfeb.append(numpy.mean(mt[mmask]))
                mask.append(True)
            else:
                tfeb.append(numpy.nan)
                mask.append(False)
        else:
            tfeb.append(numpy.nan)
            mask.append(False)

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    pc = calin.plotting.plot_camera_module_image(tfeb, stage1.run_config().camera_layout(),
                    configured_modules=stage1.run_config().configured_module_id(),
                    module_mask=mask, axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                    mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                    hatch_missing_modules=True, stats_format=u'%4.2f\u00b0C',
                    stats_fontsize=stat_label_fontsize)
    cb = axis.get_figure().colorbar(pc, label=u'FEB temperature %d [\u2103]'%temperature_set)

    if(mod_label_fontsize is not None and mod_label_fontsize>0):
        calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                configured_modules=stage1.run_config().configured_module_id(),
                                pc=pc, module_values = tfeb, fontsize=mod_label_fontsize)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    def measurement_val(m):
        m = numpy.asarray(m)
        return ' ====' if(numpy.all(m==0) or numpy.any(m<-30) or numpy.any(m>50)) else \
            ' %4.1f'%numpy.mean(m)

    if(stage1.nectarcam().ancillary_data().has_ecc_measurements() and \
            stage1.nectarcam().ancillary_data().ecc_measurements().measurement_size()>0):
        measurement_set = stage1.nectarcam().ancillary_data().ecc_measurements()

        tecc = 'ECC'
        tecc += measurement_val([measurement_set.measurement(i).temp_01() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_02() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_03() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_04() for i in range(measurement_set.measurement_size())])
        tecc += u'\n\u00b0C'
        tecc += measurement_val([measurement_set.measurement(i).temp_05() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_06() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_07() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_08() for i in range(measurement_set.measurement_size())])
        tecc += '\n'
        tecc += measurement_val([measurement_set.measurement(i).temp_09() for i in range(measurement_set.measurement_size())])[1:]
        tecc += measurement_val([measurement_set.measurement(i).temp_10() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_11() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_12() for i in range(measurement_set.measurement_size())])

        tecc4 = measurement_val([measurement_set.measurement(i).temp_13() for i in range(measurement_set.measurement_size())])[1:]
        tecc4 += measurement_val([measurement_set.measurement(i).temp_14() for i in range(measurement_set.measurement_size())])
        tecc4 += measurement_val([measurement_set.measurement(i).temp_15() for i in range(measurement_set.measurement_size())])
        tecc4 += measurement_val([measurement_set.measurement(i).temp_16() for i in range(measurement_set.measurement_size())])

        if(tecc4 != '==== ==== ==== ===='):
            tecc += '\n' + tecc4

        max_xy = max(numpy.max(numpy.abs(stage1.run_config().camera_layout().outline_polygon_vertex_x())),
                          numpy.max(numpy.abs(stage1.run_config().camera_layout().outline_polygon_vertex_y())))
        axis.text(max_xy,max_xy,tecc,ha='right',va='top',fontfamily='monospace',fontsize=stat_label_fontsize)

    return pc

def draw_nectarcam_feb_temperatures_minmax(stage1, temperature_set=1, cmap = 'inferno', axis = None,
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):

    def minmax(x):
        return numpy.max(x)-numpy.min(x)

    tfeb = []
    mask = []
    for modid in stage1.run_config().configured_module_id() :
        if(stage1.nectarcam().ancillary_data().feb_temperature_has_key(int(modid))):
            measurement_set = stage1.nectarcam().ancillary_data().feb_temperature(int(modid))
            mt = numpy.asarray([measurement_set.measurement(i).tfeb1() for i in range(measurement_set.measurement_size())]) if temperature_set==1 \
                else numpy.asarray([measurement_set.measurement(i).tfeb2() for i in range(measurement_set.measurement_size())])
            mmask = mt != 0
            if(numpy.count_nonzero(mmask)>1):
                tfeb.append(minmax(mt[mmask]))
                mask.append(True)
            else:
                tfeb.append(numpy.nan)
                mask.append(False)
        else:
            tfeb.append(numpy.nan)
            mask.append(False)

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    pc = calin.plotting.plot_camera_module_image(tfeb, stage1.run_config().camera_layout(),
                    configured_modules=stage1.run_config().configured_module_id(),
                    module_mask=mask, axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                    mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                    hatch_missing_modules=True, stats_format=u'%4.2f\u00b0C',
                    stats_fontsize=stat_label_fontsize)
    cb = axis.get_figure().colorbar(pc, label=u'FEB temperature %d range (max-min) [\u2103]'%temperature_set)

    if(mod_label_fontsize is not None and mod_label_fontsize>0):
        calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                configured_modules=stage1.run_config().configured_module_id(),
                                pc=pc, module_values = tfeb, fontsize=mod_label_fontsize)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    def measurement_val(m):
        m = numpy.asarray(m)
        return ' ====' if(numpy.all(m==0) or numpy.any(m<-30) or numpy.any(m>50)) else \
            ' %4.1f'%minmax(m)

    if(stage1.nectarcam().ancillary_data().has_ecc_measurements() and \
            stage1.nectarcam().ancillary_data().ecc_measurements().measurement_size()>0):
        measurement_set = stage1.nectarcam().ancillary_data().ecc_measurements()

        tecc = 'ECC'
        tecc += measurement_val([measurement_set.measurement(i).temp_01() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_02() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_03() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_04() for i in range(measurement_set.measurement_size())])
        tecc += '\n\u00b0C'
        tecc += measurement_val([measurement_set.measurement(i).temp_05() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_06() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_07() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_08() for i in range(measurement_set.measurement_size())])
        tecc += '\n'
        tecc += measurement_val([measurement_set.measurement(i).temp_09() for i in range(measurement_set.measurement_size())])[1:]
        tecc += measurement_val([measurement_set.measurement(i).temp_10() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_11() for i in range(measurement_set.measurement_size())])
        tecc += measurement_val([measurement_set.measurement(i).temp_12() for i in range(measurement_set.measurement_size())])

        tecc4 = measurement_val([measurement_set.measurement(i).temp_13() for i in range(measurement_set.measurement_size())])[1:]
        tecc4 += measurement_val([measurement_set.measurement(i).temp_14() for i in range(measurement_set.measurement_size())])
        tecc4 += measurement_val([measurement_set.measurement(i).temp_15() for i in range(measurement_set.measurement_size())])
        tecc4 += measurement_val([measurement_set.measurement(i).temp_16() for i in range(measurement_set.measurement_size())])

        if(tecc4 != '==== ==== ==== ===='):
            tecc += '\n' + tecc4

        max_xy = max(numpy.max(numpy.abs(stage1.run_config().camera_layout().outline_polygon_vertex_x())),
                          numpy.max(numpy.abs(stage1.run_config().camera_layout().outline_polygon_vertex_y())))
        axis.text(max_xy,max_xy,tecc,ha='right',va='top',fontfamily='monospace',fontsize=stat_label_fontsize)

    return pc

def draw_all_clock_regression(stage1,
        axis_freq = None, axis_t0 = None, axis_chi2 = None, axis_freq_spread = None, axis_t0_spread = None,
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
                        stats_fontsize=stat_label_fontsize)
        cb = axis.get_figure().colorbar(pc, label=cb_label)

        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                    configured_modules=stage1.run_config().configured_module_id(),
                                    pc=pc, module_values=mod_data,
                                    fontsize=mod_label_fontsize)

        if(stat_label_fontsize is not None and stat_label_fontsize>0):
            fmt_string = 'UCTS 10MHz : %s\nTIB 10MHz : %s\nFEB sum : %s'%(stats_format,stats_format,stats_format)
            axis.text(max_xy,max_xy,fmt_string%(cam_data[0], cam_data[1], cam_data[2]),
                fontsize=stat_label_fontsize, fontfamily='monospace',
                ha='right', va='top')

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

        return pc

    all_pc = []
    if(axis_freq is not None):
        all_pc.append(draw_it(axis_freq, freq_offset_ppm, cam_freq_offset_ppm, '%.2f ppm',
            'Mean oscillator frequency offset [ppm]'))

    if(axis_t0 is not None):
        all_pc.append(draw_it(axis_t0, time_offset_ns, cam_time_offset_ns, '%.2f ns',
            'Mean UCTS time at counter reset [ns]'))

    if(axis_chi2 is not None):
        all_pc.append(draw_it(axis_chi2, numpy.sqrt(d2_per_event), numpy.sqrt(cam_d2_per_event), '%.2f ns',
            'Linear-fit RMS residual per event [ns]'))

    if(axis_freq_spread is not None):
        all_pc.append(draw_it(axis_freq_spread, freq_spread_ppm, cam_freq_spread_ppm, '%.3f ppm',
            'Drift in oscillator frequency [ppm]'))

    if(axis_t0_spread is not None):
        all_pc.append(draw_it(axis_t0_spread, time_spread_ns, cam_time_spread_ns, '%.3f ns',
            'Drift in UCTS time at counter reset [ns]'))

    return all_pc

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

def draw_nn_failed_phy_trigger_event_fraction(stage1, cmap = 'CMRmap_r', axis = None,
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75):
    channel_count = stage1.const_charge_stats().phy_trigger_few_neighbor_channel_triggered_count()
    return draw_channel_event_fraction(stage1, channel_count, cb_label='Event fraction',
        log_scale=False, cmap=cmap, axis=axis,
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        stat_label_fontsize=stat_label_fontsize)

def draw_num_channel_triggered_hist(stage1, axis = None, phys_trigger = False, zoom = False):
    axis = axis if axis is not None else matplotlib.pyplot.gca()
    mult_hist = stage1.const_charge_stats().const_phy_trigger_num_channel_triggered_hist() \
        if phys_trigger else stage1.const_charge_stats().const_num_channel_triggered_hist()
    nn_hist = stage1.const_charge_stats().const_phy_trigger_num_contiguous_channel_triggered_hist() \
        if phys_trigger else stage1.const_charge_stats().const_num_contiguous_channel_triggered_hist()
    calin.plotting.plot_histogram(mult_hist, hatch='/', color='C0', lw=2, label='Multiplicity',
         histtype='step', axis=axis, normalise=True)
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
        subtract_pedestal = False, axis = None):
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

    has_label = False
    for ichan in range(mwf.channel_high_gain_size()):
        if(low_gain):
            chan = mwf.const_channel_low_gain(ichan)
        else:
            chan = mwf.const_channel_high_gain(ichan)
        if(chan.num_entries()):
            wf = chan.mean_waveform()
            if(subtract_pedestal):
                if(pedestals is None):
                    wf -= numpy.mean(wf)
                else:
                    wf -= pedestals[ichan]
            if(has_label):
                axis.plot(wf,'k',alpha=0.1)
            else:
                axis.plot(wf,'k',alpha=0.1,label='channel')
                has_label = True
    if(low_gain):
        cam = mwf.const_camera_low_gain()
    else:
        cam = mwf.const_camera_high_gain()
    if(cam.num_entries()):
        wf = cam.mean_waveform()
        if(subtract_pedestal):
            if(pedestals is None):
                wf -= numpy.mean(wf)
            else:
                wf -= numpy.mean(pedestals)
        axis.plot(wf,'C1',label='Camera average')
    axis.legend()
    axis.grid()
    axis.set_xlabel('Waveform sample number')
    axis.set_ylabel('Mean waveform amplitude [DC]')

def draw_mean_wf_deviation_from_camera_mean(stage1, dataset='pedestal',
        pedestals = None, low_gain=False, axis = None, cmap = 'inferno',
        draw_outline=True, pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        stat_label_fontsize=4.75, stat_format='%.2f DC'):

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
            cb_label='Waveform RMS from camera mean [DC]')

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
