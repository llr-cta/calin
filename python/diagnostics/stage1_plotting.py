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
        nsamp = stage1.config().low_gain_opt_sum().integration_n() if low_gain \
            else stage1.config().high_gain_opt_sum().integration_n()
        values = charge_stats.all_trigger_ped_win_mean()/nsamp
    else:
        nsamp = rc.num_samples()
        values = charge_stats.ped_trigger_full_wf_mean()/nsamp

    pc = calin.plotting.plot_camera_image(
        values, cl, cmap=cmap,
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
        nsamp = stage1.config().low_gain_opt_sum().integration_n() if low_gain \
            else stage1.config().high_gain_opt_sum().integration_n()
        values = charge_stats.all_trigger_ped_win_var()
    else:
        nsamp = rc.num_samples()
        values = charge_stats.ped_trigger_full_wf_var()

    pc = calin.plotting.plot_camera_image(
        numpy.sqrt(values), cl, cmap=cmap,
        configured_channels=rc.configured_channel_id(),
        draw_outline=draw_outline, pix_lw=pix_lw, outline_lw=outline_lw, outline_color=outline_color,
        axis=axis, hatch_missing_channels=True, draw_stats=True, stats_format='%.2f DC')

    cb = axis.get_figure().colorbar(pc, ax=axis, label='Pedestal %d-sample RMS [DC]'%nsamp)

    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    return pc

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

    vmin = 0.5/evts_ondisk
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

def draw_nectarcam_feb_temperatures(stage1, temperature_set=1, cmap = 'inferno', axis=None,
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):
    tfeb1 = []
    tfeb2 = []
    for modid in stage1.run_config().configured_module_id() :
        if(stage1.nectarcam().ancillary_data().feb_temperature_has_key(int(modid))):
            measurement_set = stage1.nectarcam().ancillary_data().feb_temperature(int(modid))
            tfeb1.append(numpy.mean([measurement_set.measurement(i).tfeb1() for i in range(measurement_set.measurement_size())]))
            tfeb2.append(numpy.mean([measurement_set.measurement(i).tfeb2() for i in range(measurement_set.measurement_size())]))
    tfeb = tfeb1 if temperature_set==1 else tfeb2

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    pc = calin.plotting.plot_camera_module_image(tfeb, stage1.run_config().camera_layout(),
                    configured_modules=stage1.run_config().configured_module_id(),
                    axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
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

    tfeb1 = []
    tfeb2 = []
    for modid in stage1.run_config().configured_module_id() :
        if(stage1.nectarcam().ancillary_data().feb_temperature_has_key(int(modid))):
            measurement_set = stage1.nectarcam().ancillary_data().feb_temperature(int(modid))
            tfeb1.append(minmax([measurement_set.measurement(i).tfeb1() for i in range(measurement_set.measurement_size())]))
            tfeb2.append(minmax([measurement_set.measurement(i).tfeb2() for i in range(measurement_set.measurement_size())]))
    tfeb = tfeb1 if temperature_set==1 else tfeb2

    if(axis is None):
        axis = matplotlib.pyplot.gca()

    pc = calin.plotting.plot_camera_module_image(tfeb, stage1.run_config().camera_layout(),
                    configured_modules=stage1.run_config().configured_module_id(),
                    axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
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
        axis_freq = None, axis_t0 = None, axis_chi2 = None, clockid=0, cmap = 'inferno',
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, stat_label_fontsize=4.75):

    freq_offset_ppm, time_offset_ns, d2_per_event, n_problem_bins = \
        calin.diagnostics.stage1_analysis.summarize_module_clock_regression(stage1, clockid)

    if(axis_freq is not None):
        axis = axis_freq
        data = freq_offset_ppm
        pc = calin.plotting.plot_camera_module_image(data, stage1.run_config().camera_layout(),
                        configured_modules=stage1.run_config().configured_module_id(),
                        axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                        mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_modules=True, stats_format=u'%.2f ppm',
                        stats_fontsize=stat_label_fontsize)
        cb = axis.get_figure().colorbar(pc, label='FEB oscillator frequency offset [PPM]')

        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                    configured_modules=stage1.run_config().configured_module_id(),
                                    pc=pc, module_values = data,
                                    fontsize=mod_label_fontsize)

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

    if(axis_t0 is not None):
        axis = axis_t0
        data = time_offset_ns
        pc = calin.plotting.plot_camera_module_image(data, stage1.run_config().camera_layout(),
                        configured_modules=stage1.run_config().configured_module_id(),
                        axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                        mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_modules=True, stats_format=u'%.2f ns',
                        stats_fontsize=stat_label_fontsize)
        cb = axis.get_figure().colorbar(pc, label='UCTS time at FEB clock reset [ns]')

        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                    configured_modules=stage1.run_config().configured_module_id(),
                                    pc=pc, module_values = data,
                                    fontsize=mod_label_fontsize)

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

    if(axis_chi2 is not None):
        axis = axis_chi2
        data = d2_per_event
        pc = calin.plotting.plot_camera_module_image(data, stage1.run_config().camera_layout(),
                        configured_modules=stage1.run_config().configured_module_id(),
                        axis=axis, cmap=cmap, draw_outline=True, draw_stats=True,
                        mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
                        hatch_missing_modules=True, stats_format=u'%.2f ns^2',
                        stats_fontsize=stat_label_fontsize)
        cb = axis.get_figure().colorbar(pc, label='Linrear-fit residual per event [ns$^2$]')

        if(mod_label_fontsize is not None and mod_label_fontsize>0):
            calin.plotting.add_module_numbers(axis, stage1.run_config().camera_layout(),
                                    configured_modules=stage1.run_config().configured_module_id(),
                                    pc=pc, module_values = data,
                                    fontsize=mod_label_fontsize)

        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)

    return pc
