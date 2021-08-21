# calin/python/plotting.py -- Stephen Fegan -- 2017-03-27
#
# Plotting of camera faces and histograms
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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
import matplotlib.path
import matplotlib.patches
import colorsys
import numpy
import copy

import calin.math.histogram
import calin.ix.iact_data.instrument_layout
import calin.iact_data.instrument_layout
import calin.math.regular_grid

def obsolete_plot_camera(pix_data, camera_layout, configured_channels = None, ax_in = None,
        cbar_label = None):
    if(configured_channels is not None and len(configured_channels) != len(pix_data)):
        raise ValueError('configured_channels must either be None or have same length as pix_data')
    pix = []
    max_xy = 0
    for pix_index in range(len(pix_data)):
        pix_id = int(configured_channels[pix_index]) if configured_channels is not None else  pix_index
        vx = camera_layout.channel(pix_id).outline_polygon_vertex_x_view()
        vy = camera_layout.channel(pix_id).outline_polygon_vertex_y_view()
        vv = numpy.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
        pix.append(matplotlib.pyplot.Polygon(vv,closed=True))
    ax = ax_in if ax_in is not None else matplotlib.pyplot.gca()
    pc = matplotlib.collections.PatchCollection(pix, cmap=matplotlib.cm.jet)
    pc.set_array(numpy.asarray(pix_data))
    pc.set_linewidths(0)
    ax.add_collection(pc)
    ax.axis('square')
    ax.axis(numpy.asarray([-1,1,-1,1])*1.05*max_xy)
    ax.set_xlabel('X coordinate [cm]')
    ax.set_ylabel('Y coordinate [cm]')
    cbar = matplotlib.pyplot.colorbar(pc, ax=ax)
    if cbar_label is not None:
        cbar.set_label(cbar_label)
    return pc

def obsolete_plot_module_camera(mod_data, camera_layout, configured_modules = None, ax_in = None,
        cbar_label = None, module_label_color = None, module_label_font_size = 8):
    if(configured_modules is not None and len(configured_modules) != len(mod_data)):
        raise ValueError('configured_modules must either be None or have same length as mod_data')
    mod = []
    max_xy = 0
    for mod_index in range(len(mod_data)):
        mod_id = int(configured_modules[mod_index]) if configured_modules is not None else mod_index
        vx = camera_layout.module(mod_id).outline_polygon_vertex_x_view()
        vy = camera_layout.module(mod_id).outline_polygon_vertex_y_view()
        vv = numpy.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
        mod.append(matplotlib.pyplot.Polygon(vv,closed=True))
    ax = ax_in if ax_in is not None else matplotlib.pyplot.gca()
    pc = matplotlib.collections.PatchCollection(mod, cmap=matplotlib.cm.jet)
    pc.set_array(numpy.asarray(mod_data))
    pc.set_linewidths(0)
    ax.add_collection(pc)
    if(module_label_color is not None):
        for mod_index in range(len(mod_data)):
            mod_id = int(configured_modules[mod_index]) if configured_modules is not None else mod_index
            ax.text(camera_layout.module(mod_id).x(), camera_layout.module(mod_id).y(),
                '%d'%mod_id, ha='center', va='center', fontsize=module_label_font_size,
                color=module_label_color)
    ax.axis('square')
    ax.axis(numpy.asarray([-1,1,-1,1])*1.05*max_xy)
    ax.set_xlabel('X coordinate [cm]')
    ax.set_ylabel('Y coordinate [cm]')
    cbar = matplotlib.pyplot.colorbar(pc, ax=ax)
    if cbar_label is not None:
        cbar.set_label(cbar_label)
    return pc

################################### OBSOLETE ###################################

def layout_to_polygon_vxy(layout, plate_scale = 1.0, rotation = 0.0):
    crot = numpy.cos(rotation/180.0*numpy.pi)
    srot = numpy.sin(rotation/180.0*numpy.pi)
    all_vxy = []
    ibegin = 0;
    for iend in layout.outline_polygon_vertex_index():
        vx = layout.outline_polygon_vertex_x()[ibegin:iend]*plate_scale
        vy = layout.outline_polygon_vertex_y()[ibegin:iend]*plate_scale
        vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
        all_vxy.append(numpy.column_stack([vx, vy]))
        ibegin=iend
    return all_vxy

def layout_to_polygon(layout, plate_scale = 1.0, rotation = 0.0, **args):
    all_p = []
    for vxy in layout_to_polygon_vxy(layout, plate_scale, rotation):
        all_p.append(matplotlib.pyplot.Polygon(vxy, **args))
    if len(all_p) == 1:
        return all_p[0]
    return all_p

################################################################################

def layout_to_path(layout, plate_scale = 1.0, rotation = 0.0):
    crot = numpy.cos(rotation/180.0*numpy.pi)
    srot = numpy.sin(rotation/180.0*numpy.pi)
    vx = layout.outline_polygon_vertex_x()*plate_scale
    vy = layout.outline_polygon_vertex_y()*plate_scale
    vx, vy = vx*crot - vy*srot, vy*crot + vx*srot

    nv = nv = layout.outline_polygon_vertex_x_size()+layout.outline_polygon_vertex_index_size()
    v = numpy.zeros((nv,2))
    c = numpy.zeros(nv,dtype=numpy.uint8)

    iv = 0
    for ic in range(layout.outline_polygon_vertex_index_size()):
        nv = layout.outline_polygon_vertex_index(ic) - iv
        v[iv+ic:iv+ic+nv,0] = vx[iv:iv+nv]
        v[iv+ic:iv+ic+nv,1] = vy[iv:iv+nv]
        c[iv+ic]            = matplotlib.path.Path.MOVETO
        c[iv+ic+1:iv+ic+nv] = matplotlib.path.Path.LINETO
        c[iv+ic+nv]         = matplotlib.path.Path.CLOSEPOLY
        iv = layout.outline_polygon_vertex_index(ic)

    return matplotlib.path.Path(v,c)

def layout_to_path_patch(layout, plate_scale = 1.0, rotation = 0.0, **args):
    p = layout_to_path(layout, plate_scale, rotation)
    return matplotlib.patches.PathPatch(p, **args)

def add_outline(axis, layout, plate_scale = 1.0, rotation = 0.0, fill=False,
        outline_lw = 0.5, outline_ls = '-', outline_color = '#888888', **args):
    return axis.add_patch(layout_to_path_patch(layout, plate_scale, rotation,
            fill=fill, lw=outline_lw, ls=outline_ls, edgecolor=outline_color, **args))

def add_stats(axis, max_xy, values, ids, mask=None, stats_fontsize=4.75, stats_format='%.3f',
        draw_top12 = True):

    if(mask is not None):
        if(len(mask) != len(values)):
            raise ValueError('Mask must either be None or have same length as values')

        mask = asarray(mask, dtype=bool)
        values = numpy.asarray(values)[mask]
        ids = numpyt.asarray(ids)[mask]

    if(len(values) == 0):
        return

    axis.text(-max_xy,max_xy,('Median : '+stats_format+u'\nScale : '+stats_format+
                '\nMin : '+stats_format+'\nMax : '+stats_format)%(
            numpy.median(values), 1.4826*numpy.median(abs(values - numpy.median(values))),
            numpy.min(values), numpy.max(values)),
        fontsize=stats_fontsize, fontfamily='monospace',
        ha='left', va='top')
    if(draw_top12):
        nshow = min(len(values),12)
        isort = numpy.argsort(values)
        id_len = numpy.max(list(map(lambda i: len(str(i)), ids)))
        id_fmt = '%%0%dd'%id_len
        topN = 'Top %d\n'%nshow
        for i in range(nshow):
            if(i):
                topN += '\n' if(i%4 == 0) else ' '
            topN += id_fmt%ids[isort[-(i+1)]]
        axis.text(-max_xy,-max_xy,topN,
            fontsize=stats_fontsize, fontfamily='monospace',
            ha='left', va='bottom')
        botN = 'Bottom %d\n'%nshow
        for i in range(nshow):
            if(i):
                botN += '\n' if(i%4 == 0) else ' '
            botN += id_fmt%ids[isort[i]]
        axis.text(max_xy,-max_xy,botN,
            fontsize=stats_fontsize, fontfamily='monospace',
            ha='right', va='bottom')

def add_module_numbers(axis, camera_layout, configured_modules = None, dx = 0, dy = 0,
        plate_scale = 1.0, rotation = 0.0, pc = None, module_values = None,
        brightness_hi_threshold=0.45, brightness_hi_color = 'k', brightness_lo_color = 'w',
        nonconfigured_color = None, **args):

    if(configured_modules is None):
        configured_modules = list(range(camera_layout.module_size()))

    if(module_values is not None and len(module_values) != len(configured_modules)):
        raise ValueError('module_values must be Null or have same size as configured_modules')

    modules_labelled = numpy.zeros(camera_layout.module_size(), dtype=bool)

    def rgba_to_brightness(rgba):
        return colorsys.rgb_to_hls(*rgba[:3])[1]

    crot = numpy.cos(rotation/180.0*numpy.pi)
    srot = numpy.sin(rotation/180.0*numpy.pi)
    for i, modid in enumerate(configured_modules):
        m = camera_layout.module(int(modid))
        label_color = brightness_hi_color
        if(pc is not None and module_values is not None):
            if(module_values[i] is numpy.nan):
                label_color = brightness_hi_color
            else:
                rgba = pc.cmap(pc.norm(module_values[i]))
                brightness = rgba_to_brightness(rgba)
                if(brightness < brightness_hi_threshold):
                    label_color = brightness_lo_color
        x = (m.x()*crot - m.y()*srot)*plate_scale
        y = (m.y()*crot + m.x()*srot)*plate_scale
        axis.text(x+dx, y+dy, '%d'%modid, ha='center', va='center',
            color=label_color, **args)
        modules_labelled[modid] = True

    if(nonconfigured_color is not None):
        for modid, is_labelled in enumerate(modules_labelled):
            if(not is_labelled):
                x = (m.x()*crot - m.y()*srot)*plate_scale
                y = (m.y()*crot + m.x()*srot)*plate_scale
                axis.text(x+dx, y+dy, '%d'%modid, ha='center', va='center',
                    color=nonconfigured_color, **args)

def add_colorbar_and_clipping(axis, pc, data, mask=None,
        camera_layout=None, plate_scale = 1.0, rotation = 0.0, configured_channels=None,
        cb_label=None, percentile=100, percentile_factor=2.0,
        under_color=None, over_color=None, clip_highlight='#0044ff', clip_highlight_lw=1.5):
    axis = axis if axis is not None else matplotlib.pyplot.gca()
    mask = mask if mask is not None else ones_like(data, dtype=bool)
    dmax = numpy.max(data[mask])
    dmin = numpy.min(data[mask])
    dmed = numpy.median(data[mask])
    dpch = numpy.percentile(data[mask],min(100,percentile))
    dpcl = numpy.percentile(data[mask],max(0,100-percentile))
    liml = max(dmin, dmed + (dpcl-dmed)*percentile_factor)
    limh = min(dmax, dmed + (dpch-dmed)*percentile_factor)
    pc.cmap = copy.copy(pc.cmap)
    pc.set_clim(liml,limh)
    cb_clip = 'neither'
    if(liml>dmin):
        if(under_color is not None):
            pc.cmap.set_under(under_color)
        if(limh<dmax):
            if(over_color is not None):
                pc.cmap.set_over(over_color)
            cb_clip = 'both'
        else:
            cb_clip = 'min'
    elif(limh<dmax):
        if(over_color is not None):
            pc.cmap.set_over(over_color)
        cb_clip = 'max'

    cb = axis.get_figure().colorbar(pc, ax=axis, label=cb_label, extend=cb_clip)

    if(cb_clip != 'neither' and camera_layout is not None and clip_highlight is not None):
        if(configured_channels is None):
            configured_channels = numpy.arange(camera_layout.channel_size())
        clipped_mask = numpy.bitwise_and(mask, numpy.bitwise_or(data>limh, data<liml))
        clipped_camera_layout = calin.iact_data.instrument_layout.reorder_camera_channels(
            camera_layout, configured_channels[clipped_mask])
        add_outline(axis,clipped_camera_layout, plate_scale=plate_scale, rotation=rotation,
            outline_color=clip_highlight, outline_lw=clip_highlight_lw)

    return cb

def plot_camera_image(channel_data, camera_layout, channel_mask = None,
        configured_channels = None, zero_suppression = None,
        plate_scale = None, rotation = 0.0, R = None,
        cmap = matplotlib.cm.CMRmap_r, axis = None, draw_outline = False,
        pix_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        hatch_missing_channels = False,
        draw_stats = False, stats_fontsize = 4.75, stats_format='%.3f', draw_top12 = True,
        additional_polygons = [], additional_polygon_data = []):
    if(channel_mask is None and zero_suppression is not None):
        channel_mask = numpy.asarray(channel_data)>zero_suppression
    if(channel_mask is not None and len(channel_mask) != len(channel_data)):
        raise ValueError('channel_mask must either be None or have same length as channel_data')
    if(configured_channels is not None and len(configured_channels) != len(channel_data)):
        raise ValueError('configured_channels must either be None or have same length as pix_data')
    poly = additional_polygons.copy()
    poly_data = additional_polygon_data.copy()
    pix_data = []
    pix_ids = []
    pix_gridid = []
    if plate_scale is None:
        plate_scale = 1
        if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
            plate_scale = 180/numpy.pi/camera_layout.effective_focal_length()
    crot = numpy.cos(rotation/180.0*numpy.pi)
    srot = numpy.sin(rotation/180.0*numpy.pi)
    if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
        camera_layout = camera_layout.camera()
    max_xy = max(numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_x())),
        numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_y())))*plate_scale

    for chan_index in range(len(channel_data)):
        chan_id = int(configured_channels[chan_index]) if configured_channels is not None else chan_index
        chan = camera_layout.channel(chan_id)
        if(chan.pixel_index() != -1):
            pix_gridid.append(chan.pixel_grid_index())
        if((channel_mask is None or channel_mask[chan_index]) and chan.pixel_index() != -1):
            vx = chan.outline_polygon_vertex_x_view()*plate_scale
            vy = chan.outline_polygon_vertex_y_view()*plate_scale
            vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
            vv = numpy.column_stack([vx, vy])
            poly.append(matplotlib.pyplot.Polygon(vv,closed=True))
            poly_data.append(channel_data[chan_index])
            pix_data.append(channel_data[chan_index])
            pix_ids.append(chan_id)

    axis = axis or matplotlib.pyplot.gca()

    cl_conf = camera_layout
    if draw_outline and hatch_missing_channels and configured_channels is not None and \
            len(configured_channels) != camera_layout.channel_size():
        cl_conf = calin.iact_data.instrument_layout.reorder_camera_channels(
            camera_layout, numpy.asarray(configured_channels,dtype=numpy.int32))
        add_outline(axis, camera_layout, hatch='XXXXX',
            plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color, zorder=-2)
        add_outline(axis, cl_conf,
            plate_scale=plate_scale, rotation=rotation, fill=True, facecolor='w',
            outline_lw=outline_lw, outline_color=outline_color, zorder=-1)

    pc = matplotlib.collections.PatchCollection(poly, cmap=cmap)
    pc.set_array(numpy.asarray(poly_data))
    pc.set_linewidths(pix_lw)
    axis.add_collection(pc)

    if draw_outline:
        add_outline(axis, camera_layout, plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    if draw_outline and hatch_missing_channels and configured_channels is not None and \
            len(configured_channels) != camera_layout.channel_size():
        calin.plotting.add_outline(axis, cl_conf,
            plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    if(draw_stats):
        add_stats(axis, max_xy, pix_data, pix_ids,
           stats_fontsize=stats_fontsize, stats_format=stats_format,
           draw_top12=draw_top12)

    axis.axis('square')
    axis.axis(numpy.asarray([-1,1,-1,1])*(R or 1.05*max_xy))
    return pc

def plot_camera_module_image(module_data, camera_layout, module_mask = None,
        configured_modules = None, zero_suppression = None,
        plate_scale = None, rotation = 0.0, R = None,
        cmap = matplotlib.cm.CMRmap_r, axis = None, draw_outline = False,
        mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        hatch_missing_modules = False,
        draw_stats = False, stats_fontsize = 4.75, stats_format='%.3f', draw_top12 = True,
        additional_polygons = [], additional_polygon_data = []):
    if(module_mask is None and zero_suppression is not None):
        module_mask = numpy.asarray(module_data)>zero_suppression
    if(module_mask is not None and len(module_mask) != len(module_data)):
        raise ValueError('module_mask must either be None or have same length as module_data')
    if(configured_modules is not None and len(configured_modules) != len(module_data)):
        raise ValueError('configured_modules must either be None or have same length as pix_data')
    poly = additional_polygons.copy()
    poly_data = additional_polygon_data.copy()
    mod_data = []
    mod_ids = []
    max_xy = 0
    if plate_scale is None:
        plate_scale = 1
        if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
            plate_scale = 180/numpy.pi/camera_layout.effective_focal_length()
    crot = numpy.cos(rotation/180.0*numpy.pi)
    srot = numpy.sin(rotation/180.0*numpy.pi)
    if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
        camera_layout = camera_layout.camera()
    max_xy = max(numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_x())),
        numpy.max(numpy.abs(camera_layout.outline_polygon_vertex_y())))*plate_scale

    for mod_index in range(len(module_data)):
        mod_id = int(configured_modules[mod_index]) if configured_modules is not None else mod_index
        chan = camera_layout.module(mod_id)
        if((module_mask is None or module_mask[mod_index]) and chan.module_index() != -1):
            vx = chan.outline_polygon_vertex_x_view()*plate_scale
            vy = chan.outline_polygon_vertex_y_view()*plate_scale
            vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
            vv = numpy.column_stack([vx, vy])
            poly.append(matplotlib.patches.Polygon(vv,closed=True))
            poly_data.append(module_data[mod_index])
            mod_data.append(module_data[mod_index])
            mod_ids.append(mod_id)

    axis = axis or matplotlib.pyplot.gca()

    cl_conf = camera_layout
    if draw_outline and hatch_missing_modules and configured_modules is not None and \
            len(configured_modules) != camera_layout.module_size():
        cl_conf = calin.iact_data.instrument_layout.reorder_camera_modules(
            camera_layout, numpy.array(configured_modules, dtype=numpy.int32))
        add_outline(axis, camera_layout, hatch='XXXXX',
            plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color, zorder=-2)
        add_outline(axis, cl_conf,
            plate_scale=plate_scale, rotation=rotation, fill=True, facecolor='w',
            outline_lw=outline_lw, outline_color=outline_color, zorder=-1)

    pc = matplotlib.collections.PatchCollection(poly, cmap=cmap)
    pc.set_array(numpy.asarray(poly_data))
    pc.set_linewidths(mod_lw)
    axis.add_collection(pc)

    if draw_outline:
        add_outline(axis, camera_layout, plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    if draw_outline and hatch_missing_modules and configured_modules is not None and \
            len(configured_modules) != camera_layout.module_size():
        calin.plotting.add_outline(axis, cl_conf,
            plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    if(draw_stats):
        add_stats(axis, max_xy, mod_data, mod_ids,
            stats_fontsize=stats_fontsize, stats_format=stats_format,
            draw_top12=draw_top12)

    axis.axis('square')
    axis.axis(numpy.asarray([-1,1,-1,1])*(R or 1.05*max_xy))
    return pc

def plot_histogram(h, density = False, normalise = False,
        weight_scale_function=None,
        xscale = 1, xoffset = 0, yscale = 1, yoffset = 0,
        xleft = -numpy.inf, xright = numpy.inf,
        xscale_as_log10 = False, draw_poisson_errors = False,
        histtype='steps', ecolor=None, axis = None, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
        hdy = numpy.sqrt(h.all_weight())
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        if(h.sparse_bins_size() > 0):
            h = calin.math.histogram.densify(h)
        hx = h.xval0()+h.dxval()*numpy.arange(0,h.bins_size())
        hy = h.bins()
        hdy = numpy.sqrt(h.bins())
    else:
        raise Exception('Unknown histogram type: '+str(type(h)))

    if(len(hx) == 0):
        return

    if(weight_scale_function is not None):
        hy = weight_scale_function(hy)
        hdy = weight_scale_function(hdy) # This probably doesn't make sense

    hx = numpy.append(hx, hx[-1]+h.dxval()) * xscale + xoffset

    if(hx[0] < xleft):
        if(hx[1] < xleft):
            raise ValueError("xleft outside of first bin")
        hx[0] = xleft
    if(hx[-1] > xright):
        if(hx[-2] > xright):
            raise ValueError("xright outside of final bin")
        hx[-1] = xright

    if density:
        hy /= numpy.abs(numpy.diff(hx))
        hdy /= numpy.abs(numpy.diff(hx))
    if normalise:
        hy /= h.sum_w()
        hdy /= h.sum_w()

    hy = numpy.append(hy, hy[-1]) * yscale + yoffset
    hdy = hdy * yscale + yoffset
    if(xscale_as_log10):
        hx = 10**hx
    axis = axis or matplotlib.pyplot.gca()

    histtype='step' if histtype=='steps' else histtype
    histtype='stepfilled' if histtype=='stepsfilled' else histtype
    histtype='bar' if histtype=='bars' else histtype

    if(histtype == 'floating'):
        so = axis.step(hx,hy, *args, where='post', **nargs)

        if(draw_poisson_errors):
            so1 = axis.vlines(0.5*(hx[:-1]+hx[1:]), hy[:-1]-hdy, hy[:-1]+hdy)
            so1.set_linestyles(so[0].get_linestyle())
            so1.set_color(so[0].get_color() if ecolor is None else ecolor)
            so1.set_linewidth(so[0].get_linewidth())
            so1.set_zorder(so[0].get_zorder())
            so = [ *so, so1 ]
    else:
        _, _, so = axis.hist(hx[:-1],hx,weights=hy[:-1],histtype=histtype,**nargs)

        if(draw_poisson_errors):
            so1 = axis.vlines(0.5*(hx[:-1]+hx[1:]), hy[:-1]-hdy, hy[:-1]+hdy)
            so1.set_linestyles(so[0].get_linestyle())
            so1.set_color(so[0].get_edgecolor() if ecolor is None else ecolor)
            so1.set_linewidth(so[0].get_linewidth())
            so1.set_zorder(so[0].get_zorder())
            so = [ *so, so1 ]

    if(xscale_as_log10):
        axis.set_xscale('log')

    return so

def plot_histogram_cumulative(h, plot_as_cdf = False, plot_as_cmf = False,
        right_to_left = False,
        xscale = 1, xoffset = 0, yscale = 1, yoffset = 0,
        xscale_as_log10 = False, axis = None, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        if(h.sparse_bins_size() > 0):
            h = calin.math.histogram.densify(h)
        hx = h.xval0()+h.dxval()*numpy.arange(0,h.bins_size())
        hy = h.bins()
    else:
        raise Exception('Unknown histogram type: '+str(type(h)))
    if plot_as_cdf or plot_as_cmf:
        hy /= h.sum_w()
    hx = numpy.append(hx, hx[-1]+h.dxval()) * xscale + xoffset
    hy = numpy.append(0, hy)
    if right_to_left:
        hy = numpy.sum(hy)-numpy.cumsum(hy)
    else:
        hy = numpy.cumsum(hy)
    hy =  hy * yscale + yoffset
    if(xscale_as_log10):
        hx = 10**hx
    axis = axis or matplotlib.pyplot.gca()
    so = axis.plot(hx, hy, *args, **nargs)
    if(xscale_as_log10):
        so[0].axes.set_xscale('log')
    return so
