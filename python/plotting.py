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

import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np

import calin.math.histogram
import calin.ix.iact_data.instrument_layout
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
        vv = np.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
        pix.append(plt.Polygon(vv,closed=True))
    ax = ax_in if ax_in is not None else plt.gca()
    pc = matplotlib.collections.PatchCollection(pix, cmap=matplotlib.cm.jet)
    pc.set_array(np.asarray(pix_data))
    pc.set_linewidths(0)
    ax.add_collection(pc)
    ax.axis('square')
    ax.axis(np.asarray([-1,1,-1,1])*1.05*max_xy)
    ax.set_xlabel('X coordinate [cm]')
    ax.set_ylabel('Y coordinate [cm]')
    cbar = plt.colorbar(pc, ax=ax)
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
        vv = np.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
        mod.append(plt.Polygon(vv,closed=True))
    ax = ax_in if ax_in is not None else plt.gca()
    pc = matplotlib.collections.PatchCollection(mod, cmap=matplotlib.cm.jet)
    pc.set_array(np.asarray(mod_data))
    pc.set_linewidths(0)
    ax.add_collection(pc)
    if(module_label_color is not None):
        for mod_index in range(len(mod_data)):
            mod_id = int(configured_modules[mod_index]) if configured_modules is not None else mod_index
            ax.text(camera_layout.module(mod_id).x(), camera_layout.module(mod_id).y(),
                '%d'%mod_id, ha='center', va='center', fontsize=module_label_font_size,
                color=module_label_color)
    ax.axis('square')
    ax.axis(np.asarray([-1,1,-1,1])*1.05*max_xy)
    ax.set_xlabel('X coordinate [cm]')
    ax.set_ylabel('Y coordinate [cm]')
    cbar = plt.colorbar(pc, ax=ax)
    if cbar_label is not None:
        cbar.set_label(cbar_label)
    return pc

def layout_to_polygon_vxy(layout, plate_scale = 1.0, rotation = 0.0):
    crot = np.cos(rotation/180.0*np.pi)
    srot = np.sin(rotation/180.0*np.pi)
    all_vxy = []
    ibegin = 0;
    for iend in layout.outline_polygon_vertex_index():
        vx = layout.outline_polygon_vertex_x()[ibegin:iend]*plate_scale
        vy = layout.outline_polygon_vertex_y()[ibegin:iend]*plate_scale
        vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
        all_vxy.append(np.column_stack([vx, vy]))
        ibegin=iend
    return all_vxy

def layout_to_polygon(layout, plate_scale = 1.0, rotation = 0.0, **args):
    all_p = []
    for vxy in layout_to_polygon_vxy(layout, plate_scale, rotation):
        all_p.append(plt.Polygon(vxy, **args))
    if len(all_p) == 1:
        return all_p[0]
    return all_p

def add_outline(axis, layout, plate_scale = 1.0, rotation = 0.0, fill=False,
        outline_lw = 0.5, outline_ls = '-', outline_color = '#888888', **args):
    all_p = []
    for vxy in layout_to_polygon_vxy(layout, plate_scale, rotation):
        p = axis.add_patch(plt.Polygon(vxy,
            fill=fill, lw=outline_lw, ls=outline_ls, edgecolor=outline_color, **args))
        all_p.append(p)
    return all_p

def plot_camera_image(channel_data, camera_layout, channel_mask = None,
        configured_channels = None, zero_suppression = None,
        plate_scale = None, rotation = 0.0, R = None,
        cmap = matplotlib.cm.CMRmap_r, axis = None, draw_outline = False,
        pix_lw = 0, outline_lw = 0.5, outline_color = '#888888'):
    if(channel_mask is None and zero_suppression is not None):
        channel_mask = np.asarray(channel_data)>zero_suppression
    if(channel_mask is not None and len(channel_mask) != len(channel_data)):
        raise ValueError('channel_mask must either be None or have same length as channel_data')
    if(configured_channels is not None and len(configured_channels) != len(channel_data)):
        raise ValueError('configured_channels must either be None or have same length as pix_data')
    pix = []
    pix_data = []
    pix_gridid = []
    max_xy = 0
    if plate_scale is None:
        plate_scale = 1
        if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
            plate_scale = 180/np.pi/camera_layout.effective_focal_length()
    crot = np.cos(rotation/180.0*np.pi)
    srot = np.sin(rotation/180.0*np.pi)
    if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
        camera_layout = camera_layout.camera()
    for chan_index in range(len(channel_data)):
        chan_id = int(configured_channels[chan_index]) if configured_channels is not None else chan_index
        chan = camera_layout.channel(chan_id)
        if(chan.pixel_index() != -1):
            pix_gridid.append(chan.pixel_grid_index())
        if((channel_mask is None or channel_mask[chan_index]) and chan.pixel_index() != -1):
            vx = chan.outline_polygon_vertex_x_view()*plate_scale
            vy = chan.outline_polygon_vertex_y_view()*plate_scale
            vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
            vv = np.column_stack([vx, vy])
            max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
            pix.append(plt.Polygon(vv,closed=True))
            pix_data.append(channel_data[chan_index])
    pc = matplotlib.collections.PatchCollection(pix, cmap=cmap)
    pc.set_array(np.asarray(pix_data))
    pc.set_linewidths(pix_lw)
    axis = axis or plt.gca()
    axis.add_collection(pc)

    if draw_outline:
        add_outline(axis, camera_layout, plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    axis.axis('square')
    axis.axis(np.asarray([-1,1,-1,1])*(R or 1.05*max_xy))
    return pc

def plot_camera_module_image(module_data, camera_layout, module_mask = None,
        configured_modules = None, zero_suppression = None,
        plate_scale = None, rotation = 0.0, R = None,
        cmap = matplotlib.cm.CMRmap_r, axis = None, draw_outline = False,
        mod_lw = 0, outline_lw = 0.5, outline_color = '#888888'):
    if(module_mask is None and zero_suppression is not None):
        module_mask = np.asarray(module_data)>zero_suppression
    if(module_mask is not None and len(module_mask) != len(module_data)):
        raise ValueError('module_mask must either be None or have same length as module_data')
    if(configured_modules is not None and len(configured_modules) != len(module_data)):
        raise ValueError('configured_modules must either be None or have same length as pix_data')
    poly = []
    poly_data = []
    max_xy = 0
    if plate_scale is None:
        plate_scale = 1
        if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
            plate_scale = 180/np.pi/camera_layout.effective_focal_length()
    crot = np.cos(rotation/180.0*np.pi)
    srot = np.sin(rotation/180.0*np.pi)
    if type(camera_layout) is calin.ix.iact_data.instrument_layout.TelescopeLayout:
        camera_layout = camera_layout.camera()
    for mod_index in range(len(module_data)):
        mod_id = int(configured_modules[mod_index]) if configured_modules is not None else mod_index
        chan = camera_layout.module(mod_id)
        if((module_mask is None or module_mask[mod_index]) and chan.module_index() != -1):
            vx = chan.outline_polygon_vertex_x_view()*plate_scale
            vy = chan.outline_polygon_vertex_y_view()*plate_scale
            vx, vy = vx*crot - vy*srot, vy*crot + vx*srot
            vv = np.column_stack([vx, vy])
            max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
            poly.append(plt.Polygon(vv,closed=True))
            poly_data.append(module_data[mod_index])
    pc = matplotlib.collections.PatchCollection(poly, cmap=cmap)
    pc.set_array(np.asarray(poly_data))
    pc.set_linewidths(mod_lw)
    axis = axis or plt.gca()
    axis.add_collection(pc)

    if draw_outline:
        add_outline(axis, camera_layout, plate_scale=plate_scale, rotation=rotation,
            outline_lw=outline_lw, outline_color=outline_color)

    axis.axis('square')
    axis.axis(np.asarray([-1,1,-1,1])*(R or 1.05*max_xy))
    return pc

def plot_histogram(h, density = False, normalise = False,
        xscale = 1, xoffset = 0, yscale = 1, yoffset = 0,
        xscale_as_log10 = False, draw_poisson_errors = False,
        histtype='steps', ecolor=None, axis = None, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
        hdy = np.sqrt(h.all_weight())
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        if(h.sparse_bins_size() > 0):
            h = calin.math.histogram.densify(h)
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins()
        hdy = np.sqrt(h.bins())
    else:
        raise Exception('Unknown histogram type: '+str(type(h)))
    if density:
        hy /= abs(h.dxval()*xscale)
        hdy /= abs(h.dxval()*xscale)
    if normalise:
        hy /= h.sum_w()
        hdy /= h.sum_w()
    hx = np.append(hx, hx[-1]+h.dxval()) * xscale + xoffset
    hy = np.append(hy, hy[-1]) * yscale + yoffset
    hdy = hdy * yscale + yoffset
    if(xscale_as_log10):
        hx = 10**hx
    axis = axis or plt.gca()
    if(histtype=='step' or histtype=='steps'):
        so = axis.step(hx,hy, *args, where='post', **nargs)
        if(draw_poisson_errors):
            so1 = axis.vlines(0.5*(hx[:-1]+hx[1:]), hy[:-1]-hdy, hy[:-1]+hdy)
            so1.set_linestyles(so[0].get_linestyle())
            so1.set_color(so[0].get_color() if ecolor is None else ecolor)
            so1.set_linewidth(so[0].get_linewidth())
            so1.set_zorder(so[0].get_zorder())
            so = [ *so, so1 ]
    elif(histtype=='bar' or histtype=='bars'):
        yerr = None
        if(draw_poisson_errors):
            yerr = hdy
        so = axis.bar(hx[:-1], hy[:-1], *args, width=hx[1:]-hx[:-1], align='edge',
            yerr=yerr, ecolor='black' if ecolor is None else ecolor, **nargs)
    else:
        raise Exception('Unknown histogram plotting type: '+histtype)
    if(xscale_as_log10):
        so[0].axes.set_xscale('log')
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
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins()
    else:
        raise Exception('Unknown histogram type: '+str(type(h)))
    if plot_as_cdf or plot_as_cmf:
        hy /= h.sum_w()
    hx = np.append(hx, hx[-1]+h.dxval()) * xscale + xoffset
    hy = np.append(0, hy)
    if right_to_left:
        hy = np.sum(hy)-np.cumsum(hy)
    else:
        hy = np.cumsum(hy)
    hy =  hy * yscale + yoffset
    if(xscale_as_log10):
        hx = 10**hx
    axis = axis or plt.gca()
    so = axis.plot(hx, hy, *args, **nargs)
    if(xscale_as_log10):
        so[0].axes.set_xscale('log')
    return so
