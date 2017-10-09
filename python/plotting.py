# calin/python/plotting.py -- Stephen Fegan -- 2017-03-27
#
# Plotting of camera faces and histograms
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay
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

def plot_camera(pix_data, camera_layout, configured_channels = None, ax_in = None,
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

def layout_to_polygon_vxy(layout, plate_scale = 1.0):
    all_vxy = []
    ibegin = 0;
    for iend in layout.outline_polygon_vertex_index():
        vx = layout.outline_polygon_vertex_x()[ibegin:iend]*plate_scale
        vy = layout.outline_polygon_vertex_y()[ibegin:iend]*plate_scale
        all_vxy.append(np.column_stack([vx, vy]))
    return all_vxy

def layout_to_polygon(layout, plate_scale = 1.0, **args):
    all_p = []
    for vxy in layout_to_polygon_vxy(layout, plate_scale):
        all_p.append(plt.Polygon(vxy, **args))
    if len(all_p) == 1:
        return all_p[0]
    return all_p

def add_outline(axis, layout, plate_scale = 1.0,
        outline_lw = 0.5, outline_color = '#888888'):
    for vxy in layout_to_polygon_vxy(layout, plate_scale):
        axis.add_patch(plt.Polygon(vxy,
            fill=False, lw=outline_lw, edgecolor=outline_color))

def plot_camera_image(channel_data, camera_layout, channel_mask = None,
        configured_channels = None, zero_suppression = None,
        plate_scale = None, R = None,
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
        if type(camera_layout) is calin.ix.iact_data.instrument_layout.calin.ix.iact_data.instrument_layout.TelescopeLayout:
            plate_scale = 180/np.pi/camera_layout.effective_focal_length()
    if type(camera_layout) is calin.ix.iact_data.instrument_layout.calin.ix.iact_data.instrument_layout.TelescopeLayout:
        camera_layout = camera_layout.camera()
    for chan_index in range(len(channel_data)):
        chan_id = int(configured_channels[chan_index]) if configured_channels is not None else chan_index
        chan = camera_layout.channel(chan_id)
        if(chan.pixel_index() != -1):
            pix_gridid.append(chan.pixel_grid_index())
        if((channel_mask is None or channel_mask[chan_index]) and chan.pixel_index() != -1):
            vx = chan.outline_polygon_vertex_x_view()*plate_scale
            vy = chan.outline_polygon_vertex_y_view()*plate_scale
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
        add_outline(axis, camera_layout, plate_scale=plate_scale,
            outline_lw=outline_lw, outline_color=outline_color)

    axis.axis('square')
    axis.axis(np.asarray([-1,1,-1,1])*(R or 1.05*max_xy))
    return pc

def plot_histogram(h, plot_as_pdf = False, plot_as_pmf = False, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins()
    else:
        raise Exception('Unknown histogram type: '+type(h))
    if plot_as_pdf:
        hy /= h.sum_w()/h.dxval()
    elif plot_as_pmf:
        hy /= h.sum_w()
    hx = np.append(hx, hx[-1]+h.dxval())
    hy = np.append(hy, hy[-1])
    so = plt.step(hx,hy, where='post', *args, **nargs)
    return so

def plot_histogram_cumulative(h,  plot_as_pdf = False, plot_as_pmf = False, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins()
    else:
        raise Exception('Unknown histogram type: '+type(h))
    if plot_as_pdf or plot_as_pmf:
        hy /= h.sum_w()
    hx = np.append(hx, hx[-1]+h.dxval())
    hy = np.append(0, hy)
    so = plt.plot(hx,np.cumsum(hy), *args, **nargs)
    return so
