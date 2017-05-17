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

def plot_camera(pix_data, camera_layout, configured_channels = None, ax_in = None,
        cbar_label = None):
    if(configured_channels is not None and len(configured_channels) != len(pix_data)):
        raise ValueError('configured_channels must either be None or have same length as pix_data')
    pix = []
    max_xy = 0
    for pix_index in range(len(pix_data)):
        pix_id = int(configured_channels[pix_index]) if configured_channels is not None else  pix_index
        vx = camera_layout.channel(pix_id).pixel_polygon_vertex_x_view()
        vy = camera_layout.channel(pix_id).pixel_polygon_vertex_y_view()
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

def plot_histogram(h, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins_view()
    else:
        raise Exception('Unknown histogram type: '+type(h))
    so = plt.step(hx,hy, where='post', *args, **nargs)
    return so

def plot_histogram_cumulative(h, *args, **nargs):
    if type(h) is calin.math.histogram.SimpleHist:
        hx = h.all_xval_left()
        hy = h.all_weight()
    elif type(h) is calin.ix.math.histogram.Histogram1DData:
        hx = h.xval0()+h.dxval()*np.arange(0,h.bins_size())
        hy = h.bins_view()
    else:
        raise Exception('Unknown histogram type: '+type(h))
    so = plt.step(hx,np.cumsum(hy), where='post', *args, **nargs)
    return so
