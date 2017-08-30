# calin/python/plotting.py -- Stephen Fegan -- 2017-07-15
#
# Plotting of simulation elements
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
import calin.math.hex_array
import calin.math.regular_grid
import numpy as np

def plot_mirrors(scope, label_hex_id=False, scale=1, scale_units='cm', ax_in=None, fs=None):
    max_xy = 0
    pix = []
    for imirror in range(scope.numMirrors()):
        mirror = scope.mirror(imirror)
        vx = []
        vy = []
        for i in range(6):
            v = mirror.cornerInReflectorCoords(i,scope.facetSize())
            vx.append(v[0])
            vy.append(v[2])
        vx = np.asarray(vx) * scale
        vy = np.asarray(vy) * scale
        vv = np.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        max_xy = max(max_xy, max(abs(vx)), max(abs(vy)))
        pix.append(plt.Polygon(vv,closed=True))
    ax = ax_in if ax_in is not None else plt.gca()
    pc = matplotlib.collections.PatchCollection(pix)
    #pc.set_linewidths(0)
    ax.add_collection(pc)
    ax.axis('square')
    ax.axis(np.asarray([-1,1,-1,1])*1.05*max_xy)
    ax.set_xlabel('X coordinate [%s]'%scale_units)
    ax.set_ylabel('Y coordinate [%s]'%scale_units)
    if(label_hex_id):
        for imirror in range(scope.numMirrors()):
            mirror = scope.mirror(imirror)
            args = dict()
            if fs is not None:
                args['fontsize'] = fs
            ax.text(mirror.pos()[0]*scale, mirror.pos()[2]*scale, '%d'%(mirror.hexID()+1),\
                ha='center',va='center',**args)
    return pc

def plot_image(scope, pix_data, cmap=None, clim=None, draw_outline=True, \
        plate_scale=None, ax_in=None, R=None, zero_suppress=None):
    pix = []
    idx = []
    plate_scale = scope.pixelSpacing()*(plate_scale or 1/scope.focalPlanePosition()[1]/np.pi*180.0)
    for pix_id in range(len(pix_data)):
        if(zero_suppress is not None and pix_data[pix_id]<=zero_suppress):
            continue
        pix_hexid = scope.pixel(pix_id).hexID()
        vx,vy = calin.math.hex_array.hexid_to_vertexes_xy_trans(pix_hexid,
            scope.cosPixelRotation(), scope.sinPixelRotation(), plate_scale)
        vv = np.zeros((len(vx),2))
        vv[:,0] = vx
        vv[:,1] = vy
        pix.append(plt.Polygon(vv,closed=True))
        idx.append(pix_id)

    ax = ax_in if ax_in is not None else plt.gca()
    pc = matplotlib.collections.PatchCollection(pix, cmap=cmap or matplotlib.cm.CMRmap_r)
    pc.set_array(np.asarray(pix_data)[idx])
    pc.set_linewidths(0)
    if(clim is not None):
        pc.set_clim(clim[0],clim[1])
    ax.add_collection(pc)

    if draw_outline:
        cam_hexids = list(map(lambda p: p.hexID(), scope.all_pixels()))
        grid = calin.math.regular_grid.HexGrid(plate_scale,scope.pixelRotation())
        for icurve in range(grid.num_bounday_curves(v)):
            vx,vy = grid.extract_bounday_curve(v,icurve,False)
            ax.add_patch(plt.Polygon(np.column_stack([vx, vy]),
                                fill=False,lw=0.5,edgecolor='#888888'))

    ax.axis('square')
    R = R or max(abs(np.asarray(ax.axis())))
    if R:
        ax.axis(np.asarray([-1,1,-1,1])*R)

    return pc
