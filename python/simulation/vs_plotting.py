# calin/python/plotting.py -- Stephen Fegan -- 2017-07-15
#
# Plotting of simulation elements
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole Polytechnique, CNRS/IN2P3
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
import calin.plotting
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
    if type(scope) is calin.simulation.vs_optics.VSOTelescope:
        scope = scope.convert_to_telescope_layout()
    elif type(scope) is calin.ix.simulation.vs_optics.IsotropicDCArrayParameters:
        scope = calin.simulation.vs_optics.dc_parameters_to_telescope_layout(scope)
    pc = calin.plotting.plot_camera_image(pix_data, scope,
                cmap=cmap, draw_outline=draw_outline,
                plate_scale=plate_scale, axis=ax_in, R=R,
                zero_suppression=zero_suppress)
    if clim is not None:
        pc.set_clim(clim)
    return pc
