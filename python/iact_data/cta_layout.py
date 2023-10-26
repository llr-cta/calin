# calin/python/simulation/cta_layout.py -- Stephen Fegan -- 2022-10-13
#
# Functions for CTA instriument layouts
#
# Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
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

import calin.ix.iact_data.instrument_layout
import calin.math.regular_grid

def configure_analog_trigger_patches_3(camera_layout):
    camera_layout = camera_layout.Clone()

    module_index = dict()
    for im in range(camera_layout.module_size()):
        u = camera_layout.module(im).module_grid_u()
        v = camera_layout.module(im).module_grid_v()
        module_index[(u,v)] = im

    patch_channel = []
    channel_patch = [ [] for _ in range(camera_layout.channel_size())]
    for im in range(camera_layout.module_size()):
        u = camera_layout.module(im).module_grid_u()
        v = camera_layout.module(im).module_grid_v()
        if((u+1,v) in module_index and (u,v+1) in module_index):
            channels = []
            for jm in [ im, module_index[(u+1,v)], module_index[(u,v+1)] ]:
                for c in camera_layout.module(jm).channels_in_module():
                    channels.append(c)
                    channel_patch[c].append(len(patch_channel))
            patch_channel.append(sorted(channels))
        if((u+1,v) in module_index and (u+1,v-1) in module_index):
            channels = []
            for jm in [ im, module_index[(u+1,v)], module_index[(u+1,v-1)] ]:
                for c in camera_layout.module(jm).channels_in_module():
                    channels.append(c)
                    channel_patch[c].append(len(patch_channel))
            patch_channel.append(sorted(channels))

    grid = calin.iact_data.instrument_layout.make_grid_from_instrument_layout(camera_layout)
    camera_layout.clear_trigger_patches()
    for ipc, pc in enumerate(patch_channel):
        patch = camera_layout.add_trigger_patches()
        patch.set_patch_index(ipc)
        patch.set_channels_in_patch(pc)
        bcr = grid.compute_region_boundary(
            [camera_layout.channel(int(ic)).pixel_grid_index() for ic in patch.channels_in_patch()])
        bcx, bcy = grid.extract_boundary_curve(bcr)
        patch.set_outline_polygon_vertex_x(bcx)
        patch.set_outline_polygon_vertex_y(bcy)
        patch.set_outline_polygon_vertex_index([len(bcx)])
    for ic, cp in enumerate(channel_patch):
        channel = camera_layout.mutable_channel(ic)
        channel.set_trigger_patch_indexes(cp)

    return camera_layout
