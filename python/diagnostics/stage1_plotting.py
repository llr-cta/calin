# calin/python/diagnostics/stage1_plotting.py -- Stephen Fegan -- 2021-06-23
#
# Rendering of custom stage 1 diagnostics plots
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

def draw_missing_components_fraction(stage1, cmap = matplotlib.cm.CMRmap_r, axis = None,
        draw_outline = True, mod_lw = 0, outline_lw = 0.5, outline_color = '#888888',
        mod_label_fontsize=4, aux_label_fontsize=5.5):
    SQRT3_2 = numpy.sqrt(3)/2
    RAD_30 = 30.0/180.0*numpy.pi

    rc = stage1.const_run_config()
    cl = rc.const_camera_layout()
    ri = stage1.const_run_info()

    xmin = numpy.min(cl.outline_polygon_vertex_x())
    xmax = numpy.max(cl.outline_polygon_vertex_x())
    ymin = numpy.min(cl.outline_polygon_vertex_y())
    ymax = numpy.max(cl.outline_polygon_vertex_y())

    # Delete these if you want the additional hardware to fit within the
    # camera rectangle, rather than the square
    xmin = min(xmin,ymin)
    ymin = xmin
    xmax = max(xmax,ymax)
    ymax = xmax

    num_event_missing = int(ri.event_numbers_found().end_index()[-1])-ri.num_events_found()-1+int(sum(ri.duplicate_event_numbers().count()))
    frac_event_missing = float(num_event_missing)/float(ri.event_numbers_found().end_index()[-1])

    additional_polygons = [ ]
    additional_polygon_data = [ ri.num_events_missing_tib()/ri.num_events_found(),
                                ri.num_events_missing_cdts()/ri.num_events_found(),
                                ri.num_events_missing_tib_and_cdts()/ri.num_events_found(),
                                frac_event_missing]

    polygon_name = [ 'tib', 'cdts', 'both', 'no\nevent' ]

    rad = 10
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(xmin+rad,ymax-rad*SQRT3_2),numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(xmin+rad,ymax-rad*SQRT3_2-2*rad),numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.RegularPolygon(xy=(xmin+rad+2*rad*SQRT3_2,ymax-rad*SQRT3_2-rad),numVertices=6, radius=rad, orientation=RAD_30))
    additional_polygons.append(
        matplotlib.patches.Circle(xy=(xmax-1.5*rad,ymin+1.5*rad), radius=1.5*rad))

    axis = axis or matplotlib.pyplot.gca()

    calin.plotting.add_outline(axis, cl, hatch='//', zorder=-1)

    modmissing = numpy.asarray([max(float(ri.num_events_found()-ri.module(i).num_events_present()),1e-10)/ri.num_events_found() \
        for i in range(ri.module_size())],dtype=float)

    pc = calin.plotting.plot_camera_module_image(modmissing, cl,
        configured_modules=rc.configured_module_id(),
        additional_polygons=additional_polygons,
        additional_polygon_data=additional_polygon_data,
        mod_lw=mod_lw, outline_lw=outline_lw, outline_color=outline_color,
        draw_outline=draw_outline, axis=axis)

    vmin = 0.5/ri.num_events_found()
    vmin_color_change = vmin**0.4
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

    pc.set_norm(matplotlib.colors.LogNorm(vmin=vmin,vmax=1.0))
    axis.axis(numpy.asarray([-1,1,-1,1])*1.05*max([ymax,xmax,abs(xmin),abs(ymin)]))
    axis.get_xaxis().set_visible(False)
    axis.get_yaxis().set_visible(False)

    cb = matplotlib.pyplot.colorbar(pc, ax=axis, label='Missing event fraction')

    return pc
