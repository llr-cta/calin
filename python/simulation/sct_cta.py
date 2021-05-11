# calin/python/simulation/vs_cta.py -- Stephen Fegan -- 2021-04-26
#
# Functions for returning instances of SCT for CTA arrays
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

import numpy
import calin.math.rng
import calin.math.hex_array
import calin.simulation.sct_optics

def dms(d,m,s):
    # Note that "negative" d=0 (e.g. -00:30:00) must be specified as 00:-30:00 or 00:00:-30
    sign = 1
    if(d<0):
        sign = -1
        d = abs(d)
    elif(d==0 and m<0):
        sign = -1
        m = abs(m)
    elif(d==0 and m==0 and s<0):
        sign = -1
        s = abs(s)
    return sign * (d + m/60.0 + s/3600.0)

def sct1_config(obscure = True, scope_x=0, scope_y=0, include_window = False):
    from calin.simulation.sct_optics import COS_PI_8
    from calin.simulation.sct_optics import COS_PI_16
    from calin.simulation.sct_optics import COS_PI_32
    from calin.simulation.sct_optics import SIN_PI_8
    from calin.simulation.sct_optics import SIN_PI_16
    from calin.simulation.sct_optics import SIN_PI_32

    sct = calin.ix.simulation.sct_optics.SCTRandomArrayParameters()
    sct.mutable_array_origin().set_latitude(dms(28, 45, 47.36))
    sct.mutable_array_origin().set_longitude(dms(-17, 53, 23.93))
    sct.mutable_array_origin().set_elevation(2147 * 100.0)
    try:
        for i in range(max(len(scope_x), len(scope_y))):
            scope = sct.mutable_array_layout().add_scope_positions();
            scope.set_x(scope_x[i])
            scope.set_y(scope_y[i])
            scope.set_z(mst.array_origin().elevation())
    except:
        scope = sct.mutable_array_layout().add_scope_positions();
        scope.set_x(scope_x)
        scope.set_y(scope_y)
        scope.set_z(sct.array_origin().elevation())

    F = 558.630
    q = 2/3
    alpha = 2/3

    sct.set_focal_length(F)

    pp = numpy.asarray([
        1.4237e-6,  0.11107,  -6.4869e-3,  -6.0375e-3,  1.4912e-2,
        -5.6809e-2, 0.13774,  -0.24615,    0.30847,     -0.26204,
        0.13550,    -3.3061e-2])
    pp *= F * F**(-2*numpy.arange(0,len(pp)))

    sct.set_primary_sag_polynomial(pp)

    primary_facet_scheme = sct.mutable_primary_facet_scheme()
    # primary_facet_scheme.set_inner_ring_inner_radius(219.350*COS_PI_16);
    # primary_facet_scheme.set_inner_ring_outer_radius(340.000*COS_PI_32 - 0.7)
    # primary_facet_scheme.set_outer_ring_inner_radius(340.000*COS_PI_32 + 0.7)
    # primary_facet_scheme.set_outer_ring_outer_radius(483.1875*COS_PI_32)
    primary_facet_scheme.set_inner_ring_inner_radius(215.136);
    primary_facet_scheme.set_inner_ring_outer_radius(339.305*COS_PI_32)
    primary_facet_scheme.set_outer_ring_inner_radius(339.122)
    primary_facet_scheme.set_outer_ring_outer_radius(480.929)
    primary_facet_scheme.set_long_edge_half_gap(0.7)

    primary_facet_scheme_loose = sct.mutable_primary_facet_scheme_loose()
    primary_facet_scheme_loose.set_inner_ring_inner_radius(214.436);
    primary_facet_scheme_loose.set_inner_ring_outer_radius(338.397)
    primary_facet_scheme_loose.set_outer_ring_inner_radius(338.397)
    primary_facet_scheme_loose.set_outer_ring_outer_radius(481.629)
    primary_facet_scheme_loose.set_long_edge_half_gap(0.0)

    ps = numpy.asarray([
        4.5094360e-7,  -0.4167062710,  -0.1442213045,  0.67307608,   -3.539924707,
        15.8620298,    -59.11580194,   170.5778772,    -350.8596568, 452.9519853,
        -274.1880908])
    ps *= F * F**(-2*numpy.arange(0,len(ps)))

    rsec2o = 270.816

    sct.set_secondary_distance(F/q)
    sct.set_secondary_sag_polynomial(ps)

    secondary_facet_scheme = sct.mutable_secondary_facet_scheme()
    secondary_facet_scheme.set_inner_ring_inner_radius(38.867);
    secondary_facet_scheme.set_inner_ring_outer_radius(158.969*COS_PI_16)
    secondary_facet_scheme.set_outer_ring_inner_radius(157.383)
    secondary_facet_scheme.set_outer_ring_outer_radius(265.748)
    secondary_facet_scheme.set_long_edge_half_gap(0.7)

    secondary_facet_scheme_loose = sct.mutable_secondary_facet_scheme_loose()
    secondary_facet_scheme_loose.set_inner_ring_inner_radius(38.167);
    secondary_facet_scheme_loose.set_inner_ring_outer_radius(156.649)
    secondary_facet_scheme_loose.set_outer_ring_inner_radius(156.649)
    secondary_facet_scheme_loose.set_outer_ring_outer_radius(266.448)
    secondary_facet_scheme_loose.set_long_edge_half_gap(0.0)

    pc = numpy.asarray([0, -0.8327, 4.9950]) # Values from SCT-OPTMO/121108
    pc *= F * F**(-2*numpy.arange(0,len(pc)))
    # pc = numpy.asarray([-0.2499999999998863, -1.33436e-03, 2.86525e-8]) # Values from confluence

    sct.set_camera_distance(F*(1/q - (1-alpha)))
    sct.set_camera_sag_polynomial(pc)
    sct.set_camera_radius(81*0.5)

    s_disk_r = rsec2o
    s_disk_z = sct.secondary_distance() + numpy.polyval(numpy.flipud(ps), s_disk_r**2)

    s_baffle_r = 284.48
    s_baffle_zin = 780.5
    s_baffle_zout = s_baffle_zin - 100.33

    p_baffle_r = 497.84
    p_baffle_zin = 43.38
    p_baffle_zout = p_baffle_zin + 152.4

    if(obscure):
        # Incoming rays : disk to approximate secondary
        obs = sct.add_primary_obscuration()
        obs.mutable_circular_aperture().mutable_center_pos().set_y(s_disk_z)
        obs.mutable_circular_aperture().set_diameter(s_disk_r * 2)
        obs.mutable_circular_aperture().set_invert(True)

        # Incoming rays : tube to approximate secondary baffle
        obs = sct.add_primary_obscuration()
        obs.mutable_tube().mutable_end1_pos().set_y(s_baffle_zin)
        obs.mutable_tube().mutable_end2_pos().set_y(s_baffle_zout)
        obs.mutable_tube().set_diameter(s_baffle_r * 2)
        # sct.add_secondary_obscuration().CopyFrom(obs)

        # Incoming rays : circular aperture to approximate primary baffle
        obs = sct.add_primary_obscuration()
        obs.mutable_circular_aperture().mutable_center_pos().set_y(p_baffle_zout)
        obs.mutable_circular_aperture().set_diameter(p_baffle_r * 2)

        # Post-primary rays : circular aperture to approximate secondary baffle
        obs = sct.add_secondary_obscuration()
        obs.mutable_circular_aperture().mutable_center_pos().set_y(s_baffle_zout)
        obs.mutable_circular_aperture().set_diameter(s_baffle_r * 2)
    else:
        # Incoming rays : circular aperture to approximate primary baffle
        obs = sct.add_primary_obscuration()
        obs.mutable_circular_aperture().mutable_center_pos().set_y(p_baffle_zin)
        obs.mutable_circular_aperture().set_diameter(p_baffle_r * 2)

        # Post-primary rays : circular aperture to approximate secondary baffle
        obs = sct.add_secondary_obscuration()
        obs.mutable_circular_aperture().mutable_center_pos().set_y(s_baffle_zin)
        obs.mutable_circular_aperture().set_diameter(s_baffle_r * 2)

    if(include_window):
        win = sct.mutable_window()
        win.set_front_y_coord(sct.camera_distance() + 5)
        win.set_outer_radius(0.0) # plane window
        win.set_thickness(0.6)
        win.set_refractive_index(1.49)

    return sct
