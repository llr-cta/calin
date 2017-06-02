# calin/python/simulation/vs_cta.py -- Stephen Fegan -- 2017-06-01
#
# Functions for returning instances of VSO for CTA arrays
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

import numpy
import calin.math.hex_array
import calin.simulation.vs_optics

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

def mstn1_config(obscure_camera = True):
    mst = calin.ix.simulation.vs_optics.IsotropicDCArrayParameters()
    mst.mutable_array_origin().set_latitude(dms(28, 45, 47.36))
    mst.mutable_array_origin().set_longitude(dms(-17, 53, 23.93))
    mst.mutable_array_origin().set_elevation(2147 * 100.0)
    scope = mst.mutable_prescribed_array_layout().add_scope_positions();
    scope.set_z(mst.array_origin().elevation())
    mst.mutable_reflector_frame().set_optic_axis_rotation(-90);
    dc = mst.mutable_reflector()
    dc.set_curvature_radius(1920)
    dc.set_aperture(1230)
    dc.set_facet_num_hex_rings(5)
    dc.mutable_psf_align().set_object_plane(numpy.inf) # 10 * 1e5);
    dc.set_alignment_image_plane(1600)
    dc.set_facet_spacing(122)
    dc.set_facet_size(120)
    dc.set_facet_focal_length(1607)
    dc.set_facet_focal_length_dispersion(1)
    dc.set_facet_spot_size_probability(0.8)
    dc.set_facet_spot_size(0.5 * 2.8) # Spot size of 28mm at 2F
    dc.set_facet_spot_size_dispersion(0.5 * 0.02)
    dc.set_facet_labeling_parity(True)
    dc.set_weathering_factor(1.0)
    for id in [1,62,67,72,77,82,87]: dc.add_facet_missing_list(id-1)
    mst.mutable_focal_plane().set_camera_diameter(235)
    mst.mutable_focal_plane().mutable_translation().set_y(1/(1.0/dc.alignment_image_plane()-1/(10 * 1e5)))
    mst.mutable_pixel().set_spacing(5)
    mst.mutable_pixel().set_cone_inner_diameter(5)
    mst.mutable_pixel().set_cone_survival_prob(1)
    mst.mutable_pixel().set_hex_module_size(1)
    mst.mutable_pixel().set_module_num_hex_rings(9)
    u1,v1 = calin.math.hex_array.cluster_hexid_to_center_uv(1,1)
    x1,y1 = calin.math.hex_array.uv_to_xy(u1,v1)
    rot = numpy.arctan2(-y1,x1)/numpy.pi*180 - 30
    mst.mutable_pixel().set_grid_rotation(rot)

    if(obscure_camera):
        obs_camera_box = mst.add_obscurations()
        obs_camera_box.aligned_box().max_corner().set_x(150)
        obs_camera_box.aligned_box().max_corner().set_y(mst.focal_plane().translation().y()+150)
        obs_camera_box.aligned_box().max_corner().set_z(150)
        obs_camera_box.aligned_box().min_corner().set_x(-150)
        obs_camera_box.aligned_box().min_corner().set_y(mst.focal_plane().translation().y())
        obs_camera_box.aligned_box().min_corner().set_z(-150)
        obs_camera_box.aligned_box().set_incoming_only(True)

    return mst

def make_array(cfg, rng = calin.math.rng.RNG()):
    array = calin.simulation.vs_optics.VSOArray()
    array.generateFromArrayParameters(cfg, rng)
    return array
