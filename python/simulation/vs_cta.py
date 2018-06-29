# calin/python/simulation/vs_cta.py -- Stephen Fegan -- 2017-06-01
#
# Functions for returning instances of VSO for CTA arrays
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

def mstn1_config(obscure_camera = True, scope_x=0, scope_y=0):
    mst = calin.ix.simulation.vs_optics.IsotropicDCArrayParameters()
    mst.mutable_array_origin().set_latitude(dms(28, 45, 47.36))
    mst.mutable_array_origin().set_longitude(dms(-17, 53, 23.93))
    mst.mutable_array_origin().set_elevation(2147 * 100.0)
    try:
        for i in range(max(len(scope_x), len(scope_y))):
            scope = mst.mutable_prescribed_array_layout().add_scope_positions();
            scope.set_x(scope_x[i])
            scope.set_y(scope_y[i])
            scope.set_z(mst.array_origin().elevation())
    except:
        scope = mst.mutable_prescribed_array_layout().add_scope_positions();
        scope.set_x(scope_x)
        scope.set_y(scope_y)
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
    #for id in [1,62,67,72,77,82,87]: dc.add_facet_missing_list(id-1) # 84 mirror
    for id in [1,67,72,82,87]: dc.add_facet_missing_list(id-1) # 84 mirror
    mst.mutable_focal_plane().set_camera_diameter(235)
    mst.mutable_focal_plane().mutable_translation().set_y(1/(1.0/dc.alignment_image_plane()-1/(10 * 1e5)))
    mst.mutable_pixel().set_spacing(5)
    mst.mutable_pixel().set_cone_inner_diameter(5)
    mst.mutable_pixel().set_cone_survival_prob(1)
    mst.mutable_pixel().set_hex_module_size(1)
    mst.mutable_pixel().set_hex_module_layout_use_b_configuration(True)
    mst.mutable_pixel().set_module_num_hex_rings(9)
    u1,v1 = calin.math.hex_array.cluster_hexid_to_center_uv(1,1,False)
    x1,y1 = calin.math.hex_array.uv_to_xy(u1,v1)
    rot = numpy.arctan2(-y1,x1)/numpy.pi*180 + 30
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

def lst1_config(obscure_camera = True, scope_x=0, scope_y=0):
    lst = calin.ix.simulation.vs_optics.IsotropicDCArrayParameters()
    lst.mutable_array_origin().set_latitude(dms(28, 45, 47.36))
    lst.mutable_array_origin().set_longitude(dms(-17, 53, 23.93))
    lst.mutable_array_origin().set_elevation(2147 * 100.0)
    try:
        for i in range(max(len(scope_x), len(scope_y))):
            scope = lst.mutable_prescribed_array_layout().add_scope_positions();
            scope.set_x(scope_x[i])
            scope.set_y(scope_y[i])
            scope.set_z(lst.array_origin().elevation())
    except:
        scope = lst.mutable_prescribed_array_layout().add_scope_positions();
        scope.set_x(scope_x)
        scope.set_y(scope_y)
        scope.set_z(lst.array_origin().elevation())
    lst.mutable_reflector_frame().set_optic_axis_rotation(0);
    dc = lst.mutable_reflector()
    dc.set_curvature_radius(5644) # Least square fit of sphere to mirror pos. Could be 5600 either
    dc.mutable_psf_align().set_object_plane(numpy.inf) # 10 * 1e5);
    dc.set_alignment_image_plane(2800)
    dc.set_aperture(2*1153.3)
    #dc.set_facet_num_hex_rings(8)
    dc.set_facet_spacing(154)
    dc.set_facet_size(151)
    dc.set_facet_grid_shift_z(44.46)
    dc.set_facet_focal_length(0)
    dc.set_facet_focal_length_dispersion(10)
    dc.set_facet_spot_size_probability(0.8)
    dc.set_facet_spot_size(0.5 * 1.4) # Spot size of 14mm at 2F
    dc.set_facet_spot_size_dispersion(0.5 * 0.02)
    dc.set_facet_labeling_parity(True)
    dc.set_weathering_factor(1.0)
    for id in [0,127,148,171,191]: dc.add_facet_missing_list(id)
    lst.mutable_focal_plane().set_camera_diameter(235)
    lst.mutable_focal_plane().mutable_translation().set_y(1/(1.0/dc.alignment_image_plane()-1/(10 * 1e5)))
    lst.mutable_pixel().set_spacing(5)
    lst.mutable_pixel().set_cone_inner_diameter(5)
    lst.mutable_pixel().set_cone_survival_prob(1)
    lst.mutable_pixel().set_hex_module_size(1)
    lst.mutable_pixel().set_hex_module_layout_use_b_configuration(True)
    lst.mutable_pixel().set_module_num_hex_rings(9)
    u1,v1 = calin.math.hex_array.cluster_hexid_to_center_uv(1,1,False)
    x1,y1 = calin.math.hex_array.uv_to_xy(u1,v1)
    rot = numpy.arctan2(-y1,x1)/numpy.pi*180 + 30
    lst.mutable_pixel().set_grid_rotation(rot)

    if(obscure_camera):
        obs_camera_box = lst.add_obscurations()
        obs_camera_box.aligned_box().max_corner().set_x(150)
        obs_camera_box.aligned_box().max_corner().set_y(lst.focal_plane().translation().y()+150)
        obs_camera_box.aligned_box().max_corner().set_z(150)
        obs_camera_box.aligned_box().min_corner().set_x(-150)
        obs_camera_box.aligned_box().min_corner().set_y(lst.focal_plane().translation().y())
        obs_camera_box.aligned_box().min_corner().set_z(-150)
        obs_camera_box.aligned_box().set_incoming_only(True)

    return lst

def make_array(cfg, rng = calin.math.rng.RNG('calin.simulation.vs_cta.make_array')):
    array = calin.simulation.vs_optics.VSOArray()
    if(type(cfg) is list):
        for icfg in cfg:
            array.generateFromArrayParameters(icfg, rng)
    else:
        array.generateFromArrayParameters(cfg, rng)
    return array
