# calin/python/simulation/vs_cta.py -- Stephen Fegan -- 2017-06-01
#
# Functions for returning instances of VSO for CTA arrays
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

import numpy
import calin.math.hex_array
import calin.simulation.vs_optics
import calin.simulation.detector_efficiency
import calin.provenance.system_info

def mstn_detection_efficiency(qe = 'qe_R12992-100-05b.dat',
        mirror = 'ref_MST-North-MLT_2022_06_28.dat',
        window = 'transmission_lst_window_No7-10_ave.dat',
        degradation_factor = 0.9, quiet = False):
    data_dir = calin.provenance.system_info.build_info().data_install_dir() + "/simulation/"
    det_eff = calin.simulation.detector_efficiency.DetectionEfficiency()

    # PMT quantum efficiency curve
    if(qe):
        det_eff.scaleEffFromFile(data_dir + qe)
        if(not quiet):
            print('Loading PMT QE curve :',qe,'[bandwidth = %.3f]'%det_eff.integrate())

    # Mirror reflectivity curve
    if(mirror):
        det_eff.scaleEffFromFile(data_dir + mirror)
        if(not quiet):
            print('Scaling by mirror reflectivity :',mirror,'[bandwidth = %.3f]'%det_eff.integrate())

    # Window transmission curve
    if(window):
        det_eff.scaleEffFromFile(data_dir + window)
        if(not quiet):
            print('Scaling by window transmission :',window,'[bandwidth = %.3f]'%det_eff.integrate())

    # Optional constant degradation factor
    if(degradation_factor < 1.0):
        det_eff.scaleEffByConst(degradation_factor)
        if(not quiet):
            print('Scaling by degradation factor :',degradation_factor,'[bandwidth = %.3f]'%det_eff.integrate())

    return det_eff

def mstn_cone_efficiency(cone = 'NectarCAM_lightguide_efficiency_POP_131019.dat',
        quiet = False):
    data_dir = calin.provenance.system_info.build_info().data_install_dir() + "/simulation/"

    if(cone):
        cone_eff = calin.simulation.detector_efficiency.AngularEfficiency(data_dir + cone)
        if(not quiet):
            print('Loading light-cone curve :',cone)

    return cone_eff

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

def mstn1_config(obscure_camera = True, scope_x=0, scope_y=0, include_window = False):
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
        obs_camera_box = mst.add_pre_reflection_obscuration()
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_x(2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_z(2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_y(mst.focal_plane().translation().y()+(1532.25-513.0)/10)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_x(-2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_z(-2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_y(mst.focal_plane().translation().y()-513.0/10)

        obs_outer_aperture = mst.add_camera_obscuration()
        obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_x(0)
        obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_z(-50.0/20)
        obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_y(-513.0/10)
        obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_x(2714.0/10)
        obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_z((2585.0-50.0)/10)

        obs_inner_aperture = mst.add_camera_obscuration()
        obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_x(0)
        obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_z(0)
        obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_y(-222.5/10)
        obs_inner_aperture.mutable_circular_aperture().set_diameter(2304.0/10)

    if include_window:
        win = mst.mutable_spherical_window()
        win.set_front_y_coord(mst.focal_plane().translation().y() - 427.5/10)
        win.set_outer_radius(3443.0/10)
        win.set_thickness(5.25/10.0)
        win.set_refractive_index(1.5)

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
        obs_camera_box = lst.add_post_reflection_obscuration()
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_x(150)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_y(lst.focal_plane().translation().y()+150)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_z(150)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_x(-150)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_y(lst.focal_plane().translation().y())
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_z(-150)

    return lst

def make_array(cfg, rng = calin.math.rng.RNG('calin.simulation.vs_cta.make_array')):
    array = calin.simulation.vs_optics.VSOArray()
    if(type(cfg) is list):
        for icfg in cfg:
            array.generateFromArrayParameters(icfg, rng)
    else:
        array.generateFromArrayParameters(cfg, rng)
    return array
