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
import scipy.integrate

import calin.math.hex_array
import calin.simulation.vs_optics
import calin.simulation.atmosphere
import calin.simulation.detector_efficiency
import calin.provenance.system_info
import calin.provenance.chronicle
import calin.util.utm

def ds_filename(filename):
    if(filename.find('/') >= 0):
        return filename
    else:
        data_dir = calin.provenance.system_info.build_info().data_install_dir() + "/simulation/"
        return data_dir + filename

def load_assets(filename, utm_zone, utm_hemi,
                desired_asset_types = set(['LSTN', 'MSTN', 'LSTS', 'MSTS', 'SSTS']),
                reproject_positions = True, apply_corrections = True,
                demand_ref_lon=None, demand_ref_lat=None,
                demand_ref_alt=None, set_demand_ref_from_header=True):
    ellipse = calin.util.utm.wgs84_ellipse()
    assets = []
    with open(filename, 'r') as file:
        file_record = calin.provenance.chronicle.register_file_open(filename,
            calin.ix.provenance.chronicle.AT_READ, 'calin.simulation.vs_cta.load_assets')
        comment = ""
        for line in file.readlines():
            if(line[0]=='#'):
                if(set_demand_ref_from_header):
                    if(line.find("center_lon:")>0):
                        i0 = line.find(": \"")+3
                        i1 = line.find("deg\"")
                        demand_ref_lon = float(line[i0:i1])/180.0*numpy.pi
                    if(line.find("center_lat:")>0):
                        i0 = line.find(": \"")+3
                        i1 = line.find("deg\"")
                        demand_ref_lat = float(line[i0:i1])/180.0*numpy.pi
                    if(line.find("center_alt:")>0):
                        i0 = line.find(": \"")+3
                        i1 = line.find("m\"")
                        demand_ref_alt = float(line[i0:i1])
                comment += line
            else:
                asset_data = line.split()
                if(asset_data[0] in desired_asset_types):
                    _,lat,lon = calin.util.utm.grid_to_geographic(ellipse.a,ellipse.e2,utm_zone,utm_hemi,
                        float(asset_data[3]),float(asset_data[2]))
                    assets.append([asset_data[0], asset_data[1],
                        lon, lat, float(asset_data[4]), float(asset_data[2]), float(asset_data[3])])
        file_record.set_comment(comment)
        calin.provenance.chronicle.register_file_close(file_record)

    if(demand_ref_alt is None):
        ref_alt = numpy.mean([a[4] for a in assets])
    else:
        ref_alt = demand_ref_alt

    if(demand_ref_lon is None or demand_ref_lat is None):
        ref_E = 0.5*(numpy.min([a[5] for a in assets]) + numpy.max([a[5] for a in assets]))
        ref_N = 0.5*(numpy.min([a[6] for a in assets]) + numpy.max([a[6] for a in assets]))
        _,lat_ref,lon_ref = calin.util.utm.grid_to_geographic(ellipse.a,ellipse.e2,utm_zone,utm_hemi,ref_N,ref_E)
    else:
        lat_ref = demand_ref_lat
        lon_ref = demand_ref_lon
        _, _, _, ref_N, ref_E = calin.util.utm.geographic_to_grid(ellipse.a,ellipse.e2,lat_ref,lon_ref,utm_zone,utm_hemi)

    if(reproject_positions):
        for i in range(3 if demand_ref_lon is None or demand_ref_lat is None else 1):
            assets_refit_NE = []
            for a in assets:
                N, E = calin.util.utm.geographic_to_tm(ellipse.a+ref_alt,ellipse.e2,1.0,lon_ref,0,0,a[3],a[2])
                assets_refit_NE.append([E,N])

            if(demand_ref_lon is None or demand_ref_lat is None):
                ref_E = 0.5*(numpy.min([a[0] for a in assets_refit_NE]) + numpy.max([a[0] for a in assets_refit_NE]))
                ref_N = 0.5*(numpy.min([a[1] for a in assets_refit_NE]) + numpy.max([a[1] for a in assets_refit_NE]))
                lat_ref,lon_ref = calin.util.utm.tm_to_geographic(ellipse.a+ref_alt,ellipse.e2,1.0,lon_ref,0,0,ref_N,ref_E)
            else:
                ref_N, ref_E = calin.util.utm.geographic_to_tm(ellipse.a+ref_alt,ellipse.e2,1.0,lon_ref,0,0,lat_ref,lon_ref)

        assets = [[a[0],a[1],a[2],a[3],a[4],a[5],a[6],ne[0]-ref_E,ne[1]-ref_N] for (a,ne) in zip(assets, assets_refit_NE) ]
    elif apply_corrections:
        _,_,_,_,_,grid_convergence_rad,scale = calin.util.utm.geographic_to_grid_with_convergence_and_scale(
            ellipse.a,ellipse.e2,lat_ref,lon_ref,utm_zone,utm_hemi)

        scale = (1 + ref_alt/ellipse.a)/scale

        ct = numpy.cos(grid_convergence_rad)
        st = numpy.sin(grid_convergence_rad)

        X = numpy.asarray([a[5] for a in assets]) - ref_E
        Y = numpy.asarray([a[6] for a in assets]) - ref_N

        X, Y = ct*X+st*Y, ct*Y-st*X

        if(demand_ref_lon is None or demand_ref_lat is None):
            ref_X = 0.5*(numpy.min(X) + numpy.max(X))
            ref_Y = 0.5*(numpy.min(Y) + numpy.max(Y))
            assets = [[a[0],a[1],a[2],a[3],a[4],a[5],a[6],(x-ref_X)*scale,(y-ref_Y)*scale] for (a,x,y) in zip(assets,X,Y)]
            ref_E, ref_N = ref_E+ct*ref_X-st*ref_Y, ref_N+ct*ref_Y+st*ref_X
            _,lat_ref,lon_ref = calin.util.utm.grid_to_geographic(ellipse.a,ellipse.e2,utm_zone,utm_hemi,ref_N,ref_E)
        else:
            ref_X = 0
            ref_Y = 0

        assets = [[a[0],a[1],a[2],a[3],a[4],a[5],a[6],(x-ref_X)*scale,(y-ref_Y)*scale] for (a,x,y) in zip(assets,X,Y)]
    else:
        assets = [[a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[5]-ref_E,a[6]-ref_N] for a in assets]

    return lat_ref,lon_ref,ref_alt,assets

def ctan_default_assets_filename():
    return 'CTAN_ArrayElements_Positions.ecsv'

def ctan_assets(filename = ctan_default_assets_filename(),
        utm_zone = calin.util.utm.UTM_ZONE_28, utm_hemi = calin.util.utm.HEMI_NORTH,
        **args):
    return load_assets(ds_filename(filename), utm_zone, utm_hemi, **args)

def ctas_assets(filename = 'CTAS_ArrayElements_Positions.ecsv',
        utm_zone = calin.util.utm.UTM_ZONE_19, utm_hemi = calin.util.utm.HEMI_SOUTH,
        **args):
    return load_assets(ds_filename(filename), utm_zone, utm_hemi, **args)

def ctan_observation_level(level_km = 2.156):
    return level_km * 1e5

def ctan_atmosphere(profile = 'ecmwf_intermediate', standard_profiles = {
            'ecmwf_intermediate': 'atmprof_ecmwf_north_intermediate_fixed.dat',
            'ecmwf_summer': 'atmprof_ecmwf_north_summer_fixed.dat',
            'ecmwf_winter': 'atmprof_ecmwf_north_winter_fixed.dat' },
        zobs = ctan_observation_level(), quiet=False):
    cfg = calin.simulation.atmosphere.LayeredRefractiveAtmosphere.default_config()
    cfg.set_angular_model_optimization_altitude(18e5)
    cfg.set_zn_reference(45)
    cfg.set_zn_optimize([5,10,15,20,25,30,35,40,50,55,60])
    cfg.set_high_accuracy_mode(True)
    profile = standard_profiles[profile] if profile in standard_profiles else profile
    atm = calin.simulation.atmosphere.LayeredRefractiveAtmosphere(ds_filename(profile), zobs, cfg)

    if(not quiet):
        print('Loading atmospheric profile :',profile)
        print('Observation level : %.3f km, thickness %.1f g/cm^2'%(atm.zobs(0)*1e-5,
            atm.thickness(atm.zobs(0))))
        print('Cherenkov angle at %.3f, 10, 20 km : %.3f, %.3f, %.3f deg'%(atm.zobs(0)*1e-5,
            numpy.arccos(1/(atm.n_minus_one(atm.zobs(0))+1))/numpy.pi*180,
            numpy.arccos(1/(atm.n_minus_one(10e5)+1))/numpy.pi*180,
            numpy.arccos(1/(atm.n_minus_one(20e5)+1))/numpy.pi*180))
        prop_ct = atm.propagation_ct_correction(atm.top_of_atmosphere()) - atm.propagation_ct_correction(atm.zobs(0))
        print('Vertical propagation delay from %d to %.3f km : %.2f ns (%.1f cm)'%(
            atm.top_of_atmosphere()*1e-5, atm.zobs(0)*1e-5, prop_ct*0.03335641, prop_ct))

    return atm

def ctan_atmospheric_absorption(absorption_model = 'navy_maritime', standard_models = {
            'low_extinction': 'atm_trans_2156_1_3_2_0_0_0.1_0.1.dat',
            'navy_maritime': 'atm_trans_2156_1_3_0_0_0.dat' },
        quiet=False):
    absorption_model = standard_models[absorption_model] if absorption_model in standard_models else absorption_model
    atm_abs = calin.simulation.detector_efficiency.AtmosphericAbsorption(ds_filename(absorption_model))

    if(not quiet):
        print('Loading atmospheric absoprtion model :',absorption_model)

    return atm_abs

def mstn_detection_efficiency(qe = 'qe_R12992-100-05b.dat',
        mirror = 'ref_MST-North-MLT_2022_06_28.dat',
        window = 'transmission_lst_window_No7-10_ave.dat',
        degradation_factor = 0.9, quiet = False):
    det_eff = calin.simulation.detector_efficiency.DetectionEfficiency()

    # PMT quantum efficiency curve
    if(qe):
        det_eff.scaleEffFromFile(ds_filename(qe))
        if(not quiet):
            print('Loading PMT QE curve :',qe,'[bandwidth = %.3f]'%det_eff.integrate())

    # Mirror reflectivity curve
    if(mirror):
        det_eff.scaleEffFromFile(ds_filename(mirror))
        if(not quiet):
            print('Scaling by mirror reflectivity :',mirror,'[bandwidth = %.3f]'%det_eff.integrate())

    # Window transmission curve
    if(window):
        det_eff.scaleEffFromFile(ds_filename(window))
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

    cone_eff = calin.simulation.detector_efficiency.AngularEfficiency(ds_filename(cone))
    if(not quiet):
        print('Loading light-cone curve :',cone)

    return cone_eff

def mstn_spe_amplitude_generator(spe = "spe_nectarcam_lmp_run1512.dat", dpdq_min=0,
        spline_mode = calin.simulation.detector_efficiency.SM_SQRT_LOG,
        spline_ninterval = 0, regularize_spline = True, extend_spline = False,
        spline_normalization = None, rescale_gain_to_unity = True, quiet = False):
    q = []
    dpdq = []
    with open(spe, 'r') as file:
        file_record = calin.provenance.chronicle.register_file_open(spe,
            calin.ix.provenance.chronicle.AT_READ, 'calin.simulation.vs_cta.mstn_spe_amplitude_generator')
        comment = ""
        for line in file.readlines():
            if(line[0]=='#'):
                comment += line
            else:
                spe_data = line.split()
                q.append(float(spe_data[0]))
                dpdq.append(float(spe_data[1]))
        file_record.set_comment(comment)
        calin.provenance.chronicle.register_file_close(file_record)

    q = numpy.asarray(q)
    dpdq = numpy.asarray(dpdq)
    m = dpdq>0

    spl = None
    if spline_normalization is None:
        spl = calin.simulation.detector_efficiency.SplinePEAmplitudeGenerator.make_spline(
            q[m], dpdq[m], spline_mode, regularize_spline, extend_spline, int(spline_ninterval))
    else:
        spl = calin.simulation.detector_efficiency.SplinePEAmplitudeGenerator.make_spline(
            q[m], dpdq[m], spline_mode, regularize_spline, extend_spline, int(spline_ninterval),
            spline_normalization)

    if rescale_gain_to_unity:
        if(spline_mode == calin.simulation.detector_efficiency.SM_LINEAR):
            def integrand(x):
                return spl.value(x)
            (I,Ierr) = scipy.integrate.quad(integrand,0,1,limit=1000,epsabs=-1,epsrel=1e-6)
        elif(spline_mode == calin.simulation.detector_efficiency.SM_LOG):
            def integrand(x):
                return spl.value(x)*numpy.exp(-x)
            (I,Ierr) = scipy.integrate.quad(integrand,0,numpy.inf,limit=1000,epsabs=-1,epsrel=1e-6)
        elif(spline_mode == calin.simulation.detector_efficiency.SM_SQRT_LOG):
            def integrand(x):
                return spl.value(x)*2*x*numpy.exp(-x**2)
            (I,Ierr) = scipy.integrate.quad(integrand,0,numpy.inf,limit=1000,epsabs=-1,epsrel=1e-6)
        spl.rescale(1/I)

    return calin.simulation.detector_efficiency.SplinePEAmplitudeGenerator(spl, spline_mode)

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

def mstn_generic_config(scope_x, scope_y, scope_z, array_lat, array_lon, array_alt,
        obscure_camera = True, include_window = False):
    mst = calin.ix.simulation.vs_optics.IsotropicDCArrayParameters()
    mst.mutable_array_origin().set_latitude(array_lon/180*numpy.pi)
    mst.mutable_array_origin().set_longitude(array_lon/180*numpy.pi)
    mst.mutable_array_origin().set_elevation(array_alt)
    try:
        for i in range(numpy.max([len(scope_x), len(scope_y), len(scope_z)])):
            scope = mst.mutable_prescribed_array_layout().add_scope_positions();
            scope.set_x(scope_x[i])
            scope.set_y(scope_y[i])
            # 900 cm offset from base to elevation axis : see https://gitlab.cta-observatory.org/cta-science/simulations/simulation-model/verification/verification-process/mst-structure/-/blob/146f5e633c7344fd3dddcd1b2f88028a75ebab47/meetings/Feb28_2022_call.md
            scope.set_z(scope_z[i] + 900)
    except:
        scope = mst.mutable_prescribed_array_layout().add_scope_positions();
        scope.set_x(scope_x)
        scope.set_y(scope_y)
        scope.set_z(mst.array_origin().elevation())
    mst.mutable_reflector_frame().set_optic_axis_rotation(-90);
    # See https://gitlab.cta-observatory.org/cta-science/simulations/simulation-model/verification/verification-process/mst-structure/-/blob/146f5e633c7344fd3dddcd1b2f88028a75ebab47/meetings/Feb28_2022_call.md
    mst.mutable_reflector_frame().mutable_translation().set_y(-175);
    # Taken from SCT (verified on MST Prod6 model)
    mst.mutable_reflector_frame().set_azimuth_elevation_axes_separation(160)
    dc = mst.mutable_reflector()
    dc.set_curvature_radius(1920)
    dc.set_aperture(1245)
    dc.set_facet_num_hex_rings(5)
    dc.mutable_psf_align().set_object_plane(numpy.inf) # 10 * 1e5);
    dc.set_alignment_image_plane(1600)
    dc.set_facet_spacing(124)
    dc.set_facet_size(120)
    dc.set_facet_focal_length(1608.3)
    dc.set_facet_focal_length_dispersion(3.9)
    dc.set_facet_spot_size_probability(0.8)
    # Updated 2022-08-10 so on-axis PSF could match value presented by DESY :
    # https://indico.cta-observatory.org/event/3956/contributions/32824/attachments/20968/29502/220228_MST-PSF-MC-Feedback.pdf
    # dc.set_facet_spot_size(0.5 * 2.8) # Spot size of 28mm at 2F
    dc.set_facet_spot_size(0.5 * 2.0)
    dc.set_facet_spot_size_dispersion(0.5 * 0.02)
    dc.set_facet_labeling_parity(True)
    dc.set_weathering_factor(1.0)
    #for id in [1,62,67,72,77,82,87]: dc.add_facet_missing_list(id-1) # 84 mirror
    for id in [1,67,72,82,87]: dc.add_facet_missing_list(id-1) # 86 mirror
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

    if obscure_camera:
        obs_camera_box = mst.add_pre_reflection_obscuration()
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_x(2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_z(2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_max_corner().set_y(mst.focal_plane().translation().y()+(1532.25-513.0)/10)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_x(-2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_z(-2918.0/20)
        obs_camera_box.mutable_aligned_box().mutable_min_corner().set_y(mst.focal_plane().translation().y()-513.0/10)

        # obs_outer_aperture = mst.add_camera_obscuration()
        # obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_x(0)
        # obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_z(-50.0/20)
        # obs_outer_aperture.mutable_rectangular_aperture().mutable_center_pos().set_y(-513.0/10)
        # obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_x(2714.0/10)
        # obs_outer_aperture.mutable_rectangular_aperture().set_flat_to_flat_z((2585.0-50.0)/10)
        #
        # obs_inner_aperture = mst.add_camera_obscuration()
        # obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_x(0)
        # obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_z(0)
        # obs_inner_aperture.mutable_circular_aperture().mutable_center_pos().set_y(-222.5/10)
        # obs_inner_aperture.mutable_circular_aperture().set_diameter(2304.0/10)

    if include_window:
        win = mst.mutable_spherical_window()
        win.set_front_y_coord(mst.focal_plane().translation().y() - 427.5/10)
        win.set_outer_radius(3443.0/10)
        win.set_thickness(5.25/10.0)
        win.set_refractive_index(1.5)

    return mst

def mstn1_config(obscure_camera = True, assets_file = ctan_default_assets_filename(),
        include_window = False):
    array_lat,array_lon,array_alt,all_assets = ctan_assets(filename = assets_file)
    scope_x = []
    scope_y = []
    scope_z = []
    for a in all_assets:
        if(a[0]=='MSTN' and a[1]=='03'):
            scope_x.append(a[7]*100)
            scope_y.append(a[8]*100)
            scope_z.append(a[4]*100)
    return mstn_generic_config(scope_x, scope_y, scope_z, array_lat, array_lon, array_alt*100,
        obscure_camera = obscure_camera, include_window = include_window)

def mstn_config(obscure_camera = True, assets_file = ctan_default_assets_filename(),
        include_window = False):
    array_lat,array_lon,array_alt,all_assets = ctan_assets(filename = assets_file)
    scope_x = []
    scope_y = []
    scope_z = []
    for a in all_assets:
        if(a[0]=='MSTN'):
            scope_x.append(a[7]*100)
            scope_y.append(a[8]*100)
            scope_z.append(a[4]*100)
    return mstn_generic_config(scope_x, scope_y, scope_z, array_lat, array_lon, array_alt*100,
        obscure_camera = obscure_camera, include_window = include_window)

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

def make_array(cfg, rng = None):
    rng = calin.math.rng.RNG('calin.simulation.vs_cta.make_array') if rng is None else rng
    array = calin.simulation.vs_optics.VSOArray()
    if(type(cfg) is list):
        for icfg in cfg:
            array.generateFromArrayParameters(icfg, rng)
    else:
        array.generateFromArrayParameters(cfg, rng)
    return array
