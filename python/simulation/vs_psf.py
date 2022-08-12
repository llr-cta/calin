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
import calin.simulation.ray_processor
import calin.simulation.vs_optics
import calin.simulation.vs_cta
import calin.math.geometry
import scipy.optimize
import calin.provenance.system_info

def calc_psf(vs_scope, theta_deg, phi_deg, refocus_at_infinity = False, dx=0.1, image_resolution=0.01, d80_frac=0.8):
    scope_cfg = vs_scope.dump_as_proto()
#     print(scope_cfg.DebugString())
    if(refocus_at_infinity):
        raise NotImplementedError("refocus_at_infinity not yet implemented")
        pass
    array_cfg = calin.ix.simulation.vs_optics.VSOArrayData()
    array_cfg.add_telescope().CopyFrom(scope_cfg)
    array = calin.simulation.vs_optics.VSOArray.create_from_proto(array_cfg)
    array.pointTelescopesAzEl(0/180*numpy.pi,90/180*numpy.pi)
    scope = array.telescope(0)

    ns_cfg = calin.simulation.ray_processor.NSpacePEProcessor.default_config()
    ns_cfg.set_axis_variables(calin.ix.simulation.ray_processor.XY)
    ns_cfg.set_xy_diameter(max(scope_cfg.camera_diameter(), 2*theta_deg*scope_cfg.fp_translation().y()/180*numpy.pi)*1.2)
    ns_cfg.set_xy_num_bins(int(ns_cfg.xy_diameter()/(image_resolution*scope.pixelSpacing())))

    pe_imager = calin.simulation.ray_processor.NSpacePEProcessor(ns_cfg)

    if(calin.provenance.system_info.vcl_uses_avx512f()):
        processor = calin.simulation.ray_processor.VCLRayTracerRayProcessorDouble512(array, pe_imager)
    elif(calin.provenance.system_info.vcl_uses_avx2()):
        processor = calin.simulation.ray_processor.VCLRayTracerRayProcessorDouble256(array, pe_imager)
    elif(True):
        processor = calin.simulation.ray_processor.VCLRayTracerRayProcessorDouble128(array, pe_imager)
    else:
        processor = calin.simulation.ray_processor.VCVSORayProcessor(cta, pe_imager)

    R = scope.reflectorIP()/2.0

    x00 = scope.reflectorIPCenter()
    scope.reflectorToGlobal_pos(x00)

    pos_gen = calin.math.ray_generator.HexGridPlanePositionGenerator(R, dx)
    dir_gen = calin.math.ray_generator.SingleDirectionGenerator()
    ray_gen = calin.math.ray_generator.PositionNestedByDirectionRayGenerator(x00,
        calin.math.geometry.rotation_theta_phi((180-theta_deg)/180*numpy.pi,phi_deg/180*numpy.pi),
        numpy.asarray([0,0,-100e5]), dir_gen, pos_gen)

    nray = processor.process_all_from_ray_generator(ray_gen)
    nhit = processor.nhit()

    ns = pe_imager.nspace()
    m2,m1,m0 = ns.covar_mean_and_total_weight()

    xc = ns.axis_bin_centers(0)
    yc = ns.axis_bin_centers(1)
    ns_vals = ns.as_matrix().transpose()

    ns_vals_px = numpy.sum(ns_vals, axis=0)
    ix0 = numpy.min(numpy.where(ns_vals_px>0))
    ix1 = numpy.max(numpy.where(ns_vals_px>0))

    ns_vals_py = numpy.sum(ns_vals, axis=1)
    iy0 = numpy.min(numpy.where(ns_vals_py>0))
    iy1 = numpy.max(numpy.where(ns_vals_py>0))

    xc = xc[ix0:ix1+1]
    yc = yc[iy0:iy1+1]
    Xc, Yc = numpy.meshgrid(xc,yc)

    ns_vals = ns_vals[iy0:iy1+1,ix0:ix1+1]
    ns_vals = ns_vals.flatten()

    def d80_for_xy(xy0):
        x0,y0 = xy0

        X = Xc - x0
        Y = Yc - y0

        d80x = numpy.sqrt(X.flatten()**2 + Y.flatten()**2)
        isort = numpy.argsort(d80x)
        d80x = d80x[isort]
        d80y = numpy.cumsum(ns_vals[isort])/numpy.sum(ns_vals)

        d80 = numpy.interp(d80_frac, d80y, d80x)

        return d80

    opt = scipy.optimize.minimize(d80_for_xy, m1, method='Nelder-Mead')

    return m0,m1,m2,opt.fun*2,opt.x,nray,nhit
