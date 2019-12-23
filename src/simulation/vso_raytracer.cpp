/*

   calin/simulation/vso_raytracer.cpp -- Stephen Fegan -- 2015-11-30

   Class for raytracing on a VSOTelescope

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

//-*-mode:c++; mode:font-lock;-*-

/*! \file RayTracer.cpp

  Ray tracing class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
           sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz             \n
           UCLA                        \n
           nicewicz@physics.ucla.edu   \n

  \date    12/05/2004
  \version 0.2
  \note
*/

#include <cmath>
#include <iostream>
#include <algorithm>

#include <math/hex_array.hpp>
#include <math/vector3d_util.hpp>
#include <simulation/vso_raytracer.hpp>

using namespace calin::simulation::vs_optics;

void VSOTraceInfo::reset()
{
  *this = VSOTraceInfo();
}

std::ostream& VSOTraceInfo::write(std::ostream& stream, bool cpu, bool eol) const
{
  stream
      << status << ' '                                                    // $1  --  AWK column
      << ground_x << ' '                                                  // $2
      << ground_y << ' '                                                  // $3
      << ground_dx << ' '                                                 // $4
      << ground_dy << ' '                                                 // $5
      << scope_id << ' '                                                  // $6
      << scope_id << ' '                                                  // $7
      << reflec_x << ' '                                                  // $8
      << reflec_z << ' '                                                  // $9
      << reflec_dx << ' '                                                 // $12
      << reflec_dz << ' '                                                 // $13
      << mirror_hexid << ' '                                              // $14
      << ( mirror ? mirror->hexID() : 0 ) << ' '                          // $15
      << mirror_x << ' '                                                  // $16
      << mirror_y << ' '                                                  // $17
      << mirror_z << ' '                                                  // $18
      << mirror_normal.transpose() << ' '                 // $19 $20 $21 $22 $23
      << mirror_normal_dispersion << ' '                                  // $24
      << mirror_scattered.transpose() << ' '              // $25 $26 $27 $28 $29
      << mirror_reflection_angle << ' '                                   // $30
      << fplane_x*(scope&&cpu?scope->pixelSpacing():1.0) << ' '           // $31
      << fplane_z*(scope&&cpu?scope->pixelSpacing():1.0) << ' '           // $32
      << fplane_dx*(scope&&cpu?scope->pixelSpacing():1.0) << ' '          // $33
      << fplane_dz*(scope&&cpu?scope->pixelSpacing():1.0) << ' '          // $34
      << fplane_t << ' '                                                  // $35
      << fplane_uy << ' '                                                  // $36
      << pixel_hexid << ' '                                               // $37
      << ( pixel ? pixel->hexID() : 0 ) << ' '                            // $38
      << pixel_dist*(scope&&cpu?scope->pixelSpacing():1.0) << ' '         // $39
      << concentrator_hit;                                                // $40
  if(eol)stream << '\n';
  return stream;
}


std::ostream& operator << (std::ostream& stream, const VSOTraceInfo& o)
{
  return o.write(stream,false);
}


VSORayTracer::~VSORayTracer()
{
  // nothing to see here
}

const VSOPixel* VSORayTracer::
trace(math::ray::Ray& ray, TraceInfo& info, const VSOTelescope* scope_hint)
{
  // Initialize array
  info.reset();
  info.array = fArray;
  info.scope = scope_hint;
  return(scope_trace(ray,info));
}

const VSOPixel* VSORayTracer::
trace(math::ray::Ray& ray, TraceInfo& info)
{
  // Initialize array
  info.reset();
  info.array = fArray;

  // Require photon to be down going
  if(ray.direction().z()>0)
  {
    info.status = TS_DOES_INTERSECT_GROUND;
    return 0;
  }

  // Propagate to ground
  math::ray::Ray ray_copy(ray);
  bool good;
  good = ray_copy.propagate_to_plane(Eigen::Vector3d::UnitZ(),
                                      -fArray->altitude(), true);
  if(!good)
  {
    info.status = TS_DOES_INTERSECT_GROUND;
    return 0;
  }

  //#warning Need to implement a real algorithm for finding which telescope was hit

  // Array is assumed to be hexagonal use VVV look-up function to find site
  info.ground_x = ray_copy.position().x();
  info.ground_y = ray_copy.position().y();
  info.ground_dx = 0;
  info.ground_dy = 0;
  info.scope_id  = 0;

  // Find telescope (if there is a real telescope at that site)
  info.scope = fArray->telescope(info.scope_id);
  if(info.scope==0)
  {
    info.status = TS_NO_SCOPE;
    return 0;
  }

#if 0
  // Propagate to reflector impact sphere -- test whether ray can hit scope
  // This probably does not speed things up much as the propagation to
  // the reflector (two steps down) does a similar thing
  good = ray_copy.PropagateFreeToSphere(info.scope->position(),
					info.scope->reflectorIP(),
					IP_EARLIEST,true);
  if(!good)
  {
    info.status = TS_OUTSIDE_REFLECTOR_IP;
    return 0;
  }
#endif

  return(scope_trace(ray,info));
}

const VSOPixel*
VSORayTracer::scope_trace(math::ray::Ray& ray, TraceInfo& info)
{
  bool good;
  double ref_index = 1.0;

  // Transform to reflector co-ordinate system
  info.scope->globalToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "A: " << ray.direction().transpose() << std::endl;
#endif
  // **************************************************************************
  // ****************** RAY IS NOW IN RELECTOR COORDINATES *******************
  // **************************************************************************

  // Test for obscuration of incoming ray
  unsigned nobs = info.scope->numPreReflectionObscurations();
  unsigned obs_ihit  = nobs;
  double   obs_time  = ray.ct();
  for(unsigned iobs=0;iobs<nobs;iobs++)
  {
    math::ray::Ray r_out;
    if(info.scope->pre_reflection_obscuration(iobs)->doesObscure(ray, r_out, ref_index))
    {
      if((obs_ihit==nobs)||(r_out.ct()<obs_time))
        obs_ihit = iobs, obs_time = r_out.ct();
    }
  }

  // Propagate to intersection with the reflector sphere
  good = ray.propagate_to_standard_sphere_2nd_interaction_fwd_only(
    info.scope->curvatureRadius(), ref_index);

  if(!good)
  {
    info.scope->reflectorToGlobal(ray);
    info.status = TS_MISSED_REFLECTOR_SPHERE;
    return 0;
  }

  // Assume mirrors on hexagonal grid - use hex_array routines to find which hit
  info.reflec_x     = ray.position().x();
  info.reflec_z     = ray.position().z();
  info.reflec_dx    = info.reflec_x;
  info.reflec_dz    = info.reflec_z;
  info.mirror_hexid =
    math::hex_array::xy_trans_to_hexid_with_remainder(
      info.reflec_dx, info.reflec_dz, info.scope->mirrorParity(),
      info.scope->cosReflectorRotation(), info.scope->sinReflectorRotation(),
      info.scope->facetSpacing(),
      info.scope->facetGridShiftX(), info.scope->facetGridShiftZ());

  // Find mirror - searching neighbors if desired
  good = findMirror(ray, info, ref_index);
  if(!good)
  {
    // info.status should be set appropriately by findMirror
    info.scope->reflectorToGlobal(ray);
    return 0;
  }

  // **************************************************************************
  // ******************** RAY IS NOW IN MIRROR COORDINATES ********************
  // **************************************************************************

  // Test for interaction with obscuration before mirror was hit
  if((obs_ihit!=nobs) && (obs_time<ray.ct()))
  {
    info.status = TS_OBSCURED_BEFORE_MIRROR;
    info.mirror->mirrorToReflector(ray);
    info.scope->reflectorToGlobal(ray);
    info.obscuration = info.scope->pre_reflection_obscuration(obs_ihit);
    info.obscuration_id = obs_ihit;
    return 0;
  }

  // Check to see if photon is absorbed at the mirror. Would be faster to
  // do this check before the edge check, but this way gets info.status
  // correct .. is this important ?
  if(info.mirror->degradingFactor() < 1.0 and
    fRNG->uniform() > info.mirror->degradingFactor())
  {
    info.status = TS_ABSORBED_AT_MIRROR;
    info.mirror->mirrorToReflector(ray);
    info.scope->reflectorToGlobal(ray);
    return 0;
  }

  // Find normal at the point of intersection
  double mirror_radius = info.mirror->focalLength()*2.0;
  info.mirror_normal =
    (Eigen::Vector3d(0,mirror_radius,0) - ray.position()).normalized();

  // Scatter the normal to account for the spot size ot the focal length of the
  // radius. The spot size is given as the DIAMETER at the focal distance.
  // Must divide by 2 (for reflection) and another 2 for diameter -> radius
  info.mirror_normal_dispersion =
    0.25*info.mirror->spotSize()/info.mirror->focalLength();

  info.mirror_scattered = info.mirror_normal;
  calin::math::vector3d_util::scatter_direction(info.mirror_scattered,
    info.mirror_normal_dispersion,*fRNG);

  // Reflect ray
  ray.reflect_from_surface(info.mirror_scattered);

  // Back to reflector coordinates
  info.mirror->mirrorToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "C: " << ray.direction().transpose() << std::endl;
#endif

  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************


  // Bending at window, if defined
  if(info.scope->windowThickness() > 0)
  {
    const double n = info.scope->windowRefractiveIndex()/ref_index;
    Eigen::Vector3d er;
    if(info.scope->windowOuterRadius() > 0) {
      // Spherical window
      good = ray.propagate_to_y_sphere_1st_interaction_fwd_only(
        info.scope->windowOuterRadius(), info.scope->windowFront(), ref_index);
      er = ray.position();
      er.y() -= info.scope->windowOuterRadius()+info.scope->windowFront();
      er *= 1.0/info.scope->windowOuterRadius();
    } else {
      good = ray.propagate_to_y_plane(-info.scope->windowFront(), ref_index);
      er << 0, -1, 0;
    }

    if(!good)
    {
      info.status = TS_MISSED_WINDOW;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

    ray.refract_at_surface_in(er, n);

    if(info.scope->windowOuterRadius() > 0) {
      // Spherical window
      good = ray.propagate_to_y_sphere_1st_interaction_fwd_only(
        info.scope->windowOuterRadius()-info.scope->windowThickness(),
        info.scope->windowFront()+info.scope->windowThickness(),
        info.scope->windowRefractiveIndex());
      er = ray.position();
      er.y() -= info.scope->windowOuterRadius() + info.scope->windowFront();
      er *= -1.0/(info.scope->windowOuterRadius()-info.scope->windowThickness());
    } else {
      good = ray.propagate_to_y_plane(-info.scope->windowFront(), ref_index);
      er << 0, 1, 0;
    }

    if(!good)
    {
      // Case of internal surface not encountered - ray exits outer again
      info.status = TS_MISSED_WINDOW;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

    good = ray.refract_at_surface_out(er, n);

    if(!good)
    {
      // Total internal reflection not supported
      info.status = TS_MISSED_WINDOW;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }
  }

  // Test for interaction with obscuration in reflector coordinates
  nobs = info.scope->numPostReflectionObscurations() + info.scope->numCameraObscurations();
  obs_ihit  = nobs;
  obs_time  = std::numeric_limits<double>::infinity();
  for(unsigned iobs=0;iobs<info.scope->numPostReflectionObscurations();iobs++)
  {
    math::ray::Ray r_out;
    if(info.scope->post_reflection_obscuration(iobs)->doesObscure(ray, r_out, ref_index))
    {
      if((obs_ihit==nobs)||(r_out.ct()<obs_time))
        obs_ihit = iobs, obs_time = r_out.ct();
    }
  }

  // Translate to focal plane coordinates
  info.scope->reflectorToFocalPlane(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "D: " << ray.direction().transpose() << std::endl;
#endif

  // **************************************************************************
  // ***************** RAY IS NOW IN FOCAL PLANE COORDINATES ******************
  // **************************************************************************

  // Test for interaction with obscuration in camera coordinates
  for(unsigned iobs=0;iobs<info.scope->numCameraObscurations();iobs++)
  {
    math::ray::Ray r_out;
    if(info.scope->camera_obscuration(iobs)->doesObscure(ray, r_out, ref_index))
    {
      if((obs_ihit==nobs)||(r_out.ct()<obs_time))
        obs_ihit = iobs+info.scope->numPostReflectionObscurations(), obs_time = r_out.ct();
    }
  }

  // Propagate back to camera plane
  good = ray.propagate_to_y_plane(0,false,ref_index);
  if(!good)
  {
    info.status = TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE;
    info.scope->focalPlaneToGlobal(ray);
    return 0;
  }


  if((obs_ihit!=nobs) && (obs_time<ray.ct()))
  {
    info.status = TS_OBSCURED_BEFORE_FOCAL_PLANE;
    info.scope->focalPlaneToReflector(ray);
    info.scope->reflectorToGlobal(ray);
    if(obs_ihit < info.scope->numPostReflectionObscurations())
      info.scope->post_reflection_obscuration(obs_ihit);
    else
      info.scope->post_reflection_obscuration(obs_ihit - info.scope->numPostReflectionObscurations());
    info.obscuration = info.scope->post_reflection_obscuration(obs_ihit);
    info.obscuration_id = obs_ihit;
    return 0;
  }

  // We good, record position on focal plane etc
  info.fplane_x = ray.position().x();
  info.fplane_z = ray.position().z();
  info.fplane_dx = info.fplane_x;
  info.fplane_dz = info.fplane_z;
  info.fplane_t = ray.ct() / math::constants::cgs_c;
  info.fplane_uy = ray.direction().y();

  info.pixel_hexid =
    math::hex_array::xy_trans_to_hexid_with_remainder(
      info.fplane_dx, info.fplane_dz, info.scope->pixelParity(),
      info.scope->cosPixelRotation(), info.scope->sinPixelRotation(),
      info.scope->pixelSpacing(),
      info.scope->pixelGridShiftX(), info.scope->pixelGridShiftZ());

  // Find pixel (if there is a real pixel at that site)
  info.pixel = info.scope->pixelByHexID(info.pixel_hexid);
  if(info.pixel==0)
  {
    info.status = TS_NO_PIXEL;
    info.scope->focalPlaneToGlobal(ray);
    return 0;
  }

  info.pixel_dist =
    sqrt(info.fplane_dx*info.fplane_dx + info.fplane_dz*info.fplane_dz);

  info.concentrator_hit =
      (info.pixel_dist > info.scope->cathodeDiameter()/2.0);
  if(info.concentrator_hit)
  {
    if(info.scope->concentratorSurvivalProb() < 1.0 and
      fRNG->uniform() > info.scope->concentratorSurvivalProb())
    {
      info.status = TS_ABSORBED_AT_CONCENTRATOR;
      info.scope->focalPlaneToGlobal(ray);
      return 0;
    }
  }

  // Translate to global coordinates
  info.scope->focalPlaneToGlobal(ray);

  // **************************************************************************
  // ******************** RAY IS NOW IN GLOBAL COORDINATES ********************
  // **************************************************************************

  info.status = TS_PE_GENERATED;
  return info.pixel;
}

// Find the mirror (allowing that it might not actually be one that
// neighbors the optimal position). Starts with ray at the reflector
// sphere, in the reflector coordinate system and info.mirror_hexid
// set to "nominal" mirror given position on reflector. Ends with
// either the ray in the mirror system (with info.mirror,
// mirror_hexid, mirror_x, mirror_y, mirror_z set correctly), in which
// case "true" is returned, or the mirror still in reflector
// coordinates, in the case where no mirror is encountered, where
// false is returned

bool VSORayTracer::
findMirror(math::ray::Ray& ray, TraceInfo& info, double ref_index)
{
  // **************************************************************************
  // ****************** RAY STARTS IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  constexpr unsigned search_radius = 0;
  constexpr unsigned nsearch = 3*search_radius*(search_radius+1)+1;

  math::ray::Ray ray_in(ray);

#if 0
  int               info_hexid  = info.mirror_hexid;
  Status            info_status = info.status;
  const VSOMirror*  info_mirror;
  bool              info_mirror_x;
  bool              info_mirror_y;
  bool              info_mirror_z;
  Particle          info_ray;
#endif
  bool              impinging_ray_found = false;

  unsigned hex_i0 = info.mirror_hexid;
  int hex_u0 = 0;
  int hex_v0 = 0;
  if(nsearch > 1)math::hex_array::hexid_to_uv(hex_i0, hex_u0, hex_v0);

  for(unsigned isearch = 0; isearch<nsearch; isearch++)
  {
    int test_hexid;
    if(isearch == 0)
    {
      test_hexid = info.mirror_hexid;
    }
    else
    {
      int hex_i1 = isearch;
      int hex_u1 = 0;
      int hex_v1 = 0;
      math::hex_array::hexid_to_uv(hex_i1, hex_u1, hex_v1);
      hex_u1 += hex_u0;
      hex_v1 += hex_v0;
      test_hexid = math::hex_array::uv_to_hexid(hex_u1, hex_v1);
    }

    // Find mirror (if there is a real mirror at that site)
    const VSOMirror* test_mirror = info.scope->mirrorByHexID(test_hexid);
    if(isearch==0)info.mirror = test_mirror;

    if(test_mirror==0)
    {
      if(isearch==0)info.status = TS_NO_MIRROR;
      continue;
    }

    if(test_mirror->removed())
    {
      if(isearch==0)info.status = TS_MIRROR_REMOVED;
      continue;
    }

    // Convert to mirror coordinates
    math::ray::Ray test_ray(ray_in);
    test_mirror->reflectorToMirror(test_ray);

    // **********************************************************************
    // ****************** RAY IS NOW IN MIRROR COORDINATES ******************
    // **********************************************************************

    // Propagate to intersection with the mirror sphere
    double mirror_radius = test_mirror->focalLength()*2.0;

    bool good;
    good = test_ray.
      propagate_to_standard_sphere_2nd_interaction_fwd_bwd(mirror_radius, ref_index);
    if(!good)
    {
      if(isearch==0)info.status = TS_MISSED_MIRROR_SPHERE;
      continue;
    }

    double test_mirror_x = test_ray.position().x();
    double test_mirror_y = test_ray.position().y();
    double test_mirror_z = test_ray.position().z();

    if(isearch==0)
    {
      info.mirror_x  = test_mirror_x;
      info.mirror_y  = test_mirror_y;
      info.mirror_z  = test_mirror_z;
    }

    // Check if ray impacted beyond the edge of this mirror
    constexpr double cos60 = 0.5;
    constexpr double sin60 = 0.5*CALIN_HEX_ARRAY_SQRT3;
    const double edge = 0.5*info.scope->facetSize();

    const double x_pos60 = std::fabs(cos60*test_mirror_x - sin60*test_mirror_z);
    const double x_neg60 = std::fabs(cos60*test_mirror_x + sin60*test_mirror_z);

    double xmax = std::max(x_neg60, x_pos60);
    xmax = std::max(xmax, std::fabs(test_mirror_x));

    if(xmax > edge)
    {
      if(isearch==0)info.status = TS_MISSED_MIRROR_EDGE;
      continue;
    }

    // We have a good ray that impinges on a mirror - test that it hits
    // earlier than previously found one (if any)
    if((isearch>0)
       &&((!impinging_ray_found)
          ||(test_ray.ct()<ray.ct())))
    {
      info.mirror_hexid = test_hexid;
      info.status       = info.status;
      info.mirror       = test_mirror;
      info.mirror_x     = test_mirror_x;
      info.mirror_y     = test_mirror_y;
      info.mirror_z     = test_mirror_z;
      ray               = test_ray;
    }

    if(isearch == 0)ray = test_ray;
    impinging_ray_found = true;
  }

  return impinging_ray_found;
}

bool VSORayTracer::beam(math::ray::Ray& photon,
                        const Eigen::Vector3d& origin,
                        const Eigen::Vector3d& direction,
                        double beam_start, double beam_stop,
                        double beam_radius_in, double beam_radius_out,
                        double beam_angle_lo, double beam_angle_hi,
                        double energy_ev)
{
  Eigen::Vector3d d_hat = direction.normalized();

  Eigen::Matrix3d rot =
    Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(),d_hat).toRotationMatrix();

  Eigen::Vector3d tangent_a = rot * Eigen::Vector3d::UnitX();
  Eigen::Vector3d tangent_b = rot * Eigen::Vector3d::UnitY();

  // CHOOSE PHOTON EMISSION POINT
  Eigen::Vector3d emission_point(origin);
  emission_point += d_hat*beam_start;

  // SAMPLE PHOTON EMISSION POINT FROM BEAM LENGTH
  if(beam_stop != beam_start)
  {
    emission_point += d_hat*((beam_stop-beam_start)*fRNG->uniform());
  }

  // SAMPLE PHOTON EMISSION POINT FROM BEAM AREA
  if((beam_radius_in != 0)||(beam_radius_out != 0))
  {
    double theta = 2.0*M_PI*fRNG->uniform();
    double rho2 = beam_radius_in*beam_radius_in;
    if(beam_radius_in != beam_radius_out)
      rho2 += (beam_radius_out*beam_radius_out-rho2) * fRNG->uniform();
    double rho = std::sqrt(rho2);
    emission_point += tangent_a*rho*std::cos(theta);
    emission_point += tangent_b*rho*std::sin(theta);
  }

  // SAMPLE PHOTON EMISSION DIRECTION
  double costheta = cos(beam_angle_lo);
  if(beam_angle_lo != beam_angle_hi)
  {
    costheta += fRNG->uniform() * (std::cos(beam_angle_hi)-costheta);
  }
  double theta = std::acos(costheta);

  if(theta != 0)
  {
    const double phi = fRNG->uniform() * 2.0*M_PI;
    const double sintheta = std::sin(theta);
    d_hat = rot *
      Eigen::Vector3d(sintheta*std::cos(phi),sintheta*std::sin(phi), costheta);
  }

  // MAKE PHOTON
  photon = math::ray::Ray(emission_point, d_hat, 0.0, energy_ev);
  return true;
}

bool VSORayTracer::laserBeam(math::ray::Ray& photon,
                             const Eigen::Vector3d& origin,
                             const Eigen::Vector3d& direction,
                             double d0, double sampling_radius,
                             double energy_ev)
{
  return beam(photon, origin, direction, d0, d0, 0, sampling_radius, 0, 0,
	      energy_ev);
}

bool VSORayTracer::fanBeam(math::ray::Ray& photon,
                           const Eigen::Vector3d& origin,
                           const Eigen::Vector3d& direction,
                           double half_angle_spread,
                           double energy_ev)
{
  return beam(photon, origin, direction, 0, 0, 0, 0, 0, half_angle_spread,
	      energy_ev);
}

bool VSORayTracer::muonBeam(math::ray::Ray& photon,
                            const Eigen::Vector3d& origin,
                            const Eigen::Vector3d& direction,
                            double muon_travel_distance, double opening_angle,
                            double energy_ev)
{
  return beam(photon, origin, direction, 0, muon_travel_distance, 0, 0,
	      opening_angle, opening_angle, energy_ev);
}

bool VSORayTracer::testBeam(math::ray::Ray& photon,
			    const VSOTelescope* scope,
			    double theta, double phi,
			    double U, double energy_ev)
{
  if(std::isinf(U))
  {
    // Parallel beam - must set sampling radius
    Eigen::Vector3d
      beam_dir(-sin(theta)*sin(phi) ,-cos(theta), -sin(theta)*cos(phi));
    Eigen::Vector3d beam_cen(0.0,0.0,0.0);
    beam_cen += scope->reflectorIPCenter();
    math::ray::Ray ray(beam_cen, beam_dir);
    scope->reflectorToGlobal(ray);
    return laserBeam(photon, ray.position(), ray.direction(),
                     -2.0*scope->focalPlanePosition().y(),
                     0.5*scope->reflectorIP(), energy_ev);
  }
  else
  {
    double Utantheta = U*tan(theta);
    Eigen::Vector3d
        beam_cen(Utantheta*sin(phi), U, Utantheta*cos(phi));
    double dist = beam_cen.norm();
    Eigen::Vector3d beam_dir(beam_cen*(-1.0/dist));
    beam_cen += scope->reflectorIPCenter();
    math::ray::Ray ray(beam_cen, beam_dir);
    scope->reflectorToGlobal(ray);
    return fanBeam(photon, ray.position(), ray.direction(),
                   asin(0.5*scope->reflectorIP()/dist), energy_ev);
  }
}

void VSORayTracer::calcPSF(class VSOPSFInfo& psf, const VSOTelescope* scope,
			   double theta, double phi, double U, unsigned nsim,
			   bool save_image)
{
  double PS = scope->pixelSpacing()/scope->focalPlanePosition().y()/M_PI*180.0;
  psf.reset();
  math::ray::Ray ph;
  VSOTraceInfo info;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> t;
  std::vector<double> r2;
  x.reserve(nsim);
  y.reserve(nsim);
  t.reserve(nsim);
  double cphi = std::cos(phi);
  double sphi = std::sin(phi);
  unsigned nhit = 0;
  double sumx   = 0;
  double sumy   = 0;
  double sumt   = 0;
  double sumxx  = 0;
  double sumyy  = 0;
  double sumxy  = 0;
  double sumtt  = 0;
  for(unsigned isim=0;isim<nsim;isim++)
  {
    testBeam(ph, scope, theta, phi, U);
    trace(ph, info, scope);
    info.write(std::cerr);
    std::cerr << std::endl;
    if(!info.rayHitFocalPlane())continue;
    double _x = -PS*(cphi*info.fplane_z+sphi*info.fplane_x);
    double _y = -PS*(cphi*info.fplane_x-sphi*info.fplane_z);
    double _t = info.fplane_t;
    double x2 = _x*_x;
    double y2 = _y*_y;
    nhit++;
    sumx  += _x;
    sumy  += _y;
    sumt  += _t;
    sumxx += x2;
    sumxy += _x*_y;
    sumyy += y2;
    sumtt += _t*_t;
    x.push_back(_x);
    y.push_back(_y);
    t.push_back(_t);
  }

  psf.nhit        = nhit;

  if(save_image)
  {
    psf.r_tan   = x;
    psf.r_sag   = y;
    psf.t       = t;
  }

  psf.mean_tan    = sumx/double(nhit);
  psf.mean_sag    = sumy/double(nhit);
  psf.mean_t      = sumt/double(nhit);
  psf.rms_tan     = std::sqrt(sumxx/double(nhit) - psf.mean_tan*psf.mean_tan);
  psf.rms_sag     = std::sqrt(sumyy/double(nhit) - psf.mean_sag*psf.mean_sag);
  psf.cov_tan_sag = sumxy/double(nhit) - psf.mean_tan*psf.mean_sag;
  psf.rms_t       = std::sqrt(sumtt/double(nhit) - psf.mean_t*psf.mean_t);

  r2.resize(x.size());
  for(unsigned i=0;i<x.size(); i++)
  {
    double _x = x[i] - psf.mean_tan;
    double _y = y[i] - psf.mean_sag;
    r2[i] = _x*_x + _y*_y;
  }

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  std::sort(t.begin(), t.end());

  unsigned imed = x.size()/2;
  psf.median_tan  = x[imed];
  psf.median_sag  = y[imed];
  psf.median_t    = t[imed];

  for(unsigned i=0;i<x.size(); i++)
  {
    t[i] = std::fabs(t[i] - psf.mean_t);
  }

  std::sort(t.begin(), t.end());
  std::sort(r2.begin(), r2.end());

  unsigned i80 = (x.size()*8)/10;
  psf.r80         = std::sqrt(r2[i80]);
  psf.t80         = t[i80];
}
