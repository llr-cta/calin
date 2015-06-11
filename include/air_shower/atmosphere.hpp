/* 

   calin/air_shower/atmosphere.hpp -- Stephen Fegan -- 2015-06-11

   Classes to model density and refractive index of atmosphere

*/

// Originally: Atmposphere.hpp - Classes to handle atmosphere
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// Id: Atmosphere.hpp 5416 2013-06-24 13:46:00Z sfegan

#pragma once

#include<vector>
#include<cassert>

#include"package_wide_definitions.hpp"

// Units:
// Height, z:    cm
// Density, rho: g/cm^3
// Thickness, t: g/cm^2

namespace calin { namespace air_shower { namespace atmosphere {

class AtmSlice
{
public:
  AtmSlice(double _v=0): zb(_v), zt(_v), tb(_v), tt(_v), rho(_v), 
			 n_minus_one(_v) { }
  double zb;
  double zt;
  double tb;
  double tt;
  double rho;
  double n_minus_one;

  bool operator< (const AtmSlice& o) const { return zt < o.zt; }

#ifndef SWIG
  class CmpZAsc
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b) 
    { return a.zt < b.zt; } };
  
  class CmpTDec
  { public: bool operator() (const AtmSlice& a, const AtmSlice& b) 
    { return b.tb < a.tb; } };
#endif
};

// Base class for all atmospheric models
class Atmosphere
{
 public:
  virtual ~Atmosphere();
  virtual double rho(double z) = 0;
  virtual double thickness(double z) = 0;
  virtual double n_minus_one(double z) = 0;
  virtual double propagation_time_correction(double z) = 0;
  virtual double z_for_thickness(double t) = 0;
  virtual double top_of_atmosphere() = 0;
  std::vector<AtmSlice> make_atm_slices(unsigned nslice, 
                                        double zmax, double zmin);
  std::vector<AtmSlice> make_atm_slices(unsigned nslice);
};

// Simple (i.e. unrealistic) isothermal atmosphere
class IsothermalAtmosphere: public Atmosphere
{
 public:
  IsothermalAtmosphere(double rho0=1.2e-3, double zs = 8.5e5, 
		       double zmax = 1.2e7, double nmo0=2.75e-4);
  virtual ~IsothermalAtmosphere();
  virtual double rho(double z) override;
  virtual double thickness(double z) override;
  virtual double n_minus_one(double z) override;
  virtual double propagation_time_correction(double z) override;
  virtual double z_for_thickness(double t) override;
  virtual double top_of_atmosphere() override;
 private:
  double m_ttoa;
  double m_rho0;
  double m_zs;
  double m_zmax;
  double m_nmo0;
};  

struct LayeredAtmosphereLevel
{
  double z;
  double rho;
  double t;
  double nmo;

#ifndef SWIG
  class CmpZAsc { public: bool operator() (const LayeredAtmosphereLevel& a,
                                           const LayeredAtmosphereLevel& b) {
    return a.z < b.z; } };

  class CmpTDec { public: bool operator() (const LayeredAtmosphereLevel& a,
                                           const LayeredAtmosphereLevel& b) {
    return b.t < a.t; } };
#endif // ifndef SWIG
};

#ifndef SWIG
struct LayeredAtmosphereLayer
{
  LayeredAtmosphereLayer(double _v = 0): 
      zb(_v), zt(_v), rho0(_v), rhozs(_v), t0(_v), tzs(_v), tb(_v), tt(_v),
      nmo0(_v), nmozs(_v), ptc0(_v) { }
  double zb;
  double zt;
  double rho0;
  double rhozs;
  double t0;
  double tzs;
  double tb;
  double tt;
  double nmo0;
  double nmozs;
  double ptc0;
  bool operator< (const LayeredAtmosphereLayer& o) const { return zt<o.zt; }
  
  class CmpTDec { public: bool operator() (const LayeredAtmosphereLayer& a,
                                           const LayeredAtmosphereLayer& b) {
    return b.tt<a.tt; } };
};
#endif // ifndef SWIG

// Parameterized, layered model
class LayeredAtmosphere: public Atmosphere
{
 public:
  CALIN_TYPEALIAS(Level, LayeredAtmosphereLevel);
  
  LayeredAtmosphere(const std::string& filename);
  LayeredAtmosphere(const std::vector<Level> levels);

  virtual ~LayeredAtmosphere();
  double rho(double z) override;
  double thickness(double z) override;
  double n_minus_one(double z) override;
  double propagation_time_correction(double z) override;
  double z_for_thickness(double t) override;
  double top_of_atmosphere() override;

  const std::vector<Level>& getLevels() const { return m_levels; }

 private:
  CALIN_TYPEALIAS(Layer, LayeredAtmosphereLayer);

  void initialize();
  inline std::vector<Layer>::const_iterator findZ(double z) const;

  double m_ztoa;
  double m_ttoa;
  std::vector<Level> m_levels;
  std::vector<Layer> m_layers;
  mutable std::vector<Layer>::const_iterator m_ilayer;
};

} } } // namespace calin::air_shower::atmosphere
