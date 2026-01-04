#ifndef NEUTRON_HPP
#define NEUTRON_HPP
#include "../XS/CrossSection.hpp"
#include "../Common/Random.hpp"
#include "../Geometry/Vec3.hpp"
#include <algorithm>
#include <vector>
#include <cmath>

struct Neutron
{
  double energy; 
};

struct ScatterResult
{
  double energy = 0.0; 
  Vec3 dir{};          
};

inline Vec3 sample_isotropic_direction()
{
  double mu = 2.0 * Random::instance().uniform() - 1.0;
  double phi = 2.0 * M_PI * Random::instance().uniform();
  double s = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return Vec3{s * std::cos(phi), s * std::sin(phi), mu};
}

inline Vec3 cross(const Vec3 &a, const Vec3 &b)
{
  return Vec3{
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x};
}

inline Vec3 rotate_from_axis(const Vec3 &axis_dir, double mu, double phi)
{
  Vec3 w = axis_dir.normalized();
  Vec3 helper = (std::abs(w.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
  Vec3 u = cross(helper, w).normalized();
  Vec3 v = cross(w, u);

  double s = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return u * (s * std::cos(phi)) + v * (s * std::sin(phi)) + w * mu;
}

inline ScatterResult scatter_elastic_stationary(double E_in, const Vec3 &dir_in, double A)
{
  
  double mu_cm = 2.0 * Random::instance().uniform() - 1.0;
  double phi = 2.0 * M_PI * Random::instance().uniform();

  double factor = (A * A + 1.0 + 2.0 * A * mu_cm) / ((A + 1.0) * (A + 1.0));
  factor = std::max(0.0, factor);
  double E_out = E_in * factor;

  
  double denom = std::sqrt(std::max(1e-30, A * A + 1.0 + 2.0 * A * mu_cm));
  double mu_lab = (1.0 + A * mu_cm) / denom;
  mu_lab = std::clamp(mu_lab, -1.0, 1.0);

  Vec3 out_dir = rotate_from_axis(dir_in, mu_lab, phi).normalized();
  return {E_out, out_dir};
}

inline double scatter_energy(double E_in, double A)
{
  
  return scatter_elastic_stationary(E_in, Vec3{0.0, 0.0, 1.0}, A).energy;
}

inline ScatterResult scatter_elastic_free_gas(double E_in, const Vec3 &dir_in, double target_atomic_weight_amu, double temperature_K)
{
  if (E_in <= 0.0)
    return {0.0, dir_in.normalized()};
  constexpr double neutron_mass_amu = 1.00866491588;
  double A_ratio = (target_atomic_weight_amu > 0.0) ? (target_atomic_weight_amu / neutron_mass_amu) : 0.0;
  if (temperature_K <= 0.0 || A_ratio <= 0.0)
    return scatter_elastic_stationary(E_in, dir_in, A_ratio);

  
  
  constexpr double kB_J_per_K = 1.380649e-23;
  constexpr double amu_kg = 1.66053906660e-27;
  constexpr double neutron_mass_kg = 1.67492749804e-27;
  constexpr double MeV_to_J = 1.602176634e-13;

  double target_mass_kg = target_atomic_weight_amu * amu_kg;
  double sigma_v = std::sqrt(kB_J_per_K * temperature_K / std::max(1e-30, target_mass_kg)); 

  Vec3 v_target{
      Random::instance().sample_normal(0.0, sigma_v),
      Random::instance().sample_normal(0.0, sigma_v),
      Random::instance().sample_normal(0.0, sigma_v)};

  double v_n = std::sqrt(2.0 * E_in * MeV_to_J / neutron_mass_kg); 
  Vec3 v_neutron = dir_in.normalized() * v_n;

  double m = neutron_mass_kg;
  double M = target_mass_kg;
  Vec3 V_cm = (v_neutron * m + v_target * M) / (m + M);
  Vec3 v_cm = v_neutron - V_cm;
  double v_cm_mag = v_cm.norm();
  if (v_cm_mag <= 0.0)
    return {E_in, dir_in.normalized()};

  
  double mu = 2.0 * Random::instance().uniform() - 1.0;
  double phi = 2.0 * M_PI * Random::instance().uniform();
  Vec3 v_cm_dir_out = rotate_from_axis(v_cm, mu, phi).normalized();
  Vec3 v_cm_out = v_cm_dir_out * v_cm_mag;

  Vec3 v_out = v_cm_out + V_cm;
  double E_out = 0.5 * neutron_mass_kg * v_out.dot(v_out) / MeV_to_J;
  return {std::max(0.0, E_out), v_out.normalized()};
}

inline int sample_fission_multiplicity(double nubar)
{
  if (nubar <= 0.0)
    return 0;
  return std::max(0, Random::instance().sample_poisson(nubar));
}

inline double sample_fission_energy()
{
  
  return Random::instance().sample_maxwell(1.29);
}


struct MultiplicationResult
{
  double average_total = 0.0;
  double stddev = 0.0;
};

MultiplicationResult simulate_multiplication(const Mixture &mix, double source_energy, int histories, bool include_inelastic);

struct SlowingDownResult
{
  std::vector<int> collision_index;
  std::vector<double> average_energy;
};

SlowingDownResult simulate_slowing_down(const Mixture &mix, double source_energy, int histories, int max_collisions);

struct ReactionTally
{
  std::string nuclide;
  std::string reaction;
  double mean = 0.0;
  double rel_err = 0.0;
};

std::vector<ReactionTally> tally_reactions(const Mixture &mix, double source_energy, int histories, bool include_inelastic);

struct TimeSeries
{
  std::vector<double> bin_edges;
  std::vector<double> counts;
};

TimeSeries simulate_time_dependent(const Mixture &mix, double source_energy, int histories, double t_max, int bins, bool include_inelastic);

#endif
