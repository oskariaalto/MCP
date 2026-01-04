#ifndef CROSS_SECTION_HPP
#define CROSS_SECTION_HPP
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include "../Common/Interpolation.hpp"

struct Table1D
{
  std::vector<double> energy;
  std::vector<double> value;

  double eval(double e) const
  {
    return log_interp(energy, value, e);
  }
};

struct Reaction
{
  int mt = 0;
  double q_value = 0.0; 
  Table1D table;
};

struct Nuclide
{
  std::string name;
  int Z = 0;
  int A = 0;
  double atomic_weight = 0.0;
  double temperature = 0.0;
  Table1D nubar;                              
  std::unordered_map<int, Reaction> reactions; 

  double sigma_mt(int mt, double E) const
  {
    auto it = reactions.find(mt);
    if (it == reactions.end())
      return 0.0;
    return it->second.table.eval(E);
  }

  
  
  
  double sigma_elastic(double E) const
  {
    return sigma_mt(2, E);
  }

  double sigma_inelastic(double E) const
  {
    double sum = 0.0;
    for (const auto &kv : reactions)
    {
      int mt = kv.first;
      if (mt >= 50 && mt < 90)
        sum += kv.second.table.eval(E);
    }
    return sum;
  }

  double sigma_total(double E) const
  {
    double sum = 0.0;
    for (const auto &kv : reactions)
    {
      sum += kv.second.table.eval(E);
    }
    return sum;
  }

  double sigma_capture(double E) const
  {
    return sigma_mt(102, E);
  }

  double sigma_fission(double E) const
  {
    return sigma_mt(18, E);
  }

  double sigma_absorption(double E) const
  {
    
    double sum = sigma_fission(E);
    for (const auto &kv : reactions)
    {
      int mt = kv.first;
      if (mt >= 101)
        sum += kv.second.table.eval(E);
    }
    return sum;
  }

  double nubar_value(double E) const
  {
    if (nubar.energy.empty())
      return 0.0;
    return nubar.eval(E);
  }
};

struct Species
{
  const Nuclide *nuclide = nullptr;
  double number_density = 0.0; 
};

struct Mixture
{
  std::vector<Species> species;

  double macro_sigma_total(double E) const;
  double macro_sigma_reaction(int mt, double E) const;
};

class CrossSectionLibrary
{
public:
  bool load_file(const std::string &path);
  const Nuclide *get(const std::string &name) const;
  const std::unordered_map<std::string, Nuclide> &all() const { return nuclides_; }

private:
  Nuclide parse_file(const std::string &path);
  std::unordered_map<std::string, Nuclide> nuclides_;
};

#endif
