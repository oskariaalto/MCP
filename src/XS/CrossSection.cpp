#include "CrossSection.hpp"
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <iostream>

static Table1D read_table(std::ifstream &in, int n)
{
  Table1D t;
  t.energy.reserve(n);
  t.value.reserve(n);
  for (int i = 0; i < n; ++i)
  {
    double e, v;
    if (!(in >> e >> v))
      throw std::runtime_error("Unexpected end of file while reading table");
    t.energy.push_back(e);
    t.value.push_back(v);
  }
  return t;
}

Nuclide CrossSectionLibrary::parse_file(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Could not open xs file " + path);

  Nuclide n;
  if (!(in >> n.name >> n.Z >> n.A >> n.atomic_weight >> n.temperature))
    throw std::runtime_error("Bad header in " + path);

  int nubar_points = 0;
  in >> nubar_points;
  if (nubar_points > 0)
  {
    n.nubar = read_table(in, nubar_points);
  }

  while (true)
  {
    int mt, num;
    double q;
    if (!(in >> mt >> q >> num))
      break;
    Table1D t = read_table(in, num);
    Reaction r;
    r.mt = mt;
    r.q_value = q;
    r.table = std::move(t);
    n.reactions[mt] = std::move(r);
  }

  return n;
}

bool CrossSectionLibrary::load_file(const std::string &path)
{
  try
  {
    Nuclide n = parse_file(path);
    nuclides_[n.name] = std::move(n);
    return true;
  }
  catch (const std::exception &e)
  {
    std::cerr << "Cross section load failed: " << e.what() << "\n";
    return false;
  }
}

const Nuclide *CrossSectionLibrary::get(const std::string &name) const
{
  auto it = nuclides_.find(name);
  if (it == nuclides_.end())
    return nullptr;
  return &it->second;
}

double Mixture::macro_sigma_total(double E) const
{
  double sum = 0.0;
  for (const auto &s : species)
  {
    if (s.nuclide)
      sum += s.number_density * s.nuclide->sigma_total(E); 
  }
  return sum; 
}

double Mixture::macro_sigma_reaction(int mt, double E) const
{
  double sum = 0.0;
  for (const auto &s : species)
  {
    if (s.nuclide)
      sum += s.number_density * s.nuclide->sigma_mt(mt, E);
  }
  return sum;
}
