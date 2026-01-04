#include "exercise2.hpp"
#include "../XS/CrossSection.hpp"
#include "../Physics/Neutron.hpp"
#include "../General/helperFunctions.hpp"
#include "../Common/Paths.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <stdexcept>

namespace fs = std::filesystem;

static void ensure_dir(const fs::path &p)
{
  fs::create_directories(p);
}

static const Species *sample_species(const Mixture &mix, double E, double &sigma_total_macro, double &sigma_total_micro)
{
  sigma_total_macro = mix.macro_sigma_total(E);
  if (sigma_total_macro <= 0.0)
    return nullptr;
  double xi = Random::instance().uniform() * sigma_total_macro;
  double accum = 0.0;
  for (const auto &s : mix.species)
  {
    if (!s.nuclide)
      continue;
    double contrib = s.number_density * s.nuclide->sigma_total(E);
    accum += contrib;
    if (xi <= accum)
    {
      sigma_total_micro = s.nuclide->sigma_total(E);
      return &s;
    }
  }
  
  for (auto it = mix.species.rbegin(); it != mix.species.rend(); ++it)
  {
    if (it->nuclide)
    {
      sigma_total_micro = it->nuclide->sigma_total(E);
      return &(*it);
    }
  }
  return nullptr;
}

static std::vector<double> logspace(double a, double b, int n)
{
  std::vector<double> xs;
  xs.reserve(n);
  double loga = std::log10(a);
  double logb = std::log10(b);
  double step = (logb - loga) / (n - 1);
  for (int i = 0; i < n; ++i)
  {
    xs.push_back(std::pow(10.0, loga + step * i));
  }
  return xs;
}

static void write_csv(const fs::path &file, const std::vector<std::string> &headers, const std::vector<std::vector<double>> &cols)
{
  std::ofstream out(file);
  for (size_t i = 0; i < headers.size(); ++i)
  {
    out << headers[i];
    if (i + 1 != headers.size())
      out << ",";
  }
  out << "\n";
  if (cols.empty())
    return;
  size_t rows = cols[0].size();
  for (size_t r = 0; r < rows; ++r)
  {
    for (size_t c = 0; c < cols.size(); ++c)
    {
      out << cols[c][r];
      if (c + 1 != cols.size())
        out << ",";
    }
    out << "\n";
  }
}

static Mixture make_water(const CrossSectionLibrary &lib)
{
  const Nuclide *H = lib.get("H-1");
  const Nuclide *O = lib.get("O-16");
  if (!H || !O)
    throw std::runtime_error("Missing H-1 or O-16 cross sections");
  Mixture mix;
  mix.species.push_back({H, 6.68723e-2}); 
  mix.species.push_back({O, 3.34362e-2});
  return mix;
}

static Mixture make_natural_uranium(const CrossSectionLibrary &lib)
{
  const Nuclide *U235 = lib.get("U-235");
  const Nuclide *U238 = lib.get("U-238");
  if (!U235 || !U238)
    throw std::runtime_error("Missing uranium cross sections");
  
  double total = 4.83e-2;
  Mixture mix;
  mix.species.push_back({U235, total * 0.0072});
  mix.species.push_back({U238, total * 0.9928});
  return mix;
}

static Mixture make_UO2_natural(const CrossSectionLibrary &lib)
{
  const Nuclide *U235 = lib.get("U-235");
  const Nuclide *U238 = lib.get("U-238");
  const Nuclide *O = lib.get("O-16");
  if (!U235 || !U238 || !O)
    throw std::runtime_error("Missing U/O cross sections");
  
  double u_total = 2.45e-2;
  double o_total = u_total * 2.0;
  Mixture mix;
  mix.species.push_back({U235, u_total * 0.0072});
  mix.species.push_back({U238, u_total * 0.9928});
  mix.species.push_back({O, o_total});
  return mix;
}

static Mixture make_UO2_pureU238(const CrossSectionLibrary &lib)
{
  const Nuclide *U238 = lib.get("U-238");
  const Nuclide *O = lib.get("O-16");
  if (!U238 || !O)
    throw std::runtime_error("Missing U-238/O-16 cross sections");
  double u_total = 2.45e-2;
  double o_total = u_total * 2.0;
  Mixture mix;
  mix.species.push_back({U238, u_total});
  mix.species.push_back({O, o_total});
  return mix;
}

[[maybe_unused]] static Mixture make_water_uranium_mix(const CrossSectionLibrary &lib)
{
  Mixture water = make_water(lib);
  Mixture uran = make_natural_uranium(lib);
  Mixture mix;
  for (auto s : water.species)
    mix.species.push_back(s);
  for (auto s : uran.species)
    mix.species.push_back(s);
  return mix;
}

static Nuclide make_synthetic_deuterium(const CrossSectionLibrary &lib)
{
  const Nuclide *H = lib.get("H-1");
  if (!H)
    throw std::runtime_error("Missing H-1 for synthetic D-2");
  Nuclide D = *H;
  D.name = "D-2";
  D.A = 2;
  D.atomic_weight = 2.01410177812;
  return D;
}

static void export_cross_sections(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  auto energy = logspace(1e-11, 20.0, 500);

  const Nuclide *H = lib.get("H-1");
  const Nuclide *O = lib.get("O-16");
  std::vector<double> h_tot, o_tot;
  for (double e : energy)
  {
    h_tot.push_back(H ? H->sigma_total(e) : 0.0);
    o_tot.push_back(O ? O->sigma_total(e) : 0.0);
  }
  write_csv(resultsDir / "ex2_micro_total_H_O.csv", {"energy_MeV", "H1_sigmaT", "O16_sigmaT"}, {energy, h_tot, o_tot});

  const Nuclide *U5 = lib.get("U-235");
  const Nuclide *U8 = lib.get("U-238");
  std::vector<double> u5_f, u5_c, u8_f, u8_c, u8_inl;
  for (double e : energy)
  {
    u5_f.push_back(U5 ? U5->sigma_fission(e) : 0.0);
    u5_c.push_back(U5 ? U5->sigma_capture(e) : 0.0);
    u8_f.push_back(U8 ? U8->sigma_fission(e) : 0.0);
    u8_c.push_back(U8 ? U8->sigma_capture(e) : 0.0);
    
    
    if (U8)
    {
      double inl = U8->sigma_mt(4, e);
      if (inl <= 0.0)
      {
        
        for (int mt = 50; mt < 90; ++mt)
          inl += U8->sigma_mt(mt, e);
      }
      u8_inl.push_back(inl);
    }
    else
    {
      u8_inl.push_back(0.0);
    }
  }
  write_csv(resultsDir / "ex2_micro_U_fission_capture.csv",
            {"energy_MeV", "U235_fission", "U235_capture", "U238_fission", "U238_capture"},
            {energy, u5_f, u5_c, u8_f, u8_c});
  write_csv(resultsDir / "ex2_micro_U238_inelastic.csv", {"energy_MeV", "U238_inelastic"}, {energy, u8_inl});

  Mixture water = make_water(lib);
  Mixture natU = make_natural_uranium(lib);
  std::vector<double> water_tot, nat_tot;
  for (double e : energy)
  {
    water_tot.push_back(water.macro_sigma_total(e));
    nat_tot.push_back(natU.macro_sigma_total(e));
  }
  write_csv(resultsDir / "ex2_macro_total.csv", {"energy_MeV", "water_SigT_cm^-1", "natU_SigT_cm^-1"}, {energy, water_tot, nat_tot});
}

static void run_slowdown(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  
  const Nuclide *H = lib.get("H-1");
  if (!H)
    throw std::runtime_error("Missing H-1 cross sections");
  Mixture hydrogen;
  hydrogen.species.push_back({H, 6.68723e-2});

  auto res = simulate_slowing_down(hydrogen, 2.0, 500, 40);
  std::vector<double> idx(res.collision_index.begin(), res.collision_index.end());
  write_csv(resultsDir / "ex2_slowing_hydrogen.csv", {"collision_index", "avg_energy_MeV"}, {idx, res.average_energy});

  
  const Nuclide *D_ptr = lib.get("H-2");
  Nuclide D_synth;
  if (!D_ptr)
  {
    D_synth = make_synthetic_deuterium(lib);
    D_ptr = &D_synth;
  }
  Mixture deuterium;
  deuterium.species.push_back({D_ptr, 6.68723e-2});
  auto resD = simulate_slowing_down(deuterium, 2.0, 500, 40);
  std::vector<double> idxD(resD.collision_index.begin(), resD.collision_index.end());
  write_csv(resultsDir / "ex2_slowing_deuterium.csv", {"collision_index", "avg_energy_MeV"}, {idxD, resD.average_energy});
}

static void run_multiplication(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  
  Mixture uo2_u238 = make_UO2_pureU238(lib);
  Mixture uo2_nat = make_UO2_natural(lib);
  Mixture natU = make_natural_uranium(lib);
  Mixture water = make_water(lib);

  
  for (auto &s : natU.species)
    s.number_density *= 0.5;
  for (auto &s : water.species)
    s.number_density *= 0.5;
  Mixture mix_50_50;
  for (auto s : natU.species)
    mix_50_50.species.push_back(s);
  for (auto s : water.species)
    mix_50_50.species.push_back(s);

  std::vector<std::string> labels = {"UO2_pureU238", "UO2_naturalU", "naturalU_water_50_50"};
  std::vector<Mixture *> mixes = {&uo2_u238, &uo2_nat, &mix_50_50};
  std::vector<double> avg, stddev, avg_inel;
  for (auto m : mixes)
  {
    auto res = simulate_multiplication(*m, 1.0, 10000, false);
    auto res_inel = simulate_multiplication(*m, 1.0, 10000, true);
    avg.push_back(res.average_total);
    stddev.push_back(res.stddev);
    avg_inel.push_back(res_inel.average_total);
  }
  std::vector<double> idx(labels.size());
  for (size_t i = 0; i < labels.size(); ++i)
    idx[i] = static_cast<double>(i);
  write_csv(resultsDir / "ex2_multiplication.csv", {"case_index", "avg_fission_neutrons_elastic", "stddev", "avg_fission_neutrons_with_inelastic"}, {idx, avg, stddev, avg_inel});
}

static void run_reaction_stats(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  
  Mixture uo2_u238 = make_UO2_pureU238(lib);
  Mixture uo2_nat = make_UO2_natural(lib);

  Mixture natU = make_natural_uranium(lib);
  Mixture water = make_water(lib);
  for (auto &s : natU.species)
    s.number_density *= 0.5;
  for (auto &s : water.species)
    s.number_density *= 0.5;
  Mixture mix_50_50;
  for (auto s : natU.species)
    mix_50_50.species.push_back(s);
  for (auto s : water.species)
    mix_50_50.species.push_back(s);

  int histories = 10000;
  const double source_E = 1.0;

  auto tallies_u238_no = tally_reactions(uo2_u238, source_E, histories, false);
  auto tallies_u238_in = tally_reactions(uo2_u238, source_E, histories, true);
  auto tallies_nat_no = tally_reactions(uo2_nat, source_E, histories, false);
  auto tallies_nat_in = tally_reactions(uo2_nat, source_E, histories, true);
  auto tallies_mix_no = tally_reactions(mix_50_50, source_E, histories, false);
  auto tallies_mix_in = tally_reactions(mix_50_50, source_E, histories, true);

  auto write = [&](const std::vector<ReactionTally> &t, const std::string &name)
  {
    std::ofstream out(resultsDir / name);
    out << "index,nuclide,reaction,mean,rel_err\n";
    for (size_t i = 0; i < t.size(); ++i)
    {
      out << i << "," << t[i].nuclide << "," << t[i].reaction << "," << t[i].mean << "," << t[i].rel_err << "\n";
    }
  };

  write(tallies_u238_no, "ex2_reaction_stats_UO2_pureU238_no_inelastic.csv");
  write(tallies_u238_in, "ex2_reaction_stats_UO2_pureU238_with_inelastic.csv");
  write(tallies_nat_no, "ex2_reaction_stats_UO2_naturalU_no_inelastic.csv");
  write(tallies_nat_in, "ex2_reaction_stats_UO2_naturalU_with_inelastic.csv");
  write(tallies_mix_no, "ex2_reaction_stats_natU_water_50_50_no_inelastic.csv");
  write(tallies_mix_in, "ex2_reaction_stats_natU_water_50_50_with_inelastic.csv");
}

static Mixture make_enriched_uranium(const CrossSectionLibrary &lib, double wt235)
{
  const Nuclide *U235 = lib.get("U-235");
  const Nuclide *U238 = lib.get("U-238");
  if (!U235 || !U238)
    throw std::runtime_error("Missing uranium cross sections");
  double total = 4.83e-2;
  Mixture mix;
  mix.species.push_back({U235, total * wt235});
  mix.species.push_back({U238, total * (1.0 - wt235)});
  return mix;
}

static Mixture make_water_enriched_mix(const CrossSectionLibrary &lib, double wt235)
{
  Mixture water = make_water(lib);
  for (auto &s : water.species)
    s.number_density *= 0.5; 
  Mixture u = make_enriched_uranium(lib, wt235);
  for (auto &s : u.species)
    s.number_density *= 0.5;
  Mixture m;
  for (auto s : water.species)
    m.species.push_back(s);
  for (auto s : u.species)
    m.species.push_back(s);
  return m;
}

static void run_critical_scan(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  std::vector<double> enrich = {0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 0.7, 0.9};
  std::vector<double> avgU, avgMix;
  for (double e : enrich)
  {
    auto resU = simulate_multiplication(make_enriched_uranium(lib, e), 1.0, 300, true);
    auto resMix = simulate_multiplication(make_water_enriched_mix(lib, e), 1.0, 300, true);
    avgU.push_back(resU.average_total);
    avgMix.push_back(resMix.average_total);
  }
  write_csv(resultsDir / "ex2_critical_scan_uranium.csv", {"wt235", "avg_fission_neutrons"}, {enrich, avgU});
  write_csv(resultsDir / "ex2_critical_scan_water_mix.csv", {"wt235", "avg_fission_neutrons"}, {enrich, avgMix});
}

static void run_time_dependence(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  Mixture natU = make_natural_uranium(lib);
  Mixture enriched50 = make_enriched_uranium(lib, 0.5);

  
  
  int histories = 2000;
  double t_max = 1e-6;
  int bins = 400;
  double source_E = 5.0;

  auto write = [&](const Mixture &m, const std::string &fname)
  {
    TimeSeries ts = simulate_time_dependent(m, source_E, histories, t_max, bins, true);
    std::ofstream out(resultsDir / fname);
    out << "t_mid_s,counts\n";
    for (size_t i = 0; i < ts.counts.size(); ++i)
    {
      double mid = 0.5 * (ts.bin_edges[i] + ts.bin_edges[i + 1]);
      out << mid << "," << ts.counts[i] << "\n";
    }
  };

  write(natU, "ex2_time_natU.csv");
  write(enriched50, "ex2_time_enriched50.csv");
}

static void run_four_factor(const CrossSectionLibrary &lib)
{
  const fs::path resultsDir = Paths::resultsDir("ex2");
  ensure_dir(resultsDir);
  
  Mixture mix = make_water_enriched_mix(lib, 0.0072); 
  int histories = 100;
  double total_fission_neutrons = 0.0;
  double absorptions_fuel = 0.0;
  double absorptions_total = 0.0;

  for (int h = 0; h < histories; ++h)
  {
    std::vector<Neutron> stack;
    stack.push_back({1.0});
    int steps = 0;
    const int max_steps = 500;
    while (!stack.empty())
    {
      if (++steps > max_steps)
        break;
      Neutron n = stack.back();
      stack.pop_back();
      double sigma_macro = 0.0;
      double sigma_micro = 0.0;
      const Species *spec = sample_species(mix, n.energy, sigma_macro, sigma_micro);
      if (!spec || sigma_macro <= 0.0)
        continue;

      const Nuclide *nuclide = spec->nuclide;
      double sf = nuclide->sigma_fission(n.energy);
      double sa = nuclide->sigma_capture(n.energy);
      
      
      double sinel_total = nuclide->sigma_mt(4, n.energy);
      if (sinel_total <= 0.0)
      {
        for (int mt = 50; mt < 90; ++mt)
          sinel_total += nuclide->sigma_mt(mt, n.energy);
      }
      double s_total = std::max(0.0, sigma_micro - sinel_total);
      double xi = Random::instance().uniform() * s_total;
      bool is_fuel = nuclide->Z == 92;
      double sabs_other = std::max(0.0, nuclide->sigma_absorption(n.energy) - sa - sf);
      if (xi < sa + sabs_other)
      {
        absorptions_total += 1.0;
        if (is_fuel)
          absorptions_fuel += 1.0;
        continue;
      }
      else if (xi < sa + sabs_other + sf)
      {
        absorptions_total += 1.0;
        if (is_fuel)
        {
          absorptions_fuel += 1.0;
        }
        double nubar = nuclide->nubar_value(n.energy);
        int k = sample_fission_multiplicity(nubar);
        total_fission_neutrons += k;
        int max_children = 10;
        for (int i = 0; i < k && i < max_children; ++i)
          stack.push_back({sample_fission_energy()});
      }
      else
      {
        
        double s_elastic = std::max(0.0, s_total - (sa + sabs_other + sf));
        if (s_elastic <= 0.0)
          continue;

        double e_scatter = scatter_energy(n.energy, nuclide->A);
        if (e_scatter > 0.0)
          stack.push_back({e_scatter});
      }
    }
  }

  double eta = (absorptions_fuel > 0) ? total_fission_neutrons / absorptions_fuel : 0.0;
  double f = (absorptions_total > 0) ? absorptions_fuel / absorptions_total : 0.0;
  double p = 1.0;       
  double epsilon = 1.0; 
  double keff = eta * f * p * epsilon;

  std::ofstream out(resultsDir / "ex2_four_factor.csv");
  out << "eta,f,p,epsilon,keff\n";
  out << eta << "," << f << "," << p << "," << epsilon << "," << keff << "\n";
}

void runExercise2()
{
  std::cout << "Running Exercise 2 routines...\n";
  fs::path particleDir;
  try
  {
    particleDir = Paths::particleDataDir();
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 2 error: " << e.what() << "\n";
    return;
  }

  CrossSectionLibrary lib;
  bool ok = true;
  ok &= lib.load_file((particleDir / "H1.dat").string());
  ok &= lib.load_file((particleDir / "O16.dat").string());
  ok &= lib.load_file((particleDir / "U235.dat").string());
  ok &= lib.load_file((particleDir / "U238.dat").string());
  if (fs::exists(particleDir / "H2.dat"))
    ok &= lib.load_file((particleDir / "H2.dat").string());
  if (!ok)
  {
    std::cout << "Failed to load one or more cross section files.\n";
    return;
  }

  try
  {
    std::cout << "[ex2] exporting cross sections..." << std::flush;
    export_cross_sections(lib);
    std::cout << "done\n";

    std::cout << "[ex2] slowing-down curves..." << std::flush;
    run_slowdown(lib);
    std::cout << "done\n";

    std::cout << "[ex2] multiplication studies..." << std::flush;
    run_multiplication(lib);
    std::cout << "done\n";

    std::cout << "[ex2] reaction stats..." << std::flush;
    run_reaction_stats(lib);
    std::cout << "done\n";

    std::cout << "[ex2] critical enrichment scan..." << std::flush;
    run_critical_scan(lib);
    std::cout << "done\n";

    std::cout << "[ex2] time-dependent tally..." << std::flush;
    run_time_dependence(lib);
    std::cout << "done\n";

    std::cout << "[ex2] four-factor approximation..." << std::flush;
    run_four_factor(lib);
    std::cout << "done\n";
    std::cout << "Exercise 2 results written to results/ex2/.\n";
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 2 error: " << e.what() << "\n";
  }
}
