#include "exercise5.hpp"
#include "../Common/Paths.hpp"
#include "../Transport/Tracking.hpp"
#include "../XS/CrossSection.hpp"
#include <omp.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <cmath>

namespace fs = std::filesystem;

static void ensure_dir(const fs::path &p)
{
  fs::create_directories(p);
}

static Material make_water(const CrossSectionLibrary &lib)
{
  const Nuclide *H = lib.get("H-1");
  const Nuclide *O = lib.get("O-16");
  Material m;
  m.name = "water";
  m.mix.species.push_back({H, 6.68723e-2});
  m.mix.species.push_back({O, 3.34362e-2});
  return m;
}

static Material make_air(const CrossSectionLibrary &lib)
{
  Material m;
  m.name = "air";
  const Nuclide *N = lib.get("N-14");
  if (!N)
    throw std::runtime_error("Missing N-14 cross sections");
  m.mix.species.push_back({N, 5.16062e-5});
  return m;
}

static Material make_steel(const CrossSectionLibrary &lib)
{
  Material m;
  m.name = "steel";
  const Nuclide *Fe = lib.get("Fe-56");
  if (!Fe)
    throw std::runtime_error("Missing Fe-56 cross sections");
  m.mix.species.push_back({Fe, 8.39767e-2});
  return m;
}

static GeometryContext build_ex4_cylinder_geo(double water_height_cm, const MaterialLibrary &matlib)
{
  double r_in = 20.0;
  double r_out = 20.5;
  double z_bot = 0.0;
  double z_top = 70.0;

  auto inner = std::make_shared<CylinderZ>("inner", Vec3{0, 0, 0}, r_in);
  auto outer = std::make_shared<CylinderZ>("outer", Vec3{0, 0, 0}, r_out);
  auto zmin = std::make_shared<PlaneZ>("zmin", z_bot);
  auto zsurf = std::make_shared<PlaneZ>("zsurf", water_height_cm);
  auto zmax = std::make_shared<PlaneZ>("zmax", z_top);

  auto in_fluid = surf_node(inner);
  auto z_water = and_node(not_node(surf_node(zmin)), surf_node(zsurf));
  auto z_air = and_node(not_node(surf_node(zsurf)), surf_node(zmax));

  Cell water_cell{"water", and_node(in_fluid, z_water)};
  Cell air_cell{"air", and_node(in_fluid, z_air)};
  auto z_between = and_node(not_node(surf_node(zmin)), surf_node(zmax));
  Cell steel_cell{"steel", and_node(and_node(not_node(surf_node(inner)), surf_node(outer)), z_between)};
  Cell outside{"outside", not_node(surf_node(outer))};

  Universe u;
  u.cells = {water_cell, air_cell, steel_cell, outside};

  GeometryContext geo;
  geo.universe = u;
  geo.bounds = BoundingBox{{-r_out - 1.0, -r_out - 1.0, z_bot - 1.0}, {r_out + 1.0, r_out + 1.0, z_top + 1.0}};
  geo.surfaces = {inner, outer, zmin, zsurf, zmax};
  geo.cell_material["water"] = matlib.get("water");
  geo.cell_material["air"] = matlib.get("air");
  geo.cell_material["steel"] = matlib.get("steel");
  return geo;
}

static std::vector<int> default_thread_counts()
{
  const int procs = std::max(1, omp_get_num_procs());
  const int max_threads = std::max(1, omp_get_max_threads());
  const int limit = std::max(1, std::min(procs, max_threads));

  std::vector<int> threads;
  for (int t = 1; t < limit; t *= 2)
    threads.push_back(t);
  if (threads.empty() || threads.back() != limit)
    threads.push_back(limit);
  return threads;
}

struct BenchStats
{
  double mean_s = 0.0;
  double std_s = 0.0;
};

static BenchStats bench_external_source(const GeometryContext &geo, const MaterialLibrary &matlib, bool surface_tracking, int histories, int repeats)
{
  repeats = std::max(1, repeats);
  histories = std::max(1, histories);

  double sum = 0.0;
  double sum2 = 0.0;
  volatile double sink = 0.0;

  
  {
    auto t = simulate_external_source(geo, matlib, surface_tracking, std::min(5000, histories), nullptr);
    sink += t.leaks;
  }

  for (int r = 0; r < repeats; ++r)
  {
    const auto t0 = std::chrono::steady_clock::now();
    auto t = simulate_external_source(geo, matlib, surface_tracking, histories, nullptr);
    const auto t1 = std::chrono::steady_clock::now();
    const double dt = std::chrono::duration<double>(t1 - t0).count();

    
    sink += t.leaks;

    sum += dt;
    sum2 += dt * dt;
  }
  (void)sink;

  const double mean = sum / repeats;
  const double var = std::max(0.0, (sum2 / repeats) - mean * mean);

  BenchStats s;
  s.mean_s = mean;
  s.std_s = std::sqrt(var);
  return s;
}

static void run_ex4_parallel_scaling()
{
  const fs::path particleDir = Paths::particleDataDir();
  const fs::path resultsDir = Paths::resultsDir("ex5");

  CrossSectionLibrary xs;
  xs.load_file((particleDir / "H1.dat").string());
  xs.load_file((particleDir / "O16.dat").string());
  xs.load_file((particleDir / "N14.dat").string());
  xs.load_file((particleDir / "Fe56.dat").string());

  MaterialLibrary matlib;
  matlib.add(make_water(xs));
  matlib.add(make_air(xs));
  matlib.add(make_steel(xs));

  const double water_height_cm = 35.0;
  GeometryContext geo = build_ex4_cylinder_geo(water_height_cm, matlib);

  const int histories = 50000;
  const int repeats = 3;

  const int procs = std::max(1, omp_get_num_procs());
  const int max_threads = std::max(1, omp_get_max_threads());
  const std::vector<int> threads = default_thread_counts();

  ensure_dir(resultsDir);

  auto run_one = [&](bool surface_tracking, const fs::path &outFile)
  {
    omp_set_dynamic(0);

    std::vector<BenchStats> stats;
    stats.reserve(threads.size());
    for (int t : threads)
    {
      omp_set_num_threads(t);
      stats.push_back(bench_external_source(geo, matlib, surface_tracking, histories, repeats));
    }

    const double t1 = stats.empty() ? 0.0 : stats.front().mean_s;
    std::ofstream out(outFile);
    out << "threads,num_procs,max_threads,histories,repeats,time_s_mean,time_s_std,speedup,efficiency\n";
    for (size_t i = 0; i < threads.size(); ++i)
    {
      const int t = threads[i];
      const double time_s = stats[i].mean_s;
      const double sp = (time_s > 0.0 && t1 > 0.0) ? (t1 / time_s) : 0.0;
      const double eff = (t > 0) ? (sp / t) : 0.0;
      out << t << "," << procs << "," << max_threads << ","
          << histories << "," << repeats << ","
          << time_s << "," << stats[i].std_s << ","
          << sp << "," << eff << "\n";
    }
  };

  run_one(true, resultsDir / "ex5_ex4_parallel_scaling_surface.csv");
  run_one(false, resultsDir / "ex5_ex4_parallel_scaling_delta.csv");
}

void runExercise5()
{
  std::cout << "Running Exercise 5 bonuses (parallel scaling)...\n";
  try
  {
    run_ex4_parallel_scaling();
    std::cout << "Scaling results written to results/ex5/ex5_ex4_parallel_scaling_surface.csv and results/ex5/ex5_ex4_parallel_scaling_delta.csv.\n";
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 5 error: " << e.what() << "\n";
  }
}
