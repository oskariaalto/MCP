#include "exercise4.hpp"
#include "../Transport/Tracking.hpp"
#include "../Transport/Materials.hpp"
#include "../XS/CrossSection.hpp"
#include "../General/helperFunctions.hpp"
#include "../Common/Paths.hpp"
#include "../Geometry/VolumeEstimator.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <chrono>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <stdexcept>
#include <cstdlib>
#include <omp.h>

namespace fs = std::filesystem;

static constexpr double kPi = 3.141592653589793238462643383279502884;

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

static Material make_enriched_solution(const CrossSectionLibrary &lib, double wt235, double mol_frac_u)
{
  
  
  
  
  
  mol_frac_u = std::clamp(mol_frac_u, 0.0, 1.0);

  const Nuclide *U5 = lib.get("U-235");
  const Nuclide *U8 = lib.get("U-238");
  const Nuclide *H = lib.get("H-1");
  const Nuclide *O = lib.get("O-16");
  if (!U5 || !U8 || !H || !O)
    throw std::runtime_error("Missing nuclide cross sections for enriched solution");

  Material m;
  m.name = "solution";

  
  const double N_U_total = 4.83e-2;    
  const double N_H_water = 6.68723e-2; 
  const double N_O_water = 3.34362e-2; 

  const double fU = mol_frac_u;
  const double fW = 1.0 - mol_frac_u;

  const double N_U = N_U_total * fU;
  const double N_H = N_H_water * fW;
  const double N_O = N_O_water * fW;

  m.mix.species.push_back({U5, N_U * wt235});
  m.mix.species.push_back({U8, N_U * (1.0 - wt235)});
  m.mix.species.push_back({H, N_H});
  m.mix.species.push_back({O, N_O});
  return m;
}

static GeometryContext build_cylinder_geo(double water_height, const MaterialLibrary &matlib)
{
  
  double r_in = 20.0;
  double r_out = 20.5;
  double z_bot = 0.0;
  double z_top = 70.0;

  auto inner = std::make_shared<CylinderZ>("inner", Vec3{0, 0, 0}, r_in);
  auto outer = std::make_shared<CylinderZ>("outer", Vec3{0, 0, 0}, r_out);
  auto zmin = std::make_shared<PlaneZ>("zmin", z_bot);
  auto zsurf = std::make_shared<PlaneZ>("zsurf", water_height);
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

static std::unordered_map<std::string, double> cylinder_region_volumes_cm3(double water_height_cm)
{
  double r_in = 20.0;
  double r_out = 20.5;
  double z_top = 70.0;

  std::unordered_map<std::string, double> v;
  v["water"] = kPi * r_in * r_in * water_height_cm;
  v["air"] = kPi * r_in * r_in * std::max(0.0, z_top - water_height_cm);
  v["steel"] = kPi * (r_out * r_out - r_in * r_in) * z_top;
  return v;
}

static void write_tallies(const TrackTallies &t, const std::unordered_map<std::string, double> &vol_cm3, int histories, double source_rate_per_s, const fs::path &file)
{
  ensure_dir(file.parent_path());
  std::ofstream out(file);
  out << "material,volume_cm3,"
         "collisions_analog_per_source,absorptions_analog_per_source,"
         "collisions_tle_per_source,absorptions_tle_per_source,"
         "flux_tle_per_source,flux_cfe_per_source,"
         "collisions_analog_per_s,absorptions_analog_per_s,"
         "collisions_tle_per_s,absorptions_tle_per_s,"
         "flux_tle_per_s,flux_cfe_per_s\n";

  double inv_hist = (histories > 0) ? (1.0 / histories) : 0.0;

  for (const auto &kv : vol_cm3)
  {
    const std::string &mat = kv.first;
    double V = kv.second;
    TrackTallies::Entry e{};
    auto it = t.by_material.find(mat);
    if (it != t.by_material.end())
      e = it->second;

    double c_a = e.collisions_analog * inv_hist;
    double a_a = e.absorptions_analog * inv_hist;
    double c_tle = e.collisions_tle * inv_hist;
    double a_tle = e.absorptions_tle * inv_hist;
    double tl = e.track_length * inv_hist;
    double phi_tle = (V > 0.0) ? (tl / V) : 0.0;
    double phi_cfe = (V > 0.0) ? ((e.flux_cfe * inv_hist) / V) : 0.0;

    out << mat << "," << V << ","
        << c_a << "," << a_a << ","
        << c_tle << "," << a_tle << ","
        << phi_tle << "," << phi_cfe << ","
        << c_a * source_rate_per_s << "," << a_a * source_rate_per_s << ","
        << c_tle * source_rate_per_s << "," << a_tle * source_rate_per_s << ","
        << phi_tle * source_rate_per_s << "," << phi_cfe * source_rate_per_s << "\n";
  }

  double leak_per_source = t.leaks * inv_hist;
  out << "leaks,,"
      << leak_per_source << ",,,,,"
      << leak_per_source * source_rate_per_s << ",,,,,\n";
}




static void write_mesh(const MeshTally &m, const BoundingBox &bounds, int histories, double source_rate_per_s, const fs::path &file)
{
  ensure_dir(file.parent_path());
  std::ofstream out(file);

  
  const double dx = (bounds.max.x - bounds.min.x) / m.nx;
  const double dy = (bounds.max.y - bounds.min.y) / m.ny;
  const double dz = (bounds.max.z - bounds.min.z) / m.nz;
  const double voxel_vol = dx * dy * dz;

  const double inv_hist = (histories > 0) ? (1.0 / histories) : 0.0;

  out << "ix,iy,iz,x_mid_cm,y_mid_cm,z_mid_cm,collisions,coll_rate_per_source_per_cm3,coll_rate_per_s_per_cm3\n";
  for (int iz = 0; iz < m.nz; ++iz)
    for (int iy = 0; iy < m.ny; ++iy)
      for (int ix = 0; ix < m.nx; ++ix)
      {
        size_t idx = (iz * m.ny + iy) * m.nx + ix;
        const double xmid = bounds.min.x + (ix + 0.5) * dx;
        const double ymid = bounds.min.y + (iy + 0.5) * dy;
        const double zmid = bounds.min.z + (iz + 0.5) * dz;

        const double c = m.collisions[idx];
        const double rate_per_source = (voxel_vol > 0.0) ? ((c * inv_hist) / voxel_vol) : 0.0;
        const double rate_per_s = rate_per_source * source_rate_per_s;

        out << ix << "," << iy << "," << iz << ","
            << xmid << "," << ymid << "," << zmid << ","
            << c << "," << rate_per_source << "," << rate_per_s << "\n";
      }
}

static void external_demo(const MaterialLibrary &matlib, bool surface_tracking)
{
  const fs::path resultsDir = Paths::resultsDir("ex4");
  const double water_height_cm = 35.0;
  const double source_rate_per_s = 1e6;

  const int omp_num_procs = std::max(1, omp_get_num_procs());
  const int omp_max_threads = std::max(1, omp_get_max_threads());
  const int omp_dynamic = omp_get_dynamic();
  int omp_env_threads = 0;
  if (const char *env = std::getenv("OMP_NUM_THREADS"))
    omp_env_threads = std::max(0, std::atoi(env));

  GeometryContext geo = build_cylinder_geo(water_height_cm, matlib);
  const auto volumes = cylinder_region_volumes_cm3(water_height_cm);

  const std::string tag = surface_tracking ? "surface" : "delta";

  
  const int batches = 20;
  const int histories_per_batch = 5000;
  const int total_histories = batches * histories_per_batch;

  struct RunningStats
  {
    int n = 0;
    double mean = 0.0;
    double m2 = 0.0;
    void add(double x)
    {
      ++n;
      double d = x - mean;
      mean += d / n;
      double d2 = x - mean;
      m2 += d * d2;
    }
    double sample_std() const { return (n > 1) ? std::sqrt(m2 / (n - 1)) : 0.0; }
    double stderr() const { return (n > 0) ? sample_std() / std::sqrt((double)n) : 0.0; }
    double rel_err() const
    {
      if (n == 0)
        return 0.0;
      const double denom = std::max(1e-30, std::abs(mean));
      return stderr() / denom;
    }
  };

  
  std::unordered_map<std::string, RunningStats> stats;
  auto key = [](const std::string &q, const std::string &m)
  { return q + "|" + m; };

  auto add_stats = [&](const std::string &q, const std::string &m, double x)
  { stats[key(q, m)].add(x); };

  
  TrackTallies total_tallies;
  auto add_tallies = [&](TrackTallies &dst, const TrackTallies &src)
  {
    dst.leaks += src.leaks;
    for (const auto &kv : src.by_material)
    {
      auto &d = dst.by_material[kv.first];
      d.collisions_analog += kv.second.collisions_analog;
      d.absorptions_analog += kv.second.absorptions_analog;
      d.track_length += kv.second.track_length;
      d.collisions_tle += kv.second.collisions_tle;
      d.absorptions_tle += kv.second.absorptions_tle;
      d.flux_cfe += kv.second.flux_cfe;
    }
  };

  
  MeshTally mesh;
  mesh.init(10, 10, 10, geo.bounds);

  const auto t0_all = std::chrono::high_resolution_clock::now();

  for (int b = 0; b < batches; ++b)
  {
    const auto t0 = std::chrono::high_resolution_clock::now();
    auto t = simulate_external_source(geo, matlib, surface_tracking, histories_per_batch, &mesh);
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double dt = std::chrono::duration<double>(t1 - t0).count();

    add_tallies(total_tallies, t);

    const double inv_hist = 1.0 / histories_per_batch;

    
    for (const auto &kv : volumes)
    {
      const std::string &mat = kv.first;
      const double V = kv.second;

      using Entry = std::decay_t<decltype(t.by_material.begin()->second)>;
      Entry e{};
      auto it = t.by_material.find(mat);
      if (it != t.by_material.end())
        e = it->second;

      const double c_a = e.collisions_analog * inv_hist;
      const double a_a = e.absorptions_analog * inv_hist;
      const double c_tle = e.collisions_tle * inv_hist;
      const double a_tle = e.absorptions_tle * inv_hist;
      const double phi_tle = (V > 0.0) ? ((e.track_length * inv_hist) / V) : 0.0;
      const double phi_cfe = (V > 0.0) ? ((e.flux_cfe * inv_hist) / V) : 0.0;

      add_stats("collisions_analog_per_source", mat, c_a);
      add_stats("absorptions_analog_per_source", mat, a_a);
      add_stats("collisions_tle_per_source", mat, c_tle);
      add_stats("absorptions_tle_per_source", mat, a_tle);
      add_stats("flux_tle_per_source", mat, phi_tle);
      add_stats("flux_cfe_per_source", mat, phi_cfe);
    }

    add_stats("leaks_per_source", "all", t.leaks * inv_hist);
    add_stats("batch_time_s", "all", dt);
  }

  const auto t1_all = std::chrono::high_resolution_clock::now();
  const double time_total_s = std::chrono::duration<double>(t1_all - t0_all).count();

  
  write_tallies(total_tallies, volumes, total_histories, source_rate_per_s, resultsDir / ("ex4_external_" + tag + ".csv"));
  write_mesh(mesh, geo.bounds, total_histories, source_rate_per_s, resultsDir / ("ex4_mesh_" + tag + ".csv"));

	  
	  {
	    std::ofstream out(resultsDir / ("ex4_external_" + tag + "_stats.csv"));
	    out << "quantity,material,mean_per_source,stderr_per_source,rel_error,mean_per_s,stderr_per_s,"
	           "time_total_s,time_per_history_s,omp_num_procs,omp_max_threads,omp_dynamic,omp_env_threads,fom\n";
	    const double time_per_history_s = (total_histories > 0) ? (time_total_s / total_histories) : 0.0;

	    auto emit = [&](const std::string &q, const std::string &mat)
	    {
	      const auto it = stats.find(key(q, mat));
	      if (it == stats.end())
        return;
      const RunningStats &s = it->second;

      double mean_ps = s.mean;
      double se_ps = s.stderr();
      double rel = s.rel_err();
      double mean_per_s = mean_ps * source_rate_per_s;
      double se_per_s = se_ps * source_rate_per_s;

      
      double fom = 0.0;
      if (rel > 0.0)
        fom = 1.0 / (rel * rel * std::max(1e-12, time_total_s));

	      out << q << "," << mat << ","
	          << mean_ps << "," << se_ps << "," << rel << ","
	          << mean_per_s << "," << se_per_s << ","
	          << time_total_s << "," << time_per_history_s << ","
	          << omp_num_procs << "," << omp_max_threads << "," << omp_dynamic << "," << omp_env_threads << ","
	          << fom << "\n";
	    };

    for (const auto &kv : volumes)
    {
      const std::string &mat = kv.first;
      emit("collisions_analog_per_source", mat);
      emit("absorptions_analog_per_source", mat);
      emit("collisions_tle_per_source", mat);
      emit("absorptions_tle_per_source", mat);
      emit("flux_tle_per_source", mat);
      emit("flux_cfe_per_source", mat);
    }
    emit("leaks_per_source", "all");

    
	    {
	      const auto it = stats.find(key("batch_time_s", "all"));
	      if (it != stats.end())
	        out << "avg_batch_time_s,all," << it->second.mean << "," << it->second.stderr() << "," << it->second.rel_err()
	            << ",,,"
	            << time_total_s << "," << time_per_history_s << ","
	            << omp_num_procs << "," << omp_max_threads << "," << omp_dynamic << "," << omp_env_threads << ",\n";
	    }
	  }
	}

static double critical_case_keff(const MaterialLibrary &base_matlib,
                                 const CrossSectionLibrary &xs,
                                 double water_height_cm,
                                 double wt235,
                                 double mol_frac_u,
                                 const std::string &tag,
                                 int histories = 3000,
                                 int generations = 40,
                                 int inactive = 10)
{
  const fs::path resultsDir = Paths::resultsDir("ex4");
  MaterialLibrary ml = base_matlib;

  Material sol = make_enriched_solution(xs, wt235, mol_frac_u);
  ml.add(sol);

  GeometryContext geo = build_cylinder_geo(water_height_cm, ml);
  
  geo.cell_material["water"] = ml.get("solution");

  auto res = simulate_criticality(geo, ml, histories, generations);

  ensure_dir(resultsDir);
  std::ofstream out(resultsDir / ("ex4_critical_" + tag + ".csv"));
  out << "gen,keff,water_height_cm,wt235,mol_frac_u,histories,generations,inactive\n";
  for (size_t i = 0; i < res.gen_keff.size(); ++i)
    out << i << "," << res.gen_keff[i] << ","
        << water_height_cm << "," << wt235 << "," << mol_frac_u << ","
        << histories << "," << generations << "," << inactive << "\n";

  
  const int g0 = std::min<int>(inactive, (int)res.gen_keff.size());
  double sum = 0.0;
  int n = 0;
  for (size_t i = g0; i < res.gen_keff.size(); ++i)
  {
    sum += res.gen_keff[i];
    ++n;
  }
  const double keff_active = (n > 0) ? (sum / n) : res.keff;

  
  std::ofstream summary(resultsDir / ("ex4_critical_" + tag + "_summary.csv"));
  summary << "water_height_cm,wt235,mol_frac_u,histories,generations,inactive,keff_active\n";
  summary << water_height_cm << "," << wt235 << "," << mol_frac_u << ","
          << histories << "," << generations << "," << inactive << ","
          << keff_active << "\n";
  return keff_active;
}

static double estimate_critical_height_50_50(const MaterialLibrary &matlib, const CrossSectionLibrary &xs)
{
  
  const fs::path resultsDir = Paths::resultsDir("ex4");
  ensure_dir(resultsDir);

  std::ofstream scan(resultsDir / "ex4_critical_scan_50_50.csv");
  scan << "water_height_cm,keff,wt235,mol_frac_u,histories,generations,inactive,tag\n";

  struct ScanRow
  {
    double water_height_cm = 0.0;
    double keff = 0.0;
    double wt235 = 0.0;
    double mol_frac_u = 0.0;
    int histories = 0;
    int generations = 0;
    int inactive = 0;
    std::string tag;
  };
  std::vector<ScanRow> rows;
  rows.reserve(32);

  auto keff_at = [&](double h, const std::string &tag)
  {
    const double wt235 = 0.20;
    const double mol_frac_u = 0.50;
    const int histories = 2000;
    const int generations = 30;
    const int inactive = 5;
    const double k = critical_case_keff(matlib, xs, h, wt235, mol_frac_u, tag, histories, generations, inactive);
    scan << h << "," << k << "," << wt235 << "," << mol_frac_u << ","
         << histories << "," << generations << "," << inactive << "," << tag << "\n";
    scan.flush();

    rows.push_back({h, k, wt235, mol_frac_u, histories, generations, inactive, tag});
    return k;
  };

  
  double best_h = 35.0;
  double best_k = 0.0;
  double best_diff = 1e9;

  double h_lo = 10.0, h_hi = 60.0;
  double k_lo = keff_at(h_lo, "50_50_h" + std::to_string((int)h_lo));
  double k_hi = keff_at(h_hi, "50_50_h" + std::to_string((int)h_hi));

  auto consider = [&](double h, double k)
  {
    const double d = std::abs(k - 1.0);
    if (d < best_diff)
    {
      best_diff = d;
      best_h = h;
      best_k = k;
    }
  };
  consider(h_lo, k_lo);
  consider(h_hi, k_hi);

  
  for (double h : {20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0})
  {
    double k = keff_at(h, "50_50_h" + std::to_string((int)h));
    consider(h, k);
  }

  
  if ((k_lo - 1.0) * (k_hi - 1.0) < 0.0)
  {
    for (int it = 0; it < 6; ++it)
    {
      const double h_mid = 0.5 * (h_lo + h_hi);
      const double k_mid = keff_at(h_mid, "50_50_mid" + std::to_string(it));

      consider(h_mid, k_mid);

      if ((k_lo - 1.0) * (k_mid - 1.0) < 0.0)
      {
        h_hi = h_mid;
        k_hi = k_mid;
      }
      else
      {
        h_lo = h_mid;
        k_lo = k_mid;
      }
    }
  }

  std::ofstream best(resultsDir / "ex4_critical_best_50_50.csv");
  best << "crit_height_cm,keff,wt235,mol_frac_u\n";
  best << best_h << "," << best_k << "," << 0.20 << "," << 0.50 << "\n";

  
  scan.close();
  std::sort(rows.begin(), rows.end(), [](const ScanRow &a, const ScanRow &b)
            { return a.water_height_cm < b.water_height_cm; });
  {
    std::ofstream sorted(resultsDir / "ex4_critical_scan_50_50.csv");
    sorted << "water_height_cm,keff,wt235,mol_frac_u,histories,generations,inactive,tag\n";
    for (const auto &r : rows)
    {
      sorted << r.water_height_cm << "," << r.keff << "," << r.wt235 << "," << r.mol_frac_u << ","
             << r.histories << "," << r.generations << "," << r.inactive << "," << r.tag << "\n";
    }
  }

  return best_h;
}

void runExercise4()
{
  std::cout << "Running Exercise 4 transport simulation...\n";
  fs::path particleDir;
  fs::path resultsDir;
  try
  {
    particleDir = Paths::particleDataDir();
    resultsDir = Paths::resultsDir("ex4");
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 4 error: " << e.what() << "\n";
    return;
  }

  CrossSectionLibrary xs;
  bool ok = true;
  ok &= xs.load_file((particleDir / "H1.dat").string());
  ok &= xs.load_file((particleDir / "O16.dat").string());
  ok &= xs.load_file((particleDir / "U235.dat").string());
  ok &= xs.load_file((particleDir / "U238.dat").string());
  ok &= xs.load_file((particleDir / "N14.dat").string());
  ok &= xs.load_file((particleDir / "Fe56.dat").string());
  if (!ok)
  {
    std::cout << "Failed to load one or more cross section files from " << particleDir.string() << ".\n";
    return;
  }

  MaterialLibrary mats;
  mats.add(make_water(xs));
  mats.add(make_air(xs));
  mats.add(make_steel(xs));

  
  external_demo(mats, true);
  external_demo(mats, false);

  
  
  
  const double crit_height = estimate_critical_height_50_50(mats, xs);

  const double k50 = critical_case_keff(mats, xs, crit_height, 0.20, 0.50, "50_50_at_critical", 5000, 50, 10);
  const double k25 = critical_case_keff(mats, xs, crit_height, 0.20, 0.25, "25_75_at_critical", 5000, 50, 10);
  const double k75 = critical_case_keff(mats, xs, crit_height, 0.20, 0.75, "75_25_at_critical", 5000, 50, 10);

  ensure_dir(resultsDir);
  std::ofstream kout(resultsDir / "ex4_keff_at_critical_height.csv");
  kout << "height_cm,wt235,mol_frac_u_50_50,keff_50_50,mol_frac_u_25_75,keff_25_75,mol_frac_u_75_25,keff_75_25\n";
  kout << crit_height << "," << 0.20 << ","
       << 0.50 << "," << k50 << ","
       << 0.25 << "," << k25 << ","
       << 0.75 << "," << k75 << "\n";

  
  
  
  const double power_W = 1.0;
  const double fission_energy_J = 200.0 * 1.602e-13; 
  const double fission_rate = power_W / fission_energy_J;
  const double nubar = 2.4;

  
  
  const double gen_time_s = 1e-4;

  const double neutron_rate = fission_rate * nubar;
  const double r_in = 20.0;
  const double volume_cm3 = kPi * r_in * r_in * crit_height;

  const double neutron_density_cm3 = (volume_cm3 > 0.0) ? ((neutron_rate * gen_time_s) / volume_cm3) : 0.0;

  
  const Material sol = make_enriched_solution(xs, 0.20, 0.50);
  double atomic_density_1e24 = 0.0;
  for (const auto &sp : sol.mix.species)
    atomic_density_1e24 += sp.number_density;
  const double atomic_density_cm3 = atomic_density_1e24 * 1e24;

  const double ratio = (atomic_density_cm3 > 0.0) ? (neutron_density_cm3 / atomic_density_cm3) : 0.0;

  std::ofstream nd(resultsDir / "ex4_neutron_density.csv");
  nd << "height_cm,keff_50_50,power_W,gen_time_s,neutron_density_cm3,atomic_density_cm3,ratio\n";
  nd << crit_height << "," << k50 << "," << power_W << "," << gen_time_s << ","
     << neutron_density_cm3 << "," << atomic_density_cm3 << "," << ratio << "\n";

  std::cout << "Exercise 4 outputs written to results/ex4/.\n";
}
