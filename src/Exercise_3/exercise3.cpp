#include "exercise3.hpp"
#include "../Geometry/CSG.hpp"
#include "../Geometry/VolumeEstimator.hpp"
#include "../Geometry/Surface.hpp"
#include "../Common/Random.hpp"
#include "../Common/Paths.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>

namespace fs = std::filesystem;

static void ensure_dir(const fs::path &p)
{
  fs::create_directories(p);
}

struct DemoResult
{
  std::string name;

  
  double vol_points = 0.0;
  double vol_lines = 0.0;

  
  double vol_points_stderr = 0.0;
  double vol_lines_stderr = 0.0;
  double rel_points = 0.0;
  double rel_lines = 0.0;

  
  double time_points_s = 0.0;
  double time_lines_s = 0.0;
  double fom_points = 0.0;
  double fom_lines = 0.0;

  double analytic = 0.0;
};

struct VolumeSimConfig
{
  int batches = 20;
  int total_point_samples = 1000000; 
  int total_lines = 20000;           
  int segs_per_line = 1000;          
};

static const VolumeSimConfig kVolumeSimConfig{};

static Universe build_fuel_pin()
{
  
  
  
  
  const double r_fuel = 0.824 / 2.0;
  const double clad_thickness = 0.063;
  const double r_clad = r_fuel + clad_thickness;
  const double pitch = 1.330;
  const double half_pitch = pitch / 2.0;
  const double z0 = 0.0;
  const double z1 = 10.0;

  auto fuel = std::make_shared<CylinderZ>("fuel_cyl", Vec3{0, 0, 0}, r_fuel);
  auto clad = std::make_shared<CylinderZ>("clad_cyl", Vec3{0, 0, 0}, r_clad);
  auto zmin = std::make_shared<PlaneZ>("zmin", z0);
  auto zmax = std::make_shared<PlaneZ>("zmax", z1);
  auto pin_box = std::make_shared<Cuboid>("pin_box", Vec3{-half_pitch, -half_pitch, z0}, Vec3{half_pitch, half_pitch, z1});

  auto z_between = and_node(not_node(surf_node(zmin)), surf_node(zmax)); 
  Cell fuel_cell{"fuel", and_node(surf_node(fuel), z_between)};
  Cell clad_cell{"clad", and_node(and_node(not_node(surf_node(fuel)), surf_node(clad)), z_between)};
  Cell coolant_cell{"coolant", and_node(and_node(not_node(surf_node(clad)), surf_node(pin_box)), z_between)};

  Universe u;
  u.cells = {fuel_cell, clad_cell, coolant_cell};
  return u;
}

static Universe build_hollow_cylinder()
{
  double r_inner = 2.0;
  double r_outer = 2.05;
  double z_bot = 0.0;
  double z_mid = 3.5;
  double z_top = 7.0;
  auto inner = std::make_shared<CylinderZ>("inner", Vec3{0, 0, 0}, r_inner);
  auto outer = std::make_shared<CylinderZ>("outer", Vec3{0, 0, 0}, r_outer);
  auto zmin = std::make_shared<PlaneZ>("zmin", z_bot);
  auto zmid = std::make_shared<PlaneZ>("zmid", z_mid);
  auto zmax = std::make_shared<PlaneZ>("zmax", z_top);
  auto bbox = std::make_shared<Cuboid>("hollow_bbox", Vec3{-2.2, -2.2, z_bot}, Vec3{2.2, 2.2, z_top});

  auto in_fluid = surf_node(inner);                                    
  auto z_water = and_node(not_node(surf_node(zmin)), surf_node(zmid)); 
  auto z_air = and_node(not_node(surf_node(zmid)), surf_node(zmax));   

  Cell water{"water", and_node(in_fluid, z_water)};
  Cell air{"air", and_node(in_fluid, z_air)};
  
  auto z_all = and_node(not_node(surf_node(zmin)), surf_node(zmax));
  Cell steel{"steel", and_node(and_node(not_node(surf_node(inner)), surf_node(outer)), z_all)};
  Cell outside{"outside", and_node(not_node(surf_node(outer)), surf_node(bbox))};
  Universe u;
  u.cells = {water, air, steel, outside};
  return u;
}

static Universe build_hex_prism_demo()
{
  auto hex = std::make_shared<HexPrism>("hex", 1.0, 0.0, 2.0);
  Cell hex_cell{"hex", surf_node(hex)};
  Cell outside{"outside", not_node(surf_node(hex))};
  Universe u;
  u.cells = {hex_cell, outside};
  return u;
}

static void ascii_plot(const Universe &u, const BoundingBox &b, int nx, int ny, double z, const std::string &filename)
{
  std::ofstream out(filename);
  for (int j = ny - 1; j >= 0; --j)
  {
    for (int i = 0; i < nx; ++i)
    {
      double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
      double y = b.min.y + (b.max.y - b.min.y) * (j + 0.5) / ny;
      Vec3 p{x, y, z};
      auto *cell = u.find_cell(p);
      char c = '.';
      if (cell)
      {
        if (cell->name == "fuel")
          c = 'F';
        else if (cell->name == "clad")
          c = 'C';
        else if (cell->name == "void")
          c = 'V';
        else if (cell->name == "water")
          c = 'W';
        else if (cell->name == "air")
          c = 'A';
        else if (cell->name == "steel")
          c = 'S';
        else if (cell->name == "coolant")
          c = 'M';
        else if (cell->name == "outside")
          c = 'o';
        else if (cell->name == "hex")
          c = 'H';
        else if (cell->name == "sphere_base")
          c = 'S';
        else if (cell->name == "sphere_tx")
          c = 'T';
        else if (cell->name == "box_base")
          c = 'B';
        else if (cell->name == "box_rot")
          c = 'R';
        else if (cell->name == "torus")
          c = 'O';
        else if (cell->name == "etor")
          c = 'E';
        else if (cell->name == "gplane")
          c = 'P';
      }
      out << c;
    }
    out << "\n";
  }
}

static void ascii_plot_xz(const Universe &u, const BoundingBox &b, int nx, int nz, double y, const std::string &filename)
{
  
  std::ofstream out(filename);
  for (int k = nz - 1; k >= 0; --k)
  {
    for (int i = 0; i < nx; ++i)
    {
      double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
      double z = b.min.z + (b.max.z - b.min.z) * (k + 0.5) / nz;
      Vec3 p{x, y, z};
      auto *cell = u.find_cell(p);
      char c = '.';
      if (cell)
      {
        if (cell->name == "fuel")
          c = 'F';
        else if (cell->name == "clad")
          c = 'C';
        else if (cell->name == "void")
          c = 'V';
        else if (cell->name == "water")
          c = 'W';
        else if (cell->name == "air")
          c = 'A';
        else if (cell->name == "steel")
          c = 'S';
        else if (cell->name == "coolant")
          c = 'M';
        else if (cell->name == "outside")
          c = 'o';
        else if (cell->name == "hex")
          c = 'H';
        else if (cell->name == "sphere_base")
          c = 'S';
        else if (cell->name == "sphere_tx")
          c = 'T';
        else if (cell->name == "box_base")
          c = 'B';
        else if (cell->name == "box_rot")
          c = 'R';
        else if (cell->name == "torus")
          c = 'O';
        else if (cell->name == "etor")
          c = 'E';
        else if (cell->name == "gplane")
          c = 'P';
      }
      out << c;
    }
    out << "\n";
  }
}

static void ascii_plot_compare(const Cell &original,
                               const Cell &transformed,
                               const BoundingBox &b,
                               int nx,
                               int ny,
                               double z,
                               const std::string &filename)
{
  
  
  std::ofstream out(filename);
  for (int j = ny - 1; j >= 0; --j)
  {
    for (int i = 0; i < nx; ++i)
    {
      double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
      double y = b.min.y + (b.max.y - b.min.y) * (j + 0.5) / ny;
      Vec3 p{x, y, z};
      const bool in_orig = eval_node(*original.root, p);
      const bool in_new = eval_node(*transformed.root, p);
      char c = '.';
      if (in_orig && in_new)
        c = 'X';
      else if (in_orig)
        c = 'O';
      else if (in_new)
        c = 'N';
      out << c;
    }
    out << "\n";
  }
}

static double sample_stddev(const std::vector<double> &x)
{
  if (x.size() < 2)
    return 0.0;
  const double mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
  double s2 = 0.0;
  for (double v : x)
    s2 += (v - mean) * (v - mean);
  s2 /= static_cast<double>(x.size() - 1); 
  return std::sqrt(std::max(s2, 0.0));
}

struct MethodStats
{
  double mean = 0.0;
  double stderr = 0.0; 
  double time_s = 0.0;
};

static MethodStats batched_points(const Cell &cell, const BoundingBox &bbox, int total_samples, int batches)
{
  MethodStats st;
  if (batches <= 0)
    return st;
  const int n_per = std::max(1, total_samples / batches);

  std::vector<double> batch_means;
  batch_means.reserve(batches);

  const auto t0 = std::chrono::high_resolution_clock::now();
  for (int b = 0; b < batches; ++b)
  {
    batch_means.push_back(estimate_volume_points(cell, bbox, n_per));
  }
  const auto t1 = std::chrono::high_resolution_clock::now();
  st.time_s = std::chrono::duration<double>(t1 - t0).count();

  st.mean = std::accumulate(batch_means.begin(), batch_means.end(), 0.0) / static_cast<double>(batch_means.size());
  const double s = sample_stddev(batch_means);
  st.stderr = s / std::sqrt(static_cast<double>(batch_means.size()));
  return st;
}

static MethodStats batched_lines(const Cell &cell, const BoundingBox &bbox, int total_lines, int segs_per_line, int batches)
{
  MethodStats st;
  if (batches <= 0)
    return st;
  const int n_per = std::max(1, total_lines / batches);

  std::vector<double> batch_means;
  batch_means.reserve(batches);

  const auto t0 = std::chrono::high_resolution_clock::now();
  for (int b = 0; b < batches; ++b)
  {
    batch_means.push_back(estimate_volume_lines(cell, bbox, n_per, segs_per_line));
  }
  const auto t1 = std::chrono::high_resolution_clock::now();
  st.time_s = std::chrono::duration<double>(t1 - t0).count();

  st.mean = std::accumulate(batch_means.begin(), batch_means.end(), 0.0) / static_cast<double>(batch_means.size());
  const double s = sample_stddev(batch_means);
  st.stderr = s / std::sqrt(static_cast<double>(batch_means.size()));
  return st;
}

static double safe_rel(double mean, double stderr)
{
  const double denom = std::abs(mean);
  if (denom <= 0.0)
    return 0.0;
  return stderr / denom;
}

static double safe_fom(double rel, double time_s)
{
  if (time_s <= 0.0)
    return 0.0;
  if (rel <= 0.0)
    return 0.0;
  return 1.0 / (rel * rel * time_s);
}

static void run_volume_case(const std::string &name, const Cell &cell, const BoundingBox &bbox, double analytic, std::vector<DemoResult> &results)
{
  
  
  MethodStats pts = batched_points(cell, bbox, kVolumeSimConfig.total_point_samples, kVolumeSimConfig.batches);
  MethodStats lns = batched_lines(cell, bbox, kVolumeSimConfig.total_lines, kVolumeSimConfig.segs_per_line, kVolumeSimConfig.batches);

  const double rel_p = safe_rel(pts.mean, pts.stderr);
  const double rel_l = safe_rel(lns.mean, lns.stderr);

  const double fom_p = safe_fom(rel_p, pts.time_s);
  const double fom_l = safe_fom(rel_l, lns.time_s);

  std::cout << name << " | points: " << pts.mean << " ± " << pts.stderr << " (rel " << rel_p
            << ") | lines: " << lns.mean << " ± " << lns.stderr << " (rel " << rel_l
            << ") | analytic: " << analytic << "\n";

  DemoResult r;
  r.name = name;
  r.vol_points = pts.mean;
  r.vol_lines = lns.mean;
  r.analytic = analytic;

  r.vol_points_stderr = pts.stderr;
  r.vol_lines_stderr = lns.stderr;
  r.rel_points = rel_p;
  r.rel_lines = rel_l;

  r.time_points_s = pts.time_s;
  r.time_lines_s = lns.time_s;
  r.fom_points = fom_p;
  r.fom_lines = fom_l;

  results.push_back(r);
}

static void write_results(const std::vector<DemoResult> &res)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  std::ofstream out(resultsDir / "ex3_volumes.csv");
  out << "name,"
      << "vol_points,vol_points_stderr,rel_points,time_points_s,fom_points,"
      << "vol_lines,vol_lines_stderr,rel_lines,time_lines_s,fom_lines,"
      << "analytic\n";

  for (const auto &r : res)
  {
    out << r.name << ","
        << r.vol_points << "," << r.vol_points_stderr << "," << r.rel_points << "," << r.time_points_s << "," << r.fom_points << ","
        << r.vol_lines << "," << r.vol_lines_stderr << "," << r.rel_lines << "," << r.time_lines_s << "," << r.fom_lines << ","
        << r.analytic << "\n";
  }
}

static void write_split_results(const std::vector<DemoResult> &res)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  
  {
    std::ofstream out(resultsDir / "ex3_volumes_points.csv");
    out << "name,vol,stderr,rel,time_s,fom,analytic\n";
    for (const auto &r : res)
      out << r.name << "," << r.vol_points << "," << r.vol_points_stderr << "," << r.rel_points << "," << r.time_points_s << "," << r.fom_points
          << "," << r.analytic << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_volumes_lines.csv");
    out << "name,vol,stderr,rel,time_s,fom,analytic\n";
    for (const auto &r : res)
      out << r.name << "," << r.vol_lines << "," << r.vol_lines_stderr << "," << r.rel_lines << "," << r.time_lines_s << "," << r.fom_lines
          << "," << r.analytic << "\n";
  }

  
  const fs::path perCaseDir = resultsDir / "ex3_volumes_by_case";
  ensure_dir(perCaseDir);
  for (const auto &r : res)
  {
    std::ofstream out(perCaseDir / ("ex3_volume_" + r.name + ".csv"));
    out << "method,vol,stderr,rel,time_s,fom,analytic\n";
    out << "points," << r.vol_points << "," << r.vol_points_stderr << "," << r.rel_points << "," << r.time_points_s << "," << r.fom_points << "," << r.analytic << "\n";
    out << "lines," << r.vol_lines << "," << r.vol_lines_stderr << "," << r.rel_lines << "," << r.time_lines_s << "," << r.fom_lines << "," << r.analytic << "\n";
  }
}

static void write_volume_estimator_settings()
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  std::ofstream out(resultsDir / "ex3_volume_estimator_settings.csv");
  out << "batches,total_point_samples,point_samples_per_batch,total_lines,lines_per_batch,segs_per_line\n";
  out << kVolumeSimConfig.batches << "," << kVolumeSimConfig.total_point_samples << "," << (kVolumeSimConfig.total_point_samples / kVolumeSimConfig.batches) << ","
      << kVolumeSimConfig.total_lines << "," << (kVolumeSimConfig.total_lines / kVolumeSimConfig.batches) << "," << kVolumeSimConfig.segs_per_line << "\n";
}

static void lattice_demo()
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  Universe pin = build_fuel_pin();

  
  
  Universe pure_coolant;
  {
    const double pitch = 1.330;
    const double half_pitch = pitch / 2.0;
    const double height = 10.0;
    auto box = std::make_shared<Cuboid>("pure_coolant_box", Vec3{-half_pitch, -half_pitch, 0.0}, Vec3{half_pitch, half_pitch, height});
    pure_coolant.cells = {Cell{"coolant", surf_node(box)}};
  }

  Lattice lat;
  
  const double pitch = 1.330;
  lat.pitch = {pitch, pitch, 0.0};
  lat.nx = 3;
  lat.ny = 3;
  lat.origin = {-(lat.nx * pitch) / 2.0, -(lat.ny * pitch) / 2.0, 0.0};
  lat.universe = &pin;

  const double height = 10.0;
  BoundingBox b{{lat.origin.x, lat.origin.y, 0.0}, {lat.origin.x + lat.nx * pitch, lat.origin.y + lat.ny * pitch, height}};

  
  const double r_fuel = 0.824 / 2.0;
  const double r_clad = r_fuel + 0.063;
  const double v_fuel_pin = M_PI * r_fuel * r_fuel * height;
  const double v_clad_pin = M_PI * (r_clad * r_clad - r_fuel * r_fuel) * height;
  const double v_coolant_pin = (pitch * pitch - M_PI * r_clad * r_clad) * height;

  auto estimate_lattice_volumes_points = [&](auto &&find_cell_fn, double &time_s)
  {
    const int batches = 20;
    const int total_samples = std::max(100000, kVolumeSimConfig.total_point_samples);
    const int n_per = std::max(1, total_samples / batches);

    std::vector<double> fuel_batch, clad_batch, coolant_batch;
    fuel_batch.reserve(batches);
    clad_batch.reserve(batches);
    coolant_batch.reserve(batches);

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int bidx = 0; bidx < batches; ++bidx)
    {
      int hits_fuel = 0;
      int hits_clad = 0;
      int hits_coolant = 0;
#pragma omp parallel
      {
        int lf = 0, lc = 0, lm = 0;
#pragma omp for nowait
        for (int i = 0; i < n_per; ++i)
        {
          Vec3 p = sample_point(b);
          auto *cell = find_cell_fn(p);
          if (!cell)
            continue;
          if (cell->name == "fuel")
            lf += 1;
          else if (cell->name == "clad")
            lc += 1;
          else if (cell->name == "coolant")
            lm += 1;
        }
#pragma omp atomic
        hits_fuel += lf;
#pragma omp atomic
        hits_clad += lc;
#pragma omp atomic
        hits_coolant += lm;
      }

      const double V = box_volume(b);
      fuel_batch.push_back((static_cast<double>(hits_fuel) / n_per) * V);
      clad_batch.push_back((static_cast<double>(hits_clad) / n_per) * V);
      coolant_batch.push_back((static_cast<double>(hits_coolant) / n_per) * V);
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    time_s = std::chrono::duration<double>(t1 - t0).count();

    auto mean = [](const std::vector<double> &x)
    {
      return std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
    };
    auto stderr = [&](const std::vector<double> &x)
    {
      const double s = sample_stddev(x);
      return s / std::sqrt(static_cast<double>(x.size()));
    };

    struct Out
    {
      double fuel = 0.0, clad = 0.0, coolant = 0.0;
      double fuel_se = 0.0, clad_se = 0.0, coolant_se = 0.0;
    };
    Out out;
    out.fuel = mean(fuel_batch);
    out.clad = mean(clad_batch);
    out.coolant = mean(coolant_batch);
    out.fuel_se = stderr(fuel_batch);
    out.clad_se = stderr(clad_batch);
    out.coolant_se = stderr(coolant_batch);
    return out;
  };

  
  double time_9 = 0.0;
  auto res9 = estimate_lattice_volumes_points([&](const Vec3 &p) { return lat.find_cell(p); }, time_9);
  const int n_pins_9 = lat.nx * lat.ny;
  const double fuel_analytic_9 = n_pins_9 * v_fuel_pin;
  const double clad_analytic_9 = n_pins_9 * v_clad_pin;
  const double coolant_analytic_9 = n_pins_9 * v_coolant_pin;

  std::cout << "Square lattice demo (" << lat.nx << "x" << lat.ny << ", pitch=" << pitch << " cm)\n";
  std::cout << "  fuel   ~ " << res9.fuel << " ± " << res9.fuel_se << " analytic " << fuel_analytic_9 << "\n";
  std::cout << "  clad   ~ " << res9.clad << " ± " << res9.clad_se << " analytic " << clad_analytic_9 << "\n";
  std::cout << "  coolant~ " << res9.coolant << " ± " << res9.coolant_se << " analytic " << coolant_analytic_9 << "\n";

  {
    std::ofstream out(resultsDir / "ex3_square_lattice_volumes.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "fuel," << res9.fuel << "," << res9.fuel_se << "," << safe_rel(res9.fuel, res9.fuel_se) << "," << time_9 << "," << fuel_analytic_9 << "\n";
    out << "clad," << res9.clad << "," << res9.clad_se << "," << safe_rel(res9.clad, res9.clad_se) << "," << time_9 << "," << clad_analytic_9 << "\n";
    out << "coolant," << res9.coolant << "," << res9.coolant_se << "," << safe_rel(res9.coolant, res9.coolant_se) << "," << time_9 << "," << coolant_analytic_9 << "\n";
  }
  
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_fuel.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "fuel," << res9.fuel << "," << res9.fuel_se << "," << safe_rel(res9.fuel, res9.fuel_se) << "," << time_9 << "," << fuel_analytic_9 << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_clad.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "clad," << res9.clad << "," << res9.clad_se << "," << safe_rel(res9.clad, res9.clad_se) << "," << time_9 << "," << clad_analytic_9 << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_coolant.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "coolant," << res9.coolant << "," << res9.coolant_se << "," << safe_rel(res9.coolant, res9.coolant_se) << "," << time_9 << "," << coolant_analytic_9 << "\n";
  }

  
  auto find_cell_8pins = [&](const Vec3 &p) -> const Cell *
  {
    if (lat.pitch.x == 0.0 || lat.pitch.y == 0.0)
      return nullptr;
    const double lx = (p.x - lat.origin.x) / lat.pitch.x;
    const double ly = (p.y - lat.origin.y) / lat.pitch.y;
    const int ix = static_cast<int>(std::floor(lx));
    const int iy = static_cast<int>(std::floor(ly));
    if (ix < 0 || iy < 0 || ix >= lat.nx || iy >= lat.ny)
      return nullptr;

    Vec3 local{
        p.x - (lat.origin.x + (ix + 0.5) * lat.pitch.x),
        p.y - (lat.origin.y + (iy + 0.5) * lat.pitch.y),
        p.z};

    const bool center = (ix == lat.nx / 2) && (iy == lat.ny / 2);
    return center ? pure_coolant.find_cell(local) : pin.find_cell(local);
  };

  double time_8 = 0.0;
  auto res8 = estimate_lattice_volumes_points(find_cell_8pins, time_8);

  const int n_pins_8 = lat.nx * lat.ny - 1;
  const double fuel_analytic_8 = n_pins_8 * v_fuel_pin;
  const double clad_analytic_8 = n_pins_8 * v_clad_pin;
  const double coolant_analytic_8 = n_pins_8 * v_coolant_pin + (pitch * pitch * height); 

  std::cout << "Square lattice (8 pins + center coolant)\n";
  std::cout << "  fuel   ~ " << res8.fuel << " ± " << res8.fuel_se << " analytic " << fuel_analytic_8 << "\n";
  std::cout << "  clad   ~ " << res8.clad << " ± " << res8.clad_se << " analytic " << clad_analytic_8 << "\n";
  std::cout << "  coolant~ " << res8.coolant << " ± " << res8.coolant_se << " analytic " << coolant_analytic_8 << "\n";

  {
    std::ofstream out(resultsDir / "ex3_square_lattice_8pins_volumes.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "fuel," << res8.fuel << "," << res8.fuel_se << "," << safe_rel(res8.fuel, res8.fuel_se) << "," << time_8 << "," << fuel_analytic_8 << "\n";
    out << "clad," << res8.clad << "," << res8.clad_se << "," << safe_rel(res8.clad, res8.clad_se) << "," << time_8 << "," << clad_analytic_8 << "\n";
    out << "coolant," << res8.coolant << "," << res8.coolant_se << "," << safe_rel(res8.coolant, res8.coolant_se) << "," << time_8 << "," << coolant_analytic_8 << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_8pins_fuel.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "fuel," << res8.fuel << "," << res8.fuel_se << "," << safe_rel(res8.fuel, res8.fuel_se) << "," << time_8 << "," << fuel_analytic_8 << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_8pins_clad.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "clad," << res8.clad << "," << res8.clad_se << "," << safe_rel(res8.clad, res8.clad_se) << "," << time_8 << "," << clad_analytic_8 << "\n";
  }
  {
    std::ofstream out(resultsDir / "ex3_square_lattice_8pins_coolant.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "coolant," << res8.coolant << "," << res8.coolant_se << "," << safe_rel(res8.coolant, res8.coolant_se) << "," << time_8 << "," << coolant_analytic_8 << "\n";
  }

  
  const int nx = 140;
  const int ny = 140;

  auto write_ascii = [&](const fs::path &file, auto &&find_cell_fn)
  {
    std::ofstream out(file);
    for (int j = ny - 1; j >= 0; --j)
    {
      for (int i = 0; i < nx; ++i)
      {
        double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
        double y = b.min.y + (b.max.y - b.min.y) * (j + 0.5) / ny;
        Vec3 p{x, y, 5.0};
        auto *cell = find_cell_fn(p);
        char c = '.';
        if (cell)
        {
          if (cell->name == "fuel")
            c = 'F';
          else if (cell->name == "clad")
            c = 'C';
          else
            c = 'M';
        }
        out << c;
      }
      out << "\n";
    }
  };

  
  write_ascii(resultsDir / "ex3_plot_square_lattice_z5.txt", [&](const Vec3 &p) { return lat.find_cell(p); });
  
  write_ascii(resultsDir / "ex3_plot_square_lattice_8pins_z5.txt", find_cell_8pins);

  
  Universe fuel_only;
  {
    const double r_fuel = 0.824 / 2.0;
    const double z0 = 0.0;
    const double z1 = height;
    auto fuel = std::make_shared<CylinderZ>("fuel_only_cyl", Vec3{0, 0, 0}, r_fuel);
    auto zmin = std::make_shared<PlaneZ>("fuel_only_zmin", z0);
    auto zmax = std::make_shared<PlaneZ>("fuel_only_zmax", z1);
    auto z_between = and_node(not_node(surf_node(zmin)), surf_node(zmax));
    fuel_only.cells = {Cell{"fuel", and_node(surf_node(fuel), z_between)}};
  }
  Lattice lat_fuel = lat;
  lat_fuel.universe = &fuel_only;

  
  {
    const int batches = 20;
    const int total_samples = std::max(100000, kVolumeSimConfig.total_point_samples);
    const int n_per = std::max(1, total_samples / batches);
    std::vector<double> batch;
    batch.reserve(batches);

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int bidx = 0; bidx < batches; ++bidx)
    {
      int hits = 0;
#pragma omp parallel for reduction(+ : hits)
      for (int i = 0; i < n_per; ++i)
      {
        Vec3 p = sample_point(b);
        auto *cell = lat_fuel.find_cell(p);
        if (cell && cell->name == "fuel")
          hits += 1;
      }
      batch.push_back((static_cast<double>(hits) / n_per) * box_volume(b));
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double time_s = std::chrono::duration<double>(t1 - t0).count();

    const double vol = std::accumulate(batch.begin(), batch.end(), 0.0) / static_cast<double>(batch.size());
    const double se = sample_stddev(batch) / std::sqrt(static_cast<double>(batch.size()));

    const double v_fuel_pin = M_PI * (0.824 / 2.0) * (0.824 / 2.0) * height;
    const double analytic = (lat.nx * lat.ny) * v_fuel_pin;

    std::ofstream out(resultsDir / "ex3_square_lattice_fuel_only.csv");
    out << "region,vol_points,stderr,rel,time_s,analytic\n";
    out << "fuel," << vol << "," << se << "," << safe_rel(vol, se) << "," << time_s << "," << analytic << "\n";

    std::cout << "Square lattice (3x3 fuel-only)\n";
    std::cout << "  fuel ~ " << vol << " ± " << se << " analytic " << analytic << "\n";
  }

  
  {
    std::ofstream out(resultsDir / "ex3_plot_square_lattice_fuel_only_z5.txt");
    for (int j = ny - 1; j >= 0; --j)
    {
      for (int i = 0; i < nx; ++i)
      {
        double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
        double y = b.min.y + (b.max.y - b.min.y) * (j + 0.5) / ny;
        Vec3 p{x, y, 5.0};
        auto *cell = lat_fuel.find_cell(p);
        out << ((cell && cell->name == "fuel") ? 'F' : '.');
      }
      out << "\n";
    }
  }
}

static void hex_lattice_demo()
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);
  Universe pin = build_fuel_pin();
  HexLattice hlat;
  hlat.pitch = 1.2;
  hlat.origin = {0.0, 0.0, 0.0};
  hlat.radius = 2; 
  hlat.universe = &pin;

  double extent = hlat.pitch * (hlat.radius * 2 + 1) * 1.8;
  BoundingBox b{{-extent, -extent, 0.0}, {extent, extent, 10.0}};

  int hits_fuel = 0;
  int samples = 50000;
#pragma omp parallel for reduction(+ : hits_fuel)
  for (int i = 0; i < samples; ++i)
  {
    Vec3 p = sample_point(b);
    auto *cell = hlat.find_cell(p);
    if (cell && cell->name == "fuel")
      hits_fuel += 1;
  }
  double frac = static_cast<double>(hits_fuel) / samples;
  double vol_est = frac * box_volume(b);
  std::cout << "Hex lattice demo fuel volume fraction ~ " << frac << " volume ~ " << vol_est << "\n";

  
  std::ofstream out(resultsDir / "ex3_plot_hex_lattice_z5.txt");
  int nx = 140, ny = 140;
  for (int j = ny - 1; j >= 0; --j)
  {
    for (int i = 0; i < nx; ++i)
    {
      double x = b.min.x + (b.max.x - b.min.x) * (i + 0.5) / nx;
      double y = b.min.y + (b.max.y - b.min.y) * (j + 0.5) / ny;
      Vec3 p{x, y, 5.0};
      auto *cell = hlat.find_cell(p);
      char c = '.';
      if (cell)
      {
        if (cell->name == "fuel")
          c = 'F';
        else if (cell->name == "clad")
          c = 'C';
        else
          c = 'V';
      }
      out << c;
    }
    out << "\n";
  }
}

static void transform_demo(std::vector<DemoResult> &results)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  

  
  auto sphere = std::make_shared<Sphere>("sphere_base", Vec3{0, 0, 0}, 1.0);
  auto moved_sphere = std::make_shared<TransformSurface>("sphere_tx", sphere, Vec3{1.0, 0.5, 0.0}, 0.0, 0.0, 0.0);
  Cell sphere_base_cell{"sphere_base", surf_node(sphere)};
  Cell sphere_cell{"sphere_tx", surf_node(moved_sphere)};
  Universe us_base;
  us_base.cells = {sphere_base_cell};
  Universe us;
  us.cells = {sphere_cell};
  BoundingBox bs{{-1.0, -1.0, -1.0}, {3.0, 2.0, 1.0}};
  double sphere_analytic = 4.0 / 3.0 * M_PI; 
  run_volume_case("translated_sphere", sphere_cell, bs, sphere_analytic, results);
  ascii_plot(us_base, bs, 140, 100, 0.0, (resultsDir / "ex3_plot_sphere_base_z0.txt").string());
  ascii_plot(us, bs, 140, 100, 0.0, (resultsDir / "ex3_plot_translated_sphere_z0.txt").string());
  ascii_plot_compare(sphere_base_cell, sphere_cell, bs, 140, 100, 0.0, (resultsDir / "ex3_plot_translated_sphere_compare_z0.txt").string());

  
  auto box = std::make_shared<Cuboid>("box_base", Vec3{-0.5, -0.5, -0.5}, Vec3{0.5, 0.5, 0.5});
  
  
  double angle = 45.0 * M_PI / 180.0;
  auto rotated_box = std::make_shared<TransformSurface>("box_rot", box, Vec3{0, 0, 0}, 0.0, 0.0, angle);
  Cell box_base_cell{"box_base", surf_node(box)};
  Cell box_cell{"box_rot", surf_node(rotated_box)};
  Universe ub_base;
  ub_base.cells = {box_base_cell};
  Universe ub;
  ub.cells = {box_cell};
  BoundingBox bb{{-1.5, -1.5, -1.5}, {1.5, 1.5, 1.5}};
  double box_analytic = 1.0;
  run_volume_case("rotated_cube", box_cell, bb, box_analytic, results);
  ascii_plot(ub_base, bb, 140, 140, 0.0, (resultsDir / "ex3_plot_cube_base_z0.txt").string());
  ascii_plot(ub, bb, 140, 140, 0.0, (resultsDir / "ex3_plot_rotated_cube_z0.txt").string());
  ascii_plot_compare(box_base_cell, box_cell, bb, 140, 140, 0.0, (resultsDir / "ex3_plot_rotated_cube_compare_z0.txt").string());
}

static void general_plane_demo(std::vector<DemoResult> &results)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  
  
  
  auto hex = std::make_shared<HexPrism>("gp_hex", 1.0, 0.0, 2.0);
  auto plane = std::make_shared<GeneralPlane>("gp", 0.0, 1.0, 0.0, 0.0); 
  Cell gp{"gplane", and_node(surf_node(hex), surf_node(plane))};

  Universe u;
  u.cells = {gp};
  BoundingBox b{{-1.5, -1.5, 0.0}, {1.5, 1.5, 2.0}};

  double hex_vol = 2.0 * std::sqrt(3.0) * 1.0 * 1.0 * 2.0;
  double analytic = 0.5 * hex_vol;
  run_volume_case("hex_halfspace_y_le_0", gp, b, analytic, results);
  ascii_plot(u, b, 140, 140, 1.0, (resultsDir / "ex3_plot_hex_halfplane_z1.txt").string());
}

static void torus_demo(std::vector<DemoResult> &results)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  double R = 1.5;
  double r = 0.5;
  auto torus = std::make_shared<Torus>("torus", R, r);
  Cell torus_cell{"torus", surf_node(torus)};
  Universe u;
  u.cells = {torus_cell};

  BoundingBox b{{-R - r, -R - r, -r}, {R + r, R + r, r}};
  double analytic = 2.0 * M_PI * M_PI * R * r * r;
  run_volume_case("torus", torus_cell, b, analytic, results);
  ascii_plot(u, b, 180, 180, 0.0, (resultsDir / "ex3_plot_torus_z0.txt").string());
  ascii_plot_xz(u, b, 220, 120, 0.0, (resultsDir / "ex3_plot_torus_xz_y0.txt").string());
}

static void elliptical_torus_demo(std::vector<DemoResult> &results)
{
  const fs::path resultsDir = Paths::resultsDir("ex3");
  ensure_dir(resultsDir);

  
  double R = 1.5;
  
  double a = 0.8;
  double b = 0.2;
  auto et = std::make_shared<EllipticalTorus>("etor", R, a, b);
  Cell et_cell{"etor", surf_node(et)};
  Universe u;
  u.cells = {et_cell};

  BoundingBox bb{{-(R + a), -(R + a), -b}, {R + a, R + a, b}};
  
  double analytic = 2.0 * M_PI * M_PI * R * a * b;
  run_volume_case("elliptical_torus", et_cell, bb, analytic, results);
  ascii_plot(u, bb, 180, 180, 0.0, (resultsDir / "ex3_plot_elliptical_torus_z0.txt").string());
  ascii_plot_xz(u, bb, 220, 120, 0.0, (resultsDir / "ex3_plot_elliptical_torus_xz_y0.txt").string());
}

void runExercise3()
{
  try
  {
    const fs::path resultsDir = Paths::resultsDir("ex3");
    std::cout << "Running Exercise 3 geometry routines...\n";
    ensure_dir(resultsDir);
    std::vector<DemoResult> results;

    Universe pin = build_fuel_pin();
    const double pitch = 1.330;
    const double half_pitch = pitch / 2.0;
    BoundingBox pin_box{{-half_pitch, -half_pitch, 0.0}, {half_pitch, half_pitch, 10.0}};
    const double r_fuel = 0.824 / 2.0;
    const double r_clad = r_fuel + 0.063;
    const double h_pin = 10.0;
    double pin_fuel_analytic = M_PI * r_fuel * r_fuel * h_pin;
    run_volume_case("fuel_pin_fuel", pin.cells[0], pin_box, pin_fuel_analytic, results);
    double pin_clad_analytic = M_PI * (r_clad * r_clad - r_fuel * r_fuel) * h_pin;
    run_volume_case("fuel_pin_clad", pin.cells[1], pin_box, pin_clad_analytic, results);
    double pin_coolant_analytic = (pitch * pitch - M_PI * r_clad * r_clad) * h_pin;
    run_volume_case("fuel_pin_coolant", pin.cells[2], pin_box, pin_coolant_analytic, results);
    ascii_plot(pin, pin_box, 120, 120, 5.0, (resultsDir / "ex3_plot_pin_z5.txt").string());

    Universe hollow = build_hollow_cylinder();
    BoundingBox hollow_box{{-2.2, -2.2, 0.0}, {2.2, 2.2, 7.0}};
    double steel_vol = M_PI * (2.05 * 2.05 - 2.0 * 2.0) * 7.0;
    double water_vol = M_PI * 2.0 * 2.0 * 3.5;
    double air_vol = M_PI * 2.0 * 2.0 * 3.5;
    run_volume_case("hollow_steel", hollow.cells[2], hollow_box, steel_vol, results);
    run_volume_case("hollow_water", hollow.cells[0], hollow_box, water_vol, results);
    run_volume_case("hollow_air", hollow.cells[1], hollow_box, air_vol, results);
    double outside_vol = box_volume(hollow_box) - M_PI * 2.05 * 2.05 * 7.0;
    run_volume_case("hollow_outside", hollow.cells[3], hollow_box, outside_vol, results);
    ascii_plot(hollow, hollow_box, 160, 160, 1.0, (resultsDir / "ex3_plot_hollow_z1.txt").string());
    ascii_plot(hollow, hollow_box, 160, 160, 5.0, (resultsDir / "ex3_plot_hollow_z5.txt").string());

    Universe hex = build_hex_prism_demo();
    BoundingBox hex_box{{-1.5, -1.5, 0.0}, {1.5, 1.5, 2.0}};
    
    double hex_analytic = 2.0 * std::sqrt(3.0) * 1.0 * 1.0 * 2.0;
    run_volume_case("hex_prism", hex.cells[0], hex_box, hex_analytic, results);
    ascii_plot(hex, hex_box, 120, 120, 1.0, (resultsDir / "ex3_plot_hex.txt").string());

    
    general_plane_demo(results);
    transform_demo(results);
    torus_demo(results);
    elliptical_torus_demo(results);

    write_results(results);
    write_split_results(results);
    write_volume_estimator_settings();

    
    lattice_demo();
    hex_lattice_demo();
    std::cout << "Exercise 3 outputs written to results/ex3/.\n";
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 3 error: " << e.what() << "\n";
  }
}
