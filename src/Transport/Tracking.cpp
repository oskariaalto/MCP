#include "Tracking.hpp"
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>

static void collect_surfaces(const std::shared_ptr<CsgNode> &node, std::vector<SurfacePtr> &out)
{
  if (!node)
    return;
  if (node->op == CsgOp::SURFACE && node->surface)
  {
    out.push_back(node->surface);
  }
  else
  {
    collect_surfaces(node->left, out);
    collect_surfaces(node->right, out);
  }
}

double distance_to_exit_cell(const Cell &cell, const Ray &ray)
{
  double min_t = std::numeric_limits<double>::infinity();
  std::vector<SurfacePtr> surfs;
  collect_surfaces(cell.root, surfs);
  for (const auto &s : surfs)
  {
    double t = s->distance(ray);
    if (!std::isfinite(t))
      continue;
    Vec3 p_new = ray.origin + ray.dir * (t + 1e-7);
    bool inside_now = eval_node(*cell.root, ray.origin);
    bool inside_next = eval_node(*cell.root, p_new);
    if (inside_now && !inside_next)
    {
      if (t > 1e-8 && t < min_t)
        min_t = t;
    }
  }
  return min_t;
}

static const Cell *find_cell(const Universe &u, const Vec3 &p)
{
  return u.find_cell(p);
}

static const Material *material_for(const GeometryContext &geo, const Cell *cell)
{
  if (!cell)
    return nullptr;
  auto it = geo.cell_material.find(cell->name);
  if (it == geo.cell_material.end())
    return nullptr;
  return it->second;
}

static double majorant_sigma_total(const MaterialLibrary &matlib, double E)
{
  double maj = 0.0;
  for (const auto &kv : matlib.all())
  {
    maj = std::max(maj, kv.second.mix.macro_sigma_total(E));
  }
  return maj;
}

static double macro_sigma_absorption(const Mixture &mix, double E)
{
  double sum = 0.0;
  for (const auto &s : mix.species)
  {
    if (s.nuclide)
      sum += s.number_density * s.nuclide->sigma_absorption(E);
  }
  return sum;
}

static void score_track(TrackTallies &t, const std::string &mat, double tl, double st_macro, double sa_macro)
{
  auto &e = t.by_material[mat];
  e.track_length += tl;
  e.collisions_tle += tl * st_macro;
  e.absorptions_tle += tl * sa_macro;
}

TrackTallies simulate_external_source(const GeometryContext &geo, const MaterialLibrary &matlib, bool surface_tracking, int histories, MeshTally *mesh)
{
  TrackTallies tallies;
#pragma omp parallel
  {
    TrackTallies local;
#pragma omp for nowait
    for (int h = 0; h < histories; ++h)
    {
      struct Particle
      {
        Vec3 pos;
        Vec3 dir;
        double E = 0.0;
      };
      std::vector<Particle> stack;
      Vec3 src{0.0, 0.0, 10.0}; 
      stack.push_back({src, sample_isotropic_direction(), 1.0});
      while (!stack.empty())
      {
        Particle curr = stack.back();
        stack.pop_back();
        Vec3 pos = curr.pos;
        Vec3 dir = curr.dir;
        double E = curr.E;

        Ray ray{pos, dir.normalized()};
        const Cell *cell = find_cell(geo.universe, pos);
        if (!cell)
        {
          local.leaks += 1.0;
          continue;
        }
        const Material *mat = material_for(geo, cell);
        if (!mat)
        {
          local.leaks += 1.0;
          continue;
        }

        while (true)
        {
          double sigma_tot = mat->mix.macro_sigma_total(E);
          double sigma_abs = macro_sigma_absorption(mat->mix, E);
          if (sigma_tot <= 0.0)
          {
            local.leaks += 1.0;
            break;
          }

          double dist_col;
          if (surface_tracking)
          {
            double xi = std::max(1e-16, Random::instance().uniform());
            dist_col = -std::log(xi) / sigma_tot;
          }
          else
          {
            
            double sigma_maj = std::max(majorant_sigma_total(matlib, E), 1e-12);
            double xi = std::max(1e-16, Random::instance().uniform());
            dist_col = -std::log(xi) / sigma_maj;
          }

          double dist_surf = distance_to_exit_cell(*cell, ray);
          double step = dist_col;
          bool hit_boundary = false;
          if (dist_surf < dist_col)
          {
            step = dist_surf;
            hit_boundary = true;
          }

          score_track(local, mat->name, step, sigma_tot, sigma_abs);
          Vec3 new_pos = ray.origin + ray.dir * step;

          if (hit_boundary)
          {
            ray.origin = new_pos + ray.dir * 1e-6;
            const Cell *new_cell = find_cell(geo.universe, ray.origin);
            if (!new_cell)
            {
              local.leaks += 1.0;
              break;
            }
            cell = new_cell;
            mat = material_for(geo, cell);
            if (!mat)
            {
              local.leaks += 1.0;
              break;
            }
            continue;
          }

          
          if (!surface_tracking)
          {
            double sigma_tot_local = mat->mix.macro_sigma_total(E);
            double sigma_maj = std::max(majorant_sigma_total(matlib, E), 1e-12);

            
            
            
            if (sigma_maj > 0.0)
              local.by_material[mat->name].flux_cfe += 1.0 / sigma_maj;

            double p_accept = sigma_tot_local / sigma_maj;
            if (Random::instance().uniform() > p_accept)
            {
              ray.origin = new_pos;
              continue; 
            }
          }

          
          local.by_material[mat->name].collisions_analog += 1.0;
          if (surface_tracking && sigma_tot > 0.0)
            local.by_material[mat->name].flux_cfe += 1.0 / sigma_tot;
          if (mesh)
            mesh->score(new_pos);

          
          double sigma_macro = 0.0;
          double sigma_micro = 0.0;
          const Nuclide *nuclide = nullptr;
          {
            sigma_macro = mat->mix.macro_sigma_total(E);
            if (sigma_macro > 0.0)
            {
              double xi = Random::instance().uniform() * sigma_macro;
              double accum = 0.0;
              for (const auto &s : mat->mix.species)
              {
                if (!s.nuclide)
                  continue;
                double contrib = s.number_density * s.nuclide->sigma_total(E);
                accum += contrib;
                if (xi <= accum)
                {
                  sigma_micro = s.nuclide->sigma_total(E);
                  nuclide = s.nuclide;
                  break;
                }
              }
            }
          }
          if (!nuclide)
          {
            local.leaks += 1.0;
            break;
          }

          double sf = nuclide->sigma_fission(E);
          double sa = nuclide->sigma_capture(E);
          double s_total = sigma_micro;
          double xi_reac = Random::instance().uniform() * s_total;
          double sabs_other = std::max(0.0, nuclide->sigma_absorption(E) - sa - sf);
          if (xi_reac < sa + sabs_other)
          {
            local.by_material[mat->name].absorptions_analog += 1.0;
            break; 
          }
          else if (xi_reac < sa + sabs_other + sf)
          {
            local.by_material[mat->name].absorptions_analog += 1.0;
            int k = sample_fission_multiplicity(nuclide->nubar_value(E));
            for (int i = 0; i < k; ++i)
            {
              stack.push_back({new_pos, sample_isotropic_direction(), sample_fission_energy()});
            }
            break; 
          }
          else
          {
            double aw = (nuclide->atomic_weight > 0.0) ? nuclide->atomic_weight : static_cast<double>(nuclide->A);
            auto scat = scatter_elastic_free_gas(E, ray.dir, aw, nuclide->temperature);
            ray.origin = new_pos;
            ray.dir = scat.dir;
            E = scat.energy;
            continue;
          }
        }
      }
    }
#pragma omp critical
    {
      tallies.leaks += local.leaks;
      for (const auto &kv : local.by_material)
      {
        auto &dst = tallies.by_material[kv.first];
        dst.collisions_analog += kv.second.collisions_analog;
        dst.absorptions_analog += kv.second.absorptions_analog;
        dst.track_length += kv.second.track_length;
        dst.collisions_tle += kv.second.collisions_tle;
        dst.absorptions_tle += kv.second.absorptions_tle;
        dst.flux_cfe += kv.second.flux_cfe;
      }
    }
  }
  return tallies;
}

CriticalityResult simulate_criticality(const GeometryContext &geo, const MaterialLibrary &matlib, int histories, int generations)
{
  (void)matlib;
  CriticalityResult res;
  struct Particle
  {
    Vec3 pos;
    Vec3 dir;
    double E = 0.0;
  };
  std::vector<Particle> bank;
  for (int i = 0; i < histories; ++i)
    bank.push_back({{0.0, 0.0, 10.0}, sample_isotropic_direction(), 1.0});

  for (int g = 0; g < generations; ++g)
  {
    std::vector<Particle> next_bank;
    double produced = 0.0;
#pragma omp parallel
    {
      std::vector<Particle> local_bank;
      double local_produced = 0.0;
#pragma omp for nowait
      for (int h = 0; h < static_cast<int>(bank.size()); ++h)
      {
        Vec3 pos = bank[h].pos;
        Vec3 dir = bank[h].dir;
        double E = bank[h].E;
        Ray ray{pos, dir.normalized()};

        const Cell *cell = find_cell(geo.universe, pos);
        if (!cell)
          continue;
        const Material *mat = material_for(geo, cell);
        if (!mat)
          continue;

        while (true)
        {
          double sigma_tot = mat->mix.macro_sigma_total(E);
          if (sigma_tot <= 0.0)
            break;
          double xi = std::max(1e-16, Random::instance().uniform());
          double dist_col = -std::log(xi) / sigma_tot;
          double dist_surf = distance_to_exit_cell(*cell, ray);
          if (dist_surf < dist_col)
          {
            ray.origin = ray.origin + ray.dir * (dist_surf + 1e-6);
            cell = find_cell(geo.universe, ray.origin);
            if (!cell)
              break;
            mat = material_for(geo, cell);
            if (!mat)
              break;
            continue;
          }
          Vec3 new_pos = ray.origin + ray.dir * dist_col;
          
          double sigma_macro = mat->mix.macro_sigma_total(E);
          double sigma_micro = 0.0;
          const Nuclide *nuclide = nullptr;
          double xi_n = Random::instance().uniform() * sigma_macro;
          double accum = 0.0;
          for (const auto &s : mat->mix.species)
          {
            if (!s.nuclide)
              continue;
            double contrib = s.number_density * s.nuclide->sigma_total(E);
            accum += contrib;
            if (xi_n <= accum)
            {
              sigma_micro = s.nuclide->sigma_total(E);
              nuclide = s.nuclide;
              break;
            }
          }
          if (!nuclide)
            break;
          double sf = nuclide->sigma_fission(E);
          double sa = nuclide->sigma_capture(E);
          double s_total = sigma_micro;
          double xi_reac = Random::instance().uniform() * s_total;
          double sabs_other = std::max(0.0, nuclide->sigma_absorption(E) - sa - sf);
          if (xi_reac < sa + sabs_other)
          {
            break;
          }
          else if (xi_reac < sa + sabs_other + sf)
          {
            int k = sample_fission_multiplicity(nuclide->nubar_value(E));
            local_produced += k;
            for (int i = 0; i < k; ++i)
              local_bank.push_back({new_pos, sample_isotropic_direction(), sample_fission_energy()});
            break;
          }
          else
          {
            double aw = (nuclide->atomic_weight > 0.0) ? nuclide->atomic_weight : static_cast<double>(nuclide->A);
            auto scat = scatter_elastic_free_gas(E, ray.dir, aw, nuclide->temperature);
            E = scat.energy;
            ray.dir = scat.dir;
            ray.origin = new_pos;
            continue;
          }
        }
      }
#pragma omp critical
      {
        produced += local_produced;
        next_bank.insert(next_bank.end(), local_bank.begin(), local_bank.end());
      }
    }
    double keff_gen = (bank.empty() ? 0.0 : produced / bank.size());
    res.gen_keff.push_back(keff_gen);
    if (next_bank.empty())
      break;
    
    if (static_cast<int>(next_bank.size()) > histories)
      next_bank.resize(histories);
    bank.swap(next_bank);
  }
  if (!res.gen_keff.empty())
  {
    double sum = 0.0;
    for (double k : res.gen_keff)
      sum += k;
    res.keff = sum / res.gen_keff.size();
  }
  return res;
}
