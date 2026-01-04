#include "Neutron.hpp"
#include <numeric>

static double scatter_energy_free_gas(double E_in, const Nuclide *nuclide)
{
  if (!nuclide)
    return E_in;
  double aw = (nuclide->atomic_weight > 0.0) ? nuclide->atomic_weight : static_cast<double>(nuclide->A);
  return scatter_elastic_free_gas(E_in, Vec3{0.0, 0.0, 1.0}, aw, nuclide->temperature).energy;
}

static const Nuclide *sample_nuclide(const Mixture &mix, double E, double &sigma_total_macro, double &sigma_total_micro)
{
  sigma_total_macro = mix.macro_sigma_total(E);
  if (sigma_total_macro <= 0.0)
    return nullptr;
  double xi = Random::instance().uniform() * sigma_total_macro;
  double accum = 0.0;
  for (const auto &s : mix.species)
  {
    double contrib = s.number_density * s.nuclide->sigma_total(E);
    accum += contrib;
    if (xi <= accum)
    {
      sigma_total_micro = s.nuclide->sigma_total(E);
      return s.nuclide;
    }
  }
  sigma_total_micro = mix.species.back().nuclide->sigma_total(E);
  return mix.species.back().nuclide;
}

static double sample_path_length(double sigma_tot_macro)
{
  double xi = std::max(1e-16, Random::instance().uniform());
  return -std::log(xi) / sigma_tot_macro;
}

static double sum_inelastic(const Nuclide *n, double E)
{
  double sum = 0.0;
  for (const auto &kv : n->reactions)
  {
    int mt = kv.first;
    if (mt >= 50 && mt < 90) 
    {
      sum += kv.second.table.eval(E);
    }
  }
  return sum;
}

static double inelastic_q(const Nuclide *n, int mt)
{
  auto it = n->reactions.find(mt);
  if (it == n->reactions.end())
    return 0.0;
  return std::abs(it->second.q_value);
}

static double sigma_abs_other(const Nuclide *n, double E)
{
  if (!n)
    return 0.0;
  double sa = n->sigma_capture(E);
  double sf = n->sigma_fission(E);
  double sabs = n->sigma_absorption(E);
  return std::max(0.0, sabs - sa - sf);
}

MultiplicationResult simulate_multiplication(const Mixture &mix, double source_energy, int histories, bool include_inelastic)
{
  std::vector<double> produced(histories, 0.0);

#pragma omp parallel for
  for (int h = 0; h < histories; ++h)
  {
    std::vector<Neutron> stack;
    stack.push_back({source_energy});
    int total_children = 0;
    int steps = 0;
    const int max_steps = 2000;

    while (!stack.empty())
    {
      if (++steps > max_steps)
        break;
      Neutron n = stack.back();
      stack.pop_back();
      double sigma_macro = 0.0;
      double sigma_micro = 0.0;
      const Nuclide *nuclide = sample_nuclide(mix, n.energy, sigma_macro, sigma_micro);
      if (!nuclide || sigma_macro <= 0.0)
        continue;

      
      (void)sample_path_length(sigma_macro);

      double sf = nuclide->sigma_fission(n.energy);
      double sa = nuclide->sigma_capture(n.energy);
      double sabs_other = sigma_abs_other(nuclide, n.energy);
      
      
      double sinel_total = sum_inelastic(nuclide, n.energy);
      double sinel = include_inelastic ? sinel_total : 0.0;
      double s_total = include_inelastic ? sigma_micro : std::max(0.0, sigma_micro - sinel_total);
      if (s_total <= 0.0)
        continue;

      
      double s_elastic = std::max(0.0, s_total - (sa + sabs_other + sf + sinel));

      double xi = Random::instance().uniform() * s_total;
      if (xi < sa + sabs_other)
      {
        continue; 
      }
      else if (xi < sa + sabs_other + sf)
      {
        double nubar = nuclide->nubar_value(n.energy);
        int k = sample_fission_multiplicity(nubar);
        total_children += k;
        for (int i = 0; i < k; ++i)
        {
          stack.push_back({sample_fission_energy()});
        }
      }
      else if (include_inelastic && xi < sa + sabs_other + sf + sinel)
      {
        
        
        double xi2 = Random::instance().uniform() * sinel;
        double acc = 0.0;
        int chosen_mt = 0;
        for (const auto &kv : nuclide->reactions)
        {
          int mt = kv.first;
          if (mt >= 50 && mt < 90)
          {
            double xs = kv.second.table.eval(n.energy);
            acc += xs;
            if (xi2 <= acc)
            {
              chosen_mt = mt;
              break;
            }
          }
        }
        double e_scatter = scatter_energy_free_gas(n.energy, nuclide) - inelastic_q(nuclide, chosen_mt);
        if (e_scatter > 1e-8)
          stack.push_back({e_scatter});
      }
      else
      {
        
        if (s_elastic <= 0.0)
          continue;
        double e_scatter = scatter_energy_free_gas(n.energy, nuclide);
        if (e_scatter > 1e-8)
          stack.push_back({e_scatter});
      }
    }

    produced[h] = static_cast<double>(total_children);
  }

  double mean = std::accumulate(produced.begin(), produced.end(), 0.0) / produced.size();
  double var = 0.0;
  for (double x : produced)
    var += (x - mean) * (x - mean);
  
  if (produced.size() > 1)
    var /= (produced.size() - 1);
  else
    var = 0.0;
  return {mean, std::sqrt(var)};
}

SlowingDownResult simulate_slowing_down(const Mixture &mix, double source_energy, int histories, int max_collisions)
{
  std::vector<double> energy_sums(max_collisions, 0.0);
  std::vector<int> counts(max_collisions, 0);

#pragma omp parallel for
  for (int h = 0; h < histories; ++h)
  {
    std::vector<double> local_sums(max_collisions, 0.0);
    std::vector<int> local_counts(max_collisions, 0);

    Neutron n{source_energy};
    int steps = 0;
    const int max_steps = max_collisions * 2;
    for (int c = 0; c < max_collisions; ++c)
    {
      if (++steps > max_steps)
        break;
      double sigma_macro = 0.0;
      double sigma_micro = 0.0;
      const Nuclide *nuclide = sample_nuclide(mix, n.energy, sigma_macro, sigma_micro);
      if (!nuclide || sigma_macro <= 0.0)
        break;

      (void)sample_path_length(sigma_macro);

      double sf = nuclide->sigma_fission(n.energy);
      double sa = nuclide->sigma_capture(n.energy);
      double s_total = sigma_micro;
      double sabs_other = sigma_abs_other(nuclide, n.energy);
      double s_scat = std::max(0.0, s_total - sa - sabs_other - sf);

      double xi = Random::instance().uniform() * s_total;
      if (xi < sa + sabs_other || xi >= sa + sabs_other + sf + s_scat)
        break; 
      if (xi < sa + sabs_other + sf)
      {
        
        break;
      }

      n.energy = scatter_energy_free_gas(n.energy, nuclide);
      local_sums[c] += n.energy;
      local_counts[c] += 1;
    }

    for (int c = 0; c < max_collisions; ++c)
    {
      if (local_counts[c] > 0)
      {
#pragma omp atomic
        energy_sums[c] += local_sums[c];
#pragma omp atomic
        counts[c] += local_counts[c];
      }
    }
  }

  SlowingDownResult res;
  for (int i = 0; i < max_collisions; ++i)
  {
    if (counts[i] == 0)
      break;
    res.collision_index.push_back(i + 1);
    res.average_energy.push_back(energy_sums[i] / counts[i]);
  }
  return res;
}

std::vector<ReactionTally> tally_reactions(const Mixture &mix, double source_energy, int histories, bool include_inelastic)
{
  
  
  std::vector<std::unordered_map<std::string, double>> per_history(histories);

#pragma omp parallel for
  for (int h = 0; h < histories; ++h)
  {
    std::unordered_map<std::string, double> local;
    std::vector<Neutron> stack;
    stack.push_back({source_energy});
    int steps = 0;
    const int max_steps = 2000;
    while (!stack.empty())
    {
      if (++steps > max_steps)
        break;
      Neutron n = stack.back();
      stack.pop_back();
      double sigma_macro = 0.0;
      double sigma_micro = 0.0;
      const Nuclide *nuclide = sample_nuclide(mix, n.energy, sigma_macro, sigma_micro);
      if (!nuclide || sigma_macro <= 0.0)
        continue;

      (void)sample_path_length(sigma_macro);

      double sf = nuclide->sigma_fission(n.energy);
      double sa = nuclide->sigma_capture(n.energy);
      double sabs_other = sigma_abs_other(nuclide, n.energy);
      double sinel_total = sum_inelastic(nuclide, n.energy);
      double sinel = include_inelastic ? sinel_total : 0.0;
      double s_total = include_inelastic ? sigma_micro : std::max(0.0, sigma_micro - sinel_total);
      if (s_total <= 0.0)
        continue;
      double s_elastic = std::max(0.0, s_total - (sa + sabs_other + sf + sinel));

      double xi = Random::instance().uniform() * s_total;
      if (xi < sa)
      {
        local[nuclide->name + ":capture"] += 1.0;
        continue;
      }
      else if (xi < sa + sabs_other)
      {
        local[nuclide->name + ":abs_other"] += 1.0;
        continue;
      }
      else if (xi < sa + sabs_other + sf)
      {
        local[nuclide->name + ":fission"] += 1.0;
        int k = sample_fission_multiplicity(nuclide->nubar_value(n.energy));
        for (int i = 0; i < k; ++i)
          stack.push_back({sample_fission_energy()});
      }
      else if (include_inelastic && xi < sa + sabs_other + sf + sinel)
      {
        local[nuclide->name + ":inelastic"] += 1.0;
        double xi2 = Random::instance().uniform() * sinel;
        double acc = 0.0;
        int chosen_mt = 0;
        for (const auto &kv : nuclide->reactions)
        {
          int mt = kv.first;
          if (mt >= 50 && mt < 90)
          {
            double xs = kv.second.table.eval(n.energy);
            acc += xs;
            if (xi2 <= acc)
            {
              chosen_mt = mt;
              break;
            }
          }
        }
        double e_scatter = scatter_energy_free_gas(n.energy, nuclide) - inelastic_q(nuclide, chosen_mt);
        if (e_scatter > 1e-8)
          stack.push_back({e_scatter});
      }
      else
      {
        if (s_elastic <= 0.0)
          continue;
        local[nuclide->name + ":elastic"] += 1.0;
        double e_scatter = scatter_energy_free_gas(n.energy, nuclide);
        if (e_scatter > 1e-8)
          stack.push_back({e_scatter});
      }
    }
    per_history[h] = std::move(local);
  }

  
  std::unordered_map<std::string, std::vector<double>> counts;
  counts.reserve(mix.species.size() * 4);
  for (int h = 0; h < histories; ++h)
  {
    for (const auto &kv : per_history[h])
    {
      auto &vec = counts[kv.first];
      if (vec.empty())
        vec.assign(histories, 0.0);
      vec[h] = kv.second;
    }
  }

  std::vector<ReactionTally> out;
  out.reserve(counts.size());
  for (auto &kv : counts)
  {
    const auto &vals = kv.second;
    const int N = static_cast<int>(vals.size());
    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / N;

    
    double ss = 0.0;
    for (double v : vals)
      ss += (v - mean) * (v - mean);
    double var = (N > 1) ? (ss / (N - 1)) : 0.0;
    double std_err = (N > 0) ? std::sqrt(var / N) : 0.0;
    double rel_err = mean > 0 ? (std_err / mean) : 0.0;

    auto pos = kv.first.find(':');
    ReactionTally r;
    r.nuclide = kv.first.substr(0, pos);
    r.reaction = kv.first.substr(pos + 1);
    r.mean = mean;
    r.rel_err = rel_err;
    out.push_back(r);
  }
  return out;
}

TimeSeries simulate_time_dependent(const Mixture &mix, double source_energy, int histories, double t_max, int bins, bool include_inelastic)
{
  TimeSeries ts;
  ts.bin_edges.resize(bins + 1);
  ts.counts.assign(bins, 0.0);
  double dt = t_max / bins;
  for (int i = 0; i <= bins; ++i)
    ts.bin_edges[i] = i * dt;

#pragma omp parallel
  {
    std::vector<double> local_counts(bins, 0.0);
#pragma omp for nowait
    for (int h = 0; h < histories; ++h)
    {
      std::vector<std::pair<Neutron, double>> stack; 
      stack.push_back({{source_energy}, 0.0});
      int steps = 0;
      const int max_steps = 4000;
      while (!stack.empty())
      {
        if (++steps > max_steps)
          break;
        auto current = stack.back();
        stack.pop_back();
        Neutron n = current.first;
        double t = current.second;

        double sigma_macro = 0.0;
        double sigma_micro = 0.0;
        const Nuclide *nuclide = sample_nuclide(mix, n.energy, sigma_macro, sigma_micro);
        if (!nuclide || sigma_macro <= 0.0)
          continue;

        double path = sample_path_length(sigma_macro);
        
        double mn_MeV = 939.565; 
        double beta = std::sqrt(std::max(0.0, 2.0 * n.energy / mn_MeV));
        double c_cm_s = 2.99792458e10;
        double dt_step = path / (beta * c_cm_s);
        t += dt_step;
        if (t > t_max)
          continue;

        double sf = nuclide->sigma_fission(n.energy);
        double sa = nuclide->sigma_capture(n.energy);
        double sabs_other = sigma_abs_other(nuclide, n.energy);
        double sinel_total = sum_inelastic(nuclide, n.energy);
        double sinel = include_inelastic ? sinel_total : 0.0;
        double s_total = include_inelastic ? sigma_micro : std::max(0.0, sigma_micro - sinel_total);
        if (s_total <= 0.0)
          continue;

        double s_elastic = std::max(0.0, s_total - (sa + sabs_other + sf + sinel));

        double xi = Random::instance().uniform() * s_total;
        if (xi < sa + sabs_other)
        {
          continue;
        }
        else if (xi < sa + sabs_other + sf)
        {
          int k = sample_fission_multiplicity(nuclide->nubar_value(n.energy));
          if (k > 0)
          {
            int bin = static_cast<int>(t / dt);
            if (bin >= 0 && bin < bins)
              local_counts[bin] += k;
            for (int i = 0; i < k; ++i)
              stack.push_back({{sample_fission_energy()}, t});
          }
        }
        else if (include_inelastic && xi < sa + sabs_other + sf + sinel)
        {
          double xi2 = Random::instance().uniform() * sinel;
          double acc = 0.0;
          int chosen_mt = 0;
          for (const auto &kv : nuclide->reactions)
          {
            int mt = kv.first;
            if (mt >= 50 && mt < 90)
            {
              double xs = kv.second.table.eval(n.energy);
              acc += xs;
              if (xi2 <= acc)
              {
                chosen_mt = mt;
                break;
              }
            }
          }
          double e_scatter = scatter_energy_free_gas(n.energy, nuclide) - inelastic_q(nuclide, chosen_mt);
          if (e_scatter > 1e-8)
            stack.push_back({{e_scatter}, t});
        }
        else
        {
          if (s_elastic <= 0.0)
            continue;
          double e_scatter = scatter_energy_free_gas(n.energy, nuclide);
          if (e_scatter > 1e-8)
            stack.push_back({{e_scatter}, t});
        }
      }
    }
#pragma omp critical
    {
      for (int i = 0; i < bins; ++i)
        ts.counts[i] += local_counts[i];
    }
  }
  return ts;
}
