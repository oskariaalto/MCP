#ifndef VOLUME_ESTIMATOR_HPP
#define VOLUME_ESTIMATOR_HPP
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include "../Common/Random.hpp"
#include "CSG.hpp"

struct BoundingBox
{
  Vec3 min;
  Vec3 max;
};

struct EstimateResult
{
  double estimate = 0.0; 
  double std_err = 0.0;  
  double rel_err() const { return (estimate != 0.0) ? (std_err / std::abs(estimate)) : 0.0; }
};

inline Vec3 sample_point(const BoundingBox &b)
{
  auto &rng = Random::instance();
  double x = b.min.x + (b.max.x - b.min.x) * rng.uniform();
  double y = b.min.y + (b.max.y - b.min.y) * rng.uniform();
  double z = b.min.z + (b.max.z - b.min.z) * rng.uniform();
  return {x, y, z};
}

inline double box_volume(const BoundingBox &b)
{
  return (b.max.x - b.min.x) * (b.max.y - b.min.y) * (b.max.z - b.min.z);
}

inline double estimate_volume_points(const Cell &cell, const BoundingBox &b, int samples)
{
  int inside = 0;
#pragma omp parallel for reduction(+ : inside)
  for (int i = 0; i < samples; ++i)
  {
    Vec3 p = sample_point(b);
    if (eval_node(*cell.root, p))
      inside += 1;
  }
  return static_cast<double>(inside) / samples * box_volume(b);
}

inline double estimate_volume_lines(const Cell &cell, const BoundingBox &b, int lines, int segments_per_line)
{
  double total_len = 0.0;
  double inside_len = 0.0;
#pragma omp parallel for reduction(+ : total_len, inside_len)
  for (int i = 0; i < lines; ++i)
  {
    
    
    Vec3 p0 = sample_point(b);
    const double y = p0.y;
    const double z = p0.z;
    const double L = (b.max.x - b.min.x);
    if (L <= 0.0)
      continue;

    const double seg_len = L / segments_per_line;
    total_len += L;
    for (int s = 0; s < segments_per_line; ++s)
    {
      const double x = b.min.x + (s + Random::instance().uniform()) * seg_len;
      Vec3 p{x, y, z};
      if (eval_node(*cell.root, p))
        inside_len += seg_len;
    }
  }
  return (inside_len / total_len) * box_volume(b);
}

inline EstimateResult estimate_volume_points_stats(const Cell &cell, const BoundingBox &b, int samples)
{
  int inside = 0;
  for (int i = 0; i < samples; ++i)
  {
    Vec3 p = sample_point(b);
    if (eval_node(*cell.root, p))
      inside++;
  }
  double V = box_volume(b);
  double p_in = (samples > 0) ? (double)inside / (double)samples : 0.0;
  double est = p_in * V;

  
  double var_p = (samples > 0) ? (p_in * (1.0 - p_in) / (double)samples) : 0.0;
  double std_err = std::sqrt(std::max(0.0, var_p)) * V;

  return {est, std_err};
}

inline EstimateResult estimate_volume_lines_stats(const Cell &cell, const BoundingBox &b, int lines, int segments_per_line)
{
  
  
  
  

  std::vector<double> w; 
  std::vector<double> f; 
  w.reserve(lines);
  f.reserve(lines);

  for (int i = 0; i < lines; ++i)
  {
    Vec3 start = sample_point(b);
    Vec3 end = sample_point(b);
    double L = (end - start).norm();
    if (L <= 0.0)
      continue;

    Vec3 dir = (end - start) / segments_per_line;
    double seg_len = L / segments_per_line;
    double inside_len = 0.0;
    for (int s = 0; s < segments_per_line; ++s)
    {
      Vec3 p = start + dir * (s + 0.5);
      if (eval_node(*cell.root, p))
        inside_len += seg_len;
    }

    w.push_back(L);
    f.push_back(inside_len / L);
  }

  double W = 0.0;
  double WF = 0.0;
  for (size_t i = 0; i < w.size(); ++i)
  {
    W += w[i];
    WF += w[i] * f[i];
  }
  double frac = (W > 0.0) ? (WF / W) : 0.0;

  
  double var_w = 0.0;
  double W2 = 0.0;
  for (size_t i = 0; i < w.size(); ++i)
  {
    var_w += w[i] * (f[i] - frac) * (f[i] - frac);
    W2 += w[i] * w[i];
  }
  var_w = (W > 0.0) ? (var_w / W) : 0.0;
  double n_eff = (W2 > 0.0) ? (W * W / W2) : 0.0;

  double std_err_frac = (n_eff > 1.0) ? std::sqrt(std::max(0.0, var_w / n_eff)) : 0.0;

  double V = box_volume(b);
  return {frac * V, std_err_frac * V};
}

#endif
