#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP
#include <vector>
#include <cmath>
#include <algorithm>


inline double log_interp(const std::vector<double> &x, const std::vector<double> &y, double xq)
{
  if (xq <= x.front())
    return y.front();
  if (xq >= x.back())
    return y.back();
  auto it = std::lower_bound(x.begin(), x.end(), xq);
  size_t idx = static_cast<size_t>(std::distance(x.begin(), it));
  if (idx == 0)
    return y.front();
  size_t i0 = idx - 1;
  size_t i1 = idx;
  double x0 = x[i0], x1 = x[i1];
  double y0 = y[i0], y1 = y[i1];
  double lx0 = std::log(x0), lx1 = std::log(x1), lxq = std::log(xq);
  double ly0 = std::log(std::max(y0, 1e-30));
  double ly1 = std::log(std::max(y1, 1e-30));
  double t = (lxq - lx0) / (lx1 - lx0);
  double ly = ly0 + t * (ly1 - ly0);
  return std::exp(ly);
}

#endif
