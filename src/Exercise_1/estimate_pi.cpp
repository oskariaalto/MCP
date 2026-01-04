#include "estimate_pi.hpp"
#include <omp.h>
#include <random>
#include <vector>
#include <chrono>
#include <cstdint>
#include <limits>


static inline uint64_t splitmix64(uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

typedef unsigned u4_t __attribute__((vector_size(4 * sizeof(unsigned))));

constexpr u4_t u4zero{0, 0, 0, 0};

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));

constexpr double4_t d4zero{0.0, 0.0, 0.0, 0.0};

Statistics random_point_method(double M, double N)
{
  const int n = static_cast<int>(N);
  const int m = static_cast<int>(M);

  std::vector<long double> values(n, 0.0L);

  auto start_time = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
    
    std::random_device rd;
    uint64_t seed = splitmix64(((uint64_t)rd() << 32) ^ (uint64_t)omp_get_thread_num());
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

#pragma omp for schedule(dynamic)
    for (int i = 0; i < n; ++i)
    {
      u4_t hits = u4zero;
      int j = 0;
      for (; (j + 3) < m; j += 4)
      {
        double4_t x = d4zero;
        double4_t y = d4zero;
        x[0] = dis(gen);
        x[1] = dis(gen);
        x[2] = dis(gen);
        x[3] = dis(gen);

        y[0] = dis(gen);
        y[1] = dis(gen);
        y[2] = dis(gen);
        y[3] = dis(gen);

        const double4_t r2 = x * x + y * y;

        hits[0] += static_cast<unsigned>(r2[0] < 1.0);
        hits[1] += static_cast<unsigned>(r2[1] < 1.0);
        hits[2] += static_cast<unsigned>(r2[2] < 1.0);
        hits[3] += static_cast<unsigned>(r2[3] < 1.0);
      }
      for (; j < m; ++j)
      {
        const double x = dis(gen);
        const double y = dis(gen);
        hits[0] += (x * x + y * y < 1.0);
      }

      const unsigned total_hits = hits[0] + hits[1] + hits[2] + hits[3];
      values[i] = (4.0L * static_cast<long double>(total_hits)) / static_cast<long double>(m);
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  double time_elapsed = std::chrono::duration<double>(end_time - start_time).count();

  Statistics stats(values);
  stats.set_time_elapsed(time_elapsed);
  return stats;
}

Statistics buffons_method(double M, double N, double t, double l)
{
  const int n = static_cast<int>(N);
  const int m = static_cast<int>(M);

  std::vector<long double> values(n, 0.0L);

  const long double num_factor = 2.0L * static_cast<long double>(l) * static_cast<long double>(m);
  const long double half_l = static_cast<long double>(l) / 2.0L;

  auto start_time = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
    std::random_device rd;
    uint64_t seed_base = splitmix64(((uint64_t)rd() << 32) ^ (uint64_t)omp_get_thread_num());

    std::mt19937_64 gen_theta(seed_base);
    std::mt19937_64 gen_y(splitmix64(seed_base ^ 0x12345678ULL));

    std::uniform_real_distribution<long double> dis_theta(0.0L, static_cast<long double>(PI));
    std::uniform_real_distribution<long double> dis_y(0.0L, static_cast<long double>(t) / 2.0L);

#pragma omp for schedule(dynamic)
    for (int i = 0; i < n; ++i)
    {
      u4_t hits = u4zero;
      int j = 0;
      for (; (j + 3) < m; j += 4)
      {
        const long double th0 = dis_theta(gen_theta);
        const long double th1 = dis_theta(gen_theta);
        const long double th2 = dis_theta(gen_theta);
        const long double th3 = dis_theta(gen_theta);

        const long double y0 = dis_y(gen_y);
        const long double y1 = dis_y(gen_y);
        const long double y2 = dis_y(gen_y);
        const long double y3 = dis_y(gen_y);

        
        if (y0 <= half_l * std::sin((double)th0))
          ++hits[0];
        if (y1 <= half_l * std::sin((double)th1))
          ++hits[1];
        if (y2 <= half_l * std::sin((double)th2))
          ++hits[2];
        if (y3 <= half_l * std::sin((double)th3))
          ++hits[3];
      }
      for (; j < m; ++j)
      {
        const double theta = dis_theta(gen_theta);
        const double y = dis_y(gen_y);
        if (y <= half_l * std::sin((double)theta))
          ++hits[0];
      }

      const unsigned total_hits = hits[0] + hits[1] + hits[2] + hits[3];
      if (total_hits == 0)
      {
        values[i] = std::numeric_limits<long double>::quiet_NaN();
      }
      else
      {
        values[i] = (num_factor / static_cast<long double>(t)) / static_cast<long double>(total_hits);
      }
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  double time_elapsed = std::chrono::duration<double>(end_time - start_time).count();

  Statistics stats(values);
  stats.set_time_elapsed(time_elapsed);
  return stats;
}