#ifndef RANDOM_HPP
#define RANDOM_HPP
#include <random>

class Random
{
public:
  static Random &instance()
  {
    static Random rng;
    return rng;
  }

  double uniform()
  {
    return dist_(engine_);
  }

  double sample_normal(double mean, double stddev)
  {
    std::normal_distribution<double> n(mean, stddev);
    return n(engine_);
  }

  
  
  double sample_maxwell(double a)
  {
    std::normal_distribution<double> n(0.0, std::sqrt(a));
    double x = n(engine_);
    double y = n(engine_);
    double z = n(engine_);
    return 0.5 * (x * x + y * y + z * z);
  }

  int sample_poisson(double mean)
  {
    std::poisson_distribution<int> p(mean);
    return p(engine_);
  }

private:
  Random() = default;
  static thread_local std::mt19937_64 engine_;
  static thread_local std::uniform_real_distribution<double> dist_;
};

#endif
