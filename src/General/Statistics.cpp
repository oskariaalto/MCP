#include "Statistics.hpp"
#include <numeric>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>


double chi_squared_cdf(double x)
{
  
  return 1.0 - std::exp(-x / 2.0);
}


Statistics::Statistics(const std::vector<long double> &input_data) : data(input_data), time_elapsed(0), fom(0)
{
  update_statistics();
}


void Statistics::set_time_elapsed(double time)
{
  time_elapsed = time;
  calculate_fom();
}


void Statistics::calculate_fom()
{
  fom = 1.0 / (std_dev * std_dev * time_elapsed);
}


long double Statistics::calculate_mean()
{
  return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}


long double Statistics::calculate_welford_variance()
{
  long double M2 = 0.0;
  long double N = data.size();

  for (const auto &value : data)
  {
    M2 += std::pow(value, 2);
  }
  return 1 / (N * (N - 1)) * (M2 - 1 / N * mean * mean);
}


long double Statistics::calculate_skewness()
{
  long double m3_sum = 0.0;
  long double N = data.size();
  for (const auto &x : data)
  {
    m3_sum += std::pow(x - mean, 3);
  }
  long double m3 = m3_sum / N;

  
  long double skewness = std::sqrt((N * (N - 1))) / (N - 2) * (m3 / std::pow(variance, 3 / 2));

  return skewness;
}


long double Statistics::calculate_kurtosis()
{
  long double m4_sum = 0.0;
  long double N = data.size();

  
  for (const auto &x : data)
  {
    m4_sum += std::pow(x - mean, 4);
  }

  long double m4 = m4_sum / N; 

  
  long double numerator = (N * (N + 1)) / ((N - 1) * (N - 2) * (N - 3)) * (m4 / std::pow(variance, 2));
  long double correction = 3 * std::pow(N - 1, 2) / ((N - 2) * (N - 3));

  return correction - numerator;
}


void Statistics::normalize_data()
{
  normalized_data.clear();
  for (const auto &val : data)
  {
    normalized_data.push_back((val - mean) / std_dev);
  }
}


void Statistics::update_statistics()
{
  mean = calculate_mean();
  variance = calculate_welford_variance();
  std_dev = std::sqrt(variance);
  normalize_data(); 
  skewness = calculate_skewness();
  kurtosis = calculate_kurtosis();
}


long double Statistics::jarque_bera_test() const
{
  long double S = skewness;
  long double K = kurtosis;
  long double n = static_cast<long double>(data.size());

  return n / 6.0 * (std::pow((double)S, 2) + 0.25L * std::pow((double)(K - 3), 2));
}


double Statistics::get_p_value() const
{
  long double jb_stat = jarque_bera_test();
  
  return 1.0 - chi_squared_cdf(jb_stat);
}


bool Statistics::check_confidence_interval(double expected_value) const
{
  double z = 1.96; 

  
  double margin_of_error = z * (std_dev);

  
  double lower_bound = mean - margin_of_error;
  double upper_bound = mean + margin_of_error;

  std::cout << "Confidence Interval: [" << lower_bound << ", " << upper_bound << "]\n";

  
  return expected_value >= lower_bound && expected_value <= upper_bound;
}


double Statistics::get_time_elapsed() const
{
  return time_elapsed;
}


double Statistics::get_fom() const
{
  return fom;
}


long double Statistics::get_mean() const
{
  return mean;
}


long double Statistics::get_variance() const
{
  return variance;
}


long double Statistics::get_standard_deviation() const
{
  return std_dev;
}


long double Statistics::get_skewness() const
{
  return skewness;
}


long double Statistics::get_kurtosis() const
{
  return kurtosis;
}


std::vector<long double> Statistics::get_normalized_data() const
{
  return normalized_data;
}


std::ostream &operator<<(std::ostream &os, const Statistics &stats)
{
  os << "Mean: " << stats.mean << "\n"
     << "Variance: " << stats.variance << "\n"
     << "Standard Deviation: " << stats.std_dev << "\n"
     << "Skewness: " << stats.skewness << "\n"
     << "Kurtosis: " << stats.kurtosis << "\n"
     << "Jarque-Bera Test Statistic: " << stats.jarque_bera_test() << "\n"
     << "P-Value: " << stats.get_p_value() << "\n" 
     << "Time Elapsed: " << stats.time_elapsed << " seconds\n"
     << "Figure of Merit (FOM): " << stats.fom << "\n\n";
  return os;
}

void Statistics::write_to_bin_file(const std::string &filename) const
{
  std::string full_filename = expand_tilde(filename);
  
  std::ofstream out_file(full_filename, std::ios::binary | std::ios::trunc);

  if (!out_file)
  {
    std::cerr << "Unable to open file " << full_filename << " for writing." << std::endl;
    return;
  }

  
  size_t size = data.size();
  out_file.write(reinterpret_cast<const char *>(&size), sizeof(size));

  
  out_file.write(reinterpret_cast<const char *>(data.data()), size * sizeof(long double));

  out_file.close(); 
}

Row Statistics::to_row(int M, int N, long double expected_value)
{
  Row r{};
  r.M = M;
  r.N = N;
  r.mean = mean;
  r.variance = variance;
  r.stddev = std_dev;
  r.ci_lower = mean - 1.96L * std_dev / std::sqrt((long double)data.size());
  r.ci_upper = mean + 1.96L * std_dev / std::sqrt((long double)data.size());
  r.abs_err = fabsl(mean - expected_value);
  r.ci_ok = (r.ci_lower <= expected_value && expected_value <= r.ci_upper);
  r.time_s = time_elapsed;
  r.skewness = skewness;
  r.kurtosis = kurtosis;
  r.fom = fom;
  return r;
}
