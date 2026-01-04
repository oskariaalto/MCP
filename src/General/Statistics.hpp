#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <vector>
#include <string>
#include <ostream>
#include <chrono>
#include "helperFunctions.hpp"

class Statistics
{
private:
  std::vector<long double> data;            
  std::vector<long double> normalized_data; 
  long double mean;
  long double variance;
  long double std_dev;
  long double skewness;
  long double kurtosis;
  double time_elapsed; 
  double fom;          

  
  long double calculate_mean();
  long double calculate_welford_variance();
  long double calculate_skewness();
  long double calculate_kurtosis();

  
  void normalize_data();

  
  void update_statistics();

  
  void calculate_fom();

public:
  
  Statistics(const std::vector<long double> &input_data);

  
  void set_time_elapsed(double time);

  
  long double get_mean() const;
  long double get_variance() const;
  long double get_standard_deviation() const;
  long double get_skewness() const;
  long double get_kurtosis() const;
  double get_time_elapsed() const;
  double get_fom() const;

  
  long double jarque_bera_test() const;

  
  double get_p_value() const;

  
  bool check_confidence_interval(double expected_value) const;

  
  std::vector<long double> get_normalized_data() const;

  void write_to_bin_file(const std::string &filename) const;

  
  friend std::ostream &operator<<(std::ostream &os, const Statistics &stats);

  Row to_row(int M, int N, long double expected_value);
};

#endif 
