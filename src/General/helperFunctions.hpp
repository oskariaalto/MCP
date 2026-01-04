#ifndef HELPER_FUNCTIONS_HPP
#define HELPER_FUNCTIONS_HPP
#include <string>
#include <vector>

struct Row
{
  int M, N;
  long double mean;
  long double variance;
  long double stddev;
  long double ci_lower;
  long double ci_upper;
  long double abs_err;
  bool ci_ok;
  double time_s;
  long double skewness;
  long double kurtosis;
  double fom;
};

std::string expand_tilde(const std::string &path);
void writeVectorToCSV(const std::vector<long double> &x, const std::vector<long double> &y, const std::string &filename);
void write_rows_to_csv(const std::vector<Row> &rows, const std::string &filename);

#endif