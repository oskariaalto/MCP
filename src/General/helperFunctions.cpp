#include "helperFunctions.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <iomanip>

std::string expand_tilde(const std::string &path)
{
  if (path.empty() || path[0] != '~')
  {
    return path; 
  }

  
  const char *home = std::getenv("HOME");
  if (!home)
  {
    std::cerr << "Error: HOME environment variable is not set." << std::endl;
    return path; 
  }

  
  std::string full_path = std::string(home) + path.substr(1); 
  return full_path;
}

void writeVectorToCSV(const std::vector<long double> &x, const std::vector<long double> &y, const std::string &filename)
{
  if (x.size() != y.size())
  {
    std::cerr << "Error: Vectors x and y must have the same length." << std::endl;
    return;
  }

  
  std::string expanded_filename = expand_tilde(filename);

  std::ofstream file(expanded_filename, std::ios::out); 

  if (!file)
  {
    std::cerr << "Error: Could not open or create file " << expanded_filename << std::endl;
    return;
  }

  
  file << "x,y\n"; 

  
  for (size_t i = 0; i < x.size(); ++i)
  {
    file << x[i] << "," << y[i] << "\n"; 
  }

  file.close();
  std::cout << "Data successfully written to " << expanded_filename << std::endl;
}

void write_rows_to_csv(const std::vector<Row> &rows, const std::string &filename)
{
  const std::string expanded = expand_tilde(filename);

  std::ofstream csv(expanded, std::ios::out | std::ios::trunc);
  if (!csv)
  {
    std::cerr << "Error: could not open " << expanded << " for writing.\n";
    return;
  }
  csv << "M,N,mean,variance,stddev,ci_lower,ci_upper,abs_err,ci_contains_pi,time_s,skewness,kurtosis,fom\n";
  csv << std::setprecision(17);
  for (const auto &r : rows)
  {
    csv << r.M << ',' << r.N << ','
        << (double)r.mean << ','
        << (double)r.variance << ','
        << (double)r.stddev << ','
        << (double)r.ci_lower << ','
        << (double)r.ci_upper << ','
        << (double)r.abs_err << ','
        << (r.ci_ok ? 1 : 0) << ','
        << r.time_s << ','
        << (double)r.skewness << ','
        << (double)r.kurtosis << ','
        << r.fom << '\n';
  }
}
