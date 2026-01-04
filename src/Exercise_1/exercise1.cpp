#include "exercise1.hpp"
#include "../Common/Paths.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <cmath>

namespace fs = std::filesystem;

static void ensure_dir(const fs::path &p)
{
  fs::create_directories(p);
}

static int nearest_divisor(long long S, long long target)
{
  if (target < 1)
    target = 1;
  for (long long d = 0;; ++d)
  {
    if (target - d >= 1 && S % (target - d) == 0)
      return (int)(target - d);
    if (S % (target + d) == 0)
      return (int)(target + d);
  }
}

void runExercise1()
{
  fs::path root;
  try
  {
    root = Paths::repoRoot();
  }
  catch (const std::exception &e)
  {
    std::cout << "Exercise 1 error: " << e.what() << "\n";
    return;
  }

  InputReader input_reader((root / "values.txt").string());
  std::unordered_map<std::string, double> values;

  
  if (!input_reader.readInputs(values))
  {
    return; 
  }

  
  const int M0 = static_cast<int>(values["M"]);
  const int N0 = static_cast<int>(values["N"]);

  const long long S = M0 * N0 * 1LL;
  double t = values["t"];
  double l = values["l"];

  
  const long double ratios[4] = {1, 2, 3, 4};

  std::vector<Row> rows_point;
  rows_point.reserve(4);
  std::vector<Row> rows_buffon;
  rows_buffon.reserve(4);

  const fs::path outDir = Paths::estimatePiDataDir();
  ensure_dir(outDir);

  for (auto r : ratios)
  {
    const long long M_target = std::llround(std::sqrt(S / r));
    const int M = nearest_divisor(S, M_target);
    const int N = S / M;
    
    Statistics sp = random_point_method(M, N);
    Statistics sb = buffons_method(M, N, t, l);

    
    rows_point.push_back(sp.to_row(M, N, M_PI));
    rows_buffon.push_back(sb.to_row(M, N, M_PI));

    sb.write_to_bin_file(
        (outDir / ("buffons_" + std::to_string(M) + "x" + std::to_string(N) + ".bin")).string());

    sp.write_to_bin_file(
        (outDir / ("point_method_" + std::to_string(M) + "x" + std::to_string(N) + ".bin")).string());

    std::cout << "Done. For M=" << M << "and N=" << N << "\n";
  }
  write_rows_to_csv(rows_point, (outDir / "point_method_ratios.csv").string());
  write_rows_to_csv(rows_buffon, (outDir / "buffons_ratios.csv").string());

  std::cout << "Done. S=" << S << " held constant across four N/M ratios.\n";

  return;
}
