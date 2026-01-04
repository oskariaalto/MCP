#ifndef PATHS_HPP
#define PATHS_HPP

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>

namespace Paths
{
namespace fs = std::filesystem;

inline std::optional<fs::path> findRepoRoot()
{
  fs::path p = fs::current_path();
  for (int i = 0; i < 10; ++i)
  {
    const bool ok = fs::exists(p / "particle_data") && fs::is_directory(p / "particle_data") &&
                    fs::exists(p / "src") && fs::is_directory(p / "src") &&
                    (fs::exists(p / "Makefile") || fs::exists(p / "CMakeLists.txt"));
    if (ok)
      return p;

    if (!p.has_parent_path())
      break;
    const fs::path parent = p.parent_path();
    if (parent == p)
      break;
    p = parent;
  }
  return std::nullopt;
}

inline fs::path repoRoot()
{
  static const fs::path root = []() -> fs::path {
    auto found = findRepoRoot();
    if (!found)
    {
      throw std::runtime_error(
          "Could not locate repository root (expected `particle_data/`, `src/`, and a build file like `Makefile`). "
          "Run from the repo root or a build subdirectory.");
    }
    return *found;
  }();
  return root;
}

inline fs::path particleDataDir()
{
  return repoRoot() / "particle_data";
}

inline fs::path resultsDir()
{
  return repoRoot() / "results";
}

inline fs::path resultsDir(const std::string &subdir)
{
  if (subdir.empty())
    return resultsDir();
  return resultsDir() / subdir;
}

inline fs::path estimatePiDataDir()
{
  return repoRoot() / "estimate_pi_data";
}
} 

#endif
