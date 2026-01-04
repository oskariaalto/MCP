#ifndef MATERIALS_HPP
#define MATERIALS_HPP
#include "../XS/CrossSection.hpp"
#include <unordered_map>
#include <string>

struct Material
{
  std::string name;
  Mixture mix;
};

class MaterialLibrary
{
public:
  void add(const Material &m) { mats_[m.name] = m; }
  const Material *get(const std::string &name) const
  {
    auto it = mats_.find(name);
    if (it == mats_.end())
      return nullptr;
    return &it->second;
  }

  std::unordered_map<std::string, Material> &all() { return mats_; }
  const std::unordered_map<std::string, Material> &all() const { return mats_; }

private:
  std::unordered_map<std::string, Material> mats_;
};

#endif
