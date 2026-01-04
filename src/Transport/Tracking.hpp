#ifndef TRACKING_HPP
#define TRACKING_HPP
#include "../Geometry/CSG.hpp"
#include "../Transport/Materials.hpp"
#include "../Physics/Neutron.hpp"
#include "../Common/Random.hpp"
#include "../Geometry/VolumeEstimator.hpp"
#include <unordered_map>
#include <vector>

struct MaterialCell
{
  const Cell *cell = nullptr;
  const Material *material = nullptr;
};

struct GeometryContext
{
  Universe universe;
  std::unordered_map<std::string, const Material *> cell_material;
  std::vector<SurfacePtr> surfaces; 
  BoundingBox bounds;
};

struct TrackTallies
{
  struct Entry
  {
    double collisions_analog = 0.0;  
    double absorptions_analog = 0.0; 
    double track_length = 0.0;       
    double collisions_tle = 0.0;      
    double absorptions_tle = 0.0;     
    double flux_cfe = 0.0;            
  };
  std::unordered_map<std::string, Entry> by_material;
  double leaks = 0.0;
};

struct MeshTally
{
  int nx = 0, ny = 0, nz = 0;
  BoundingBox bbox;
  std::vector<double> collisions; 

  void init(int nx_, int ny_, int nz_, const BoundingBox &b)
  {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    bbox = b;
    collisions.assign(nx * ny * nz, 0.0);
  }

  void score(const Vec3 &p)
  {
    int ix = static_cast<int>(nx * (p.x - bbox.min.x) / (bbox.max.x - bbox.min.x));
    int iy = static_cast<int>(ny * (p.y - bbox.min.y) / (bbox.max.y - bbox.min.y));
    int iz = static_cast<int>(nz * (p.z - bbox.min.z) / (bbox.max.z - bbox.min.z));
    if (ix >= 0 && iy >= 0 && iz >= 0 && ix < nx && iy < ny && iz < nz)
      collisions[(iz * ny + iy) * nx + ix] += 1.0;
  }
};

double distance_to_exit_cell(const Cell &cell, const Ray &ray);

TrackTallies simulate_external_source(const GeometryContext &geo, const MaterialLibrary &matlib, bool surface_tracking, int histories, MeshTally *mesh);

struct CriticalityResult
{
  double keff = 0.0;
  std::vector<double> gen_keff;
};

CriticalityResult simulate_criticality(const GeometryContext &geo, const MaterialLibrary &matlib, int histories, int generations);

#endif
