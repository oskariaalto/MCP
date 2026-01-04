#ifndef CSG_HPP
#define CSG_HPP
#include <vector>
#include <memory>
#include <functional>
#include "Surface.hpp"
#include <cmath>

enum class CsgOp
{
  SURFACE,
  NOT,
  AND,
  OR
};

struct CsgNode
{
  CsgOp op = CsgOp::SURFACE;
  SurfacePtr surface;
  std::shared_ptr<CsgNode> left;
  std::shared_ptr<CsgNode> right;
};

inline std::shared_ptr<CsgNode> surf_node(const SurfacePtr &s)
{
  auto n = std::make_shared<CsgNode>();
  n->op = CsgOp::SURFACE;
  n->surface = s;
  return n;
}

inline std::shared_ptr<CsgNode> not_node(const std::shared_ptr<CsgNode> &a)
{
  auto n = std::make_shared<CsgNode>();
  n->op = CsgOp::NOT;
  n->left = a;
  return n;
}

inline std::shared_ptr<CsgNode> and_node(const std::shared_ptr<CsgNode> &a, const std::shared_ptr<CsgNode> &b)
{
  auto n = std::make_shared<CsgNode>();
  n->op = CsgOp::AND;
  n->left = a;
  n->right = b;
  return n;
}

inline std::shared_ptr<CsgNode> or_node(const std::shared_ptr<CsgNode> &a, const std::shared_ptr<CsgNode> &b)
{
  auto n = std::make_shared<CsgNode>();
  n->op = CsgOp::OR;
  n->left = a;
  n->right = b;
  return n;
}

inline bool eval_node(const CsgNode &n, const Vec3 &p)
{
  switch (n.op)
  {
  case CsgOp::SURFACE:
    return n.surface->evaluate(p) <= 0.0;
  case CsgOp::NOT:
    return !eval_node(*n.left, p);
  case CsgOp::AND:
    return eval_node(*n.left, p) && eval_node(*n.right, p);
  case CsgOp::OR:
    return eval_node(*n.left, p) || eval_node(*n.right, p);
  }
  return false;
}

struct Cell
{
  std::string name;
  std::shared_ptr<CsgNode> root;
};

struct Universe
{
  std::vector<Cell> cells;

  const Cell *find_cell(const Vec3 &p) const
  {
    for (const auto &c : cells)
    {
      if (eval_node(*c.root, p))
        return &c;
    }
    return nullptr;
  }
};

struct Lattice
{
  Vec3 pitch;      
  Vec3 origin;     
  int nx = 1, ny = 1;
  const Universe *universe = nullptr;

  const Cell *find_cell(const Vec3 &p) const
  {
    if (pitch.x == 0.0 || pitch.y == 0.0)
      return nullptr;
    double lx = (p.x - origin.x) / pitch.x;
    double ly = (p.y - origin.y) / pitch.y;
    int ix = static_cast<int>(std::floor(lx));
    int iy = static_cast<int>(std::floor(ly));
    if (ix < 0 || iy < 0 || ix >= nx || iy >= ny)
      return nullptr;
    
    Vec3 local{
        p.x - (origin.x + (ix + 0.5) * pitch.x),
        p.y - (origin.y + (iy + 0.5) * pitch.y),
        p.z};
    return universe ? universe->find_cell(local) : nullptr;
  }
};

struct HexLattice
{
  double pitch = 1.0; 
  Vec3 origin;        
  int radius = 1;     
  const Universe *universe = nullptr;

  struct Axial
  {
    int q;
    int r;
  };

  static Axial round_axial(double qf, double rf)
  {
    double xf = qf;
    double zf = rf;
    double yf = -xf - zf;
    int xr = (int)std::round(xf);
    int yr = (int)std::round(yf);
    int zr = (int)std::round(zf);

    double dx = std::abs(xr - xf);
    double dy = std::abs(yr - yf);
    double dz = std::abs(zr - zf);

    if (dx > dy && dx > dz)
      xr = -yr - zr;
    else if (dy > dz)
      yr = -xr - zr;
    else
      zr = -xr - yr;
    return {xr, zr};
  }

  const Cell *find_cell(const Vec3 &p) const
  {
    Vec3 rel = p - origin;
    double qf = (std::sqrt(3.0) / 3.0 * rel.x - 1.0 / 3.0 * rel.y) / pitch;
    double rf = (2.0 / 3.0 * rel.y) / pitch;
    Axial ax = round_axial(qf, rf);
    int s = -ax.q - ax.r;
    int dist = std::max({std::abs(ax.q), std::abs(ax.r), std::abs(s)});
    if (dist > radius)
      return nullptr;

    double cx = pitch * (std::sqrt(3.0) * ax.q + std::sqrt(3.0) / 2.0 * ax.r);
    double cy = pitch * (1.5 * ax.r);
    Vec3 local{rel.x - cx, rel.y - cy, rel.z};
    return universe ? universe->find_cell(local) : nullptr;
  }
};

#endif
