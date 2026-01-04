#ifndef SURFACE_HPP
#define SURFACE_HPP
#include <string>
#include <memory>
#include "Vec3.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

class Surface
{
public:
  explicit Surface(std::string name_) : name(std::move(name_)) {}
  virtual ~Surface() = default;

  
  virtual double evaluate(const Vec3 &p) const = 0;
  
  virtual double distance(const Ray &ray) const = 0;

  std::string name;
};

using SurfacePtr = std::shared_ptr<Surface>;

class PlaneX : public Surface
{
public:
  PlaneX(std::string n, double x0) : Surface(std::move(n)), x(x0) {}
  double evaluate(const Vec3 &p) const override { return p.x - x; }
  double distance(const Ray &r) const override
  {
    if (std::abs(r.dir.x) < 1e-12)
      return std::numeric_limits<double>::infinity();
    double t = (x - r.origin.x) / r.dir.x;
    return t > 1e-12 ? t : std::numeric_limits<double>::infinity();
  }
  double x;
};

class PlaneY : public Surface
{
public:
  PlaneY(std::string n, double y0) : Surface(std::move(n)), y(y0) {}
  double evaluate(const Vec3 &p) const override { return p.y - y; }
  double distance(const Ray &r) const override
  {
    if (std::abs(r.dir.y) < 1e-12)
      return std::numeric_limits<double>::infinity();
    double t = (y - r.origin.y) / r.dir.y;
    return t > 1e-12 ? t : std::numeric_limits<double>::infinity();
  }
  double y;
};

class PlaneZ : public Surface
{
public:
  PlaneZ(std::string n, double z0) : Surface(std::move(n)), z(z0) {}
  double evaluate(const Vec3 &p) const override { return p.z - z; }
  double distance(const Ray &r) const override
  {
    if (std::abs(r.dir.z) < 1e-12)
      return std::numeric_limits<double>::infinity();
    double t = (z - r.origin.z) / r.dir.z;
    return t > 1e-12 ? t : std::numeric_limits<double>::infinity();
  }
  double z;
};

class Sphere : public Surface
{
public:
  Sphere(std::string n, Vec3 c, double r) : Surface(std::move(n)), center(c), radius(r) {}
  double evaluate(const Vec3 &p) const override
  {
    Vec3 d = p - center;
    return d.dot(d) - radius * radius;
  }
  double distance(const Ray &ray) const override
  {
    Vec3 oc = ray.origin - center;
    double b = oc.dot(ray.dir);
    double c = oc.dot(oc) - radius * radius;
    double disc = b * b - c;
    if (disc < 0)
      return std::numeric_limits<double>::infinity();
    double t = -b - std::sqrt(disc);
    if (t > 1e-8)
      return t;
    t = -b + std::sqrt(disc);
    return t > 1e-8 ? t : std::numeric_limits<double>::infinity();
  }
  Vec3 center;
  double radius;
};

class CylinderZ : public Surface
{
public:
  CylinderZ(std::string n, Vec3 c, double r) : Surface(std::move(n)), center(c), radius(r) {}
  double evaluate(const Vec3 &p) const override
  {
    double dx = p.x - center.x;
    double dy = p.y - center.y;
    return dx * dx + dy * dy - radius * radius;
  }
  double distance(const Ray &ray) const override
  {
    Vec3 d = ray.dir;
    Vec3 oc = ray.origin - center;
    double a = d.x * d.x + d.y * d.y;
    double b = 2.0 * (oc.x * d.x + oc.y * d.y);
    double c = oc.x * oc.x + oc.y * oc.y - radius * radius;
    double disc = b * b - 4 * a * c;
    if (disc < 0 || a == 0.0)
      return std::numeric_limits<double>::infinity();
    double t0 = (-b - std::sqrt(disc)) / (2 * a);
    double t1 = (-b + std::sqrt(disc)) / (2 * a);
    double t = std::numeric_limits<double>::infinity();
    if (t0 > 1e-8)
      t = t0;
    else if (t1 > 1e-8)
      t = t1;
    return t;
  }
  Vec3 center;
  double radius;
};

class Cuboid : public Surface
{
public:
  Cuboid(std::string n, Vec3 min_, Vec3 max_) : Surface(std::move(n)), min(min_), max(max_) {}
  double evaluate(const Vec3 &p) const override
  {
    
    
    
    double ox = std::max(std::max(min.x - p.x, 0.0), p.x - max.x);
    double oy = std::max(std::max(min.y - p.y, 0.0), p.y - max.y);
    double oz = std::max(std::max(min.z - p.z, 0.0), p.z - max.z);
    double outside = ox + oy + oz;
    if (outside > 0.0)
      return outside;

    double dx_in = std::min(p.x - min.x, max.x - p.x);
    double dy_in = std::min(p.y - min.y, max.y - p.y);
    double dz_in = std::min(p.z - min.z, max.z - p.z);
    double dmin = std::min(dx_in, std::min(dy_in, dz_in));
    return -dmin;
  }
  double distance(const Ray &r) const override
  {
    
    const double INF = std::numeric_limits<double>::infinity();
    double t0 = -INF;
    double t1 = INF;

    auto update_interval = [&](double origin, double dir, double mn, double mx)
    {
      if (std::abs(dir) < 1e-14)
      {
        
        if (origin < mn || origin > mx)
          return false;
        return true;
      }
      double ta = (mn - origin) / dir;
      double tb = (mx - origin) / dir;
      if (ta > tb)
        std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      return t0 <= t1;
    };

    if (!update_interval(r.origin.x, r.dir.x, min.x, max.x))
      return INF;
    if (!update_interval(r.origin.y, r.dir.y, min.y, max.y))
      return INF;
    if (!update_interval(r.origin.z, r.dir.z, min.z, max.z))
      return INF;

    
    if (t0 > 1e-8)
      return t0;
    if (t1 > 1e-8)
      return t1;
    return INF;
  }
  Vec3 min, max;
};



class SquarePrismZ : public Surface
{
public:
  SquarePrismZ(std::string n, Vec3 c, double half_width_)
      : Surface(std::move(n)), center(c), half_width(half_width_) {}

  double evaluate(const Vec3 &p) const override
  {
    double dx = std::abs(p.x - center.x);
    double dy = std::abs(p.y - center.y);
    return std::max(dx, dy) - half_width; 
  }

  double distance(const Ray &r) const override
  {
    
    const double INF = std::numeric_limits<double>::infinity();
    double best = INF;
    auto try_plane = [&](double t)
    {
      if (t <= 1e-8)
        return;
      Vec3 p = r.origin + r.dir * t;
      if (evaluate(p) <= 1e-10)
        best = std::min(best, t);
    };

    if (std::abs(r.dir.x) > 1e-14)
    {
      try_plane((center.x - half_width - r.origin.x) / r.dir.x);
      try_plane((center.x + half_width - r.origin.x) / r.dir.x);
    }
    if (std::abs(r.dir.y) > 1e-14)
    {
      try_plane((center.y - half_width - r.origin.y) / r.dir.y);
      try_plane((center.y + half_width - r.origin.y) / r.dir.y);
    }
    return best;
  }

  Vec3 center;
  double half_width;
};



class TruncatedCylinderZ : public Surface
{
public:
  TruncatedCylinderZ(std::string n, Vec3 c, double r, double zmin_, double zmax_)
      : Surface(std::move(n)), center(c), radius(r), zmin(zmin_), zmax(zmax_) {}

  double evaluate(const Vec3 &p) const override
  {
    double dx = p.x - center.x;
    double dy = p.y - center.y;
    double radial = dx * dx + dy * dy - radius * radius;
    double below = zmin - p.z;
    double above = p.z - zmax;
    
    return std::max(radial, std::max(below, above));
  }

  double distance(const Ray &ray) const override
  {
    
    const double INF = std::numeric_limits<double>::infinity();
    double best = INF;

    auto accept = [&](double t)
    {
      if (t <= 1e-8)
        return;
      Vec3 p = ray.origin + ray.dir * t;
      if (evaluate(p) <= 1e-10)
        best = std::min(best, t);
    };

    
    Vec3 d = ray.dir;
    Vec3 oc = ray.origin - center;
    double a = d.x * d.x + d.y * d.y;
    double b = 2.0 * (oc.x * d.x + oc.y * d.y);
    double c = oc.x * oc.x + oc.y * oc.y - radius * radius;
    double disc = b * b - 4 * a * c;
    if (disc >= 0.0 && a > 0.0)
    {
      double sdisc = std::sqrt(disc);
      double t0 = (-b - sdisc) / (2 * a);
      double t1 = (-b + sdisc) / (2 * a);
      auto check_side = [&](double t)
      {
        if (t <= 1e-8)
          return;
        Vec3 p = ray.origin + ray.dir * t;
        if (p.z >= zmin - 1e-10 && p.z <= zmax + 1e-10)
          accept(t);
      };
      check_side(t0);
      check_side(t1);
    }

    
    if (std::abs(ray.dir.z) > 1e-14)
    {
      auto check_cap = [&](double zcap)
      {
        double t = (zcap - ray.origin.z) / ray.dir.z;
        if (t <= 1e-8)
          return;
        Vec3 p = ray.origin + ray.dir * t;
        double dx = p.x - center.x;
        double dy = p.y - center.y;
        if (dx * dx + dy * dy <= radius * radius + 1e-10)
          accept(t);
      };
      check_cap(zmin);
      check_cap(zmax);
    }
    return best;
  }

  Vec3 center;
  double radius;
  double zmin, zmax;
};

class GeneralPlane : public Surface
{
public:
  GeneralPlane(std::string n, double a_, double b_, double c_, double d_) : Surface(std::move(n)), a(a_), b(b_), c(c_), d(d_) {}
  double evaluate(const Vec3 &p) const override { return a * p.x + b * p.y + c * p.z + d; }
  double distance(const Ray &r) const override
  {
    double denom = a * r.dir.x + b * r.dir.y + c * r.dir.z;
    if (std::abs(denom) < 1e-12)
      return std::numeric_limits<double>::infinity();
    double t = -(a * r.origin.x + b * r.origin.y + c * r.origin.z + d) / denom;
    return t > 1e-8 ? t : std::numeric_limits<double>::infinity();
  }
  double a, b, c, d;
};

class HexPrism : public Surface
{
public:
  HexPrism(std::string n, double r_flat, double zmin_, double zmax_) : Surface(std::move(n)), R(r_flat), zmin(zmin_), zmax(zmax_) {}
  double evaluate(const Vec3 &p) const override
  {
    
    
    const double s3 = 0.8660254037844386; 
    const double h = 0.5;

    double max_violation = -std::numeric_limits<double>::infinity();
    auto update = [&](double nx, double ny)
    {
      max_violation = std::max(max_violation, nx * p.x + ny * p.y - R);
    };

    update(1.0, 0.0);
    update(-1.0, 0.0);
    update(h, s3);
    update(-h, s3);
    update(h, -s3);
    update(-h, -s3);

    double dz = std::max(std::max(zmin - p.z, 0.0), p.z - zmax);
    return std::max(max_violation, dz);
  }
  double distance(const Ray &ray) const override
  {
    
    
    const double INF = std::numeric_limits<double>::infinity();
    double t_enter = 0.0;
    double t_exit = INF;

    auto update_halfspace = [&](double nx, double ny, double nz, double rhs) -> bool
    {
      double num = rhs - (nx * ray.origin.x + ny * ray.origin.y + nz * ray.origin.z);
      double den = nx * ray.dir.x + ny * ray.dir.y + nz * ray.dir.z;

      if (std::abs(den) < 1e-14)
      {
        
        return num >= 0.0;
      }

      double t = num / den;
      if (den > 0.0)
      {
        
        t_exit = std::min(t_exit, t);
      }
      else
      {
        
        t_enter = std::max(t_enter, t);
      }
      return t_enter <= t_exit;
    };

    
    const double s3 = 0.8660254037844386; 
    const double h = 0.5;

    if (!update_halfspace(1.0, 0.0, 0.0, R))
      return INF;
    if (!update_halfspace(-1.0, 0.0, 0.0, R))
      return INF;

    if (!update_halfspace(h, s3, 0.0, R))
      return INF;
    if (!update_halfspace(-h, s3, 0.0, R))
      return INF;
    if (!update_halfspace(h, -s3, 0.0, R))
      return INF;
    if (!update_halfspace(-h, -s3, 0.0, R))
      return INF;

    
    if (!update_halfspace(0.0, 0.0, 1.0, zmax))
      return INF;
    if (!update_halfspace(0.0, 0.0, -1.0, -zmin))
      return INF;

    
    bool inside0 = evaluate(ray.origin) <= 0.0;
    double t = inside0 ? t_exit : t_enter;
    if (t > 1e-8)
      return t;
    
    double t2 = inside0 ? t_enter : t_exit;
    if (t2 > 1e-8)
      return t2;
    return INF;
  }
  double R, zmin, zmax;
};

class Torus : public Surface
{
public:
  Torus(std::string n, double R_, double r_) : Surface(std::move(n)), R(R_), r(r_) {}
  double evaluate(const Vec3 &p) const override
  {
    double q = std::sqrt(p.x * p.x + p.y * p.y) - R;
    return q * q + p.z * p.z - r * r;
  }
  double distance(const Ray &ray) const override
  {
    
    
    
    const double INF = std::numeric_limits<double>::infinity();
    auto f = [&](double t)
    {
      Vec3 p = ray.origin + ray.dir * t;
      return evaluate(p);
    };

    double t0 = 0.0;
    double f0 = f(t0);
    
    if (std::abs(f0) < 1e-12)
    {
      t0 = 1e-6;
      f0 = f(t0);
    }

    
    double t1 = 1e-3;
    double f1 = f(t1);
    int it = 0;
    while (f0 * f1 > 0.0 && it < 60 && t1 < 1e6)
    {
      t1 *= 2.0;
      f1 = f(t1);
      it++;
    }
    if (f0 * f1 > 0.0)
      return INF;

    
    for (int k = 0; k < 80; ++k)
    {
      double tm = 0.5 * (t0 + t1);
      double fm = f(tm);
      if (f0 * fm <= 0.0)
      {
        t1 = tm;
        f1 = fm;
      }
      else
      {
        t0 = tm;
        f0 = fm;
      }
    }
    double t = 0.5 * (t0 + t1);
    return (t > 1e-8) ? t : INF;
  }
  double R, r;
};



class EllipticalTorus : public Surface
{
public:
  EllipticalTorus(std::string n, double R_, double a_, double b_)
      : Surface(std::move(n)), R(R_), a(a_), b(b_) {}

  double evaluate(const Vec3 &p) const override
  {
    double rho = std::sqrt(p.x * p.x + p.y * p.y);
    double q = rho - R;
    return (q * q) / (a * a) + (p.z * p.z) / (b * b) - 1.0;
  }

  double distance(const Ray &ray) const override
  {
    
    const double INF = std::numeric_limits<double>::infinity();
    auto f = [&](double t)
    {
      Vec3 p = ray.origin + ray.dir * t;
      return evaluate(p);
    };

    double t0 = 0.0;
    double f0 = f(t0);
    if (std::abs(f0) < 1e-12)
    {
      t0 = 1e-6;
      f0 = f(t0);
    }

    double t1 = 1e-3;
    double f1 = f(t1);
    int it = 0;
    while (f0 * f1 > 0.0 && it < 70 && t1 < 1e6)
    {
      t1 *= 2.0;
      f1 = f(t1);
      it++;
    }
    if (f0 * f1 > 0.0)
      return INF;

    for (int k = 0; k < 90; ++k)
    {
      double tm = 0.5 * (t0 + t1);
      double fm = f(tm);
      if (f0 * fm <= 0.0)
      {
        t1 = tm;
        f1 = fm;
      }
      else
      {
        t0 = tm;
        f0 = fm;
      }
    }
    double t = 0.5 * (t0 + t1);
    return (t > 1e-8) ? t : INF;
  }

  double R, a, b;
};

struct Mat3
{
  std::array<double, 9> m{};
  Vec3 apply(const Vec3 &v) const
  {
    return {
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z};
  }
  Mat3 transpose() const
  {
    Mat3 t;
    t.m = {m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]};
    return t;
  }
};

inline Mat3 rotation_matrix(double rx, double ry, double rz)
{
  double cx = std::cos(rx), sx = std::sin(rx);
  double cy = std::cos(ry), sy = std::sin(ry);
  double cz = std::cos(rz), sz = std::sin(rz);
  Mat3 Rx{{1, 0, 0, 0, cx, -sx, 0, sx, cx}};
  Mat3 Ry{{cy, 0, sy, 0, 1, 0, -sy, 0, cy}};
  Mat3 Rz{{cz, -sz, 0, sz, cz, 0, 0, 0, 1}};
  auto mul = [](const Mat3 &A, const Mat3 &B)
  {
    Mat3 R;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        R.m[i * 3 + j] = 0.0;
        for (int k = 0; k < 3; ++k)
          R.m[i * 3 + j] += A.m[i * 3 + k] * B.m[k * 3 + j];
      }
    return R;
  };
  return mul(mul(Rz, Ry), Rx);
}

class TransformSurface : public Surface
{
public:
  TransformSurface(std::string n, SurfacePtr base_, Vec3 translation_, double rx_, double ry_, double rz_)
      : Surface(std::move(n)), base(std::move(base_)), translation(translation_), R(rotation_matrix(rx_, ry_, rz_)), Rinv(R.transpose()) {}

  double evaluate(const Vec3 &p) const override
  {
    Vec3 local = to_local(p);
    return base->evaluate(local);
  }

  double distance(const Ray &ray) const override
  {
    Ray local;
    local.origin = to_local(ray.origin);
    local.dir = Rinv.apply(ray.dir);
    double t = base->distance(local);
    return t;
  }

private:
  Vec3 to_local(const Vec3 &p) const
  {
    Vec3 shifted = p - translation;
    return Rinv.apply(shifted);
  }

  SurfacePtr base;
  Vec3 translation;
  Mat3 R;
  Mat3 Rinv;
};

#endif
