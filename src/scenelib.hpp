#pragma once
#include "types.hpp"
#include <functional>
#include <map>
#include <stdexcept>

template <typename FT> class Scenes {
  using T = Types<FT>;
  const static std::map<std::string, typename T::SurfaceFunction> surfaces;

public:
  static T::SurfaceFunction GetSceneByName(std::string name) {
    auto it = surfaces.find(name);
    if (it == surfaces.end())
      throw std::domain_error("The requested surface is not found");
    return it->second;
  }

private:
  // Quadratic surfaces:
  constexpr static T::SurfaceFunction sphere() {
    return [](FT x, FT y, FT z) {
      FT R = 0.5;
      return x * x + y * y + z * z - R * R;
    };
  }

  // Cubic surfaces:
  // source: https://wims.univ-cotedazur.fr/gallery/
  constexpr static T::SurfaceFunction cubic_sheet() {
    return [](FT x, FT y, FT z) {
      return x * x * y + y * y * z + z * z * x - static_cast<FT>(0.1);
    };
  }

  // source: https://wims.univ-cotedazur.fr/gallery/
  constexpr static T::SurfaceFunction cubic_holes() {
    return [](FT x, FT y, FT z) {
      return x * x * x - x + y * y * y - y + z * z * z - z;
    };
  }

  // Quartic surfaces:
  // source: https://wims.univ-cotedazur.fr/gallery/
  constexpr static T::SurfaceFunction riemann() {
    return [](FT x, FT y, FT z) {
      return z * z * x * x + (z * z + 1) * y * y - 5 * (z * z * z * z + z * z);
    };
  }

  // source: https://wims.univ-cotedazur.fr/gallery/
  constexpr static T::SurfaceFunction kummer() {
    return [](FT x, FT y, FT z) {
      return glm::pow(x, 4) + glm::pow(y, 4) + glm::pow(z, 4) -
             glm::pow(0.5f * x + 1, 2) * (x * x + y * y + z * z) -
             (x * x * y * y + x * x * z * z + y * y * z * z) +
             glm::pow(0.5f * x + 1, 4);
    };
  }

  constexpr static T::SurfaceFunction torus_deg_4() {
    return [](FT x, FT y, FT z) {
      FT R = 2.0;
      FT r = 0.5;
      FT lhs = 4 * R * R * (x * x + y * y);
      FT rhs = (x * x + y * y + z * z + R * R - r * r);
      return rhs * rhs - lhs;
    };
  }

  // Sextic surfaces:
  // source: https://wims.univ-cotedazur.fr/gallery/
  constexpr static T::SurfaceFunction barth_sextic() {
    return [](FT x, FT y, FT z) {
      return 4 * (2.618f * x * x - y * x) * (2.618f * y * y - z * z) *
                 (2.618f * z * z - x * x) -
             4.236f * glm::pow(x * x + y * y + z * z - 1, 2);
    };
  }

  // SDFs:
  // source: https://iquilezles.org/articles/distfunctions/
  constexpr static T::SurfaceFunction sphere_sdf() {
    return [](FT x, FT y, FT z) {
      FT R = 0.5;
      return glm::sqrt(x * x + y * y + z * z) - R;
    };
  }

  constexpr static T::SurfaceFunction box() {
    return [](FT x, FT y, FT z) {
      typename T::gv3 b = typename T::gv3(0.5, 1, 1.5); // half bounds
      typename T::gv3 p = typename T::gv3(x, y, z);
      typename T::gv3 q = glm::abs(p) - b;
      return glm::length(glm::max(q, static_cast<FT>(0))) +
             glm::min(glm::max(q.x, glm::max(q.y, q.z)), static_cast<FT>(0));
    };
  }

  constexpr static T::SurfaceFunction torus_sdf() {
    return [](FT x, FT y, FT z) {
      typename T::gv3 p = typename T::gv3(x, y, z);
      typename T::gv2 t = typename T::gv2(2.0, 0.5);
      typename T::gv2 q(glm::length(typename T::gv2(p.x,p.y))-t.x,p.y);
      return glm::length(q)-t.y;
    };
  }
};

extern template class Scenes<float>;
extern template class Scenes<double>;
