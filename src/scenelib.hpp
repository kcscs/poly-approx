#pragma once
#include "types.hpp"
#include <functional>
#include <map>
#include <stdexcept>

template<typename FT>
class Scenes {
  using T = Types<FT>;
  const static std::map<std::string, typename T::SurfaceFunction> surfaces;
  
public:
    static T::SurfaceFunction GetSceneByName(std::string name) {
      auto it = surfaces.find(name);
      if(it == surfaces.end())
        throw std::domain_error("The requested surface is not found");
      return it->second;
    }

  private:


    constexpr static T::SurfaceFunction torus_deg_4() {
      return [](FT x, FT y, FT z) {
        FT R = 2.0;
        FT r = 0.5;
        FT lhs = 4*R*R*(x*x+y*y);
        FT rhs = (x*x+y*y+z*z+R*R-r*r);
        return rhs*rhs-lhs;
      };
    }

    constexpr static T::SurfaceFunction sphere() {
    return [](FT x, FT y, FT z) {
      FT R = 0.5;
      return x*x+y*y+z*z-R*R;
    };
  }

};
