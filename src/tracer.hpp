#pragma once

#include "chebseg.hpp"
#include "exceptions.hpp"
#include "fun.hpp"
#include "glm_parsers.hpp"
#include "monseg.hpp"
#include "types.hpp"
#include <limits>

template <typename FT> class TraceMethod {

public:
  using T = Types<FT>;
  struct TraceResult {
    bool hit = false;
    FT distance;
    json metadata;
  };

  virtual TraceResult trace(T::Ray ray, T::SurfaceFunction f) const = 0;
};

/// Interpolates a chebyshev polynomial over the ray which is then converted to
/// monomial basis. Then Cem Yuksel's rootfinding is applied.
template <typename FT> class PolynomialTracer : public TraceMethod<FT> {
public:
  using typename TraceMethod<FT>::TraceResult;
  using typename TraceMethod<FT>::T;
  PolynomialTracer(json settings) : TraceMethod<FT>() {
    clip_distances = settings["clip"];
    this->settings = settings;
  }

  TraceResult trace(T::Ray ray, T::SurfaceFunction f) const override {

    TraceResult res;
    typename T::RRFunction func = [&](FT t) {
      typename T::gv3 p = ray.start + t * ray.dir;
      return f(p.x, p.y, p.z);
    };

    json cheb_approx_metadata;
    SegFunApproximator<FT, ChebFun<FT>, ChebSeg<FT>> approximator(settings);
    ChebFun<FT> chebfun =
        approximator(func, clip_distances.x, clip_distances.y, cheb_approx_metadata);
    res.metadata["chebyshev_interpolation"] = cheb_approx_metadata;

    json monomial_conversion_metadata;
    SegFun<FT, MonSeg<FT>> monfun = chebfun.template Convert<MonSeg<FT>>(monomial_conversion_metadata);
    res.metadata["conversion_to_monomial"] = monomial_conversion_metadata;

    json root_finding_metadata;
    auto roots = monfun.FindRoots(root_finding_metadata);
    res.metadata["rootfinding"] = root_finding_metadata;

    if (!roots.empty()) {
      res.hit = true;
      res.distance = roots[0];
    } else {
      res.hit = false;
      res.distance = std::numeric_limits<FT>::quiet_NaN();
    }
    return res;
  }

private:
  T::gv2 clip_distances;
  json settings;
};
