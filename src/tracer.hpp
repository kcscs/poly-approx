#pragma once

#include "exceptions.hpp"
#include "types.hpp"
#include "chebseg.hpp"
#include "monseg.hpp"
#include "fun.hpp"
#include "glm_parsers.hpp"

template <typename FT> class TraceMethod {

public:
  using T = Types<FT>;
  struct TraceResult {
    bool hit = false;
    FT distance;
    json metadata;
  };

  virtual TraceResult trace(T::Ray ray, T::SurfaceFunction f, json params) const = 0;
};

/// Interpolates a chebyshev polynomial over the ray which is then converted to monomial basis.
/// Then Cem Yuksel's rootfinding is applied.
template <typename FT> class PolynomialTracer : public TraceMethod<FT> {
public:
  using typename TraceMethod<FT>::TraceResult;
  using typename TraceMethod<FT>::T;
  TraceResult trace(T::Ray ray, T::SurfaceFunction f,
                    json params) const override {
                      glm::vec2 clip_distances = params.value("clip", glm::vec2(0.0f, 20.0f));

                      typename T::RRFunction func = [&](FT t) {
                        typename T::gv3 p = ray.start + t * ray.dir;
                        return f(p.x, p.y, p.z);
                      };

                      SegFunApproximator<FT, ChebFun<FT>, ChebSeg<FT>> approximator;
                      ChebFun<FT> chebfun = approximator(func, clip_distances.x, clip_distances.y);


                    }
};
