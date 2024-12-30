#pragma once

#include "chebseg.hpp"
#include "exceptions.hpp"
#include "fun.hpp"
#include "glm_parsers.hpp"
#include "logging.hpp"
#include "monseg.hpp"
#include "types.hpp"
#include "split_strategy.hpp"
#include "stop_strategy.hpp"
#include <limits>
#include <memory>

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
    ChebFun<FT> chebfun = approximator(func, clip_distances.x, clip_distances.y,
                                       cheb_approx_metadata);
    res.metadata["chebyshev_interpolation"] = cheb_approx_metadata;

    json monomial_conversion_metadata;
    SegFun<FT, MonSeg<FT>> monfun =
        chebfun.template Convert<MonSeg<FT>>(monomial_conversion_metadata);
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

protected:
  T::gv2 clip_distances;
  json settings;
};

template <typename FT>
class FirstRootPolynomialTracer : public PolynomialTracer<FT> {
  using typename TraceMethod<FT>::TraceResult;
  using typename TraceMethod<FT>::T;

public:
  FirstRootPolynomialTracer(json settings) : PolynomialTracer<FT>(settings) {
    max_degree = settings["max_degree"];
    target_precision = settings["target_precision"];
    split_strategy = CreateSplitStrategy<FT>(settings);
    stop_strategy = CreateStopStrategy<FT>(settings);
  }

  TraceResult trace(T::Ray ray, T::SurfaceFunction f) const override {

    Logger &log = Logger::Get();
    TraceResult res;
    typename T::RRFunction func = [&](FT t) {
      typename T::gv3 p = ray.start + t * ray.dir;
      return f(p.x, p.y, p.z);
    };

    FT cur_begin = this->clip_distances.x;
    FT cur_end = this->clip_distances.y;

    std::vector<ChebSeg<FT>> computed_cheb_segments;
    std::vector<MonSeg<FT>> computed_mon_segments;

    stop_strategy->reset();

    while (cur_begin < this->clip_distances.y) {
      ChebSeg<FT> seg =
          ChebSeg<FT>::Interpolate(func, max_degree, cur_begin, cur_end);

      auto err = seg.Error(func);
      FT err_val = std::get<0>(err);
      FT max_err_place = std::get<3>(err);
      seg.metadata["interstitial_error"] = err_val;
      seg.metadata["interstitial_error_max_loc"] = max_err_place;

      if (!stop_strategy->stop(err_val)) {
        SplitContext<FT> ctx;
        ctx.max_error_place = max_err_place;
        log << "split"_cat << "split: " << cur_begin << "-" << cur_end
            << " e:" << err_val << "\n";
        cur_end = split_strategy->split(cur_begin, cur_end, ctx);

      } else {
        log << "segment"_cat << "segment: " << cur_begin << "-" << cur_end
            << " err: " << err_val << "\n";

        computed_cheb_segments.push_back(seg);
        MonSeg<FT> monseg = MonSeg<FT>::FitAtChebPoints(seg);
        computed_mon_segments.push_back(monseg);
        json rootdata;
        std::vector<FT> roots = monseg.FindRoots(rootdata);
        if (roots.size() > 0) {
          res.hit = true;
          res.distance = roots[0];
          res.metadata["rootfinding"] = rootdata;
          res.metadata["chebyshev_segments"] = computed_cheb_segments;
          res.metadata["power_segments"] = computed_mon_segments;
          return res;
        } else {
          cur_begin = cur_end;
          cur_end = this->clip_distances.y;
          stop_strategy->reset();
        }
      }
    }
    res.hit = false;
    res.metadata["chebyshev_segments"] = computed_cheb_segments;
    res.metadata["power_segments"] = computed_mon_segments;
    return res;
  }

protected:
  int max_degree;
  FT target_precision;

  std::unique_ptr<SplitStrategy<FT>> split_strategy;
  std::unique_ptr<StopStrategy<FT>> stop_strategy;
};
