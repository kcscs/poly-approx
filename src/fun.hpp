#pragma once
#include "chebseg.hpp"
#include "exceptions.hpp"
#include "logging.hpp"
#include "monseg.hpp"
#include "seg.hpp"
#include "types.hpp"
#include <Eigen/QR>
#include <concepts>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>

template <typename FT> class Fun : public Types<FT> {
public:
  using typename Types<FT>::RRFunction;
  virtual FT Evaluate(FT x) const = 0;
};

template <typename FT, typename SegT>
  requires std::derived_from<SegT, Seg<FT>>
class SegFun : public Fun<FT> { // TODO: optimize
public:
  using typename Types<FT>::RRFunction;

  SegFun(std::vector<SegT> segments) : segments(segments) {}

  FT Evaluate(FT x) const override {
    if (x < segments[0].begin)
      throw std::out_of_range("evaluated point is before first segment");
    int i = 0;
    while (i < segments.size() && segments[i].end < x)
      ++i;
    if (i == segments.size())
      throw std::out_of_range("evaluated point is after last segment");

    return segments[i].Eval(x);
  }

  template <typename OtherSegT>
    requires std::derived_from<OtherSegT, Seg<FT>>
  SegFun<FT, OtherSegT> Convert() const {
    std::vector<OtherSegT> conv_segments;
    conv_segments.reserve(segments.size());
    for (int i = 0; i < segments.size(); ++i) {
      conv_segments.push_back(OtherSegT::FitAtChebPoints(segments[i]));
      // std::cout<<"asd: "<<segments[i].coeffs.size()<<" -> "<<conv_segments.back().coeffs.size()<<"\n";
    }
    return SegFun<FT, OtherSegT>(conv_segments);
  }

  const std::vector<SegT> &GetSegments() const { return segments; }

  std::vector<FT> FindRoots() const {
    Logger& log = Logger::Get();
    std::vector<FT> roots;
    for (const auto &seg : segments) {
      log<<Logger::cat("rootfinder")<<"coeffs: "<<seg.coeffs.size()<<"\n";
      std::vector<FT> seg_roots = seg.FindRoots();
      roots.insert(roots.end(), seg_roots.begin(), seg_roots.end());
    }
    return roots;
  }

protected:
  std::vector<SegT> segments;
};

template <typename FT, typename FunT, typename SegT>
  requires std::derived_from<FunT, SegFun<FT, SegT>>
class SegFunApproximator : public Types<FT> {
public:
  using T = Types<FT>;
  using typename T::em;
  using typename T::ev;
  using typename T::RRFunction;
  using ET = FT;

  RRFunction f;
  int split_degree;
  FT target_precision;

public:
  SegFunApproximator(json params) {
    split_degree = params["max_degree"];
    target_precision = params["target_precision"];
  }

  virtual FunT operator()(RRFunction f, FT x_begin, FT x_end) {
    Logger &log = Logger::Get();
    std::vector<SegT> segments;
    FT cur_begin = x_begin;
    FT cur_end = x_end;

    ET last_err = std::numeric_limits<ET>::infinity();
    log << Logger::cat("split") << "initial range: " << x_begin << " - "
        << x_end << "\n";
    while (cur_begin + T::eps < x_end) {
      int cur_degree = 1;
      double err = std::numeric_limits<double>::infinity();
      ev coeffs;
      double prev_err = std::numeric_limits<double>::infinity();
      bool first = true;
      bool split = false;
      int best_degree = -1;
      double best_error = std::numeric_limits<double>::infinity();
      FT best_error_place = cur_end;
      while (err > target_precision) {
        prev_err = err;

        ChebSeg<FT> chebseg =
            ChebSeg<FT>::Interpolate(f, cur_degree, cur_begin, cur_end);
        // auto interp_res =
        //     chebyshev_interpolate(f, cur_degree, cur_begin, cur_end);

        auto err_res = chebseg.Error(f);
        log << Logger::cat("cheb_res");
        err = std::get<0>(err_res);
        log << "d:" << cur_degree << "  e:" << err << "\n";

        if (err < best_error) {
          best_error = err;
          best_degree = cur_degree;
          best_error_place = std::get<3>(err_res);
        }

        cur_degree *= 2;
        if (cur_degree > split_degree && err > target_precision) {
          split = true;
          break;
        }
      }

      if (split /*&& cur_end-cur_begin > (x_end-x_begin)/100*/) {
        // cur_end = (cur_begin + cur_end) / 2;
        cur_end = best_error_place;
        log << Logger::cat("split") << "split: " << cur_begin << "-" << cur_end
            << " e:" << err << "\n";
      } else {
        ChebSeg<FT> chebseg =
            ChebSeg<FT>::Interpolate(f, best_degree, cur_begin, cur_end);
        auto err_res = chebseg.Error(f);
        err = std::get<0>(err_res);
        log << Logger::cat("segment") << "segment: " << cur_begin << "-"
            << cur_end << " deg: " << best_degree << " err: " << err << "\n";
        // std::vector<FT> coeff_vec(coeffs.size());
        // for (int i = 0; i < coeffs.size(); i++) {
        //   coeff_vec[i] = coeffs(i);
        // }
        // segments.emplace_back(coeff_vec, cur_begin, cur_end);
        segments.push_back(chebseg);
        cur_begin = cur_end;
        cur_end = x_end;
      }
    }

    return {segments};
  }
};

template <typename FT> class ChebFun : public SegFun<FT, ChebSeg<FT>> {
public:
  using T = Types<FT>;
  using typename Types<FT>::ev;
  using typename Types<FT>::em;

  using typename Types<FT>::RRFunction;

  ChebFun(std::vector<ChebSeg<FT>> _segments)
      : SegFun<FT, ChebSeg<FT>>(_segments) {}

private:
};

template <typename FT> class MonFun : public SegFun<FT, MonSeg<FT>> {
public:
  using T = Types<FT>;
  using typename Types<FT>::ev;
  using typename Types<FT>::em;

  MonFun(std::vector<MonSeg<FT>> _segments)
  : SegFun<FT, MonSeg<FT>>(_segments) {}

private:
};
