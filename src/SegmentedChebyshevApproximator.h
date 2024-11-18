#pragma once

#include "FunctionApproximator.h"
#include "monom.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <cmath>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <matplot/matplot.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

using namespace matplot;

using Eigen::MatrixXd;
using Eigen::Vector2;
using Eigen::VectorXd;

// coeffs[0] is constant, coeffs[i] is deg i
template <typename funcval_T>
funcval_T eval_cheb(const Eigen::Matrix<funcval_T, Eigen::Dynamic, 1> &coeffs,
                    funcval_T x, funcval_T a, funcval_T b) {
  // funcval_T res = 0;
  // for (int i = 0; i < coeffs.size(); ++i) {
  //   res += coeffs(i) * cos(i * acos((2 * x - (b + a)) / (b - a)));
  // }

  x = (2 * x - (b + a)) / (b - a); // transform to [-1;1]

  int deg = coeffs.size() - 1;
  if (deg == 0)
    return coeffs(0);
  if (deg == 1)
    return coeffs(0) + coeffs(1) * x;

  funcval_T Tprev = x;
  funcval_T Tcur = 2 * x * x - 1;
  funcval_T res = coeffs(0) + coeffs(1) * x;
  for (int i = 2; i < deg; ++i) {
    res += coeffs(i) * Tcur;
    funcval_T Ttemp = Tcur;
    Tcur = 2 * x * Tcur - Tprev;
    Tprev = Ttemp;
  }
  res += coeffs(deg) * Tcur;

  return res;
}

template <typename funcval_T> struct ChebyshevSegment {
  using VectorT = Eigen::Matrix<funcval_T, Eigen::Dynamic, 1>;
  using MatrixT = Eigen::Matrix<funcval_T, Eigen::Dynamic, Eigen::Dynamic>;
  VectorT coeffs;
  funcval_T begin, end;

  ChebyshevSegment<funcval_T>() {}
  ChebyshevSegment<funcval_T>(std::initializer_list<funcval_T> coefficients,
                              funcval_T begin, funcval_T end)
      : begin(begin), end(end), coeffs(coefficients.size()) {
    auto iter = coefficients.begin();
    for (int i = 0; i < coefficients.size(); ++i) {
      coeffs[i] = *iter;
      iter++;
    }
  }

  funcval_T Evaluate(funcval_T x) const {
    return eval_cheb(coeffs, x, begin, end);
  }

  void Show(bool open_window = true) const {
    std::vector<double> x = linspace(begin, end);
    std::vector<double> y(x.size());
    for (int i = 0; i < y.size(); ++i) {
      y[i] = Evaluate(x[i]);
    }
    plot(x, y)->line_width(2);
    // if(open_window)
    //     show();
  }

  MonomSegment<funcval_T> ConvertToMonomial() const {
    const int deg = coeffs.size() - 1;
    std::vector<funcval_T> xs(deg + 1);
    VectorT ys(deg + 1);
    for (int i = 0; i <= deg; ++i) {
      // xs[i] = (static_cast<funcval_T>(i) / deg)*2-1; //
      xs[i] = - cos(static_cast<funcval_T>(i) / deg * glm::pi<funcval_T>());
      ys(i) = eval_cheb(coeffs, xs[i], static_cast<funcval_T>(-1.0),static_cast<funcval_T>(1.0));
    }

    Eigen::Matrix<funcval_T, Eigen::Dynamic, Eigen::Dynamic> A(deg + 1,
                                                               deg + 1);
    for (int i = 0; i <= deg; ++i) {
      for (int j = 0; j <= deg; ++j) {
        A(i, j) = pow(xs[i], j);
      }
    }

    VectorT mon_coeffs = A.inverse() * ys;
    return MonomSegment<funcval_T>(mon_coeffs, begin, end);
  }
};

template <typename funcval_T>
class SegmentedChebyshevApproximator : FunctionApproximator<funcval_T> {
  using VectorT = Eigen::Matrix<funcval_T, Eigen::Dynamic, 1>;
  using MatrixT = Eigen::Matrix<funcval_T, Eigen::Dynamic, Eigen::Dynamic>;
  using ChebyshevSegment = ChebyshevSegment<funcval_T>;

private:
  std::vector<ChebyshevSegment> segments;

  int split_degree;
  double target_precision;

  funcval_T &x_begin = FunctionApproximator<funcval_T>::x_begin;
  funcval_T &x_end = FunctionApproximator<funcval_T>::x_end;
  FunctionApproximator<funcval_T>::funcRR &f =
      FunctionApproximator<funcval_T>::f;

      /// Solving transcendental equations 3.2
  std::tuple<VectorT, VectorT, VectorT>
  chebyshev_interpolate(std::function<funcval_T(funcval_T)> func, int degree,
                        funcval_T a, funcval_T b) {
    VectorT x(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
      // x[i] = (b - a) / 2 * cos((pi * i) / degree) + (b + a) / 2;
      x[i] = cos((pi*i)/degree);
    // std::cout << "Chebysev points:\n" << x << "\n\n";

    VectorT vals(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
      vals(i) = func((b - a) / 2 * x(i) + (b + a) / 2);
    // std::cout << "Function values:\n" << vals << "\n\n";

    MatrixT J(degree + 1, degree + 1);
    for (int j = 0; j < degree + 1; ++j) {
      for (int k = 0; k < degree + 1; ++k) {
        int pj = j == 0 || j == degree ? 2 : 1;
        int pk = k == 0 || k == degree ? 2 : 1;
        J(j, k) = 2.0 / (pj * pk * degree) * cos((j * pi * k) / degree);
      }
    }
    // std::cout << "J:\n" << J << "\n\n";
    //
    Eigen::JacobiSVD<MatrixT> svd(J);
    double cond = svd.singularValues()(0) / (svd.singularValues()(svd.singularValues().size()-1));
    if(cond > 3)
      std::cout<<"cond: "<<cond<<"\n";

    VectorT coeffs = J * vals;
    // std::cout << "Chebysev coefficients:\n" << coeffs << "\n\n";
    return std::make_tuple(coeffs, x, vals);
  }

  std::tuple<double, VectorT, VectorT, funcval_T>
  calc_interstitial_error(std::function<funcval_T(funcval_T)> func,
                          VectorT coeffs, funcval_T a, funcval_T b) {
    int degree = coeffs.size() - 1;
    VectorT x(degree);
    for (int i = 1; i < 2 * degree + 1; i += 2)
      x(i / 2) = (b - a) / 2 * cos((pi * i) / (2 * degree)) + (b + a) / 2;

    VectorT y(degree);
    double interstitial_error = -1;
    funcval_T interstitial_error_place = -1;
    for (int i = 0; i < x.size(); ++i) {
      y(i) = func(x(i));
      double approx = eval_cheb(coeffs, x(i), a, b);
      double err = abs(approx - y(i));
      if (err > interstitial_error){
        interstitial_error = err;
        interstitial_error_place = x(i);
      }
    }
    return std::make_tuple(interstitial_error, x, y, interstitial_error_place);
  }

public:
  SegmentedChebyshevApproximator(FunctionApproximator<funcval_T>::funcRR f,
                                 funcval_T x_begin, funcval_T x_end,
                                 int split_degree = 64,
                                 double target_precision = 5e-10)
      : FunctionApproximator<funcval_T>(f, x_begin, x_end),
        split_degree(split_degree), target_precision(target_precision) {
    Fit();
  };

  virtual void Fit() override {

    funcval_T cur_begin = x_begin;
    funcval_T cur_end = x_end;
    segments.clear();

    double last_err = std::numeric_limits<double>::infinity();

    while (cur_begin + eps<funcval_T>() < x_end) {
      int cur_degree = 1;
      double err = std::numeric_limits<double>::infinity();
      VectorT coeffs;
      double prev_err = std::numeric_limits<double>::infinity();
      bool first = true;
      bool split = false;
      int best_degree = -1;
      double best_error = std::numeric_limits<double>::infinity();
      funcval_T best_error_place = cur_end;
      while (err > target_precision) {
        prev_err = err;

        auto interp_res =
            chebyshev_interpolate(f, cur_degree, cur_begin, cur_end);
        coeffs = std::get<0>(interp_res);
        auto err_res = calc_interstitial_error(f, coeffs, cur_begin, cur_end);
        err = std::get<0>(err_res);
        // std::cout << "d:" << cur_degree << "  e:" << err << "\n";

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
        std::cout << "split: " << cur_begin << "-" << cur_end << " e:" << err
                  << "\n";
      } else {
        auto interp_res =
            chebyshev_interpolate(f, best_degree, cur_begin, cur_end);
        coeffs = std::get<0>(interp_res);
        auto err_res = calc_interstitial_error(f, coeffs, cur_begin, cur_end);
        err = std::get<0>(err_res);
        std::cout << "segment: " << cur_begin << "-" << cur_end
                  << " deg: " << best_degree << " err: " << err << "\n";
        segments.emplace_back();
        auto &seg = segments.back();
        seg.begin = cur_begin;
        seg.end = cur_end;
        seg.coeffs = coeffs;
        cur_begin = cur_end;
        cur_end = x_end;
      }
      // if (cur_degree >= split_degree || prev_err * 2 < err) {
      //	/*if (err > last_err + eps<double>()) {
      //		cur_end = cur_end * 2 - cur_begin;
      //		cur_degree = split_degree;
      //		std::cout << "segment (fail): " << cur_begin << "-" <<
      // cur_end << " deg: " << cur_degree << " err: " << err << "\n";
      // auto interp_res = chebyshev_interpolate(f, cur_degree, cur_begin,
      // cur_end); 		coeffs = std::get<0>(interp_res);
      // segments.emplace_back();
      //		auto& seg = segments.back();
      //		seg.begin = cur_begin;
      //		seg.end = cur_end;
      //		seg.coeffs = coeffs;
      //		cur_begin = cur_end;
      //		cur_end = x_end;
      //		last_err = std::numeric_limits<double>::infinity();
      //	}*/
      //	//else {
      //		cur_end = (cur_begin + cur_end) / 2;
      //		last_err = err;
      //		std::cout << "split: " << cur_begin << "-" << cur_end <<
      //" e:"<<err<< "\n";
      //	//}
      // }
      // else {
      //	std::cout << "segment: " << cur_begin << "-" << cur_end << "
      // deg: " << cur_degree << " err: " << err << "\n";
      //	segments.emplace_back();
      //	auto& seg = segments.back();
      //	seg.begin = cur_begin;
      //	seg.end = cur_end;
      //	seg.coeffs = coeffs;
      //	cur_begin = cur_end;
      //	cur_end = x_end;
      // }
    }
  }

  virtual funcval_T Evaluate(funcval_T x) const override {

    if (x < segments[0].begin)
      throw std::out_of_range("evaluated point is before first segment");
    int i = 0;
    while (i < segments.size() && segments[i].end < x)
      ++i;
    if (i == segments.size())
      throw std::out_of_range("evaluated point is after last segment");

    return segments[i].Evaluate(x);
  };

  virtual void Show() const override {
    hold(true);
    std::vector<double> all_x;
    std::vector<double> correct_y;
    for (const ChebyshevSegment &seg : segments) {
      std::vector<double> x = linspace(seg.begin, seg.end);
      std::vector<double> y(x.size());
      for (int i = 0; i < y.size(); ++i) {
        y[i] = seg.Evaluate(x[i]);
        all_x.push_back(x[i]);
        correct_y.push_back(f(x[i]));
      }
      plot(x, y, "-")->line_width(2);
    }
    auto p = plot(all_x, correct_y, ":");
    p->line_width(1);
    p->display_name("Ground truth");
    // p->visible(false);

    // legend(true);
    show();
  }

  funcval_T FindFirstRoot() const {}

  /// Make sure the chebyshev segments are fitted first
  std::vector<MonomSegment<funcval_T>> GetMonomSegments() const {


    std::vector<MonomSegment<funcval_T>> mon_segments;
    mon_segments.reserve(segments.size());
    for (const auto &seg : segments) {
      // if(seg.coeffs.size() > 5)
      assert(seg.coeffs.size() <= 5);
      mon_segments.push_back(seg.ConvertToMonomial());
    }

    return mon_segments;
  }
};

template <typename funcval_T>
std::pair<funcval_T, funcval_T>
FindQuadraticChebyshevRoot(const ChebyshevSegment<funcval_T> &segment) {
  const auto &a = segment.coeffs;
  int degree = a.size() - 1;
  constexpr funcval_T Tnan = std::numeric_limits<funcval_T>::signaling_NaN();
  constexpr funcval_T Tinf = std::numeric_limits<funcval_T>::infinity();
  if (degree == 0)
    return a[0] < eps<funcval_T>() ? std::make_pair(Tinf, Tinf)
                                   : std::make_pair(Tnan, Tnan);
  if (degree == 1 || degree == 2 && abs(2 * a[2]) < eps<funcval_T>())
    return std::make_pair(-a[0] / a[1], Tnan);
  if (degree == 2) {
    funcval_T D = a[1] * a[1] + 8 * a[2] * a[2] - 8 * a[0] * a[2];
    if (D < 0)
      return std::make_pair(Tnan, Tnan);
    if (D < eps<funcval_T>())
      return std::make_pair(-a[1] / (4 * a[2]), Tnan);
    funcval_T s = -(a[1] + (a[1] >= 0 ? 1 : -1) * sqrt(D));
    funcval_T r1 = s / (4 * a[2]);
    funcval_T r2 =
        (a[0] - a[2]) / (2 * a[2] * r1); // Might be unstable if a0 ~ a2
    return std::make_pair(min(r1, r2), max(r1, r2));
  }
  throw std::invalid_argument("input segment degree must be at most 2");
}
