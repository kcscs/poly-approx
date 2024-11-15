#pragma once
#include "seg.hpp"
#include <glm/gtc/constants.hpp>

template <typename FT> struct ChebSeg : public Seg<FT> {
  using T = Types<FT>;
  using Seg<FT>::coeffs;
  using Seg<FT>::begin;
  using Seg<FT>::end;
  using typename T::em;
  using typename T::ev;

public:
  ChebSeg(std::vector<FT> coeffs, FT begin, FT end)
      : Seg<FT>(coeffs, begin, end) {}

  FT EvalNorm(FT x) const override {

    int deg = coeffs.size() - 1;
    if (deg == 0)
      return coeffs[0];
    if (deg == 1)
      return coeffs[0] + coeffs[1] * x;

    FT Tprev = x;
    FT Tcur = 2 * x * x - 1;
    FT res = coeffs[0] + coeffs[1] * x;
    for (int i = 2; i < deg; ++i) {
      res += coeffs[i] * Tcur;
      FT Ttemp = Tcur;
      Tcur = 2 * x * Tcur - Tprev;
      Tprev = Ttemp;
    }
    res += coeffs[deg] * Tcur;

    return res;
  }

  std::tuple<FT, ev, ev, FT> Error(T::RRFunction ground_truth)
      const override { // return type should have ET for the error
    using ET = FT;
    int degree = this->deg();
    ev x(degree);
    for (int i = 1; i < 2 * degree + 1; i += 2)
      x(i / 2) = (end - begin) / 2 * cos((T::pi * i) / (2 * degree)) +
                 (end + begin) / 2;

    ev y(degree);
    ET interstitial_error = -1;
    FT interstitial_error_place = -1;
    for (int i = 0; i < x.size(); ++i) {
      y(i) = ground_truth(x(i));
      FT approx = this->Eval(x(i));
      ET err = abs(static_cast<ET>(approx) - y(i));
      if (err > interstitial_error) {
        interstitial_error = err;
        interstitial_error_place = x(i);
      }
    }
    return std::make_tuple(interstitial_error, x, y, interstitial_error_place);
  }

  /// Solving transcendental equations 3.2
  static ChebSeg<FT> Interpolate(T::RRFunction func, int degree, FT a, FT b) {
    ev x(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
      // x[i] = (b - a) / 2 * cos((pi * i) / degree) + (b + a) / 2;
      x[i] = cos((T::pi * i) / degree);
    // std::cout << "Chebysev points:\n" << x << "\n\n";

    ev vals(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
      vals(i) = func((b - a) / 2 * x(i) + (b + a) / 2);
    // std::cout << "Function values:\n" << vals << "\n\n";

    em J(degree + 1, degree + 1);
    for (int j = 0; j < degree + 1; ++j) {
      for (int k = 0; k < degree + 1; ++k) {
        int pj = j == 0 || j == degree ? 2 : 1;
        int pk = k == 0 || k == degree ? 2 : 1;
        J(j, k) = 2.0 / (pj * pk * degree) * cos((j * T::pi * k) / degree);
      }
    }
    // std::cout << "J:\n" << J << "\n\n";
    //
    Eigen::JacobiSVD<em> svd(J);
    double cond = svd.singularValues()(0) /
                  (svd.singularValues()(svd.singularValues().size() - 1));
    // if (cond > 3)
    //   std::cout << "cond: " << cond << "\n";

    ev coeffs = J * vals;
    // std::cout << "Chebysev coefficients:\n" << coeffs << "\n\n";
    return std::make_tuple(coeffs, x, vals);
  }
};
