#pragma once
#include "logging.hpp"
#include "seg.hpp"
#include <glm/gtc/constants.hpp>
#include <iostream>
#include <limits>
#include <cmath>
#include <format>

template <typename FT> struct ChebSeg : public Seg<FT> {
  using T = Types<FT>;
  using Seg<FT>::coeffs;
  using Seg<FT>::begin;
  using Seg<FT>::end;
  using Seg<FT>::min_val;
  using Seg<FT>::max_val;
  using typename T::em;
  using typename T::ev;

public:
  ChebSeg(std::vector<FT> coeffs, FT begin, FT end, FT min_val = 0, FT max_val = 1)
      : Seg<FT>(coeffs, begin, end, min_val, max_val) {}

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

  

  /// Solving transcendental equations 3.2
  static ChebSeg<FT> Interpolate(T::RRFunction func, int degree, FT a, FT b) {
    Logger &log = Logger::Get();
    log << Logger::cat("cheb_interp");
    ev x(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
      // x[i] = (b - a) / 2 * cos((pi * i) / degree) + (b + a) / 2;
      x[i] = cos((T::pi * i) / degree);
    log << "Chebysev points:\n" << x << "\n\n";

    for (int i = 0; i < degree + 1; ++i)
      x[i] = (b - a) / 2 * x(i) + (b + a) / 2;
    log << "Transformed points:\n" << x << "\n\n";

    ev vals(degree + 1);
    FT min_val = std::numeric_limits<FT>::max();
    FT max_val = std::numeric_limits<FT>::min();
    for (int i = 0; i < degree + 1; ++i){
      vals(i) = func(x[i]);
      if(vals(i) < min_val)
        min_val = vals(i);
      if(vals(i) > max_val)
        max_val = vals(i);
    }

    FT range = max_val - min_val;
    log << "Function values:\n";
    for (int i = 0; i < degree + 1; ++i){
      log << vals(i);
      FT ulp = std::nextafter(vals(i), std::numeric_limits<FT>::infinity()) - vals(i);
      log << " ulp: "<<std::format("{:.0e}", ulp)<< "   ";
      vals(i) = (vals(i)-min_val)/range;
      log << " (" << vals(i);
      ulp = std::nextafter(vals(i), std::numeric_limits<FT>::infinity()) - vals(i);
      log << " ulp: "<<std::format("{:.0e}", ulp)<< ")\n";
    }

    log << "Function range: "<<min_val<<"-"<<max_val<<"\n";
    log << "Normalizing to 0-1\n";

    

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
    log << "Chebysev coefficients:\n" << coeffs << "\n\n";

    std::vector<FT> coeff_vec(coeffs.size());
    for (int i = 0; i < coeffs.size(); i++) {
      coeff_vec[i] = coeffs(i);
    }
    return {coeff_vec, a, b, min_val, max_val};
  }
};
