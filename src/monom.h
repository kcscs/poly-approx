#pragma once

#include <Eigen/Dense>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <matplot/matplot.h>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>

using namespace matplot;

// coeffs[0] is constant, coeffs[i] is deg i
template <typename funcval_T>
funcval_T eval_mon(const Eigen::Matrix<funcval_T, Eigen::Dynamic, 1> &coeffs,
                   funcval_T x, funcval_T a, funcval_T b) {

  x = (2 * x - (b + a)) / (b - a); // transform to [-1;1]

  int deg = coeffs.size() - 1;
  if (deg == 0)
    return coeffs(0);
  if (deg == 1)
    return coeffs(0) + coeffs(1) * x;
  funcval_T res = coeffs(deg);
  for (int i = deg - 1; i >= 0; --i) {
    res = res * x + coeffs(i);
  }

  return res;
}

template <typename funcval_T> struct MonomSegment {
  using VectorT = Eigen::Matrix<funcval_T, Eigen::Dynamic, 1>;
  VectorT coeffs;
  funcval_T begin, end;
  funcval_T begin_val = NAN, end_val = NAN;

  MonomSegment<funcval_T>() {}
  MonomSegment<funcval_T>(VectorT coeffs, funcval_T begin, funcval_T end)
      : coeffs(coeffs), begin(begin), end(end) {}
  MonomSegment<funcval_T>(std::initializer_list<funcval_T> coefficients,
                          funcval_T begin, funcval_T end)
      : begin(begin), end(end), coeffs(coefficients.size()) {
    auto iter = coefficients.begin();
    for (int i = 0; i < coefficients.size(); ++i) {
      coeffs[i] = *iter;
      iter++;
    }
  }

  funcval_T Evaluate(funcval_T x) const {
    return eval_mon(coeffs, x, begin, end);
  }

  void Show(bool open_window = true) const {
    std::vector<double> x = linspace(begin, end);
    std::vector<double> y(x.size());
    for (int i = 0; i < y.size(); ++i) {
      y[i] = Evaluate(x[i]);
    }
    auto p = plot(x, y);
    p->line_width(2);
    // p->color("red");
    // if(open_window)
    //     show();
  }

  MonomSegment<funcval_T> GetDerivative() const {
    const auto &c = coeffs;
    VectorT dc(c.size() - 1);
    int ddeg = c.size() - 2;
    for (int i = 0; i <= ddeg; ++i) {
      dc(i) = c(i + 1) * (i + 1);
    }

    return MonomSegment<funcval_T>(dc, begin, end);
  }

  // works on [0;1]
  std::vector<funcval_T> FindRoots01() const {
    constexpr funcval_T eps = std::numeric_limits<funcval_T>::epsilon();
    const auto &c = coeffs;
    const int deg = c.size() - 1;
    // std::cout<<"fr: "<<begin<<"-"<<end<<" d: "<<deg<<"\n";
    if (deg == 0)
      throw std::invalid_argument("Infinite roots");

    if (deg == 1) {
      funcval_T r = -c(0) / c(1);
      return abs(r) <= 1 ? std::vector<funcval_T>({r})
                         : std::vector<funcval_T>();
    }

    if (deg == 2) {
      funcval_T D = c(1) * c(1) - 4 * c(2) * c(0);

      if (D < -eps)
        return std::vector<funcval_T>();
      if (abs(D) < eps) {
        funcval_T r = -c(1) / (2 * c(2));
        return abs(r) <= 1 ? std::vector<funcval_T>({r})
                           : std::vector<funcval_T>();
      } else {
        funcval_T sqrtD = sqrt(D);
        funcval_T sgnb = c(1) < 0 ? -1 : 1;
        funcval_T r1 = -2 * c(0) / (c(1) + sgnb * sqrtD);
        funcval_T r2 = -(c(1) + sgnb * sqrtD) / (2 * c(2));

        std::vector<funcval_T> roots;
        if (abs(r1) <= 1)
          roots.push_back(r1);
        if (abs(r2) <= 1)
          roots.push_back(r2);
        if (roots.size() == 2 && roots[0] > roots[1])
          std::swap(roots[0], roots[1]);
        return roots;
      }
    } else {
      MonomSegment<funcval_T> derivative = GetDerivative();
      std::vector<funcval_T> critical_points = derivative.FindRoots01();
      // std::cout<<"d: "<<deg<<" crits: "<<critical_points.size()<<"\n";
      std::vector<funcval_T> borders = {-1};
      borders.insert(borders.end(), critical_points.begin(),
                     critical_points.end());
      borders.push_back(1);

      std::vector<funcval_T> border_values(borders.size());
      for (int i = 0; i < borders.size(); ++i) {
        border_values[i] = eval_mon(c, borders[i], static_cast<funcval_T>(-1.0),static_cast<funcval_T>(1.0));
      }

      std::vector<funcval_T> roots;
      for (int i = 0; i < borders.size() - 1; ++i) {
        // process segment
        funcval_T xl = borders[i], xr = borders[i + 1], vl = border_values[i],
                  vr = border_values[i + 1];
        if (vl * vr > eps)
          continue;

        funcval_T xc = (borders[i] + borders[i + 1]) / 2;
        funcval_T vc = eval_mon(c, xc, static_cast<funcval_T>(-1.0),static_cast<funcval_T>(1.0));
        if (vl * vc > eps) {
          vl = vc;
          xl = xc;
        } else {
          vr = vc;
          xr = xc;
        }

        // xc = (xl + xr) / 2;
        funcval_T xn = (xl + xr) / 2, xn_prev = -2;
        int max_steps = 2000;
        int step = 0;
        while (abs(xn - xn_prev) > eps*50 && step < max_steps) { // newton-bisection hybrid iterations
          while (abs(xn - xn_prev) > eps*50 && xn >= xl && xn <= xr && step < max_steps) {
            xn_prev = xn;
            xn = xn - eval_mon(c, xn, static_cast<funcval_T>(-1.0),static_cast<funcval_T>(1.0)) /
                          eval_mon(derivative.coeffs, xn, static_cast<funcval_T>(-1.0),static_cast<funcval_T>(1.0));
            ++step;
          }
          if (xn < xl){
            xn = (xl+xn_prev)/2;
            xr = xn_prev;
            ++step;
          }
          else if (xn > xr) {
            xn = (xn_prev+xr)/2;
            xl = xn_prev;
            ++step;
          }
        }
        // if(step == max_steps){
        //   std::cout<<"warning: early stop with last change: "<<abs(xn-xn_prev)<<"\n";
        // }
        roots.push_back(xn);
      }
      return roots;
    }
  }

  std::vector<funcval_T> FindRoots() const { // With Yuksel's method (High-Performance
                                    // Polynomial Root Finding for Graphics)
    std::vector<funcval_T> roots = FindRoots01();
    for (funcval_T& r : roots) {
      r = (r+1)/2*(end-begin)+begin;
    }
    return roots;
  }
};

template<typename funcval_T>
std::ostream& operator<<(std::ostream& os, const MonomSegment<funcval_T>& ms) {
  os << "---\nPolynomial of degree "<<ms.coeffs.size()-1<<" representing "<<ms.begin<<" to "<<ms.end<<"\n";
  os <<"  "<< ms.coeffs(0)<<"\n";
  for(int i = 1; i < ms.coeffs.size(); ++i){
    os <<"+ "<< ms.coeffs(i)<<" * x^"<<i<<"\n";
  }
  return os<<"---\n";
}
