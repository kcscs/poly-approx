#pragma once
#include "chebseg.hpp"
#include "logging.hpp"
#include "seg.hpp"
#include "types.hpp"
#include <cassert>
#include <glm/gtc/constants.hpp>
#include <vector>

template <typename FT> struct MonSeg : public Seg<FT> {
  using typename Seg<FT>::T::em;
  using typename Seg<FT>::T::ev;
  using Seg<FT>::coeffs;
  using Seg<FT>::begin;
  using Seg<FT>::end;
  using Seg<FT>::min_val;
  using Seg<FT>::max_val;

public:
  MonSeg(std::vector<FT> coeffs, FT begin, FT end, FT min_val = 0,
         FT max_val = 1)
      : Seg<FT>(coeffs, begin, end, min_val, max_val) {}

  static MonSeg<FT> FitAtChebPoints(const Seg<FT> &other) {
    const int deg = other.coeffs.size() - 1;
    std::vector<FT> xs(deg + 1);
    ev ys(deg + 1);
    for (int i = 0; i <= deg; ++i) {
      // xs[i] = (static_cast<funcval_T>(i) / deg)*2-1; //
      xs[i] = -cos(static_cast<FT>(i) / deg * glm::pi<FT>());
      ys(i) = other.EvalNorm(xs[i]);
    }

    em A(deg + 1, deg + 1);
    for (int i = 0; i <= deg; ++i) {
      for (int j = 0; j <= deg; ++j) {
        A(i, j) = pow(xs[i], j);
      }
    }

    ev mon_coeffs = A.inverse() * ys;

    std::vector<FT> mon_coeffs_vec(mon_coeffs.size());
    for (int i = 0; i < mon_coeffs.size(); ++i) {
      mon_coeffs_vec[i] = mon_coeffs(i);
    }

    return MonSeg<FT>(mon_coeffs_vec, other.begin, other.end, other.min_val, other.max_val);
  }

  Seg<FT> Differentiate() const override { return TDifferentiate(); }

  MonSeg<FT> TDifferentiate() const {
    const auto &c = coeffs;
    std::vector<FT> dc(c.size() - 1);
    int ddeg = c.size() - 2;
    for (int i = 0; i <= ddeg; ++i) {
      dc[i] = c[i + 1] * (i + 1) * (max_val-min_val);
    }

    return MonSeg<FT>(dc, begin, end, 0, 1);
  }

  FT EvalNorm(FT x) const override { return eval_mon(coeffs, x); }

private:
  static FT eval_mon(const std::vector<FT> &coeffs, FT x) {
    int deg = coeffs.size() - 1;
    if (deg == 0)
      return coeffs[0];
    if (deg == 1)
      return coeffs[0] + coeffs[1] * x;
    FT res = coeffs[deg];
    for (int i = deg - 1; i >= 0; --i) {
      res = res * x + coeffs[i];
    }

    return res;
  }

  virtual std::vector<FT> FindRootsNorm() const override {
    Logger &log = Logger::Get();
    constexpr FT eps = std::numeric_limits<FT>::epsilon();
    const auto &c = coeffs;
    assert(c.size() > 0);
    const int deg = c.size() - 1;
    log << Logger::cat("rootfinder");
    log << "Note: FindRootsNorm works without normalized range\n";
    log << "fr " << begin << "-" << end << " d: " << deg << "\n";
    if (deg == 0)
      throw std::invalid_argument("Infinite roots");

    if (deg == 1) {
      // FT r = (-min_val/(max_val-min_val)  -c[0]) / c[1];
      FT r = (-min_val/(max_val-min_val)  -c[0]) / c[1];
      return abs(r) <= 1 ? std::vector<FT>({r}) : std::vector<FT>();
    }

    if (deg == 2) {
      // FT D = c[1] * c[1] - 4 * c[2] * c[0];
      FT ran = max_val-min_val;
      FT sqrtran = glm::sqrt(ran);
      FT D = ran*(c[1]*c[1]-4*c[2]*c[0])-min_val*4*c[2];

      if (D < -eps)
        return std::vector<FT>();
      if (abs(D) < eps) {
        // FT r = -c[1] / (2 * c[2]);
        
        FT r = -c[1]*sqrtran/(2*c[2]*sqrtran);
        return abs(r) <= 1 ? std::vector<FT>({r}) : std::vector<FT>();
      } else {
        FT sqrtD = sqrt(D);
        FT sgnb = c[1] < 0 ? -1 : 1;
        // FT r1 = -2 * c[0] / (c[1] + sgnb * sqrtD);
        // FT r2 = -(c[1] + sgnb * sqrtD) / (2 * c[2]);

        FT r1 = (-2*c[0]*ran-2*min_val)/(c[1]*ran+sgnb*sqrtran*sqrtD);
        FT r2 = -(c[1]*sqrtran+sgnb*sqrtD)/(2*c[2]*sqrtran);

        std::vector<FT> roots;
        if (abs(r1) <= 1)
          roots.push_back(r1);
        if (abs(r2) <= 1)
          roots.push_back(r2);
        if (roots.size() == 2 && roots[0] > roots[1])
          std::swap(roots[0], roots[1]);
        return roots;
      }
    } else {
      MonSeg<FT> derivative = TDifferentiate();
      std::vector<FT> critical_points = derivative.FindRootsNorm();
      log << "d: " << deg << " crits: " << critical_points.size() << "\n";
      std::vector<FT> borders = {-1};
      borders.insert(borders.end(), critical_points.begin(),
                     critical_points.end());
      borders.push_back(1);

      std::vector<FT> border_values(borders.size());
      for (int i = 0; i < borders.size(); ++i) {
        border_values[i] = MonSeg<FT>(c, static_cast<FT>(-1.0),
                                      static_cast<FT>(1.0), min_val, max_val)
                               .Eval(borders[i]);
      }

      std::vector<FT> roots;
      for (int i = 0; i < borders.size() - 1; ++i) {
        // process segment
        FT xl = borders[i], xr = borders[i + 1], vl = border_values[i],
           vr = border_values[i + 1];
        if (vl * vr > eps)
          continue;

        FT xc = (borders[i] + borders[i + 1]) / 2;
        // FT vc = Eval(c, xc, static_cast<FT>(-1.0), static_cast<FT>(1.0));
        FT vc = MonSeg<FT>(c, static_cast<FT>(-1.0), static_cast<FT>(1.0),
                           min_val, max_val)
                    .Eval(xc);
        if (vl * vc > eps) {
          vl = vc;
          xl = xc;
        } else {
          vr = vc;
          xr = xc;
        }

        // xc = (xl + xr) / 2;
        FT xn = (xl + xr) / 2, xn_prev = -2;
        int max_steps = 2000;
        int step = 0;
        while (abs(xn - xn_prev) > eps * 50 &&
               step < max_steps) { // newton-bisection hybrid iterations
          while (abs(xn - xn_prev) > eps * 50 && xn >= xl && xn <= xr &&
                 step < max_steps) {
            xn_prev = xn;
            // xn = xn -
            //      Eval(c, xn, static_cast<FT>(-1.0), static_cast<FT>(1.0)) /
            //          Eval(derivative.coeffs, xn, static_cast<FT>(-1.0),
            //                   static_cast<FT>(1.0));
            xn = xn - MonSeg<FT>(c, static_cast<FT>(-1.0), static_cast<FT>(1.0),
                                 min_val, max_val)
                              .Eval(xn) /
                          MonSeg<FT>(derivative.coeffs, static_cast<FT>(-1.0),
                                     static_cast<FT>(1.0), min_val, max_val)
                              .Eval(xc);

            ++step;
          }
          if (xn < xl) {
            xn = (xl + xn_prev) / 2;
            xr = xn_prev;
            ++step;
          } else if (xn > xr) {
            xn = (xn_prev + xr) / 2;
            xl = xn_prev;
            ++step;
          }
        }
        // if(step == max_steps){
        //   std::cout<<"warning: early stop with last change:
        //   "<<abs(xn-xn_prev)<<"\n";
        // }
        roots.push_back(xn);
      }
      return roots;
    }
  }
};
