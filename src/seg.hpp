#pragma once

#include "exceptions.hpp"
#include "run_config.hpp"
#include "types.hpp"
#include <vector>

template <typename FT> struct Seg : Types<FT> {
  using T = Types<FT>;
  using typename T::ev;

public:
  /// Contains the coefficients in increasing order by degree
  std::vector<FT> coeffs;
  FT begin, end, min_val, max_val;
  json metadata;

  Seg(std::vector<FT> coeffs, FT begin, FT end, FT min_val = 0, FT max_val = 1)
      : coeffs(coeffs), begin(begin), end(end), min_val(min_val), max_val(max_val) {}

  virtual FT Eval(FT x) const {
    x = (2 * x - (begin + end)) / (end - begin);
    FT y = EvalNorm(x);
    y = y*(max_val-min_val) + min_val;
    return y;
  }

  virtual Seg<FT> Differentiate() const { throw NotImplemented(); }

  virtual Seg<FT> Integrate() const { throw NotImplemented(); }

  virtual std::vector<FT> FindRoots(json& metadata) const {
    json norm_rootfind_data;
    std::vector<FT> roots = FindRootsNorm(norm_rootfind_data);

    metadata["normalized_rootfinder"] = norm_rootfind_data;
    metadata["normalized_roots"] = roots;

    for (auto &x : roots) { // Transform from [-1;1] to [begin;end]
      x = (x / static_cast<FT>(2.0) + static_cast<FT>(0.5)) * (end - begin) + begin;
    }

    metadata["roots"] = roots;
    return roots;
  }

  virtual FT EvalNorm(FT x) const { throw NotImplemented(); }

  virtual std::tuple<FT, ev, ev, FT> InterstitialError(T::RRFunction ground_truth) const { // return type should have ET for the error
      using ET = FT;
      Logger& log = Logger::Get();
      log << Logger::cat("interstitial");

      int degree = this->deg();
      ev x(degree);
      for (int i = 1; i < 2 * degree + 1; i += 2)
          x(i / 2) = cos((T::pi * i) / (2 * degree));

      ev y(degree);
      ET interstitial_error = -1;
      FT interstitial_error_place = -1;
      log << "calculating error\ncoeffs: " << coeffs << "\n";
      log << "scaling min and max values: " << min_val << " " << max_val << "\n";
      for (int i = 0; i < x.size(); ++i) {
          FT x2 = (end - begin) / 2 * x(i) + (end + begin) / 2;
          y(i) = (ground_truth(x2) - min_val) / (max_val - min_val);
          FT approx = this->EvalNorm(x(i));
          ET err = abs(static_cast<ET>(approx) - y(i));
          log << "at " << x(i) << " gt:" << y(i) << " approx:" << approx
              << " err:" << err << "\n";
          if (err > interstitial_error) {
              interstitial_error = err;
              interstitial_error_place = x(i);
          }
      }
      return std::make_tuple(interstitial_error, x, y, interstitial_error_place);
  }

  inline size_t deg() const { return coeffs.size() - 1; }

protected:
  /// Evaluate on [-1;1] interval

  virtual std::vector<FT> FindRootsNorm(json& metadata) const { throw NotImplemented(); }
};

template<typename T>
void to_json(json& j, const Seg<T>& s){
  j["domain"] = typename Types<T>::gv2(s.begin, s.end);
  j["coeffs"] = s.coeffs;
  j["degree"] = s.coeffs.size()-1;
  j["range"] = typename Types<T>::gv2(s.min_val, s.max_val);
  j["metadata"] = s.metadata;
}

template<typename T>
void from_json(const json &j, Seg<T>& s) {
  typename Types<T>::gv2 domain = j["domain"];
  s.begin = domain.begin;
  s.end = domain.end;
  s.coeffs = j["coeffs"];
  assert(s.coeffs.size()-1 == j["degree"]);
  typename Types<T>::gv2 range = j["range"];
  s.min_val = range.x;
  s.max_val = range.y;
  s.metadata = j["metadata"];
}
