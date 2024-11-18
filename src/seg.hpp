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
  FT begin, end;

  Seg(std::vector<FT> coeffs, FT begin, FT end)
      : coeffs(coeffs), begin(begin), end(end) {}

  virtual FT Eval(FT x) const {
    x = (2 * x - (begin + end)) / (end - begin);

    return EvalNorm(x);
  }

  virtual Seg<FT> Differentiate() const { throw NotImplemented(); }

  virtual Seg<FT> Integrate() const { throw NotImplemented(); }

  virtual std::vector<FT> FindRoots(json& metadata) const {
    json norm_rootfind_data;
    norm_rootfind_data["info"] = "TODO";
    std::vector<FT> roots = FindRootsNorm();

    for (auto &x : roots) { // Transform from [-1;1] to [begin;end]
      x = (x / 2 + static_cast<FT>(0.5)) * (end - begin) + begin;
    }

    metadata["normalized_rootfinder"] = norm_rootfind_data;
    metadata["roots"] = roots;
    return roots;
  }

  virtual FT EvalNorm(FT x) const { throw NotImplemented(); }

  virtual std::tuple<FT, ev, ev, FT> Error(T::RRFunction ground_truth) const { throw NotImplemented(); }

  inline size_t deg() const { return coeffs.size() - 1; }

protected:
  /// Evaluate on [-1;1] interval

  virtual std::vector<FT> FindRootsNorm() const { throw NotImplemented(); }
};

template<typename T>
void to_json(json& j, const Seg<T>& s){
  j["domain"] = typename Types<T>::gv2(s.begin, s.end);
  j["coeffs"] = s.coeffs;
  j["degree"] = s.coeffs.size()-1;
}

template<typename T>
void from_json(const json &j, Seg<T>& s) {
  typename Types<T>::gv2 domain = j["domain"];
  s.begin = domain.begin;
  s.end = domain.end;
  s.coeffs = j["coeffs"];
  assert(s.coeffs.size()-1 == j["degree"]);
}
