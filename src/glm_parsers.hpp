#pragma once
#include <glm/glm.hpp>
#include <nlohmann/json.hpp>
#include <vector>

namespace nlohmann {
template <glm::length_t L, typename float_T> struct adl_serializer<glm::vec<L, float_T>> {
  using gvec = glm::vec<L, float_T>;

  static void to_json(json &j, const gvec &vec) {
    std::vector<float_T> d(L);
    for (int i = 0; i < L; ++i) {
      d[i] = vec[i];
    }
    j = json{d};
  }

  static void from_json(const json &j, gvec &vec) {
    std::vector<float> d;
    j.get_to(d);
    gvec v;
    for (int i = 0; i < L; ++i) {
      v[i] = d[i];
    }
    vec = v;
  }
};
} // namespace nlohmann
