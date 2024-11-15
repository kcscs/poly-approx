#pragma once
#include <functional>
#include <glm/glm.hpp>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <glm/gtc/constants.hpp>
using json = nlohmann::json;

template<typename FT>
struct Types{
  public:
  using gv2 = glm::vec<2, FT>;
  using gv3 = glm::vec<3, FT>;
  using gv4 = glm::vec<4, FT>;

  using gm3 = glm::mat<3,3,FT>;
  using gm4 = glm::mat<4,4,FT>;

  using ev = Eigen::Vector<FT, Eigen::Dynamic>;
  using ev3 = Eigen::Vector<FT, 3>;

  using em = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;

  struct Ray{
    gv3 start;
    gv3 dir;
  };

  using SurfaceFunction = std::function<FT(FT,FT,FT)>;
  using RRFunction = std::function<FT(FT)>;


  static constexpr FT pi = glm::pi<FT>();
  static constexpr FT eps = glm::epsilon<FT>();
};
