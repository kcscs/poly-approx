#include "scenelib.hpp"



template<typename FT>
const std::map<std::string, typename Types<FT>::SurfaceFunction> Scenes<FT>::surfaces = {
  {"torus", torus_deg_4()},
  {"sphere", sphere()}
};

template class Scenes<float>;
template class Scenes<double>;

