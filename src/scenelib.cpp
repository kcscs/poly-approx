#include "scenelib.hpp"

template <typename FT>
const std::map<std::string, typename Types<FT>::SurfaceFunction>
    Scenes<FT>::surfaces = {{"torus", torus_deg_4()},
                            {"sphere", sphere()},
                            {"cubic_sheet", cubic_sheet()},
                            {"cubic_holes", cubic_holes()},
                            {"riemann", riemann()},
                            {"kummer", kummer()},
                            {"barth_sextic", barth_sextic()},
                            {"sphere_sdf", sphere_sdf()},
                            {"box", box()},
                            {"torus_sdf", torus_sdf()}};

template class Scenes<float>;
template class Scenes<double>;
