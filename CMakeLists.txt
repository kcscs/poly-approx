cmake_minimum_required(VERSION 3.14)

project(poly-approx)

set(CMAKE_BUILD_TYPE Debug)
# find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_subdirectory(dependencies/matplotplusplus)
# add_subdirectory(dependencies/eigen)

file(GLOB src CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")
include_directories(dependencies/eigen)

include_directories(dependencies/stb)

add_executable(polyapprox ${src})
target_compile_features(polyapprox PRIVATE cxx_std_20)

include(FetchContent)

FetchContent_Declare(
 glm
 GIT_REPOSITORY	https://github.com/g-truc/glm.git
 GIT_TAG 	bf71a834948186f4097caa076cd2663c69a10e1e #refs/tags/1.0.1
)

FetchContent_MakeAvailable(glm)

target_link_libraries(polyapprox PRIVATE glm::glm)
target_link_libraries(polyapprox PUBLIC matplot)
# target_link_libraries (polyapprox Eigen3::Eigen)
# target_link_libraries(my_target PUBLIC Eigen3)
