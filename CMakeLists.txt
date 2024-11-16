cmake_minimum_required(VERSION 3.15)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# SET (CMAKE_C_COMPILER             "/usr/bin/clang")
# SET (CMAKE_C_FLAGS                "-Wall -std=c99")
# SET (CMAKE_C_FLAGS_DEBUG          "-g")
# SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
# SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
# SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

# SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
# SET (CMAKE_CXX_FLAGS                "-Wall")
# SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
# SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
# SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
# SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# SET (CMAKE_AR      "/usr/bin/llvm-ar")
# SET (CMAKE_LINKER  "/usr/bin/llvm-ld")
# SET (CMAKE_NM      "/usr/bin/llvm-nm")
# SET (CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
# SET (CMAKE_RANLIB  "/usr/bin/llvm-ranlib")

project(cheb_renderer)

set(CMAKE_BUILD_TYPE Debug)
add_subdirectory(dependencies/matplotplusplus)

file(GLOB src CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")
include_directories(dependencies/eigen)
include_directories(dependencies/stb)

add_executable(chebrenderer ${src})
target_compile_features(chebrenderer PRIVATE cxx_std_20)
target_compile_options(chebrenderer PRIVATE -g)
# target_compile_options(chebrenderer PRIVATE -Wall -Wextra -Wpedantic)

find_package(OpenMP REQUIRED)
target_link_libraries(chebrenderer PRIVATE OpenMP::OpenMP_CXX)

include(FetchContent)

FetchContent_Declare(
 glm
 GIT_REPOSITORY	https://github.com/g-truc/glm.git
 GIT_TAG 	bf71a834948186f4097caa076cd2663c69a10e1e #refs/tags/1.0.1
)
FetchContent_MakeAvailable(glm)

FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v3.1
)
FetchContent_MakeAvailable(argparse)

FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

target_link_libraries(chebrenderer PRIVATE glm::glm)
target_link_libraries(chebrenderer PUBLIC matplot)
target_link_libraries(chebrenderer PUBLIC argparse)
target_link_libraries(chebrenderer PRIVATE nlohmann_json::nlohmann_json)
