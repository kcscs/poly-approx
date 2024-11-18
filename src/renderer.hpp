#pragma once
#include "experiment_config.hpp"
#include "fun.hpp"
#include "monseg.hpp"
#include "tracer.hpp"
#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ios>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <Eigen/Dense>
#include <chrono>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <iostream>
#include <string>

template <typename FT> class Renderer {
  using ns_duration = std::chrono::nanoseconds;
  using T = Types<FT>;
  // using typename T::gv3;
  // using typename Types<FT>::ev3;
  // using typename Types<FT>::gm4;
  // using typename Types<FT>::gv3;
  // using typename Types<FT>::gv4;
  // using typename Types<FT>::Ray;
  // using typename Types<FT>::RRFunction;
  // using typename Types<FT>::SurfaceFunction;

public:
  json render(T::SurfaceFunction f, CameraConfig cam_conf, ViewConfig view_conf,
              T::gv3 light_dir, std::string outfile, TraceMethod<FT> *tracer) {
    Logger &log = Logger::Get();
    float yedge = tanf(glm::radians(cam_conf.fovy) / 2);
    int width = cam_conf.resolution.x;
    int height = cam_conf.resolution.y;
    float xedge = static_cast<float>(width) / height * yedge;
    std::cout << "edges: " << xedge << " " << yedge << "\n";

    json results;

    typename T::gv3 eye(view_conf.eye);
    typename T::gv3 at(view_conf.at);

    typename T::gm4 invView =
        glm::inverse(glm::lookAt(eye, at, typename T::gv3(0.0f, 1.0f, 0.0f)));

    std::vector<unsigned char> pixels;
    pixels.reserve(width * height * 3);

    std::vector<json> pixels_metadata;
    pixels_metadata.reserve(width * height * 3);

    double all_pixels = width * height;

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        float progress = (x + 1 + y * height) / all_pixels * 100;
        std::cout << std::fixed << "\r" << (x + 1 + y * height) << "/"
                  << (int)all_pixels << " " << progress << "%    ";
        std::cout.flush();
        // std::cout<<x<<";"<<y<<"\n";

        json pixel_data;
        pixel_data["coordinate"] =
            glm::uvec2(x, y); // x: horizontal rightwards, y: vertical downwards

        float yf = (static_cast<float>(y) / (height - 1) * 2 - 1) * -1;
        float xf = static_cast<float>(x) / (width - 1) * 2 - 1;

        yf *= yedge;
        xf *= xedge;

        typename T::gv4 ahead = {xf, yf, -1, 1};
        typename T::gv3 dir = typename T::gv3(invView * ahead) - eye;

        typename T::Ray r = {eye, glm::normalize(dir)};
        log << "render"_cat << "Ray: " << "d:" << r.dir.x << " " << r.dir.y
            << " " << r.dir.z << " | o:" << r.start.x << " " << r.start.y << " "
            << r.start.z << "\n";
        typename TraceMethod<FT>::TraceResult hit =
            tracer->trace(r, f); // TraceMethod<FT>::TraceResult
        pixel_data["hit"] = hit.hit;
        pixel_data["distance"] = hit.distance;
        pixel_data["trace_data"] = hit.metadata;

        if (hit.hit) {

          typename T::gv3 wp = r.start + r.dir * hit.distance;
          typename T::gv3 norm;
          norm.x = f(wp.x + 0.01, wp.y, wp.z) - f(wp.x - 0.01, wp.y, wp.z);
          norm.y = f(wp.x, wp.y + 0.01, wp.z) - f(wp.x, wp.y - 0.01, wp.z);
          norm.z = f(wp.x, wp.y, wp.z + 0.01) - f(wp.x, wp.y, wp.z - 0.01);

          typename T::gv3 to_eye = glm::normalize(eye - wp);
          norm = glm::normalize(norm);
          if (glm::dot(norm, to_eye) < 0)
            norm *= -1;

          if (y == height / 4 && x == width / 2) {
            std::cout << "worldpos: " << wp.x << " " << wp.y << " " << wp.z
                      << "\n";
            std::cout << "normal: " << norm.x << " " << norm.y << " " << norm.z
                      << "\n";
          }

          float diffuse = glm::dot(norm, -light_dir);
          if (diffuse < 0)
            diffuse = 0;

          float specular = glm::pow(
              glm::clamp(glm::dot(to_eye, glm::reflect(light_dir, norm)),
                         static_cast<FT>(0.0), static_cast<FT>(1.0)),
              27);
          float col = glm::clamp(diffuse + specular, 0.0f, 1.0f);
          pixels.push_back(col * 255);
          pixels.push_back(col * 255);
          pixels.push_back(col * 255);
        } else {
          pixels.push_back(255);
          pixels.push_back(0);
          pixels.push_back(0);
        }

        pixels_metadata.push_back(pixel_data);
      }
    }
    std::cout << "Finished, writing to " << outfile << "\n";

    stbi_write_png(outfile.c_str(), width, height, 3, pixels.data(), width * 3);

    results["pixels"] = pixels_metadata;
    return results;
  }
};
