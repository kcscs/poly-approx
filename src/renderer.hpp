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
  void render(T::SurfaceFunction f, CameraConfig cam_conf, ViewConfig view_conf,
              T::gv3 light_dir, std::string outfile, TraceMethod<FT> *tracer) {
    float yedge = tanf(glm::radians(cam_conf.fovy) / 2);
    int width = cam_conf.resolution.x;
    int height = cam_conf.resolution.y;
    float xedge = static_cast<float>(width) / height * yedge;
    std::cout << "edges: " << xedge << " " << yedge << "\n";

    typename T::gm4 invView = glm::inverse(
        glm::lookAt(view_conf.eye, view_conf.at, typename T::gv3(0.0f, 1.0f, 0.0f)));

    std::vector<unsigned char> pixels;
    pixels.reserve(width * height * 3);

    double all_pixels = width * height;

    ns_duration time_spent_on_cheb = ns_duration::zero();
    ns_duration time_spent_on_mon = ns_duration::zero();
    ns_duration time_spent_on_roots = ns_duration::zero();

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        float progress = (x + 1 + y * height) / all_pixels * 100;
        double avg_ns_cheb =
            static_cast<double>(time_spent_on_cheb.count()) / (x + y * height);
        double avg_ns_mon =
            static_cast<double>(time_spent_on_mon.count()) / (x + y * height);
        double avg_ns_roots =
            static_cast<double>(time_spent_on_roots.count()) / (x + y * height);
        int prec = std::cout.precision();
        std::cout.precision(2);
        std::cout << std::fixed << "\r" << (x + 1 + y * height) << "/"
                  << (int)all_pixels << " " << progress << "%"
                  << "  avg_times_ms (c/m/r): " << avg_ns_cheb / 1e6 << " "
                  << avg_ns_mon / 1e6 << " " << avg_ns_roots / 1e6 << "     "
                  << std::defaultfloat;
        std::cout.flush();
        std::cout.precision(prec);
        // std::cout<<x<<";"<<y<<"\n";
        float yf = (static_cast<float>(y) / (height - 1) * 2 - 1) * -1;
        float xf = static_cast<float>(x) / (width - 1) * 2 - 1;

        yf *= yedge;
        xf *= xedge;

        typename T::gv4 ahead = {xf, yf, -1, 1};
        typename T::gv3 dir = typename T::gv3(invView * ahead) - view_conf.eye;

        typename T::Ray r = {view_conf.eye, dir};
        typename TraceMethod<FT>::TraceResult hit =
            tracer->trace(r, f, json()); // TraceMethod<FT>::TraceResult

        // std::cout << "Finished, writing to " << outfile << "\n";
        //     std::cout << "average times per pixel (ms): \n"
        //               << "\tchebyshev interpolation: "
        //               << time_spent_on_cheb.count() / all_pixels / 1e6 <<
        //               "\n";
        //     std::cout << "\tconversion to monomial: "
        //               << time_spent_on_mon.count() / all_pixels / 1e6 <<
        //               "\n";
        //     std::cout << "\troot finding: "
        //               << time_spent_on_roots.count() / all_pixels / 1e6 <<
        //               "\n";

        //     stbi_write_png(outfile.c_str(), width, height, 3, pixels.data(),
        //     width * 3);v
      }
    }
  }
};
