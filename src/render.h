#include "SegmentedChebyshevApproximator.h"
#include "monom.h"
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

template <typename funcval_T> class Renderer {
  using ns_duration = std::chrono::nanoseconds;
  using glmvec2 = glm::vec<2, funcval_T>;
  using glmvec3 = glm::vec<3, funcval_T>;
  using glmvec4 = glm::vec<4, funcval_T>;
  using glmmat4 = glm::mat<4,4,funcval_T>;

public:
  using Vec3 = Eigen::Matrix<funcval_T, 3, 1>;

  void render(std::function<funcval_T(funcval_T, funcval_T, funcval_T)> f,
              glmvec3 eye, glmvec3 center, int width, int height,
              float fovy_rad, glmvec3 light_dir, std::string outfile) {
    float yedge = tanf(fovy_rad / 2);
    float xedge = static_cast<float>(width) / height * yedge;
    std::cout << "edges: " << xedge << " " << yedge << "\n";

    glmmat4 invView =
        glm::inverse(glm::lookAt(eye, center, glmvec3(0.0f, 1.0f, 0.0f)));

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

        glmvec4 ahead = {xf, yf, -1, 1};
        glmvec3 dir = glmvec3(invView * ahead) - eye;

        std::function<funcval_T(funcval_T)> func = [&](funcval_T t) {
          glmvec3 p = eye + t * dir;
          return f(p.x, p.y, p.z);
        };

        // std::cout<<"Fitting chebyshev\n";

        auto t_cheb_begin = std::chrono::high_resolution_clock::now();
        // SegmentedChebyshevApproximator<funcval_T> cheb(func, 0, 20, 64, 5e-14);
        // SegmentedChebyshevApproximator<funcval_T> cheb(func, 0, 20, 64, 5e-14);
        SegmentedChebyshevApproximator<funcval_T> cheb(func, 0, 20, 16, 5e-13);
        // std::cout<<"Approximated\n";
        auto t_mon_begin = std::chrono::high_resolution_clock::now();
        auto MonSegments = cheb.GetMonomSegments();

        auto t_roots_begin = std::chrono::high_resolution_clock::now();
        std::vector<funcval_T> roots;
        for (const auto &s : MonSegments) {
          // std::cout<<"root finding "<<s;
          auto r = s.FindRoots();
          roots.insert(roots.end(), r.begin(), r.end());
        }
        auto t_roots_end = std::chrono::high_resolution_clock::now();
        time_spent_on_cheb += t_mon_begin - t_cheb_begin;
        time_spent_on_mon += t_roots_begin - t_mon_begin;
        time_spent_on_roots += t_roots_end - t_roots_begin;

        if (y == height / 4 && x == width / 2) {
          std::cout << "dir: " << dir.x << " " << dir.y << " " << dir.z << "\n";
          figure();
          hold(on);
          for (const auto &s : MonSegments) {
            s.Show();
          }
          for (auto r : roots) {
            std::cout << "root: " << r << "\n";
          }

          std::vector<double> rootzeros(roots.size());
          auto p = plot(roots, rootzeros, "x");
          p->line_width(3);
          show();
        }

        if (roots.size() == 0) {
          pixels.push_back(0);
          pixels.push_back(0);
          pixels.push_back(255);
        } else {
          glmvec3 wp = eye + roots[0] * dir;
          glmvec3 norm;
          norm.x = f(wp.x + 0.01, wp.y, wp.z) - f(wp.x - 0.01, wp.y, wp.z);
          norm.y = f(wp.x, wp.y + 0.01, wp.z) - f(wp.x, wp.y - 0.01, wp.z);
          norm.z = f(wp.x, wp.y, wp.z + 0.01) - f(wp.x, wp.y, wp.z - 0.01);

          glmvec3 to_eye = glm::normalize(eye - wp);
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
              glm::clamp(glm::dot(to_eye, glm::reflect(light_dir, norm)), static_cast<funcval_T>(0.0),static_cast<funcval_T>(1.0)),
              27);
          float col = glm::clamp(diffuse + specular, 0.0f, 1.0f);
          pixels.push_back(roots.size());
          pixels.push_back(col * 255);
          pixels.push_back(col * 255);
        }
      }
    }

    std::cout << "Finished, writing to " << outfile << "\n";
    std::cout << "average times per pixel (ms): \n"
              << "\tchebyshev interpolation: "
              << time_spent_on_cheb.count() / all_pixels / 1e6 << "\n";
    std::cout << "\tconversion to monomial: "
              << time_spent_on_mon.count() / all_pixels / 1e6 << "\n";
    std::cout << "\troot finding: "
              << time_spent_on_roots.count() / all_pixels / 1e6 << "\n";

    stbi_write_png(outfile.c_str(), width, height, 3, pixels.data(), width * 3);
  }

  void
  raymarch_test_ray(std::function<funcval_T(funcval_T, funcval_T, funcval_T)> f,
                    glmvec3 start, glmvec3 dir, glmvec2 interp_range,
                    glmvec2 test_range, int resolution, int max_degree = 32,
                    funcval_T interp_max_error = 5e-2) {
    std::cout << "raymarch test\n";
    std::function<funcval_T(funcval_T)> func = [&](funcval_T t) {
      glmvec3 p = start + t * dir;
      return f(p.x, p.y, p.z);
    };

    std::cout << "interpolating\n";
    SegmentedChebyshevApproximator<funcval_T> cheb(
        func, interp_range.x, interp_range.y, max_degree, interp_max_error);

    auto monom = cheb.GetMonomSegments();
    std::cout<<"segments: "<<monom.size()<<"\n";
    for(auto& m : monom)
      std::cout<<m.coeffs.size()<<"\n";
    std::cout.flush();

    std::vector<funcval_T> errors(resolution);
    for (int i = 0; i < resolution; ++i) {
      float progress = i / (resolution - 1.0f) * 100;
      int prec = std::cout.precision();
      std::cout.precision(2);
      std::cout << std::fixed << "\r" << progress << "%    " << std::defaultfloat;
      std::cout.flush();
      std::cout.precision(prec);

      funcval_T x = interp_range.x +
                    (interp_range.y - interp_range.x) / (resolution - 1) * i;
      funcval_T y = func(x);
      funcval_T y_est = cheb.Evaluate(x);

      MonomSegment<funcval_T> search_val;
      search_val.begin = search_val.end = x;
      auto mon_segment = std::lower_bound(
          monom.begin(), monom.end(), search_val,
          [&](const MonomSegment<funcval_T> &a,
              const MonomSegment<funcval_T> &b) { return a.end < b.end; });

      funcval_T mon_est = (*mon_segment).Evaluate(x);


      funcval_T err = abs(y - y_est);
      if (err > interp_max_error)
        std::cout << x << " cheb err " << err << "\n";
      funcval_T mon_err = abs(y - mon_est);
      if (mon_err > interp_max_error)
        std::cout << x << " mon err " << mon_err << "\n";
    }
  }

  void raymarch_test_pixel(
      std::function<funcval_T(funcval_T, funcval_T, funcval_T)> f,
      glmvec3 eye, glmvec3 center, glm::ivec2 resolution,
      glm::ivec2 pixel_coordinate, float fovy_rad, glmvec2 interp_range,
      glmvec2 test_range, int march_resolution, int max_degree = 32,
      funcval_T interp_max_error = 5e-2) {
    float yedge = tanf(fovy_rad / 2);
    float xedge = static_cast<float>(resolution.x) / resolution.y * yedge;
    float yf =
        (static_cast<float>(pixel_coordinate.y) / (resolution.y - 1) * 2 - 1) *
        -1;
    float xf =
        static_cast<float>(pixel_coordinate.x) / (resolution.x - 1) * 2 - 1;

    yf *= yedge;
    xf *= xedge;

    glmmat4 invView =
        glm::inverse(glm::lookAt(eye, center, glmvec3(0.0f, 1.0f, 0.0f)));

    glmvec4 ahead = {xf, yf, -1, 1};
    glmvec3 dir = glmvec3(invView * ahead) - eye;

    raymarch_test_ray(f, eye, dir, interp_range, test_range, march_resolution,
                      max_degree, interp_max_error);
  }
};
