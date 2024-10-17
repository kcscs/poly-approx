#include "SegmentedChebyshevApproximator.h"
#include "monom.h"
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <Eigen/Dense>
#include <functional>
#include <string>
#include <iostream>

template <typename funcval_T> class Renderer {
public:

  using Vec3 = Eigen::Matrix<funcval_T, 3, 1>;

  void render(std::function<funcval_T(funcval_T, funcval_T, funcval_T)> f,
              Vec3 eye, Vec3 view, int width, int height, float fovy_rad, Vec3 light_dir,
              std::string outfile) {
    float yedge = tanf(fovy_rad);
    float xedge = static_cast<float>(width)/height * yedge;

    std::vector<unsigned char> pixels;
    pixels.reserve(width*height*3);

    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        std::cout<<x<<";"<<y<<"\n";
        float yf = static_cast<float>(y) / (height-1)*2-1;
        float xf = static_cast<float>(x) / (width-1)*2-1;

        yf *= yedge;
        xf *= xedge;

        Vec3 dir;
        dir(0) = xf;
        dir(1) = yf;
        dir(2) = 1;

        dir = dir-eye;


        std::function<funcval_T(funcval_T)> func = [&](funcval_T t) {
          Vec3 p = eye + t*dir;
          return f(p(0),p(1),p(2));
        };

        SegmentedChebyshevApproximator<funcval_T> cheb(func, 0, 20, 32, 5e-10);
        auto MonSegments = cheb.GetMonomSegments();
        std::vector<double> roots;
        for(const auto& s : MonSegments){
          auto r = s.FindRoots();
          roots.insert(roots.end(),r.begin(),r.end());
        }

        if(roots.size() == 0){
          pixels.push_back(0);
          pixels.push_back(0);
          pixels.push_back(255);
        } else {
          Vec3 wp = eye + roots[0] * dir;
          Vec3 norm;
          norm(0) = f(wp(0)+0.01,wp(1),wp(2)) - f(wp(0)-0.01,wp(1),wp(2));
          norm(1) = f(wp(0),wp(1)+0.01,wp(2)) - f(wp(0),wp(1)-0.01,wp(2));
          norm(2) = f(wp(0),wp(1),wp(2)+0.01) - f(wp(0),wp(1),wp(2)-0.01);

          if(norm(1) < 0)
            norm = -norm;

          float diffuse = norm.transpose() * (-light_dir);
          if(diffuse < 0)
            diffuse = 0;
          pixels.push_back(diffuse);
          pixels.push_back(diffuse);
          pixels.push_back(diffuse);
        }
      }
    }

    stbi_write_png(outfile.c_str(), width, height, 3, pixels.data(), width*3);
  }
};
