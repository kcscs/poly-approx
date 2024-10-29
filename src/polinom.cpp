// polinom.cpp : This file contains the 'main' function. Program execution
// begins and ends there.
//

#include "SegmentedChebyshevApproximator.h"
#include "render.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <ios>
#include <iostream>
#include <limits>
#include <math.h>
#include <matplot/matplot.h>
#include <random>
#include <set>
// #include <matplot/backend/opengl.h>

using Eigen::MatrixXd;
using Eigen::Vector2;
using Eigen::VectorXd;
using namespace matplot;

constexpr int SAMPLE_COUNT = 12;
constexpr int DEGREE = 10;
constexpr int ESTIMATION_DEGREE = 3;
constexpr double RANGE_MIN = -20;
constexpr double RANGE_MAX = 20;
constexpr int MAX_DEGREE =
    DEGREE > ESTIMATION_DEGREE ? DEGREE : ESTIMATION_DEGREE;
// constexpr double TARGET_PRECISION = 5e-16;
constexpr double TARGET_PRECISION = 5e-10;
constexpr int SPLIT_DEGREE = 50;

// consteval size_t choose(size_t n, size_t k) {
//	if (k == 0)
//		return 1;
//	if (n == k)
//		return 1;
//	return choose(n - 1, k) + choose(n-1, k - 1);
// }

struct nCrTable {
  size_t table[MAX_DEGREE][MAX_DEGREE] = {0};

  consteval nCrTable() {
    table[0][0] = 1;
    for (size_t i = 1; i < MAX_DEGREE; ++i) {
      table[i][0] = table[i][i] = 1;
      for (size_t j = 1; j < i; ++j) {
        table[i][j] = table[i - 1][j] + table[i - 1][j - 1];
      }
    }
  }
};

constexpr nCrTable NCHOOSEK;

double eval_monomial(const std::vector<double> &coeffs, double x) {
  double cur_base = 1;
  double res = 0;
  for (auto c : coeffs) {
    res += c * cur_base;
    cur_base *= x;
  }
  return res;
}

MatrixXd construct_matrix_monomial(const std::vector<double> &points_x,
                                   size_t degree) {
  MatrixXd mat(points_x.size(), degree + 1);
  for (size_t i = 0; i < points_x.size(); ++i) {
    double cur_base = 1;
    for (size_t j = 0; j <= degree; ++j) {
      mat(i, j) = cur_base;
      cur_base *= points_x[i];
    }
  }
  return mat;
}

double eval_cheb(const std::vector<double> &coeffs, double x, int kind = 1) {
  double cur_base = kind == 1 ? x : 2 * x;
  double prev_base = 1;
  double res = 1 * coeffs[0];
  for (size_t i = 1; i < coeffs.size(); ++i) {
    res += coeffs[i] * cur_base;
    double tmp = cur_base;
    cur_base = 2 * x * cur_base - prev_base;
    prev_base = tmp;
  }
  return res;
}

MatrixXd construct_matrix_cheb(const std::vector<double> &points_x,
                               size_t degree, int kind = 1) {
  MatrixXd mat(points_x.size(), degree + 1);
  for (size_t i = 0; i < points_x.size(); ++i) {
    mat(i, 0) = 1;
    mat(i, 1) = kind == 1 ? points_x[i] : 2 * points_x[i];
    for (size_t j = 2; j <= degree; ++j) {
      mat(i, j) = 2 * points_x[i] * mat(i, j - 1) - mat(i, j - 2);
    }
  }
  return mat;
}

double eval_bernstein(const std::vector<double> &coeffs, double x) {
  double res = 0;
  double xp = 1;
  double mxp = std::pow(1 - x, coeffs.size() - 1);
  for (size_t i = 0; i < coeffs.size(); ++i) {
    res += NCHOOSEK.table[coeffs.size() - 1][i] * xp * mxp;
    xp *= x;
    mxp /= (1 - x);
  }
  return res;
}

MatrixXd construct_matrix_bernstein(const std::vector<double> &points_x,
                                    size_t degree) {
  MatrixXd mat(points_x.size(), degree + 1);
  for (size_t i = 0; i < points_x.size(); ++i) {
    double xp = 1;
    double mxp = std::pow(1 - points_x[i], degree);
    for (size_t j = 0; j <= degree; ++j) {
      mat(i, j) = NCHOOSEK.table[degree][j] * xp * mxp;
      xp *= points_x[i];
      mxp /= (1 - points_x[i]);
    }
  }
  return mat;
}

std::tuple<VectorXd, VectorXd, VectorXd>
chebyshev_interpolate(std::function<double(double)> func, int degree, double a,
                      double b) {
  VectorXd x(degree + 1);
  for (int i = 0; i < degree + 1; ++i)
    x[i] = (b - a) / 2 * cos((pi * i) / degree) + (b + a) / 2;
  // std::cout << "Chebysev points:\n" << x << "\n\n";

  VectorXd vals(degree + 1);
  for (int i = 0; i < degree + 1; ++i)
    vals(i) = func(x(i));
  // std::cout << "Function values:\n" << vals << "\n\n";

  MatrixXd J(degree + 1, degree + 1);
  for (int j = 0; j < degree + 1; ++j) {
    for (int k = 0; k < degree + 1; ++k) {
      int pj = j == 0 || j == degree ? 2 : 1;
      int pk = k == 0 || k == degree ? 2 : 1;
      J(j, k) = 2.0 / (pj * pk * degree) * cos((j * pi * k) / degree);
    }
  }
  // std::cout << "J:\n" << J << "\n\n";

  VectorXd coeffs = J * vals;
  // std::cout << "Chebysev coefficients:\n" << coeffs << "\n\n";
  return std::make_tuple(coeffs, x, vals);
}

double approx_cheb(const VectorXd &coeffs, double x, double a, double b) {
  double res = 0;
  for (int i = 0; i < coeffs.size(); ++i) {
    res += coeffs(i) * cos(i * acos((2 * x - (b + a)) / (b - a)));
  }
  return res;
}

std::tuple<double, VectorXd, VectorXd>
calc_interstitial_error(std::function<double(double)> func, VectorXd coeffs,
                        double a, double b) {
  int degree = coeffs.size() - 1;
  VectorXd x(degree);
  for (int i = 1; i < 2 * degree + 1; i += 2)
    x(i / 2) = (b - a) / 2 * cos((pi * i) / (2 * degree)) + (b + a) / 2;

  VectorXd y(degree);
  double interstitial_error = -1;
  for (int i = 0; i < x.size(); ++i) {
    y(i) = func(x(i));
    double approx = approx_cheb(coeffs, x(i), RANGE_MIN, RANGE_MAX);
    double err = abs(approx - y(i));
    if (err > interstitial_error)
      interstitial_error = err;
  }
  return std::make_tuple(interstitial_error, x, y);
}

// using ChebyshevSegment = ChebyshevSegment<double>;
// struct ChebyshevSegment {
// 	VectorXd coeffs;
// 	double begin, end;
// };

std::vector<ChebyshevSegment<double>>
approximate_adaptive_chebyshev(std::function<double(double)> func, double begin,
                               double end) {
  std::vector<ChebyshevSegment<double>> segments;

  double cur_begin = begin;
  double cur_end = end;
  /*segments.emplace_back();
  auto& seg = segments.back();
  seg.begin = cur_begin;
  seg.end = cur_end;*/

  while (cur_begin + eps<double>() < end) {
    int cur_degree = 1;
    double err = std::numeric_limits<double>::infinity();
    VectorXd coeffs;
    while (err > TARGET_PRECISION && cur_degree < SPLIT_DEGREE) {
      cur_degree *= 2;
      auto interp_res =
          chebyshev_interpolate(func, cur_degree, cur_begin, cur_end);
      coeffs = std::get<0>(interp_res);
      auto err_res = calc_interstitial_error(func, coeffs, cur_begin, cur_end);
      err = std::get<0>(err_res);
      std::cout << "e:" << err << "\n";
    }
    if (cur_degree >= SPLIT_DEGREE) {
      cur_end = (cur_begin + cur_end) / 2;
    } else {
      std::cout << "segment: " << cur_begin << "-" << cur_end
                << " deg: " << cur_degree << " err: " << err << "\n";
      segments.emplace_back();
      auto &seg = segments.back();
      seg.begin = cur_begin;
      seg.end = cur_end;
      seg.coeffs = coeffs;
      cur_begin = cur_end;
      cur_end = end;
    }
  }

  return segments;
}
/*
int main()
{
        MatrixXd m(2, 2);
        m(0, 0) = 3;
        m(1, 0) = 2.5;
        m(0, 1) = -1;
        m(1, 1) = m(1, 0) + m(0, 1);
        std::cout << m << std::endl;

        std::default_random_engine rng(std::time(nullptr));degree
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        auto rand = std::bind(unif, rng);

        std::vector<double> coeffs(DEGREE+1);
        for (auto& c : coeffs)
                c = rand()*2-1;

        std::vector<double> x = linspace(RANGE_MIN, RANGE_MAX);
        std::vector<double> y = transform(x, [&](auto x) {return
eval_monomial(coeffs, x); });


        std::vector<double> points_x = linspace(RANGE_MIN, RANGE_MAX,
SAMPLE_COUNT); std::vector<double> points_y = transform(points_x, [&](auto x)
{return eval_monomial(coeffs, x); });

         //LSQ interpolation in different bases
        //MatrixXd A = construct_matrix_monomial(points_x, ESTIMATION_DEGREE+1);
        MatrixXd A = construct_matrix_cheb(points_x, ESTIMATION_DEGREE + 1, 2);
        //MatrixXd A = construct_matrix_bernstein(points_x, ESTIMATION_DEGREE +
1); MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();
        //MatrixXd pinv = (A.transpose() * A).inverse() * A.transpose();
        VectorXd b(points_y.size());
        for (size_t i = 0; i < points_y.size(); ++i)
                b(i,0) = points_y[i];
        VectorXd res_lsq = pinv * b;
        std::cout << "Coeffs lsq:\n" << res_lsq << "\n\n";

        // Chebyshev interpolation
        auto f = [&](double x) {return eval_monomial(coeffs, x); };
        auto cheb_res = chebyshev_interpolate(f, ESTIMATION_DEGREE, RANGE_MIN,
RANGE_MAX); VectorXd res_cheb = std::get<0>(cheb_res); for (int i = 0; i <
ESTIMATION_DEGREE + 1; ++i) { auto& cheb_points = std::get<1>(cheb_res); auto&
cheb_vals = std::get<2>(cheb_res); points_x[i] = cheb_points[i]; points_y[i] =
cheb_vals[i];
        }
        points_x.resize(ESTIMATION_DEGREE + 1);
        points_y.resize(ESTIMATION_DEGREE + 1);

        std::vector<double> coeffs_est(ESTIMATION_DEGREE+1);
        for (size_t i = 0; i < coeffs_est.size(); ++i) {
                coeffs_est[i] = res_cheb(i);
        }

        //std::vector<double> y_est = transform(x, [&](auto x) {return
eval_monomial(coeffs_est, x); });
        //std::vector<double> y_est = transform(x, [&](auto x) {return
eval_cheb(coeffs_est, x, 2); });
        //std::vector<double> y_est = transform(x, [&](auto x) {return
eval_bernstein(coeffs_est, x); });

        std::vector<double> y_est = transform(x, [&](auto x) {return
approx_cheb(res_cheb, x, RANGE_MIN, RANGE_MAX); });


        auto err_res = calc_interstitial_error(f, res_cheb, RANGE_MIN,
RANGE_MAX); std::cout << "Interstitial error: " << std::get<0>(err_res) <<
"\n\n";

        std::vector<double> inter_x(std::get<1>(err_res).size());
        std::vector<double> inter_y(std::get<1>(err_res).size());
        for (int i = 0; i < inter_x.size(); ++i) {
                const auto& inter_points = std::get<1>(err_res);
                const auto& inter_vals = std::get<2>(err_res);
                inter_x[i] = inter_points(i);
                inter_y[i] = inter_vals(i);
        }

        plot(x, y, "--", x, y_est, points_x, points_y, "x", inter_x, inter_y,
"o");
        //std::set<std::vector<double>> Y = {
        //	{16, 5, 9, 4}, {2, 11, 7, 14}, {3, 10, 6, 15}, {13, 8, 12, 1} };
        //plot(Y);

        show();



        auto segments = approximate_adaptive_chebyshev(f, RANGE_MIN, RANGE_MAX);

}*/

template <typename T> T wilkinson(T x) {
  T res = 1;
  for (int i = 1; i <= 20; ++i)
    res *= (x - i);
  return res / 1000.0f;
}

template <typename T, T a, T b> T linear(T x) { return a * x + b; }

double sind2(double d) { return sin(d); }

int main() {
  // std::ios_base::sync_with_stdio(false);

  std::default_random_engine rng(std::time(nullptr));
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  auto rand = std::bind(unif, rng);

  std::vector<double> coeffs(50 + 1);
  std::cout << "coeffs: ";
  for (auto &c : coeffs) {
    c = rand() * 2 - 1;
    std::cout << ", " << c;
  }
  std::cout << "\n";

  // figure(); //<backend::opengl>

  // SegmentedChebyshevApproximator<double> cheb(wilkinson<double>, 0.0, 1.0,
  // 32, 5e-8); SegmentedChebyshevApproximator<double> cheb(sind2, 0.0, 10.0,
  // 32, 5e-16);

  // SegmentedChebyshevApproximator<double> cheb(
  //     [&](double x) { return eval_monomial(coeffs, x); }, -1, 1, 32, 1e-12);
  // // SegmentedChebyshevApproximator<double> cheb(linear<double, 1.0, 0.0>,
  // // 0.0, 5.0, 128, 5e-10);
  // // cheb.Show();
  // figure();
  // hold(on);
  // auto MonSegments = cheb.GetMonomSegments();
  // for (const auto &s : MonSegments) {
  //   s.Show();
  // }

  // std::vector<double> roots;
  // for (const auto &s : MonSegments) {
  //   auto rs = s.FindRoots();
  //   roots.insert(roots.end(), rs.begin(), rs.end());
  // }
  // for (auto r : roots) {
  //   std::cout << "root: " << r << "\n";
  // }

  // std::vector<double> rootzeros(roots.size());
  // auto p = plot(roots, rootzeros, "x");
  // p->line_width(3);
  // show();

  // ChebyshevSegment<double> root_test_segment({1,1.6,0.5},0,1);
  // // auto roots = FindQuadraticChebyshevRoot(root_test_segment);

  // figure();
  // hold(on);
  // root_test_segment.Show();
  // figure();
  // root_test_segment.ConvertToMonomial().Show();
  // show();

  // std::vector<double> rootx = {roots.first,roots.second};
  // std::vector<double> rooty = {0,0};
  // plot(rootx,rooty,"x");
  // show();
  // std::cout<<"Roots: "<<rootx[0]<<"("<<root_test_segment.Evaluate(rootx[0])
  // << "), "<<rootx[1]<<"("<<root_test_segment.Evaluate(rootx[1])<<")\n";

  // Renderer
  Renderer<float> renderer;
  glm::dvec3 eye = {0, 0, -5};
  glm::dvec3 at = {0, 0, 0};
  glm::dvec3 light = glm::normalize(glm::dvec3(0.5, -1, -0.2));

  // Sphere
  // renderer.render(
  //   [](double x, double y, double z) { return sqrt(x*x+y*y+z*z)-3; },
  //     eye, at, 64, 64, 3.14/2,
  //     light, "Test6.png");

  // Box
  // renderer.render(
  //   [](double x, double y, double z) {
  //     glm::dvec3 p(x,y,z);
  //     glm::dvec3 b(1,0.5,1.5);
  //     glm::dvec3 q = glm::abs(p)-b;
  //     return glm::length(glm::max(q,0.0)) +
  //     glm::min(glm::max(q.x,glm::max(q.y,q.z)),0.0);
  //   },
  //     glm::dvec3(2,3,-4), at, 64, 64, 3.14/2,
  // glm::normalize(glm::dvec3(1.0f,-3.0f,1.0f)), "Test15.png");

  // Torus
  // renderer.render(
  //     [](double x, double y, double z) {
  //       glm::dvec3 p(x, y, z);
  //       glm::dvec2 t(2, 0.5);
  //       glm::dvec2 q = glm::dvec2(glm::length(glm::dvec2(p.x, p.z)) - t.x, p.y);
  //       return glm::length(q) - t.y;
  //     },
  //     glm::dvec3(2, 3, -4), at, 16, 16, 3.14 / 2,
  //     glm::normalize(glm::dvec3(1.0f, -3.0f, 1.0f)), "Test17.png");
  // test 11: max deg 128 (64), max err 1e-13
  // test 12: max deg 128 (actually 128), max err 1e-13
  // test 13: max deg 16, max err 1e-13
  // test 14: max deg 8, max err 5e-15

  // Torus stable (degree 4)
  // renderer.render(
  //     [](double x, double y, double z) {
  //       double R = 2.0;
  //       double r = 0.5;
  //       double lhs = 4*R*R*(x*x+y*y);
  //       double rhs = (x*x+y*y+z*z+R*R-r*r);
  //       return rhs*rhs-lhs;
  //     },
  //     glm::dvec3(2, 3, -4), at, 64, 64, 3.14 / 2,
  //     glm::normalize(glm::dvec3(1.0f, -3.0f, 1.0f)), "Test16.png");


  renderer.raymarch_test_ray(
      [](double x, double y, double z) {
        glm::dvec3 p(x, y, z);
        glm::dvec2 t(2, 0.5);
        glm::dvec2 q = glm::dvec2(glm::length(glm::dvec2(p.x, p.z)) - t.x, p.y);
        return glm::length(q) - t.y;
      },
      glm::dvec3(2, 3, -4), glm::dvec3(0.0, 2.0, 0.0) - glm::normalize(glm::dvec3(2, 3, -4)),
      {0.0, 20.0}, {4.0, 8.0}, 1000, 8, 1e-10);

  // renderer.render(
  //   [](double x, double y, double z) { return sin(x)-y*y*y-z*z; },
  //     glm::dvec3(-1,3,-4), at, 64, 64, 3.14/2,
  //     glm::normalize(glm::dvec3(1.0f,-3.0f,1.0f)), "Test7.png");
}
