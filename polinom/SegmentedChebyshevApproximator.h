#pragma once

#include "FunctionApproximator.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <vector>
#include <matplot/matplot.h>
#include <iostream>
#include <stdexcept>

using namespace matplot;


using Eigen::MatrixXd;
using Eigen::Vector2;
using Eigen::VectorXd;



template<typename funcval_T>
class SegmentedChebyshevApproximator : FunctionApproximator<funcval_T> {
	using VectorT = Eigen::Matrix<funcval_T, Eigen::Dynamic, 1>;
private:
	struct ChebyshevSegment {
		VectorT coeffs;
		funcval_T begin, end;
	};
	std::vector<ChebyshevSegment> segments;

	int split_degree;
	double target_precision;

	funcval_T& x_begin = FunctionApproximator<funcval_T>::x_begin;
	funcval_T& x_end = FunctionApproximator<funcval_T>::x_end;
	FunctionApproximator<funcval_T>::funcRR& f = FunctionApproximator<funcval_T>::f;

	
	std::tuple<VectorT, VectorT, VectorT> chebyshev_interpolate(std::function<funcval_T(funcval_T)> func, int degree, funcval_T a, funcval_T b) {
		VectorXd x(degree + 1);
		for (int i = 0; i < degree + 1; ++i)
			x[i] = (b - a) / 2 * cos((pi * i) / degree) + (b + a) / 2;
		//std::cout << "Chebysev points:\n" << x << "\n\n";

		VectorXd vals(degree + 1);
		for (int i = 0; i < degree + 1; ++i)
			vals(i) = func(x(i));
		//std::cout << "Function values:\n" << vals << "\n\n";

		MatrixXd J(degree + 1, degree + 1);
		for (int j = 0; j < degree + 1; ++j) {
			for (int k = 0; k < degree + 1; ++k) {
				int pj = j == 0 || j == degree ? 2 : 1;
				int pk = k == 0 || k == degree ? 2 : 1;
				J(j, k) = 2.0 / (pj * pk * degree) * cos((j * pi * k) / degree);
			}
		}
		//std::cout << "J:\n" << J << "\n\n";


		VectorXd coeffs = J * vals;
		//std::cout << "Chebysev coefficients:\n" << coeffs << "\n\n";
		return std::make_tuple(coeffs, x, vals);
	}
	
	std::tuple<double, VectorT, VectorT> calc_interstitial_error(std::function<funcval_T(funcval_T)> func, VectorT coeffs, funcval_T a, funcval_T b) {
		int degree = coeffs.size() - 1;
		VectorT x(degree);
		for (int i = 1; i < 2 * degree + 1; i += 2)
			x(i / 2) = (b - a) / 2 * cos((pi * i) / (2 * degree)) + (b + a) / 2;

		VectorT y(degree);
		double interstitial_error = -1;
		for (int i = 0; i < x.size(); ++i) {
			y(i) = func(x(i));
			double approx = approx_cheb(coeffs, x(i), a, b);
			double err = abs(approx - y(i));
			if (err > interstitial_error)
				interstitial_error = err;
		}
		return std::make_tuple(interstitial_error, x, y);
	}
	
	funcval_T approx_cheb(const VectorT& coeffs, funcval_T x, funcval_T a, funcval_T b) const {
		funcval_T res = 0;
		for (int i = 0; i < coeffs.size(); ++i) {
			res += coeffs(i) * cos(i * acos((2 * x - (b + a)) / (b - a)));
		}
		return res;
	}

public:
	SegmentedChebyshevApproximator(FunctionApproximator<funcval_T>::funcRR f, funcval_T x_begin, funcval_T x_end, int split_degree = 64, double target_precision = 5e-10)
		: FunctionApproximator<funcval_T>(f, x_begin, x_end), split_degree(split_degree), target_precision(target_precision) {
		Fit();
	};

	virtual void Fit() override {

		funcval_T cur_begin = x_begin;
		funcval_T cur_end = x_end;
		segments.clear();

		double last_err = std::numeric_limits<double>::infinity();

		while (cur_begin + eps<funcval_T>() < x_end) {
			int cur_degree = 1;
			double err = std::numeric_limits<double>::infinity();
			VectorXd coeffs;
			double prev_err = std::numeric_limits<double>::infinity();
			bool first = true;
			bool split = false;
			int best_degree = -1;
			double best_error = std::numeric_limits<double>::infinity();
			while (err > target_precision) {
				prev_err = err;

				auto interp_res = chebyshev_interpolate(f, cur_degree, cur_begin, cur_end);
				coeffs = std::get<0>(interp_res);
				auto err_res = calc_interstitial_error(f, coeffs, cur_begin, cur_end);
				err = std::get<0>(err_res);
				std::cout <<"d:" << cur_degree << "  e:" << err << "\n";

				if (err < best_error) {
					best_error = err;
					best_degree = cur_degree;
				}

				cur_degree *= 2;
				if (cur_degree > split_degree) {
					split = true;
					break;
				}
			}

			if (split /*&& cur_end-cur_begin > (x_end-x_begin)/100*/) {
				cur_end = (cur_begin + cur_end) / 2;
				std::cout << "split: " << cur_begin << "-" << cur_end << " e:"<<err<< "\n";
			}
			else {
				auto interp_res = chebyshev_interpolate(f, best_degree, cur_begin, cur_end);
				coeffs = std::get<0>(interp_res);
				auto err_res = calc_interstitial_error(f, coeffs, cur_begin, cur_end);
				err = std::get<0>(err_res);
				std::cout << "segment: " << cur_begin << "-" << cur_end << " deg: " << best_degree << " err: " << err << "\n";
				segments.emplace_back();
				auto& seg = segments.back();
				seg.begin = cur_begin;
				seg.end = cur_end;
				seg.coeffs = coeffs;
				cur_begin = cur_end;
				cur_end = x_end;
			}
			//if (cur_degree >= split_degree || prev_err * 2 < err) {
			//	/*if (err > last_err + eps<double>()) {
			//		cur_end = cur_end * 2 - cur_begin;
			//		cur_degree = split_degree;
			//		std::cout << "segment (fail): " << cur_begin << "-" << cur_end << " deg: " << cur_degree << " err: " << err << "\n";
			//		auto interp_res = chebyshev_interpolate(f, cur_degree, cur_begin, cur_end);
			//		coeffs = std::get<0>(interp_res);
			//		segments.emplace_back();
			//		auto& seg = segments.back();
			//		seg.begin = cur_begin;
			//		seg.end = cur_end;
			//		seg.coeffs = coeffs;
			//		cur_begin = cur_end;
			//		cur_end = x_end;
			//		last_err = std::numeric_limits<double>::infinity();
			//	}*/
			//	//else {
			//		cur_end = (cur_begin + cur_end) / 2;
			//		last_err = err;
			//		std::cout << "split: " << cur_begin << "-" << cur_end << " e:"<<err<< "\n";
			//	//}
			//}
			//else {
			//	std::cout << "segment: " << cur_begin << "-" << cur_end << " deg: " << cur_degree << " err: " << err << "\n";
			//	segments.emplace_back();
			//	auto& seg = segments.back();
			//	seg.begin = cur_begin;
			//	seg.end = cur_end;
			//	seg.coeffs = coeffs;
			//	cur_begin = cur_end;
			//	cur_end = x_end;
			//}
		}
	}

	virtual funcval_T Evaluate(funcval_T x) const override {

		if (x < segments[0].begin)
			throw std::out_of_range("evaluated point is before first segment");
		int i = 0;
		while (i < segments.size() && segments[i].end < x)
			++i;
		if (i == segments.size())
			throw std::out_of_range("evaluated point is after last segment");

		return approx_cheb(segments[i].coeffs, x, segments[i].begin, segments[i].end);
	};

	virtual void Show() const {
		hold(true);
		std::vector<double> all_x;
		std::vector<double> correct_y;
		for (const ChebyshevSegment& seg : segments) {
			std::vector<double> x = linspace(seg.begin, seg.end);
			std::vector<double> y(x.size());
			for (int i = 0; i < y.size(); ++i) {
				y[i] = approx_cheb(seg.coeffs, x[i], seg.begin, seg.end);
				all_x.push_back(x[i]);
				correct_y.push_back(f(x[i]));
			}
			plot(x, y, "-");
		}
		plot(all_x, correct_y);
		show();
	}
};
