#pragma once

#include <functional>

template<typename funcval_T>
class FunctionApproximator {
protected:
	using funcRR = std::function<funcval_T(funcval_T)>;
	funcRR f;
	funcval_T x_begin, x_end;

public:
	FunctionApproximator(funcRR f, funcval_T x_begin, funcval_T x_end) :f(f), x_begin(x_begin), x_end(x_end) {}

	virtual void Fit() = 0;
	virtual funcval_T Evaluate(funcval_T x) const = 0;

	virtual void Show() const = 0;

};