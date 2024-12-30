import numpy as np
import numpy.polynomial as poly

def power_eval(x, coeffs, dom, ran = [0,1]):
    p = poly.Polynomial(coeffs, dom)
    y = p(x)
    y = y*(ran[1]-ran[0])+ran[0]
    return y

def cheb_eval(x, coeffs, dom, ran = [0,1]):
    p = poly.Chebyshev(coeffs, dom)
    y = p(x)
    y = y*(ran[1]-ran[0])+ran[0]
    return y