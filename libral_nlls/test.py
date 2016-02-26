#!/usr/bin/python
#
# Attempts to fit the model y_i = x_1 e^(x_2 t_i)
# For parameters x_1 and x_2, and input data (t_i, y_i)

import numpy
import ral_nlls

# Calculate r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
def r(x, data):
    x1 = x[0]; x2 = x[1]
    t = data[:,0]; y = data[:,1]

    return x1 * numpy.exp(x2 * t) - y

# Calculate:
# J_i1 = e^(x_2 * t_i)
# J_i2 = t_i x_1 e^(x_2 * t_i)
def J(x, data):
    x1 = x[0] ;x2 = x[1]
    t = data[:,0]; y = data[:,1]

    return [
        numpy.exp(x2*t),
        t * x1 * numpy.exp(x2*t)
        ]

# Calculate:
# HF = sum_i r_i H_i
# Where H_i = [ 1                t_i e^(x_2 t_i)    ]
#             [ t_i e^(x_2 t_i)  t_i^2 e^(x_2 t_i)  ]
def HF(x, r, data):
    x1 = x[1]
    x2 = x[2]
    HF = numpy.zeros((2,2))
    for (r, (t, y)) in zip(r,data):
        HF[1,1] += r                    # H_11
        HF[2,1] += r * t * exp(x2*t)    # H_21
        HF[2,2] += r * t*t * exp(x2*t)  # H_22

    return HF

# Data to be fitted
data = numpy.array([
    [1.0,  3.0],
    [2.0,  4.0],
    [4.0,  6.0],
    [5.0, 11.0],
    [8.0, 20.0]
    ])

# Starting guess
x0 = numpy.array([2.5, 0.25])

# Call fitting routine
x = ral_nlls.solve(x0, r, J, params=data)

# Print result
print "Found a local optimum at x = ", x
