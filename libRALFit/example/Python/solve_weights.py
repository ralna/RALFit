#!/usr/bin/python
#
# Attempts to fit the model y_i = x_1 e^(x_2 t_i)
# For parameters x_1 and x_2, and input data (t_i, y_i)

import numpy
import ral_nlls

# Calculate r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
def r(x, t, y):
    x1 = x[0]; x2 = x[1]

    return x1 * numpy.exp(x2 * t) - y

# Calculate:
# J_i1 = e^(x_2 * t_i)
# J_i2 = t_i x_1 e^(x_2 * t_i)
def J(x, t, y):
    x1 = x[0]; x2 = x[1]

    return numpy.column_stack((
        numpy.exp(x2*t),
        t * x1 * numpy.exp(x2*t)
        ))

# Calculate:
# Hr = sum_i r_i H_i
# Where H_i = [ 0                     t_i x_1 e^(x_2 t_i)    ]
#             [ t_i x_1 e^(x_2 t_i)   x_1 t_i^2 e^(x_2 t_i)  ]
def Hr(x, r, t, y):
    x1 = x[0]; x2 = x[1]

    Hr = numpy.zeros((2,2))
    Hr[0,0] = 0.0                       # H_11
    v = t * numpy.exp(x2*t)
    Hr[1,0] = numpy.dot(r, v)           # H_21
    Hr[1,1] = numpy.dot(r, (t*x1)*v)    # H_22

    return Hr

# Data to be fitted
t = numpy.array([1.0, 2.0, 4.0,  5.0,  8.0])
y = numpy.array([3.0, 4.0, 6.0, 11.0, 20.0])

# weights
w = numpy.array([1.0, 10.0, 5.0, 0.01, 0.5])

# Starting guess
x0 = numpy.array([2.5, 0.25])

# Call fitting routine
(x, inform) = ral_nlls.solve(x0, r, J, Hr=Hr, params=(t,y), weights=w)

# Print result
print("Found a local optimum at x = [ {0:.8f}  {1:.8f} ]".format(x[0],x[1]))
print("RALFit converged in {:d} iterations".format(inform["iter"]))
print("The cost function is {:.4f} at the minimum, and ||J^Tr||={:.2e}"
      .format(inform["obj"],inform["norm_g"])
)


