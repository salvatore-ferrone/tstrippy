import numpy as np

def riemann(f, a, b, n=100):
    """Approximate the integral of f(x) from a to b by the Riemann sum.
    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a , b : numbers
        Endpoints of the interval [a,b]
    n : integer
        Number of subintervals of equal length in the partition of [a,b]
    Returns
    -------
    float
        Approximation of the integral of f(x) from a to b using the
        Riemann sum with n subintervals of equal length.
    """
    return np.sum(f(np.linspace(a, b, n + 1)) * (b - a) / n)