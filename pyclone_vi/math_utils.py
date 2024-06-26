from numba import njit
import numpy as np
from math import lgamma


# Special functions


@njit(cache=True)
def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)


@njit(cache=True)
def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(n - x) - log_factorial(x)


@njit(cache=True)
def log_factorial(x):
    return log_gamma(x + 1)


@njit(cache=True)
def log_gamma(x):
    return lgamma(x)


@njit(cache=True)
def log_normalize(x):
    return x - log_sum_exp(x)


@njit(cache=True)
def log_sum_exp(log_X):
    max_exp = np.max(log_X)

    if np.isinf(max_exp):
        return max_exp

    total = 0

    for x in log_X:
        total += np.exp(x - max_exp)

    return np.log(total) + max_exp


# Probability densities


@njit(cache=True)
def log_beta_binomial_pdf(n, x, a, b):
    return log_binomial_coefficient(n, x) + log_beta(a + x, b + n - x) - log_beta(a, b)


@njit(cache=True)
def log_binomial_pdf(n, x, p):
    return log_binomial_coefficient(n, x) + x * np.log(p) + (n - x) * np.log1p(-p)
