import numpy as np
import scipy.stats as stats


def simulate_binomial_data_point(n, p, rng, grid_size):
    p = np.atleast_1d(p)

    data = []

    for p_i in p:
        x = stats.binom.rvs(n, p_i, random_state=rng)

        data.append(log_binomial_likelihood(n, x, grid_size=grid_size))

    data = np.atleast_2d(data)

    return data


def log_binomial_likelihood(n, x, eps=1e-10, grid_size=101):
    grid = np.linspace(0 + eps, 1 - eps, grid_size)

    return stats.binom.logpmf(x, n, grid)
