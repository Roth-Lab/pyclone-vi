from __future__ import annotations
from scipy.special import gammaln as log_gamma, logsumexp, psi
# from numba import njit, prange
import numpy as np


def fit_annealed(
    log_p_data: np.ndarray,
    priors: Priors,
    var_params: VariationalParameters,
    annealing_power=1.0,
    convergence_threshold=1e-6,
    max_iters=int(1e4),
    num_annealing_steps=10,
    print_freq=100,
):
    if num_annealing_steps == 1:
        annealing_ladder = [1.0]

    else:
        annealing_ladder = np.linspace(0, 1.0, num_annealing_steps) ** annealing_power

    elbo_trace = []

    for t in annealing_ladder:
        print("Setting annealing factor to : {}".format(t))
        print()

        log_p_data_annealed = t * log_p_data

        if t == 1.0:
            convergence_threshold_t = convergence_threshold

        else:
            convergence_threshold_t = convergence_threshold * 1e-2

        elbo_trace = fit(
            log_p_data_annealed,
            priors,
            var_params,
            convergence_threshold=convergence_threshold_t,
            max_iters=max_iters,
            print_freq=print_freq,
        )

    return elbo_trace


def fit(
    log_p_data: np.ndarray,
    priors: Priors,
    var_params: VariationalParameters,
    convergence_threshold=1e-6,
    max_iters=int(1e4),
    print_freq=100,
):
    elbo_trace = [compute_elbo(log_p_data, priors, var_params)]

    for i in range(max_iters):
        if i % print_freq == 0:
            num_clusters = len(set(var_params.z.argmax(axis=1)))
            print("Iteration: {}".format(i))
            print("ELBO: {}".format(elbo_trace[-1]))
            print("Number of clusters used: {}".format(num_clusters))
            print()

        var_params.update_z(log_p_data)

        var_params.update_pi(priors)

        var_params.update_theta(log_p_data, priors)

        curr_elbo = compute_elbo(log_p_data, priors, var_params)

        prev_elbo = elbo_trace[-1]

        elbo_trace.append(curr_elbo)

        diff = (curr_elbo - prev_elbo) / np.abs(curr_elbo)

        if diff < convergence_threshold:
            break

    return elbo_trace


class Priors(object):
    __slots__ = "pi", "theta", "log_theta", "pi_log_gamma"

    def __init__(self, num_clusters: int, num_grid_points: int, mix_weight_prior: float):
        self.pi = np.full(num_clusters, mix_weight_prior, dtype=np.float64, order="C")

        theta_fill_val = 1 / num_grid_points
        self.theta = np.full(num_grid_points, theta_fill_val, dtype=np.float64, order="C")

        self.log_theta = np.log(self.theta)

        self.pi_log_gamma = log_gamma(self.pi.sum()) - log_gamma(self.pi).sum()


class VariationalParameters(object):
    __slots__ = "pi", "theta", "z"

    def __init__(
        self,
        num_clusters: int,
        num_data_points: int,
        num_dims: int,
        num_grid_points: int,
        rng: np.random.Generator,
    ):
        ones_arr = np.ones(num_clusters)

        self.pi = rng.dirichlet(ones_arr)

        pre_theta = rng.gamma(1, 1, size=(num_clusters, num_dims, num_grid_points))
        pre_theta /= pre_theta.sum(axis=2, keepdims=True)

        self.theta = pre_theta

        self.z = rng.dirichlet(ones_arr, size=num_data_points)

    def update_pi(self, priors: Priors):
        self.pi = np.add(priors.pi, self.z.sum(axis=0), out=self.pi)

    def update_z(self, log_p_data):
        new_z = get_log_p_data_theta(log_p_data, self.theta)

        psi_term = psi(self.pi)
        psi_term -= psi(self.pi.sum())

        new_z += psi_term

        new_z -= logsumexp(new_z, axis=1, keepdims=True)

        self.z = np.exp(new_z, order="C", out=self.z)

    def update_theta(self, log_p_data, priors: Priors):

        # log_p_data_z = np.zeros((self.z.shape[1], log_p_data.shape[1], log_p_data.shape[2]), order="C")
        # compute_log_p_data_z(log_p_data, self.z, log_p_data_z)

        log_p_data_z = np.tensordot(self.z, log_p_data, axes=([0], [0]))

        log_p_data_z += priors.log_theta

        log_p_data_z -= logsumexp(log_p_data_z, axis=2, keepdims=True)
        self.theta = np.exp(log_p_data_z, order="C", out=self.theta)


def compute_elbo(log_p_data, priors: Priors, var_params: VariationalParameters):
    return compute_e_log_p(log_p_data, priors, var_params) - compute_e_log_q(var_params)


def compute_e_log_p(log_p_data, priors: Priors, var_params: VariationalParameters):
    log_p = priors.pi_log_gamma

    p_pi_z_term = priors.pi + var_params.z.sum(axis=0)
    p_pi_z_term -= 1

    pi_psi_term = psi(var_params.pi)
    pi_psi_term -= psi(var_params.pi.sum())
    pi_psi_term *= p_pi_z_term
    pi_psi_term = np.asarray(pi_psi_term)

    log_p += pi_psi_term.sum()

    log_p += (var_params.theta * priors.log_theta).sum()

    log_p_data_theta = get_log_p_data_theta(log_p_data, var_params.theta)

    log_p_data_theta *= var_params.z

    log_p += log_p_data_theta.sum()

    return log_p


def get_log_p_data_theta(log_p_data, theta):
    # log_p_data_theta = np.zeros((log_p_data.shape[0], theta.shape[0]), order="C")
    # compute_log_p_data_theta(log_p_data, theta, log_p_data_theta)
    log_p_data_theta = np.tensordot(log_p_data, theta, axes=([2, 1], [2, 1]))

    return log_p_data_theta


def compute_e_log_q(var_params: VariationalParameters):
    log_p = 0.0

    pi_sum = var_params.pi.sum()

    log_p += log_gamma(pi_sum) - log_gamma(var_params.pi).sum()

    pi_psi_term = psi(var_params.pi)
    pi_psi_term -= psi(pi_sum)
    pi_psi_term *= (var_params.pi - 1)
    pi_psi_term = np.asarray(pi_psi_term)

    log_p += pi_psi_term.sum()

    theta_term = np.log(var_params.theta + 1e-6)
    theta_term *= var_params.theta

    log_p += theta_term.sum()

    z_term = np.log(var_params.z + 1e-6)
    z_term *= var_params.z

    log_p += z_term.sum()

    return log_p


# @njit(parallel=True)
# def compute_log_p_data_z(log_p_data, z, result):
#     """Equivalent to np.sum(var_params.z[:, :, np.newaxis, np.newaxis] * log_p_data[:, np.newaxis, :, :], axis=0)"""
#     N, D, G = log_p_data.shape
#
#     K = z.shape[1]
#
#     for cluster in prange(K):
#         for mut in range(N):
#             for sample in range(D):
#                 for grid_point in range(G):
#                     result[cluster, sample, grid_point] += log_p_data[mut, sample, grid_point] * z[mut, cluster]
#
#
# @njit(parallel=True, fastmath=True)
# def compute_log_p_data_theta(log_p_data, theta, result):
#     """Equivalent to np.sum(var_params.theta[np.newaxis, :, :, :] * log_p_data[:, np.newaxis, :, :], axis=(2, 3))"""
#     N, D, G = log_p_data.shape
#
#     K = theta.shape[0]
#
#     for mut in prange(N):
#         for cluster in range(K):
#             for sample in range(D):
#                 for grid_point in range(G):
#                     result[mut, cluster] += log_p_data[mut, sample, grid_point] * theta[cluster, sample, grid_point]
