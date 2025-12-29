from __future__ import annotations
from scipy.special import gammaln as log_gamma, logsumexp, psi
import numpy as np
import click


def fit_pyclone_model(
    priors: Priors,
    var_params: VariationalParameters,
    data_preproc: DataPreprocessor,
    convergence_threshold=1e-6,
    max_iters=int(1e4),
    print_freq=100,
):

    epsilon = 1e-6
    elbo_trace = [compute_elbo(priors, var_params, data_preproc, epsilon)]

    for i in range(max_iters):
        if i % print_freq == 0:
            num_clusters = len(set(var_params.z.argmax(axis=1)))
            click.echo("Iteration: {}".format(i))
            click.echo("ELBO: {}".format(elbo_trace[-1]))
            click.echo("Number of clusters used: {}".format(num_clusters))
            click.echo()

        var_params.update_z(data_preproc)

        var_params.update_pi(priors)

        var_params.update_theta(priors, data_preproc)

        curr_elbo = compute_elbo(priors, var_params, data_preproc, epsilon)

        prev_elbo = elbo_trace[-1]

        elbo_trace.append(curr_elbo)

        diff = (curr_elbo - prev_elbo) / np.abs(curr_elbo)

        if diff < convergence_threshold:
            break

    return elbo_trace


class DataPreprocessor(object):
    __slots__ = "theta_update_data", "z_update_data", "theta_update_shape", "z_update_shape"

    def __init__(self, log_p_data: np.ndarray):
        self.theta_update_data = self._reshape_data_for_inference(log_p_data, [0, 1, 2])
        self.z_update_data = self._reshape_data_for_inference(log_p_data, [0, 2, 1])

        self.theta_update_shape = log_p_data.shape[1], log_p_data.shape[2]

        self.z_update_shape = log_p_data.shape[0]

    @staticmethod
    def _reshape_data_for_inference(log_p_data: np.ndarray, axis_order: list[int]) -> np.ndarray:
        new_axes_order = axis_order
        contraction_axis_size = log_p_data.shape[2] * log_p_data.shape[1]
        new_shape = [log_p_data.shape[0], contraction_axis_size]
        reshaped_data_arr = log_p_data.transpose(new_axes_order).reshape(new_shape)
        reshaped_data_arr.setflags(write=False)
        return reshaped_data_arr


class Priors(object):
    __slots__ = "pi", "theta", "log_theta", "pi_log_gamma"

    def __init__(self, num_clusters: int, num_grid_points: int, mix_weight_prior: float):
        self.pi = np.full(num_clusters, mix_weight_prior, dtype=np.float64, order="C")

        theta_fill_val = 1 / num_grid_points
        self.theta = np.full(num_grid_points, theta_fill_val, dtype=np.float64, order="C")

        self.log_theta = np.log(self.theta)

        self.pi_log_gamma = log_gamma(self.pi.sum()) - log_gamma(self.pi).sum()

        self._make_prior_arrays_read_only()

    def _make_prior_arrays_read_only(self):
        self.pi.setflags(write=False)
        self.theta.setflags(write=False)
        self.log_theta.setflags(write=False)


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

        self.theta = self._draw_initial_theta_value(num_clusters, num_dims, num_grid_points, rng)

        self.z = rng.dirichlet(ones_arr, size=num_data_points)

    @staticmethod
    def _draw_initial_theta_value(num_clusters: int, num_dims: int, num_grid_points: int, rng: np.random.Generator):
        pre_theta = rng.gamma(1, 1, size=(num_clusters, num_dims, num_grid_points))
        pre_theta /= pre_theta.sum(axis=2, keepdims=True)
        return pre_theta

    def update_pi(self, priors: Priors):
        self.pi = np.add(priors.pi, self.z.sum(axis=0), out=self.pi)

    def update_z(self, data_preproc: DataPreprocessor):
        new_z = get_log_p_data_theta(self.theta, data_preproc)

        psi_term = psi(self.pi)
        psi_term -= psi(self.pi.sum())

        new_z += psi_term

        new_z -= logsumexp(new_z, axis=1, keepdims=True)

        self.z = np.exp(new_z, order="C", out=self.z)

    def update_theta(self, priors: Priors, data_preproc: DataPreprocessor):

        reshaped_z = self.z.transpose([1, 0]).reshape(self.z.shape[1], self.z.shape[0])
        log_p_data_z = np.dot(reshaped_z, data_preproc.theta_update_data)
        log_p_data_z = log_p_data_z.reshape(self.z.shape[1], *data_preproc.theta_update_shape)

        log_p_data_z += priors.log_theta

        log_p_data_z -= logsumexp(log_p_data_z, axis=2, keepdims=True)
        self.theta = np.exp(log_p_data_z, order="C", out=self.theta)


def compute_elbo(priors: Priors, var_params: VariationalParameters, data_preproc: DataPreprocessor, epsilon: float):
    return compute_e_log_p(priors, var_params, data_preproc) - compute_e_log_q(var_params, epsilon)


def compute_e_log_p(priors: Priors, var_params: VariationalParameters, data_preproc: DataPreprocessor):
    log_p = priors.pi_log_gamma

    p_pi_z_term = priors.pi + var_params.z.sum(axis=0)
    p_pi_z_term -= 1

    pi_psi_term = psi(var_params.pi)
    pi_psi_term -= psi(var_params.pi.sum())
    pi_psi_term *= p_pi_z_term
    pi_psi_term = np.asarray(pi_psi_term)

    log_p += pi_psi_term.sum()

    log_p += (var_params.theta * priors.log_theta).sum()

    log_p_data_theta = get_log_p_data_theta(var_params.theta, data_preproc)

    log_p_data_theta *= var_params.z

    log_p += log_p_data_theta.sum()

    return log_p


def get_log_p_data_theta(theta: np.ndarray, data_preproc: DataPreprocessor):

    new_axes_order = [2, 1, 0]
    contraction_axis_size = theta.shape[2] * theta.shape[1]
    new_theta_shape = [contraction_axis_size, theta.shape[0]]
    reshaped_theta_arr = theta.transpose(new_axes_order).reshape(new_theta_shape)

    log_p_data_theta = np.dot(data_preproc.z_update_data, reshaped_theta_arr)
    log_p_data_theta = log_p_data_theta.reshape(data_preproc.z_update_shape, theta.shape[0])

    return log_p_data_theta


def compute_e_log_q(var_params: VariationalParameters, epsilon: float):
    log_p = 0.0

    pi_sum = var_params.pi.sum()

    log_p += log_gamma(pi_sum) - log_gamma(var_params.pi).sum()

    pi_psi_term = psi(var_params.pi)
    pi_psi_term -= psi(pi_sum)
    pi_psi_term *= var_params.pi - 1
    pi_psi_term = np.asarray(pi_psi_term)

    log_p += pi_psi_term.sum()

    theta_term = var_params.theta.clip(min=epsilon)
    theta_term = np.log(theta_term)
    theta_term *= var_params.theta

    log_p += theta_term.sum()

    z_term = var_params.z.clip(min=epsilon)
    z_term = np.log(z_term)
    z_term *= var_params.z

    log_p += z_term.sum()

    return log_p
