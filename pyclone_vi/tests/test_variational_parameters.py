import unittest
from threadpoolctl import threadpool_limits
import numpy as np
from scipy.special import logsumexp as log_sum_exp, psi
from pyclone_vi.inference import VariationalParameters, Priors, DataPreprocessor
from pyclone_vi.tests.simulate import simulate_binomial_data_point
from numba import njit, prange, set_num_threads


class OldVariationalParameters(object):
    def __init__(self, pi, theta, z):
        self.pi = pi

        self.theta = theta

        self.z = z


@njit(parallel=True)
def compute_log_p_data_z(log_p_data, z):
    """Equivalent to np.sum(var_params.z[:, :, np.newaxis, np.newaxis] * log_p_data[:, np.newaxis, :, :], axis=0)"""
    N, D, G = log_p_data.shape

    K = z.shape[1]

    result = np.zeros((K, D, G))

    for k in prange(K):
        for d in range(D):
            for g in range(G):
                for n in range(N):
                    result[k, d, g] += log_p_data[n, d, g] * z[n, k]

    return result


@njit(parallel=True)
def compute_log_p_data_theta(log_p_data, theta):
    """Equivalent to np.sum(var_params.theta[np.newaxis, :, :, :] * log_p_data[:, np.newaxis, :, :], axis=(2, 3))"""
    N, D, G = log_p_data.shape

    K = theta.shape[0]

    result = np.zeros((N, K))

    for n in prange(N):
        for k in range(K):
            for d in range(D):
                for g in range(G):
                    result[n, k] += log_p_data[n, d, g] * theta[k, d, g]

    return result


def update_pi(priors, var_params):
    var_params.pi = priors.pi + np.sum(var_params.z, axis=0)


def update_z(log_p_data, var_params):
    var_params.z = compute_log_p_data_theta(log_p_data, var_params.theta)

    var_params.z += (psi(var_params.pi) - psi(np.sum(var_params.pi)))[np.newaxis, :]

    var_params.z = var_params.z - log_sum_exp(var_params.z, axis=1)[:, np.newaxis]

    var_params.z = np.exp(var_params.z)


def update_theta(log_p_data, priors, var_params):
    var_params.theta = np.log(priors.theta[np.newaxis, np.newaxis, :]) + compute_log_p_data_z(log_p_data, var_params.z)

    var_params.theta = var_params.theta - log_sum_exp(var_params.theta, axis=2)[:, :, np.newaxis]

    var_params.theta = np.exp(var_params.theta)


class TestVariationalParameterUpdates(unittest.TestCase):

    def setUp(self) -> None:
        self.default_grid_size = 100
        self.rng_seed = 242643578967193853558243570818064774262
        self.num_threads = 10
        set_num_threads(self.num_threads)
        self.rng = np.random.default_rng(self.rng_seed)

    def create_log_p_data(self, depth, num_data_points, num_dims, num_grid_points):
        generator_exp = (
            simulate_binomial_data_point(depth, self.rng.random(num_dims), self.rng, num_grid_points) for _ in range(num_data_points)
        )
        log_p_data = np.fromiter(
            generator_exp,
            dtype=np.dtype((np.float64, (num_dims, num_grid_points))),
            count=num_data_points,
        )
        return log_p_data

    def create_test_structures_big(self):
        num_clusters = 40
        num_data_points = 1000
        num_dims = 100
        num_grid_points = self.default_grid_size
        depth = 100

        expected_var_params, log_p_data, priors, actual_var_params = self.create_test_structures(
            depth,
            num_clusters,
            num_data_points,
            num_dims,
            num_grid_points,
        )

        return log_p_data, actual_var_params, priors, expected_var_params

    def create_test_structures_small(self):
        num_clusters = 10
        num_data_points = 10
        num_dims = 10
        num_grid_points = self.default_grid_size
        depth = 100

        expected_var_params, log_p_data, priors, actual_var_params = self.create_test_structures(
            depth,
            num_clusters,
            num_data_points,
            num_dims,
            num_grid_points,
        )

        return log_p_data, actual_var_params, priors, expected_var_params

    def create_test_structures(self, depth, num_clusters, num_data_points, num_dims, num_grid_points):
        log_p_data = self.create_log_p_data(depth, num_data_points, num_dims, num_grid_points)
        priors = Priors(num_clusters, num_grid_points, 1.0)
        actual_var_params = VariationalParameters(num_clusters, num_data_points, num_dims, num_grid_points, self.rng)
        expected_var_params = OldVariationalParameters(
            actual_var_params.pi.copy(),
            actual_var_params.theta.copy(),
            actual_var_params.z.copy(),
        )
        return expected_var_params, log_p_data, priors, actual_var_params

    def test_update_pi_small(self):

        log_p_data, actual_var_params, priors, expected_var_params = self.create_test_structures_small()

        update_pi(priors, expected_var_params)

        with threadpool_limits(limits=self.num_threads, user_api="blas"):
            actual_var_params.update_pi(priors)

        np.testing.assert_allclose(actual_var_params.pi, expected_var_params.pi)

    def test_update_z_small(self):
        log_p_data, actual_var_params, priors, expected_var_params = self.create_test_structures_small()

        data_preproc = DataPreprocessor(log_p_data)

        update_z(log_p_data, expected_var_params)

        with threadpool_limits(limits=self.num_threads, user_api="blas"):
            actual_var_params.update_z(data_preproc)

        np.testing.assert_allclose(actual_var_params.z, expected_var_params.z)

    def test_update_theta_small(self):
        log_p_data, actual_var_params, priors, expected_var_params = self.create_test_structures_small()

        data_preproc = DataPreprocessor(log_p_data)

        update_theta(log_p_data, priors, expected_var_params)

        with threadpool_limits(limits=self.num_threads, user_api="blas"):
            actual_var_params.update_theta(priors, data_preproc)

        np.testing.assert_allclose(actual_var_params.theta, expected_var_params.theta)


if __name__ == "__main__":
    unittest.main()
