import unittest
import numpy as np
from scipy.special import gammaln as log_gamma, psi
from pyclone_vi.inference import compute_e_log_q, VariationalParameters, compute_e_log_p, Priors, DataPreprocessor
from pyclone_vi.tests.simulate import simulate_binomial_data_point
from numba import njit, prange, set_num_threads


def compute_e_log_q_old(var_params):
    log_p = 0

    log_p += log_gamma(np.sum(var_params.pi)) - np.sum(log_gamma(var_params.pi))

    log_p += np.sum((var_params.pi - 1) * (psi(var_params.pi) - psi(np.sum(var_params.pi))))

    log_p += np.sum(var_params.theta * np.log(var_params.theta + 1e-6))

    log_p += np.sum(var_params.z * np.log(var_params.z + 1e-6))

    return log_p


def compute_e_log_p_old(log_p_data, priors, var_params):
    log_p = 0

    log_p += log_gamma(np.sum(priors.pi)) - np.sum(log_gamma(priors.pi))

    log_p += np.sum((priors.pi + np.sum(var_params.z, axis=0) - 1) * (psi(var_params.pi) - psi(np.sum(var_params.pi))))

    log_p += np.sum(var_params.theta * np.log(priors.theta)[np.newaxis, np.newaxis, :])

    log_p += np.sum(var_params.z * compute_log_p_data_theta(log_p_data, var_params.theta))

    return log_p


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


class TestComputeELogQ(unittest.TestCase):
    def __init__(self, method_name: str = ...):
        super().__init__(method_name)

        self.default_grid_size = 100

        self.rng_seed = 242643578967193853558243570818064774262

        self.rng = None

    def setUp(self) -> None:
        self.rng = np.random.default_rng(self.rng_seed)

    def run_test(self, num_clusters, num_data_points, num_dims, num_grid_points):
        var_params = VariationalParameters(num_clusters, num_data_points, num_dims, num_grid_points, self.rng)
        expected = compute_e_log_q_old(var_params)
        actual = compute_e_log_q(var_params)
        print("Expected = ", expected)
        print("Actual = ", actual)
        self.assertEqual(expected, actual)

    def test_compute_e_log_q_small(self):
        num_clusters = 10
        num_data_points = 200
        num_dims = 10
        num_grid_points = self.default_grid_size

        self.run_test(num_clusters, num_data_points, num_dims, num_grid_points)

    def test_compute_e_log_q_big(self):
        num_clusters = 200
        num_data_points = 2000
        num_dims = 100
        num_grid_points = self.default_grid_size

        self.run_test(num_clusters, num_data_points, num_dims, num_grid_points)


class TestComputeELogP(unittest.TestCase):
    def __init__(self, method_name: str = ...):
        super().__init__(method_name)

        self.default_grid_size = 100

        self.rng_seed = 242643578967193853558243570818064774262

        self.rng = None

        set_num_threads(10)

    def setUp(self) -> None:
        self.rng = np.random.default_rng(self.rng_seed)

    def run_test(self, num_clusters, num_data_points, num_dims, num_grid_points, log_p_data):
        priors = Priors(num_clusters, num_grid_points, 1.0)
        var_params = VariationalParameters(num_clusters, num_data_points, num_dims, num_grid_points, self.rng)
        expected = compute_e_log_p_old(log_p_data, priors, var_params)
        data_preproc = DataPreprocessor(log_p_data)
        actual = compute_e_log_p(priors, var_params, data_preproc)
        print("Expected = ", expected)
        print("Actual = ", actual)
        np.testing.assert_almost_equal(actual, expected)

    def create_log_p_data(self, depth, num_data_points, num_dims, num_grid_points):
        p_grid = np.full(num_dims, 1.0)
        generator_exp = (
            simulate_binomial_data_point(depth, p_grid, self.rng, num_grid_points) for _ in range(num_data_points)
        )
        log_p_data = np.fromiter(
            generator_exp, dtype=np.dtype((np.float64, (num_dims, num_grid_points))), count=num_data_points
        )
        return log_p_data

    def test_compute_e_log_p_small(self):
        num_clusters = 10
        num_data_points = 10
        num_dims = 10
        num_grid_points = self.default_grid_size
        depth = 100

        log_p_data = self.create_log_p_data(depth, num_data_points, num_dims, num_grid_points)

        self.run_test(num_clusters, num_data_points, num_dims, num_grid_points, log_p_data)

    def test_compute_e_log_p_big(self):
        num_clusters = 100
        num_data_points = 1000
        num_dims = 100
        num_grid_points = self.default_grid_size
        depth = 100

        log_p_data = self.create_log_p_data(depth, num_data_points, num_dims, num_grid_points)

        self.run_test(num_clusters, num_data_points, num_dims, num_grid_points, log_p_data)


if __name__ == "__main__":
    unittest.main()
