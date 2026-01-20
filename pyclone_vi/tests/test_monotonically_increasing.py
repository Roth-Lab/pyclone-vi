import unittest
from threadpoolctl import threadpool_limits
import numpy as np
from pyclone_vi.inference import VariationalParameters, Priors, DataPreprocessor, fit_pyclone_model
from pyclone_vi.tests.simulate import simulate_binomial_data_point
from numba import set_num_threads
import pandas as pd


class TestMonotonicallyIncreasingELBO(unittest.TestCase):

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
        num_grid_points = 200
        depth = 100

        log_p_data, priors, var_params, data_preproc = self.create_test_structures(
            depth,
            num_clusters,
            num_data_points,
            num_dims,
            num_grid_points,
        )

        return log_p_data, var_params, priors, data_preproc

    def create_test_structures_small(self):
        num_clusters = 10
        num_data_points = 10
        num_dims = 10
        num_grid_points = self.default_grid_size
        depth = 100

        log_p_data, priors, var_params, data_preproc = self.create_test_structures(
            depth,
            num_clusters,
            num_data_points,
            num_dims,
            num_grid_points,
        )

        return log_p_data, var_params, priors, data_preproc

    def create_test_structures(self, depth, num_clusters, num_data_points, num_dims, num_grid_points):
        log_p_data = self.create_log_p_data(depth, num_data_points, num_dims, num_grid_points)
        priors = Priors(num_clusters, num_grid_points, 1.0)
        var_params = VariationalParameters(num_clusters, num_data_points, num_dims, num_grid_points, self.rng)
        data_preproc = DataPreprocessor(log_p_data)
        return log_p_data, priors, var_params, data_preproc

    def test_monotonically_increasing_ELBO_small(self):
        log_p_data, var_params, priors, data_preproc = self.create_test_structures_small()
        with threadpool_limits(limits=self.num_threads, user_api="blas"):
            elbo_trace = fit_pyclone_model(priors,
                                           var_params,
                                           data_preproc,
                                           print_freq=1)

        elbo_series = pd.Series(elbo_trace)

        self.assertTrue(elbo_series.is_monotonic_increasing)

    def test_monotonically_increasing_ELBO_big(self):
        log_p_data, var_params, priors, data_preproc = self.create_test_structures_big()
        with threadpool_limits(limits=self.num_threads, user_api="blas"):
            elbo_trace = fit_pyclone_model(priors,
                                           var_params,
                                           data_preproc,
                                           print_freq=1)

        elbo_series = pd.Series(elbo_trace)

        self.assertTrue(elbo_series.is_monotonic_increasing)


if __name__ == "__main__":
    unittest.main()
