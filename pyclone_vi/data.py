from collections import OrderedDict
from numba import njit, prange, int64, float64

from numba.experimental import jitclass
import numpy as np
import pandas as pd

from pyclone_vi.math_utils import (
    log_beta_binomial_pdf,
    log_binomial_pdf,
    log_sum_exp,
)

from functools import lru_cache


def load_data(file_name, density="binomial", num_grid_points=100, precision=200):

    print("Parsing Input Data...\n")

    data, samples = load_pyclone_data(file_name)

    mutations = data.index.to_list()

    log_p_data = data.apply(dp_likelihood_grid, args=(density, num_grid_points, precision))

    log_p_data = log_p_data.to_numpy(dtype=np.dtype((np.float64, (len(samples), num_grid_points))))

    print("#" * 100)
    print()

    return log_p_data, mutations, samples


def dp_likelihood_grid(dp, density, num_grid_points, precision):
    return dp.to_likelihood_grid(density, num_grid_points, precision=precision)


def load_pyclone_data(file_name):
    df = pd.read_csv(file_name, sep="\t")
    df = df.drop_duplicates()

    df = _remove_cn_zero_mutations(df)

    _process_required_columns(df)

    samples = sorted(df["sample_id"].unique())

    df = _remove_duplicated_and_partially_absent_mutations(df, samples)

    data = _create_loaded_pyclone_data_dict(df, samples)

    get_major_cn_prior.cache_clear()

    print("Num Samples: {}".format(len(samples)))
    if len(samples) > 10:
        print("Samples: {}...".format(" ".join(samples[:5])))
    else:
        print("Samples: {}".format(" ".join(samples)))

    print("Num Mutations: {}".format(len(data)))
    print()

    return data, samples


def _create_loaded_pyclone_data_dict(df, samples):
    df.set_index("sample_id", inplace=True)
    df.sort_index(inplace=True)
    grouped = df.groupby("mutation_id", sort=False)
    samples = pd.Index(samples, name="sample_id")

    data = grouped.apply(make_datapoint_from_group, samples=samples, include_groups=False)

    return data


def make_datapoint_from_group(group, samples):
    sample_dp_df = group.agg(create_sample_data_point, axis=1)
    return DataPoint(samples, sample_dp_df)


def create_sample_data_point(row_series):
    major_cn = int(row_series["major_cn"])
    minor_cn = row_series["minor_cn"]
    normal_cn = row_series["normal_cn"]
    error_rate = row_series["error_rate"]
    ref_count = row_series["ref_counts"]
    alt_count = row_series["alt_counts"]
    tumour_content = row_series["tumour_content"]

    cn, mu, log_pi = get_major_cn_prior(
            major_cn,
            minor_cn,
            normal_cn,
            error_rate,
        )

    sample_dp = SampleDataPoint(ref_count, alt_count, cn, mu, log_pi, tumour_content)

    return sample_dp


def _process_required_columns(df):
    df["sample_id"] = df["sample_id"].astype(str)
    if "error_rate" not in df.columns:
        print("Error rate column not found, setting values to {}.\n".format(1e-3))
        df["error_rate"] = 1e-3

    if "tumour_content" not in df.columns:
        print("Tumour content column not found, setting values to 1.0.\n")
        df["tumour_content"] = 1.0


def _remove_cn_zero_mutations(df):
    num_dels = len(df.loc[df["major_cn"] == 0])
    if num_dels > 0:
        print("Removing {} mutations with major copy number zero".format(num_dels))
    df = df.loc[df["major_cn"] > 0]
    return df


def _remove_duplicated_and_partially_absent_mutations(df, samples):
    samples_len = len(samples)
    group_transform = df.groupby("mutation_id")["sample_id"].transform("size")
    num_not_present_in_all = len(df.loc[group_transform < samples_len]["mutation_id"].unique())
    num_duplicates = len(df.loc[group_transform > samples_len]["mutation_id"].unique())
    if num_duplicates > 0:
        if num_duplicates == 1:
            pl = ""
        else:
            pl = "s"
        print("Removing {} duplicate mutation ID{}".format(num_duplicates, pl))
    if num_not_present_in_all > 0:
        if num_not_present_in_all == 1:
            pl = ("", "is")
        else:
            pl = ("s", "are")
        print(
            "Removing {} mutation{} that {} not present in all samples".format(
                num_not_present_in_all,
                pl[0],
                pl[1],
            )
        )
    df = df.loc[group_transform == samples_len]

    if (num_duplicates > 0) or (num_not_present_in_all > 0):
        print()

    return df


@lru_cache(maxsize=4096)
def get_major_cn_prior(major_cn, minor_cn, normal_cn, error_rate=1e-3):
    total_cn = major_cn + minor_cn

    # Consider all possible mutational genotypes consistent with mutation before CN change
    cn = [(normal_cn, normal_cn, total_cn) for _ in range(1, major_cn + 1)]
    mu = [(error_rate, error_rate, min(1 - error_rate, x / total_cn)) for x in range(1, major_cn + 1)]

    # Consider mutational genotype of mutation before CN change if not already added
    if total_cn != normal_cn:
        mutation_after_cn = (normal_cn, total_cn, total_cn)
        cn.append(mutation_after_cn)
        mu.append((error_rate, error_rate, min(1 - error_rate, 1 / total_cn)))
        assert len(set(cn)) == 2

    cn = np.array(cn, dtype=np.int64)
    mu = np.array(mu, dtype=np.float64)

    log_pi_val = -np.log(len(cn))
    log_pi = np.full(len(cn), log_pi_val)

    return cn, mu, log_pi


class DataPoint(object):
    __slots__ = "samples", "sample_data_points"

    def __init__(self, samples, sample_data_points):
        self.samples = samples
        self.sample_data_points = sample_data_points

        if not samples.equals(sample_data_points.index):
            self.sample_data_points.sort_index(inplace=True)
            assert samples.equals(sample_data_points.index)

    @staticmethod
    def get_ccf_grid(grid_size, eps=1e-6):
        return np.linspace(eps, 1 - eps, grid_size)

    def to_dict(self):
        return self.sample_data_points.to_dict(into=OrderedDict)

    def to_likelihood_grid(self, density, num_grid_points, precision=200):
        grid = self.get_ccf_grid(num_grid_points)

        if density == "beta-binomial":
            grid_res = self.sample_data_points.apply(log_pyclone_beta_binomial_pdf_grid_helper, args=(grid, precision, num_grid_points),)
        elif density == "binomial":
            grid_res = self.sample_data_points.apply(log_pyclone_binomial_pdf_grid_helper, args=(grid, num_grid_points),)
        else:
            raise NotImplemented("Unknown density: {}".format(density))

        log_ll = grid_res.to_numpy(dtype=np.dtype((np.float64, num_grid_points)))
        return log_ll


@jitclass(
    [
        ("a", int64),
        ("b", int64),
        ("cn", int64[:, :]),
        ("mu", float64[:, :]),
        ("log_pi", float64[:]),
        ("t", float64)
    ]
)
class SampleDataPoint(object):
    def __init__(self, a, b, cn, mu, log_pi, t):
        self.a = a
        self.b = b
        self.cn = cn
        self.mu = mu
        self.log_pi = log_pi
        self.t = t


def log_pyclone_beta_binomial_pdf_grid_helper(data_point, grid, precision, num_grid_points):
    log_ll = np.empty(num_grid_points, dtype=np.float64, order="C")
    log_pyclone_beta_binomial_pdf_grid(data_point, grid, precision, log_ll)
    return log_ll


def log_pyclone_binomial_pdf_grid_helper(data_point, grid, num_grid_points):
    log_ll = np.empty(num_grid_points, dtype=np.float64, order="C")
    log_pyclone_binomial_pdf_grid(data_point, grid, log_ll)
    return log_ll


@njit(parallel=True)
def log_pyclone_beta_binomial_pdf_grid(data_point, grid, precision, log_ll):
    for i in prange(len(grid)):
        log_ll[i] = log_pyclone_beta_binomial_pdf(data_point, grid[i], precision)


@njit(parallel=True)
def log_pyclone_binomial_pdf_grid(data_point, grid, log_ll):
    for i in prange(len(grid)):
        log_ll[i] = log_pyclone_binomial_pdf(data_point, grid[i])


@njit
def log_pyclone_beta_binomial_pdf(data, f, s):
    t = data.t

    C = len(data.cn)

    population_prior = np.zeros(3)
    population_prior[0] = 1 - t
    population_prior[1] = t * (1 - f)
    population_prior[2] = t * f

    ll = np.full(C, -np.inf, dtype=np.float64)

    for c in range(C):
        e_vaf = 0

        norm_const = 0

        for i in range(3):
            e_cn = population_prior[i] * data.cn[c, i]

            e_vaf += e_cn * data.mu[c, i]

            norm_const += e_cn

        e_vaf /= norm_const

        a = e_vaf * s

        b = s - a

        ll[c] = data.log_pi[c] + log_beta_binomial_pdf(data.a + data.b, data.b, a, b)

    return log_sum_exp(ll)


@njit
def log_pyclone_binomial_pdf(data, f):
    t = data.t

    C = len(data.cn)

    population_prior = np.zeros(3)
    population_prior[0] = 1 - t
    population_prior[1] = t * (1 - f)
    population_prior[2] = t * f

    ll = np.full(C, -np.inf, dtype=np.float64)

    for c in range(C):
        e_vaf = 0

        norm_const = 0

        for i in range(3):
            e_cn = population_prior[i] * data.cn[c, i]

            e_vaf += e_cn * data.mu[c, i]

            norm_const += e_cn

        e_vaf /= norm_const

        ll[c] = data.log_pi[c] + log_binomial_pdf(data.a + data.b, data.b, e_vaf)

    return log_sum_exp(ll)
