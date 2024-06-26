import h5py
import numpy as np
from numba import set_num_threads

from pyclone_vi.data import load_data

import pyclone_vi.inference
import pyclone_vi.post_process


def fit(
    in_file,
    out_file,
    convergence_threshold=1e-6,
    density="binomial",
    annealing_power=1.0,
    max_iters=int(1e4),
    mix_weight_prior=1.0,
    num_annealing_steps=1,
    num_clusters=10,
    num_grid_points=100,
    num_restarts=1,
    num_threads=1,
    precision=200,
    print_freq=100,
    seed=None,
):
    set_num_threads(num_threads)

    rng = instantiate_and_seed_RNG(seed)

    log_p_data, mutations, samples = load_data(
        in_file, density, num_grid_points, precision=precision
    )

    best_elbo = float("-inf")

    result = None

    priors = None

    for i in range(num_restarts):
        print("Performing restart {}".format(i))

        priors = pyclone_vi.inference.get_priors(num_clusters, num_grid_points)

        priors.pi = np.ones(num_clusters) * mix_weight_prior

        var_params = pyclone_vi.inference.get_variational_params(len(priors.pi), log_p_data.shape[0],
                                                                 log_p_data.shape[1], log_p_data.shape[2], rng)

        elbo_trace = pyclone_vi.inference.fit_annealed(
            log_p_data,
            priors,
            var_params,
            annealing_power=annealing_power,
            convergence_threshold=convergence_threshold,
            max_iters=max_iters,
            num_annealing_steps=num_annealing_steps,
            print_freq=print_freq,
        )

        if elbo_trace[-1] > best_elbo:
            best_elbo = elbo_trace[-1]

            result = (elbo_trace, var_params)

        print("Fitting completed")
        print("ELBO: {}".format(elbo_trace[-1]))
        print(
            "Number of clusters used: {}".format(len(set(var_params.z.argmax(axis=1))))
        )
        print()

    elbo_trace, var_params = result

    print("All restarts completed")
    print("Final ELBO: {}".format(elbo_trace[-1]))
    print("Number of clusters used: {}".format(len(set(var_params.z.argmax(axis=1)))))

    _create_fit_results_file(elbo_trace, log_p_data, mutations, out_file, priors, samples, var_params)


def _create_fit_results_file(elbo_trace, log_p_data, mutations, out_file, priors, samples, var_params):
    with h5py.File(out_file, "w") as fh:
        fh.create_dataset(
            "/data/mutations",
            data=np.array(mutations, dtype=h5py.string_dtype(encoding="utf-8")),
        )

        fh.create_dataset(
            "/data/samples",
            data=np.array(samples, dtype=h5py.string_dtype(encoding="utf-8")),
        )

        fh.create_dataset("/data/log_p", data=log_p_data)

        fh.create_dataset("/priors/pi", data=priors.pi)

        fh.create_dataset("/priors/theta", data=priors.theta)

        fh.create_dataset("/var_params/pi", data=var_params.pi)

        fh.create_dataset("/var_params/theta", data=var_params.theta)

        fh.create_dataset("/var_params/z", data=var_params.z)

        fh.create_dataset("/stats/elbo", data=np.array(elbo_trace))


def write_results_file(in_file, out_file, compress=False):
    df = pyclone_vi.post_process.load_results_df(in_file)

    df = pyclone_vi.post_process.fix_cluster_ids(df)

    if compress:
        df.to_csv(
            out_file, compression="gzip", float_format="%.4f", index=False, sep="\t"
        )

    else:
        df.to_csv(out_file, float_format="%.4f", index=False, sep="\t")


def instantiate_and_seed_RNG(seed):
    if seed is not None:
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()
    return rng
