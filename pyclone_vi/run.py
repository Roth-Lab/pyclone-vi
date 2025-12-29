from threadpoolctl import threadpool_limits
import h5py
import numpy as np
from numba import set_num_threads
import click

from pyclone_vi.data import load_data
from pyclone_vi.inference import Priors, fit_pyclone_model, VariationalParameters, DataPreprocessor
from pyclone_vi.post_process import load_results_df, fix_cluster_ids
from pyclone_vi.utils import print_command_header
from pathlib import Path


def fit(
    in_file,
    out_file,
    convergence_threshold=1e-6,
    density="binomial",
    max_iters=int(1e4),
    mix_weight_prior=1.0,
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

    seed_val = print_welcome_message(
        num_restarts,
        density,
        num_threads,
        seed,
        rng,
        num_clusters,
        num_grid_points,
        mix_weight_prior,
    )

    log_p_data, mutations, samples = load_data(in_file, density, num_grid_points, precision=precision)

    priors = Priors(num_clusters, num_grid_points, mix_weight_prior)

    run_var_params_setup_dict = {
        "num_clusters": len(priors.pi),
        "num_data_points": log_p_data.shape[0],
        "num_dims": log_p_data.shape[1],
        "num_grid_points": log_p_data.shape[2],
    }

    data_preproc = DataPreprocessor(log_p_data)

    best_elbo = float("-inf")

    result = None

    click.echo("Running PyClone-VI:\n")

    with threadpool_limits(limits=num_threads, user_api="blas"):
        for i in range(num_restarts):
            click.echo("Performing restart {}".format(i))

            var_params = VariationalParameters(
                run_var_params_setup_dict["num_clusters"],
                run_var_params_setup_dict["num_data_points"],
                run_var_params_setup_dict["num_dims"],
                run_var_params_setup_dict["num_grid_points"],
                rng,
            )

            elbo_trace = fit_pyclone_model(
                priors,
                var_params,
                data_preproc,
                convergence_threshold=convergence_threshold,
                max_iters=max_iters,
                print_freq=print_freq,
            )

            if elbo_trace[-1] > best_elbo:
                best_elbo = elbo_trace[-1]

                result = (elbo_trace, var_params)

            click.echo("Fitting completed")
            click.echo("ELBO: {}".format(elbo_trace[-1]))
            click.echo("Number of clusters used: {}".format(len(set(var_params.z.argmax(axis=1)))))
            click.echo()

    elbo_trace, var_params = result

    click.echo("All restarts completed")
    click.echo("Final ELBO: {}".format(elbo_trace[-1]))
    click.echo("Number of clusters used: {}".format(len(set(var_params.z.argmax(axis=1)))))

    _create_fit_results_file(elbo_trace, log_p_data, mutations, out_file, priors, samples, var_params, seed_val)


def _create_fit_results_file(elbo_trace, log_p_data, mutations, out_file, priors, samples, var_params, seed_val):
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

        fh.attrs["seed"] = str(seed_val)


def write_results_file(in_file, out_file, compress=False):
    print_command_header("Write Results File")

    df = load_results_df(in_file)

    df = fix_cluster_ids(df)

    if compress:
        out_path = Path(out_file)
        if out_path.suffix != ".gz":
            out_file = str(out_path.with_suffix(out_path.suffix + ".gz"))
        df.to_csv(out_file, float_format="%.4f", index=False, sep="\t")
    else:
        df.to_csv(out_file, float_format="%.4f", index=False, sep="\t")

    click.echo("Results table written to:\n{}\n".format(out_file))
    click.echo("#" * 100)


def instantiate_and_seed_RNG(seed):
    if seed is not None:
        rng = np.random.default_rng(seed)
    else:
        rng = np.random.default_rng()
    return rng


def print_welcome_message(
    num_restarts,
    density,
    num_threads,
    seed,
    rng,
    num_clusters,
    num_grid_points,
    mix_weight_prior,
):
    print_command_header("Fit")
    click.echo("Running with the following parameters:\n")
    click.echo("Density: {}".format(density))
    click.echo("Max number of clusters: {}".format(num_clusters))
    click.echo("Number of random restarts: {}".format(num_restarts))
    click.echo("Number of CCF approximation grid points: {}".format(num_grid_points))
    click.echo("Mix weight prior: {}".format(mix_weight_prior))
    click.echo("Number of threads: {}".format(num_threads))
    seed_ret = rng.bit_generator.seed_seq.entropy
    if seed is not None:
        seed_msg = "(user-provided)"
        assert seed_ret == seed
    else:
        seed_msg = "(machine-entropy)"
    click.echo("Random seed: {} {}".format(seed_ret, seed_msg))
    click.echo()
    click.echo("#" * 100)
    click.echo()
    return seed_ret
