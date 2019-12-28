import h5py
import numpy as np
import pandas as pd

from pyclone_vi.data import load_data

import pyclone_vi.inference


def fit(
        in_file,
        out_file,
        convergence_threshold=1e-6,
        density='binomial',
        annealing_power=1.0,
        max_iters=int(1e4),
        mix_weight_prior=1.0,
        num_annealing_steps=None,
        num_clusters=10,
        num_grid_points=100,
        precision=200,
        print_freq=100,
        seed=None):

    if seed is not None:
        np.random.seed(seed)

    log_p_data, mutation_ids = load_data(in_file, density, num_grid_points, precision=precision)

    priors = pyclone_vi.inference.get_priors(num_clusters, num_grid_points)

    priors.pi = np.ones(num_clusters) * mix_weight_prior

    var_params = pyclone_vi.inference.get_variational_params(
        len(priors.pi),
        log_p_data.shape[0],
        log_p_data.shape[1],
        log_p_data.shape[2]
    )

    if num_annealing_steps is None:
        elbo_trace = pyclone_vi.inference.fit(
            log_p_data,
            priors,
            var_params,
            convergence_threshold=convergence_threshold,
            max_iters=max_iters,
            print_freq=print_freq
        )

    else:
        elbo_trace = pyclone_vi.inference.fit_annealed(
            log_p_data,
            priors,
            var_params,
            annealing_power=annealing_power,
            convergence_threshold=convergence_threshold,
            max_iters=max_iters,
            num_annealing_steps=num_annealing_steps,
            print_freq=print_freq
        )

    print('Fitting completed')
    print('Final ELBO: {}'.format(elbo_trace[-1]))
    print('Number of clusters used: {}'.format(len(set(var_params.z.argmax(axis=1)))))

    with h5py.File(out_file, 'w') as fh:
        fh.create_dataset(
            '/data/mutation_ids',
            data=np.array(mutation_ids, dtype=h5py.string_dtype(encoding='utf-8'))
        )

        fh.create_dataset('/priors/pi', data=priors.pi)

        fh.create_dataset('/priors/theta', data=priors.theta)

        fh.create_dataset('/var_params/pi', data=var_params.pi)

        fh.create_dataset('/var_params/theta', data=var_params.theta)

        fh.create_dataset('/var_params/z', data=var_params.z)

        fh.create_dataset('/stats/elbo', data=np.array(elbo_trace))


def write_cluster_file(in_file, out_file):
    with h5py.File(in_file, 'r') as fh:
        m = fh['/data/mutation_ids'][()]

        z = fh['/var_params/z'][()]

    clusters = np.argmax(z, axis=1)

    df = pd.DataFrame({
        'mutation_id': m,
        'cluster': clusters,
        'cluster_prob': np.max(z, axis=1)
    })

    cluster_map = dict(zip(
        df['cluster'].unique(), np.arange(df['cluster'].nunique())
    ))

    df['cluster'] = df['cluster'].map(cluster_map)

    df.to_csv(out_file, index=False, sep='\t')
