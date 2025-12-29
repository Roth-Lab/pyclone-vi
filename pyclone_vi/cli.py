import click

import pyclone_vi.run
import pathlib


def _validate_out_file(ctx, param, value):
    parent_dir = pathlib.Path(value).parent
    checker = click.Path(exists=True, file_okay=False, dir_okay=True, writable=True)
    checker.convert(parent_dir, param, ctx)
    return value


def _validate_positive_value(ctx, param, value):
    if value < 0:
        raise click.BadParameter("Value must be positive.")
    return value


def _validate_nullable_positive_value(ctx, param, value):
    if value is None:
        return value
    if value < 0:
        raise click.BadParameter("Value must be positive.")
    return value


@click.command(context_settings={"max_content_width": 120}, name="fit")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True, readable=True, file_okay=True, dir_okay=False),
    help="""Path to TSV format file with copy number and allele count information for all samples. """
    """See the examples directory in the GitHub repository for format.""",
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(resolve_path=True, writable=True, file_okay=True, dir_okay=False),
    callback=_validate_out_file,
    help="""Path to where results will be written in HDF5 format.""",
)
@click.option(
    "-c",
    "--num-clusters",
    default=10,
    type=click.IntRange(2, clamp=True),
    show_default=True,
    help="""Number of clusters to use in variational approximation distribution. """
    """Note that not all clusters may not be assigned data points, so the final number of clusters could be lower. """
    """Default is 10.""",
)
@click.option(
    "-d",
    "--density",
    default="binomial",
    type=click.Choice(["beta-binomial", "binomial"]),
    show_default=True,
    help="""Allele count density in the PyClone model. Use beta-binomial for high coverage sequencing. """
    """Default binomial.""",
)
@click.option(
    "-g",
    "--num-grid-points",
    default=100,
    type=click.IntRange(10, clamp=True),
    show_default=True,
    help="""Number of points used to approximate CCF values. Default is 100.""",
)
@click.option(
    "-r",
    "--num-restarts",
    default=1,
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="""Number of random restarts of variational inference. Default is 1.""",
)
@click.option(
    "-t",
    "--num-threads",
    default=1,
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="""Number of threads to use. Default is 1.""",
)
@click.option(
    "--convergence-threshold",
    default=1e-6,
    type=float,
    help="""Maximum relative ELBO difference between iterations to decide on convergence. Default is 10^-6.""",
)
@click.option(
    "--max-iters",
    default=int(1e4),
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="""Maximum number of ELBO optimization iterations. Default is 10,000.""",
)
@click.option(
    "--mix-weight-prior",
    default=1.0,
    type=float,
    show_default=True,
    callback=_validate_positive_value,
    help="""Parameter value of symmetric Dirichlet prior distribution on mixture weights. 
    Higher values will produce more clusters. Default is 1.0 which is the uniform prior.""",
)
@click.option(
    "--precision",
    default=200,
    type=float,
    show_default=True,
    callback=_validate_positive_value,
    help="""Precision for Beta-Binomial density. Has no effect when using Binomial. Default is 200.""",
)
@click.option(
    "--print-freq",
    default=100,
    type=click.IntRange(1, clamp=True),
    show_default=True,
    help="""How often to print information about optimization. Default is every 100 iteration.""",
)
@click.option(
    "--seed",
    default=None,
    type=int,
    callback=_validate_nullable_positive_value,
    show_default=True,
    help="""Set random seed so results can be reproduced. By default, a random seed is chosen.""",
)
def fit(**kwargs):
    """Fit PyClone-VI model to data."""
    pyclone_vi.run.fit(**kwargs)


@click.command(context_settings={"max_content_width": 120}, name="write-results-file")
@click.option(
    "-i",
    "--in-file",
    required=True,
    type=click.Path(exists=True, resolve_path=True, readable=True, file_okay=True, dir_okay=False),
    help="""Path to HDF5 format file produced by the `fit` command.""",
)
@click.option(
    "-o",
    "--out-file",
    required=True,
    type=click.Path(resolve_path=True, writable=True, file_okay=True, dir_okay=False),
    callback=_validate_out_file,
    help="""Path to where results will be written in tsv format.""",
)
@click.option(
    "-c",
    "--compress",
    is_flag=True,
    help="""If the output file should be compressed using gzip.""",
)
def write_results_file(**kwargs):
    """Write the results of a fitted model to file."""
    pyclone_vi.run.write_results_file(**kwargs)


@click.group(name="pyclone-vi")
@click.version_option()
def main():
    pass


main.add_command(fit)
main.add_command(write_results_file)
