[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pyclone-vi/README.html)

PyClone-VI
=========

PyClone-VI is a fast method for inferring clonal population structure.

Paper: [PyClone-VI: scalable inference of clonal population structures using whole genome data](https://doi.org/10.1186/s12859-020-03919-2)

--------

# Overview
1. [Installation](#installation)
2. [Usage: Quickstart Example](#quick-start)
3. [Input Format](#input-format)
4. [Output Format](#output-format)
5. [Running PyClone-VI](#running-pyclone-vi)
   * [Fitting the PyClone-VI Model](#fit-command)
   * [Writing Output](#write-results-file-command)

-------

# Installation

The recommended way to install PyClone-VI is through [conda](https://github.com/conda-forge/miniforge) and the [Bioconda](https://bioconda.github.io/index.html) package channel.

1. Ensure you have a working `conda` installation. You can do this by installing [Miniforge](https://github.com/conda-forge/miniforge#install).
2. Configure the [Bioconda channel](https://bioconda.github.io/#usage):
   ```
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```
3. Once the Bioconda channel has been configured in your conda install:

    * To install into a newly created environment **(Recommended)**:
    
    ```
    conda create --name pyclone-vi pyclone-vi
    ```
> [!NOTE]
> If the above command (step 3) fails due to conda being unable to find the PyClone-VI package, you may need to specify the channel, e.g.:
> ```conda create --name pyclone-vi bioconda::pyclone-vi```
> or
> ```conda create --name pyclone-vi -c bioconda pyclone-vi```

4. Next, test the installation. Activate the `conda` environment.
```
conda activate pyclone-vi
```
> Note: You will have to do this whenever you open a new terminal and want to run PyClone-VI.
5. If everything worked PyClone-VI should be available on the command line.
```
pyclone-vi --help
```

-------

# Quick start

1. Activate the `conda` environment
```
conda activate pyclone-vi
```

2. Fit the model to the data here we use the [TRACERx file](examples/tracerx.tsv) provided in the `examples` folder.
This assumes you are running in the base directory of the git repo.
Here we run allowing for up to 40 clusters (clones), using the Beta-Binomial distribution and performing 10 random restarts.
This should take under five minutes.
```
pyclone-vi fit -i examples/tracerx.tsv -o tracerx.h5 -c 40 -d beta-binomial -r 10
```

3. Next we output the final results from the best random restart.
```
pyclone-vi write-results-file -i tracerx.h5 -o tracerx.tsv
```

-------

# File Formats

## Input format

To run a PyClone-VI analysis you will need to prepare an input file.
The file should be in tab delimited format and have the following columns.
> [!TIP]
> There is an example file in [examples/tracerx.tsv](examples/tracerx.tsv)

1. `mutation_id` - Unique identifier for the mutation.
This is free form but should match across all samples.
> [!WARNING]
> PyClone-VI will remove any mutations without entries for all detected samples.
> If you have mutations with no data in a subset of the samples, the correct procedure is to extract ref and alt counts for these mutations from each affected sample's associated BAM file.
> Please refer to [this thread](https://groups.google.com/g/pyclone-user-group/c/wgXV7tq470Y) for further detail.

2. `sample_id` - Unique identifier for the sample.

3. `ref_counts` - Number of reads matching the reference allele.

4. `alt_counts` - Number of reads matching the alternate allele.

5. `major_cn` - Major copy number of segment overlapping mutation.

6. `minor_cn` - Minor copy number of segment overlapping mutation.

7. `normal_cn` - Total copy number of segment in healthy tissue.
For autosomes this will be two, and for male sex chromosomes one.

You can include the following optional columns.

1. `tumour_content` - The tumour content (cellularity) of the sample.
Default value is 1.0 if column is not present.
> Note: In principle this could be different for each mutation/sample.
However, in most cases it should be the same for all mutations in a sample.

2. `error_rate` - Sequencing error rate.
Default value is 0.001 if column is not present.
> Note: Most users will not need to change this value.

## Output format

The results file output by `write-results-file` is in tab delimited format.
There six columns:

1. `mutation_id` - Mutation identifier as used in the input file.

2. `sample_id` - Unique identifier for the sample as used in the input file.

3. `cluster_id` - Most probable cluster or clone the mutation was assigned to.

4. `cellular_prevalence` - Proportion of malignant cells with the mutation in the sample.
This is also called cancer cell fraction (CCF) in the literature.

5. `cellular_prevalence_std` - Standard error of the cellular_prevalence estimate.

6. `cluster_assignment_prob` - Posterior probability the mutation is assigned to the cluster.
This can be used as a confidence score to remove mutations with low probability of belonging to a cluster.

-------

# Running PyClone-VI

PyClone-VI has two sub-commands `fit` and `write-results-file`.
Typical usage is to run `fit` to perform inference and then `write-results-file` to select the best fit and post-process the results.

## `fit` command

The `fit` command is used for performing the inference step.
It supports performing multiple restarts, the best of which will be selected by the `write-results-file` command.

There are a few mandatory arguments:

* `-i, --in-file` - Path where the input file is located.
This file should be in the format specified above.

* `-o, --out-file` - Path where the output file will be written.
The output file is in HDF5 file format.
Most users will execute the `write-results-file` to extract the final results from this file.

There are several optional arguments:

* `-c, --num-clusters` - The number of clusters to use while fitting.
This should be set to a value larger than the expected number of clusters.
The software will then automatically determine how many to use.
Usually a value of 10-40 will work.
In general this value should increase if as more samples are used.

* `-d, --density` - The probability density used to model the read count data.
Choices are `beta-binomial` and `binomial`.
`binomial` is a common choice for sequencing data.
`beta-binomial` is useful when the data is over-dispersed which has been observed frequently in sequencing data.

* `-g, --num-grid-points` - Number of grid points used for approximating the posterior distribution.
Higher values should be used for deeply sequenced data.
The default value of 100 will likely work for most users.

* `-r, --num-restarts` - Number of random restarts of variational inference.
More restarts will have a higher probability of finding the optimal variational approximation.
This also increases running time.
Usually a value of 10-100 will work.

Additional arguments can be viewed by running `pyclone-vi fit --help`

## `write-results-file` command

The `write-results-file` will select the best solution found by the `fit` command and post-process the results.
The output format is tab delimited file which can be imported and manipulated using tools such as R and Python.

There are two mandatory arguments:

* `-i, --in-file` - Path to the output file generated by the `fit` command.

* `-o, --out-file` - Path where the final results will be written in tab delimited format.

There is one optional argument:

* `-c, --compress` - If set the output will be compressed using gzip.
This is useful where a large number mutations are input to reduce the size of the results file.


# License

PyClone-VI is licensed under the GPL v3, see the [LICENSE.txt](LICENSE.txt) file for details.
