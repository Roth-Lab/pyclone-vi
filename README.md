# Overview

PyClone-VI is a fast method for inferring clonal population structure.


## Installation

PyClone-Vi is currently in development so the following proceedure has a few steps.

1. Ensure you have a working `conda` installation.
You can do this by installing [Miniconda](https://conda.io/miniconda.html)

2. Install the required dependencies using conda.
We will create a new `conda` environment with the dependencies.
From inside the checked out PhyClone repository folder run the following.
```
conda create -n pyclone-vi --file requirements.txt --yes
```

3. Activate the `conda` environment.
```
source activate pyclone-vi
```
> Note: You will have to do this whenever you open a new terminal and want to run PyClone-VI.

4. Install PyClone-VI
```
pip install git+ssh://git@github.com/aroth85/pyclone-vi.git
```

5. If everything worked PyClone-VI should be available on the command line.
```
pyclone-vi --help
```

## Usage

## Quick start

1. Activate the `conda` environment
```
source activate pyclone-vi
```

2. Fit the model to the data here we use the TRACERx file provided in the `examples` folder.
This assumes you are running in the base directory of the git repo.
Here we run allowing for up to 40 clusters (clones), using the Beta-Binomial distribution and performing 10 random restarts.
This should take a under five minutes.
```
pyclone-vi fit -i examples/tracerx.tsv -o tracerx.h5 -c 40 -d beta-binomial -r 100
```

3. Next we output the final results from the best random restart.
```
pyclone-vi write-results-file -i tracerx.h5 -o tracerx.tsv
```

### Input format

To run a PyClone-VI analysis you will need to prepare an input file.
The file should be in tab delimited format and have the following columns.
> Note: There is an example file in examples/data/tracerx.tsv

1. mutation_id - Unique identifier for the mutation.
This is free form but should match across all samples.
> Note: PyClone-VI will remove any mutations without entries for all detected samples.
If you have mutations with no data in some samples set their ref/alt counts to 0 for the corresponding sample.

2. sample_id - Unique identifier for the sample.

3. ref_counts - Number of reads matching the reference allele.

4. alt_counts - Number of reads matching the alternate allele.

5. major_cn - Major copy number of segment overlapping mutation.

6. minor_cn - Minor copy number of segment overlapping mutation.

7. normal_cn - Total copy number of segment in healthy tissue.
For autosome this will be two and male sex chromosomes one.

You can include the following optional columns.

1. tumour_content - The tumour content (cellularity) of the sample.
Default value is 1.0 if column is not present.
> Note: In principle this could be different for each mutations/sample.
However it most cases it should be the same for all mutations in a sample.

2. error_rate - Sequencing error rate.
Default value is 0.001 if column is not present.
> Note: Most users will not need to change this value.

### Output format

The results file output by `write-results-file` is in tab delimited format.
There six columns:

1. mutation_id - Mutation identifier as used in the input file.

2. sample_id - Unique identifier for the sample as used in the input file.

3. cluster_id - Most probable cluster or clone the mutation was assigned to.

4. cellular_prevalence - Proportion of malignant cells with the mutation in the sample.
This is also called cancer cell fraction (CCF) in the literature.

5. cellular_prevalence_std - Standard error of the cellular_prevalence estimate.

6. cluster_assignment_prob - Posterior probability the mutation is assigned to the cluster.
This can be used as a confidence score to remove mutations with low probability of belonging to a cluster.


### Running PyClone-VI

PyClone-VI has two sub-commands `fit` and `write-results-file`.
Typical usage is to run `fit` to perform inference and then `write-results-file` to select the best fit and post-process the results.

#### `fit` command

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

#### `write-results-file` command

The `write-results-file` will select the best solution found by the `fit` command and post-process the results.
The output format is tab delimited file which an be imported and manipulated using tools such as R and Python.

The are two mandatory arguments:

* `-i, --in-file` - Path to the output file generated by the `fit` command.

* `-o, --out-file` - Path where the final results will be written in tab delimited format.

There is one optional argument:

* `-c, --compress` - If set the output will be compressed using gzip.
This is useful where a large number mutations are input to reduce the size of the results file.


# License

PyClone-VI is licensed under the GPL v3, see the LICENSE.txt file for details.

# Versions

## 0.1.0

First release.
