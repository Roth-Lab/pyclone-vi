# Overview

PyClone-VI is fast method for inferring clonal population structure.


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
> Note: You will have to do this whenver you open a new terminal and want to run PyClone-VI. 

4. Install PyClone-VI

```
pip install git+ssh://git@github.com/aroth85/pyclone-vi.git
```

5. If everything worked PyClone-VI should be available on the command line.
```
pyclone-vi --help
```

## Usage

### Input format

To run a PyClone-VI analysis you will need to prepare an input file.
The file should be in tab delimited format and have the following columns.
> Note: There is an example file in examples/data/mixing.tsv

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

### Running PyClone-VI

PyClone-VI analyses are broken into two parts.
First the model is fit using the `fit` command.
Second the fitted model is post-processed and a results file produced using the `write-results-file` command.

Sampling can be run as follows
```
phyclone run -i INPUT.tsv -t TRACE.pkl 
``` 
which will take the file `INPUT.tsv` as described above and write the trace file `TRACE.pkl` in a Python pickle format.

The `-n` command can be used to control the number of iterations of sampling to perform.

The `-b` command can be used to control the number of burnin iterations to perform.
> Note: burnin is done using a heuristic strategy of unconditional SMC.
All samples from the burnin are discarded as they will not target the posterior.

The `-d` command can be used to select the emission density.
As in PyClone the `binomial` and `beta-binomial` densities are available.
> Note: Unlike PyClone, PhyClone does not estimate the precision parameter of the Beta-Binomial.
This parameter can be set with the --precision flag.

To build the final results of PhyClone you can run the `consensus` command as follows.
```
phyclone consensus -d INPUT.tsv -t TRACE.pkl -n TREE.nwk -o TABLE.tsv
``` 
where `INPUT.tsv` and `TRACE.pkl` are the same as the previous step and `TREE.nwk` is the clone tree in newick format and `TABLE.tsv` the assignment of mutations to clones.

# License

PhyClone is licensed under the GPL v3, see the LICENSE.txt file for details.

# Versions

## 0.1.0

