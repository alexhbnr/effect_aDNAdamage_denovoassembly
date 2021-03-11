# Quantifying the effect of ancient DNA damage on the assembled sequences using simulations

This repository contains the code to both conduct the simulation experiments as well as re-generate
the plots of the analysis shown in the Supplementary Information 6 of Wibowo *et al.*:
**Reconstruction of ancient microbial genomes from the human gut (2021)**.

## Simulation experiments

The code to conduct the simulation experiment is in the sub-folder `scripts`. The simulation
pipeline was written in Snakemake v5.3 (Köster *et al.*, 2012) and uses the package manager
[conda](https://docs.conda.io/en/latest/) in order to install all subsequent bioinformatics tools.

### Prepare for running Snakemake

Snakemake is based on **Python3**. Next to the standard modules, the following additional modules
are used and are required to be installed prior to running the workflow.

  - snakemake
  - tqdm
  - numpy
  - pandas
  - scipy
  - pysam
  - pyfastx

These modules can be installed using either `pip`

```
pip install snakemake tqdm numpy pandas scipy pysam pyfastx
```

or `conda`

```
conda install snakemake tqdm numpy pandas scipy pysam pyfastx
```

After having installed these additional modules, the Snakemake workflow can be executed.

### Running the simulations

The simulation experiment consists out of four steps that are all part of the Snakefile:

  1. Download and prepare bacterial reference genomes (**reference_genomes**)
  2. Simulate short-read sequencing data with different underlying read length profiles, amounts of
     ancient DNA damage, and coverage (**simulate**)
  3. Assemble these short-read sequencing data samples using MEGAHIT (**assembly**)
  4. Compare the resulting contig sequences to the reference genomes (**evaluation**)

Each of these steps can be performed independently by specifying the target rule (highlighted in
bold after each step) when executing the workflow.

By default, all these steps are run consecutively. To execute the workflow, run:

```
snakemake -s scripts/simulation.Snakefile -d <output directory> --use-conda -j <number of parallel jobs>
```

For further questions regarding executing `snakemake`, please refer to its 
[documentation](https://snakemake.readthedocs.io/en/stable/index.html).

All results necessary to reproduce the figure can be found in the sub-folder `results`. 

### Summary of contig stats assembled in this study

To obtain **panel a** of the Extended Data Figure 9, we summarised the contigs of high and
medium quality bins assembled in this study regarding their contig length, coverage, amount of
ancient DNA damage on the 5' terminal base, GC content, and the mean read length of the short-read
sequencing data aligning against these contigs.

All parameters but the mean read length were obtained using
[PyDamage](https://github.com/maxibor/pydamage). The mean read length was determined using `bioawk`.
The data is provided in `results/empiricalbins_summarystats.tsv`.

## Generating the Extended Data Figure

The Extended Data Figure was generated using the programming language R using the R package
`ggplot2` within the larger [tidyverse](https://www.tidyverse.org/) framework. The three additional
R packages that were used were `data.table`, `santoku` and `patchwork`, which are all available via
CRAN.

The code is available within a R Markdown document in the sub-folder `analysis`. The markdown report
for viewing the figures directly on GitHub can be found [here](analysis/analysis.md)
