# binXdrives

This tool models the impact of a Bipartite Expression Drive (BED) on a population. The BED may be intended to spread 1) toxic proteins that sterilize BED-carrying females or 2) lethal seminal proteins that BED-carrying males transfer to females, killing them right after mating.

The simulations produced for the publication have been packed into gzip files acording to the figure of the paper they correspond to and can be found [here](simulations/paper/). The commands to run each of them are detailed [here](simulations.md).


## Dependencies

This pipeline makes use of some third-party software, packages, and databases. The versions provided below are the ones that were used to run the analysis for the publication.

### Software

- r-base v4.0.0

### R packages

- bettermc v1.1.1
- dplyr v1.0.7
- ggplot2 v3.3.5
- logger v0.2.0
- optparse v1.6.6
- rlist v0.4.6.1 
- stringr v1.4.0
- zeallot v0.1.0


## Prerequisites

- If you DON'T have all the above mentioned software and packages available in your environment, you could run the tool using CONDA. A minimal installer can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html).


## Installation

1. Clone the repo:

```bash
git clone https://github.com/santiagorevale/binXdrives.git
```

2. Manually install R and its dependencies or give it a go at steps 3 and 4.

3. Create a CONDA environment based on the provided `environment.yml` file:

```bash
cd binXdrives
conda env create -f environment.yml
```

4. The R package `bettermc` can't be install through CONDA, so you will have to install it from R. Use the following lines to install the package. Feel free to choose the repos that are more convenient for you, or just open `R` and plainly run `install.packages(c("bettermc"))`.

```bash
R --no-save << EOT
repos = c(
  CRAN="http://cran.ma.imperial.ac.uk/",
  CRAN_RSTUDIO="https://cran.rstudio.com/"
);
install.packages(c("bettermc"), repos=repos)
EOT
```

## Usage

Here's an example command on how to run the tool. This is just an illustrative example, because many of the below specified parameters are using the default values.

```bash
Rscript --vanilla binXdrives.R \
    --output fertility_BED_reference \
    --threads 8 \
    --simulations 6 \
    --generations 24 \
    --bed_design yes \
    --carrying_capacity 1800 \
    --fecundity 48 \
    --release_A 60 \
    --release_B 60 \
    --release_critical_difference 100 \
    --polyandry 3 \
    --unintended_reproductive_cost_A 0.025 \
    --unintended_reproductive_cost_B 0.025 \
    --unintended_viability_cost_A 0.025 \
    --unintended_viability_cost_B 0.025 \
    --terminator_efficiency 0 \
    --dominance 0.6 \
    --conversion_efficiency 0.9 \
    --intended_fecundity_cost 1 \
    --resistance_formation 0.05 \
    --resistance_functionality 0 \
    --larval_survival 0.5
```


## Output

By default, results will be stored in the folder `results` in the same path where the pipeline was executed.

The final output will be comprised by the following files:

- `all_simulations.rds`: R object containing all the simulations.
- `simulation_n.txt`: one file per simulation containing its results in a text format.
- `plot.pdf`: population size and allele frequencies simulations draft plot.
- `summary.txt`: it compiles an average value for each variable per generation.
- `suppression.txt`: it compiles for each simulation the time until population elimination.


## License

[MIT](LICENSE)


## Authors

- [@santiagorevale](https://github.com/santiagorevale)
- [@hurtadojuan](https://github.com/hurtadojuan)


## Citation

If you use **binXdrives** for your analysis, you can cite it as follows:

> **Propagation of seminal toxins through binary expression gene drives could suppress populations**
>
> Hurtado, Juan; Revale, Santiago; Matzkin, Luciano M.
>
> _Unpublished work_ DATE_TBD. doi: [TBD]().


## References

### Software tools

* [R Project](https://www.r-project.org/)
    > R Core Team (2019) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

### R packages

* [bettermc](https://cran.r-project.org/package=bettermc)
    > Kersting A, et al. (2021). bettermc: Enhanced Fork-Based Parallelization. R package version 1.1.1.

* [dplyr](https://cran.r-project.org/package=dplyr)
    > Wickham H, et al. (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7.

* [ggplot2](https://cran.r-project.org/package=ggplot2)
    > Wickham H, et al. (2021). ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics. R package version 3.3.5.

* [logger](https://cran.r-project.org/package=logger)
    > Daróczi G. (2021). logger: A Lightweight, Modern and Flexible Logging Utility. R package version 0.2.0.

* [optparse](https://CRAN.R-project.org/package=optparse)
    > Davis TL. (2019). optparse: Command Line Option Parser. R package version 1.6.6.

* [stringr](https://CRAN.R-project.org/package=stringr)
    > Wickham H. (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0.

* [rlist](https://cran.r-project.org/package=rlist)
    > Ren K. (2021). rlist: A Toolbox for Non-Tabular Data Manipulation. R package version 0.4.6.1.

* [zeallot](https://cran.r-project.org/package=zeallot)
    > Teetor N. (2018). zeallot: Multiple, Unpacking, and Destructuring Assignment. R package version 0.1.0.

### Software packaging/containerisation tools

* [Anaconda](https://anaconda.com)
    > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

* [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)
    > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.
