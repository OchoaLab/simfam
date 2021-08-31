# simfam

The goal of simfam is to simulate and model families with founders drawn from a structured population.
The main function simulates a random pedigree for many generations with realistic features.
Additional functions calculate kinship matrices, admixture matrices, and draw random genotypes across arbitrary pedigree structures starting from the corresponding founder values.

## Installation

You can install the released version of simfam from [CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("simfam")
```

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('OchoaLab/simfam', build_vignettes = TRUE)
```

You can see the package vignette, which has more detailed documentation and examples, by typing this into your R session:
```R
vignette('simfam')
```

## Examples

These are some basic ways of calling the main functions.

``` r
# load package!
library(simfam)
```

Simulate a random pedigree with a desired number of individuals per generation `n` and a number of generations `G`:
```r
data <- sim_pedigree( n, G )
# creates a plink-formatted FAM table
# (describes pedigree, most important!)
fam <- data$fam
# lists of IDs split by generation
ids <- data$ids
# and local kinship of last generation
kinship_local_G <- data$kinship_local
```

The basics of encoding a pedigree in a `fam` table (a data.frame) is that every individual in the pedigree is a row, column `id` identifies the individual with a unique number or string, columns `pat` and `mat` identify the parents of the individual (who are themselves earlier rows), and `sex` encodes the sex of the individual numerically (1=male, 2=female).
The following functions work with arbitrary pedigrees/`fam` data.frames:

Prune a given `fam`, to speed up simulations/etc, by removing individuals without descendants among set of individuals `ids` (in this example, the last generation from the output of `sim_pedigree`):
```r
fam <- prune_fam( fam, ids[[G]] )
```

Draw genotypes `X` through pedigree, starting from genotypes of founders (`X_1`):
```r
X <- geno_fam( X_1, fam )
# Version for last generation only, which uses less memory.
# (`ids` must be as from `sim_pedigree`,
# a list partitioning non-overlapping generations)
X_G <- geno_last_gen( X_1, fam, ids )
```

Calculate kinship through pedigree, starting from kinship of founders (`kinship_1`):
```r
kinship <- kinship_fam( kinship_1, fam )
# Version for last generation only, which uses less memory.
kinship_G <- kinship_last_gen( kinship_1, fam, ids )
```

Calculate expected admixture proportions through pedigree, starting from admixture of founders (`admix_proportions_1`):
```r
admix_proportions <- admix_fam( admix_proportions_1, fam )
# Version for last generation only, which uses less memory.
admix_proportions_G <- admix_last_gen( admix_proportions_1, fam, ids )
```

