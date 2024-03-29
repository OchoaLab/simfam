% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recomb_geno_inds.R
\name{recomb_geno_inds}
\alias{recomb_geno_inds}
\title{Reduce haplotype data to genotype matrix}
\usage{
recomb_geno_inds(haplos)
}
\arguments{
\item{haplos}{A list of diploid individuals, each of which is a list with two haploid individuals named "pat" and "mat", each of which is a list of chromosomes.
Each chromosome can be a list, in which case the named element "x" must give the haplotype vector (ideally with values in zero and one counting reference alleles, including NA), otherwise the chromosome must be this vector (accommodating both output formats from \code{\link[=recomb_haplo_inds]{recomb_haplo_inds()}} automatically).}
}
\value{
The genotype matrix, which is the sum of the haplotype values (with values in 0, 1, 2, and NA, counting reference alleles), with individuals along columns in same order as \code{haplos} list, and loci along rows in order of appearance concatenating chromosomes in numerical order.
}
\description{
This function accepts haplotype data, such as the output from \code{\link[=recomb_haplo_inds]{recomb_haplo_inds()}}, and reduces it to a genotype matrix.
The haplotype data is more detailed because it is phased, while phase is lost in the genotype representation.
Moreover, the haplotype data separates individuals and chromosomes into lists (the way it is simulated), but the output genotype matrix concatenates data from all chromosomes into a single matrix, as it appears in simpler simulations and real data.
}
\examples{
# Lengthy code creates individuals with recombination data to map
# The smallest pedigree, two parents and a child (minimal fam table).
library(tibble)
fam <- tibble(
  id = c('father', 'mother', 'child'),
  pat = c(NA, NA, 'father'),
  mat = c(NA, NA, 'mother')
)
# use latest human recombination map, but just first two chrs to keep this example fast
map <- recomb_map_hg38[ 1L:2L ]
# initialize parents with this other function
founders <- recomb_init_founders( c('father', 'mother'), map )
# draw recombination breaks for child
inds <- recomb_fam( founders, fam )
# now add base pair coordinates to recombination breaks
inds <- recomb_map_inds( inds, map )

# also need ancestral haplotypes
# these should be simulated carefully as needed, but for this example we make random data
haplo <- vector( 'list', length( map ) )
# names of ancestor haplotypes for this scenario
# (founders of fam$id but each with "_pat" and "_mat" suffixes)
anc_names <- c( 'father_pat', 'father_mat', 'mother_pat', 'mother_mat' )
n_ind <- length( anc_names )
# number of loci per chr, for toy test
m_loci <- 10L
for ( chr in 1L : length( map ) ) {
    # draw random positions
    pos_chr <- sample.int( max( map[[ chr ]]$pos ), m_loci )
    # draw haplotypes
    X_chr <- matrix(
        rbinom( m_loci * n_ind, 1L, 0.5 ),
        nrow = m_loci,
        ncol = n_ind
    )
    # required column names!
    colnames( X_chr ) <- anc_names
    # add to structure, in a list
    haplo[[ chr ]] <- list( X = X_chr, pos = pos_chr )
}
# determine haplotypes of descendants given ancestral haplotypes
haplos <- recomb_haplo_inds( inds, haplo )

# finally, run desired function!
# convert haplotypes structure to a plain genotype matrix
X <- recomb_geno_inds( haplos )

}
\seealso{
\code{\link[=recomb_fam]{recomb_fam()}} for drawing recombination (ancestor) blocks, defined in terms of genetic distance.

\code{\link[=recomb_map_inds]{recomb_map_inds()}} for transforming genetic to basepair coordinates given a genetic map.

\code{\link[=recomb_haplo_inds]{recomb_haplo_inds()}} for determining haplotypes of descendants given ancestral haplotypes (creates input to this function).
}
