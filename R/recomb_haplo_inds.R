#' Construct haplotypes of individuals given their ancestral blocks and the ancestral haplotype variants
#'
#' @param inds A list of individuals in the same format as the output of [recomb_fam()] after being processed with [recomb_map_inds()].
#' More specifically, each individual is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes (always identified by number, but may also be named arbitrarily), each of which is a data.frame/tibble with implicit ranges (`pos` is end coordinates in base pairs; start is the end of the previous block plus one, 1 for the first block) and ancestors `anc` as strings.
#' @param haplo The ancestral haplotypes, which is a list of chromosomes, each of which is a list with two named elements: `X` is a matrix of haplotype markers (loci along rows, ancestral individuals along columns, which must be named as in `anc` strings in `inds` above), and `pos` is a vector of locus positions in base pair coordinates.
#' @param ret_anc If `TRUE`, returns local ancestries (per position) along with haplotypes, otherwise only haplotypes are returned.
#'
#' @return A list of diploid individuals, each of which is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes.
#' If `ret_anc = FALSE` (default), each chromosome is a haplotype (vector of values copied from ancestors in `haplo`);
#' if `ret_anc = TRUE`, each chromosome is a list with named elements `x` for the haplotype vector and `anc` for the vector of ancestor name per position.
#'
#' @examples
#' # Lengthy code creates individuals with recombination data to map
#' # The smallest pedigree, two parents and a child (minimal fam table).
#' library(tibble)
#' fam <- tibble(
#'   id = c('father', 'mother', 'child'),
#'   pat = c(NA, NA, 'father'),
#'   mat = c(NA, NA, 'mother')
#' )
#' # use latest human recombination map, but just first two chrs to keep this example fast
#' map <- recomb_map_hg38[ 1L:2L ]
#' # initialize parents with this other function
#' founders <- recomb_init_founders( c('father', 'mother'), map )
#' # draw recombination breaks for child
#' inds <- recomb_fam( founders, fam )
#' # now add base pair coordinates to recombination breaks
#' inds <- recomb_map_inds( inds, map )
#'
#' # also need ancestral haplotypes
#' # these should be simulated carefully as needed, but for this example we make random data
#' haplo <- vector( 'list', length( map ) )
#' # names of ancestor haplotypes for this scenario
#' # (founders of fam$id but each with "_pat" and "_mat" suffixes)
#' anc_names <- c( 'father_pat', 'father_mat', 'mother_pat', 'mother_mat' )
#' n_ind <- length( anc_names )
#' # number of loci per chr, for toy test
#' m_loci <- 10L
#' for ( chr in 1L : length( map ) ) {
#'     # draw random positions
#'     pos_chr <- sample.int( max( map[[ chr ]]$pos ), m_loci )
#'     # draw haplotypes
#'     X_chr <- matrix(
#'         rbinom( m_loci * n_ind, 1L, 0.5 ),
#'         nrow = m_loci,
#'         ncol = n_ind
#'     )
#'     # required column names!
#'     colnames( X_chr ) <- anc_names
#'     # add to structure, in a list
#'     haplo[[ chr ]] <- list( X = X_chr, pos = pos_chr )
#' }
#'
#' # finally, run desired function!
#' # determine haplotypes of descendants given ancestral haplotypes
#' data <- recomb_haplo_inds( inds, haplo )
#' 
#'
#' @seealso
#' [recomb_fam()] for drawing recombination (ancestor) blocks, defined in terms of genetic distance.
#'
#' [recomb_map_inds()] for transforming genetic to basepair coordinates given a genetic map.
#'
#' [recomb_geno_inds()] for transforming the output of this function from haplotypes (a nested lists structure) to a plain genotype matrix.
#' 
#' @export
recomb_haplo_inds <- function( inds, haplo, ret_anc = FALSE ) {
    # validations
    if ( missing( inds ) )
        stop( '`inds` is required!' )
    if ( missing( haplo ) )
        stop( '`haplo` is required!' )
    if ( !is.list( inds ) )
        stop( '`inds` must be a list!' )
    if ( !is.list( haplo ) )
        stop( '`haplo` must be a list!' )

    # apply ind processing to each ind in inds
    return( lapply( inds, recomb_haplo_ind, haplo, ret_anc = ret_anc ) )
}
