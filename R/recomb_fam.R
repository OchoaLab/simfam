#' Draw recombination breaks for autosomes from a pedigree
#'
#' Create random recombination breaks for all autosomes of all individuals in the provided pedigree FAM table.
#' Recombination lengths follow an exponential distribution with mean of 100 centiMorgans (cM).
#' The output specifies identical-by-descent (IBD) blocks as ranges per chromosome (per individual) and the founder chromosome they arose from (are IBD with).
#' All calculations are in terms of genetic distance (not base pairs), and no genotypes are constructed/drawn in this step.
#'
#' @inheritParams geno_fam
#' @param founders The named list of founders with their chromosomes.
#' For unstructured founders, initialize with [recomb_init_founders()].
#' Each element of this list is a diploid individual, which is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes (always identified by number, but may also be named arbitrarily), each of which is a data.frame/tibble with implicit ranges (`posg` is end coordinates in cM; start is the end of the previous block, zero for the first block) and ancestors `anc` as strings.
#' For true founders each chromosome may be trivial (each chromosome is a single block with ID equal to itself but distinguishing maternal from paternal copy), but input itself can be recombined (for iterating).
#' This list must have names that identify each founder (matching codes in `fam$id`).
#' Individuals may be in a different order than `fam$id`.
#' Extra individuals in `founders` but absent in `fam$id` will be silently ignored.
#'
#' @return The list of individuals with recombined chromosomes of the entire `fam` table, in the same format as `founders` above.
#' The names of this list correspond to `fam$id` in that order.
#'
#' @examples
#' # The smallest pedigree, two parents and a child.
#' # A minimal fam table with the three required columns.
#' # Note "mother" and "father" have missing parent IDs, while "child" does not
#' library(tibble)
#' fam <- tibble(
#'   id = c('father', 'mother', 'child'),
#'   pat = c(NA, NA, 'father'),
#'   mat = c(NA, NA, 'mother')
#' )
#' 
#' # initialize parents with this other function
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' ids <- c('mother', 'father')
#' # simulate three chromosomes with these lengths in cM
#' lengs <- c(50, 100, 150)
#' founders <- recomb_init_founders( ids, lengs )
#'
#' # draw recombination breaks for the whole fam table now:
#' inds <- recomb_fam( founders, fam )
#' 
#' # This is a length-3 list with names matching fam$id.
#' # The parent data equals the input (reordered),
#' # but now there's data to the child too
#' inds
#' 
#' @seealso
#' [recomb_init_founders()] to initialize `founders` for this function.
#' 
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
recomb_fam <- function( founders, fam, missing_vals = c('', 0) ) {
    if ( missing( founders ) )
        stop( '`founders` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.list( founders ) )
        stop( '`founders` must be a list!' )
    
    # ensures that `founders` and `fam` agree, maps parent indexes
    # also loads of `fam` validations
    data <- match_fam_founders( fam, names( founders ), 'founders', 'list', missing_vals )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    n_ind_fam <- nrow( fam )
    
    # initialize output list
    inds <- vector( 'list', n_ind_fam )
    # copy names of individuals and ancestries
    names( inds ) <- fam$id
    # copy founders
    # reorder/subset founders if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    inds[ fam$founder ] <- founders[ indexes ]
    
    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, using precalculated indexes of the current inds
        mat <- inds[[ fam$mati[ i ] ]]
        pat <- inds[[ fam$pati[ i ] ]]
        
        # each parent in turn has two haploid copies of each chromosome, recombine each those!
        child_mat <- recomb_hap( mat$mat, mat$pat )
        child_pat <- recomb_hap( pat$mat, pat$pat )
        
        # store child's two sets of chromosomes to output, as named list
        inds[[ i ]] <- list( pat = child_pat, mat = child_mat )
    }
    
    # return entire matrix!
    return( inds )
}
