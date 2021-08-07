#' Calculate admixture matrix of a pedigree with known admixture of founders
#'
#' Calculates a full admixture proportions matrix (for all individuals in the provided pedigree FAM table) starting from the admixture proportions of the founders as provided.
#'
#' @inheritParams geno_fam
#' @param admix The admixture proportions matrix of the founders (individuals along rows and ancestries along columns).
#' This matrix must have row names that identify each founder (matching codes in `fam$id`).
#' Individuals may be in a different order than `fam$id`.
#' Extra individuals in `admix` but absent in `fam$id` will be silently ignored.
#' All values should be non-negative and each row of `admix` should sum to one; for speed, this code does not check that `admix` is valid, just averages data as-is.
#'
#' @return The admixture proportions matrix of the entire `fam` table, based on the admixture of the founders.
#' These are expectations, calculated for each individual as the average ancestry proportion of the parents.
#' The rows of this admixture matrix correspond to `fam$id` in that order.
#' The columns (ancestries) are the same as in the input `admix`.
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
#' # admixture proportions of the parents
#' admix <- rbind( c(0.3, 0.3, 0.4), c(0.5, 0.25, 0.25) )
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' rownames( admix ) <- c('mother', 'father')
#' # name ancestries too
#' colnames( admix ) <- c('African', 'European', 'Asian')
#'
#' # Calculate the full admixture proportions matrix
#' admix_all <- admix_fam( admix, fam )
#' 
#' # This is a 3x3 matrix with row names matching fam$id.
#' # The parent submatrix equals the input (reordered),
#' # but now there's admixture to the child too (averages of parents)
#' admix_all
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
admix_fam <- function( admix, fam, missing_vals = c('', 0) ) {
    if ( missing( admix ) )
        stop( '`admix` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( admix ) )
        stop( '`admix` must be a matrix!' )
    
    # ensures that `admix` and `fam` agree, maps parent indexes
    # also loads of `fam` validations
    data <- match_fam_founders( fam, rownames( admix ), 'admix', 'row', missing_vals )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    n_ind_fam <- nrow( fam )
    k_subpops <- ncol( admix )
    
    # initialize output matrix
    admix_fam <- matrix(
        NA,
        nrow = n_ind_fam,
        ncol = k_subpops
    )
    # copy names of individuals and ancestries
    rownames( admix_fam ) <- fam$id
    colnames( admix_fam ) <- colnames( admix )
    # copy founders
    # reorder/subset columns of X if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    admix_fam[ fam$founder, ] <- admix[ indexes, ]

    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, as indexes of the current admix matrix (precalculated)
        parents_j <- c( fam$pati[ i ], fam$mati[ i ] )
        # and extracts their admixture proportions, which are simply averaged
        # store directly at child location
        admix_fam[ i, ] <- colMeans( admix_fam[ parents_j, ] )
    }
    
    # return entire matrix!
    return( admix_fam )
}
