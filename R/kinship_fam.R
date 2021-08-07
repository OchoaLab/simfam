#' Calculate kinship matrix of a pedigree with structured founders
#'
#' Calculates a full kinship matrix (between all individuals in the provided pedigree FAM table) taking into account the relatedness of the founders as provided.
#' Output agrees with [kinship2::kinship()] but only when founders are unrelated/outbred (in other words, that function does not allow relatedness between founders).
#'
#' @inheritParams geno_fam
#' @param kinship The kinship matrix of the founders.
#' This matrix must have column and row names that identify each founder (matching codes in `fam$id`).
#' Individuals may be in a different order than `fam$id`.
#' Extra individuals in `kinship` but absent in `fam$id` will be silently ignored.
#' A traditional pedigree calculation would use `kinship = diag(n)/2` (plus appropriate column/row names), where `n` is the number of founders, to model unrelated and outbred founders.
#' However, if `kinship` measures the population kinship estimates between founders, the output is also a population kinship matrix (which combines the structural/ancestral and local/pedigree relatedness values into one).
#'
#' @return The kinship matrix of the entire `fam` table, taking the relatedness of the founders into account.
#' The rows and columns of this kinship matrix correspond to `fam$id` in that order.
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
#' # Kinship of the parents, here two unrelated/outbred individuals:
#' kinship <- diag(2)/2
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' colnames( kinship ) <- c('mother', 'father')
#' rownames( kinship ) <- c('mother', 'father')
#' # For a clearer example, make the father slightly inbred
#' # (a self-kinship value that exceeds 1/2):
#' kinship[2,2] <- 0.6
#'
#' # Calculate the full kinship matrix
#' kinship_all <- kinship_fam( kinship, fam )
#' 
#' # This is a 3x3 matrix with row/col names matching fam$id.
#' # The parent submatrix equals the input (reordered),
#' # but now there's relatedness to the child too
#' kinship_all
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
kinship_fam <- function( kinship, fam, missing_vals = c('', 0) ) {
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( kinship ) )
        stop( '`kinship` must be a matrix!' )
    if ( !isSymmetric( kinship ) )
        stop( '`kinship` must be a symmetric matrix!' )
    
    # ensures that `kinship` and `fam` agree, maps parent indexes
    # also loads of `fam` validations
    data <- match_fam_founders( fam, rownames( kinship ), 'kinship', 'row', missing_vals )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    n_ind_fam <- nrow( fam )
    
    # initialize output matrix
    kinship_fam <- matrix(
        NA,
        nrow = n_ind_fam,
        ncol = n_ind_fam
    )
    # copy names of individuals
    rownames( kinship_fam ) <- fam$id
    colnames( kinship_fam ) <- fam$id
    # copy founders
    # reorder/subset columns of X if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    kinship_fam[ fam$founder, fam$founder ] <- kinship[ indexes, indexes ]

    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, as indexes of the current kinship matrix (precalculated)
        p1 <- fam$pati[ i ]
        p2 <- fam$mati[ i ]
        # estimate kinship between child and everybody else
        # simple average works for everybody, including parents themselves
        kinship_i <- ( kinship_fam[ p1, ] + kinship_fam[ p2, ] ) / 2
        # copy back to big matrix
        kinship_fam[ i, ] <- kinship_i
        kinship_fam[ , i ] <- kinship_i
        # inbreeding is the kinship of the parents
        # encode as proper self kinship
        kinship_fam[ i, i ] <- ( 1 + kinship_fam[ p1, p2 ] ) / 2
    }
    
    # return entire matrix!
    return( kinship_fam )
}
