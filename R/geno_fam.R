#' Draw random genotypes on a pedigree with known founder genotypes
#'
#' Constructs a random genotype matrix (for all individuals in the provided pedigree FAM table) starting from the genotype matrix of the founders as provided.
#'
#' @param X The genotype matrix of the founders (loci along rows, individuals along columns).
#' This matrix must have column names that identify each founder (matching codes in `fam$id`).
#' Individuals may be in a different order than `fam$id`.
#' Extra individuals in `admix` but absent in `fam$id` will be silently ignored.
#' All values should be in `c(0L, 1L, 2L)`; for speed, this code does not check that `X` is valid (i.e. fractional values between 0 and 2 may not cause errors).
#' @param fam The pedigree data.frame, in plink FAM format.
#' Only columns `id`, `pat`, and `mat` are required.
#' `id` must be unique and non-missing.
#' Founders must be present, and their `pat` and `mat` values must be missing (see below).
#' Non-founders must have both their parents be non-missing.
#' Parents must appear earlier than their children in the table.
#' @param missing_vals The list of ID values treated as missing.
#' `NA` is always treated as missing.
#' By default, the empty string ('') and zero (0) are also treated as missing (remove values from here if this is a problem).
#'
#' @return The random genotype matrix of the entire `fam` table, starting from the genotypes of the founders.
#' The columns of this matrix correspond to `fam$id` in that order.
#' The rows (loci) are the same as in the input `X`.
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
#' # genotypes of the parents at 4 loci
#' X <- cbind( c(1, 2, 0, 2), c(0, 2, 2, 1) )
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' colnames( X ) <- c('mother', 'father')
#' # name loci too
#' rownames( X ) <- paste0( 'rs', 1:4 )
#'
#' # Draw the full genotype matrix
#' X_all <- geno_fam( X, fam )
#' 
#' # This is a 4x3 matrix with column names matching fam$id.
#' # The parent submatrix equals the input (reordered),
#' # but now there's random genotypes for the child too
#' X_all
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
geno_fam <- function( X, fam, missing_vals = c('', 0) ) {
    # validate inputs
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )

    # ensures that `X` and `fam` agree, maps parent indexes
    data <- match_fam_founders( fam, colnames( X ), 'X', 'column', missing_vals )
    fam <- data$fam # has new columns: founder, pati, mati
    
    # set variables up for passing to C++, and validate stuff
    # NOTE: these are R indexes, need to turn into zero-based C++ indexes!
    indexes <- data$indexes - 1 # this is `i_founder_in` in the C++ function, everything else is unambiguous
    i_founder_out <- which( fam$founder ) - 1
    i_child <- which( !fam$founder ) - 1
    i_pat <- fam$pati[ !fam$founder ] - 1
    i_mat <- fam$mati[ !fam$founder ] - 1
    
    # create the output matrix using C++!
    X_fam <- geno_fam_cpp( X, indexes, i_founder_out, i_child, i_pat, i_mat )
    
    # copy names of individuals and loci
    colnames( X_fam ) <- fam$id
    rownames( X_fam ) <- rownames( X ) # may be NULL, meh
    
    # done, return genotype matrix of all FAM
    return( X_fam )
}

