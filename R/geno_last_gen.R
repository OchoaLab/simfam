#' Draw random genotypes for last generation of a pedigree with known founder genotypes
#'
#' A wrapper around the more general [geno_fam()], specialized to save memory when only the last generation is desired ([geno_fam()] returns genotypes for the entire pedigree in a single matrix).
#' This function assumes that generations are non-overlapping (met by the output of [sim_pedigree()]), in which case each generation `g` can be drawn from generation `g-1` data only.
#' That way, only two consecutive generations need be in memory at any given time.
#' The partitioning of individuals into generations is given by the `ids` parameter (again matches the output of [sim_pedigree()]).
#'
#' @inheritParams geno_fam
#' @param ids A list containing vectors of IDs for each generation.
#' All these IDs must be present in `fam$id`.
#' If IDs in `fam` and `ids` do not fully agree, the code processes the IDs in the intersection, which is helpful when `fam` is pruned but `ids` is the original (larger) set.
#'
#' @return The random genotype matrix of the last generation (the intersection of `ids[ length(ids) ]` and `fam$id`).
#' The columns of this matrix are last-generation individuals in the order that they appear in `fam$id`.
#' The rows (loci) are the same as in the input `X`.
#'
#' @examples
#' # A small pedigree, two parents and two children.
#' # A minimal fam table with the three required columns.
#' # Note "mother" and "father" have missing parent IDs, while children do not
#' library(tibble)
#' fam <- tibble(
#'   id = c('father', 'mother', 'child', 'sib'),
#'   pat = c(NA, NA, 'father', 'father'),
#'   mat = c(NA, NA, 'mother', 'mother')
#' )
#' # need an `ids` list separating the generations
#' ids <- list( c('father', 'mother'), c('child', 'sib') )
#' 
#' # genotypes of the parents at 4 loci
#' X <- cbind( c(1, 2, 0, 2), c(0, 2, 2, 1) )
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' colnames( X ) <- c('mother', 'father')
#' # name loci too
#' rownames( X ) <- paste0( 'rs', 1:4 )
#'
#' # Draw the genotype matrix of the children
#' X2 <- geno_last_gen( X, fam, ids )
#' X2
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
geno_last_gen <- function( X, fam, ids, missing_vals = c('', 0) ) {
    # validate inputs
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    if ( !is.list( ids ) )
        stop( '`ids` must be a list!' )

    # draw generations sequentially
    G <- length( ids )
    for ( g in 2 : G ) {
        ids1 <- ids[[ g-1 ]] # parents (exist)
        ids2 <- ids[[ g ]] # children (to draw)
        # subset fam to contain only previous (parents) and current (children) generations
        fam_g <- fam[ fam$id %in% c(ids1, ids2), ]
        # parents must have their own parents set to NA (to be treated as founders)
        indexes1 <- fam_g$id %in% ids1
        fam_g$pat[ indexes1 ] <- NA
        fam_g$mat[ indexes1 ] <- NA
        # now run through this more general function!
        X <- geno_fam( X, fam_g, missing_vals )
        # subset to extract current generation now, overwrite genotype matrix!
        X <- X[ , !indexes1, drop = FALSE ]
    }
    # done, return the genotypes of the final generation
    return( X )
}
