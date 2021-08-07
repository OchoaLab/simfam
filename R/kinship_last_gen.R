#' Calculate kinship matrix for last generation of a pedigree with structured founders
#'
#' A wrapper around the more general [kinship_fam()], specialized to save memory when only the last generation is desired ([kinship_fam()] returns kinship for the entire pedigree in a single matrix).
#' This function assumes that generations are non-overlapping (met by the output of [sim_pedigree()]), in which case each generation `g` can be drawn from generation `g-1` data only.
#' That way, only two consecutive generations need be in memory at any given time.
#' The partitioning of individuals into generations is given by the `ids` parameter (again matches the output of [sim_pedigree()]).
#'
#' @inheritParams kinship_fam
#' @inheritParams geno_last_gen
#'
#' @return The kinship matrix of the last generation (the intersection of `ids[ length(ids) ]` and `fam$id`).
#' The columns/rows of this matrix are last-generation individuals in the order that they appear in `fam$id`.
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
#' # calculate the kinship matrix of the children
#' kinship2 <- kinship_last_gen( kinship, fam, ids )
#' kinship2
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
kinship_last_gen <- function( kinship, fam, ids, missing_vals = c('', 0) ) {
    # validate inputs
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( !is.matrix( kinship ) )
        stop( '`kinship` must be a matrix!' )
    if ( !isSymmetric( kinship ) )
        stop( '`kinship` must be a symmetric matrix!' )
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
        kinship <- kinship_fam( kinship, fam_g, missing_vals )
        # subset to extract current generation now, overwrite kinship matrix!
        kinship <- kinship[ !indexes1, !indexes1, drop = FALSE ]
    }
    # done, return the kinship of the final generation
    return( kinship )
}
