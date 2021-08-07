#' Calculate admixture matrix for last generation of a pedigree with admixture of founders
#'
#' A wrapper around the more general [admix_fam()], specialized to save memory when only the last generation is desired ([admix_fam()] returns admixture for the entire pedigree in a single matrix).
#' This function assumes that generations are non-overlapping (met by the output of [sim_pedigree()]), in which case each generation `g` can be drawn from generation `g-1` data only.
#' That way, only two consecutive generations need be in memory at any given time.
#' The partitioning of individuals into generations is given by the `ids` parameter (again matches the output of [sim_pedigree()]).
#'
#' @inheritParams admix_fam
#' @inheritParams geno_last_gen
#'
#' @return The admixture proportions matrix of the last generation (the intersection of `ids[ length(ids) ]` and `fam$id`).
#' The rows of this matrix are last-generation individuals in the order that they appear in `fam$id`.
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
#' # admixture proportions of the parents
#' admix <- rbind( c(0.3, 0.3, 0.4), c(0.5, 0.25, 0.25) )
#' # Name the parents with same codes as in `fam`
#' # (order can be different)
#' rownames( admix ) <- c('mother', 'father')
#' # name ancestries too
#' colnames( admix ) <- c('African', 'European', 'Asian')
#'
#' # calculate the admixture matrix of the children
#' admix2 <- admix_last_gen( admix, fam, ids )
#' admix2
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
admix_last_gen <- function( admix, fam, ids, missing_vals = c('', 0) ) {
    # validate inputs
    if ( missing( admix ) )
        stop( '`admix` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( !is.matrix( admix ) )
        stop( '`admix` must be a matrix!' )
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
        admix <- admix_fam( admix, fam_g, missing_vals )
        # subset to extract current generation now, overwrite admixture matrix!
        admix <- admix[ !indexes1, , drop = FALSE ]
    }
    # done, return the admixture matrix of the final generation
    return( admix )
}
