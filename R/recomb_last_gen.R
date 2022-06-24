#' Draw recombination breaks for autosomes for last generation of a pedigree
#'
#' A wrapper around the more general [recomb_fam()], specialized to save memory when only the last generation is desired ([recomb_fam()] returns recombination blocks for the entire pedigree).
#' This function assumes that generations are non-overlapping (met by the output of [sim_pedigree()]), in which case each generation `g` can be drawn from generation `g-1` data only.
#' That way, only two consecutive generations need be in memory at any given time.
#' The partitioning of individuals into generations is given by the `ids` parameter (again matches the output of [sim_pedigree()]).
#'
#' @inheritParams recomb_fam
#' @param ids A list containing vectors of IDs for each generation.
#' All these IDs must be present in `fam$id`.
#' If IDs in `fam` and `ids` do not fully agree, the code processes the IDs in the intersection, which is helpful when `fam` is pruned but `ids` is the original (larger) set.
#'
#' @return The list of individuals with recombined chromosomes of the last generation (the intersection of `ids[ length(ids) ]` and `fam$id`), in the same format as `founders` above.
#' The names of this list are last-generation individuals in the order that they appear in `fam$id`.
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
#' # initialize parents with this other function
#' # simulate three chromosomes with these lengths in cM
#' lengs <- c(50, 100, 150)
#' founders <- recomb_init_founders( ids[[1]], lengs )
#'
#' # draw recombination breaks for the children
#' inds <- recomb_last_gen( founders, fam, ids )
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
recomb_last_gen <- function( founders, fam, ids, missing_vals = c('', 0) ) {
    # validate inputs
    if ( missing( founders ) )
        stop( '`founders` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( !is.list( founders ) )
        stop( '`founders` must be a list!' )
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
        founders <- recomb_fam( founders, fam_g, missing_vals )
        # subset to extract current generation now, overwrite individuals
        founders <- founders[ !indexes1 ]
    }
    # done, return the genotypes of the final generation
    return( founders )
}
