#' Remove non-ancestors of a set of individuals from pedigree
#'
#' This function accepts an input pedigree and a list of individuals of interest, and returns the subset of the pedigree including only the individuals of interest and their direct ancestors.
#' This is useful in simulations, to avoid modeling/drawing genotypes of individuals without descendants in the last generation.
#'
#' @param fam The pedigree data.frame, in plink FAM format.
#' Only columns `id`, `pat`, and `mat` are required.
#' `id` must be unique and non-missing.
#' Founders must be present, and their `pat` and `mat` values must be missing (see below).
#' Non-founders must have both their parents be non-missing.
#' Parents must appear earlier than their children in the table.
#' @param ids The list of individuals of interest, whose ancestors we want to keep.
#' All must be present in `fam$id`.
#' @param missing_vals The list of ID values treated as missing.
#' `NA` is always treated as missing.
#' By default, the empty string ('') and zero (0) are also treated as missing (remove values from here if this is a problem).
#' 
#' @return The filtered FAM table with non-ancestors of `ids` excluded.
#' IDs that are `NA`-equivalent (see `missing_vals`) will be mapped to `NA`.
#'
#' @examples
#' # construct a family with three founders, but one "bob" has no descendants
#' library(tibble)
#' fam <- tibble(
#'     id  = c('mom', 'dad', 'bob', 'child'),
#'     pat = c(   NA,    NA,    NA,   'dad'),
#'     mat = c(   NA,    NA,    NA,   'mom')
#' )
#' # only want 'child' and its ancestors
#' ids <- 'child'
#' fam2 <- prune_fam( fam, ids )
#' # the filtered pedigree has "bob" removed:
#' fam2
#' 
#' @export
prune_fam <- function( fam, ids, missing_vals = c('', 0) ) {
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( !is.data.frame( fam ) )
        stop( '`fam` must be a data.frame!' )
    if ( length( ids ) == 0 )
        stop( '`ids` must have non-zero length!' )
    if ( is.null( fam$id ) )
        stop( '`fam$id` is required!' )
    if ( is.null( fam$pat ) )
        stop( '`fam$pat` is required!' )
    if ( is.null( fam$mat ) )
        stop( '`fam$mat` is required!' )
    
    # normalize missing values for individuals (id, pat, mat)
    # the only obvious value that is missing no matter what we chose is NA
    # these others are treated as NA-equivalent
    for ( missing_val in missing_vals ) {
        fam$id[ fam$id == missing_val ] <- NA
        fam$pat[ fam$pat == missing_val ] <- NA
        fam$mat[ fam$mat == missing_val ] <- NA
    }
    
    # in case input is weird, this will behave
    ids <- unique( ids )
    # ids should be in FAM!
    if ( !all( ids %in% fam$id ) )
        stop( 'At least some `ids` are not in `fam$id`!' )
    # now start moving backward recursively
    ids_keep <- ids # obviously want these guys
    while ( length( ids ) > 0 ) {
        # identify rows of current individuals
        indexes <- fam$id %in% ids
        # collect their parents, make unique (siblings share some)
        parents <- unique( c( fam$pat[ indexes ], fam$mat[ indexes ] ) )
        # and never include people we've already analyzed (relevant for complex pedigrees with individiuals mating across generations)
        # also exclude missing parents (NA)
        parents <- setdiff( parents, c(ids_keep, NA) )
        # add to IDs to keep (since we already subtracted shared people, simple concatenation keeps list unique)
        ids_keep <- c( ids_keep, parents )
        # now look for their parents!
        # NOTE: eventually we hit founders, whose parents are all missing, which we've removed, so list will be empty
        ids <- parents
    }
    # by construction, ids_keep is non-redundant and has no missing IDs
    # finally, find the rows of FAM to keep
    indexes <- fam$id %in% ids_keep
    # apply filter
    fam <- fam[ indexes, ]
    # done!
    return( fam )
}
