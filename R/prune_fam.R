# make a `fam` pruner that eliminates individuals that are not ancestors of a certain focal population
prune_fam <- function( fam, ids ) {
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
        # also exclude missing parents (0)
        parents <- setdiff( parents, c(ids_keep, 0) )
        # add to IDs to keep (since we already subtracted shared people, simple concatenation keeps list unique)
        ids_keep <- c( ids_keep, parents )
        # now look for their parents!
        # NOTE: eventually we hit founders, whose parents are all "0", which we've removed, so list will be empty
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
