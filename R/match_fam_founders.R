# shared by `draw_geno_fam` and `kinship_fam`
match_fam_founders <- function( fam, names_founders, name_var, name_dim ) {
    # validate inputs
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( names_founders ) )
        stop( '`names_founders` is required!' )
    if ( missing( name_var ) )
        stop( '`name` is required!' )
    if ( missing( name_dim ) )
        stop( '`name_dim` is required!' )
    
    # check here that fam has all of the required parts
    if ( !is.data.frame( fam ) )
        stop( '`fam` must be a data.frame!' )
    if ( is.null( fam$id ) )
        stop( '`fam$id` is required!' )
    if ( is.null( fam$pat ) )
        stop( '`fam$pat` is required!' )
    if ( is.null( fam$mat ) )
        stop( '`fam$mat` is required!' )
    
    # map happens via `names_founders`
    # make sure everything we need exists
    if ( is.null( names_founders ) )
        stop( '`', name_var, '` must have ', name_dim, ' names for parents!' )
    # number of founders now properly defined
    n <- length( names_founders )
    
    # let's only require that `fam` founders are present in `names_founders`
    # so we'll allow different orders, and for `names_founders` to have extra people
    # outside this function we'll want to subset those data so they agree with `fam`, let's get that mapping here too
    # founder status (i.e. won't assume founders are all in the beginning of fam, can be discontiguous)
    fam$founder <- fam$pat == 0 | fam$mat == 0 # logical
    indexes <- match( fam$id[ fam$founder ], names_founders )
    # so only error is if any of these are NA (`fam$id` of founder not in `names_founders`)
    if ( anyNA( indexes ) )
        stop( 'Founders in `fam` (missing either parent) must be present in `', name_var, '`!' )
    
    # need indexes of parents for fastest computation
    # this should be in terms of final order (i.e. fam$id)
    fam$pati <- match( fam$pat, fam$id )
    fam$mati <- match( fam$mat, fam$id )
    # make sure all non-missing parents are actually present as individuals earlier
    if ( anyNA( fam$pati[ !fam$founder ] ) )
        stop( 'All non-missing parents in `fam$pat` must be present in `fam$id`!' )
    if ( anyNA( fam$mati[ !fam$founder ] ) )
        stop( 'All non-missing parents in `fam$mat` must be present in `fam$id`!' )

    # return fam with parent indexes mapped
    # and indexes of founders in `names_founders`
    return( list( fam = fam, indexes = indexes ) )
}

