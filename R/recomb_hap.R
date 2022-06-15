# parents are assumed to be a list of data.frames, one per chr, with 2 columns: posg, anc
# (start is 0 for first chunk, end from previous row otherwise)
# mat/pat are symmetric (no distinction in practice)
recomb_hap <- function( mat, pat ) {
    # validate data
    if ( missing( mat ) )
        stop( '`mat` is required!' )
    if ( missing( pat ) )
        stop( '`pat` is required!' )
    if ( !is.list( mat ) )
        stop( '`mat` must be a list' )
    if ( !is.list( pat ) )
        stop( '`pat` must be a list' )
    # make sure chromosomes agree
    n_chr <- length( mat )
    if ( n_chr != length( pat ) )
        stop( 'Number of chromosomes in `mat` and `pat` disagree!' )
    # names too? if any are present...
    if ( !is.null( names( mat ) ) || !is.null( names( pat ) ) ) {
        if ( is.null( names( mat ) ) || is.null( names( pat ) ) )
            stop( 'One parent has names for chromosome and other does not!  Names must match or both be missing!' )
        # if we're here, then both parents have names
        if ( any( names( mat ) != names( pat ) ) )
            stop( 'Names for chromosomes of `mat` and `pat` disagree!' )
    }

    # start generating child, looping through chromosomes
    child <- vector( 'list', n_chr )
    names( child ) <- names( mat )

    # access by index, as names are not required to be present
    for ( chr in 1 : n_chr ) {
        # copy down these values
        chr_mat <- mat[[ chr ]]
        chr_pat <- pat[[ chr ]]
        # to draw breaks, need chr length
        # get from mat, will validate consistency with pat inside `recomb_chr`
        # minimal checks for this purpose
        if ( !is.data.frame( chr_mat ) )
            stop( '`mat[[', chr, ']]` must be a data.frame (including tibble)!' )
        if ( !( 'posg' %in% names( chr_mat ) ) )
            stop( '`mat[[', chr, ']]` must have column "posg"!' )
        chr_leng <- max( chr_mat$posg )
        # draw random chromosome breaks now!
        breaks <- recomb_breaks( chr_leng )
        # construct the ancestries of all chunks of child's chromosome now!
        child[[ chr ]] <- recomb_chr( breaks, chr_mat, chr_pat )
    }
    
    # all done!  return child info only
    return( child )
}
