recomb_haplo_chr <- function( chr, X, pos, ret_anc = FALSE ) {
    # validations
    if ( missing( chr ) )
        stop( '`chr` is required!' )
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( pos ) )
        stop( '`pos` is required!' )
    # chr validations in detail
    if ( !is.data.frame( chr ) )
        stop( '`chr` must be a data.frame (including tibble)!' )
    if ( !( 'anc' %in% names( chr ) ) )
        stop( '`chr` must have column `anc`!' )
    if ( !( 'pos' %in% names( chr ) ) )
        stop( '`chr` must have column `pos`!' )
    # X validations in detail
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    ## if ( !is.numeric( X ) )
    ##     stop( '`X` must be numeric!' )
    m_loci <- nrow( X )
    # pos validations in detail
    if ( !is.numeric( pos ) )
        stop( '`pos` must be numeric!' )
    if ( length( pos ) != m_loci )
        stop( 'length of `pos` must equal number of rows of `X`!' )

    # make sure all ancestor names in `chr` are also in `X`
    if ( is.null( colnames( X ) ) )
        stop( '`X` must have column names!' )
    indexes_anc <- match( chr$anc, colnames( X ) )
    if ( anyNA( indexes_anc ) )
        stop( 'These ancestors in `chr$anc` are missing from `colnames(X)`: ', toString( chr$anc[ is.na( indexes_anc ) ] ) )

    # get started copying data from ancestors (`X`, `pos`) to new individual (`chr`)!
    # inputs are haploid so output will be too
    # this is the new individual's haploid variant vector (technically not "genotype", but meh)
    # at this level values don't need to be numeric, they'll just be copied whatever they are, though more specifically if they are integers the output will be too, but will allow numeric in general and also other options
    xc <- vector( storage.mode( X ), m_loci )
    # can also be useful to return ancestry per position, calculate that now if requested
    if ( ret_anc )
        ancs <- vector( 'character', m_loci )

    # in keeping with end-only data, this tells us where we left off in pos vector
    i_pos_start <- 1L
    # navigate chr rows
    for ( i_chr in 1L : nrow( chr )) {
        # find greatest `pos <= chr$pos[ i_chr ]`
        # NOTE: this is end-inclusive!
        i_pos_end <- pos <= chr$pos[ i_chr ]
        # check for cases where there are no such cases
        # note i_pos_start is unchanged in those cases
        if ( !any( i_pos_end ) )
            next
        # turn logicals to indexes
        i_pos_end <- which( i_pos_end )
        # and because which returns them in order, we want the last one if there was more than one
        if ( length( i_pos_end ) > 1L )
            i_pos_end <- i_pos_end[ length( i_pos_end ) ]
        # another required condition is that end is later or equal than start
        if ( i_pos_start > i_pos_end )
            next
        # if all of these passed, we've found a segment (of length 1 or more) to copy from ancestor to descendant!
        
        # now copy data
        indexes_pos <- i_pos_start : i_pos_end
        # `indexes_anc[ i_chr ]` gives us the ancestral chromosome to process, given pre-mapped relationship between chr$anc and colnames(X)
        xc[ indexes_pos ] <- X[ indexes_pos, indexes_anc[ i_chr ] ]
        # if want ancestries, for this block all are the same so RHS is scalar
        if ( ret_anc )
            ancs[ indexes_pos ] <- chr$anc[ i_chr ]

        # after a successful copy, can now increment i_pos_start for next round
        i_pos_start <- i_pos_end + 1L
    }

    if ( ret_anc ) {
        return( list( x = xc, anc = ancs ) )
    } else 
        return( xc )
}
