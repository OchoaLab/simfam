# posg is a vector of positions in genetic distance
# map is a data.frame/tibble with names "pos" and "posg"
recomb_map_posg <- function( posg, map ) {
    # validations
    if ( missing( posg ) )
        stop( '`posg` is missing!' )
    if ( missing( map ) )
        stop( '`map` is missing!' )
    if ( !is.numeric( posg ) )
        stop( '`posg` must be numeric!' )
    if ( !is.data.frame( map ) )
        stop( '`map` must be a data.frame (including tibble)!' )
    if ( !( 'pos' %in% names( map ) ) )
        stop( '`map` must have column `pos`!' )
    if ( !( 'posg' %in% names( map ) ) )
        stop( '`map` must have column `posg`!' )

    # wrapper for this mostly!
    # `rule = 2` returns value at closest extreme outside of range (never NA)
    # this `ties` says posg is ordered (don't check), average ties (avoids warnings by specifying explicitly)
    out <- stats::approx( map$posg, map$pos, posg, rule = 2, ties = list( 'ordered', mean ) )
    # these are the values in basepairs
    pos <- out$y
    
    # because linear interpolation was used these are not integers!
    # I think it doesn't matter, but keep it in mind!
    # this rounds and sets to integers
    pos <- as.integer( round( pos ) )
    
    return( pos )
}
