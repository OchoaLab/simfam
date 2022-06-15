# extracts chromosome lengths in genetic distance
recomb_map_lengs <- function( map, name = 'map' ) {
    if ( missing( map ) )
        stop( '`', name, '` is required!' )
    if ( !is.list( map ) )
        stop( '`', name, '` must be a list!' )

    # make sure every chromosome has desired column, then extract the length
    lengs <- sapply( map, function(x) {
        if ( !( 'posg' %in% names( x ) ) )
            stop( '`', name, '` has a chromosome missing column "posg"!' )
        
        return( x$posg[ nrow(x) ] )
    })
    
    return( lengs )
}
