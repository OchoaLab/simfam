recomb_map_hap <- function( hap, map ) {
    # validations
    if ( missing( hap ) )
        stop( '`hap` is required!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( !is.list( hap ) )
        stop( '`hap` must be a list!' )
    if ( !is.list( map ) )
        stop( '`map` must be a list!' )

    # navigate chromosomes, apply processing to each
    n_chr <- length( hap )
    # access by index, as names are not required to be present
    for ( chr in 1 : n_chr ) {
        # make sure this map exists
        map_chr <- map[[ chr ]]
        if ( is.null( map_chr ) )
            stop( 'Recombination `map` missing for chr ', chr, '!' )
        # let other validations happen inside this function
        hap[[ chr ]] <- recomb_map_chr( hap[[ chr ]], map_chr )
    }
    return( hap )
}
