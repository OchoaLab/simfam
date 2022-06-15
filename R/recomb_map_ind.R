recomb_map_ind <- function( ind, map ) {
    # validations
    if ( missing( ind ) )
        stop( '`ind` is required!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( !is.list( ind ) )
        stop( '`ind` must be a list!' )
    if ( !is.list( map ) )
        stop( '`map` must be a list!' )
    if ( !( 'pat' %in% names( ind ) ) )
        stop( '`ind` must have named element "pat"!' )
    if ( !( 'mat' %in% names( ind ) ) )
        stop( '`ind` must have named element "mat"!' )

    # just apply to both parental copies
    ind$pat <- recomb_map_hap( ind$pat, map )
    ind$mat <- recomb_map_hap( ind$mat, map )

    return( ind )
}
