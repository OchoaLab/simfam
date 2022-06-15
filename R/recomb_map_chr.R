recomb_map_chr <- function( chr, map ) {
    # validations
    if ( missing( chr ) )
        stop( '`chr` is missing!' )
    if ( missing( map ) )
        stop( '`map` is missing!' )
    if ( !is.data.frame( chr ) )
        stop( '`chr` must be a data.frame (including tibble)!' )
    if ( !( 'posg' %in% names( chr ) ) )
        stop( '`chr` must have column `posg`!' )
    if ( !( 'anc' %in% names( chr ) ) ) # not used except at reordering step
        stop( '`chr` must have column `anc`!' )
    if ( !is.data.frame( map ) )
        stop( '`map` must be a data.frame (including tibble)!' )
    if ( !( 'pos' %in% names( map ) ) )
        stop( '`map` must have column `pos`!' )
    if ( !( 'posg' %in% names( map ) ) )
        stop( '`map` must have column `posg`!' )
    
    # actual mapping
    chr$pos <- recomb_map_posg( chr$posg, map )

    # it's more visually pleasing to show pos first, reorder columns!
    # first columns
    names_order <- c('pos', 'posg', 'anc')
    # keep all other columns but move them last (if any)
    names_order <- c( names_order, setdiff( names( chr ), names_order ) )
    # apply order!
    chr <- chr[ , names_order ]
    
    return( chr )
}
