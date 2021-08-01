draw_num_children_per_fam <- function( n_fam, pop_size, children_min = 1L ) {
    if ( missing( n_fam ) )
        stop( '`n_fam` is required!' )
    if ( missing( pop_size ) )
        stop( '`pop_size` is required!' )
    # to keep deltas as integers, make sure these are integers
    pop_size <- as.integer( pop_size )
    n_fam <- as.integer( n_fam )
    children_min <- as.integer( children_min )
    if ( children_min < 0 )
        stop( '`children_min` must be non-negative!' )
    
    # average number of children per family
    avg_children_per_fam <- pop_size / n_fam
    
    # this number must be larger than the minimum number of children
    # (Poisson doesn't do negative lambdas)
    if ( avg_children_per_fam < children_min )
        stop( '`children_min` exceeds the average number of children per family in this generation!  Either decrease `children_min` (recommended) or increase `pop_size` for this generation relative to the previous generation (tricky because it also depends on number of families, which is random).' )
    
    # draw from Poisson
    # since Poisson produces zeroes, which are very unproductive in our case, let's allow forcing a minimum number in all cases, and take that away from the mean
    # NOTE: these are integers because `rpois` returns integers and `children_min` was coerced into int
    children_per_fam <- children_min + stats::rpois( n_fam, avg_children_per_fam - children_min )
    
    # this count may be off from target by small amounts
    # fix counts iteratively respecting minimum fam sizes
    # NOTE: `delta` is also an int because everything here is an int (by coercion or otherwise)
    delta <- pop_size - sum( children_per_fam )
    # enter loop in case `delta` is too large and it takes several of these steps to fix (single increments each family)
    while ( delta != 0L ) {
        if ( delta > 0L ) {
            # `children_per_fam` should be incremented
            # select `delta` random families to increment each by 1
            if ( n_fam >= delta ) {
                indexes <- sample( n_fam, delta )
            } else {
                # unlikely that `delta > n_fam`, but just in case, increment all families and do another round
                indexes <- 1 : n_fam
            }
            children_per_fam[ indexes ] <- children_per_fam[ indexes ] + 1L
        } else if ( delta < 0L ) {
            delta <- -delta # turn positive for simplicity
            # `children_per_fam` should be decremented
            # here we must respect the minimum number of children, so we want to decrement larger families only
            indexes <- which( children_per_fam > children_min )
            # select a random subset out of those
            if ( length( indexes ) >= delta )
                indexes <- sample( indexes, delta )
            # else (if delta > length( indexes )) just use all indexes to decrement and do another round
            children_per_fam[ indexes ] <- children_per_fam[ indexes ] - 1L
        }
        # recalculate `delta` for next round if needed
        delta <- pop_size - sum( children_per_fam )
    }
    return( children_per_fam )
}
