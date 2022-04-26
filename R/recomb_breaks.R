# input is total genetic distance of one chromosome
# output is vector of recombination locations in genetic distance
# do not include trivial values (0 and end)
recomb_breaks <- function( leng, mean = 100 ) {
    # validate input
    if ( missing( leng ) )
        stop( '`leng` is required!' )
    if ( !is.numeric( leng ) )
        stop( '`leng` must be numeric!' )
    if ( length( leng ) != 1L )
        stop( '`leng` must be a scalar!' )
    if ( leng <= 0 )
        stop( '`leng` must be positive!' )
    
    # initialize desired list of breaks
    # start empty but non-NULL.
    # this gets returned if there were no recombinations
    breaks <- numeric(0)

    # total length of recombined chunks so far
    # starts at zero, updated in loop
    leng_done <- 0
    # iteratively cut chromosome into bits with expected exponential length
    while ( TRUE ) { # leng_done < leng ) { 
        # draw a length from distribution
        x <- stats::rexp( 1, rate = 1 / mean )
        # don't exceed length remaining
        if ( x >= leng - leng_done ) {
            # if we've passed end, don't update x/breaks/leng_done, just return
            return( breaks )
        } else {
            # update length for next iteration
            leng_done <- leng_done + x
            # update desired output vector
            breaks <- c( breaks, leng_done )
            # note this guarantees that on next iteration `leng_done < leng` as desired
        }
    }
}

