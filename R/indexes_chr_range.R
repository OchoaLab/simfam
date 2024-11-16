# extracts info we want from a minimal structure that doesn't make this very easy to get
indexes_chr_range <- function( indexes_chr_ends, chr ) {
    # this is the easiest to retrieve
    # note that if chr is actually missing from structure, this will return NA for the end
    end <- indexes_chr_ends[ chr ]
    # let's make it easier outside, return NA for both ends of range if chr is missing
    if ( is.na( end ) ) return( c(NA, NA) )
    
    # this one is more tricky
    # first chromosome always starts at the first position (implied)
    start <- 1L
    if ( chr > 1L ) {
        # all subsequent chromosomes start where the previous one started plus 1
        # however, if the immediately previous chromosome was missing, we want to move backward until we find something that is non-zero; if all previous chromosomes were missing, the max part returns zero, which after +1 returns 1 as desired!
        # these are the candidates to consider, but there can be missingness here
        x <- indexes_chr_ends[ 1L : ( chr - 1L ) ]
        # if all are missing, stay with start=1, otherwise this is the value we want
        if ( !all( is.na( x ) ) )
            start <- max( x, na.rm = TRUE ) + 1L
    }
    
    # return range
    return( c( start, end ) )
}
