#' Simplify recombination map of one chromosome to a desired numerical precision
#'
#' Given an input recombination map, this function iteratively removes rows that can be interpolated to less than a given error `tol`.
#' This is a heuristic that works very well in practice, resulting in average interpolation errors well below `tol`, and maximum final errors no greater than `3 * tol` in our internal benchmarks (expected in extremely concave or convex regions of the map; final errors are rarely above `tol` with few exceptions).
#' 
#' This function reduces recombination map sizes drastically, in order to include them in packages, and also makes linear interpolation faster.
#' This simplification operation can be justified as the precision of many existing maps is both limited and overstated, and a high accuracy is not needed for simulations with many other approximations in place.
#'
#' @param map A tibble with two columns: `pos` position in base pairs, and `posg` position in centiMorgans (cM).
#' @param tol Tolerance of interpolation errors, in cM.
#'
#' @return The recombination map with rows (positions) removed (if they are interpolated with errors below `tol` in most cases).
#'
#' @examples
#' library(tibble)
#' # create a toy recombination map to simplify
#' # in this case all middle rows can be interpolated from the ends with practically no error
#' map <- tibble(
#'     pos  = c(  1L, 1e6L, 2e6L, 3e6L ),
#'     posg = c( 0.0,  1.0,  2.0,  3.0 )
#' )
#' 
#' # simplify map!
#' map_simple <- recomb_map_simplify_chr( map )
#' # inspect
#' map_simple
#'
#' @seealso
#' [recomb_map_fix_ends_chr()] to shift and extrapolate recombination map to ends of chromosome.
#' 
#' @export 
recomb_map_simplify_chr <- function( map, tol = 0.1 ) {
    if ( missing( map ) )
        stop( '`map` is required!' )
    # go all the way because this isn't done often, fully validate monotonicity/etc
    if ( !is.data.frame( map ) )
        stop( '`map` must be a data.frame (including tibble)!' )
    if ( !( 'pos' %in% names( map ) ) )
        stop( '`map` must have column named `pos`!' )
    if ( !( 'posg' %in% names( map ) ) )
        stop( '`map` must have column named `posg`!' )
    # force map$pos to be integer now
    if ( !is.integer( map$pos ) )
        map$pos <- as.integer( map$pos )
    if ( !all( diff( map$pos ) > 0L ) )
        stop( '`map$pos` must be strictly monotonically increasing!' )
    if ( !all( diff( map$posg ) >= 0 ) ) # allow for ties here for now, raw data has these ties
        stop( '`map$posg` must be monotonically increasing!' )
    # once monotonicity has been established, the rest of the checks concern the ends only
    pos_min <- map$pos[ 1L ]
    posg_min <- map$posg[ 1L ]
    if ( pos_min <= 0L )
        stop( 'The minimum `map$pos` must be positive!' )
    if ( posg_min < 0 )
        stop( 'The minimum `map$posg` must be non-negative!' )

    # start an infinite loop
    while ( TRUE ) {
        # for all points except both extremes, calculate how far they are from the linear interpolation of their two neighbors
        # as we do iterative removals, `m` gets smaller each time!
        m <- nrow( map )
        # need at least three points for this to be non-trivial
        if ( m < 3L )
            return( map )
        # copy these vectors down
        x <- map$pos
        y <- map$posg
        # can compute desired errors with diff!
        dx1 <- diff( x )[ -(m-1L) ] # final value not used
        dy1 <- diff( y )[ -(m-1L) ]
        dx2 <- diff( x, lag = 2L )
        dy2 <- diff( y, lag = 2L )
        # `dy2 / dx2` is slope from extrema
        # `dy2 / dx2 * dx1` is prediction of `dy1`
        errors <- abs( dy2 / dx2 * dx1 - dy1 )
        indexes <- errors < tol
        if ( !any( indexes ) )
            return( map ) # done! and end infinite loop
        # else we have at least one removal!
        
        # NOTE: consecutive removals must not happen because final error can double! (results on paper)
        # will only need subset of these errors now
        errors <- errors[ indexes ]
        # and these as actual indexes
        indexes <- which( indexes )
        # what seems most foolproof to me is to go by order of errors, approve the smallest ones, and not bigger ones if they are next to small ones
        # default order is increasing, so perfect
        indexes_order <- indexes[ order( errors ) ]
        # now navigate and edit as needed
        # do it this awkward way because `indexes_order` might get shorter (elements will be deleted)
        # (last element of `indexes_order` doesn't delete anything else, so no need to scan it)
        k <- 1L
        while ( k < length( indexes_order ) ) {
            j <- indexes_order[ k ]
            k_rm <- c()
            # look for cases that are worse-ranking than j and are its neighbors
            for ( k2 in ( k + 1L ) : length( indexes_order ) ) {
                if ( abs( indexes_order[ k2 ] - j ) == 1L )
                    k_rm <- c( k_rm, k2 )
            }
            # apply removals if any
            if ( length( k_rm ) > 0 )
                indexes_order <- indexes_order[ -k_rm ]
            # increment for next iteration
            k <- k + 1L
        }
        # indexes left are the ones we want to remove
        # they were all off by 1 because they came from `errors`, which starts at 2, so fix that now
        indexes_order <- indexes_order + 1L
        ## # DEBUGGING
        ## message( 'Removals: ', toString( indexes_order ) )
        # apply removal (there is still at least one)
        map <- map[ -indexes_order, ]
        # now attempt another iteration!
    }
    
    # this should never happen but meh
    return( map )
}
