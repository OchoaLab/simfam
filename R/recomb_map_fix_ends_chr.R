#' Extrapolate and shift recombination map of one chromosome to ends
#'
#' Given an existing recombination map and a chromosome length in base pairs, extrapolates the map to ensure all positions are covered, and shifts to ensure position one in basepairs corresponds to position 0 in genetic position.
#' Recombination rates are extrapolated from the first and last 10Mb of data by default (separately per end).
#' Therefore fixes the fact that common maps start genetic position zero at base pair position `>> 1` and do not extend to ends (some SNPs from modern projects fall out of range without fixes).
#'
#' @param map A tibble with two columns: `pos` position in base pairs, and `posg` position in centiMorgans (cM).
#' @param pos_length The length of the chromosome in base pairs.
#' @param pos_delta The size of the window used to extrapolate recombination rates.
#'
#' @return The extrapolated recombination map, shifted so the first non-trivial position maps to the genetic distance expected from the extrapolated rate at the beginning, then added a first trivial position (`pos=1, posg=0`) and final basepair position at length of chromosome and expected genetic position from end extrapolation.
#'
#' @examples
#' library(tibble)
#' # create a toy recombination map with at least 10Mb at each end
#' map <- tibble(
#'     pos  = c( 3L,  15L, 100L, 120L ) * 1e6L,
#'     posg = c(  0, 10.4, 90.1,  110 )
#' )
#' # and length
#' pos_length <- 150L * 1e6L
#' 
#' # apply function!
#' map_fixed <- recomb_map_fix_ends_chr( map, pos_length )
#' # inspect
#' map_fixed
#' 
#' @seealso
#' [recomb_map_simplify_chr()] to simplify recombination maps to a desired numerical accuracy.
#'
#' @export
recomb_map_fix_ends_chr <- function( map, pos_length, pos_delta = 10000000L ) {
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( missing( pos_length ) )
        stop( '`pos_length` is required!' )
    # force pos_length to be integer now
    if ( !is.integer( pos_length ) )
        pos_length <- as.integer( pos_length )
    # further validate inputs
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
    if ( pos_length <= 0L )
        stop( '`pos_length` must be positive!' )
    # once monotonicity has been established, the rest of the checks concern the ends only
    pos_min <- map$pos[ 1L ]
    posg_min <- map$posg[ 1L ]
    m <- nrow( map )
    pos_max <- map$pos[ m ]
    posg_max <- map$posg[ m ]
    if ( pos_min <= 0L )
        stop( 'The minimum `map$pos` must be positive!' )
    if ( posg_min < 0 )
        stop( 'The minimum `map$posg` must be non-negative!' )
    if ( pos_max > pos_length )
        stop( 'The maximum `map$pos` must be equal or less than `pos_length`!' )
    
    # fix lower end first
    # our data always has pos_min > 1, but often posg_min == 0 (all but a few chrs)
    # the positions are such that those posg_min shouldn't be zero (far greater than numeric precision)
    # sensible thing to do here is to shift posg
    # let's do even better and use the first pos_delta window to see the rate and extrapolate with that
    #pos_min_delta <- pos_min + delta
    posg_min_delta <- stats::approx( map$pos, map$posg, pos_min + pos_delta )$y
    # rate in cM / base
    rate0 <- ( posg_min_delta - posg_min ) / pos_delta # (pos_min_delta - pos_min)
    # rate0 >= 0 guaranteed by monotonicity
    # so the observed min position should be at this cM:
    #posg_min_exp <- pos_min * rate0
    # lastly, shift depending on observed minimum in cM
    #posg_shift <- posg_min_exp - posg_min
    # same but all in one step
    posg_shift <- pos_min * rate0 - posg_min
    # apply shift to posg column
    map$posg <- map$posg + posg_shift
    # (the way this is coded, handles all posg_min cases, no need to treat posg_min==0 differently

    # now fix upper end
    # extract $posg again because shift modified the numbers! (m is unchanged for now)
    posg_max <- map$posg[m]
    # confirmed end of chromosome is never reached exactly, so always have to add this
    # calculate rate as before, but at other end
    posg_max_delta <- stats::approx( map$pos, map$posg, pos_max - pos_delta )$y
    # rate in cM / base
    ratem <- ( posg_max - posg_max_delta ) / pos_delta # (pos_max - pos_max_delta)
    # ratem >= 0 guaranteed by monotonicity
    # no shifting here!  just extrapolation
    posg_length <- posg_max + (pos_length - pos_max) * ratem
    
    # now add first and last rows as desired
    # map$pos will stay integer given previous data castings
    map <- rbind(
        tibble::tibble( pos = 1L, posg = 0 ),
        map,
        tibble::tibble( pos = pos_length, posg = posg_length )
    )
    # all done with both ends!
    return( map )
}
