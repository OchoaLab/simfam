#' Calculate genetic positions from base pair positions and a genetic map
#'
#' Given a table of base pair positions (a data frame with chromosome and position values), and a genetic map (see [`recomb_map_hg`]), this function calculates genetic positions.
#' If genetic positions existed in input, they are overwritten.
#'
#' Base pair positions are converted to genetic positions from the provided map using linear interpolation, using [stats::approx()] with options `rule = 2` (out of range cases are set to nearest end's value) and `ties = list( 'ordered', mean )` (assume data is ordered, interpolate ties in base pair positions in map using mean of genetic positions).
#' Output will be incorrect, without throwing errors, if genetic map is not ordered.
#'
#' @param bim The table of variants, which is a data.frame/tibble with at least two columns: `chr` (must be numeric between 1 and the maximum chromosome in `map` below for map to work, otherwise ignored with a warning) and `pos` (base pair position, usually an integer).
#' If column `posg` is present, it will be overwritten in the output, otherwise it is created.
#' Additional columns may be present and will be unedited.
#' @param map The genetic map, a list of chromosomes each of which is a data.frame/tibble with columns `pos` for base pair position and `posg` for genetic position.
#'
#' @return The `bim` input with new or overwritten column `posg` of genetic positions in cM.
#' Rows with values of `chr` that are not numeric or are out of range (for given `map`) are unedited if the `posg` column was present, or assigned `NA` otherwise.
#'
#' @examples
#' # let's define a very simple table of base pair positions, with minimal information
#' library(tibble)
#' bim <- tibble(
#'   chr = c(1, 1, 2, 2),
#'   pos = c(50, 200, 30, 123) * 1000000
#' )
#' # use latest human recombination map
#' map <- recomb_map_hg38
#' 
#' # now use this function to add genetic positions column to `bim`!
#' bim <- bim_add_posg( bim, map )
#'
#' @seealso
#' [`recomb_map_hg`] for simplified human recombination maps included in this package.
#' 
#' @export
bim_add_posg <- function( bim, map ) {
    # validate inputs
    if ( missing( bim ) )
        stop( '`bim` is missing!' )
    if ( missing( map ) )
        stop( '`map` is missing!' )
    if ( !is.data.frame( bim ) )
        stop( '`bim` must be a data.frame (including tibble)!' )
    if ( !( 'chr' %in% names( bim ) ) )
        stop( '`bim` must have column `chr`!' )
    if ( !( 'pos' %in% names( bim ) ) )
        stop( '`bim` must have column `pos`!' )
    if ( !is.list( map ) )
        stop( '`map` must be a list!' )

    # initialize column to edit, if it's missing
    if ( !( 'posg' %in% names( bim ) ) )
        bim$posg <- NA
        
    # have to linearly interpolate each chromosome separately

    # map is a list of tibbles, each chromosome is one element of the list
    # this will tell us if data is in range or not
    n_chr <- length( map )
    # expected chromosomes
    chrs_exp <- 1L : n_chr
    
    # navigate only chromosomes on bim table
    chrs <- unique( bim$chr )
    for ( chr in chrs ) {
        # skip with warning if not numeric in desired range
        if ( ! chr %in% chrs_exp ) {
            warning( 'Chromosome `', chr, '` in `bim` input is not expected numeric between 1 and ', n_chr, ' (length of `map`)!  Skipping!' )
            next
        }
        # cast chr as numeric or lookup fails (map[[1]] works but map[['1']] fails because elements aren't named like that!)
        # this should be a tibble
        map_chr <- map[[ as.numeric( chr ) ]]
        if ( !is.data.frame( map_chr ) )
            stop( '`map[[', chr, ']]` must be a data.frame (including tibble)!' )
        if ( !( 'pos' %in% names( map_chr ) ) )
            stop( '`map[[', chr, ']]` must have column `pos`!' )
        if ( !( 'posg' %in% names( map_chr ) ) )
            stop( '`map[[', chr, ']]` must have column `posg`!' )

        # get table subset on bim
        indexes <- bim$chr == chr
        chr_pos <- bim$pos[ indexes ]
        # perform desired interpolation
        out <- stats::approx( map_chr$pos, map_chr$posg, chr_pos, rule = 2, ties = list( "ordered", mean ) )
        # overwrite previous posg for this chromosome!
        bim$posg[ indexes ] <- out$y
    }

    # all done, return edited `bim` tibble
    return( bim )
}
