#' Map recombination breaks from genetic positions to base pair coordinates
#'
#' Given a list of individuals with recombination breaks given in genetic distance (such as the output of [recomb_fam()]), and a genetic map (see [`recomb_map_hg`]), this function determines all positions in base pair coordinates.
#' If base pair positions existed in input, they are overwritten.
#'
#' Genetic positions are converted to base pair positions from the provided map using linear interpolation, using [stats::approx()] with options `rule = 2` (out of range cases are set to nearest end's value) and `ties = list( 'ordered', mean )` (assume data is ordered, interpolate ties in genetic distance in map using mean of base pair positions).
#' Output will be incorrect, without throwing errors, if genetic map is not ordered.
#' Base pair positions are rounded to integers.
#'
#' @param inds The list of individuals, each of which is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes (always identified by number, but may also be named arbitrarily), each of which is a data.frame/tibble with implicit ranges (`posg` is end coordinates in cM; start is the end of the previous block, zero for the first block) and ancestors `anc` as strings.
#' @param map The genetic map, a list of chromosomes each of which is a data.frame/tibble with columns `pos` for base pair position and `posg` for genetic position.
#'
#' @return The input list of individuals, with each chromosome added column `pos` corresponding to end coordinate in base pairs.
#' Each chromosome has columns reordered so `pos`, `posg`, and `anc` appear first, and any additional columns appear afterwards.
#'
#' @examples
#' # Lengthy code creates individuals with recombination data to map
#' # The smallest pedigree, two parents and a child (minimal fam table).
#' library(tibble)
#' fam <- tibble(
#'   id = c('father', 'mother', 'child'),
#'   pat = c(NA, NA, 'father'),
#'   mat = c(NA, NA, 'mother')
#' )
#' # use latest human recombination map, but just first two chrs to keep this example fast
#' map <- recomb_map_hg38[ 1:2 ]
#' # initialize parents with this other function
#' founders <- recomb_init_founders( c('father', 'mother'), map )
#' # draw recombination breaks for child
#' inds <- recomb_fam( founders, fam )
#' 
#' # now use this function to add base pair coordinates for recombination breaks!
#' inds <- recomb_map_inds( inds, map )
#'
#' @seealso
#' [recomb_fam()] for drawing recombination breaks of individuals from a pedigree.
#'
#' [`recomb_map_hg`] for simplified human recombination maps included in this package.
#' 
#' @export
recomb_map_inds <- function( inds, map ) {
    # validations
    if ( missing( inds ) )
        stop( '`inds` is required!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( !is.list( inds ) )
        stop( '`inds` must be a list!' )
    if ( !is.list( map ) )
        stop( '`map` must be a list!' )

    # just apply to every individual
    inds <- lapply( inds, recomb_map_ind, map )

    return( inds )
}
