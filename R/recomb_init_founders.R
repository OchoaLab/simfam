#' Initialize chromosome strctures for founders
#'
#' This function initializes what is otherwise a tedious structure for founders, to be used for simulating recombination in a pedigree.
#' The genetic structure is trivial, in that these "founder" chromosomes are each of a single ancestral individual (none are recombined).
#'
#' @param ids The list of IDs to use for each individual
#' @param lengs The lengths of each chromosome in centi-Morgans (cM).
#' If this vector is named, the output inherits these chromosome names.
#'
#' @return A named list of diploid individuals, each of which is a list with two haploid individuals named "pat" and "mat", each of which is a list of chromosomes (inherits names of `lengs` if present), each of which is a tibble with a single row and two columns: "end" equals the chromosome length, and "anc" equals the ID of the individual (from `ids`) concatenated with either "_pat" or "_mat" depending on which parent it is.
#'
#' @examples
#' ancs <- recomb_init_founders( c('a', 'b'), c(100, 200) )
#' ancs
#'
#' @seealso
#' [recomb_fam()] to simulate recombination across a pedigree using the founders initialized here.
#' 
#' @export
recomb_init_founders <- function( ids, lengs ) {
    # validate inputs
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( missing( lengs ) )
        stop( '`lengs` is required!' )
    # `ids` can be almost anything, no need to test further
    # validate lengs
    if ( !is.numeric( lengs ) )
        stop( '`lengs` must be numeric!' )
    if ( any( lengs <= 0 ) )
        stop( '`lengs` must be positive!' )
    n_chr <- length( lengs )
    
    # initialize list
    inds <- vector('list', length( ids ) )
    names( inds ) <- ids
    parents <- c('pat', 'mat')
    # loop through list, initializing each parent
    for ( id in ids ) {
        ind <- vector( 'list', 2 )
        names( ind ) <- parents
        for ( parent in parents ) {
            # ID to use for ancestor
            id_par <- paste0( id, '_', parent )
            # output list of chromosomes
            hap <- vector( 'list', n_chr )
            for ( chr in 1 : n_chr ) {
                # output tibble
                hap[[ chr ]] <- tibble::tibble(
                                               end = lengs[ chr ],
                                               anc = id_par
                                           )
            }
            # pass chr names if available
            if ( !is.null( names( lengs ) ) )
                names( hap ) <- names( lengs )
            # store in bigger list
            ind[[ parent ]] <- hap
        }
        # store in outermost list
        inds[[ id ]] <- ind
    }
    # all done, return big list!
    return( inds )
}
