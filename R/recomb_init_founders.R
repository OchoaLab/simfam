#' Initialize chromosome structures for founders
#'
#' This function initializes what is otherwise a tedious structure for founders, to be used for simulating recombination in a pedigree.
#' The genetic structure is trivial, in that these "founder" chromosomes are each of a single ancestral individual (none are recombined).
#'
#' @param ids The list of IDs to use for each individual
#' @param lengs The lengths of each chromosome in centiMorgans (cM).
#' If this vector is named, the output inherits these chromosome names.
#' If it is a list, it is assumed to be a recombination map (see [`recomb_map_hg`] for examples) and the desired lengths extracted automatically (taken as the last value of column `posg` of each chromosome).
#'
#' @return A named list of diploid individuals, each of which is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes (inherits names of `lengs` if present), each of which is a tibble with a single row and two columns: `posg` equals the chromosome length, and `anc` equals the ID of the individual (from `ids`) concatenated with either `_pat` or `_mat` depending on which parent it is.
#'
#' @examples
#' # version with explicit recombination lengths
#' ancs <- recomb_init_founders( c('a', 'b'), c(100, 200) )
#' ancs
#'
#' # version using genetic map (uses provided human map) from which lengths are extracted
#' ancs <- recomb_init_founders( c('a', 'b'), recomb_map_hg38 )
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
    if ( is.list( lengs ) ) {
        # extract from recombination map if that was provided
        lengs <- recomb_map_lengs( lengs, name = 'lengs' )
    }
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
                                               posg = lengs[ chr ],
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
