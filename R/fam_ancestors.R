#' Construct an ancestors-only pedigree for one person `G`-generations deep
#'
#' Creates an idealized pedigree listing all ancestors of one individual starting from `G` generations ago, without inbreeding (a binary tree).
#' IDs are automatically generated strings indicating generation and individual number within generation.
#' Useful for simple simulations of individuals with explicit ancestors.
#'
#' @param G The desired number of generations.
#' `G=1` returns a trivial pedigree with a single individual; `G=2` an individual and its two parents; `G=3` an individual, its parents and grandparents, etc.
#'
#' @return A list with two named elements:
#' - `fam`: a tibble describing the pedigree, with the following columns
#'   - `id`: The ID of each individual, a string in the format "g-i" joining with a dash the generation number ("g", numbered backward in time) and the individual number within the generation ("i").
#'   - `pat`: The paternal ID.  For individual "g-i" parent is (g+1)"-"(2*i-1), except for last generation it is `NA` (their parents are missing).
#'   - `mat`: The maternal ID.  For individual "g-i" parent is (g+1)"-"(2*i), except for last generation it is `NA` (their parents are missing).
#'   - `sex`: 1 (male) for all odd-numbered individuals, 2 (female) for even-numbered individuals, consistent with pedigree structure.  Side-effect is first-generation individual ("1-1") is always male (edit afterwards as desired).
#' - `ids`: A list containing vectors of IDs separated by generation, but here starting from the last generation (highest "g"), to be consistent with output of [sim_pedigree()] and the expected input of all `*_last_gen` functions.
#'
#' @examples
#' # construct the 8-generation ancestor tree of one individual:
#' data <- fam_ancestors( 8 )
#' # this is the pedigree
#' fam <- data$fam
#' # and this is the handy list of IDs by discrete generation,
#' # used by `*_last_gen` functions to reduce memory usage
#' ids <- data$ids
#'
#' @seealso
#' [sim_pedigree()] to simulate a random pedigree with a given number of generations, generation sizes, and other parameters.
#'
#' @export
fam_ancestors <- function( G ) {
    if ( missing( G ) )
        stop( 'The number of generations `G` is required!' )
    if ( G < 1L )
        stop( '`G` must be at least 1!' )
    
    # initialize vectors of interest
    # this is the known length for the sum of powers of two, adjusting for ending at G-1 instead of G
    # https://math.stackexchange.com/questions/1990137/the-idea-behind-the-sum-of-powers-of-2
    length_G <- as.integer( 2L^G ) - 1L
    id <- vector( 'character', length_G )
    pat <- vector( 'character', length_G )
    mat <- vector( 'character', length_G )
    sex <- vector( 'integer', length_G )
    # create discrete generations vector common in other functions
    ids <- vector( 'list', G )
    
    # here g counts generations backwards in time
    for ( g in 1L : G ) {
        # number of individuals in this generation
        n_g <- as.integer( 2L^( g - 1L ) )
        inds <- 1L : n_g
        # also set indexes for where to place them
        # again, since the sums of powers of two is a power of two, this odd trick works for the shift in indexes:
        indexes <- n_g - 1L + inds
        # IDs of individuals
        # copy from list of generations and fam table
        ids[[ g ]] <- paste0( g, '-', inds )
        id[ indexes ] <- ids[[ g ]]
        # IDs of parents, by sex
        if ( g == G ) {
            # for last generation parents must be missing
            pat[ indexes ] <- NA
            mat[ indexes ] <- NA
        } else {
            # usual case, for all but the last generation
            pat[ indexes ] <- paste0( g + 1L, '-', 2L * inds - 1L )
            mat[ indexes ] <- paste0( g + 1L, '-', 2L * inds )
        }
        # for completeness, let's assign sex
        # all odd parents are male, even are female
        # only final child may be anything, keeping the pattern below forces to male
        sex[ indexes ] <- rep_len( c(1L, 2L), n_g )
    }
    
    # put it all together in a tibble
    fam <- tibble::tibble( id = id, pat = pat, mat = mat, sex = sex )
    # NOTE: above code listed children before their parents, reverse so the usual ordering/convention holds
    fam <- fam[ nrow( fam ) : 1L, ]
    # same for IDs, though code is more straightforward
    ids <- rev( ids )
    # return data
    return( list( fam = fam, ids = ids ) )
}

