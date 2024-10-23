#' Tidy recombination block data inherited in individuals from founders
#'
#' This function takes a compact structure of nested lists and tables describing the founder blocks of identity by descent (IBD) inherited in a pedigree with recombination, and outputs the same data organized as a single table following tidy conventions, which is easy to manipulate externally although it is also larger and more redundant.
#'
#' @param inds The output value from `recomb_map_inds`.
#' In particular, a list of named individuals, each of which is a list with two parents labeled "pat" and "mat", each of which is a list of chromosomes identified by index, each of which is a table with one row per IBD block, and at least two columns, labeled "anc" is the label of the founder individual from which this block was derived, and "pos" which is the end coordinate in basepairs of the block (the start of the block is given implicitly by the "pos" of the previous row plus 1, or 1 for the first block).
#'
#' @return A table with one row per IBD block, and the following columns: "ind" containing the label of the individual as given in the input object, "parent" equal to either "pat" or "mat", "chr" the chromosome, "start" and "end" the range of the block in basepairs, and "anc" the label of the founder individual.
#'
#' @examples
#' # simulate the ancestors of one person to 3 generations
#' obj <- fam_ancestors( 3 )
#' fam <- obj$fam
#' ids <- obj$ids
#' # initialize founders
#' founders <- recomb_init_founders( ids[[ 1 ]], recomb_map_hg38 )
#' # draw recombination breaks along pedigree, with coordinates in genetic distance (centiMorgans),
#' # with information for last generation only
#' inds <- recomb_last_gen( founders, fam, ids )
#' # map recombination break coordinates to base pairs
#' inds <- recomb_map_inds( inds, recomb_map_hg38 )
#'
#' # now that the input structure is ready, this function returns a tidy table version!
#' inds_tidy <- tidy_recomb_map_inds( inds )
#'
#' @seealso
#' [recomb_map_inds()]
#' 
#' @export
tidy_recomb_map_inds <- function( inds ) {
    # in the first pass just reorganize this data into a bigger single table that is easier to manipulate
    data <- NULL

    # https://cran.r-project.org/web/packages/dplyr/vignettes/in-packages.html
    
    # navigate structure
    for ( name_i in names( inds ) ) {
        ind <- inds[[ name_i ]]
        for ( parent in c('pat', 'mat') ) {
            par <- ind[[ parent ]]
            n_chr <- length( par )
            for ( chr_i in 1 : n_chr ) {
                # finally, this is a tibble describing the founders of origin here and their chunks
                chr <- par[[ chr_i ]]
                # we don't care about genetic positions here anymore, and rename pos to end because that's what this is
                chr <- dplyr::select( chr, end = 'pos', 'anc' )
                # initialize new column start, whose value is 1 for very first row so meh
                chr$start <- 1L
                # then calculate start positions of 2nd and on rows, as they are implied to be, shifting to be previous row's values
                if ( nrow( chr ) > 1 )
                    chr$start[ 2L : nrow( chr ) ] <- chr$end[ 1L : ( nrow( chr ) - 1L) ] + 1L
                # also add new column that marks chromosome
                chr$chr <- chr_i
                # for completeness, mark subject that contained this chunk
                chr$ind <- name_i
                chr$parent <- parent
                # append to table containing data
                data <- dplyr::bind_rows( data, chr )
            }
        }
    }

    # reorder columns now that we're here
    data <- dplyr::select( data, 'ind', 'parent', 'chr', 'start', 'end', 'anc' )
    
    return( data )
}

