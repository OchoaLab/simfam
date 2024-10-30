#' Identify IBD blocks of founders that were inherited by at least one individual in this data
#'
#' Transforms a tidy table describing the founder blocks inherited by each individual into a more summarized table describing the contiguous blocks of each founder inherited by at least one individual.
#' This is especially useful when the focal individuals are several generations removed from the founders so only a small fraction of founder chromosomes are inherited.
#'
#' @param inds The tidy table of recombination block data inherited in individuals from founders, such as the return value of [tidy_recomb_map_inds()].
#' More broadly, it is a table with a row for every IBD block, and at least these columns: "anc" is the label of the founder individual, "chr" is the chromosome number, "start" and "end" are the range of the block.
#' Although this data usually also has columns "ind" and "parent" identifying the individual and parental haplotype that inherited the block, this information is not used by this function and is removed if present.
#'
#' @return A table in a similar format as the input (columns "anc", "chr", "start", and "end") where each row is a contiguous IBD block whose each basepair is inherited by at least one individual.
#'
#' @examples
#' # manually construct a toy sample input,
#' # with individuals marked for clarity although they are not required
#' # note first two individuals
#' library(tibble)
#' inds <- tibble(
#'     ind = letters[1:3],
#'     parent = c('pat', 'mat', 'mat'),
#'     chr = c(1, 1, 2),
#'     start = c(   1,  800,   5),
#'     end   = c(1000, 2000, 100),
#'     anc   = c('f1', 'f1', 'f2')
#' )
#'
#' # the new table merges the first two rows,
#' # because they overlapped from the same ancestor,
#' # while the third row stays unchanged (after removing individual info)
#' founder_blocks <- recomb_founder_blocks_inherited( inds )
#'
#' @seealso
#' [tidy_recomb_map_inds()]
#'
#' @export
recomb_founder_blocks_inherited <- function( inds ) {
    # first erase individual info, if present, since we only care about founders
    inds <- dplyr::select( inds, -tidyselect::any_of( c( 'ind', 'parent' ) ) )
    # sort by founder (column `anc`), then by position
    inds <- dplyr::arrange( inds, 'anc', 'chr', 'start' )
    
    # the hardest part is merging overlapping regions if any
    # take advantage of data already being sorted
    # loop this way because the number of rows will get reduced if there are merges
    # always start with the second row
    i <- 2L
    while ( i < nrow( inds ) ) {
        # for mergings to happen, it has to be the same ancestor and chromosome, then check for overlap in the stand and end respectively
        if ( inds$anc[ i ] == inds$anc[ i - 1L ] && inds$chr[ i ] == inds$chr[ i - 1L ] && inds$start[ i ] <= inds$end[ i - 1L ] + 1L ) {
            # it's as easy as just extending the end of the previous row to be the end of the current row!
            inds$end[ i - 1L ] <- inds$end[ i ]
            # now delete the current row
            inds <- inds[ -i, ]
            # don't increment here because we want to process the new row `i` again
        } else {
            # else move on, have to increment manually
            i <- i + 1L
        }
    }
    # done, return the edited tibble
    return ( inds )
}

