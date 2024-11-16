# simulate a single haplotype (one half individual) given a population of haplotypes, using a Li-Stephens-like model of haplotype copying
# this is organized differently than other recomb_* series because here we have no need to keep track of recombination locations, we just want the haplotypes
# @param haps Matrix of zeroes and ones, one row per location, one column per haplotype (half individual), transposed if BEDMatrix or if `loci_on_cols` is TRUE, for this chromosome only (or subset to be so with `indexes_loci`).  BEDMatrix also works!
# @param pos BP positions within given chromosome
# @param map Recombination map for given chromosome
# @param G Number of generations (to multiply standard recombination rate)
# @param loci_on_cols transpose haps
# @param indexes_loci subset hap loci, a range of indexes to keep.  Applied to haps only, `pos` has to be subset already!  (haps can be much bigger and can be read from disk, pos has no such advantage)
pop_recomb_chr <- function( haps, pos, map, G, loci_on_cols = FALSE, indexes_loci = NULL ) {
    # most of these are redundant with checks outside, but they are all low-computation so meh
    # validations
    if ( missing( haps ) )
        stop( '`haps` is required!' )
    if ( missing( pos ) )
        stop( '`pos` is required!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( missing( G ) )
        stop( '`G` is required!' )
    # detailed validations
    if ( !is.logical( loci_on_cols ) )
        stop( '`loci_on_cols` must be logical!' )
    # haps validations in detail
    if ( !is.matrix( haps ) )
        stop( '`haps` must be a matrix!' )
    # same as general matrix but transposed
    # this is always imposed for this particular format!
    if ( 'BEDMatrix' %in% class( haps ) ) 
        loci_on_cols <- TRUE
    if ( loci_on_cols ) {
        n_ind <- nrow( haps )
        m_loci <- ncol( haps )
    } else {
        m_loci <- nrow( haps )
        n_ind <- ncol( haps )
    }
    # by default output has same length as input, but indexes_loci changes that!
    m_loci_out <- m_loci
    # indexes_loci validations in detail
    if ( !is.null( indexes_loci ) ) {
        if ( length( indexes_loci ) != 2L )
            stop( '`indexes_loci` must have two elements only (it is a range)!' )
        if ( !is.numeric( indexes_loci ) )
            stop( '`indexes_loci` must be numeric!' )
        if ( indexes_loci[ 1L ] > indexes_loci[ 2L ] )
            stop( '`indexes_loci` must have the first value be equal or smaller than the second one!' )
        if ( indexes_loci[ 1L ] < 1L )
            stop( '`indexes_loci` must have first value >= 1!' )
        if ( indexes_loci[ 2L ] > m_loci )
            stop( '`indexes_loci` must have second value <= the number of loci in `haps`!' )
        m_loci_out <- indexes_loci[ 2L ] - indexes_loci[ 1L ] + 1L
    }
    # pos validations in detail
    if ( !is.numeric( pos ) )
        stop( '`pos` must be numeric!' )
    if ( length( pos ) != m_loci_out )
        stop( 'length of `pos` must equal either the number of rows of `haps` or the length implied by `indexes_loci` if not NULL!' )
    # G validation in detail
    if ( !is.numeric( G ) )
        stop( '`G` must be numeric!' )
    if ( length( G ) != 1L )
        stop( '`G` must be scalar!' )
    # map validations in detail
    if ( !is.data.frame( map ) )
        stop( '`map` must be a data.frame (including tibble)!' )
    if ( !( 'pos' %in% names( map ) ) )
        stop( '`map` must have column `pos`!' )
    if ( !( 'posg' %in% names( map ) ) )
        stop( '`map` must have column `posg`!' )

    # infer chromosome length in genetic positions (cM) from map
    chr_leng <- max( map$posg )
    
    # draw random recombination breaks for a single chromosome
    # here we speed up recombination rate for multiple generations
    breaks_posg <- recomb_breaks( chr_leng, mean = 100 / G )
    # add trivial length to end (it is missing otherwise) so its processed in the data
    breaks_posg <- c( breaks_posg, chr_leng )
    # convert recombination breaks from genetic to BP positions
    breaks_pos <- recomb_map_posg( breaks_posg, map )
    
    # start generating new haplotype
    # (closest analog is `recomb_haplo_chr`, from which much of the following code was copied)

    # get started copying data from ancestors (`X`, `pos`) to new individual (`chr`)!
    # inputs are haploid so output will be too
    # this is the new individual's haploid variant vector (technically not "genotype", but meh)
    # at this level values don't need to be numeric, they'll just be copied whatever they are, though more specifically if they are integers the output will be too, but will allow numeric in general and also other options
    # BEDMatrix: storage.mode returns "S4", not what we want, so here we hack that case
    mode <- if ( 'BEDMatrix' %in% class( haps ) ) 'integer' else storage.mode( haps )
    hap_new <- vector( mode, m_loci_out )

    # in keeping with end-only data, this tells us where we left off in pos vector
    i_pos_start <- 1L
    # navigate chr rows
    for ( i_breaks in 1L : length( breaks_pos ) ) {
        # find greatest `pos <= breaks_pos[ i_breaks ]`
        # NOTE: this is end-inclusive!
        i_pos_end <- pos <= breaks_pos[ i_breaks ]
        # check for cases where there are no such cases
        # note i_pos_start is unchanged in those cases
        if ( !any( i_pos_end ) )
            next
        # turn logicals to indexes
        i_pos_end <- which( i_pos_end )
        # and because which returns them in order, we want the last one if there was more than one
        if ( length( i_pos_end ) > 1L )
            i_pos_end <- i_pos_end[ length( i_pos_end ) ]
        # another required condition is that end is later or equal than start
        if ( i_pos_start > i_pos_end )
            next
        # if all of these passed, we've found a segment (of length 1 or more) to copy from ancestor to descendant!
        
        # now copy data
        # these are destination indexes
        indexes_pos <- i_pos_start : i_pos_end
        # source indexes (in `haps` matrix) are the same unless we're subsetting, in which case map those indexes now too
        indexes_pos_haps <- indexes_pos
        if ( !is.null( indexes_loci ) ) {
            i_pos_haps_start <- indexes_loci[1L] - 1L + i_pos_start
            i_pos_haps_end <- indexes_loci[1L] - 1L + i_pos_end 
            indexes_pos_haps <- i_pos_haps_start : i_pos_haps_end
            # old way of doing this
            indexes_pos_haps2 <- ( indexes_loci[1L] : indexes_loci[2L] )[ indexes_pos ]
            # confirm that they match
            stopifnot( all( indexes_pos_haps == indexes_pos_haps2 ) )
        }
        # select random haplotype to copy (ok to rarely select same one as before in a row)
        index_hap <- sample.int( n_ind, 1L )
        hap_new[ indexes_pos ] <- if ( loci_on_cols ) {
                                      haps[ index_hap, indexes_pos_haps ]
                                  } else {
                                      haps[ indexes_pos_haps, index_hap ]
                                  }
        
        # after a successful copy, can now increment i_pos_start for next round
        i_pos_start <- i_pos_end + 1L
    }

    return( hap_new )
}
