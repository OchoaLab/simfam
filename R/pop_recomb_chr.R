# simulate a single haplotype (one half individual) given a population of haplotypes, using a Li-Stephens-like model of haplotype copying
# this is organized differently than other recomb_* series because here we have no need to keep track of recombination locations, we just want the haplotypes
# @param haps Matrix of zeroes and ones, one row per location, one column per haplotype (half individual), for this chromosome only
# @param pos BP positions within given chromosome
# @param map Recombination map for given chromosome
# @param G Number of generations (to multiply standard recombination rate)
pop_recomb_chr <- function( haps, pos, map, G ) {
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
    # haps validations in detail
    if ( !is.matrix( haps ) )
        stop( '`haps` must be a matrix!' )
    m_loci <- nrow( haps )
    n_ind <- ncol( haps )
    # pos validations in detail
    if ( !is.numeric( pos ) )
        stop( '`pos` must be numeric!' )
    if ( length( pos ) != m_loci )
        stop( 'length of `pos` must equal number of rows of `haps`!' )
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
    hap_new <- vector( storage.mode( haps ), m_loci )

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
        indexes_pos <- i_pos_start : i_pos_end
        # select random haplotype to copy (ok to rarely select same one twice)
        index_hap <- sample.int( n_ind, 1L )
        hap_new[ indexes_pos ] <- haps[ indexes_pos, index_hap ]

        # after a successful copy, can now increment i_pos_start for next round
        i_pos_start <- i_pos_end + 1L
    }

    return( hap_new )
}
