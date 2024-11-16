#' Simulate genotypes with linkage disequilibrium (LD) given a population of haplotypes, using a Li-Stephens-like model of haplotype copying
#'
#' Each new genome has breaks drawn from the recombination model accelerated for the desired number of generations `G`, and haplotypes are drawn and copied randomly from the population.
#' This results in a population with individuals drawn independently and identically distributed but with LD within individuals, since the ancestral LD is preserved to some extent (attenuated by recombination after `G` generations).
#' If there is no LD in `haps`, then the output will not have LD either except when the number of columns of `haps` is very small (which resembles a bottleneck).
#' This model does not introduce mutations (unlike the original Li-Stephens).
#' Genotypes, when requested, are simply sums of independently drawn haplotype values.
#'
#' @inheritParams recomb_map_inds
#' @param haps Regular matrix or `BEDMatrix` object of haplotype values, one row per locus, one column per haplotype (half individual), or transposed if `loci_on_cols = TRUE` and for `BEDMatrix` objects.
#' If `geno = TRUE` (default), these values should be numeric (recommended are zeroes and ones, indicating absence or presence of reference allele), but if `geno = FALSE` code will work with any values, including strings, which are just copied to outputs in blocks.
#' @param bim The table of variants of `haps`, which is a data.frame/tibble with at least two columns: `chr` (must be numeric between 1 and the maximum chromosome in `map` below for map to work) and `pos` (base pair position, usually an integer).
#' @param G Number of generations since most recent common ancestor of population (to multiply standard recombination rate)
#' @param n_ind Number of individuals (if geno = TRUE) or haplotypes (half individuals, if geno = FALSE) desired in output
#' @param geno If `TRUE` (default) returns matrix of genotypes (values in 0,1,2 if `haps` is binary, otherwise double whatever the range of values in `haps` is), otherwise returns matrix of haplotypes (half individuals, same values of input `haps`)
#' @param loci_on_cols If `TRUE`, `haps` has loci on columns and individuals on rows; if `FALSE` (default), loci are on rows and individuals on columns.
#' If `haps` is a `BEDMatrix` object, `loci_on_cols` is ignored (set automatically to `TRUE` internally).
#' @param indexes_loci Vector of index ranges (start, end) of loci to simulate, to request to simulate only a subset of the loci in `haps`/`bim`.
#' Default `NULL` simulates all loci in `haps`/`bim`.
#' @param indexes_chr_ends Optional vector mapping each chromosome (index in vector) to the last index in the `bim` table; missing chromosomes must be mapped to `NA` values.
#' If `NULL` (default), this vector is calculated from the provided `bim` table.
#' However, if this function is called repeatedly on the same `bim` data, it can improve efficiency to calculate it outside just once (use [indexes_chr()] for this).
#'
#' @return A matrix with the same number of rows as `haps` and `n_ind` columns, with values copied from `haps` in (recombination) blocks if `geno = FALSE`, or sums of two such values drawn independently when `geno = TRUE`.
#'
#' @examples
#' # simulate a tiny population with few SNPs for example
#' library(tibble)
#' bim <- tibble(
#'   chr = c(2, 2, 3, 3, 22),
#'   pos = c(100, 121, 53, 154, 66) * 1e6
#' )
#' m_loci <- nrow( bim )
#' # Most often, haplotypes are binary data as simulated here.
#' # Here haplotypes will be totally unstructured, but to have LD in the output use real human data
#' # or data simulated to have LD
#' n_ind_haps <- 5
#' haps <- matrix(
#'   rbinom( m_loci * n_ind_haps, 1, 0.5 ),
#'   nrow = m_loci,
#'   ncol = n_ind_haps
#' )
#' # makes sense to have a lot of recombination at the population level
#' G <- 500
#' # ask for small output for example
#' n_ind <- 7
#' # use the recombination map for the same genome build as your data!
#' map <- recomb_map_hg38
#'
#' # simulate genotypes!  (Usually more convenient, but phase information is lost)
#' X <- pop_recomb( haps, bim, map, G, n_ind )
#'
#' # simulate haplotypes instead (preserves true phase)
#' H <- pop_recomb( haps, bim, map, G, n_ind, geno = FALSE )
#'
#' @export
pop_recomb <- function( haps, bim, map, G, n_ind, geno = TRUE, loci_on_cols = FALSE, indexes_loci = NULL, indexes_chr_ends = NULL ) {
    # validations
    if ( missing( haps ) )
        stop( '`haps` is required!' )
    if ( missing( bim ) )
        stop( '`bim` is missing!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( missing( G ) )
        stop( '`G` is required!' )
    if ( missing( n_ind ) )
        stop( '`n_ind` is required!' )
    # haps validations in detail
    if ( !is.matrix( haps ) )
        stop( '`haps` must be a matrix!' )
    # same as general matrix but transposed
    # this is always imposed for this particular format!
    if ( 'BEDMatrix' %in% class( haps ) ) {
        loci_on_cols <- TRUE
        # BEDMatrix is always numeric, no need to test, but BEDMatrix object returns false to is.numeric, so we have to hack it this way
    } else if ( geno && !is.numeric( haps ) )
        stop( '`haps` must be numeric if `geno = TRUE`!' )
    # bim validations in detail
    if ( !is.data.frame( bim ) )
        stop( '`bim` must be a data.frame (including tibble)!' )
    if ( !( 'chr' %in% names( bim ) ) )
        stop( '`bim` must have column `chr`!' )
    if ( !( 'pos' %in% names( bim ) ) )
        stop( '`bim` must have column `pos`!' )
    m_loci <- nrow( bim )
    if ( loci_on_cols ) {
        if ( m_loci != ncol( haps ) )
            stop( 'Number of rows of `bim` and columns of `haps` must be equal!' )
    } else {
        if ( m_loci != nrow( haps ) )
            stop( 'Number of rows of `bim` and `haps` must be equal!' )
    }
    # (in some parts chr could be non-numeric, but the genetic map does require them as indexes)
    if ( !is.numeric( bim$chr ) )
        stop( '`bim$chr` must be numeric!' )
    chrs <- unique( bim$chr )
    if ( !is.numeric( bim$pos ) )
        stop( '`bim$pos` must be numeric!' )
    # map validation in detail
    if ( !is.list( map ) )
        stop( '`map` must be a list!' )
    if ( any( chrs > length( map ) ) )
        stop( '`map` has ', length( map ), ' chromosomes but `chr` has these values that exceed that: ', chrs[ chrs > length( map ) ] )
    # G, n_ind validation in detail
    if ( !is.numeric( G ) )
        stop( '`G` must be numeric!' )
    if ( length( G ) != 1L )
        stop( '`G` must be scalar!' )
    if ( !is.numeric( n_ind ) )
        stop( '`n_ind` must be numeric!' )
    if ( length( n_ind ) != 1L )
        stop( '`n_ind` must be scalar!' )
    # indexes_loci validation
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
    }

    # calculate this if it wasn't provided, but providing it can be more efficient!
    if ( is.null( indexes_chr_ends ) )
        indexes_chr_ends <- indexes_chr( bim$chr )
    
    # final output to concatenate to
    X <- NULL

    # simulate each chromosome in turn
    for ( chr_i in chrs ) {
        # subset haplotype data
        # first retrieve range of chromosome
        index_range <- indexes_chr_range( indexes_chr_ends, chr_i )
        
        # subset ranges further if requested
        if ( !is.null( indexes_loci ) ) {
            # first check is that end of one is lower than start of the other, in these cases there's no overlap and we skip quietly
            if ( index_range[2] < indexes_loci[1] || indexes_loci[2] < index_range[1] )
                next
            # there is an overlap!  Reduce range if needed
            if ( index_range[1] < indexes_loci[1] )
                index_range[1] <- indexes_loci[1]
            if ( index_range[2] > indexes_loci[2] )
                index_range[2] <- indexes_loci[2]
        }
        pos_i <- bim$pos[ index_range[1] : index_range[2] ]
        map_i <- map[[ chr_i ]]
        m_loci_i <- length( pos_i )
        # NOTE: never subset `haps`, which is in general inefficient but extremely so when it is a BEDMatrix object.  It gets subset by internal `pop_recomb_chr` efficiently via `indexes_loci = indexes` option below

        # map_i validations in detail
        if ( !is.data.frame( map_i ) )
            stop( '`map[[', chr_i , ']]` must be a data.frame (including tibble)!' )
        if ( !( 'pos' %in% names( map_i ) ) )
            stop( '`map[[', chr_i , ']]` must have column `pos`!' )
        if ( !( 'posg' %in% names( map_i ) ) )
            stop( '`map[[', chr_i , ']]` must have column `posg`!' )
        
        # initialize output matrix for this chromosome
        X_i <- matrix( NA, nrow = m_loci_i, ncol = n_ind )

        # create data for every individual, which are IID from this distribution
        for ( j in 1L : n_ind ) {
            # draw one haplotype
            hap_new <- pop_recomb_chr( haps, pos_i, map_i, G, indexes_loci = index_range )

            if ( geno ) {
                # draw again and add to create genotypes
                hap_new <- hap_new + pop_recomb_chr( haps, pos_i, map_i, G, indexes_loci = index_range )
            }

            # store in matrix as desired
            X_i[ , j ] <- hap_new
        }

        # concatenate to output
        X <- rbind( X, X_i )
    }

    return( X )
}
