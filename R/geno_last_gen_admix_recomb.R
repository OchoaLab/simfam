#' Simulate an admixed family efficiently with founders with LD
#'
#' This function in essence combines [pop_recomb()] to simulate founders of known ancestries with LD (following a Li-Stephens-like model), draws recombination breaks of focal last-generation descendants from the specified pedigree using [recomb_last_gen()], and their genomes from the founders variants using [recomb_haplo_inds()].
#' However, since a limited portion of founder sequences is actually inherited, the simulation is made much more efficient by simulating only those subsequences that were inherited, which saves time, and utilizing sparse matrices, which saves memory too.
#' See below for a more detailed algorithm.
#'
#' This function wraps around several exported package functions to achieve its objectives, which are roughly grouped into the following 4 phases.
#' Phase 1 simulates recombination in the family without explicit sequences.
#' In particular, it initializes the founder haplotype structure (without variants yet) using [recomb_init_founders()], then simulates recombination breaks along the pedigree and identifies all the founder haplotype blocks in the focal individuals using [recomb_last_gen()], and maps recombination breaks in cM to basepairs using [recomb_map_inds()].
#' Phase 2 reorganizes this data to identify the unique founder blocks that were inherited, first by making the data tidy with [tidy_recomb_map_inds()], then applying [recomb_founder_blocks_inherited()].
#' Phase 3 initializes founder haplotypes using sparse matrices from the package `Matrix`, and draws inherited founder subsequences according to their known ancestries and using the Li-Stephens-like haplotype model of [pop_recomb()].
#' Phase 4 constructs the genotype matrices of focal individuals using the haplotypes of the founders drawn in phase 3 and the known origin of focal blocks from founders from phase 1, first constructing this data at the phased haplotype level with [recomb_haplo_inds()], reencoding as unphased genotypes using [recomb_geno_inds()], and constructing the corresponding local ancestry dosages using [recomb_admix_inds()].
#'
#' @param anc_haps A named list that maps the code used for each ancestry to its haplotype matrix.
#' Each of the haplotype matrices the argument `haps` passed to [pop_recomb()], namely is a regular matrix or `BEDMatrix` object of haplotype values, one row per locus, one column per haplotype (half individual), or transposed if `loci_on_cols = TRUE` and for `BEDMatrix` objects.
#' Here, these values must be numeric (recommended are zeroes and ones, indicating absence or presence of reference allele).
#' @inheritParams pop_recomb
#' @inheritParams recomb_last_gen
#' @param founders_anc a named vector that maps every founder haplotype (the names of this vector) to its ancestry code.
#' Ancestry codes must match the codes used in `anc_haps` above.
#' Founder haplotypes are the founder individual IDs from the pedigree (values in `ids[[1]]`) appearing twice, suffixed with "_pat" and "_mat", respectively (so the parents of the founders are unadmixed, though founders be first generation admixed this way).
#'
#' @return A named list with three elements:
#' - `X`: the genotype matrix of the focal individuals, as returned by [recomb_geno_inds()].
#' - `Ls`: a list, mapping each ancestry to its matrix of local ancestry dosages, as returned by [recomb_admix_inds()].
#' - `haplos`: a phased version of the haplotypes and local ancestries of the focal individuals, structured as nested lists, as returned by [recomb_haplo_inds()].
#'
#' @examples
#' library(tibble)
#' 
#' # simulate random haplotypes for example
#' # this toy data has 10 SNPs per chromosome, in fixed positions for simplicity
#' bim <- tibble( chr = rep( 1 : 22, each = 10 ), pos = rep( (1:10) * 1e6, 22 ) )
#' # and random haplotype data to go with this
#' n_ind_hap <- 10
#' m_loci <- nrow( bim )
#' # NOTE ancestry labels can be anything but must match `founders_anc` below
#' anc_haps <- list(
#'     'AFR' = matrix( rbinom( m_loci * n_ind_hap, 1L, 0.5 ), nrow = m_loci, ncol = n_ind_hap ),
#'     'EUR' = matrix( rbinom( m_loci * n_ind_hap, 1L, 0.2 ), nrow = m_loci, ncol = n_ind_hap )
#' )
#'
#' # now simulate a very small family with one individual, 2 parents, 4 implicit grandparents
#' data <- fam_ancestors( 2 )
#' fam <- data$fam
#' ids <- data$ids
#' # select ancestries for each of the 4 grandparents / founder haplotypes (unadmixed)
#' founders_anc <- c('AFR', 'AFR', 'AFR', 'EUR')
#' # set names of founders with _pat/mat, needed to match recombination structure
#' # order is odd but choices were random so that doesn't matter
#' names( founders_anc ) <- c(
#'     paste0( ids[[1]], '_pat' ),
#'     paste0( ids[[1]], '_mat' )
#' )
#'
#' # this performs the simulation!
#' data <- geno_last_gen_admix_recomb( anc_haps, bim, recomb_map_hg38, 10, fam, ids, founders_anc )
#' # this is the genotype matrix for the one admixed individual
#' data$X
#' # the corresponding local ancestry dosage matrices
#' # names match input labels
#' data$Ls$AFR
#' data$Ls$EUR
#' # if desired, a more complete but more complicated structure holding phased haplotypes
#' # and phased local ancestry information
#' data$haplos
#' 
#' @seealso
#' [recomb_init_founders()], [recomb_last_gen()], [recomb_map_inds()], [tidy_recomb_map_inds()], [recomb_founder_blocks_inherited()], [pop_recomb()], [recomb_haplo_inds()], [recomb_geno_inds()], [recomb_admix_inds()].
#'
#' @export
#' @importFrom rlang .data
geno_last_gen_admix_recomb <- function( anc_haps, bim, map, G, fam, ids, founders_anc, indexes_chr_ends = NULL, loci_on_cols = FALSE, missing_vals = c('', 0) ) {
    # other pop_recomb options left as defaults: indexes_loci = NULL
    
    # require everything that isn't default
    if ( missing( anc_haps ) )
        stop( '`anc_haps` is required!' )
    if ( missing( bim ) )
        stop( '`bim` is required!' )
    if ( missing( map ) )
        stop( '`map` is required!' )
    if ( missing( G ) )
        stop( '`G` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( missing( ids ) )
        stop( '`ids` is required!' )
    if ( missing( founders_anc ) )
        stop( '`founders_anc` is required!' )
    
    # validate new parameters, but the rest get validated by the other functions that use them
    if ( !is.list( anc_haps ) )
        stop( '`anc_haps` must be a list!' )
    anc_haps_names <- names( anc_haps )
    if ( is.null( anc_haps_names ) )
        stop( '`anc_haps` list must have named elements!' )
    if ( !is.vector( founders_anc ) )
        stop( '`founders_anc` must be a vector!' )
    if ( is.null( names( founders_anc ) ) )
        stop( '`founders_anc` vector must have named elements!' ) 
    if ( !all( founders_anc %in% anc_haps_names ) )
        stop( 'All ancestry labels in `founders_anc` must be present in the names of `anc_haps`!' )
    
    # initialize trivial founders chromosomes (unrecombined chromosomes with unique labels)
    # must use IDs of first generation from pedigree as input
    founders <- recomb_init_founders( ids[[ 1 ]], map )
    # draw recombination breaks along pedigree, with coordinates in genetic distance (centiMorgans)
    inds <- recomb_last_gen( founders, fam, ids, missing_vals = missing_vals )
    # map recombination break coordinates to base pairs
    inds <- recomb_map_inds( inds, map )

    # to be as efficient as possible, given that these haplotype-based genotype simulations are expensive, let's only simulate the haplotype segments of ancestors that actually had descendants
    # first tidy up the structure so it's easier to manipulate
    inds_tidy <- tidy_recomb_map_inds( inds )
    # then determine the segments of each ancestor that are present in at least one individual, with merged overlaps
    founder_blocks <- recomb_founder_blocks_inherited( inds_tidy )
    # now use this info to actually carry out sparse data simulated
    # sinced it's simulated by chromosome, sort by chr:start now
    founder_blocks <- dplyr::arrange( founder_blocks, .data$chr, .data$start )

    # precalculate this so these things aren't recalculated in every call to `pop_recomb` inside the loop
    # calculate this if it wasn't provided, but providing it can be more efficient!
    if ( is.null( indexes_chr_ends ) )
        indexes_chr_ends <- indexes_chr( bim$chr )

    # this awkward loop simulates exactly the ancestral haplotypes needed for last generation, and no more
    # a value we access repeatedly
    n_founders <- length( founders_anc )
    n_chr <- length( map )
    haplos_anc <- vector( 'list', n_chr )
    for ( chr_i in 1L : n_chr ) {
        # extract ranges for this chromosome
        founder_block_chr_i <- dplyr::filter( founder_blocks, .data$chr == chr_i )
        # map between indexes in this chromosome and indexes in global matrix
        index_range <- indexes_chr_range( indexes_chr_ends, chr_i )
        # skip if we don't have data for this chromosome (test data is like this)
        if ( all( is.na( index_range ) ) ) next
        # initialize sparse matrix of founder haplotypes
        X_chr <- Matrix::Matrix( nrow = index_range[2] - index_range[1] + 1L, ncol = n_founders, data = 0, sparse = TRUE )
        # copy names of founder haplotypes here, required!
        colnames( X_chr ) <- names( founders_anc )
        # navigate each of these rows
        for ( i in 1L : nrow( founder_block_chr_i ) ) {
            # find rows from original bim table to simulate
            indexes_sim_range <- indexes_chr_pos( bim$pos, index_range[1], index_range[2], founder_block_chr_i$start[ i ], founder_block_chr_i$end[ i ] )
            # skip if nothing was found
            if ( all( is.na( indexes_sim_range ) ) ) next
            # destination in chr-specific matrix is shifted down this way
            indexes_in_chr_range <- indexes_sim_range - index_range[1] + 1L
            # now get founder
            founder_i <- founder_block_chr_i$anc[ i ]
            # use appropriate haplotype matrix depending on ancestry of this founder
            X_in <- anc_haps[[ founders_anc[ founder_i ] ]]
            # simulate this segment only, store where we desire it
            X_chr[ indexes_in_chr_range[1] : indexes_in_chr_range[2], founder_i ] <- pop_recomb( X_in, bim, map, G, 1, geno = FALSE, indexes_loci = indexes_sim_range, loci_on_cols = loci_on_cols )
        }
        # add to structure, in a list
        haplos_anc[[ chr_i ]] <- list( X = X_chr, pos = bim$pos[ index_range[1] : index_range[2] ] )
    }
    
    # determine haplotypes of descendants given ancestral haplotypes
    # `ret_anc = TRUE` ensures ancestry labels are given at every position too! (needed for local ancestry calculation)
    haplos <- recomb_haplo_inds( inds, haplos_anc, ret_anc = TRUE )
    # and reduce haplotype data to genotypes, same standard matrix format as previous examples
    X <- recomb_geno_inds( haplos )
    # calculate local ancestry matrices (dosages of population ancestry)
    anc_map <- tibble::tibble(
                           anc = names( founders_anc ),
                           pop = founders_anc
                       )
    Ls <- recomb_admix_inds( haplos, anc_map )
    
    # return the desired data
    return( list( X = X, Ls = Ls, haplos = haplos ) )
}
