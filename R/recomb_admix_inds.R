#' Reduce haplotype ancestry data to population ancestry dosage matrices
#'
#' This function accepts haplotype data, such as the output from [recomb_haplo_inds()] with `ret_anc = TRUE` (required), and reduces it to a list of population ancestry dosage matrices.
#' In this context, "ancestors/ancestry" refer to haplotype blocks from specific ancestor individuals, whereas "population ancestry" groups these ancestors into populations (such as African, European, etc.).
#' Although the haplotype data separates individuals and chromosomes into lists (the way it is simulated), the output matrices concatenates data from all chromosomes into a single matrix, as it appears in simpler simulations and real data, and matching the format of [recomb_geno_inds()].
#'
#' @param haplos A list of diploid individuals, each of which is a list with two haploid individuals named `pat` and `mat`, each of which is a list of chromosomes, each of which must be a list with a named element `anc` must give the vector of ancestor names per position (the output format from [recomb_haplo_inds()] with `ret_anc = TRUE`).
#' @param anc_map A data.frame or tibble with two columns: `anc` lists every ancestor haplotype name present in `haplos`, and `pop` the population assignment of that haplotype.
#' @param pops Optional order of populations in output, by default sorted alphabetically from `anc_map$pop`.
#'
#' @return A named list of population ancestry dosage matrices, ordered as in `pops`, each of which counts populations in both alleles (in 0, 1, 2), with individuals along columns in same order as `haplos` list, and loci along rows in order of appearance concatenating chromosomes in numerical order.
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
#' map <- recomb_map_hg38[ 1L:2L ]
#' # initialize parents with this other function
#' founders <- recomb_init_founders( c('father', 'mother'), map )
#' # draw recombination breaks for child
#' inds <- recomb_fam( founders, fam )
#' # now add base pair coordinates to recombination breaks
#' inds <- recomb_map_inds( inds, map )
#'
#' # also need ancestral haplotypes
#' # these should be simulated carefully as needed, but for this example we make random data
#' haplo <- vector( 'list', length( map ) )
#' # names of ancestor haplotypes for this scenario
#' # (founders of fam$id but each with "_pat" and "_mat" suffixes)
#' anc_names <- c( 'father_pat', 'father_mat', 'mother_pat', 'mother_mat' )
#' n_ind <- length( anc_names )
#' # number of loci per chr, for toy test
#' m_loci <- 10L
#' for ( chr in 1L : length( map ) ) {
#'     # draw random positions
#'     pos_chr <- sample.int( max( map[[ chr ]]$pos ), m_loci )
#'     # draw haplotypes
#'     X_chr <- matrix(
#'         rbinom( m_loci * n_ind, 1L, 0.5 ),
#'         nrow = m_loci,
#'         ncol = n_ind
#'     )
#'     # required column names!
#'     colnames( X_chr ) <- anc_names
#'     # add to structure, in a list
#'     haplo[[ chr ]] <- list( X = X_chr, pos = pos_chr )
#' }
#' # determine haplotypes and per-position ancestries of descendants given ancestral haplotypes
#' haplos <- recomb_haplo_inds( inds, haplo, ret_anc = TRUE )
#'
#' # define individual to population ancestry map
#' # take four ancestral haplotypes from above, assign them population labels
#' anc_map <- tibble(
#'     anc = anc_names,
#'     pop = c('African', 'European', 'African', 'African')
#' )
#'
#' # finally, run desired function!
#' # convert haplotypes structure to list of population ancestry dosage matrices
#' Xs <- recomb_admix_inds( haplos, anc_map )
#'
#' @seealso
#' [recomb_fam()] for drawing recombination (ancestor) blocks, defined in terms of genetic distance.
#'
#' [recomb_map_inds()] for transforming genetic to basepair coordinates given a genetic map.
#'
#' [recomb_haplo_inds()] for determining haplotypes of descendants given ancestral haplotypes (creates input to this function).
#' 
#' @export
recomb_admix_inds <- function( haplos, anc_map, pops = sort( unique( anc_map$pop ) ) ) {
    if ( missing( haplos ) )
        stop( '`haplos` is required!' )
    if ( missing( anc_map ) )
        stop( '`anc_map` is required!' )
    if ( !is.list( haplos ) )
        stop( '`haplos` must be a list!' )
    if ( !is.data.frame( anc_map ) )
        stop( '`anc_map` must be a data.frame (including tibble)!' )
    if ( !('anc' %in% names( anc_map ) ) )
        stop( '`anc_map` must have column "anc"!' )
    if ( !('pop' %in% names( anc_map ) ) )
        stop( '`anc_map` must have column "pop"!' )
    if ( !missing( pops ) )
        if ( !all( anc_map$pop %in% pops ) )
            stop( 'Populations in `anc_map$pop` are missing in `pops`: ', toString( unique( anc_map$pop[ !(anc_map$pop %in% pops) ] ) ) )
    K <- length( pops )
    if ( K <= 1L )
        stop( 'The number of `pops` should be 2 or more!' )
    
    n_ind <- length( haplos )
    # process each individual
    for ( i in 1L : n_ind ) {
        # get each parental haplotype, validate extensively
        haplo_i <- haplos[[ i ]]
        if ( !is.list( haplo_i ) )
            stop( 'Individual ', i, ' is not a list!' )
        if ( !( 'pat' %in% names( haplo_i ) ) )
            stop( 'Individual ', i, ' must have element "pat"!' )
        if ( !( 'mat' %in% names( haplo_i ) ) )
            stop( 'Individual ', i, ' must have element "mat"!' )
        pat_i <- haplo_i$pat
        mat_i <- haplo_i$mat
        if ( !is.list( pat_i ) )
            stop( 'Individual ', i, ' element "pat" must be a list!' )
        if ( !is.list( mat_i ) )
            stop( 'Individual ', i, ' element "mat" must be a list!' )

        # record number of chromosomes only first time, require all individuals to have the same values!
        if ( i == 1L ) {
            n_chr <- length( pat_i )
        } else {
            if ( length( pat_i ) != n_chr )
                stop( 'Individual ', i, ' number of chromosomes in paternal haplotype (', length( pat_i ), ') does not match that of previous individuals (', n_chr, ')!' )
        }
        if ( length( mat_i ) != n_chr )
            stop( 'Individual ', i, ' number of chromosomes in paternal (', n_chr, ') and maternal (', length( mat_i ), ') haplotypes differs!' )

        # a growing vector of ancestry dosages for this individual
        # have to process each chr separately, then concatenate that
        # this initializes each ancestry with NULL
        x_i <- vector( 'list', K )
        names( x_i ) <- pops
        
        # navigate chromosomes
        for ( chr in 1L : n_chr ) {
            pat_i_chr <- pat_i[[ chr ]]
            mat_i_chr <- mat_i[[ chr ]]
            # ancestry version requires data to be lists with "anc" (not just haplotype vectors)
            if ( !is.list( pat_i_chr ) )
                stop( 'Individual ', i, ' element "pat" chr ', chr, ' must be a list!' )
            if ( !( "anc" %in% names( pat_i_chr ) ) )
                stop( 'Individual ', i, ' element "pat" chr ', chr, ' list is missing element "anc"!' )
            pat_i_chr <- pat_i_chr$anc
            # repeat for mat copy
            if ( !is.list( mat_i_chr ) )
                stop( 'Individual ', i, ' element "mat" chr ', chr, ' must be a list!' )
            if ( !( "anc" %in% names( mat_i_chr ) ) )
                stop( 'Individual ', i, ' element "mat" chr ', chr, ' list is missing element "anc"!' )
            mat_i_chr <- mat_i_chr$anc

            # translate each (individual founder) ancestry vector to populations
            indexes <- match( pat_i_chr, anc_map$anc )
            if ( anyNA( indexes ) )
                stop( 'Individual ', i, ' element "pat" chr ', chr, ' has ancestries missing in `anc_map`: ', toString( unique( pat_i_chr[ is.na( indexes ) ] ) ) )
            pat_i_chr <- anc_map$pop[ indexes ]
            # repeat for mat copy
            indexes <- match( mat_i_chr, anc_map$anc )
            if ( anyNA( indexes ) )
                stop( 'Individual ', i, ' element "mat" chr ', chr, ' has ancestries missing in `anc_map`: ', toString( unique( mat_i_chr[ is.na( indexes ) ] ) ) )
            mat_i_chr <- anc_map$pop[ indexes ]

            # now calculate dosages for every population
            for ( pop in pops ) {
                # this is the desired dosage data
                x_i_chr <- ( pat_i_chr == pop ) + ( mat_i_chr == pop )
                # concatenate this data now!
                x_i[[ pop ]] <- c( x_i[[ pop ]], x_i_chr )
            }
        }
        
        # now that we're done adding this individual, need to add to bigger output matrix
        if ( i == 1L ) {
            # if this is the first individual, this is the first time we know the number of loci, so we initialize the matrix
            # NOTE: same for all ancestries because they all came from the same x_i_chr
            m_loci <- length( x_i[[1]] )
            # initialize each population
            X <- vector( 'list', K )
            names( X ) <- pops
            for ( pop in pops )
                X[[ pop ]] <- matrix(
                    0L,
                    ncol = n_ind,
                    nrow = m_loci
                )
        } else {
            if ( length( x_i[[1]] ) != m_loci )
                stop( 'Individual ', i, ' number of loci (', length( x_i[[1]] ), ') does not equal number of loci for previous individuals (', m_loci, ')!' )
        }
        # now add data
        for ( pop in pops ) {
            X[[ pop ]][ , i ] <- x_i[[ pop ]]
        }
    }
    
    return( X )
}
