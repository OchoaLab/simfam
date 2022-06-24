#' Reduce haplotype data to genotype matrix
#'
#' This function accepts haplotype data, such as the output from [recomb_haplo_inds()], and reduces it to a genotype matrix.
#' The haplotype data is more detailed because it is phased, while phase is lost in the genotype representation.
#' Moreover, the haplotype data separates individuals and chromosomes into lists (the way it is simulated), but the output genotype matrix concatenates data from all chromosomes into a single matrix, as it appears in simpler simulations and real data.
#'
#' @param haplos A list of diploid individuals, each of which is a list with two haploid individuals named "pat" and "mat", each of which is a list of chromosomes.
#' Each chromosome can be a list, in which case the named element "x" must give the haplotype vector (ideally with values in zero and one counting reference alleles, including NA), otherwise the chromosome must be this vector (accommodating both output formats from [recomb_haplo_inds()] automatically).
#'
#' @return The genotype matrix, which is the sum of the haplotype values (with values in 0, 1, 2, and NA, counting reference alleles), with individuals along columns in same order as `haplos` list, and loci along rows in order of appearance concatenating chromosomes in numerical order.
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
#' # determine haplotypes of descendants given ancestral haplotypes
#' haplos <- recomb_haplo_inds( inds, haplo )
#'
#' # finally, run desired function!
#' # convert haplotypes structure to a plain genotype matrix
#' X <- recomb_geno_inds( haplos )
#'
#' @seealso
#' [recomb_fam()] for drawing recombination (ancestor) blocks, defined in terms of genetic distance.
#'
#' [recomb_map_inds()] for transforming genetic to basepair coordinates given a genetic map.
#'
#' [recomb_haplo_inds()] for determining haplotypes of descendants given ancestral haplotypes (creates input to this function).
#' 
#' @export
recomb_geno_inds <- function( haplos ) {
    if ( missing( haplos ) )
        stop( '`haplos` is required!' )
    if ( !is.list( haplos ) )
        stop( '`haplos` must be a list!' )
    
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

        # a growing vector of genotypes for this individual
        # have to process each chr separately, then concatenate that
        x_i <- NULL
        
        # navigate chromosomes
        for ( chr in 1L : n_chr ) {
            pat_i_chr <- pat_i[[ chr ]]
            mat_i_chr <- mat_i[[ chr ]]
            # two versions: if ancestry is missing, these are vectors, otherwise they are lists
            if ( is.list( pat_i_chr ) ) {
                pat_i_chr <- pat_i_chr$x
                if ( is.null( pat_i_chr ) )
                    stop( 'Individual ', i, ' element "pat" chr ', chr, ' was a list but is missing element "x"!' )
            }
            if ( is.list( mat_i_chr ) ) {
                mat_i_chr <- mat_i_chr$x
                if ( is.null( mat_i_chr ) )
                    stop( 'Individual ', i, ' element "mat" chr ', chr, ' was a list but is missing element "x"!' )
            }
            # now we should have vectors, validate that
            # because we will add them up, they must be numeric
            if ( !is.numeric( pat_i_chr ) )
                stop( 'Individual ', i, ' element "pat" chr ', chr, ' must be numeric!' )
            if ( !is.numeric( mat_i_chr ) )
                stop( 'Individual ', i, ' element "mat" chr ', chr, ' must be numeric!' )
            # now we can add them up!
            # this is the desired genotype data
            x_i_chr <- pat_i_chr + mat_i_chr
            # concatenate this data now!
            x_i <- c( x_i, x_i_chr )
        }
        
        # now that we're done adding this individual, need to add to bigger output matrix
        if ( i == 1L ) {
            # if this is the first individual, this is the first time we know the number of loci, so we initialize the matrix
            m_loci <- length( x_i )
            X <- matrix(
                0L,
                ncol = n_ind,
                nrow = m_loci
            )
        } else {
            if ( length( x_i ) != m_loci )
                stop( 'Individual ', i, ' number of loci (', length( x_i ), ') does not equal number of loci for previous individuals (', m_loci, ')!' )
        }
        # now add data
        X[ , i ] <- x_i
    }
    
    return( X )
}
