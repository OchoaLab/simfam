# replaces old `update_admix_proportions_children` and `sim_children_generations_admix_proportions`
# NOTE: doesn't check that `admix` is valid (i.e that rows are non-negative and sum to one), just averages data
admix_fam <- function(admix, fam) {
    if ( missing( admix ) )
        stop( '`admix` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( admix ) )
        stop( '`admix` must be a matrix!' )
    
    # ensures that `admix` and `fam` agree, maps parent indexes
    # also loads of `fam` validations
    data <- match_fam_founders( fam, rownames( admix ), 'admix', 'row' )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    n_ind_fam <- nrow( fam )
    k_subpops <- ncol( admix )
    
    # initialize output matrix
    admix_fam <- matrix(
        NA,
        nrow = n_ind_fam,
        ncol = k_subpops
    )
    # copy names of individuals and ancestries
    rownames( admix_fam ) <- fam$id
    colnames( admix_fam ) <- colnames( admix )
    # copy founders
    # reorder/subset columns of X if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    admix_fam[ fam$founder, ] <- admix[ indexes, ]

    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, as indexes of the current admix matrix (precalculated)
        parents_j <- c( fam$pati[ i ], fam$mati[ i ] )
        # and extracts their admixture proportions, which are simply averaged
        # store directly at child location
        admix_fam[ i, ] <- colMeans( admix_fam[ parents_j, ] )
    }
    
    # return entire matrix!
    return( admix_fam )
}
