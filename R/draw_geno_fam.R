# multiplexes draw_geno_child in an efficient manner, for large output matrices
# NOTE: this generalized version of original `draw_geno_children` also replaces old `sim_children_generations_genotypes`!
draw_geno_fam <- function(X, fam) {
    # validate inputs
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )

    # ensures that `X` and `fam` agree, maps parent indexes
    data <- match_fam_founders( fam, colnames( X ), 'X', 'column' )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    m_loci <- nrow( X )
    n_ind_fam <- nrow( fam )
    
    # initialize output matrix
    X_fam <- matrix(
        NA,
        nrow = m_loci,
        ncol = n_ind_fam
    )
    # copy names of individuals and loci
    colnames( X_fam ) <- fam$id
    rownames( X_fam ) <- rownames( X ) # may be NULL, meh
    # copy founders
    # reorder/subset columns of X if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    X_fam[ , fam$founder ] <- X[ , indexes ]
    
    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, as indexes of the current genotype matrix (precalculated)
        parents_j <- c( fam$pati[ i ], fam$mati[ i ] )
        # and extracts their genotypes
        X_parents_j <- X_fam[ , parents_j, drop = FALSE ]
        # draw child, add to output matrix
        X_fam[, i] <- draw_geno_child( X_parents_j )
    }

    # done, return genotype matrix of all FAM
    return( X_fam )
}

