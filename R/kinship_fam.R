# output agrees with [kinship2::kinship()] when founders are unrelated (`kinship = diag(n)/2` for the right number of founders and column/row names matching the founders in `fam`).
# This function is more general since it allows founders to be related in arbitrary ways
kinship_fam <- function(kinship, fam) {
    if ( missing( kinship ) )
        stop( '`kinship` is required!' )
    if ( missing( fam ) )
        stop( '`fam` is required!' )
    if ( !is.matrix( kinship ) )
        stop( '`kinship` must be a matrix!' )
    if ( !isSymmetric( kinship ) )
        stop( '`kinship` must be a symmetric matrix!' )
    
    # ensures that `kinship` and `fam` agree, maps parent indexes
    # also loads of `fam` validations
    data <- match_fam_founders( fam, rownames( kinship ), 'kinship', 'row' )
    fam <- data$fam # has new columns: founder, pati, mati
    indexes <- data$indexes # to reorder founders if needed
    
    # dimensions of output data
    n_ind_fam <- nrow( fam )
    
    # initialize output matrix
    kinship_fam <- matrix(
        NA,
        nrow = n_ind_fam,
        ncol = n_ind_fam
    )
    # copy names of individuals
    rownames( kinship_fam ) <- fam$id
    colnames( kinship_fam ) <- fam$id
    # copy founders
    # reorder/subset columns of X if needed!
    # place them where those founders are in `fam` (`fam$founder` is logical)
    kinship_fam[ fam$founder, fam$founder ] <- kinship[ indexes, indexes ]

    # navigate individuals
    for ( i in 1 : n_ind_fam ) {
        # skip founders
        if ( fam$founder[i] ) next
        # get parents, as indexes of the current kinship matrix (precalculated)
        p1 <- fam$pati[ i ]
        p2 <- fam$mati[ i ]
        # estimate kinship between child and everybody else
        # simple average works for everybody, including parents themselves
        kinship_i <- ( kinship_fam[ p1, ] + kinship_fam[ p2, ] ) / 2
        # copy back to big matrix
        kinship_fam[ i, ] <- kinship_i
        kinship_fam[ , i ] <- kinship_i
        # inbreeding is the kinship of the parents
        # encode as proper self kinship
        kinship_fam[ i, i ] <- ( 1 + kinship_fam[ p1, p2 ] ) / 2
    }
    
    # return entire matrix!
    return( kinship_fam )
}
