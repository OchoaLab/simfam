library(tibble)

test_that( "draw_couples_nearest works", {
    # create kinship for unrelated founders first
    n <- 10
    kinship_local <- diag( n ) / 2
    # even sex distribution to ensure prefect pairing
    sex <- rep.int( c(1, 2), n / 2 )

    # cause errors on purpose
    # first two arguments are mandatory
    expect_error( draw_couples_nearest( ) )
    expect_error( draw_couples_nearest( kinship_local = kinship_local ) )
    expect_error( draw_couples_nearest( sex = sex ) )
    expect_error( draw_couples_nearest( 1:n, sex = sex ) ) # must be matrix
    expect_error( draw_couples_nearest( cbind( 1:n ), sex = sex ) ) # must be square
    expect_error( draw_couples_nearest( kinship_local, sex[-1] ) ) # dimensions must agree

    # a successful run
    expect_silent(
        parents <- draw_couples_nearest( kinship_local, sex )
    )
    expect_true( is.matrix( parents ) )
    expect_equal( nrow( parents ), 2 )
    expect_equal( ncol( parents ), n / 2 ) # equality only because everybody is unrelated in this case
    expect_true( !anyNA( parents ) )
    # parents are a subset of all individuals available (as indexes)
    expect_true( all( parents %in% 1 : n ) )
    # and no parent appears twice
    expect_equal( length( parents ), length( unique( parents ) ) )
    # verify parent sexes
    expect_true( all( sex[ parents[ 1, ] ] == 1 ) )
    expect_true( all( sex[ parents[ 2, ] ] == 2 ) )
    # no need to verify kinship here because it was all trivial

    # create a smaller case where a unique solution exists
    # start with two nuclear families (can't be paired within because of thresholds)
    # the numbers are hacked so that across only one pairing works for each sibling (not realistic numbers but good for test)
    # pairings must have local kinship strictly < 1 / 4^3 = 0.015625
    n <- 4
    sex <- c(1, 2, 2, 1)
    kinship_local <- matrix(
        c(
            0.50, 0.25, 0.00, 0.06,
            0.25, 0.50, 0.06, 0.00,
            0.00, 0.06, 0.50, 0.25,
            0.06, 0.00, 0.25, 0.50
        ),
        nrow = n,
        ncol = n
    )
    expect_silent(
        parents <- draw_couples_nearest( kinship_local, sex )
    )
    expect_true( is.matrix( parents ) )
    expect_equal( nrow( parents ), 2 )
    expect_equal( ncol( parents ), n / 2 ) # exact solution also ensures equality here
    expect_true( !anyNA( parents ) )
    # parents are a subset of all individuals available (as indexes)
    expect_true( all( parents %in% 1 : n ) )
    # and no parent appears twice
    expect_equal( length( parents ), length( unique( parents ) ) )
    # unfortunately, though the pairings are unique, there are two column permutations in which they can appear in output
    # NOTE: sex means fathers (first row) can only be individuals 1 and 4
    parents_exp = cbind( c(1,3), c(4,2) )
    expect_true( parents == parents_exp || parents == parents_exp[, 2:1 ] )
    # verify parent sexes
    expect_true( all( sex[ parents[ 1, ] ] == 1 ) )
    expect_true( all( sex[ parents[ 2, ] ] == 2 ) )
    # the restricted parent pairings confirms kinship filter was correct (no need for additional kinship validation)

    # now a more challenging random example
    cutoff <- 1 / 4^3 # have out here for tests (same as default)
    # also make it odd, should work now!
    n <- 11
    # a random positive definite matrix with all values between 0 and 1
    # divided by extra 1/n to make numbers tiny, and as it is there are frequently between 0-2 pairs, great test for that edge case
    kinship_local <- crossprod( matrix( runif( n*n ), nrow = n, ncol = n ) ) / n^2
    expect_true( min( kinship_local ) >= 0 )
    expect_true( max( kinship_local ) <= 1 )
    sex <- sample( c(1, 2), n, replace = TRUE ) # completely random sex too
    # start tests
    expect_silent(
        parents <- draw_couples_nearest( kinship_local, sex, cutoff = cutoff )
    )
    expect_true( is.matrix( parents ) )
    expect_equal( nrow( parents ), 2 )
    expect_true( ncol( parents ) < n / 2 ) # must be fewer here because n is odd!
    #message( 'pairs: ', ncol( parents ) )
    expect_true( !anyNA( parents ) )
    # parents are a subset of all individuals available (as indexes)
    expect_true( all( parents %in% 1 : n ) )
    # and no parent appears twice
    expect_equal( length( parents ), length( unique( parents ) ) )
    # verify parent sexes
    expect_true( all( sex[ parents[ 1, ] ] == 1 ) )
    expect_true( all( sex[ parents[ 2, ] ] == 2 ) )
    # verify kinship (only non-trivial test here)
    # (sometimes there's no parents at all)
    n_fam <- ncol( parents )
    if ( n_fam > 0 ) {
        kinship_parents <- vector( 'numeric', n_fam )
        for ( j in 1 : n_fam ) {
            kinship_parents[ j ] <- kinship_local[ parents[ 1, j ], parents[ 2, j ] ]
        }
        expect_true( all( kinship_parents < cutoff ) )
    }

})

test_that( "draw_num_children_per_fam works", {
    pop_size <- 30
    n_fam <- 10
    children_min <- 1
    
    expect_silent( 
        children_per_fam <- draw_num_children_per_fam( n_fam, pop_size, children_min )
    )
    expect_equal( length( children_per_fam ), n_fam )
    expect_equal( sum( children_per_fam ), pop_size )
    expect_true( min( children_per_fam ) >= children_min )
})

test_that( "match_fam_founders", {
    # smallest toy example
    fam <- tibble(
        id = c('a', 'b', 'c'),
        pat = c(0, 0, 'a'),
        mat = c(0, 0, 'b')
    )
    names_founders <- c('a', 'b')
    # start with case that works
    expect_silent( 
        data <- match_fam_founders( fam, names_founders, 'X', 'column' )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'indexes') )
    # check that new columns are as they should be
    fam2 <- data$fam
    expect_equal( fam2$founder, c(TRUE, TRUE, FALSE) )
    expect_equal( fam2$pati, c(NA, NA, 1) )
    expect_equal( fam2$mati, c(NA, NA, 2) )
    expect_equal( data$indexes, 1:2 )
    # ultimately we want this to be true
    expect_equal( fam2$id[ fam2$founder ], names_founders[ data$indexes ] )

    # works with founders in reversed order!
    expect_silent( 
        data <- match_fam_founders( fam, rev(names_founders), 'X', 'column' )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'indexes') )
    # check that new columns are as they should be
    fam2 <- data$fam
    expect_equal( fam2$founder, c(TRUE, TRUE, FALSE) )
    expect_equal( fam2$pati, c(NA, NA, 1) )
    expect_equal( fam2$mati, c(NA, NA, 2) )
    expect_equal( data$indexes, 2:1 )
    # ultimately we want this to be true
    expect_equal( fam2$id[ fam2$founder ], rev(names_founders)[ data$indexes ] )
    
    # works with founders with extra IDs
    names_founders <- c('a', 'd', 'b')
    expect_silent( 
        data <- match_fam_founders( fam, names_founders, 'X', 'column' )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'indexes') )
    # check that new columns are as they should be
    fam2 <- data$fam
    expect_equal( fam2$founder, c(TRUE, TRUE, FALSE) )
    expect_equal( fam2$pati, c(NA, NA, 1) )
    expect_equal( fam2$mati, c(NA, NA, 2) )
    expect_equal( data$indexes, c(1, 3) )
    # ultimately we want this to be true
    expect_equal( fam2$id[ fam2$founder ], names_founders[ data$indexes ] )
    
    # cause errors on purpose
    # missing one founder
    expect_error( match_fam_founders( fam, c('a'), 'X', 'column' ) )
    expect_error( match_fam_founders( fam, c('b'), 'X', 'column' ) )

    # larger toy example
    # founders intermingled with non-founders
    fam <- tibble(
        id  = c('a', 'b', 'c', 'd', 'e', 'f'),
        pat = c(  0,   0, 'a',   0,   0, 'e'),
        mat = c(  0,   0, 'b',   0,   0, 'd')
    )
    names_founders <- c('d', 'e', 'a', 'b')
    # expect to work
    expect_silent( 
        data <- match_fam_founders( fam, names_founders, 'X', 'column' )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'indexes') )
    # check that new columns are as they should be
    fam2 <- data$fam
    expect_equal( fam2$founder, c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) )
    expect_equal( fam2$pati, c(NA, NA, 1, NA, NA, 5) )
    expect_equal( fam2$mati, c(NA, NA, 2, NA, NA, 4) )
    expect_equal( data$indexes, c(3, 4, 1, 2) )
    # ultimately we want this to be true
    expect_equal( fam2$id[ fam2$founder ], names_founders[ data$indexes ] )
})


# tests for recurrent data
validate_kinship_proper <- function( kinship, ids ) {
    expect_true( is.matrix( kinship ) )
    expect_true( isSymmetric( kinship ) )
    expect_equal( nrow( kinship ), length( ids ) )
    expect_true( !anyNA( kinship ) )
    expect_true( min( kinship ) >= 0 )
    expect_true( max( kinship ) <= 1 )
    expect_equal( rownames( kinship ), ids )
}

# kinship calculated from FAM table using external package `kinship2`, for validation
kin2 <- function( fam ) {
    # need to massage data so kinship2 likes it
    # need to transform missing parents into empty strings (zeroes only work when other IDs are numeric)
    fam$pat[ fam$pat == 0 ] <- ''
    fam$mat[ fam$mat == 0 ] <- ''
    # create the pedigree object
    pedigree <- kinship2::pedigree(
                              id = fam$id,
                              dadid = fam$pat,
                              momid = fam$mat,
                              sex = fam$sex,
                              affected = fam$pheno
                          )
    
    # calculate full kinship matrix using this package, which assumes (as we do) that founders are unrelated
    return( kinship2::kinship( pedigree ) )
}

test_that( "kinship_fam works", {
    # here any kinship matrix will work
    n <- 9
    n2 <- floor( n / 2 )
    # a different number for more fun
    n_kids <- 11
    # a random kinship matrix for founders
    # a random positive definite matrix with all values between 0 and 1
    # divided by extra 1/n to make numbers tiny, and as it is there are frequently between 0-2 pairs, great test for that edge case
    kinship <- crossprod( matrix( runif( n*n ), nrow = n, ncol = n ) ) / n^2
    # new code needs names, use letters to be extra non-trivial
    rownames( kinship ) <- letters[ 1 : n ]
    colnames( kinship ) <- letters[ 1 : n ]
    # individuals just paired consecutively
    # (parents matrix in old format)
    parents <- matrix( letters[ 1 : ( 2 * n2 ) ], nrow = 2, ncol = n2 )
    # draw family sizes as usual
    children_per_fam <- draw_num_children_per_fam( n2, n_kids )
    # minimal plink FAM table
    # contains founders and their children! (complete pedigree)
    fam <- tibble(
        # kid names should be letters after that, non-overlapping with parents
        id = c( letters[ 1 : n ], letters[ n + (1 : n_kids) ] ),
        # founders have missing parents
        # for children, expand parents * children_per_fam to make FAM pat/mat vectors
        pat = c( rep.int( 0, n ), rep.int( parents[ 1, ], children_per_fam ) ),
        mat = c( rep.int( 0, n ), rep.int( parents[ 2, ], children_per_fam ) )
    )
    
    # make kinship for all
    expect_silent(
        kinship_all <- kinship_fam( kinship, fam )
    )
    validate_kinship_proper( kinship_all, fam$id )
    # and check that kinship of founders appears in output as it should (here trivial because orders agree)
    expect_equal( kinship_all[ 1:n, 1:n ], kinship )

    # a trivial case for unrelated founders results in agreement with kinship2
    if ( suppressMessages(suppressWarnings(require(kinship2))) ) {
        # recalculate our way first, with unrelated founders
        kinship <- diag( n ) / 2
        # new code needs names, use letters to be extra non-trivial
        rownames( kinship ) <- letters[ 1 : n ]
        colnames( kinship ) <- letters[ 1 : n ]
        expect_silent(
            kinship_all <- kinship_fam( kinship, fam )
        )
        
        # calculate the kinship matrix with their code
        # first fix `fam` (kinship2 triggers checks that are irrelevant here, but meh)
        fam$sex <- rep_len( c(1L, 2L), nrow( fam ) ) # this agrees with our parent setup
        fam$pheno <- 0
        kinship_kinship2 <- kin2( fam )
        
        # these should be identical!
        expect_equal( kinship_all, kinship_kinship2 )
    }

    # a smaller toy example with unrelated founders, but with interlaced founders and non-founders, and different kinship order, to check that case explicitly (and also not dependent on kinship2 being available for testing)
    fam <- tibble(
        id  = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
        pat = c(  0,   0, 'a',   0,   0, 'e', 'f'),
        mat = c(  0,   0, 'b',   0,   0, 'd', 'c')
    )
    n <- 4
    # since all are unrelated, order is unimportant for numeric values, though column names will change that
    names_founders <- c('d', 'e', 'a', 'b')
    kinship <- diag( n ) / 2
    colnames( kinship ) <- names_founders
    rownames( kinship ) <- names_founders
    # expected output
    # note order is that of `fam`
    kinship_all_exp <- matrix(
        c(
            1/2,   0, 1/4,   0,   0,   0, 1/8,
              0, 1/2, 1/4,   0,   0,   0, 1/8,
            1/4, 1/4, 1/2,   0,   0,   0, 1/4,
              0,   0,   0, 1/2,   0, 1/4, 1/8,
              0,   0,   0,   0, 1/2, 1/4, 1/8,
              0,   0,   0, 1/4, 1/4, 1/2, 1/4,
            1/8, 1/8, 1/4, 1/8, 1/8, 1/4, 1/2
        ),
        nrow = nrow(fam),
        ncol = nrow(fam)
    )
    colnames( kinship_all_exp ) <- fam$id
    rownames( kinship_all_exp ) <- fam$id
    expect_silent(
        kinship_all <- kinship_fam( kinship, fam )
    )
    expect_equal( kinship_all, kinship_all_exp )
})

validate_fam <- function( fam, n, G ) {
    # fix n if needed
    if ( length( n ) == 1 )
        n <- rep.int( n, G )
    expect_true( is.data.frame( fam ) )
    expect_equal( names( fam ), c('fam', 'id', 'pat', 'mat', 'sex', 'pheno') )
    expect_equal( nrow( fam ), sum( n ) )
    # check values broadly
    expect_true( all( fam$fam == 'fam1' ) ) # dummy default
    expect_true( all( fam$pheno == 0 ) ) # dummy default
    expect_true( all( fam$sex %in% c(1L, 2L) ) )
    # check that all IDs we expect are indeed in FAM file, in the expected order!
    # need a loop to build vector
    ids_exp <- NULL
    for ( g in 1 : G ) {
        ids_exp <- c( ids_exp, paste0( 'g', g, '-', 1 : n[g] ) )
    }
    expect_equal( fam$id, ids_exp )
    # check that everybody's parents are the sex the should be
    # NOTE: list of parents exclude founders (who have no parents)
    expect_true( all( fam$sex[ fam$id %in% fam$pat[ -(1:n[1]) ] ] == 1L ) ) # all fathers are male
    expect_true( all( fam$sex[ fam$id %in% fam$mat[ -(1:n[1]) ] ] == 2L ) ) # all mothers are female
}

test_that( "sim_pedigree works", {
    G <- 3
    n <- 16

    # cause errors on purpose
    expect_error( sim_pedigree() )
    expect_error( sim_pedigree( G = G ) )
    expect_error( sim_pedigree( n = n ) )

    # a minimal, successful run
    # NOTE: had to set sex of founders to be exactly half male/female because otherwise they are set randomly, and there is a good chance (for small `n`) that they are all the same sex, in which case there's no solution!
    expect_silent(
        data <- sim_pedigree( G, n, sex = rep_len( c(1L, 2L), n ), verbose = FALSE )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'kinship_local') )
    # check fam
    validate_fam( data$fam, n, G )
    # check kinship_local
    names_kinship_exp <- paste0( 'g', G, '-', 1 : n )
    validate_kinship_proper( data$kinship_local, names_kinship_exp )

    # it'd be nice to compare kinship calculated by `kinship2` package!
    if ( suppressMessages(suppressWarnings(require(kinship2))) ) {
        # calculate the kinship matrix directly from FAM
        kinship_kinship2 <- kin2( data$fam )
        
        # subset because we only have last generation (though that suffices, since it couldn't be correct if the previous steps weren't also correct due to the recursivity of calculations)
        indexes <- (G-1) * n + 1:n # these are the indexes we want
        kinship_kinship2 <- kinship_kinship2[ indexes, indexes ]
        # finally, compare!
        expect_equal( data$kinship_local, kinship_kinship2 )
    }

    # do a version with variable `n` per generation and full output
    # better to have population grow since there's a minimum number of children
    G <- 3
    n <- c(16, 19, 21)

    # another type of error if G and n are mismatched in dimensions
    expect_error( sim_pedigree( G, n[1:2], verbose = FALSE ) )

    # and now the proper run
    # again set sex to ensure at least one pair in first generation
    expect_silent(
        data <- sim_pedigree( G, n, sex = rep_len( c(1L, 2L), n[1] ), verbose = FALSE, full = TRUE )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('fam', 'kinship_local') )
    # check fam
    validate_fam( data$fam, n, G )
    # check kinship_local
    # this time it's a list because `full = TRUE`
    expect_true( is.list( data$kinship_local ) )
    expect_equal( length( data$kinship_local ), G )
    for ( g in 1 : G ) {
        names_kinship_exp <- paste0( 'g', g, '-', 1 : n[g] )
        validate_kinship_proper( data$kinship_local[[g]], names_kinship_exp )
    }
    # compare to `kinship2` package again!
    if ( suppressMessages(suppressWarnings(require(kinship2))) ) {
        # calculate the kinship matrix directly from FAM
        kinship_kinship2 <- kin2( data$fam )

        # here we subset each generation and compare
        for ( g in 1 : G ) {
            # subset the current generation
            # these are the indexes we want
            indexes <- sum( n[ seq_len(g-1) ] ) + 1 : n[g]
            kinship_kinship2_g <- kinship_kinship2[ indexes, indexes ]
            # finally, compare!
            expect_equal( data$kinship_local[[ g ]], kinship_kinship2_g )
        }
    }
})

test_that( "draw_geno_child works", {
    # cause errors on purpose
    expect_error( draw_geno_child() )
    expect_error( draw_geno_child( 1:10 ) ) # must be matrix
    expect_error( draw_geno_child( matrix(0, 3, 3) ) ) # must have 2 columns
    
    # though normally output is random, some cases are deterministic
    # here parents are opposite homozygous
    m <- 3
    X <- cbind( rep.int( 0L, m ), rep.int( 2L, m ) )
    expect_silent(
        x <- draw_geno_child( X )
    )
    expect_equal( x, rep.int( 1L, m ) )
    
    # same homozygous
    x_exp <- c(0L, 2L, 0L)
    X <- cbind( x_exp, x_exp )
    expect_silent(
        x <- draw_geno_child( X )
    )
    expect_equal( x, x_exp )

    # heterozygous, result is random but certainly have some minimal expectations for genotypes
    X <- matrix( 1L, nrow = m, ncol = 2 )
    expect_silent(
        x <- draw_geno_child( X )
    )
    expect_equal( length( x ), m )
    expect_true( all( x %in% 0L:2L ) )

    # and completely random data
    X <- matrix( rbinom( m * 2, 2, 0.5 ), nrow = m, ncol = 2 )
    expect_silent(
        x <- draw_geno_child( X )
    )
    expect_equal( length( x ), m )
    expect_true( all( x %in% 0L:2L ) )
})

test_that( "draw_geno_fam works", {
    # construct genotypes for parents
    G <- 2
    n <- c( 4, 5 )
    m <- 10
    # random genotype data
    X <- matrix( rbinom( m * n[1], 2, 0.5 ), nrow = m, ncol = n[1] )
    
    # and FAM table
    # NOTE: had to set sex of founders to be exactly half male/female because otherwise they are set randomly, and there is a good chance (for small `n`) that they are all the same sex, in which case there's no solution!
    expect_silent(
        data <- sim_pedigree( G, n, sex = rep_len( c(1L, 2L), n[1] ), verbose = FALSE )
    )
    # extract table including both generations
    fam <- data$fam
    
    # cause errors on purpose
    # stuff missing
    expect_error( draw_geno_fam() )
    expect_error( draw_geno_fam( X = X ) )
    expect_error( draw_geno_fam( fam = fam ) )
    # X is missing colnames!
    expect_error( draw_geno_fam( X, fam ) )

    # add colnames to X now
    colnames( X ) <- fam$id[ 1 : n[1] ]

    # successful run
    expect_silent(
        X_all <- draw_geno_fam( X, fam )
    )
    expect_true( is.matrix( X_all ) )
    expect_equal( ncol( X_all ), sum(n) )
    expect_equal( nrow( X_all ), m )
    expect_true( !anyNA( X_all ) )
    expect_true( all( X_all %in% 0L:2L ) )
    expect_equal( colnames( X_all ), fam$id )
    
    # a more explicit toy case with deterministic solution
    # also tests for edge case of a single child
    X <- cbind( rep.int( 0L, m ), rep.int( 2L, m ) )
    colnames( X ) <- letters[1:2]
    X_all_exp <- cbind( X, rep.int( 1L, m ) )
    colnames( X_all_exp ) <- letters[1:3]
    # minimal FAM table for this
    fam <- tibble(
        id = c('a', 'b', 'c'),
        pat = c(0, 0, 'a'),
        mat = c(0, 0, 'b')
    )
    expect_silent(
        X_all <- draw_geno_fam( X, fam )
    )
    expect_equal( X_all, X_all_exp )

    # same homozygous
    x_exp <- c(0L, 2L, 0L)
    X <- cbind( x_exp, x_exp )
    colnames( X ) <- letters[1:2]
    X_all_exp <- cbind( X, x_exp )
    colnames( X_all_exp ) <- letters[1:3]
    expect_silent(
        X_all <- draw_geno_fam( X, fam )
    )
    expect_equal( X_all, X_all_exp )

    # cause a new error by including a parent missing in X
    fam <- tibble(
        id = c('a', 'b', 'c'),
        pat = c(0, 0, 'a'),
        mat = c(0, 0, 'd') # not "b"
    )
    expect_error( draw_geno_fam( X, fam_g2 ) )
})

test_that( "prune_fam works", {
    # create a pedigree with a founder without descendants
    fam <- tibble(
        id = letters[1:4],
        pat = c(0, 0, 0, 'a'),
        mat = c(0, 0, 0, 'b')
    )
    # only want 'd' and its ancestors
    ids <- 'd'
    
    # cause errors on purpose
    expect_error( prune_fam() )
    expect_error( prune_fam( fam ) )
    expect_error( prune_fam( ids = ids ) )

    # successful run
    expect_silent(
        fam_out <- prune_fam( fam, ids )
    )
    expect_equal( fam_out, fam[ -3, ] ) # exclude 3rd row (id="c")

    # more complicated
    # here 'c' mated with both 'd' and 'e'
    fam <- tibble(
        id  = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
        pat = c(  0,   0,   0,   0, 'a', 'c', 'e'),
        mat = c(  0,   0,   0,   0, 'b', 'd', 'c')
    )
    # listed both 'g' (leaf node) and one of its ancestors (to try to trip it up)
    ids <- c('g', 'c')
    expect_silent(
        fam_out <- prune_fam( fam, ids )
    )
    expect_equal( fam_out, fam[ -c(4, 6), ] ) # exclude c('d', 'f')
    
})

test_that( "admix_fam works", {
    # toy example, with interspersed founders and reordered founders in `admix`
    fam <- tibble(
        id  = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
        pat = c(  0,   0, 'a',   0,   0, 'e', 'f'),
        mat = c(  0,   0, 'b',   0,   0, 'd', 'c')
    )
    k_subpops <- 3
    names_subpops <- paste0( 'S', 1:3 )
    admix <- matrix(
        c(
            1.0, 0.0, 0.0, # a
            0.0, 1.0, 0.0, # b
            0.0, 0.0, 1.0, # e
            0.2, 0.3, 0.5  # d
        ),
        byrow = TRUE,
        ncol = k_subpops
    )
    rownames( admix ) <- c('a', 'b', 'e', 'd') # scrambled from FAM
    # expected output
    admix2_exp <- matrix(
        c(
            1.0, 0.0, 0.0, # a
            0.0, 1.0, 0.0, # b
            0.5, 0.5, 0.0, # c
            0.2, 0.3, 0.5, # d
            0.0, 0.0, 1.0, # e
            0.1, 0.15, 0.75, # f
            0.3, 0.325, 0.375 # g
        ),
        byrow = TRUE,
        ncol = k_subpops
    )
    rownames( admix2_exp ) <- fam$id
    # a succcessful run
    expect_silent(
        admix2 <- admix_fam( admix, fam )
    )
    expect_equal( admix2, admix2_exp )
})