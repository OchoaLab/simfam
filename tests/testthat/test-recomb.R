library(tibble)

# validates structure of a single chromosome (tibble)
# flexible for cases where breaks are known and cases where they aren't
validate_chr <- function( chr, leng, breaks = NULL, ancs1 = NULL, ancs2 = NULL, ancs_set = NULL ) {
    # simple tests
    expect_true( is.data.frame( chr ) )
    expect_equal( names( chr ), c('end', 'anc') )
    
    # always test total length, other basic properties
    # in our tests, we only know number of rows when we have all breaks, so those go together
    if ( is.null( breaks ) ) {
        # there's always at least one row
        expect_true( nrow( chr ) > 0 )
        # end are non-negative, strictly increasing, fit in chr length
        expect_true( all( chr$end > 0 ) )
        expect_true( all( diff( chr$end ) > 0 ) )
        expect_true( all( chr$end <= leng ) )
        # last value should always be length
        expect_true( chr$end[ nrow( chr ) ] == leng )
    } else {
        # some expected values
        n_chunks <- length( breaks ) + 1
        # require exact number of rows
        expect_equal( nrow( chr ), n_chunks )
        # require exact breaks
        expect_equal( chr$end, c( breaks, leng ) )
    }

    # test anc labels
    if ( !is.null( ancs1 ) && !is.null( ancs2 ) ) {
        expect_true( all( chr$anc == ancs1 ) || all( chr$anc == ancs2 ) )
    } else {
        # if there isn't a precise order to compare against, at least make sure adjacent ancestors don't repeat
        # need more than one row for this to make sense
        # NOTE: this can fail under inbreeding, though we don't have that in our testing examples
        if ( nrow( chr ) > 1 )
            for ( i in 2 : nrow( chr ) )
                expect_true( chr$anc[ i ] != chr$anc[ i-1 ] )
        
        # also test set if provided
        if ( !is.null( ancs_set ) )
            expect_true( all( chr$anc %in% ancs_set ) )
    }
}

# validate structure of a halploid individual (all chromosomes)
# at this level breaks will be unknown, so ancs are also unknown in precision (but set can be known)
validate_hap <- function( hap, lengs, ancs_set = NULL ) {
    n_chr <- length( lengs )
    expect_true( is.list( hap ) )
    expect_equal( length( hap ), n_chr )
    # test each chromosome now
    for ( chr in 1 : n_chr ) {
        # number of breaks is unknown, but just test broader expectations
        validate_chr( hap[[ chr ]], lengs[ chr ], ancs_set = ancs_set )
    }
}

# validate structure of a diploid individual (both sets of chromosomes)
validate_ind <- function( ind, lengs, ancs_set = NULL ) {
    expect_true( is.list( ind ) )
    expect_equal( names( ind ), c('pat', 'mat') )
    validate_hap( ind$pat, lengs, ancs_set = ancs_set )
    validate_hap( ind$mat, lengs, ancs_set = ancs_set )
}

# validate structure of a list of diploid individuals
validate_inds <- function( inds, ids, lengs, ancs_set = NULL ) {
    expect_true( is.list( inds ) )
    expect_equal( names( inds ), ids )
    for ( id in ids ) {
        validate_ind( inds[[ id ]], lengs, ancs_set )
    }
}

test_that( "recomb_breaks works", {
    # try invalid inputs, expect errors
    expect_error( recomb_breaks() )
    expect_error( recomb_breaks( 'a' ) )
    expect_error( recomb_breaks( 1:10 ) )
    expect_error( recomb_breaks( -1 ) )

    # now test valid inputs
    for ( i in 1 : 10 ) {
        # draw a random genetic length (in cM) with two recombinations on average
        leng <- rexp( 1, rate = 1 / 200 )
        expect_silent( 
            breaks <- recomb_breaks( leng )
        )
        expect_true( is.numeric( breaks ) )
        expect_true( all( breaks > 0 ) )
        expect_true( all( breaks < leng ) )
        # breaks are monotonically increasing
        expect_true( all( diff( breaks ) > 0 ) )
    }

    # repeat for a much shorter chromosome to ensure zero-recomb edge case works
    # (NOTE: `all(logical(0)) == TRUE` )
    leng <- 1
    expect_silent( 
        breaks <- recomb_breaks( leng )
    )
    expect_true( is.numeric( breaks ) )
    expect_true( all( breaks > 0 ) )
    expect_true( all( breaks < leng ) )
    # breaks are monotonically increasing
    expect_true( all( diff( breaks ) > 0 ) )
})

test_that( "recomb_chr works", {
    # construct parents who are founders (trivial ancestries)
    leng <- 200 # a larger chromosome
    mother <- tibble( end = leng, anc = 'm' )
    father <- tibble( end = leng, anc = 'f' )

    # test errors
    expect_error( recomb_chr() )
    # add more...

    # test working cases
    
    # explicitly test all small cases (0 to 3 recombinations)
    breaks_fixed <- c(50.3, 100.1, 149.99)
    for ( i in 0 : 3 ) {
        # subset values
        breaks <- breaks_fixed[ seq_len( i ) ]
        # some expected values
        n_chunks <- length( breaks ) + 1
        # which parent is first is random, but otherwise there are only two possible patterns
        ancs1 <- rep_len( c('f', 'm'), n_chunks )
        ancs2 <- rep_len( c('m', 'f'), n_chunks )
        # tests
        expect_silent(
            child <- recomb_chr( breaks, mother, father )
        )
        validate_chr( child, leng, breaks, ancs1, ancs2 )
    }
    
    # now pick random recombination breaks, also more direct joint test of these two functions
    expect_silent(
        breaks <- recomb_breaks( leng )
    )
    # some expected values
    n_chunks <- length( breaks ) + 1
    ancs1 <- rep_len( c('f', 'm'), n_chunks )
    ancs2 <- rep_len( c('m', 'f'), n_chunks )
    # tests
    expect_silent(
        child <- recomb_chr( breaks, mother, father )
    )
    validate_chr( child, leng, breaks, ancs1, ancs2 )

    # a more complicated example, with parents with multiple ancestors (2 each)
    mother <- tibble( end = c( leng * 1/3, leng ), anc = paste0('m', 1:2) )
    father <- tibble( end = c( leng * 2/3, leng ), anc = paste0('f', 1:2) )
    # one recombination right in the middle
    breaks <- leng * 1/2
    # here number of breaks is variable because of which parent may be chosen first
    # if mother is first
    child1 <- tibble( end = c( leng * 1/3, leng * 1/2, leng * 2/3, leng ), anc = c('m1', 'm2', 'f1', 'f2') )
    # if father is first
    child2 <- tibble( end = c( leng * 1/2, leng ), anc = c('f1', 'm2') )
    # tests!
    expect_silent(
        child <- recomb_chr( breaks, mother, father )
    )
    expect_true( is.data.frame( child ) )
    expect_equal( names( child ), c('end', 'anc') )
    expect_true( nrow( child ) %in% c(2L, 4L) )
    # explicitly compare both possiblities in full detail
    if ( nrow( child ) == 4L ) {
        expect_equal( child, child1 )
    } else {
        expect_equal( child, child2 )
    }

    # complicate even more, with more ancestors, but with regular spacing for simplicity
    mother <- tibble( end = leng * (1:10)/10, anc = paste0( 'm', 1:10 ) )
    father <- tibble( end = leng * (1:10)/10, anc = paste0( 'f', 1:10 ) )
    # new breaks that don't overlap with ancestral breaks
    breaks <- leng * (1:2)/3
    # if mother is picked first, this is what we expect
    child1 <- tibble(
        end = sort( c( mother$end, breaks ) ),
        anc = c( paste0('m', 1:4 ), paste0('f', 4:7 ), paste0('m', 7:10 ) )
    )
    # if father is picked first, answer is the same but with parent roles reversed
    child2 <- child1
    # reverse letters in several steps (wish there was a tr/// like perl's here)
    child2$anc <- sub( 'm', 'x', child2$anc )
    child2$anc <- sub( 'f', 'm', child2$anc )
    child2$anc <- sub( 'x', 'f', child2$anc )
    # test starts
    expect_silent(
        child <- recomb_chr( breaks, mother, father )
    )
    expect_true( is.data.frame( child ) )
    expect_equal( names( child ), c('end', 'anc') )
    expect_equal( nrow( child ), nrow( child1 ) )
    # explicitly compare both possiblities in full detail
    if ( child$anc[1] == 'm1' ) {
        expect_equal( child, child1 )
    } else {
        expect_equal( child, child2 )
    }
    
})

test_that( "recomb_hap works", {
    # construct trivial parents for first test
    n_chr <- 22
    mother <- vector( 'list', n_chr )
    father <- vector( 'list', n_chr )
    # draw random chromosome lengths, mean 200 cM
    lengs <- rexp( n_chr, rate = 1 / 200 )
    for ( chr in 1 : n_chr ) {
        # each parent has a single ancestor (itself)
        mother[[ chr ]] <- tibble( end = lengs[ chr ], anc = 'm' )
        father[[ chr ]] <- tibble( end = lengs[ chr ], anc = 'f' )
    }

    # successful test
    expect_silent(
        child <- recomb_hap( mother, father )
    )
    validate_hap( child, lengs, ancs_set = c('f', 'm') )
})

test_that( "recomb_init_founders works", {
    # just try one case that works
    ids <- letters[1:2]
    lengs <- 100 * (1:2)
    expect_silent(
        ancs <- recomb_init_founders( ids, lengs )
    )
    ancs_set <- c(
        paste0( ids, '_pat' ),
        paste0( ids, '_mat' )
    )
    validate_inds( ancs, ids, lengs, ancs_set )
})

test_that( "recomb_fam works", {
    # initialize fam table
    fam <- tibble(
        id = c('father', 'mother', 'child'),
        pat = c(NA, NA, 'father'),
        mat = c(NA, NA, 'mother')
    )
    # initialize founders
    ids <- c( 'mother', 'father' )
    lengs <- c( 50, 100, 150 )
    expect_silent(
        ancs <- recomb_init_founders( ids, lengs )
    )
    # expected ancestral IDs
    ancs_set <- c(
        paste0( ids, '_pat' ),
        paste0( ids, '_mat' )
    )
    # actual run
    expect_silent(
        inds <- recomb_fam( ancs, fam )
    )
    # full validation!
    validate_inds( inds, fam$id, lengs, ancs_set )

    # add another generation for a more challenging test
    fam <- tibble(
        id  = c('g1', 'g2', 'g3', 'g4', 'p1', 'p2', 'c1', 'c2'),
        pat = c(  NA,   NA,   NA,   NA, 'g1', 'g3', 'p1', 'p1'),
        mat = c(  NA,   NA,   NA,   NA, 'g2', 'g4', 'p2', 'p2')
    )
    # initialize founders
    ids <- paste0( 'g', 4:1 ) # list them backwards
    lengs <- c( 50, 100, 150 )
    expect_silent(
        ancs <- recomb_init_founders( ids, lengs )
    )
    # expected ancestral IDs
    ancs_set <- c(
        paste0( ids, '_pat' ),
        paste0( ids, '_mat' )
    )
    # actual run
    expect_silent(
        inds <- recomb_fam( ancs, fam )
    )
    # full validation!
    validate_inds( inds, fam$id, lengs, ancs_set )
})
