library(tibble)

# some tests require additional packages that aren't strictly required
test_BEDMatrix <- FALSE
if (suppressMessages(suppressWarnings(require(BEDMatrix)))) {
    if (suppressMessages(suppressWarnings(require(genio)))) {
        test_BEDMatrix <- TRUE
    }
}

# validates structure of a single chromosome (tibble)
# flexible for cases where breaks are known and cases where they aren't
validate_chr <- function( chr, leng, breaks = NULL, ancs1 = NULL, ancs2 = NULL, ancs_set = NULL, mapped = FALSE ) {
    # decide on which columns are expected
    chr_cols <- c( 'posg', 'anc' )
    # if mapped, "pos" will be added and be first
    if ( mapped )
        chr_cols <- c( 'pos', chr_cols )
    
    # simple tests
    expect_true( is.data.frame( chr ) )
    expect_equal( names( chr ), chr_cols )
    expect_true( is.numeric( chr$posg ) )
    
    # always test total length, other basic properties
    # in our tests, we only know number of rows when we have all breaks, so those go together
    if ( is.null( breaks ) ) {
        # there's always at least one row
        expect_true( nrow( chr ) > 0 )
        # posg are non-negative, strictly increasing, fit in chr length
        expect_true( all( chr$posg > 0 ) )
        expect_true( all( diff( chr$posg ) > 0 ) )
        expect_true( all( chr$posg <= leng ) )
        # last value should always be length
        expect_true( chr$posg[ nrow( chr ) ] == leng )

        if ( mapped ) {
            # repeat tests on basepair position "pos"
            expect_true( is.integer( chr$pos ) )
            expect_true( all( chr$pos >= 1L ) )
            expect_true( all( diff( chr$pos ) > 0L ) )
        }
    } else {
        # some expected values
        n_chunks <- length( breaks ) + 1
        # require exact number of rows
        expect_equal( nrow( chr ), n_chunks )
        # require exact breaks
        expect_equal( chr$posg, c( breaks, leng ) )
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
validate_hap <- function( hap, lengs, ancs_set = NULL, mapped = FALSE ) {
    n_chr <- length( lengs )
    expect_true( is.list( hap ) )
    expect_equal( length( hap ), n_chr )
    # test each chromosome now
    for ( chr in 1 : n_chr ) {
        # number of breaks is unknown, but just test broader expectations
        validate_chr( hap[[ chr ]], lengs[ chr ], ancs_set = ancs_set, mapped = mapped )
    }
}

# validate structure of a diploid individual (both sets of chromosomes)
validate_ind <- function( ind, lengs, ancs_set = NULL, mapped = FALSE ) {
    expect_true( is.list( ind ) )
    expect_equal( names( ind ), c('pat', 'mat') )
    validate_hap( ind$pat, lengs, ancs_set = ancs_set, mapped = mapped )
    validate_hap( ind$mat, lengs, ancs_set = ancs_set, mapped = mapped )
}

# validate structure of a list of diploid individuals
validate_inds <- function( inds, ids, lengs, ancs_set = NULL, mapped = FALSE ) {
    expect_true( is.list( inds ) )
    expect_equal( names( inds ), ids )
    for ( id in ids ) {
        validate_ind( inds[[ id ]], lengs, ancs_set, mapped = mapped )
    }
}

validate_geno <- function( X, m_loci, n_ind, int = TRUE, vals = 0L : 2L ) {
    expect_true( is.matrix( X ) )
    expect_equal( nrow( X ), m_loci )
    expect_equal( ncol( X ), n_ind )
    expect_true( !anyNA( X ) )
    expect_true( all( X %in% vals ) )
    # not true in all cases
    if ( int ) 
        expect_true( is.integer( X ) )
}

test_that( 'indexes_chr, indexes_chr_range work', {
    expect_silent(
        indexes <- indexes_chr( 1:10 )
    )
    # outputs are 1-based
    expect_equal( indexes, 1:10 )
    # retrieve coordinates of one chromosome
    expect_silent( 
        r <- indexes_chr_range( indexes, 9 )
    )
    # in this case start and end are the same (just one index)
    expect_equal( r, c(9, 9) )

    # a more complicated case with a known answer
    indexes_in <- c( 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 5, 6 )
    indexes_out <- c( 3, 7, 8, 11, 12, 13 )
    expect_silent(
        indexes_obs <- indexes_chr( indexes_in )
    )
    expect_equal( indexes_obs, indexes_out )
    # test retrieval too
    expect_silent( 
        r <- indexes_chr_range( indexes_obs, 4 )
    )
    expect_equal( r, c(9, 11) )
    
    # skipped chromosomes have their ends set to zero
    indexes_in <- c( 1, 2, 4 )
    indexes_out <- c( 1, 2, NA, 3 )
    expect_silent(
        indexes_obs <- indexes_chr( indexes_in )
    )
    expect_equal( indexes_obs, indexes_out )
    # test retrieval too, this is actually the trickest case too
    expect_silent( 
        r <- indexes_chr_range( indexes_obs, 4 )
    )
    expect_equal( r, c(3, 3) )
    # and confirm that retrieval returns NA for missing chromosomes
    expect_silent( 
        r <- indexes_chr_range( indexes_obs, 3 )
    )
    expect_true( all( is.na( r ) ) )

    # harder case where first chr is missing
    indexes_in <- 2
    indexes_out <- c( NA, 1 )
    expect_silent(
        indexes_obs <- indexes_chr( indexes_in )
    )
    expect_equal( indexes_obs, indexes_out )
    # test retrieval too
    expect_silent( 
        r <- indexes_chr_range( indexes_obs, 2 )
    )
    expect_equal( r, c(1, 1) )
    
    # cause errors now
    # inputs have to be sorted!
    expect_error( indexes_chr( 2:1 ) )
    # and contiguous
    expect_error( indexes_chr( c( 1, 1, 1, 2, 2, 2, 1 ) ) )
})

test_that( 'indexes_chr_pos works', {
    # come up with a big enough example where we know if we're looking beyond where we should
    pos <- rep.int( NA, 20 )
    # only these have positions, we will only look at those
    chr_start <- 5L
    chr_end <- 10L
    pos[ chr_start : chr_end ] <- ( chr_start : chr_end ) * 100L
    
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 630L, 880L )
    )
    # here we want to see these indexes
    expect_equal( r, c(7, 8) )

    # an example with ties
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 600L, 800L )
    )
    expect_equal( r, c(6, 8) )
    # and related edge cases
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 599L, 800L )
    )
    expect_equal( r, c(6, 8) )
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 601L, 800L )
    )
    expect_equal( r, c(7, 8) )
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 600L, 801L )
    )
    expect_equal( r, c(6, 8) )
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 600L, 799L )
    )
    expect_equal( r, c(6, 7) )

    # an example that doesn't find anything
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 601L, 699L )
    )
    expect_equal( r, as.integer( c(NA, NA) ) )

    # example that should fail because query ends before data does
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 1L, 10L )
    )
    expect_equal( r, as.integer( c(NA, NA) ) )

    # example that should fail because query starts after data does
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, 2000L, 3000L )
    )
    expect_equal( r, as.integer( c(NA, NA) ) )

    # some examples with random data
    pos <- sort( sample.int( 100, 10 ) )
    pos_range <- sort( sample.int( 100, 2 ) )
    chr_start <- 1L
    chr_end <- 10L
    expect_silent(
        r <- indexes_chr_pos( pos, chr_start, chr_end, pos_range[1], pos_range[2] )
    )
    # test all cases explicitly
    indexes <- pos_range[1] <= pos & pos <= pos_range[2]
    if ( all( is.na( r ) ) ) {
        # in this case there should have been no positions in between the two desired ones, let's see if that's true
        expect_true( !any( indexes ) )
    } else {
        expect_true( any( indexes ) )
        # confirm ranges otherwise
        indexes <- which( indexes )
        expect_equal( r[1], min( indexes ) )
        expect_equal( r[2], max( indexes ) )
    }

    # make sure errors are caught when inputs are out of order
    expect_error( indexes_chr_pos( pos, chr_end, chr_start, pos_range[1], pos_range[2] ) )
    expect_error( indexes_chr_pos( pos, chr_start, chr_end, pos_range[2], pos_range[1] ) )
})

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
    mother <- tibble( posg = leng, anc = 'm' )
    father <- tibble( posg = leng, anc = 'f' )

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
    mother <- tibble( posg = c( leng * 1/3, leng ), anc = paste0('m', 1:2) )
    father <- tibble( posg = c( leng * 2/3, leng ), anc = paste0('f', 1:2) )
    # one recombination right in the middle
    breaks <- leng * 1/2
    # here number of breaks is variable because of which parent may be chosen first
    # if mother is first
    child1 <- tibble( posg = c( leng * 1/3, leng * 1/2, leng * 2/3, leng ), anc = c('m1', 'm2', 'f1', 'f2') )
    # if father is first
    child2 <- tibble( posg = c( leng * 1/2, leng ), anc = c('f1', 'm2') )
    # tests!
    expect_silent(
        child <- recomb_chr( breaks, mother, father )
    )
    expect_true( is.data.frame( child ) )
    expect_equal( names( child ), c('posg', 'anc') )
    expect_true( nrow( child ) %in% c(2L, 4L) )
    # explicitly compare both possiblities in full detail
    if ( nrow( child ) == 4L ) {
        expect_equal( child, child1 )
    } else {
        expect_equal( child, child2 )
    }

    # complicate even more, with more ancestors, but with regular spacing for simplicity
    mother <- tibble( posg = leng * (1:10)/10, anc = paste0( 'm', 1:10 ) )
    father <- tibble( posg = leng * (1:10)/10, anc = paste0( 'f', 1:10 ) )
    # new breaks that don't overlap with ancestral breaks
    breaks <- leng * (1:2)/3
    # if mother is picked first, this is what we expect
    child1 <- tibble(
        posg = sort( c( mother$posg, breaks ) ),
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
    expect_equal( names( child ), c('posg', 'anc') )
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
        mother[[ chr ]] <- tibble( posg = lengs[ chr ], anc = 'm' )
        father[[ chr ]] <- tibble( posg = lengs[ chr ], anc = 'f' )
    }

    # successful test
    expect_silent(
        child <- recomb_hap( mother, father )
    )
    validate_hap( child, lengs, ancs_set = c('f', 'm') )
})

test_that( "recomb_map_lengs works", {
    # just apply directly to latest genome version
    map <- recomb_map_hg38
    expect_silent(
        lengs <- recomb_map_lengs( map, name = 'test' )
    )
    expect_true( is.numeric( lengs ) )
    expect_equal( length( lengs ), length( map ) )
    expect_true( all( lengs > 0 ) )
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

test_that( "recomb_last_gen works", {
    # need a 3-generation pedigree
    fam <- tibble(
        id  = c('g1', 'g2', 'g3', 'g4', 'p1', 'p2', 'c1', 'c2'),
        pat = c(  NA,   NA,   NA,   NA, 'g1', 'g3', 'p1', 'p1'),
        mat = c(  NA,   NA,   NA,   NA, 'g2', 'g4', 'p2', 'p2')
    )
    # list IDs of each generation
    ids <- list(
        paste0( 'g', 1:4 ),
        c('p1', 'p2'),
        c('c1', 'c2')
    )
    # initialize founders
    lengs <- c( 50, 100, 150 )
    expect_silent(
        founders <- recomb_init_founders( ids[[1]], lengs )
    )
    # expected ancestral IDs
    ancs_set <- c(
        paste0( ids[[1]], '_pat' ),
        paste0( ids[[1]], '_mat' )
    )
    
    # actual run
    expect_silent(
        inds <- recomb_last_gen( founders, fam, ids )
    )
    # full validation!
    validate_inds( inds, ids[[3]], lengs, ancs_set )
})

test_that( "recomb_map_fix_ends_chr works", {
    # define an incomplete map, toy example only
    # because of default window size, ensure there's at least 10 Mb at each end!
    map_in <- tibble(
        pos  = c( 4000L, 11000000L, 80000000L, 100000000L ),
        posg = c(     0,        10,        70,         80 )
    )
    pos_length <- 101000000L
    # apply code!
    expect_silent( 
        map_out <- recomb_map_fix_ends_chr( map_in, pos_length )
    )
    # make sure it looks good!
    expect_true( is.data.frame( map_out ) )
    expect_equal( names( map_out ), c('pos', 'posg') )
    expect_true( is.integer( map_out$pos ) )
    expect_true( is.numeric( map_out$posg ) )
    # positions are non-negative, strictly increasing, fit in chr length
    # but actually have strict requirements for ends, particularly for $pos
    # pos
    expect_equal( map_out$pos[1], 1L )
    expect_true( all( diff( map_out$pos ) > 0 ) )
    expect_equal( map_out$pos[ nrow( map_out ) ], pos_length )
    # posg
    expect_equal( map_out$posg[1], 0 )
    expect_true( all( diff( map_out$posg ) >= 0 ) ) # these have ties in real data
    # no expectation for end of posg
})

test_that( "recomb_map_simplify_chr works", {
    # create a simple case where all data falls on the perfect line
    # (1Mb to 10Mb)
    map_in <- tibble( pos = ( 1L : 10L ) * 1000000L )
    # use average recombination rate (rule of thumb)
    map_in$posg <- map_in$pos * 1e-6
    # apply!
    expect_silent( 
        map_out <- recomb_map_simplify_chr( map_in )
    )
    # expect the output to be just the first and last rows
    # (good test for handling of ties too, since all errors are zero at least theoretically!)
    expect_equal( map_out, map_in[ c(1L, 10L), ] )

    # a more random example, and also calculate actual final errors and make sure they are controlled
    # make sure this is the same in function and in validation
    tol <- 1e-6 / 2
    # 1Mb to 100Mb in 1Mb increments, and values will be approximately 1e-6 but random (sorted though)
    # this weight reduces noise further, turns out otherwise I wasn't eliminating anything under the desired tolerance.  In this case the weight equal to the tolerance gives good results (only about 1/3 of positions are not eliminated!)
    w <- tol
    map_in <- tibble(
        pos = ( 1L : 100L ) * 1000000L,
        posg = ( ( 1 - w ) * ( 1L : 100L ) + w * 100 * sort( runif( 100L ) ) )
    )
    # apply!
    expect_silent( 
        map_out <- recomb_map_simplify_chr( map_in, tol = tol )
    )
    # results are random but can expect dimensions to be less or equal
    expect_true( nrow( map_out ) <= nrow( map_in ) )
    #message( nrow( map_in ), ' -> ', nrow( map_out ) ) # debugging, making sure problem is non-trivial enough
    # estimate all original positions using new map, to confirm error is small as expected
    pos <- map_in$pos
    posg <- map_in$posg
    expect_silent(
        posg_est <- stats::approx( map_out$pos, map_out$posg, pos )$y
    )
    # make sure there are no NAs (all should be interpolations)
    expect_true( !anyNA( posg_est ) )
    # now actually look at errors
    # a few errors are greater than the tolerance but not by more than 2x
    expect_true( all( abs( posg_est - posg ) < 2 * tol ) )
    
})

test_that( "recomb_map_posg works", {
    # here test all maps we have (all chromosomes, both hg versions provided), to make sure none produce warnings
    for ( chr in 1 : 22 ) {
        # test both maps provided with package 
        for ( map in list( recomb_map_hg38[[ chr ]], recomb_map_hg37[[ chr ]]) ) {
            leng <- max( map$posg )
            # draw some breaks in genetic distance
            posg <- recomb_breaks( leng )
            # to avoid trivial case (zero recombinations), add length to end
            posg <- c( posg, leng )
            # map the breaks to bp positions!
            expect_silent(
                pos <- recomb_map_posg( posg, map )
            )
            expect_equal( length( pos ), length( posg ) )
            expect_true( all( pos > 1 ) )
            expect_true( all( diff( pos ) > 0 ) )
        }
    }
})

test_that( "recomb_map_chr works", {
    # since all maps were thoroughly tested already, here a single random choice suffices
    # pick a random human chromosome
    chr <- sample( 22, 1 )
    # use its map from hg38 (data has been made at this point)
    map <- recomb_map_hg38[[ chr ]]
    leng <- max( map$posg )
    # draw some breaks in genetic distance
    # to avoid trivial case (zero recombinations), add length to end (as proper chr inputs should be anyway)
    chr <- tibble(
        posg = c( recomb_breaks( leng ), leng )
    )
    # add arbitrary ancestors, just for code to work
    ancs_set <- c('a', 'b')
    chr$anc <- rep_len( ancs_set, nrow( chr ) )
    # map the breaks to bp positions!
    expect_silent(
        chr <- recomb_map_chr( chr, map )
    )
    validate_chr( chr, leng, ancs_set = ancs_set, mapped = TRUE )
})

test_that( "recomb_map_hap, recomb_map_ind, recomb_map_inds work", {
    # come up with toy test inputs
    # use real recombination maps, but only a few random chromosomes to speed up tests
    map <- recomb_map_hg38[ sample( 22, 2 ) ]
    # so it's not entirely trivial, construct a child with recombined chromosomes of parents
    # initialize fam table
    fam <- tibble(
        id = c('father', 'mother', 'child'),
        pat = c(NA, NA, 'father'),
        mat = c(NA, NA, 'mother')
    )
    # initialize founder structures
    lengs <- recomb_map_lengs( map )
    expect_silent(
        ancs <- recomb_init_founders( fam$id[1:2], lengs )
    )
    # draw recombined child (returns parents too)
    expect_silent(
        inds <- recomb_fam( ancs, fam )
    )
    # test child only
    ind <- inds$child
    # get one of the two sets of haploid chromosomes (pat, or mat) of individual
    hap <- ind$pat
    
    # expect errors when things are missing
    expect_error( recomb_map_hap( ) )
    expect_error( recomb_map_hap( hap ) )
    expect_error( recomb_map_hap( map = map ) )
    expect_error( recomb_map_ind( ) )
    expect_error( recomb_map_ind( ind ) )
    expect_error( recomb_map_ind( map = map ) )
    expect_error( recomb_map_inds( ) )
    expect_error( recomb_map_inds( inds ) )
    expect_error( recomb_map_inds( map = map ) )

    # proper runs
    expect_silent(
        hap_out <- recomb_map_hap( hap, map )
    )
    # generic validator
    validate_hap( hap_out, lengs, c('father_pat', 'father_mat'), mapped = TRUE )
    # make sure only $pos data was added, i.e. if that were removed the object would be the same as input
    hap_out_cleaned <- lapply( hap_out, function(x) { x$pos <- NULL; return( x ) } )
    expect_equal( hap_out_cleaned, hap )

    # run diploid individual test
    expect_silent(
        ind_out <- recomb_map_ind( ind, map )
    )
    # generic validator
    validate_ind( ind_out, lengs, c('father_pat', 'father_mat', 'mother_pat', 'mother_mat'), mapped = TRUE )
    # make sure only $pos data was added, i.e. if that were removed the object would be the same as input
    ind_out_cleaned <- ind_out
    ind_out_cleaned$pat <- lapply( ind_out_cleaned$pat, function(x) { x$pos <- NULL; return( x ) } )
    ind_out_cleaned$mat <- lapply( ind_out_cleaned$mat, function(x) { x$pos <- NULL; return( x ) } )
    expect_equal( ind_out_cleaned, ind )
    # and make sure haploid $pat processing was the same (it is not random conditioning on breaks)
    expect_equal( ind_out$pat, hap_out )

    # run multiple individuals test
    expect_silent(
        inds_out <- recomb_map_inds( inds, map )
    )
    # generic validator
    validate_inds( inds_out, fam$id, lengs, c('father_pat', 'father_mat', 'mother_pat', 'mother_mat'), mapped = TRUE )
    # parents are quite trivial (will skip additional validations), but make sure child matches that of previous step
    expect_equal( inds_out$child, ind_out )
    # make sure only $pos data was added, i.e. if that were removed the object would be the same as input
    inds_out_cleaned <- lapply( inds_out, function( ind ) {
        ind$pat <- lapply( ind$pat, function(x) { x$pos <- NULL; return( x ) } )
        ind$mat <- lapply( ind$mat, function(x) { x$pos <- NULL; return( x ) } )
        return( ind )
    })
    expect_equal( inds_out_cleaned, inds )

})

test_that( "recomb_haplo_chr works", {
    # make toy chr with breaks
    # nothing random so output is deterministic
    # only columns required are `anc` and `pos`
    chr <- tibble(
        pos = c(1000L, 3000L, 20000L, 50000L, 100000L),
        anc = c('c', 'aaa', 'bb', 'z', 'd')
    )
    # now ancestor data
    # some positions match breaks to test tie behavior (end inclusive)
    pos <- c(2000L, 2500L, 3000L, 4000L, 5000L, 19000L, 30000L, 40000L, 50000L, 60000L, 90000L)
    # ancestries in alphabetical order for clarity
    ancs <- c('aaa', 'bb', 'c', 'd', 'z')
    # filled by column! (default for R but not the way I normally think of it)
    # make sure neighboring ancestries differ at breakpoints to test that correctly
    X <- matrix(
        c(
            1L, 0L, NA, 0L, 0L, 1L, 1L, 0L, NA, 1L, 0L,
            0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L, NA, 1L,
            1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L,
            NA, 1L, 0L, 0L, 1L, 0L, 1L, NA, NA, 0L, 0L,
            1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L
        ),
        nrow = length( pos ),
        ncol = length( ancs )
    )
    colnames(X) <- ancs
    # expected output given breaks
    # NOTE: first block "c" has no SNPs! (important case to test), so start with "aaa"
    # list of ancestries (for second test, but also for my reference)
    ancs_exp <- c('aaa', 'aaa', 'aaa', 'bb', 'bb', 'bb', 'z', 'z', 'z', 'd', 'd')
    x_exp <- c(1L, 0L, NA, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L)
    # run test!
    expect_silent( 
        x_out <- recomb_haplo_chr( chr, X, pos )
    )
    expect_equal( x_out, x_exp )

    # repeat with version that returns ancestries
    expect_silent( 
        data <- recomb_haplo_chr( chr, X, pos, ret_anc = TRUE )
    )
    # gather expected output into list, as expected for single test
    data_exp = list( x = x_exp, anc = ancs_exp )
    expect_equal( data, data_exp )
})

test_that( "recomb_haplo_hap, recomb_haplo_ind, recomb_haplo_inds, recomb_geno_inds, recomb_admix_inds work", {
    # come up with toy test inputs
    # use real recombination maps, but only a few random chromosomes to speed up tests
    n_chr <- 2
    map <- recomb_map_hg38[ sample( 22, n_chr ) ]
    # so it's not entirely trivial, construct a child with recombined chromosomes of parents
    # initialize fam table
    fam <- tibble(
        id = c('father', 'mother', 'child'),
        pat = c(NA, NA, 'father'),
        mat = c(NA, NA, 'mother')
    )
    # need later...
    anc_names <- c( 'father_pat', 'father_mat', 'mother_pat', 'mother_mat' )
    # initialize founder structures
    lengs <- recomb_map_lengs( map )
    expect_silent(
        ancs <- recomb_init_founders( fam$id[1:2], lengs )
    )
    # draw recombined child (returns parents too)
    expect_silent(
        inds <- recomb_fam( ancs, fam )
    )
    # map all breaks to basepair positions, needed for this test
    expect_silent(
        inds <- recomb_map_inds( inds, map )
    )
    # test child only
    ind <- inds$child
    # get one of the two sets of haploid chromosomes (pat, or mat) of individual
    hap <- ind$pat

    # now construct random haplotype input
    haplo <- vector( 'list', n_chr )
    # number of ancestors (x2 because we do haploids here)
    n_ind <- length( ancs ) * 2L
    # number of loci per chr, for toy test
    m_loci <- 10L
    for ( chr in 1L : n_chr ) {
        # draw positions
        pos_chr <- sort( sample.int( max( map[[ chr ]]$pos ), m_loci ) )
        # draw haplotypes
        X_chr <- matrix(
            rbinom( m_loci * n_ind, 1L, 0.5 ),
            nrow = m_loci,
            ncol = n_ind
        )
        # required column names!
        colnames( X_chr ) <- anc_names
        # add to structure, in a list
        haplo[[ chr ]] <- list( X = X_chr, pos = pos_chr )
    }
    
    # run desired tests!
    expect_silent(
        data_hap <- recomb_haplo_hap( hap, haplo )
    )
    # validate
    expect_true( is.list( data_hap ) )
    expect_equal( length( data_hap ), n_chr )
    for ( chr in 1 : n_chr ) {
        x <- data_hap[[ chr ]]
        expect_true( is.integer( x ) )
        expect_equal( length( x ), m_loci )
        expect_true( all( x %in% c(0L, 1L) ) )
    }
    
    # test version with per-pos ancestry
    expect_silent(
        data_hap_anc <- recomb_haplo_hap( hap, haplo, ret_anc = TRUE )
    )
    # validate
    expect_true( is.list( data_hap_anc ) )
    expect_equal( length( data_hap_anc ), n_chr )
    for ( chr in 1 : n_chr ) {
        data_chr <- data_hap_anc[[ chr ]]
        expect_true( is.list( data_chr ) )
        expect_equal( names( data_chr ), c('x', 'anc') )
        x <- data_chr$x
        anc <- data_chr$anc
        expect_true( is.integer( x ) )
        expect_equal( length( x ), m_loci )
        expect_true( all( x %in% c(0L, 1L) ) )
        expect_true( is.character( anc ) )
        expect_equal( length( anc ), m_loci )
        expect_true( all( anc %in% anc_names ) )
    }

    # ind version
    expect_silent( 
        data_ind <- recomb_haplo_ind( ind, haplo )
    )
    expect_true( is.list( data_ind ) )
    expect_equal( names( data_ind ), c('pat', 'mat') )
    # these should be the same!
    expect_equal( data_ind$pat, data_hap )
    
    # ind version with ancestries
    expect_silent( 
        data_ind_anc <- recomb_haplo_ind( ind, haplo, ret_anc = TRUE )
    )
    expect_true( is.list( data_ind_anc ) )
    expect_equal( names( data_ind_anc ), c('pat', 'mat') )
    # these should be the same!
    expect_equal( data_ind_anc$pat, data_hap_anc )

    # inds version
    expect_silent(
        data_inds <- recomb_haplo_inds( inds, haplo )
    )
    expect_true( is.list( data_inds ) )
    expect_equal( names( data_inds ), fam$id )
    # these should be the same!
    expect_equal( data_inds$child, data_ind )

    # inds version with ancestries
    expect_silent(
        data_inds_anc <- recomb_haplo_inds( inds, haplo, ret_anc = TRUE )
    )
    expect_true( is.list( data_inds_anc ) )
    expect_equal( names( data_inds_anc ), fam$id )
    # these should be the same!
    expect_equal( data_inds_anc$child, data_ind_anc )

    ## # quick tests where inputs are missing chromosomes
    ## # only makes sense to delete one that isn't the last one, let's do the first
    ## haplo_del <- haplo
    ## haplo_del[[1]] <- NULL
    ## expect_silent(
    ##     data_inds <- recomb_haplo_inds( inds, haplo_del )
    ## )
    ## expect_true( is.list( data_inds ) )
    ## expect_equal( names( data_inds ), fam$id )
    
    # test function that converts output of previous functions into genotype matrices
    expect_silent(
        X <- recomb_geno_inds( data_inds )
    )
    validate_geno( X, m_loci * n_chr, nrow( fam ) )

    # make sure this works and output is identical when data with ancestries is passed
    expect_silent(
        X2 <- recomb_geno_inds( data_inds_anc )
    )
    expect_equal( X2, X )

    # now test admixture dosages function
    # define ancestry map
    anc_map <- tibble(
        anc = anc_names, # recall these are 4 haplotypes
        pop = c('AFR', 'EUR', 'AFR', 'AFR')
    )
    pops <- c('AFR', 'EUR') # expected output names
    # must use the data with ancestries (expect error otherwise!)
    expect_error( recomb_admix_inds( data_inds, anc_map ) )
    # successful run
    expect_silent(
        X_anc <- recomb_admix_inds( data_inds_anc, anc_map )
    )
    expect_true( is.list( X_anc ) )
    expect_equal( length( X_anc ), length( pops ) )
    expect_equal( names( X_anc ), pops )
    for ( pop in pops ) {
        X_pop <- X_anc[[ pop ]]
        validate_geno( X_pop, m_loci * n_chr, nrow( fam ) )
    }
    
    # simulate a case where each parent has a single ancestry (both of its haplotypes), so the child's dosages are 1 everywhere!
    anc_map <- tibble(
        anc = anc_names, # recall these are 4 haplotypes
        pop = c('AFR', 'AFR', 'EUR', 'EUR')
    )
    # successful run
    expect_silent(
        X_anc <- recomb_admix_inds( data_inds_anc, anc_map )
    )
    expect_true( is.list( X_anc ) )
    expect_equal( length( X_anc ), length( pops ) )
    expect_equal( names( X_anc ), pops )
    for ( pop in pops ) {
        X_pop <- X_anc[[ pop ]]
        validate_geno( X_pop, m_loci * n_chr, nrow( fam ) )
        # for parents, none are 1
        expect_true( all( X_pop[ , 1L:2L ] != 1L ) )
        # for child, all are 1
        expect_true( all( X_pop[ , 3L ] == 1L ) )
    }
    
})

test_that( "bim_add_posg works", {
    # a toy dataset that should succeed
    # chr, pos are minimal required columns
    # random is an additional column, which should be fine to have
    bim <- tibble( chr = 1:22, pos = (22:1) * 10000000, random = 'x' )
    map <- recomb_map_hg38
    expect_silent(
        bim2 <- bim_add_posg( bim, map )
    )
    # comfirm that rest of tibble is unaltered, just posg was added
    bim3 <- bim2
    bim3$posg <- NULL
    expect_equal( bim, bim3 )
    # now check new column
    # in this case we just want numbers, positive values
    expect_true( is.numeric( bim2$posg ) )
    expect_true( min( bim2$posg ) >= 0 )

    # repeat test but with a single valid chromosome, which should work without issues
    bim$chr <- 3
    expect_silent(
        bim2 <- bim_add_posg( bim, map )
    )
    # comfirm that rest of tibble is unaltered, just posg was added
    bim3 <- bim2
    bim3$posg <- NULL
    expect_equal( bim, bim3 )
    # now check new column
    # in this case we just want numbers, positive values
    expect_true( is.numeric( bim2$posg ) )
    expect_true( min( bim2$posg ) >= 0 )
    # because the positions were sorted backwards, all differences must be negative
    # (tests that pos and posg are monotonically increasing)
    expect_true( max( diff( bim2$posg ) ) <= 0 )

    # out-of-range chromosomes are just unmapped with a warning
    bim_bad <- bim
    bim_bad$chr <- 'x'
    # initialize to facilitate comparison to output
    bim_bad$posg <- NA
    expect_warning(
        bim_bad2 <- bim_add_posg( bim_bad, map )
    )
    # in this case output matches input
    expect_equal( bim_bad2, bim_bad )
})

test_that( "pop_recomb_chr works", {
    # some recombination
    G <- 10
    # need single chromosome, but let's pick a realistic one
    map <- recomb_map_hg38[[ 22L ]]
    # sample random positions allowed by the map
    m_loci <- 10
    pos <- sort( sample( max( map$pos ), m_loci ) )
    # and random haplotype data to go with this
    n_ind <- 9
    haps <- matrix( rbinom( m_loci * n_ind, 1, 0.5 ), nrow = m_loci, ncol = n_ind )
    # for subsetting tests
    indexes_loci <- c(2, 4)
    m_loci_subset <- indexes_loci[2] - indexes_loci[1] + 1L
    pos_subset <- pos[ indexes_loci[1] : indexes_loci[2] ]

    # test a trivial single haplotype as input, so answer is deterministic
    # (but must pass as matrix, hence `drop = FALSE`)
    expect_silent( 
        hap_new <- pop_recomb_chr( haps[, 1, drop = FALSE ], pos, map, G )
    )
    # equal as vectors
    expect_equal( hap_new, haps[ , 1 ] )

    # test subsetting with indexes
    expect_silent( 
        hap_new <- pop_recomb_chr( haps[, 1, drop = FALSE ], pos_subset, map, G, indexes_loci = indexes_loci )
    )
    expect_equal( hap_new, haps[ indexes_loci[1] : indexes_loci[2], 1 ] )

    # now proper full dataset, intended binary case
    expect_silent( 
        hap_new <- pop_recomb_chr( haps, pos, map, G )
    )
    expect_equal( length( hap_new ), m_loci )
    expect_true( all( hap_new %in% c(0L, 1L) ) )
    
    # repeat with subset
    expect_silent( 
        hap_new <- pop_recomb_chr( haps, pos_subset, map, G, indexes_loci = indexes_loci )
    )
    expect_equal( length( hap_new ), m_loci_subset )
    expect_true( all( hap_new %in% c(0L, 1L) ) )

    if ( test_BEDMatrix ) {
        # test version that loads the haplotypes using BEDMatrix
        # this requires data that can be written to BEDMatrix, so the weird cases below (higher numbers, characters) don't get tested this way
        
        # first write the exact same simulated haplotypes to file
        name_haps <- tempfile('delete-me-random-test') # output name without extensions!
        write_plink( name_haps, haps, verbose = FALSE )
        # reload as BEDMatrix
        haps_BM <- suppressMessages( suppressWarnings( BEDMatrix( name_haps ) ) )
        
        # now repeat test!
        # the output is of course random, so it will be different, but it should at least be correct
        expect_silent( 
            hap_new <- pop_recomb_chr( haps_BM, pos, map, G )
        )
        expect_equal( length( hap_new ), m_loci )
        expect_true( all( hap_new %in% c(0L, 1L) ) )

        # and subset version
        expect_silent( 
            hap_new <- pop_recomb_chr( haps_BM, pos_subset, map, G, indexes_loci = indexes_loci )
        )
        expect_equal( length( hap_new ), m_loci_subset )
        expect_true( all( hap_new %in% c(0L, 1L) ) )
        
        # delete all three outputs (bed/bim/fam) when done
        delete_files_plink( name_haps )
    }
    
    # a weird artificial case that identifies the ancestor only, by number
    n_ind <- 3L
    haps <- t( matrix( 1 : n_ind, nrow = n_ind, ncol = m_loci ) )
    expect_silent( 
        hap_new <- pop_recomb_chr( haps, pos, map, G )
    )
    expect_equal( length( hap_new ), m_loci )
    expect_true( all( hap_new %in% 1L : n_ind ) )

    # and by character
    haps <- matrix( letters[ haps ], nrow = m_loci, ncol = n_ind )
    expect_silent( 
        hap_new <- pop_recomb_chr( haps, pos, map, G )
    )
    expect_equal( length( hap_new ), m_loci )
    expect_true( all( hap_new %in% letters[ 1L : n_ind ] ) )
})

test_that( "pop_recomb works", {
    # some recombination
    G <- 10
    # simulate from whole human genome!
    map <- recomb_map_hg38
    # sample random positions per chr allowed by the map
    m_loci_per_chr <- 6
    # test case that skips chromosomes, which should be totally fine!
    chrs <- c(2,5,21)
    pos <- NULL
    chr <- NULL
    for ( chr_i in chrs ) {
        map_i <- map[[ chr_i ]]
        pos <- c( pos, sort( sample( max( map_i$pos ), m_loci_per_chr ) ) )
        chr <- c( chr, rep.int( chr_i, m_loci_per_chr ) )
    }
    bim <- tibble( chr = chr, pos = pos )
    # and random haplotype data to go with this
    n_ind_hap <- 10
    m_loci <- m_loci_per_chr * length( chrs )
    haps <- matrix( rbinom( m_loci * n_ind_hap, 1, 0.5 ), nrow = m_loci, ncol = n_ind_hap )
    # number of desired individuals is a different value altogether
    n_ind <- 7

    # a successful run
    expect_silent( 
        X <- pop_recomb( haps, bim, map, G, n_ind, geno = TRUE )
    )
    validate_geno( X, m_loci, n_ind )

    # version that simulates only a subset of genome available
    indexes_loci <- sort( sample.int( m_loci, 2 ) )
    m_loci_subset <- indexes_loci[2] - indexes_loci[1] + 1L
    # repeat earlier test otherwise
    expect_silent( 
        X <- pop_recomb( haps, bim, map, G, n_ind, geno = TRUE, indexes_loci = indexes_loci )
    )
    validate_geno( X, m_loci_subset, n_ind )
    
    if ( test_BEDMatrix ) {
        # test version that loads the haplotypes using BEDMatrix
        
        # first write the exact same simulated haplotypes to file
        name_haps <- tempfile('delete-me-random-test') # output name without extensions!
        write_plink( name_haps, haps, verbose = FALSE )
        # reload as BEDMatrix
        haps_BM <- suppressMessages( suppressWarnings( BEDMatrix( name_haps ) ) )
        
        # now repeat test!
        # the output is of course random, so it will be different, but it should at least be correct
        expect_silent( 
            X <- pop_recomb( haps_BM, bim, map, G, n_ind, geno = TRUE )
        )
        validate_geno( X, m_loci, n_ind )

        # repeat subset version too
        expect_silent( 
            X <- pop_recomb( haps_BM, bim, map, G, n_ind, geno = TRUE, indexes_loci = indexes_loci )
        )
        validate_geno( X, m_loci_subset, n_ind )
        
        # delete all three outputs (bed/bim/fam) when done
        delete_files_plink( name_haps )
    }
    
    # deterministic version, where the population has a single haplotype (since here `geno=TRUE`, everything is doubled)
    haps1 <- haps[ , 1L, drop = FALSE ]
    expect_silent( 
        X <- pop_recomb( haps1, bim, map, G, n_ind, geno = TRUE )
    )
    x_exp <- matrix( 2L * haps1, nrow = m_loci, ncol = n_ind )
    expect_equal( X, x_exp )
    # trick version that just copies the same haplotype twice, result should still be the same
    expect_silent( 
        X <- pop_recomb( cbind( haps1, haps1 ), bim, map, G, n_ind, geno = TRUE )
    )
    expect_equal( X, x_exp )
    
    if ( test_BEDMatrix ) {
        # test version that loads the haplotypes using BEDMatrix
        # test the trick version where there are two identical haplotypes, which should give the same deterministic result as before
        
        # first write the exact same simulated haplotypes to file
        name_haps <- tempfile('delete-me-random-test') # output name without extensions!
        write_plink( name_haps, cbind( haps1, haps1 ), verbose = FALSE )
        # reload as BEDMatrix
        haps_BM <- suppressMessages( suppressWarnings( BEDMatrix( name_haps ) ) )
        
        # now repeat test!
        expect_silent( 
            X <- pop_recomb( haps_BM, bim, map, G, n_ind, geno = TRUE )
        )
        expect_equal( X, x_exp )

        # delete all three outputs (bed/bim/fam) when done
        delete_files_plink( name_haps )
    }
    
    # version that requests haplotypes (non-default)
    # no BEDMatrix version of test since this isn't expected to fail if BEDMatrix geno=TRUE succeeded (given that the non-BEDMatrix geno=FALSE below also succeeds)
    expect_silent( 
        X <- pop_recomb( haps, bim, map, G, n_ind, geno = FALSE )
    )
    validate_geno( X, m_loci, n_ind, vals = 0L : 1L )

    # test non-numeric version
    haps <- matrix( letters[ haps + 1L ], nrow = m_loci, ncol = n_ind_hap )
    expect_silent( 
        X <- pop_recomb( haps, bim, map, G, n_ind, geno = FALSE )
    )
    validate_geno( X, m_loci, n_ind, vals = letters[ 1:2 ], int = FALSE )

    # deterministic version, where the population has a single haplotype (only makes sense if `geno=FALSE`), otherwise everything is doubled
    haps1 <- haps[ , 1L, drop = FALSE ]
    expect_silent( 
        X <- pop_recomb( haps1, bim, map, G, n_ind, geno = FALSE )
    )
    x_exp <- matrix( haps1, nrow = m_loci, ncol = n_ind )
    expect_equal( X, x_exp )
    
    # cause one particular weird error on purpose, make sure it's caught
    bim$chr[1] <- 33 # out of range!
    expect_error( 
        pop_recomb( haps, bim, map, G, n_ind, geno = FALSE )
    )
})

validate_inds_tidy <- function( inds, ids, G, map, founders_only = FALSE ) {
    founders_haps <- c(
        paste0( ids[[1]], '_pat' ),
        paste0( ids[[1]], '_mat' )
    )
    if ( ! founders_only ) {
        expect_true( all( inds$ind %in% ids[[ G ]] ) )
        expect_true( all( inds$parent %in% c('pat', 'mat') ) )
    }
    expect_true( all( inds$anc %in% founders_haps ) )
    expect_true( all( inds$chr %in% 1 : length( map ) ) )
    expect_true( all( inds$start >= 0 ) )
    expect_true( all( inds$start <= inds$end ) )
    # end should be below chr length, but that takes more work to construct, for now meh...
}

test_that( 'tidy_recomb_map_inds, recomb_founder_blocks_inherited work', {
    # this creates toy data to use to test this function
    # simulate from whole human genome!
    map <- recomb_map_hg38
    # simulate the ancestors of one person to 3 generations
    Gp <- 3
    data <- fam_ancestors( Gp )
    fam <- data$fam
    ids <- data$ids
    # initialize founders
    founders <- recomb_init_founders( ids[[ 1 ]], map )
    # draw recombination breaks along pedigree, with coordinates in genetic distance (centiMorgans)
    inds <- recomb_last_gen( founders, fam, ids )
    # map recombination break coordinates to base pairs
    inds <- recomb_map_inds( inds, map )
    
    # this is what we actually want to test
    expect_silent(
        inds2 <- tidy_recomb_map_inds( inds )
    )
    # the number of rows is variable because each corresponds to a recombination event, and those are random too
    # and testing it exactly would amount to just rewriting the function
    # instead let's make sure some basic properties make sense
    validate_inds_tidy( inds2, ids, Gp, map )
    
    # now apply that data to this other function
    expect_silent(
        founders_with_descendants <- recomb_founder_blocks_inherited( inds2 )
    )
    # output is very similar
    validate_inds_tidy( founders_with_descendants, ids, Gp, map, founders_only = TRUE )

    # repeat tests with a more random pedigree with multiple individuals in the last generation
    n <- 20
    data <- sim_pedigree( n, G = Gp )
    fam <- data$fam
    ids <- data$ids
    founders <- recomb_init_founders( ids[[ 1 ]], map )
    inds <- recomb_last_gen( founders, fam, ids )
    inds <- recomb_map_inds( inds, map )
    expect_silent(
        inds2 <- tidy_recomb_map_inds( inds )
    )
    validate_inds_tidy( inds2, ids, Gp, map )
    expect_silent(
        founders_with_descendants <- recomb_founder_blocks_inherited( inds2 )
    )
    validate_inds_tidy( founders_with_descendants, ids, Gp, map, founders_only = TRUE )
})

test_that( "geno_last_gen_admix_recomb works", {
    # some recombination
    G <- 10L
    # simulate all chromosomes between 1 and the max!  Code is unhappy otherwise
    n_chr <- 3L
    chrs <- 1L : n_chr
    # simulate from those chromosomes only
    map <- recomb_map_hg38[ chrs ]
    # sample random positions per chr allowed by the map
    m_loci_per_chr <- 6L
    pos <- NULL
    chr <- NULL
    for ( chr_i in chrs ) {
        map_i <- map[[ chr_i ]]
        pos <- c( pos, sort( sample( max( map_i$pos ), m_loci_per_chr ) ) )
        chr <- c( chr, rep.int( chr_i, m_loci_per_chr ) )
    }
    bim <- tibble( chr = chr, pos = pos )
    # and random haplotype data to go with this
    n_ind_hap <- 10
    m_loci <- m_loci_per_chr * length( chrs )
    haps_afr <- matrix( rbinom( m_loci * n_ind_hap, 1L, 0.5 ), nrow = m_loci, ncol = n_ind_hap )
    haps_eur <- matrix( rbinom( m_loci * n_ind_hap, 1L, 0.5 ), nrow = m_loci, ncol = n_ind_hap )
    anc_haps <- list(
        'AFR' = haps_afr,
        'EUR' = haps_eur
    )
    
    # simulate a very small family with one individual, 2 parents, 4 implicit grandparents
    Gp <- 2L
    data <- fam_ancestors( Gp )
    fam <- data$fam
    ids <- data$ids
    # number of output individuals in this simulation
    n_ind <- length( ids[[ Gp ]] )
    # select ancestries for each of the 4 founder haplotypes (unadmixed)
    # make it non-random to most easily ensure that there's at least one of each ancestry
    founders_anc <- c('AFR', 'EUR', 'AFR', 'EUR')
    # set names of founders with _pat/mat, needed to match recombination structure
    # order is odd but choices were random so that doesn't matter
    names( founders_anc ) <- c(
        paste0( ids[[1]], '_pat' ),
        paste0( ids[[1]], '_mat' )
    )
    
    # a successful run
    expect_silent( 
        data <- geno_last_gen_admix_recomb( anc_haps, bim, map, G, fam, ids, founders_anc )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('X', 'Ls', 'haplos') )
    Ls <- data$Ls
    # NOTE: in none of these cases do we test `haplos` directly, though there are indirect tests through X and Ls so meh
    # genotypes and local ancestry have very similar formats
    validate_geno( data$X, m_loci, n_ind, int = FALSE )
    # local ancestry object
    expect_true( is.list( Ls ) )
    expect_equal( names( Ls ), c('AFR', 'EUR') )
    # test both now, as genotypes
    validate_geno( Ls$AFR, m_loci, n_ind )
    validate_geno( Ls$EUR, m_loci, n_ind )

    # try to break things forcing all ancestors of a single ancestry
    founders_anc_no_admix <- rep.int( 'AFR', length( founders_anc ) )
    # set names of founders with _pat/mat, needed to match recombination structure
    # order is odd but choices were random so that doesn't matter
    names( founders_anc_no_admix ) <- c(
        paste0( ids[[1]], '_pat' ),
        paste0( ids[[1]], '_mat' )
    )
    # repeat test, except for a single ancestry
    expect_silent( 
        data <- geno_last_gen_admix_recomb( anc_haps, bim, map, G, fam, ids, founders_anc_no_admix )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('X', 'Ls', 'haplos') )
    Ls <- data$Ls
    validate_geno( data$X, m_loci, n_ind, int = FALSE )
    expect_true( is.list( Ls ) )
    expect_equal( names( Ls ), 'AFR' )
    validate_geno( Ls$AFR, m_loci, n_ind, vals = 2L )

    # test cases where all haplotypes are 1, so non-1s reveal unset values
    anc_haps_ones <- list(
        AFR = matrix( 1L, nrow = m_loci, ncol = 1L ),
        EUR = matrix( 1L, nrow = m_loci, ncol = 1L )
    )
    expect_silent( 
        data_ones <- geno_last_gen_admix_recomb( anc_haps_ones, bim, map, G, fam, ids, founders_anc )
    )
    expect_true( is.list( data_ones ) )
    expect_equal( names( data_ones ), c('X', 'Ls', 'haplos') )
    Ls <- data_ones$Ls
    # for this data all genotypes must be 2, those that aren't weren't set correctly!
    validate_geno( data_ones$X, m_loci, n_ind, int = FALSE, vals = 2 )
    expect_true( is.list( Ls ) )
    expect_equal( names( Ls ), c('AFR', 'EUR') )
    validate_geno( Ls$AFR, m_loci, n_ind )
    validate_geno( Ls$EUR, m_loci, n_ind )

    # do a test with predicable fixed haplotypes by ancestry, so answer is fixed too
    anc_haps_fixed <- list(
        AFR = matrix( 1L, nrow = m_loci, ncol = 1L ),
        EUR = matrix( 0L, nrow = m_loci, ncol = 1L )
    )
    expect_silent( 
        data <- geno_last_gen_admix_recomb( anc_haps_fixed, bim, map, G, fam, ids, founders_anc )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('X', 'Ls', 'haplos') )
    Ls <- data$Ls
    # genotypes and local ancestry have very similar formats
    validate_geno( data$X, m_loci, n_ind, int = FALSE )
    # local ancestry object
    expect_true( is.list( Ls ) )
    expect_equal( names( Ls ), c('AFR', 'EUR') )
    # test both now, as genotypes
    validate_geno( Ls$AFR, m_loci, n_ind )
    validate_geno( Ls$EUR, m_loci, n_ind )
    # in this case genotypes reflect ancestry exactly, so we expect these equivalences
    expect_equal( data$X, Ls$AFR )
    expect_equal( data$X, 2 - Ls$EUR )
    
    if ( test_BEDMatrix ) {
        # test version that loads the haplotypes using BEDMatrix
        # test only one version with deterministic outputs
        
        # first write the exact same simulated haplotypes to file
        # output names without extensions!
        name_haps_afr <- tempfile('delete-me-random-test_afr')
        name_haps_eur <- tempfile('delete-me-random-test_eur')
        write_plink( name_haps_afr, anc_haps_ones$AFR, verbose = FALSE )
        write_plink( name_haps_eur, anc_haps_ones$EUR, verbose = FALSE )
        # reload as BEDMatrix
        anc_haps_ones_BM <- list(
            AFR = suppressMessages( suppressWarnings( BEDMatrix( name_haps_afr ) ) ),
            EUR = suppressMessages( suppressWarnings( BEDMatrix( name_haps_eur ) ) )
        )
        
        # now repeat test!
        expect_silent( 
            data_ones_BM <- geno_last_gen_admix_recomb( anc_haps_ones_BM, bim, map, G, fam, ids, founders_anc )
        )
        # haplos isn't as before, but the other two objects match!
        expect_equal( data_ones_BM$X, data_ones$X )
        expect_equal( data_ones_BM$Ls, data_ones$Ls )
        
        # delete all three outputs (bed/bim/fam) when done
        delete_files_plink( name_haps_afr )
        delete_files_plink( name_haps_eur )
    }
    
})

