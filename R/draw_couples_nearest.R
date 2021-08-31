# select parents for next generation randomly pairing only opposite-sex individuals while avoiding relatives, but extremely strongly biased for nearest pairings
# default is to avoid up to second cousins (so 3rd cousins and on is OK)
draw_couples_nearest <- function(
                                 kinship_local,
                                 sex,
                                 cutoff = 1 / 4^3
                                 ) {
    # validate kinship_local
    if ( missing( kinship_local ) )
        stop( '`kinship_local` is required!' )
    if ( !is.matrix( kinship_local ) )
        stop( '`kinship_local` must be a matrix!' )
    n <- nrow( kinship_local )
    if ( ncol( kinship_local ) != n )
        stop( '`kinship_local` must be a square matrix!' )
    # kinship_local should be pedigree-based so let's be strict with ranges
    if ( min( kinship_local ) < 0 )
        stop( '`kinship_local` must be non-negative!' )
    if ( max( kinship_local ) > 1 )
        stop( '`kinship_local` must not exceed 1!' )
    
    # validate sex
    if ( missing( sex ) )
        stop( '`sex` is required!' )
    if ( length( sex ) != n )
        stop( '`sex` must have length (', length( sex ), ') equal to `kinship_local` rows/columns (', n, ')!' )
    if ( !all( sex %in% c( 1, 2 ) ) )
        stop( '`sex` must take only values in: 1, 2!' )
    
    # the maximum number of parents/pairs
    n2 <- floor( n / 2 )
    
    # initially everybody is available
    available_men <- which( sex == 1 )
    available_women <- which( sex == 2 )
    # will store pairings in a 2 x (n/2) matrix, just as in other cases
    parents <- matrix( NA, nrow = 2, ncol = n2 )
    # construct each couple explicitly
    # order of priority is random
    current_pair <- 1 # current pair number (not every iteration advances counter)
    # infinite loop we must break out of
    while ( TRUE ) {
        # stop when we're out of either men or women
        # (if there's at least one of each, then we can pair them)
        if ( length( available_men ) == 0 || length( available_women ) == 0 )
            break
        # draw `current_man` randomly among men
        # unfortunately singletons such as `sample( 3, 1 )` samples from 1:3 instead of only picking 3, so we have to look for singletons and treat them differently
        current_man <- if ( length( available_men ) == 1 ) available_men else sample( available_men, 1 )
        # remove from available men
        available_men <- setdiff( available_men, current_man )
        # calculate kinship from current man to all available women, set threshold
        available_women_unrelated <- available_women[ kinship_local[ available_women, current_man ] < cutoff ]
        # here is a potential failure point, if all other remaining women are too related to current man, will have to skip to next iteration
        # there's no solution involving this man, who has already been removed from `available_men`, which ensures we don't accidentally pick him again, so just move on
        if ( length( available_women_unrelated ) == 0 ) 
            next
        # now calculate coordinate distances, select closest, this is the current woman
        current_woman <- available_women_unrelated[ which.min( abs( available_women_unrelated - current_man ) ) ]
        # remove from available women
        available_women <- setdiff( available_women, current_woman )
        # store both parents
        parents[ 1, current_pair ] <- current_man
        parents[ 2, current_pair ] <- current_woman
        # advance pair counter for next iteration
        current_pair <- current_pair + 1
    }
    # in the end, if not everybody was paired, the final columns will have NAs; if so, remove them
    if ( current_pair <= n2 ) {
        parents <- parents[ , -(current_pair : n2), drop = FALSE ]
    }
    # done, return parents!
    return( parents )
}

