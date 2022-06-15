# parents are assumed to be data.frames with 2 columns: posg, anc
# (start is 0 for first chunk, posg from previous row otherwise)
# mat/pat are symmetric (no distinction in practice)
recomb_chr <- function( breaks, mat, pat ) {
    # validate data
    if ( missing( breaks ) )
        stop( '`breaks` is required!' )
    if ( missing( mat ) )
        stop( '`mat` is required!' )
    if ( missing( pat ) )
        stop( '`pat` is required!' )
    # validate breaks
    if ( !is.numeric( breaks ) )
        stop( '`breaks` must be numeric!' )
    if ( any( breaks <= 0 ) )
       stop( '`breaks` must be positive!' ) 
    if ( any( diff( breaks ) <= 0 ) )
        stop( '`breaks` must be strictly monotonically increasing!' )
    # validate mat
    if ( !is.data.frame( mat ) )
        stop( '`mat` must be a data.frame (including tibble)' )
    if ( !all( c('posg', 'anc') %in% names( mat ) ) )
        stop( '`mat` must have column names "posg" and "anc"!' )
    if ( nrow( mat ) < 1 )
        stop( '`mat` must have at least one row!' )
    # validate pat
    if ( !is.data.frame( pat ) )
        stop( '`pat` must be a data.frame (including tibble)' )
    if ( !all( c('posg', 'anc') %in% names( pat ) ) )
        stop( '`pat` must have column names "posg" and "anc"!' )
    if ( nrow( pat ) < 1 )
        stop( '`pat` must have at least one row!' )
    # make sure chromosomes have the same length
    leng <- mat$posg[ nrow( mat ) ]
    if ( pat$posg[ nrow( pat ) ] != leng )
        stop( '`mat` and `pat` chomosomes must be same length!' )
    # and make sure breaks are consistent with this length
    if ( any( breaks > leng ) )
        stop( '`breaks` must be smaller than chromosome length!' )
    
    # randomly select first parent order
    if ( stats::runif(1) > 0.5 ) {
        p1 <- mat
        p2 <- pat
    } else {
        p1 <- pat
        p2 <- mat
    }
    
    # if no recombination happened, desired output is first parent (trivial case)
    if ( length( breaks ) == 0 )
        return( p1 )

    # else, for the code to work most concisely, add length as final break
    breaks <- c( breaks, leng )
    
    # create child by recombining parents (non-trivial case)
    break_last <- 0 # first implicit break, updated as last break as we go
    parent_curr <- 1 # indicate current parent
    # earliest row to consider for each parent
    p1_i_start <- 1
    p2_i_start <- 1
    # create child, blank at first
    child <- NULL
    for ( break_curr in breaks ) {
        # parent to work with for this chunk
        if ( parent_curr == 1 ) {
            p <- p1
            i_start <- p1_i_start
        } else {
            p <- p2
            i_start <- p2_i_start
        }
        # identify rows to copy
        # bottom range first
        while ( p$posg[ i_start ] < break_last )
            i_start <- i_start + 1
        # top of range now
        i_end <- i_start
        while ( p$posg[ i_end ] < break_curr )
            i_end <- i_end + 1

        # extract rows to pass to child for this chunk
        child_tmp <- p[ i_start:i_end, ]
        # since starts are implicit, if they get "shorter" there's no need here to edit that (they start where the previous chunk ends)
        # however, last end must be truncated to recombination position
        child_tmp$posg[ nrow( child_tmp ) ] <- break_curr
        # append to rest of child's chunks
        child <- rbind( child, child_tmp )
        
        # update for next round
        if ( parent_curr == 1 ) {
            p1_i_start <- i_end
            parent_curr <- 2
        } else {
            p2_i_start <- i_end
            parent_curr <- 1
        }
        break_last <- break_curr
    }
    
    # all done!  return child info only
    return( child )
}
