recomb_haplo_ind <- function( ind, haplo, ret_anc = FALSE ) {
    # validations
    if ( missing( ind ) )
        stop( '`ind` is required!' )
    if ( missing( haplo ) )
        stop( '`haplo` is required!' )
    if ( !is.list( ind ) )
        stop( '`ind` must be a list!' )
    if ( !is.list( haplo ) )
        stop( '`haplo` must be a list!' )
    if ( !( 'pat' %in% names( ind ) ) )
        stop( '`ind` must have named element "pat"!' )
    if ( !( 'mat' %in% names( ind ) ) )
        stop( '`ind` must have named element "mat"!' )

    # just apply to both parental copies
    data <- list(
        pat = recomb_haplo_hap( ind$pat, haplo, ret_anc = ret_anc ),
        mat = recomb_haplo_hap( ind$mat, haplo, ret_anc = ret_anc )
    )
    
    return( data )
}
