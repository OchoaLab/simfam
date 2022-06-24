recomb_haplo_hap <- function( hap, haplo, ret_anc = FALSE ) {
    if ( missing( hap ) )
        stop( '`hap` is required!' )
    if ( missing( haplo ) )
        stop( '`haplo` is required!' )
    if ( !is.list( hap ) )
        stop( '`hap` must be a list!' )
    if ( !is.list( haplo ) )
        stop( '`haplo` must be a list!' )

    # navigate chromosomes, apply processing to each
    n_chr <- length( hap )
    # create output structure, which is also a list with elements for each chr
    data <- vector( 'list', n_chr )
    # access by index, as names are not required to be present
    for ( chr in 1 : n_chr ) {
        # make sure this chr's haplo data exists
        haplo_chr <- haplo[[ chr ]]
        if ( is.null( haplo_chr ) )
            stop( '`haplo` data missing for chr ', chr, '!' )
        if ( !is.list( haplo_chr ) )
            stop( '`haplo[[chr]]` must be a list!  Failed for chr ', chr )
        # extract components, make sure they're not missing
        X_chr <- haplo_chr$X
        pos_chr <- haplo_chr$pos
        if ( is.null( X_chr ) )
            stop( 'Haplotype matrix element "X" missing from `haplo` data for chr `', chr, '`!' )
        if ( is.null( pos_chr ) )
            stop( 'Position vector element "pos" missing from `haplo` data for chr `', chr, '`!' )
        # let other validations happen inside this function
        # data will have haplotypes and if requested also local ancestry for each position
        data[[ chr ]] <- recomb_haplo_chr( hap[[ chr ]], X_chr, pos_chr, ret_anc = ret_anc )
    }
    return( data )
}
