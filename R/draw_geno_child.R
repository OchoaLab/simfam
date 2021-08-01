# most basic element for drawing genotypes in pedigrees!
# inputs are genotype vectors of parents
# returns random genotypes of one or more children (always a matrix)
# NOTES
# - does it assume no NAs?  What happens if there are NAs?
draw_geno_child <- function( X ) {
    # validations
    if ( missing( X ) )
        stop( '`X` is required!' )
    if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    if ( ncol( X ) != 2 )
        stop( '`X` must have two columns!' )
    m_loci <- nrow( X )
    # turn into frequencies before passing to rbinom
    X <- X / 2
    
    # draw each haplotype properly first (from each parent, ignoring the other), then combine
    x <- stats::rbinom( m_loci, 1, X[, 1] ) + stats::rbinom( m_loci, 1, X[, 2] )
    # this is the final genotype of the child
    return( x )
}

