#' Draw sex values randomly for a list of individuals
#'
#' Each individual has their sex drawn between male and female with equal probability.
#' Sex is encoded numerically following the convention for plink FAM files (see below).
#'
#' @param n The number of individuals.
#'
#' @return The length-`n` vector of integer sex assignments: `1L` corresponds to male, `2L` to female.
#'
#' @examples
#' draw_sex( 10 )
#' 
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
draw_sex <- function(n) {
    if ( missing( n ) )
        stop( '`n` is required!' )
    if ( n < 1 )
        stop( '`n >= 1` is required!' )
    return( sample( c(1L, 2L), n, replace = TRUE ) )
}
