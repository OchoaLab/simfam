#include <Rcpp.h>
using namespace Rcpp;

//' Extract last index of every chromosome in a vector of chromosomes
//'
//' This function tests that chromosomes appear in their numerical order (dies with an error otherwise), and returns the last index of each chromosome.
//' 
//' @param chrs Vector of integer chromosomes as it appears in a `bim` table.
//' 
//' @return Vector of end indexes per chromosome.
//' The numerical value of the chromosome is its index in this vector.
//' Chromosomes that were not observed have `NA` values as their end indexes.
//'
//' @examples
//' # the end of chr1 is at index 4, for chr2 it's 7, and for chr3 it's 9
//' indexes_chr( c(1,1,1,1,2,2,2,3,3) )
//'
//' @export
// [[Rcpp::export]]
IntegerVector indexes_chr( IntegerVector chrs ) {
  int m = chrs.length();
  // assuming data is sorted, the last value should be the last chromosome!
  int n_chr = chrs[ m - 1 ];
  // initialize output, starting values are NAs here!
  IntegerVector chr_to_end( n_chr, NA_INTEGER );
  // the last value is the last index too
  chr_to_end[ n_chr - 1 ] = m;
  
  // initialize variables with first element values
  int chr_prev = chrs[ 0 ];
  int i, chr_i;
  for ( i = 1; i < m; i++ ) {
    chr_i = chrs[ i ];
    if ( chr_i > chr_prev ) {
      // we've encountered a new chromosome!  record end of previous one
      // place record in base-0 index (back in R it will be base-1 as desired)
      // indexes are 1-based! (so they're useful in R)
      chr_to_end[ chr_prev - 1 ] = i;
      // now update the new chromosome for next round
      chr_prev = chr_i;
    } else if ( chr_i < chr_prev ) {
      stop( "Chromosomes are not sorted!  Observed '%i' before '%i'", chr_prev, chr_i );
    }
  }
  // all done!
  return chr_to_end;
}
