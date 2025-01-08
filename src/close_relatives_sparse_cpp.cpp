#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::vector<size_t> > close_relatives_sparse_cpp(
							      S4 kinship,
							      double cutoff
							      ) {
  // pull down copies of the key sparse kinship variables
  // unfortunately these have to be R objects initially, can't convert more easily
  IntegerVector dims = kinship.slot( "Dim" );
  NumericVector xs = kinship.slot( "x" );
  IntegerVector is = kinship.slot( "i" );
  IntegerVector ps = kinship.slot( "p" );
  // NOTE: `is` and `ps` are already zero based, as is C!  Keep it that way throughout

  // make native C versions of these things
  size_t n = dims[ 0 ];
  size_t ps_length = ps.length();
  size_t i, j, p, ps_end;
  
  // initialize output, which is a 2D vector
  std::vector< std::vector<size_t> > close_relatives(n, std::vector<size_t>(0));
  
  // navigate structure now
  // easiest to go by column this way
  for ( j = 0; j < ps_length - 1; j++ ) {
    // and now the rows of each column
    ps_end = ps[ j + 1 ];
    for ( p = ps[ j ]; p < ps_end; p++ ) {
      i = is[ p ];
      if ( i != j && xs[ p ] >= cutoff ) {
	// this is a closely related pair, add to structure both ways!
	// only here convert to 1-based in lists, but as indexes they stay 0-based because we're still in C++
	close_relatives[ i ].push_back( j + 1 );
	close_relatives[ j ].push_back( i + 1 );
      }
    }
  }
  
  // that's it, return structure, which gets converted to a list automatically!
  return close_relatives;
}
