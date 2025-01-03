#include <string>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

// https://en.cppreference.com/w/cpp/container/vector
// https://teuder.github.io/rcpp4everyone_en/160_s3_s4.html#s4-class

// [[Rcpp::export]]
S4 kinship_fam_sparse_cpp(
			  S4 kinship,
			  IntegerVector pars1_R,
			  IntegerVector pars2_R,
			  CharacterVector ids_new_R
			  ) {
  // pull down copies of the key sparse kinship variables
  // unfortunately these have to be R objects initially, can't convert more easily
  IntegerVector dims = kinship.slot( "Dim" );
  NumericVector xs_R = kinship.slot( "x" );
  IntegerVector is_R = kinship.slot( "i" );
  IntegerVector ps_R = kinship.slot( "p" );
  // NOTE: `is` and `ps` are already zero based, as is C!  Keep it that way throughout

  // make native C versions of these things
  size_t n = dims[ 0 ];
  size_t is_length = is_R.length();
  size_t ps_length = ps_R.length();
  std::vector<double> xs = as< std::vector<double> >( xs_R );
  std::vector<size_t> is = as< std::vector<size_t> >( is_R );
  std::vector<size_t> ps = as< std::vector<size_t> >( ps_R );
  size_t j, p, p2, indiv, ps1_start, ps1_end, ps1_length, ps1_filled, par1, par2;
  double f;
  
  // for row searches, precalculate the columns at each position `p` of `is`
  std::vector<size_t> js( is_length );
  for ( j = 0; j < ps_length - 1; j++ )
    for ( p = ps[ j ]; p < ps[ j + 1 ]; p++ )
      js[ p ] = j;
  
  // process each individual
  size_t n_indiv = pars1_R.length();
  for ( indiv = 0; indiv < n_indiv; indiv++ ) {
    par1 = pars1_R[ indiv ];
    par2 = pars2_R[ indiv ];

    // NOTE: par1 != par2 is assumed but not tested
    
    // make sure par1 < par2
    // (j reused temporarily, though it will be rewritten later as needed)
    if ( par1 > par2 ) {
      j = par2;
      par2 = par1;
      par1 = j;
    }

    // search column par1
    ps1_start = ps[ par1 ];
    ps1_end = ps[ par1 + 1 ];
    ps1_length = ps1_end - ps1_start;
    std::vector<size_t> is2( ps1_length );
    std::vector<double> xs2( ps1_length );
    for ( p = ps1_start, p2 = 0; p < ps1_end; p++, p2++ ) {
      // these vectors are new values we'll be adding in the end
      is2[ p2 ] = is[ p ];
      xs2[ p2 ] = xs[ p ];
    }
    
    // search rows for par1, starting at current position (ps1_end)
    for ( p = ps1_end; p < is_length; p++ ) {
      if ( is[ p ] == par1 ) {
	is2.push_back( js[ p ] );
	xs2.push_back( xs[ p ] );
      }
    }

    // search column par2
    ps1_start = ps[ par2 ];
    ps1_end = ps[ par2 + 1 ];
    ps1_length = ps1_end - ps1_start;
    // resize is2/xs2 instead of using push_back, since we know the length of the new data in advance
    ps1_filled = is2.size();
    is2.resize( ps1_filled + ps1_length );
    xs2.resize( ps1_filled + ps1_length );
    // also find and set inbreeding
    f = 0;
    for ( p = ps1_start, p2 = ps1_filled; p < ps1_end; p++, p2++ ) {
      // these vectors are new values we'll be adding in the end
      is2[ p2 ] = is[ p ];
      xs2[ p2 ] = xs[ p ];
      // test if parents are related, set their kinship as the individual's inbreeding too
      if ( is[ p ] == par1 ) {
	f = xs[ p ];
      }
    }

    // search rows for par2, starting at current position (ps1_end)
    for ( p = ps1_end; p < is_length; p++ ) {
      if ( is[ p ] == par2 ) {
	is2.push_back( js[ p ] );
	xs2.push_back( xs[ p ] );
      }
    }
    
    // almost done, first sort data by `i`
    // because we need parallel sorting, let's do it this way
    // https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
    // first create a simple index vector
    ps1_length = is2.size();
    std::vector<size_t> index( ps1_length );
    for ( p = 0 ; p < ps1_length ; p++ )
      index[ p ] = p;
    // this sorts index
    sort( index.begin(), index.end(),
	  [&]( const size_t& a, const size_t& b ) {
	    return ( is2[a] < is2[b] );
	  }
	  );
    // now apply index to vectors, but need to store them in new objects
    std::vector<size_t> is2_sorted( ps1_length );
    std::vector<double> xs2_sorted( ps1_length );
    for ( p = 0 ; p < ps1_length ; p++ ) {
      p2 = index[ p ];
      is2_sorted[ p ] = is2[ p2 ];
      // good time to halve all new x values, as they must be
      xs2_sorted[ p ] = xs2[ p2 ] / 2;
    }
    // done, replace objects
    is2 = is2_sorted;
    xs2 = xs2_sorted;
    
    // now add values where both parents are related to a given person, instead of listing twice
    // NOTE: there's always at least two individuals, since this person is related to both of their parents, so this won't fail
    // NOTE: it's ok to increment p even when there was a deletion, because there are never more than two overlaps (because there are only two parents). Skipping the unnecessary test should be faster!
    for ( p = 1 ; p < ps1_length ; p++ ) {
      if ( is2[ p ] == is2[ p - 1 ] ) {
	// since is2 is now sorted, repeats must be adjacent!
	// add up value here onto previous one
	xs2[ p - 1 ] += xs2[ p ];
	// now remove entry from both vectors
	// (need this awkward notation because input must be an iterator)
	is2.erase( is2.begin() + p );
	xs2.erase( xs2.begin() + p );
	// decrement this too
	ps1_length--;
      }
    }

    // the above loops never set self kinship directly in any form, since the new individual isn't even on the matrix, so we have to do it as an extra step.  However, we have already captured the kinship between parents if any.
    // the self kinship is always the last value of the new column, so it can be added after sorting and halving
    n++; // now's the best time to increment size of output matrix
    ps1_length++;
    is2.push_back( n - 1 ); // remember it's 0-based!
    xs2.push_back( ( 1 + f ) / 2 );
    
    // now apply updates to existing data, which gets used in the next round
    xs.insert( xs.end(), xs2.begin(), xs2.end() );
    is.insert( is.end(), is2.begin(), is2.end() );
    // for the new ps, we just added one more column, and the increment is the number of new values
    ps.push_back( ps[ ps_length - 1 ] + ps1_length );
    // though only used internally, update this too so it works in the next round
    std::vector<size_t> js2(ps1_length, n-1); // remember it's 0-based!
    js.insert( js.end(), js2.begin(), js2.end() );

    // grow lengths as we go
    // ps_length++;
    // is_length += is2.size();
    // calculate new lengths
    is_length = is.size();
    ps_length = ps.size();
  }
  
  // when done, create a new object, otherwise the original gets edited by reference which is bad practice!
  IntegerVector dims_new_R = IntegerVector::create( n, n );
  NumericVector xs_new_R = wrap( xs );
  IntegerVector is_new_R = wrap( is );
  IntegerVector ps_new_R = wrap( ps );

  // process IDs now
  // get original ones first
  List dimnames_R = kinship.slot( "Dimnames" );
  CharacterVector ids_R = dimnames_R[0];
  std::vector<std::string> ids = as< std::vector<std::string> >( ids_R );
  std::vector<std::string> ids_new = as< std::vector<std::string> >( ids_new_R );
  // append new ones into new object
  ids.insert( ids.end(), ids_new.begin(), ids_new.end() );
  // form new vector and list
  ids_new_R = wrap( ids );
  List dimnames_new_R = List::create( ids_new_R, ids_new_R );
  
  // finalize new object to return!
  S4 kinship_new("dsCMatrix");
  kinship_new.slot( "Dim" ) = dims_new_R;
  kinship_new.slot( "Dimnames" ) = dimnames_new_R;
  kinship_new.slot( "x" ) = xs_new_R;
  kinship_new.slot( "i" ) = is_new_R;
  kinship_new.slot( "p" ) = ps_new_R;
  return kinship_new;
}
