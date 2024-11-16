#include <Rcpp.h>
using namespace Rcpp;

// this one is internal so we won't document
// chr_start/end: indexes in `pos` that belong to desired chromosome
// pos_start/end: values of `pos` we are interested in
// [[Rcpp::export]]
IntegerVector indexes_chr_pos( IntegerVector pos, int chr_start, int chr_end, int pos_start, int pos_end ) {
  // want to find index range of overlap, in original indexes

  // validate inputs quickly
  if ( chr_start > chr_end )
    stop( "Input `chr_start > chr_end` has wrong order!" );
  if ( pos_start > pos_end )
    stop( "Input `pos_start > pos_end` has wrong order!" );
  
  // all values from outside are 1-based!
  // convert these to 0-based for simplicity, to simplify access to `pos`
  chr_start--;
  chr_end--;
  // keep `pos_start/end` as 1-based so they match actual values in `pos`

  // initialize output
  IntegerVector inds( 2, NA_INTEGER );

  // one simple edge case is when 
  
  // start from the bottom, want to find first index whose position is equal or greater than the desired position
  // note i is 0-based now
  int i;
  int i_start = -1;
  for ( i = chr_start; i <= chr_end; i++ ) {
    if ( pos[i] >= pos_start ) {
      // remember this index, then stop looking!
      i_start = i;
      break;
    }
  }
  // if no such position was found, return NA coordinates
  if ( i_start == -1 ) return inds;
  // otherwise continue

  // now look for the end of the range, start looking at the start that we found (more likely to be closer to this than to the end of the chromosome)
  // initialize this one to -2, because it's possible if the condition below is triggered at i=i_start=0 for i_end=-1, we want to distinguish that case from the condition never being triggered!
  int i_end = -2;
  for ( i = i_start; i <= chr_end; i++ ) {
    if ( pos[i] > pos_end ) {
      // if we've exceeded the position we were looking for, we want the previous one as the one we go with
      i_end = i - 1;
      break;
    }
  }
  // if the condition was never triggered, then `pos_end` exceeds all positions in the chromosome, so we want the last index to be what we keep in this case
  if ( i_end == -2 ) i_end = chr_end;
  
  // NOTE: since `pos[ i_start ]` could have exceeded pos_start by a lot, it could have exceeded pos_end too!  In that case `i_end = i_start - 1` gets set right away and the loop breaks, let's see if that happened and return the NAs in that case
  // 
  if ( i_end < i_start ) return inds;

  // otherwise it should all be good!
  // coordinates we have found are 0-based, but we need them 1-based outside!!!
  inds[0] = i_start + 1;
  inds[1] = i_end + 1;
  return inds;
}
