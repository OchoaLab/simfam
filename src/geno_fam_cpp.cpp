#include <Rcpp.h>
using namespace Rcpp;

// draw allele from a single genotype (scalar) of a parent
// maybe it should return IntegerVector so NA behaves (initial test in R space suggested otherwise)

// NOTE: this exports to R space (so I can test in testthat) but not to namespace (as desired!)
// [[Rcpp::export]]
int draw_allele( int x ) {
  // homozygote cases are deterministic
  if ( x == 2 ) {
    return 1;
  } else if ( x == 0 ) {
    return 0;
  } else if ( x == 1 ) {
    // heterozygous case is a coin toss
    // here I draw a random uniform number from R (using its RNG) in scalar mode
    // https://gallery.rcpp.org/articles/random-number-generation/
    if ( R::runif(0,1) > 0.5 ) {
      return 1;
    } else {
      return 0;
    }
  } else if ( x == NA_INTEGER ) {
    // NA always results in NA
    return NA_INTEGER;
  } else {
    // anything else must be fatal!
    stop("A genotype was neither 0, 1, 2, or NA!");
  }
}

// [[Rcpp::export]]
IntegerMatrix geno_fam_cpp(IntegerMatrix X_in, IntegerVector i_founder_in, IntegerVector i_founder_out, IntegerVector i_child, IntegerVector i_pat, IntegerVector i_mat) {
  // gather some important dimensions
  int m_loci = X_in.nrow();
  int n_founder = i_founder_in.length();
  int n_child = i_child.length();
  int n_ind_out = n_founder + n_child;
  
  // initialize output genotype matrix (contains founders and all descendants)
  IntegerMatrix X_out(m_loci, n_ind_out);

  // copy founders but in their right location
  // NOTE: i_founder_in/out have same length by construction, here we just assume it's fine
  int i, j, j_in, j_out;
  for (j = 0; j < n_founder; j++) {
    // get input and output indexes for this founder
    j_in = i_founder_in[ j ];
    j_out = i_founder_out[ j ];
    // copy its genome (works because loci are not reordered)
    X_out( _, j_out ) = X_in( _, j_in );
  }

  // now draw children (really every descendant)
  // NOTE: i_child/pat/mat have same length by construction, here we just assume it's fine
  int j_child, j_pat, j_mat, x_pat, x_mat;
  for (j = 0; j < n_child; j++) {
    // get all indexes relating to this family
    j_child = i_child[ j ];
    j_pat = i_pat[ j ];
    j_mat = i_mat[ j ];
    // navigate every locus now
    for (i = 0; i < m_loci; i++) {
      // extract genotypes and draw alleles for each parent
      // Note: although `draw_allele` handles NAs, they are returned as integers and their sum doesn't come out right
      // so it's best to check before passing values if things are NA, and assume after that that the results will be non-NA and will behave!
      x_pat = X_out(i, j_pat);
      if ( x_pat == NA_INTEGER ) {
	X_out(i, j_child) = NA_INTEGER;
      } else {
	x_mat = X_out(i, j_mat);
	if ( x_mat == NA_INTEGER ) {
	  X_out(i, j_child) = NA_INTEGER;
	} else {
	  // add and store for child
	  X_out(i, j_child) = draw_allele( x_pat ) + draw_allele( x_mat );
	}
      }
    }
  }
  
  // done! return genotype matrix
  return X_out;
}
