library(simfam)

# construct genotypes for parents
G <- 3
#n <- c( 8, 9, 11 )
#n <- c( 9, 10, 11 )
#n <- c( 10, 11, 12 )
#n <- c( 12, 13, 14 )
n <- c( 14, 15, 16 )
#n_reps <- 10000
n_reps <- 100000

# start counting

n_err <- 0
for ( i in 1 : n_reps ) {
    # and FAM table
    out <- try( sim_pedigree( n, sex = rep_len( c(1L, 2L), n[1] ) ), silent = TRUE )
    if ( class( out ) == 'try-error' )
        n_err <- n_err + 1
}

message( 'Error rate: ', n_err / n_reps )

# time Rscript error-rate.R
# n_reps <- 10000
# Error rate: 0.0041 # 1m1.359s # n <- c( 8, 9, 11 )
# Error rate: 0.0048 # 1m1.324s
# Error rate: 0.0022 # 1m1.230s # n <- c( 9, 10, 11 )
# Error rate: 0.0011 # 1m2.663s # n <- c( 10, 11, 12 )
# Error rate: 1e-04  # 1m3.152s # n <- c( 12, 13, 14 )
# Error rate: 1e-04  # 1m3.662s
# Error rate: 3e-04  # 1m3.017s # (removed unnecessary `X` construction, reduces runtime only)
# Error rate: 0      # 1m3.714s # n <- c( 14, 15, 16 )
# n_reps <- 100000
# Error rate: 2e-05  # 10m33.288s
