library(microbenchmark)

# load scripts to test
source( '../R/update_kinship_child.R' )
source( '../R/update_kinship_children.R' )
source( '../R/update_kinship_children2.R' )

# a random case, but much bigger than before
n <- 1000
# a random positive definite matrix with all values between 0 and 1
# no thresholds are applied so no careful scaling is needed here
kinship <- crossprod( matrix( runif( n*n ), nrow = n, ncol = n ) ) / n
# individuals just paired consecutively
parents <- matrix( 1:n, nrow = 2, ncol = n/2 )

# run test!
microbenchmark(
    update_kinship_children(kinship, parents),
    update_kinship_children2(kinship, parents),
    times = 100L,
    check = 'equal'
)
## Unit: milliseconds
##                                        expr      min       lq     mean   median
##   update_kinship_children(kinship, parents) 75.85423 76.89719 81.01667 78.69742
##  update_kinship_children2(kinship, parents) 63.53177 64.85965 67.29630 65.21414
##        uq       max neval cld
##  79.27206 101.41005   100   b
##  65.73540  99.94985   100  a 

# so the new code is a bit faster, but not game changing
