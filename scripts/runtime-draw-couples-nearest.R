library(microbenchmark)

# CONCLUSIONS
# - new version is always faster!
#   - expected in hard cases, as there's no attempt to pair everybody
#   - however, even in easy cases (where everybody is pairable) it's faster!
#     - could have been slower (had this been C code) because new version calculates more distances/cutoffs
#     - but it's more parallel, resulting in speed and also briefer code!

# load scripts to test
source( '../R/draw_couples_nearest.R' )
source( '../R/draw_couples_nearest2.R' )

# a random case, but much bigger than before
n <- 1000
# a random positive definite matrix with all values between 0 and 1
kinship_local <- crossprod( matrix( runif( n*n ), nrow = n, ncol = n ) ) / n
# inspection shows distribution of `kinship_local` is centered around 1/4 and low variance (in one case min was 0.2184223)
# so without changes, all individuals will be too related and classic code will fail completely
# so shrink further so at least half of the data is at threshold
kinship_local <- kinship_local / 4^2

# and make a version that is easy to pair, for a separate test
# in this case all are valid pairs (max value is below cutoff)
kinship_local_easy <- kinship_local / 2

# start with easy case
microbenchmark(
    draw_couples_nearest( kinship_local_easy ),
    draw_couples_nearest2( kinship_local_easy ),
    times = 100L
)
## Unit: milliseconds
##                                       expr      min       lq     mean   median
##   draw_couples_nearest(kinship_local_easy) 20.56078 23.06034 23.67303 23.64667
##  draw_couples_nearest2(kinship_local_easy) 17.66595 17.94360 19.59214 20.28266
##        uq      max neval cld
##  24.26743 41.17967   100   b
##  20.45335 33.93510   100  a 

# now hard case!
microbenchmark(
    draw_couples_nearest( kinship_local ),
    draw_couples_nearest2( kinship_local ),
    times = 100L
)
## Unit: milliseconds
##                                  expr       min         lq       mean
##   draw_couples_nearest(kinship_local) 2580.5268 2744.38677 2823.36547
##  draw_couples_nearest2(kinship_local)   17.6862   18.25259   21.08912
##      median         uq        max neval cld
##  2812.83142 2897.34277 3127.16434   100   b
##    21.25392   21.38493   94.32936   100  a 

# NOTE: inspection in one case actually showed original code just timed out (returned NULL) in this hard case, whereas new code returns quite a decent-looking solution with most people paired (457/500 max pairs).

