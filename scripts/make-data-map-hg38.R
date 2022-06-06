library(readr)
library(testthat)
library(usethis)
library(ochoalabtools)
#library(simfam)
source( '../R/recomb_map_fix_ends_chr.R' ) # hack while it's not exported
source( '../R/recomb_map_simplify_chr.R' )

# params
hgV <- 38
tol <- 1e-1
n_chr <- 22
chrs <- 1 : n_chr

plot_maps <- function( name, data ) {
    # figure out plot ranges
    pos_range <- range( sapply( data, function(x) range( x$pos ) ) )
    posg_range <- range( sapply( data, function(x) range( x$posg ) ) )
    fig_start( name, height = 6, width = 6 )
    # make plot area
    plot( NA, xlim = pos_range, ylim = posg_range, xlab = 'Position (bp)', ylab = 'Genetic distance (cM)' )
    for ( chr in chrs ) {
        i <- which( chr == chrs )
        data_chr <- data[[ chr ]]
        x <- data_chr$pos
        y <- data_chr$posg
        lines( x, y, col = i )
        # hard to label these things well, but marking extreme points seems good
        text( range(x), range(y), labels = chr, cex = 0.3 )
    }
    fig_end()
}

# recall where we started
dir_orig <- getwd()

# load general recombination map
# autosomes only
setwd( '~/dbs/genetic-maps/' )
setwd( paste0( 'hg', hgV ) )

data <- vector( 'list', n_chr )
for ( chr in chrs ) {
    file <- paste0( 'chr', chr, '.gz' )
    data[[ chr ]] <- read_tsv( file, col_names = c('pos', 'posg'), col_types = 'id' )
}
# load lengths in bp
pos_lengths <- read_tsv( 'lengths.txt', col_types = 'ci' )

# we know there are issues with the maps, fix ends
for ( chr in chrs ) {
    message( 'recomb_map_fix_ends_chr: chr', chr )
    # copy down data
    map <- data[[ chr ]]
    pos_length <- pos_lengths$length[ chr ]
    # apply fix
    expect_silent( 
        map <- recomb_map_fix_ends_chr( map, pos_length )
    )
    
    # standard validations
    expect_true( is.data.frame( map ) )
    expect_equal( names( map ), c('pos', 'posg') )
    expect_true( is.integer( map$pos ) )
    expect_true( is.numeric( map$posg ) )
    # positions are non-negative, strictly increasing, fit in chr length
    # but actually have strict requirements for ends, particularly for $pos
    # pos
    expect_equal( map$pos[1], 1L )
    expect_true( all( diff( map$pos ) > 0 ) )
    expect_equal( map$pos[ nrow( map ) ], pos_length )
    # posg
    expect_equal( map$posg[1], 0 )
    expect_true( all( diff( map$posg ) >= 0 ) ) # these have ties in real data
    # no expectation for end of posg

    # if all clear, overwrite in big structure
    data[[ chr ]] <- map
}

# make plot of "fixed" data, to compare to original (made by another script)
plot_maps( 'plot-fixed', data )

# second step is simplifications, as these tables are huge in a probably unnecessary way
for ( chr in chrs ) {
    message( 'recomb_map_simplify_chr: chr', chr )
    # copy down data
    map_in <- data[[ chr ]]
    pos_length <- pos_lengths$length[ chr ]
    # apply fix
    expect_silent( 
        map <- recomb_map_simplify_chr( map_in, tol = tol )
    )

    # expect a reduction
    expect_true( nrow( map ) <= nrow( map_in ) )
    message( nrow( map_in ), ' -> ', nrow( map ) )
    # estimate all original positions using new map, to confirm error is small as expected
    pos <- map_in$pos
    posg <- map_in$posg
    expect_silent(
        posg_est <- stats::approx( map$pos, map$posg, pos )$y
    )
    # make sure there are no NAs (all should be interpolations)
    expect_true( !anyNA( posg_est ) )
    # now actually look at errors
    # mean error should certainly be under the tolerance
    expect_true( mean( abs( posg_est - posg ) ) < tol )
    # a few errors are greater than the tolerance but not by more than 3x
    expect_true( all( abs( posg_est - posg ) < 3 * tol ) )
    
    # if all clear, overwrite in big structure
    data[[ chr ]] <- map
}

# make plot of "fixed" data, to compare to original (made by another script)
plot_maps( paste0( 'plot-fixed-simplified-', tol ), data )

# when all done, save as data object
# go back to project
setwd( dir_orig )
# rename object for better loading
if ( hgV == 38 ) {
    recomb_map_hg38 <- data
    use_data( recomb_map_hg38, overwrite = TRUE )
} else if ( hgV == 37 ) {
    recomb_map_hg37 <- data
    use_data( recomb_map_hg37, overwrite = TRUE )
}
# hg38 notes
# 30M   unsimplified
# 23M     simplified tol = 1e-6/2  62m17.656s ideapad
# 14M     simplified tol = 1e-5   173m35.426s viiiaR5
#  7.9M   simplified tol = 1e-4   282m20.460s viiiaR5
#  3.2M   simplified tol = 1e-3   368m03.242s viiiaR5
#  926K   simplified tol = 1e-2   377m47.968s viiiaR5
#  145K   simplified tol = 1e-1   383m01.328s viiiaR5
# hg37
#  145K   simplified tol = 1e-1   393m28.098s viiiaR5



# all below comments are hg38 only

## # report rates given this transform
## rates <- sapply( data, function( map ) map$posg[2]/map$pos[2] )
## rates_end <- sapply( data, function( map ) { m <- nrow(map); ( map$posg[m] - map$posg[m-1])/( map$pos[m] - map$pos[m-1]) } )

# 1Mb window, looks ok but some are too large
## rates
##  [1] 2.614333e-06 1.188656e-06 1.388770e-06 1.048527e-06 9.456073e-07
##  [6] 2.636357e-06 1.797883e-06 1.359231e-06 1.936096e-06 6.574553e-07
## [11] 1.724730e-06 2.378421e-06 9.457458e-07 2.493364e-06 6.340934e-06
## [16] 2.174737e-06 3.712705e-06 3.205198e-06 4.407401e-06 4.940805e-06
## [21] 0.000000e+00 1.178499e-07
# 2Mb window, max got smaller a bit, but still have a zero!
##  [1] 2.173985e-06 1.299818e-06 1.780631e-06 1.074230e-06 1.678681e-06
##  [6] 2.986391e-06 1.721501e-06 2.160381e-06 2.892618e-06 1.420254e-06
## [11] 1.924131e-06 2.306618e-06 1.284211e-06 3.439610e-06 5.346668e-06
## [16] 2.107185e-06 3.224838e-06 2.444792e-06 3.890799e-06 3.540473e-06
## [21] 0.000000e+00 1.565753e-06
# 3Mb window, no zeroes but first one is very small
##  [1] 1.911194e-06 1.694020e-06 2.341765e-06 9.820897e-07 2.469818e-06
##  [6] 2.645928e-06 1.802910e-06 2.321849e-06 2.943351e-06 1.721178e-06
## [11] 2.206627e-06 2.303406e-06 1.346848e-06 3.220490e-06 5.609154e-06
## [16] 2.112919e-06 2.911847e-06 2.615271e-06 3.780135e-06 3.254609e-06
## [21] 9.629827e-08 1.995264e-06
# 10Mb window, all much tighter around 2e-6 (min 1.5e-6, max 3.5e-6), will keep this!
##  [1] 2.173288e-06 2.601961e-06 2.474924e-06 2.576055e-06 2.273779e-06
##  [6] 2.240991e-06 1.996913e-06 2.139267e-06 2.338478e-06 2.206776e-06
## [11] 1.881684e-06 2.300819e-06 2.098755e-06 2.174835e-06 3.501910e-06
## [16] 2.359202e-06 2.586357e-06 2.959741e-06 3.054318e-06 2.765094e-06
## [21] 1.550994e-06 2.288018e-06

# end rates are tight wih 10Mb window (only one used), around 1e-6 again, min 1.5e-6, max 3e-6)
## rates_end
##  [1] 2.053891e-06 2.113561e-06 2.091141e-06 2.535335e-06 2.269811e-06
##  [6] 2.012442e-06 2.486103e-06 1.567666e-06 2.615426e-06 2.427888e-06
## [11] 2.088652e-06 2.809929e-06 2.374710e-06 1.915847e-06 2.921564e-06
## [16] 3.030498e-06 2.603084e-06 2.281555e-06 2.912121e-06 2.646369e-06
## [21] 2.171426e-06 2.387881e-06

# genetic lengths
## lengs <- sapply( data, function(x) max( x$posg ) )
## # for info
## pos_min <- sapply( data, function(x) min( x$pos ) )
## posg_min <- sapply( data, function(x) min( x$posg ) )
## pos_max <- sapply( data, function(x) max( x$pos ) )

## lengs # after (0,0) fix but no end fix
##  [1] 286.39996 268.87343 223.40920 214.72172 204.13616 192.24612 187.29104
##  [8] 168.46078 166.45500 181.24778 158.59134 174.77827 164.42222 161.59439
## [15] 211.23213 134.11805 128.91276 117.74303 108.50945 108.49133  78.21894
## [22] 107.63084
# after end fix, all small differences as expected
##  [1] 286.46492 268.96787 223.71848 214.99317 204.70236 192.37441 187.32060
##  [8] 168.55934 166.91319 181.32318 158.61430 174.98367 164.47638 161.90772
## [15] 211.30600 134.84890 129.35541 118.00766 108.60213 108.82497  78.28352
## [22] 107.69553


