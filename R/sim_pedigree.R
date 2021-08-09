#' Construct a random pedigree
#'
#' Specify the number of individuals per generation, and some other optional parameters, and a single pedigree with those properties will be simulated, where close relatives are never paired, sex is drawn randomly per individual and pairings are strictly across opposite-sex individuals, and otherwise closest individuals (on an underlying 1D geography given by their index) are paired in a random order.
#' Pairs are reordered based on the average of their indexes, where their children are placed (determines their indexes in the 1D geography).
#' The procedure may leave some individuals unpaired in the next generation, and family sizes vary randomly (with a fixed minimum family size) to achieve the desired population size in each generation.
#'
#' @param n The number of individuals per generation.
#' If scalar, the number of generations `G >= 2` must also be specified.
#' Otherwise, the length of `n` is the number of generations.
#' @param G The number of generations (optional).
#' Note `G == 1` is founders only, so it is invalid (there is no pedigree).
#' Must specify a `G >= 2` if `n` is a scalar.
#' If both `G` is specified and `length(n) > 1`, both values must agree.
#' @param sex The numeric sex values for the founders (1L for male, 2L for female).
#' By default they are drawn randomly using [draw_sex()].
#' @param kinship_local The local kinship matrix of the founder population.
#' The default value is half the identity matrix, which corresponds to locally unrelated and locally outbred founders.
#' This "local" kinship is the basis for all kinship calculations used to decide on close relative avoidance.
#' The goal is to make a decision to not pair close relatives based on the pedigree only (and not based on population structure, which otherwise increases all kinship values), so the default value is appropriate.
#' @param cutoff Local kinship values strictly less than `cutoff` are required for pairs.
#' The default value of `1/4^3` corresponds to second cousins, so those and closer relatives are forbidden pairs (but a third cousin pair is allowed).
#' @param children_min The minimum number of children per family.
#' Must be 0 or larger, but not exceed the average number of children per family in each generation (varies depending on how many individuals were left unpaired, but this upper limit is approximately `2 * n[i] / n[i-1]` for generation `i`).
#' The number of children for each given family is first chosen as `children_min` plus a Poisson random variable with parameter equal to the mean number of children per family needed to achieve the desired population size (`n`) minus `children_min`.
#' As these numbers may not exactly equal the target population size, random families are incremented or decremented (respecting the minimum family size) by single counts until the target population size is met.
#' @param full If `TRUE`, part of the return object is a list of local kinship matrices for every generation.
#' If `FALSE` (default), only the local kinship matrix of the last generation is returned.
#'
#' @return A list with these named elements:
#' - `fam`: the pedigree, a tibble in plink FAM format.  Following the column naming convention of the related `genio` package, it contains columns:
#'   - `fam`: Family ID, trivial "fam1" for all individuals
#'   - `id`: Individual ID, in this case a code of format (in regular expression) "(\\d+)-(\\d+)" where the first integer is the generation number and the second integer is the index number (1 to `n[g]` for generation `g`).
#'   - `pat`: Paternal ID.  Matches an `id` except for founders, which have fathers set to `NA`.
#'   - `mat`: Maternal ID.  Matches an `id` except for founders, which have mothers set to `NA`.
#'   - `sex`: integers 1L (male) or 2L (female) which were drawn randomly; no other values occur in these outputs.
#'   - `pheno`: Phenotype, here all 0 (missing value).
#' - `ids`: a list of IDs for each generation (indexed in the list by generation).
#' - `kinship_local`: if `full = FALSE`, the local kinship matrix of the last generation, otherwise a list of local kinship matrices for every generation.
#'
#' @examples
#' # number of individuals for each generation
#' n <- c(15, 20, 25)
#' 
#' # create random pedigree with 3 generations, etc
#' data <- sim_pedigree( n )
#' 
#' # this is the FAM table defining the entire pedigree,
#' # which is the most important piece of information desired!
#' data$fam
#'
#' # the IDs separated by generation
#' data$ids
#' 
#' # bonus: the local kinship matrix of the final generation
#' data$kinship_local
#'
#' @seealso
#' Plink FAM format reference:
#' <https://www.cog-genomics.org/plink/1.9/formats#fam>
#'
#' @export
sim_pedigree <- function(
                         n,
                         G = length( n ),
                         sex = draw_sex( n[1] ),
                         kinship_local = diag( n[1] ) / 2,
                         cutoff = 1 / 4^3,
                         children_min = 1L,
                         full = FALSE
                         ) {
    # check for missing stuff
    if ( missing( n ) )
        stop( '`n` (number of individuals) is required!' )
    
    # allow each generation to potentially have a different size
    # if a scalar is passed, turn to vector:
    if ( length( n ) == 1 ) {
        if ( G == 1 )
            stop( 'Either specify `G >= 2` or `n` with `length(n) >= 2`' )
        n <- rep.int( n, G )
    } else if ( length( n ) != G )
        stop( '`n` must either be scalar or a length-`G` vector!' )

    # code depends on this
    if ( children_min < 0 )
        stop( '`children_min` must be non-negative!' )
    # check other inputs
    # won't repeat things checked in `draw_couples_nearest`, but check against `n` doesn't happen there (and checking that it's a matrix first makes sense)
    if ( !is.matrix( kinship_local ) )
        stop( '`kinship_local` must be a matrix!' )
    if ( nrow( kinship_local ) != n[1] )
        stop( 'Number of individuals in `kinship_local` must equal `n[1]`!' )

    # initialize `fam` tibble with data for founders
    fam <- tibble::tibble(
        fam = 'fam1', # place in desired order, but there aren't families really
        # names of individuals in first generation (founders)
        id = paste0( '1-', 1 : n[1] ),
        pat = NA,
        mat = NA,
        sex = sex,
        pheno = 0
    )
    # create a list with IDs separated by generation, simplifies stuff a lot for users outside!
    ids <- vector('list', G)
    ids[[1]] <- fam$id # save this generation's IDs

    # code requires kinship matrices to have IDs as names, for maximum consistency with FAM table
    rownames( kinship_local ) <- fam$id
    colnames( kinship_local ) <- fam$id
    
    if (full) {
        # this holds all data
        LG <- vector('list', G)
        # store initial values
        LG[[1]] <- kinship_local
    }
    
    for (g in 2:G) {
        # let's pick pairs of parents
        parents <- draw_couples_nearest(kinship_local, sex, cutoff = cutoff)
        n_fam <- ncol( parents )
        # worst-case scenario is everybody is too related so there isn't a single parent and no more generations can be picked
        # just die if it's that bad
        if ( n_fam == 0 )
            stop( 'No parent pairs could be picked because everybody is too related!  To avoid this, you can increase `n` values (preferred), reduce `G`, or reduce `cutoff`.' )
        # reorder for visualization and to preserve ancestry continuity in further generations
        # reorders indexes by their average value (of each parent/column)
        parents <- parents[ , order( colMeans( parents ) ), drop = FALSE ]
        
        # draw how many children each family has
        children_per_fam <- draw_num_children_per_fam( n_fam, n[g], children_min )

        # FAM table for this generation
        fam_g <- tibble::tibble(
            fam = 'fam1',
            # names of individuals in current generation (children)
            id = paste0( g, '-', 1 : n[g] ),
            # expand parents * children_per_fam to make FAM pat/mat vectors
            # NOTE: generation of parents is `g-1` (previous of current)!
            pat = paste0( g-1, '-', rep.int( parents[ 1, ], children_per_fam ) ),
            mat = paste0( g-1, '-', rep.int( parents[ 2, ], children_per_fam ) ),
            # draw sex
            sex = draw_sex( n[g] ),
            pheno = 0
        )
        ids[[g]] <- fam_g$id # save this generation's IDs
        
        # need a different fam for `kinship_fam` below that includes not just children but also their parents (treated as founders too though)
        # extract from master table, make it match current `kinship_local` (of parents) directly!
        fam_p <- fam[ fam$id %in% rownames( kinship_local ), ]
        # set parents as missing now, necessary for code to work
        fam_p$pat <- NA
        fam_p$mat <- NA
        # all set, other parts don't matter for these purposes
        # add children to end
        fam_p <- rbind( fam_p, fam_g )
        
        # calculates kinship for entire `fam_p` provided!
        kinship_local <- kinship_fam( kinship_local, fam_p )
        # subset to keep children only now
        indexes <- rownames( kinship_local ) %in% fam_g$id
        kinship_local <- kinship_local[ indexes, indexes ]
        
        # for next round, overwrite `sex` with current generation's values
        sex <- fam_g$sex
        
        # concatenate to master FAM table
        fam <- rbind( fam, fam_g )
        # local kinship too, if needed
        if (full)
            LG[[g]] <- kinship_local
    }

    # return list of lists, but with better names
    return(
        list(
            fam = fam,
            ids = ids,
            kinship_local = if ( full ) LG else kinship_local
        )
    )
}
