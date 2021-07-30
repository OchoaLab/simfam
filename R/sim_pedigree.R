# NOTE: unlike original `sim_children_generations_kinship`, here we'll ignore population kinship (only local kinship matters)
#
# input are genotypes and kinship of founders, and number of generations required
# simulates two children per couple to have a stable population
# and at each generation there is inbreeding avoidance!
# however, other than very close relatives, sim is strongly biased for pairing individuals with the most similar ancestry (via their coordinate and repeated seriation runs to reorder upon admixture).
# this way population structure is preserved across generations
# note that G==1 is founders
# NOTE: this is kinship-only version, no genotype drawing
# NOTE: because not everybody is guaranteed to be paired, some people without children may (probably will) appear in the FAM table (excluding the trivial case of the last generation, which is left unpaired)
sim_pedigree <- function(
                         G,
                         n,
                         kinship_local = diag( n[1] ) / 2,
                         sex = draw_sex( n[1] ),
                         children_min = 1L,
                         full = FALSE,
                         cutoff = 1 / 4^3,
                         verbose = TRUE
                         ) {
    # check for missing stuff
    if ( missing( G ) )
        stop( '`G` (number of generations) is required!' )
    if ( missing( n ) )
        stop( '`n` (number of individuals) is required!' )
    # code depends on this
    if ( G < 2 )
        stop( '`G` should be at least 2!' )
    if ( children_min < 0 )
        stop( '`children_min` must be non-negative!' )
    # check other inputs
    # won't repeat things checked in `draw_couples_nearest`, but check against `n` doesn't happen there (and checking that it's a matrix first makes sense)
    if ( !is.matrix( kinship_local ) )
        stop( '`kinship_local` must be a matrix!' )
    if ( nrow( kinship_local ) != n[1] )
        stop( 'Number of individuals in `kinship_local` must equal `n[1]`!' )

    # allow each generation to potentially have a different size
    # if a scalar is passed, turn to vector:
    if ( length( n ) == 1 ) {
        n <- rep.int( n, G )
    } else if ( length( n ) != G )
        stop( '`n` must either be scalar or a length-`G` vector!' )

    # initialize `fam` tibble with data for founders
    fam <- tibble(
        fam = 'fam1', # place in desired order, but there aren't families really
        # names of individuals in first generation (founders)
        id = paste0( 'g1-', 1 : n[1] ),
        pat = 0,
        mat = 0,
        sex = sex,
        pheno = 0
    )
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
        if (verbose)
            message('g = ', g)
        # let's pick pairs of parents
        if (verbose)
            message('draw_couples_nearest')
        parents <- draw_couples_nearest(kinship_local, sex, cutoff = cutoff)
        ###DEBUG
        if ( anyNA( parents ) )
            stop('parents (1) had NAs')
        n_fam <- ncol( parents )
        # worst-case scenario is everybody is too related so there isn't a single parent and no more generations can be picked
        # just die if it's that bad
        if ( n_fam == 0 )
            stop( 'No parent pairs could be picked because everybody is too related!' )
        # reorder for visualization and to preserve ancestry continuity in further generations
        # reorders indexes by their average value (of each parent/column)
        parents <- parents[ , order( colMeans( parents ) ), drop = FALSE ]
        ###DEBUG
        if ( anyNA( parents ) )
            stop('parents (2) had NAs')
        
        # draw how many children each family has
        if (verbose)
            message('draw_num_children_per_fam')
        children_per_fam <- draw_num_children_per_fam( n_fam, n[g], children_min )
        ###DEBUG
        if ( anyNA( children_per_fam ) )
            stop('children_per_fam had NAs')

        # FAM table for this generation
        fam_g <- tibble(
            fam = 'fam1',
            # names of individuals in current generation (children)
            id = paste0( 'g', g, '-', 1 : n[g] ),
            # expand parents * children_per_fam to make FAM pat/mat vectors
            # NOTE: generation of parents is `g-1` (previous of current)!
            pat = paste0( 'g', g-1, '-', rep.int( parents[ 1, ], children_per_fam ) ),
            mat = paste0( 'g', g-1, '-', rep.int( parents[ 2, ], children_per_fam ) ),
            # draw sex
            sex = draw_sex( n[g] ),
            pheno = 0
        )
        
        # need a different fam for `kinship_fam` below that includes not just children but also their parents (treated as founders too though)
        # extract from master table, make it match current `kinship_local` (of parents) directly!
        fam_p <- fam[ fam$id %in% rownames( kinship_local ), ]
        # set parents as missing now, necessary for code to work
        fam_p$pat <- 0
        fam_p$mat <- 0
        # all set, other parts don't matter for these purposes
        # add children to end
        fam_p <- rbind( fam_p, fam_g )
        
        if (verbose)
            message('kinship_fam (local)')
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
            kinship_local = if ( full ) LG else kinship_local
        )
    )
}
