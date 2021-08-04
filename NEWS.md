# simfam 0.0.0.9000 (2021-07-30)

- First commit of redesigned family structure simulation code, meant for public (eventually).
  - Original (private repo hacky scripts, was not a proper package) was not tested systematically, assumed a rigid setup with every individual paired in the next generation and exactly two children per family, and sex retroactively assigned as needed.
  - New code draws sex first and respects those values while pairing (resulting in more realistic constraints), allows unpaired individuals (fixes cases in old setup where pairing everybody was impossible and restarts were required, which in extreme cases led to infinite loops), variable number of children per family (constrained to target population size, which may vary per generation).
  - Also, new code for calculating kinship and admixture matrices, and drawing genotypes, work for arbitrary FAM tables (i.e. arbitrary pedigrees).  Again, old code was very rigid in assuming non-overlapping generations and fixed family sizes, and accepted a more limited data structure for the same reason.
- This early version is thoroughly tested, but not documented at all.  No functions are exported.

# simfam 0.0.1.9000 (2021-07-31)

- Documented (Roxygen2) and exported main functions: `sim_pedigree`, `kinship_fam`, `admix_fam`, `draw_geno_fam`, `draw_sex`, `prune_fam`.
- Function `sim_pedigree`: made `n` first and only mandatory argument (used to be second), now `G` is second (used to be first) and defaults to `G = length(n)`.
- Fixed import namespaces, other minor changes to pass R CMD check.

# simfam 0.0.2.9000 (2021-08-03)

- Function `sim_pedigree` now assigns IDs without `g` prefix (format is just `\d+-\d+` with two integers denoting generation and index, separated by a dash).

# simfam 0.0.3.9000 (2021-08-03)

- Function `sim_pedigree` now returns parents of founders as `NA` (used to be `0`).
- All functions that accept a FAM table as input now treat `NA` parents correctly as missing (i.e., those individuals with missing parents are treated as founders), and by default the empty strings ('') and zero (0) are also treated as missing (used to be only `0` was treated as missing).

# simfam 0.0.4.9000 (2021-08-04)

- Added vignette with beautiful examples!
