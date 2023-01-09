# simfam 0.0.0.9000 (2021-07-30)

- First commit of redesigned family structure simulation code, meant for public (eventually).
  - Original (scripts in a private repository, not a proper package) was not tested systematically, assumed a rigid setup with every individual paired in the next generation and exactly two children per family, and sex retroactively assigned as needed.
  - New code draws sex first and respects those values while pairing (resulting in more realistic constraints), allows unpaired individuals (fixes cases in old setup where pairing everybody was impossible and restarts were required, which in extreme cases led to infinite loops), variable number of children per family (constrained to target population size, which may vary per generation).
  - Also, new code for calculating kinship and admixture matrices, and drawing genotypes, work for arbitrary FAM tables (i.e. arbitrary pedigrees).  Again, old code was very rigid in assuming non-overlapping generations and fixed family sizes, and accepted a more limited data structure for the same reason.
- This early version is thoroughly tested, but not documented at all.  No functions are exported.

# simfam 0.0.1.9000 (2021-07-31)

- Documented (Roxygen2) and exported main functions: `sim_pedigree`, `kinship_fam`, `admix_fam`, `draw_geno_fam`, `draw_sex`, `prune_fam`.
- Function `sim_pedigree`: made `n` first and only mandatory argument (used to be second), now `G` is second (used to be first) and defaults to `G = length(n)`.
- Fixed import namespaces, other minor changes to pass R checks.

# simfam 0.0.2.9000 (2021-08-03)

- Function `sim_pedigree` now assigns IDs without `g` prefix (format is just `\d+-\d+` with two integers denoting generation and index, separated by a dash).

# simfam 0.0.3.9000 (2021-08-03)

- Function `sim_pedigree` now returns parents of founders as `NA` (used to be `0`).
- All functions that accept a FAM table as input now treat `NA` parents correctly as missing (i.e., those individuals with missing parents are treated as founders), and by default the empty strings ('') and zero (0) are also treated as missing (used to be only `0` was treated as missing).

# simfam 0.0.4.9000 (2021-08-04)

- Added vignette with beautiful examples!

# simfam 0.0.5.9000 (2021-08-04)

- Added README
- Added documentation for package entry (`simfam-package`).
- Fixed wording in DESCRIPTION and vignette.

# simfam 0.0.6.9000 (2021-08-04)

- Function `sim_pedigree` removed `verbose` option (it was a holdout from original code, which could get stuck in some situations; the new code doesn't get stuck).

# simfam 0.0.7.9000 (2021-08-05)

- Function `sim_pedigree` now returns `ids` (ids of IDs separated by generation) among its list elements, after `fam` but before `kinship_local`.

# simfam 0.0.8.9000 (2021-08-05)

- Added function `draw_geno_last_gen` for drawing genotypes for last generation only, of a pedigree with non-overlapping generations, saving lots of memory when the number of generations is large (behavior resembles original function, though internally it's a wrapper around the more general `draw_geno_fam`).

# simfam 0.0.9.9000 (2021-08-06)

- Rewrote core of function `draw_geno_fam` in C++ (using `Rcpp`).
  New version is much faster and uses about half as much memory as the previous pure-R version!

# simfam 0.0.10.9000 (2021-08-06)

- Added function `kinship_last_gen` for calculating kinship for last generation only, of a pedigree with non-overlapping generations, saving lots of memory when the number of generations is large (behavior resembles original function, though internally it's a wrapper around the more general `kinship_fam`).
- Other functions now inherit parameters, the main function being `draw_geno_fam` (for the other `*_fam` functions, which are respectively sources for `*_last_gen` functions).
- Removed obsolete comments.

# simfam 0.0.11.9000 (2021-08-06)

- Added function `admix_last_gen`. same deal as previous `*_last_gen` functions (less coding in practice, memory savings).
- Fixed some minor typos in other functions, added `drop = FALSE` in some necessary cases.

# simfam 0.0.12.9000 (2021-08-06)

- Dropped redundant `draw_` prefix from genotype functions:
  - `draw_geno_fam` -> `geno_fam`
  - `draw_geno_last_gen` -> `geno_last_gen`

# simfam 1.0.0 (2021-08-09)

- First version publicly available on GitHub!
- Changes prepping for first CRAN submission
- Corrected spelling errors
- Clarified bnpsd dependence (older versions don't propagate names, which the vignette relies upon).

# simfam 1.0.1 (2021-08-12)

- More changes prepping for CRAN submission (based on feedback from another new package).
- Added bioRxiv paper reference to description.
- Reset `par()` in vignette examples.

# simfam 1.0.2 (2021-08-31)

- First CRAN submission!
- README
  - Uncommented CRAN installation instructions
  - Corrected GitHub installation instructions to use `build_vignettes = TRUE` instead of `build_opts = c()` (which did not build vignettes anymore).
- Function `draw_couples_nearest` removed unnecessary checks (redundant with unit tests)
- Added file `cran-comments.md`.
- DESCRIPTION corrected title for proper title case

# simfam 1.0.3 (2021-09-27)

- 2nd CRAN submission!
- Tests: adjusted sample sizes in an internal unit test to make a particular `sim_pedigree` error rate ~100x less likely than before (error due to not being able to pair everybody; previously, this unlikely error expected 0.4% of the time actually occurred on CRAN, now it is expected 0.002% of the time).

# simfam 1.1.0.9000 (2022-04-26)

First set of updates for simulating with recombination!

- Updated reference DOI to newest preprint
- Added functions `recomb_fam` and `recomb_init_founders` for simulating recombination breaks for a pedigree!
  - These breaks are on genetic distance only so far.
  - Other nice features (drawing genotypes, calculating kinship, or local ancestry) are not yet present but they are planned.

# simfam 1.1.1.9000 (2022-06-06)

- Added simplified human genome recombination maps (data) and the functions used to create them.
  - Added function `recomb_map_fix_ends_chr` to shift and extrapolate genetic map to chromosome ends.
  - Added function `recomb_map_simplify_chr` to simplify genetic maps by removing rows that can be interpolated to within a desired error.
  - Added and documented data `recomb_map_hg38` and `recomb_map_hg37`, which was created from existing maps processed by the above two functions.

# simfam 1.1.2.9000 (2022-06-15)

- Functions `recomb_fam` and `recomb_init_founders` slight change in input and output formats: each chromosome list now has column `posg` indicating end of recombination block in genetic position (the column used to be called `end`; changed to match notation in recombination map, where `pos` is position in base pairs and `posg` is in genetic distance).
- Function `recomb_init_founders` argument `lengs` may now be a recombination map for simplicity, from which the desired chromosome lengths are extracted, rather than having to extract them in a separate step.
- Added function `recomb_map_inds` to map recombination breaks from genetic positions to base pair coordinates.

# simfam 1.1.3.9000 (2022-06-24)

- Added function `recomb_last_gen`, a wrapper around `recomb_fam` that processes data in discrete generations and returns the recombination breaks/blocks of the final generation only, to reduce memory usage.
  Same analogy of previous `*_last_gen` functions and their corresponding `*_fam` versions.
- Added function `recomb_haplo_inds` to construct the haplotypes of descendant individuals given the haplotypes of the ancestors.
- Added function `recomb_geno_inds` to construct a standard genotype matrix from the haplotypes of individuals (a complex nested list structure).

# simfam 1.1.4.9000 (2022-06-27)

- Added function `fam_ancestors` to construct simple ancestor pedigrees for a single person with a desired number of generations, and automatic names.

# simfam 1.1.5.9000 (2022-06-28)

- Added function `recomb_admix_inds` to produce true population ancestry dosage matrices that parallel genotype matrices, useful for regression models that incorporate local ancestry.

# simfam 1.1.6 (2023-01-06)

- CRAN submission
- Documentation edits
  - Updated hyperlink to raw genetic maps.
  - Fixed typos
  - Escaped more code to simplify spellcheck in the future.
- Updated `cran-comments.md`
