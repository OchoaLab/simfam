#' Simplified recombination maps for human genomes
#'
#' Human genetic recombination maps for builds 38 (GRCh38/hg38) and 37 (GRCh37/hg19, below suffixed as hg37 for simplicity although technically incorrect).
#' Processed each first with [recomb_map_fix_ends_chr()] to shift and extrapolate to sequence ends, then simplified with [recomb_map_simplify_chr()] to remove all values that can be extrapolated with an error of up to `tol = 0.1`, in order to reduce their sizes and interpolation runtime.
#' Defaults were used, which resulted in extrapolated recombination rates close to and centered around the average of 1e-6 cM/base).
#' Autosomes only.
#'
#' @format
#' A list with 22 elements (autosomes, not named), each a tibble with two columns defining the recombination map at that chromosome:
#' - `pos`: position in base pairs
#' - `posg`: position in centiMorgans (cM)
#'
#' @source
#' Raw genetic maps downloaded from this location prior to above processing:
#' <https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/>
#'
#' Chromosome lengths from:
#' <https://www.ncbi.nlm.nih.gov/grc/human/data>
#' @name recomb_map_hg

#' @rdname recomb_map_hg
"recomb_map_hg38"

#' @rdname recomb_map_hg
"recomb_map_hg37"
