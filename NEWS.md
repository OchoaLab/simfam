# simfam 0.0.0.9000 (2021-07-30)

- First commit of redesigned family structure simulation code, meant for public (eventually).
  - Original (private repo hacky scripts, was not a proper package) was not tested systematically, assumed a rigid setup with every individual paired in the next generation and exactly two children per family, and sex retroactively assigned as needed.
  - New code draws sex first and respects those values while pairing (resulting in more realistic constraints), allows unpaired individuals (fixes cases in old setup where pairing everybody was impossible and restarts were required, which in extreme cases led to infinite loops), variable number of children per family (constrained to target population size, which may vary per generation).
  - Also, new code for calculating kinship and admixture matrices, and drawing genotypes, work for arbitrary FAM tables (i.e. arbitrary pedigrees).  Again, old code was very rigid in assuming non-overlapping generations and fixed family sizes, and accepted a more limited data structure for the same reason.
- This early version is thoroughly tested, but not documented at all.  No functions are exported.
