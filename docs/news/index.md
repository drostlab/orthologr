# Changelog

## `orthologr` version 0.4.3

### New Function

- [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  enables `diamond deepclust` functionalities for clustering protein
  sequences. Multiple input `.fasta` files or `.fasta.gz` are supported.
- [`deepclust_realign()`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
  enables `diamond deepclust` functionalities for realigning protein
  sequences from the output of
  [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md).
  Multiple input `.fasta` files or `.fasta.gz` are supported.
- [`deepclust_annotate()`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
  enables `diamond deepclust` functionalities for annotating protein
  sequences from the output of
  [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  based on user provided names for each proteome file (or the name of
  the proteome file if not provided).
  [`deepclust_annotate()`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
  can also be used to annotate the output of
  [`deepclust_realign()`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md).
  Multiple input `.fasta` files or `.fasta.gz` are supported.

### Miscellaneous

- Unit tests for
  [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md),
  [`deepclust_realign()`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
  and
  [`deepclust_annotate()`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
  have been added to ensure the correct functionality of the new
  features.
- Unit tests for
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md),
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  and
  [`divergence_map()`](https://drostlab.github.io/orthologr/reference/divergence_map.md)
  functions have been created.

## `orthologr` version 0.4.2

### New Function

- [`diamond()`](https://drostlab.github.io/orthologr/reference/diamond.md),
  [`set_diamond()`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
  [`diamond_best()`](https://drostlab.github.io/orthologr/reference/diamond_best.md)
  and
  [`diamond_rec()`](https://drostlab.github.io/orthologr/reference/diamond_rec.md)
  enable a massive speed-up in pairwise sequence alignment
  functionalities.

### New Features

- a logo is now set for `orthologr`.
- DIAMOND2 is now used by default in
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) and
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  unless blast is specified (`aligner = "blast"`).

## `orthologr` version 0.4.1

### New Features

- the
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  and
  [`divergence_map()`](https://drostlab.github.io/orthologr/reference/divergence_map.md)
  functions now include the parameter `n_quantile`, which enables users
  to choose the number of quantiles to generate for the
  `divergence map`. This could allow users to get a higher-resolution
  `divergence map` if `n_quantile` is greater than 10. Alternatively,
  this can resolve the issue of empty divergence strata when deciling
  the dNdS values for closely related organisms with dNdS = 0 for over
  10% of the genes.

## `orthologr` version 0.4.0

### New Functions

- new function
  [`check_annotation()`](https://drostlab.github.io/orthologr/reference/check_annotation.md)
  helps to detect corrupt GFF or GTF annotation files and removes such
  outlier lines

### New Features

- the
  [`generate_ortholog_tables()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables.md)
  and
  [`retrieve_longest_isoforms()`](https://drostlab.github.io/orthologr/reference/retrieve_longest_isoforms.md)
  now include the new
  [`check_annotation()`](https://drostlab.github.io/orthologr/reference/check_annotation.md)
  function to capture corrupt GFF or GTF files and fix them

- adding a new argument `of_path` to
  [`orthofinder2()`](https://drostlab.github.io/orthologr/reference/orthofinder2.md)
  to allow users to specify their own path to their locally installed
  `orthofinder` executable

- adding new argument `task` to
  [`map_generator_lnc()`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md)
  and
  [`orthologs_lnc()`](https://drostlab.github.io/orthologr/reference/orthologs_lnc.md)
  to allow users to use the full `blastn` range provided by
  [blast_nucleotide_to_nucleotide()](https://drostlab.github.io/metablastr/reference/blast_nucleotide_to_nucleotide.html)

- adding new argument `path` to
  [`map_generator_lnc()`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md)
  to allow users to specify their local installation path of BLAST

## `orthologr` version 0.3.0

### New Functions

- new function
  [`plot_pairwise_orthologs()`](https://drostlab.github.io/orthologr/reference/plot_pairwise_orthologs.md)
  allows users to plot pairwise orthologs tables for multiple pairwise
  comparisons

- new function
  [`retrieve_core_orthologs()`](https://drostlab.github.io/orthologr/reference/retrieve_core_orthologs.md)
  allows users to retrieve a core set of orthologous gene loci from
  several pairwise ortholog tables

- new functions
  [`generate_ortholog_tables()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables.md)
  and
  [`generate_ortholog_tables_all()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md)
  allow users to generate ortholog tables by gene locus and splice
  varaint for a set of species

- new function
  [`retrieve_longest_isoforms_all()`](https://drostlab.github.io/orthologr/reference/retrieve_longest_isoforms_all.md)
  allows users to specify folders and retrieve the longest splice
  variants for all proteomes stored in a folder

- new functions
  [`translate_cds_to_protein()`](https://drostlab.github.io/orthologr/reference/translate_cds_to_protein.md)
  and
  [`translate_cds_to_protein_all()`](https://drostlab.github.io/orthologr/reference/translate_cds_to_protein_all.md)
  which translate coding sequences into amino acid sequences for single
  or multiple files

### New Features

- in
  [`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
  the default value of `delete_corrupt_cds` changed from
  `delete_corrupt_cds = TRUE` to `delete_corrupt_cds = FALSE` to be
  consistent with
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) and
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)

- the
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  function received a new argument `dnds_est.method` which now allows
  users to select different dNdS estimation methods when running
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  (suggested by Momir Futo)

- the
  [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  function now allows to change the `eval` argument which wasn’t passed
  down to the
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
  call within the function (Many thanks to Momir Futo)

- the function `map.generator()` was renamed to
  [`map_generator_dnds()`](https://drostlab.github.io/orthologr/reference/map_generator_dnds.md)
  to be more consistent with the notation of other functions

- the function `map.generator.lnc()` was renamed to
  [`map_generator_lnc()`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md)
  to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to
  [`divergence_map()`](https://drostlab.github.io/orthologr/reference/divergence_map.md)
  to be more consistent with the notation of other functions

- the function `DivergenceMap()` was renamed to
  [`divergence_map()`](https://drostlab.github.io/orthologr/reference/divergence_map.md)
  to be more consistent with the notation of other functions

- the function `orthologs.lnc` was renamed to `orthologs_lnc` to be more
  consistent with the notation of other functions

- the function `OF2CoreOrthologs()` was renamed to
  [`orthofinder2_retrieve_core_orthologs()`](https://drostlab.github.io/orthologr/reference/orthofinder2_retrieve_core_orthologs.md)

### Removed Functions

- the function `advanced_blast()` is not supported anymore and thus is
  not available to users anymore (please consult the [metablastr
  package](https://github.com/drostlab/metablastr) in case you need this
  functionality)

- the function `advanced_makedb()` is not supported anymore and thus is
  not available to users anymore (please consult the [metablastr
  package](https://github.com/drostlab/metablastr) in case you need this
  functionality)

- the function `blast.nr()` is not supported anymore and thus is not
  available to users anymore (please consult the [metablastr
  package](https://github.com/drostlab/metablastr) in case you need this
  functionality)

- the function `delta.blast()` is not supported anymore and thus is not
  available to users anymore (please consult the [metablastr
  package](https://github.com/drostlab/metablastr) in case you need this
  functionality)

- the function `ProteinOrtho()` is not supported anymore and thus is not
  available to users anymore

## `orthologr` version 0.2.0

### New Features

- new function
  [`retrieve_longest_isoforms()`](https://drostlab.github.io/orthologr/reference/retrieve_longest_isoforms.md)
  which enables retrieval of the longest isoforms from a proteome file
  and save results as fasta file for downstream analyses
- new function `OF2CoreOrthologs()` to retrieve core orthologs across
  multiple species from Orthofinder2 output
- new function
  [`extract_features()`](https://drostlab.github.io/orthologr/reference/extract_features.md):
  Helper function to extract gene loci and splice variant IDs from GFF
  files
- new function
  [`filter_best_hits()`](https://drostlab.github.io/orthologr/reference/filter_best_hits.md):
  Helper function to select best BLAST hit based on minimum evalue
- new function
  [`generate_ortholog_tables()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables.md):
  Generate ortholog tables by gene locus and splice varaint

## `orthologr` version 0.1.0

### New Features

- [`read.cds()`](https://drostlab.github.io/orthologr/reference/read.cds.md)
  now trimms corrupted CDS (= CDS not divisible by 3) when
  `delete_corrupt_cds = FALSE` is specified
- the new default value of argument `delete_corrupt_cds` in
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) is
  now `FALSE`. Thus, given the new trimming feature in
  [`read.cds()`](https://drostlab.github.io/orthologr/reference/read.cds.md),
  corrupted CDS equences will be trimmed before being translated.

## `orthologr` version 0.0.5

### Bug fixes

- The default setting of the `BLAST` argument `max_target_seqs 1` was
  removed from
  [`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md)
  and
  [`blast_rec()`](https://drostlab.github.io/orthologr/reference/blast_rec.md)
  due to the misunderstood functionality of the `BLAST` argument (See
  details [here](https://doi.org/10.1093/bioinformatics/bty833) and
  [\#9](https://github.com/drostlab/orthologr/issues/9) ; Many thanks to
  [@armish](https://github.com/armish))

## `orthologr` version 0.0.4

### New Features

- Users can now control via the new `delete_corrupt_cds` argument in
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) and
  related downstream functions whether or not corrupted input coding
  sequences shall be removed prior to dN/dS inference. In case corrupted
  CDS exist, the
  [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
  function will now generate a separate fasta file which stores all
  corrupted CDS so that they can be investigated. See issue
  [\#8](https://github.com/drostlab/orthologr/issues/8) for details.

#### Function updates

- [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`read.cds()`](https://drostlab.github.io/orthologr/reference/read.cds.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`cds2aa()`](https://drostlab.github.io/orthologr/reference/cds2aa.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`set_blast()`](https://drostlab.github.io/orthologr/reference/set_blast.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`blast()`](https://drostlab.github.io/orthologr/reference/blast.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

- [`blast_rec()`](https://drostlab.github.io/orthologr/reference/blast_rec.md)
  receives new argument `delete_corrupt_cds` to remove corrupted input
  coding sequences (`delete_corrupt_cds` is set to `TRUE` as default)

## `orthologr` version 0.0.3

- Fixing internal path bug that caused that wrong pal2nal paths were
  generated when using multiple sequence aligners -\> see issue
  <https://github.com/drostlab/orthologr/issues/5> (Many thanks to
  Dr. Mario López-Pérez)

## `orthologr` version 0.0.2

- Fixing a major bug that caused KaKs_Calculator to not be able to
  correctly parse the kaks computation output (Many thanks to Hongyi Li
  who spotted the bug and found a solution).
