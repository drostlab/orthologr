# Package index

## Sequence Search (`DIAMOND2`)

Run DIAMOND2 homology searches and build DIAMOND databases.

- [`diamond()`](https://drostlab.github.io/orthologr/reference/diamond.md)
  : Perform a DIAMOND2 search

- [`diamond_best()`](https://drostlab.github.io/orthologr/reference/diamond_best.md)
  : Perform a DIAMOND2 best hit search

- [`diamond_rec()`](https://drostlab.github.io/orthologr/reference/diamond_rec.md)
  : Perform a DIAMOND2 reciprocal best hit (RBH) search

- [`set_diamond()`](https://drostlab.github.io/orthologr/reference/set_diamond.md)
  :

  Create a DIAMONDable database with `diamond makedb`

## Sequence Search (BLAST+)

Run BLAST+ homology searches, filter hits, and build BLAST databases.

- [`blast()`](https://drostlab.github.io/orthologr/reference/blast.md) :
  Perform a BLAST+ search
- [`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md)
  : Perform a BLAST+ best hit search
- [`blast_rec()`](https://drostlab.github.io/orthologr/reference/blast_rec.md)
  : Perform a BLAST+ reciprocal best hit (RBH) search
- [`set_blast()`](https://drostlab.github.io/orthologr/reference/set_blast.md)
  : Create a BLASTable database with makeblastdb
- [`filter_best_hits()`](https://drostlab.github.io/orthologr/reference/filter_best_hits.md)
  : Helper function to select best BLAST hit based on minimum evalue

## Sequence Alignment

Pairwise and multiple sequence alignments, including codon-aware
alignments.

- [`pairwise_aln()`](https://drostlab.github.io/orthologr/reference/pairwise_aln.md)
  : Compute Pairwise Alignments
- [`multi_aln()`](https://drostlab.github.io/orthologr/reference/multi_aln.md)
  : Compute Multiple Sequence Alignments
- [`codon_aln()`](https://drostlab.github.io/orthologr/reference/codon_aln.md)
  : Compute a Codon Alignment

## Orthology Inference

Infer orthologous genes between pairs or sets of species.

- [`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
  : Main Orthology Inference Function
- [`select_orthologs()`](https://drostlab.github.io/orthologr/reference/select_orthologs.md)
  : Select orthologs based on either gene locus or splice variant

## Sequence Clustering

Cluster sequences with `diamond deepclust` and post-process clusters.

- [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  :

  Cluster sequences using `diamond deepclust`

- [`deepclust_annotate()`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
  :

  Annotate `deepclust` output with source file metadata

- [`deepclust_realign()`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
  :

  Realign sequences within `deepclust` clusters

## dNdS / Substitution Rates

Compute, filter, and import pairwise synonymous (dS) and non-synonymous
(dN) substitution rates, and generate multi-species dNdS maps.

- [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) :
  Compute dNdS values for two organisms

- [`compute_dnds()`](https://drostlab.github.io/orthologr/reference/compute_dnds.md)
  : Compute dNdS Values For A Given Pairwise Alignment

- [`substitutionrate()`](https://drostlab.github.io/orthologr/reference/substitutionrate.md)
  : Internal function for dNdS computations

- [`filter_dNdS()`](https://drostlab.github.io/orthologr/reference/filter_dNdS.md)
  : Filter dNdS values

- [`read.dnds.tbl()`](https://drostlab.github.io/orthologr/reference/read.dnds.tbl.md)
  :

  Import a dNdS table generated with `dNdS`

- [`map_generator_dnds()`](https://drostlab.github.io/orthologr/reference/map_generator_dnds.md)
  : Generate dNdS Maps Between a Query Organism and Multiple Subject
  Organisms

## Divergence Stratigraphy

Perform divergence stratigraphy to assign genes to evolutionary
divergence strata based on dNdS values.

- [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  : Perform Divergence Stratigraphy
- [`divergence_map()`](https://drostlab.github.io/orthologr/reference/divergence_map.md)
  : Sort dNdS Values Into Divergence Strata

## Core Ortholog Tables

Generate, retrieve, and import pairwise and multi-species ortholog
tables.

- [`generate_ortholog_tables()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables.md)
  : Generate ortholog tables by gene locus and splice varaint

- [`generate_ortholog_tables_all()`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md)
  : Generate ortholog tables by gene locus and splice varaint for a set
  of species

- [`import_ortholog_tables_all()`](https://drostlab.github.io/orthologr/reference/import_ortholog_tables_all.md)
  :

  Importing output pairwise orthologs tables generated with
  `generate_ortholog_tables_all`

- [`retrieve_core_orthologs()`](https://drostlab.github.io/orthologr/reference/retrieve_core_orthologs.md)
  : Retrieve a core set of orthologs from pairwise ortholog tables

- [`filter_core_set()`](https://drostlab.github.io/orthologr/reference/filter_core_set.md)
  : Helper function to extract a core set of orthologs

## OrthoFinder2 Interface

Run OrthoFinder2 and retrieve core orthologs from its output.

- [`orthofinder2()`](https://drostlab.github.io/orthologr/reference/orthofinder2.md)
  : Interface function to Orthofinder2
- [`orthofinder2_retrieve_core_orthologs()`](https://drostlab.github.io/orthologr/reference/orthofinder2_retrieve_core_orthologs.md)
  : Retrieve core orthologs across multiple species from Orthofinder2
  output

## lncRNA Orthology

Infer orthologous long non-coding RNAs (lncRNAs) and summarise results.

- [`orthologs_lnc()`](https://drostlab.github.io/orthologr/reference/orthologs_lnc.md)
  : Orthology Inference of lncRNAs
- [`lnc_map_core_orthologs()`](https://drostlab.github.io/orthologr/reference/lnc_map_core_orthologs.md)
  : Retrieve a core set of orthologous lncRNAs (deprecated)
- [`lnc_map_counts()`](https://drostlab.github.io/orthologr/reference/lnc_map_counts.md)
  : Count number of orthologous lncRNAs per pairwise species comparison
- [`map_generator_lnc()`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md)
  : Infer orthologous lncRNAs between multiple species
- [`filter_core_set_lnc()`](https://drostlab.github.io/orthologr/reference/filter_core_set_lnc.md)
  : Helper function to extract a core set of orthologous lncRNAs
  (deprecated)

## Promoter Divergence

Estimate divergence of promoter sequences for orthologous genes.

- [`promotor_divergence_estimation()`](https://drostlab.github.io/orthologr/reference/promotor_divergence_estimation.md)
  : Estimate the DNA distance between promotor sequences
- [`promotor_divergence_of_orthologous_genes()`](https://drostlab.github.io/orthologr/reference/promotor_divergence_of_orthologous_genes.md)
  : Compute promotor sequence divergence of orthologous genes

## Visualization

Plot ortholog counts and homology-threshold curves across species
comparisons.

- [`plot_pairwise_orthologs()`](https://drostlab.github.io/orthologr/reference/plot_pairwise_orthologs.md)
  :

  A line plot visualizing the number of pairwise orthologs within a
  ortho table generated with `generate_ortholog_tables_all`

- [`plot_diverse_homology_thresholds()`](https://drostlab.github.io/orthologr/reference/plot_diverse_homology_thresholds.md)
  :

  Diverse line plots visualizing the number of pairwise orthologs within
  a ortho table generated with `generate_ortholog_tables_all` based on
  different sets of homology thresholds.

- [`plot_diverse_homology_thresholds_core_orthologs()`](https://drostlab.github.io/orthologr/reference/plot_diverse_homology_thresholds_core_orthologs.md)
  :

  Diverse line plots visualizing the number of core orthologs within a
  ortho table generated with `generate_ortholog_tables_all` based on
  different sets of homology thresholds.

- [`testCoreOrthoParamsGeneLocus()`](https://drostlab.github.io/orthologr/reference/testCoreOrthoParamsGeneLocus.md)
  :

  Helper function for `plot_diverse_homology_thresholds_core_orthologs`

- [`testCoreOrthoParamsSpliceVariant()`](https://drostlab.github.io/orthologr/reference/testCoreOrthoParamsSpliceVariant.md)
  :

  Helper function for splice variant based
  `plot_diverse_homology_thresholds_core_orthologs`

## Sequence I/O

Read and write genome, proteome, and CDS files.

- [`read.genome()`](https://drostlab.github.io/orthologr/reference/read.genome.md)
  : Read the genome of a given organism
- [`read.proteome()`](https://drostlab.github.io/orthologr/reference/read.proteome.md)
  : Read the proteome of a given organism
- [`read.cds()`](https://drostlab.github.io/orthologr/reference/read.cds.md)
  : Read the CDS of a given organism
- [`write.proteome()`](https://drostlab.github.io/orthologr/reference/write.proteome.md)
  : Save a proteome in fasta format

## Sequence Processing

Translate CDS sequences and retrieve longest isoforms from proteome
files.

- [`transl()`](https://drostlab.github.io/orthologr/reference/transl.md)
  : Translate DNA to Amino Acids
- [`cds2aa()`](https://drostlab.github.io/orthologr/reference/cds2aa.md)
  : Translate CDS file to Amino Acids file
- [`translate_cds_to_protein()`](https://drostlab.github.io/orthologr/reference/translate_cds_to_protein.md)
  : Translate coding sequences into amino acid sequences
- [`translate_cds_to_protein_all()`](https://drostlab.github.io/orthologr/reference/translate_cds_to_protein_all.md)
  : Translate coding sequences into amino acid sequences for multiple
  files
- [`retrieve_longest_isoforms()`](https://drostlab.github.io/orthologr/reference/retrieve_longest_isoforms.md)
  : Retrieve the longest isoforms from a proteome file and save results
  as fasta file
- [`retrieve_longest_isoforms_all()`](https://drostlab.github.io/orthologr/reference/retrieve_longest_isoforms_all.md)
  : Retrieve the longest isoforms of several proteome files stored in a
  folder

## Annotation Utilities

Check annotation files, extract gene loci / splice variant features from
GFF files, and manage internal folder structure.

- [`check_annotation()`](https://drostlab.github.io/orthologr/reference/check_annotation.md)
  : Check whether an annotation file contains outlier lines
- [`extract_features()`](https://drostlab.github.io/orthologr/reference/extract_features.md)
  : Helper function to extract gene loci and splice variant IDs from GFF
  files
- [`clean_all_folders()`](https://drostlab.github.io/orthologr/reference/clean_all_folders.md)
  : Delete the internal folder hierarchy
