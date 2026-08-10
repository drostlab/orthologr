# Interface function to Orthofinder2

Run Orthofinder2 from R.

## Usage

``` r
orthofinder2(
  proteome_folder,
  use_existing_output = FALSE,
  import_type = NULL,
  comp_cores = 1,
  orthogroup_only = TRUE,
  of_path = NULL
)
```

## Arguments

- proteome_folder:

  file path to a folder storing the proteome sequences of the species
  for which orthology inference shall be performed.

- use_existing_output:

  a logical value indicating whether or not an existing `Orthofinder2`
  output folder shall be used fo further import and processing. If
  `use_existing_output = TRUE` is selected then please specify the file
  path to the in `proteome_folder` where the `Orthofinder2` output
  folder can be found.

- import_type:

  type of `Orthofinder2` output that shall be imported after running
  `Orthofinder2`. Options are:

  - `import_type = "orthogroups_core"`

  - `import_type = "orthogroups_pairwise"`

  - `import_type = ""`

  - `import_type = ""`

- comp_cores:

  number of cores that shall be used for parallel processing. Default is
  `cores = 1`.

- orthogroup_only:

  logical, default `TRUE`. When `TRUE`, OrthoFinder is run with the
  `-og` flag, which stops the analysis after orthogroup inference and
  skips the more computationally expensive gene-tree and orthologue
  inference steps. Set to `FALSE` to run the full OrthoFinder pipeline.

- of_path:

  a character string specifying the path to the locally installed
  `orthofinder` executable. A possible specification could be
  `of_path = "/opt/miniconda3/bin/"` which internally will translate to
  `/opt/miniconda3/bin/orthofinder`. The default is `of_path = NULL`
  which means that orthofinder assumes users have their `orthofinder`
  executable stored at `/opt/miniconda3/bin/orthofinder`.

## Note

This function assumes that users have `OrthoFinder` installed via
`miniconda` and stored at `~/opt/miniconda3/bin/`. In addition, DIAMOND
needs to be installed as executable tool (/usr/local/bin).

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# specify species names
orgs <- c("Arabidopsis lyrata",
          "Capsella rubella", "Solanum lycopersicum")
# download proteome files for all species
biomartr::getProteomeSet(db = "refseq", organisms = orgs, path = "of_proteomes")
# download annotation files for all species
biomartr::getGFFSet(db = "refseq", organisms = orgs, path = "of_gff")
# select longest splice variant per gene locus
retrieve_longest_isoforms_all(proteome_folder = "of_proteomes",
                              annotation_folder = "of_gff",
                              annotation_format = "gff",
                              output_folder = "of_proteomes_longest_sv")
# run orthofinder2 to infer ortho groups for the specified species
orthofinder2(proteome_folder = "of_proteomes_longest_sv", comp_cores = 4)
} # }
```
