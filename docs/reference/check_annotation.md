# Check whether an annotation file contains outlier lines

Some annotation files include lines with character lengths greater than
65000. This causes problems when trying to import such annotation files
into R using
[`import`](https://rdrr.io/bioc/rtracklayer/man/export.html). To
overcome this issue, this function screens for such lines in a given
annotation file and removes these lines so that
[`import`](https://rdrr.io/bioc/rtracklayer/man/export.html) can handle
the file.

## Usage

``` r
check_annotation(annotation_file, remove_annotation_outliers = FALSE)
```

## Arguments

- annotation_file:

  a file path tp the annotation file.

- remove_annotation_outliers:

  shall outlier lines be removed from the input `annotation_file`? If
  yes, then the initial `annotation_file` will be overwritten and the
  removed outlier lines will be stored at
  [`tempdir`](https://rdrr.io/r/base/tempfile.html) for further
  exploration.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# download an example annotation file from NCBI RefSeq
Ath_path <- biomartr::getGFF(organism = "Arabidopsis thaliana")
# run annotation file check on the downloaded file
orthologr::check_annotation(Ath_path)
# several outlier lines were detected, thus we re-run the
# function using 'remove_annotation_outliers = TRUE'
# to remove the outliers and overwrite the file
check_annotation(Ath_path, remove_annotation_outliers = TRUE)
} # }
```
