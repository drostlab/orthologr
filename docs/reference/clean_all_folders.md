# Delete the internal folder hierarchy

This function deletes all internal folders that have been created during
pipeline processing. Internally this function uses
[`unlink`](https://rdrr.io/r/base/unlink.html) to delete all folders
created by the pipline.

## Usage

``` r
clean_all_folders(foldernames)
```

## Arguments

- foldernames:

  a character vector storing the folder names that shall be deleted.

## Value

This is a void function.

## Details

This function takes a vector storing the names of the folders that shall
be deleted, e.g.:
`clean_all_folders( c("_alignment", "_blast_db", "_calculation") )`.

Since orthologr is a package for pipeline processing based on interface
functions, many partial results need to be written to a hard drive to
allow subsequent programs to access the output of previous computations.
The R core conventions do not favour this behaviour of functions. This
would indicate that the `clean_folders` argument had to be `TRUE` by
default in all functions. Nevertheless, this would distract some
pipeline functions and therefore `clean_folders` = `FALSE` by default
leaving you with a folder environment that needs to be cleaned by hand,
meaning that when using a specific function or pipeline function you
must specify `clean_folders` = `TRUE` in case you want to remove all
internal folders after pipeline processing.

## See also

[`divergence_stratigraphy`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md),
[`blast`](https://drostlab.github.io/orthologr/reference/blast.md),
[`blast_best`](https://drostlab.github.io/orthologr/reference/blast_best.md),
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)

## Author

Hajk-Georg Drost
