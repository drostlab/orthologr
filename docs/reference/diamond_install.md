# Install the DIAMOND aligner into the user cache

Downloads a pre-compiled DIAMOND2 binary from the official GitHub
releases (<https://github.com/bbuchfink/diamond/releases>) that matches
the operating system of the user and stores it in a user-level cache
directory (see [`R_user_dir`](https://rdrr.io/r/tools/userdir.html)). No
administrator rights are required. Once installed,
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
and
[`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
find the cached binary automatically whenever `path = NULL` and no
`diamond` is available on the `PATH`.

## Usage

``` r
diamond_install(version = "latest", force = FALSE, quiet = FALSE)
```

## Arguments

- version:

  a DIAMOND release tag, e.g. `"v2.1.9"`. The default `"latest"` uses
  the most recent release.

- force:

  logical; re-download even if DIAMOND is already cached.

- quiet:

  logical; suppress progress messages.

## Value

Invisibly, the directory containing the `diamond` executable. This can
be passed to the `path` argument of the DIAMOND functions.

## Details

Pre-compiled binaries are provided upstream for Linux (x86_64), macOS
and Windows (x86_64). On other platforms (e.g. Linux ARM) please install
DIAMOND via conda: `conda install -c bioconda diamond`.

## See also

[`diamond_path`](https://drostlab.github.io/orthologr/reference/diamond_path.md),
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md)

## Examples

``` r
if (FALSE) { # \dontrun{
diamond_install()
diamond_path()
} # }
```
