# Locate the DIAMOND executable used by orthologr

Returns the directory holding the `diamond` executable that `orthologr`
will use. A `diamond` on the `PATH` takes precedence (returned as
`NULL`, i.e. no `path` needed); otherwise the copy installed with
[`diamond_install`](https://drostlab.github.io/orthologr/reference/diamond_install.md)
is returned.

## Usage

``` r
diamond_path()
```

## Value

`NULL` if `diamond` is on the `PATH`, the cache directory if only the
cached copy exists, and an error otherwise.

## See also

[`diamond_install`](https://drostlab.github.io/orthologr/reference/diamond_install.md)
