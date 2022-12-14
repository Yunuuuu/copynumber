
<!-- README.md is generated from README.Rmd. Please edit that file -->

# copynumber with hg38

This is an modified version of the
[copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html)
R package. It has been adjusted to support the human hg38 genome builds
(Following the tutorial in <https://github.com/aroneklund/copynumber>).
Speicifically, this is used for the `assembly` parameter in `aspcf`,
`multipcf`, `pcf`, and `winsorize` functions.

## Installation

You can install the development version of copynumber from
[GitHub](https://github.com/) with:

``` r
pak::pkg_install("Yunuuuu/copynumber")
```

Similar work also can be found in
[aroneklund](https://github.com/aroneklund/copynumber),
[igordot](https://github.com/igordot/copynumber) and
[ShixiangWang](https://github.com/ShixiangWang/copynumber). Especially,
[Irrationone](https://github.com/Irrationone/copynumber) provides a
species-agnostic approach.
