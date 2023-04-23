
<!-- README.md is generated from README.Rmd. Please edit that file -->

# copynumber with hg38

This is an modified version of the
[copynumber](http://bioconductor.org/packages/release/bioc/html/copynumber.html)
R package. It has been adjusted to support the human hg38 genome builds
(Following the tutorial in <https://github.com/aroneklund/copynumber>),
in addition, it will be okay to provide a customized cytoband data.frame
directly.

## Assembly usage

Speicifically, this is used for the “assembly” parameter in `aspcf`,
`multipcf`, `pcf`, and `winsorize` functions. We can provide a scalar
string or a data.frame directly.

- string Only “hg16”, “hg17”, “hg18”, “hg19”, “hg38”, “mm7”, “mm8” and
  “mm9” are supported.
- a data.frame with 5 columns (names don’t matter) in the order:
  1.  chromosome
  2.  start_pos
  3.  end_pos
  4.  band
  5.  gieStain

Notes: it’s easy to do this with:

``` r
# assembly should be a genome string like "hg38", "hg19", or "mm9", "mm10".
data.table::fread(
    sprintf("http://hgdownload.cse.ucsc.edu/goldenpath/%s/database/cytoBand.txt.gz",
    assembly),
    header = FALSE
)
```

## Installation

You can install the development version of copynumber from
[GitHub](https://github.com/) with:

``` r
pak::pkg_install("Yunuuuu/copynumber")
```

Similar work can also be found in
[aroneklund](https://github.com/aroneklund/copynumber),
[igordot](https://github.com/igordot/copynumber),
[ShixiangWang](https://github.com/ShixiangWang/copynumber). And
[Irrationone](https://github.com/Irrationone/copynumber) provides a
species-agnostic approach.

## Official copynumber

Since the official copynumber has been labelled as deprecated and will
be removed from Bioconductor version 3.18. This repo also kept a backup
of the official copynumber in `official` branch (You can install it with
`pak::pkg_install(Yunuuuu/copynumber@official)`).
