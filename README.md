# MaMaMIA

`MaMaMIA` (Max Mapping for Meticulous Introgression Analysis) detects interspecific introgressions from Genotyping-by-Sequencing (GBS) coverage data.

It builds a Read Count Array (`RCA`) from windowed donor and recipient coverage and GC-content tables, corrects GC bias with a zero-inflated negative binomial model, and segments pairwise coverage differences with circular binary segmentation (CBS) into an Introgression Segment Array (`ISA`). Putative introgressions can be plotted, filtered, and exported as BEDPE.

`MaMaMIA` ships with an example wheat dataset (`triticum`: *T. timopheevii* donor vs *T. aestivum* Chinese Spring T2T recipient). See `vignette("Overview", package = "MaMaMIA")` for a full walkthrough.

## Installation

You can install `MaMaMIA` from [GitHub](https://github.com/alermol/mamamia) with:

```r
# install.packages("devtools")
devtools::install_github("alermol/mamamia")
```

## Resource use

The GC-bias fit differentiates the likelihood window by window, so **peak memory
and runtime scale with the number of fitted windows, not with `cores`**: budget
roughly 0.1 MB per window and 0.85 s per 1,000 windows at `cores = 4`. The
bundled `triticum` example fits 34,000 windows and peaks near 3.6 GB, so fitting
is practical up to a few tens of thousands of windows on a typical machine
(350,000 would need roughly 36 GB).

For a batch of samples, run **one process per sample with `cores = 1`** rather
than one wide multi-threaded call, with `OMP_NUM_THREADS=1` in each worker. The
thread count shifts fitted values by around `1e-7`, so fixing it makes a sample's
result independent of the worker it ran on, while throughput comes from the
number of workers. The overview vignette has the full guidance and figures.

## Disclaimer

This package is still under active development, the content is therefore subject to change. 

## Contact

Suggestions and bug reports: please [open an issue](https://github.com/alermol/mamamia/issues).
