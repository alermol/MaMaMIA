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

## Quick start

`runMaMaMIA()` is the main entry point: it opens a locally running Shiny
dashboard that reads the input tables, runs the whole MaMaMIA pipeline (RCA,
GC-bias correction, segmentation), and then hands over to the interactive
explorer for filtering and BEDPE export. The dashboard packages are optional
dependencies:

```r
install.packages(c("shiny", "bslib", "shinyWidgets", "DT"))

library(MaMaMIA)
runMaMaMIA(
  analysis_dir = "mamamia_analysis",
  assemblies_dir = "assemblies"
)
```

The app expects an analysis directory with one subdirectory per sample and an
assemblies directory with one subdirectory per donor reference:

```
analysis_dir/<sample>/<reference>.w<N>.dcov.bed   # donor coverage windows
analysis_dir/<sample>/<reference>.w<N>.rcov.bed   # recipient coverage windows
assemblies_dir/<reference>/donor.w<N>.gc.bed      # donor GC content
assemblies_dir/<reference>/recipient.w<N>.gc.bed  # recipient GC content
assemblies_dir/<reference>/chromosomes.tsv        # chr_id, chr_name, subgenome
```

See `?readMaMaMIAInput` for the table format and
`?defaultTargetPairs` for how donor and recipient chromosomes are paired.
Called without arguments, the app starts on the packaged `triticum` example
and the directories can be entered in the sidebar.

The BioCircos-based circos view is temporarily disabled while the BioCircos
dependency is removed: the app shows the segment table with the BEDPE export.

## BAM-first workflow

The BAM-first app starts from BAM files instead of precomputed coverage tables.
It reads the BAM headers, matches them to an assembly registry, generates
windows, extracts GC content from a `gcrank` index, counts reads, corrects GC
bias, normalizes cohorts and segments every sample:

```sh
uvr sync --frozen
uvr activate --write-shim
. .uvr/activate
R CMD INSTALL .
Rscript -e 'MaMaMIA::runMaMaMIA(config = "assemblies.yml", input = "bam")'
```

`assemblies.yml` declares role-free species (GC indexes and contigs) and,
optionally, named assemblies for scripted runs; see
`inst/config/assemblies.example.yml`. In the app the user chooses a recipient
species and a donor species, and contigs are paired by chromosome group. The
same species can be used on either side; the app has no predefined-assembly
selector. A "Swap donor and recipient" button exchanges the two roles, and
auto-detection is available as a checkbox when named assemblies are configured.
The registry is validated before any BAM work starts
(`?validateMaMaMIAConfig`). Uploaded BAMs are staged into a private session
directory and removed on session end; server-side BAMs under `server.bam_root`
are only read.

The app accepts one or more BAMs, shows a cohort summary, keeps table
selections per sample, and exports a per-sample BEDPE or a combined BEDPE with
a `sample` column. Resolution (`windows_per_chromosome`) and normalization
(`none` or `median`) are selectable; the GC fit is capped by `max_fit_windows`
(50,000 by default) for high-resolution runs.

## Deployment

The reproducible image and ShinyProxy configuration live in `deployment/`.
Build from the repository root with
`docker build -f deployment/Dockerfile -t mamamia:0.99.0 .` and see
[`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) for the data
prerequisites, local smoke test, read-only mounts, upload/timeouts and resource
limits.

## Scripted workflow

The same pipeline can be driven from a script; `exploreIntrogressions()` opens
the explorer for an existing result:

```r
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
isa <- segments(rca,
  target_pairs = c("OY997261.1" = "CP155614.1"),
  verbose = FALSE
)
exploreIntrogressions(isa)
writeToBedpe(isa, "introgressions.bedpe")
```

For the BAM workflow, resolve the two species into an assembly first:

```r
config <- readMaMaMIAConfig("assemblies.yml")
pair <- resolveSpeciesAssembly(
  config, donor = "secale_cereale", recipient = "triticum_aestivum"
)
result <- analyzeBamInput(
  c("sample1.bam", "sample2.bam"), pair$config,
  assembly = pair$key, windows_per_chromosome = 1000L
)
writeCohortToBedpe(result, "cohort.bedpe")
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
