<p align="center">
  <img src="inst/www/mamamia-logo.svg" alt="MaMaMIA logo" width="96" height="96">
</p>

<h1 align="center">MaMaMIA</h1>

<p align="center">
  <strong>Max Mapping for Meticulous Introgression Analysis</strong><br>
  Version 0.1.0 (first public release) · GPL (>= 3) · R 4.6+
</p>

MaMaMIA detects interspecific introgressions from Genotyping-by-Sequencing
(GBS) coverage. It measures breadth and read counts in exact windows, compares
donor and recipient depth with a paired log-ratio across any number of genome
components, corrects GC bias, masks unreliable windows, segments with CBS, and
calibrates calls by segment length. The method is ploidy-agnostic and BAM-first:
the same pipeline powers the Shiny dashboard and the R API.

## Highlights

- Calls substitutions, additions, deletions, complex, and other events.
- Supports diploid through hexaploid designs, multi-donor samples, and
  overlapping component labels.
- Single-sample or cohort analysis with GC correction, recurrence/blacklist
  masks, and length-stratified empirical calibration.
- Interactive karyotype, sortable call table, and single-archive export bundle.
- Reproducible R 4.6.1 environment via `uvr`.
- Paired and donor-only scoring both recover 6/6 reference events on the
  24-BAM *Aegilops speltoides* validation set.

## Requirements

- R 4.6 or newer and a C++23 compiler (for Rhtslib and the vendored gcrank).
- `uvr` for the pinned package set.
- `samtools` and `gcrank-build` only for preparing BAMs and GC indexes.

## Install

```sh
git clone https://github.com/alermol/MaMaMIA.git
cd MaMaMIA
uvr sync --frozen
uvr activate --write-shim
. .uvr/activate
R CMD INSTALL .
```

## Quick start

Create a version 2 assembly registry. Start from
[`inst/config/assemblies.example.yml`](inst/config/assemblies.example.yml), or
locate the installed copy with
`system.file("config/assemblies.example.yml", package = "MaMaMIA")`.

Input BAMs must be coordinate-sorted; MaMaMIA creates missing CSI indexes.
Then analyze a BAM:

```r
library(MaMaMIA)

config <- readMaMaMIAConfig("assemblies.yml")
pair <- resolveSpeciesAssembly(
    config, donor = "aegilops_speltoides", recipient = "triticum_aestivum"
)
result <- analyzeBamInput(
    "sample.bam", pair$config,
    assembly = pair$key,
    windows_per_chromosome = 1000L,
    cores = 8L
)

writeCallsToBedpe(result, "calls.bedpe")
```

Launch the dashboard with:

```r
runMaMaMIA(config = "assemblies.yml")
```

The dashboard accepts one BAM per run; a multi-sample BAM is one file read
through its `@RG SM` tags. Oversized BAMs can be slimmed locally with
[bamhusk](https://github.com/alermol/bamhusk) before upload.

## Documentation

Detailed configuration, method, result-field, and deployment guides will be
maintained in the [MaMaMIA wiki](https://github.com/alermol/MaMaMIA/wiki).

- `vignette("Overview", package = "MaMaMIA")` — end-to-end workflow.
- `?analyzeBamInput`, `?runMaMaMIA` — function reference.
- [`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) — deployment notes.

## License

GPL (>= 3). See [`LICENSE.md`](LICENSE.md).

## Contact

Aleksey Ermolaev <alerm@al-erm.ru> ·
[Issues](https://github.com/alermol/MaMaMIA/issues)
