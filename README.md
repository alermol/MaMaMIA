# MaMaMIA

`MaMaMIA` (Max Mapping for Meticulous Introgression Analysis) detects
interspecific introgressions from Genotyping-by-Sequencing (GBS) coverage data.

It builds a Read Count Array (`RCA`) from windowed donor and recipient coverage
and GC-content, corrects GC bias with a zero-inflated negative binomial model,
and segments pairwise coverage differences with circular binary segmentation
(CBS) into an Introgression Segment Array (`ISA`). Putative introgressions can
be plotted, filtered, and exported as BEDPE.

The **BAM-first workflow** is the only dashboard workflow: `MaMaMIA` reads BAM
headers, creates CSI indexes, tiles the reference into windows, reads GC content
from a `gcrank` index, counts reads, and runs the cohort pipeline. The assembly
registry lists role-free **species**; the user assigns one species as the donor
and another as the recipient at run time. The previous directory-layout launcher
and its readers were removed; see
[Migration from the directory layout](#migration-from-the-directory-layout) if
you have precomputed tables.

## Requirements

- R 4.6 or newer and a C++23 compiler (for Rhtslib and the C++ sources).
- `uvr` for the reproducible, pinned environment. The repository commits
  `uvr.toml`, `uvr.lock` and `.r-version` (R 4.6.1); `uvr` installs the exact
  package set into the project library `.uvr/library/`.
- `samtools` and `gcrank-build` only for preparing BAMs and GC indexes; the
  analysis itself does not call them.

## Reproducible environment with uvr

The project migrated from `renv` to `uvr`. There is no `renv.lock` and no
`renv/` directory: dependencies are declared in `uvr.toml`, resolved in
`uvr.lock`, and installed into `.uvr/library/`.

```sh
git clone https://github.com/alermol/MaMaMIA.git
cd MaMaMIA
uvr sync --frozen          # install R 4.6.1 and every locked package
uvr activate --write-shim  # (re)create .uvr/activate
. .uvr/activate            # put the uvr-managed R on PATH
R CMD INSTALL .
uvr doctor                 # diagnose a broken environment
```

`uvr sync --frozen` fails instead of silently updating when `uvr.toml` and
`uvr.lock` disagree. `.Rprofile` activates `.uvr/library/` automatically in a
fresh R session. For a one-off script, `uvr run path/to/script.R -- args` uses
the same environment without activating a shell.

If R 4.6 and the Bioconductor dependencies are already installed, the package
can also be installed with `remotes::install_github("alermol/mamamia")`.

## BAM-first quick start

Write an assembly registry, for example `assemblies.yml`. Species are
role-free and declare only their label, GC indexes and contigs:

```yaml
version: 1

defaults:
  windows_per_chromosome: 1000
  allowed_windows: [100, 1000, 10000]
  mapq: 20
  exclude_flags: 2308      # unmapped + secondary + supplementary
  max_fit_windows: 50000
  normalization_by_side: true

server:
  max_bams: 20
  max_total_bytes: 5368709120   # 5 GiB

species:
  triticum_aestivum_cs_t2t:
    label: "Triticum aestivum CS T2T"
    default_subgenome: A
    gc_indexes: [/srv/mamamia/gcidx/triticum-aestivum/index.gcrank]
    contigs:
      - {id: CP155614.1, name: 1A, subgenome: A}
      - {id: CP155611.1, name: 2A, subgenome: A}
      # ... one entry per target contig
  secale_cereale:
    label: "Secale cereale (rye)"
    default_subgenome: R
    gc_indexes: [/srv/mamamia/gcidx/secale-cereale/index.gcrank]
    contigs:
      - {id: OZ284382.1, name: 1R, subgenome: R}
      - {id: OZ284383.1, name: 2R, subgenome: R}
      # ... one entry per target contig
```

`gc_indexes` points at a `*.gcrank` directory with `contigs.tsv` and
`manifest.json` (format version `0.1.0`). Contig lengths are read from the
index. Relative paths are resolved against the directory containing the YAML
file. The complete example is `inst/config/assemblies.example.yml`; validate a
production registry before using it:

```sh
Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml
```

Start the dashboard:

```r
library(MaMaMIA)
runMaMaMIA(config = "assemblies.yml")
```

With `config`, `runMaMaMIA()` opens the BAM-first app. The sidebar starts
with one BAM upload control and the **Example data** button. Samples are read
from the `@RG SM` tags and read groups from the `@RG ID` tags in the BAM header;
a single BAM can therefore contain one or several samples. There is no
multi-file upload and no manual single/multi mode switch: the interface detects
the sample count from the header and enables cohort controls (reference sample,
combined export) automatically when more than one sample is present. The
analysis controls are revealed only after preflight succeeds, as collapsible
panels (only one panel is open at a time): **Analysis** and **Assembly**. The
**Run analysis** button sits below both panels and stays disabled until
preflight and the selected resolution are ready.

Choose a **recipient species** and a **donor species**; each species can be
used on either side of different runs (one run needs two different species), and
the **Swap donor and recipient** button reverses the comparison without
reselecting anything. Subgenome selectors appear only for species with several
subgenomes. The **Assembly** panel contains only these selectors and the swap
button. There is no predefined-assembly selector and no BAM-header
auto-detection checkbox: contigs are paired by chromosome group (the optional
`group` field, otherwise the first number in `name`). A registry without a
`species:` section falls back to header-based auto-detection automatically.
BAM processing uses the MAPQ, flag mask, fit cap and core count from the
validated registry; those controls are not exposed in the interface.

Upload one BAM, choose a resolution and, for several samples, a reference
sample or median normalization. BAMs are accepted only through the upload
control; the container does not mount a BAM directory. A shared filter bar above
the results controls chromosomes, chromosome pairs and the putative-only view
for every tab. For multi-sample BAMs the same bar carries a **Sample** selector
that switches the circos, coverage plot, segment table and exports together.
The main area is split into **Introgression circos**,
**Coverage / dosage** (mirror or difference), **Segments** and **Export** tabs.
The circos is zoomable (scroll to zoom, drag to pan), labels chromosomes with
their registry names, and hides hover tooltips when the pointer leaves the
plot. Call threshold, minimum segment width and row selection sit in
**Segments**; cohort summaries and the per-sample and combined BEDPE downloads
live in **Export**. Preflight status, detected assemblies, problems and every
pipeline stage are written to the calculation log in the footer. The line above
the log reports only the number of BAM files. The app keeps table selections
per sample and exports a per-sample BEDPE or a combined BEDPE with a `sample`
column. Uploaded BAMs are staged into a private session directory and removed on
session end.

A small synthetic example ships as pre-computed tables and is loaded with the
**Example data** button. It has two species and a variable-GC index, needs no
BAM upload and no mounted BAM, fills the example species/resolution fields, and
replays the normal pipeline stage messages before showing the stored result. The
button is available in both modes; pressing it again clears the example and
restores the registry resolution choices. Set `windows_per_chromosome` to 100
or 1,000 for the example (the bundled demo contigs are too short for a useful
10,000-window GC fit).

### Scripted BAM-first run

The small `inst/extdata/bam-demo` fixture remains available for scripted runs
and documentation examples:

```r
demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
pair <- resolveSpeciesAssembly(
  config, donor = "secale_like", recipient = "triticum_like"
)

result <- analyzeBamInput(
  file.path(demo, "demo.bam"), pair$config,
  assembly = pair$key,
  windows_per_chromosome = 100L,
  normalization = "none"
)
result$samples
writeCohortToBedpe(result, "demo.bedpe")
```

`pair$config` contains the derived assembly; `pair$unpaired` explains any
excluded contigs. Swap the two species in `resolveSpeciesAssembly()` to run the
reverse comparison. The same two calls accept two or more BAMs and
`normalization = "median"` or `"reference"`.

## Configuration reference

| Key | Default | Meaning |
|---|---|---|
| `version` | required | Registry format, currently `1`. |
| `defaults.windows_per_chromosome` | `1000` | Windows per target contig. |
| `defaults.allowed_windows` | `[100, 1000, 10000]` | Resolutions offered by the app and accepted by validation. |
| `defaults.mapq` | `20` | Minimum mapping quality. |
| `defaults.exclude_flags` | `2308` | SAM flag mask; a read is excluded when **any** bit is set. |
| `defaults.max_fit_windows` | `50000` | Maximum ideal windows used for one sample's GC fit. |
| `defaults.gc_mmap` | `true` | Memory-map GC indexes. |
| `defaults.gc_validate_crc` | `false` | Verify every index payload against the manifest at open. |
| `defaults.normalization_by_side` | `true` | Normalize donor and recipient windows independently. |
| `server.max_bams` | `20` | Maximum BAMs per run. |
| `server.max_total_bytes` | 5 GiB | Maximum total upload size. |
| `species.<key>.label` | required | Display name. |
| `species.<key>.default_subgenome` | none | Preselected subgenome for a polyploid. |
| `species.<key>.gc_indexes` | required | One or more `*.gcrank` directories. |
| `species.<key>.contigs` | required | `id`, `name`, `subgenome`, optional `group`. |
| `assemblies.<key>` | optional | Named preset for scripted runs; the app does not expose a preset list. |

A species-only registry needs at least two species. An optional named assembly
may reference species (`recipient: wheat`, `donor: rye`) or declare inline
contigs. The species builder pairs contigs by chromosome group and only
includes the paired contigs; `buildSampleTables()` falls back to
`defaultTargetPairs()` when an assembly declares no `pairs`.

### Registry growth and long-lived containers

Mount the whole GC-index directory, not individual indexes. A new
`<species>/index.gcrank` becomes readable inside running containers immediately
because indexes are opened on demand and memory-mapped. The YAML registry is
different: the app reads and validates it once when the R process starts.
Running ShinyProxy containers keep their loaded registry; close the app and
open it again (or restart the container) after changing `assemblies.yml`. The
sidebar then shows the new species in both selectors. See
[`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) for the atomic-replace
procedure.

## Resolution, normalization and the fit cap

**Resolution.** `windows_per_chromosome` tiles every target contig into exactly
that many contiguous half-open windows. Donor and recipient windows are matched
by position index, not by genomic coordinates. Every contig must be at least as
long as the resolution, otherwise the app disables Run and names the contig.
1000 is a good default; 100 makes quick exploratory runs, 10,000 is the
high-resolution mode.

**Normalization.** With one BAM, `normalization = "none"` compares
GC-corrected coverage (`cor.gc`) between the paired contigs. With a cohort:

- `"median"` (default): `dose = (C + eps) / (median(C) + eps)` across samples at
  each window;
- `"reference"`: `dose = (C + eps) / (C_ref + eps)` using `reference_sample`;
- `eps` is `pseudocount` (default `0.5`), so a zero side stays finite and an
  all-zero window gets dose `1`.

By default the donor and recipient sides are normalized independently
(`normalization_by_side = true`), because they come from different genomes and
mappability regimes. CBS then segments `dose.don - dose.rec`; without
normalization it segments `cor.gc.don - cor.gc.rec`. Normalization reads
`cor.gc` and never overwrites it.

**Fit cap.** The GC-bias fit is the memory bottleneck. `max_fit_windows` caps
the number of ideal windows fitted per sample. Larger sets are subsampled
deterministically by subgenome and GC decile (at least ten windows per non-empty
stratum); prediction still covers every window. A cohort whose total window
count exceeds the cap also prints a performance warning. `Inf` disables the
cap.

## Resource requirements and high-resolution runs

The GC-bias fit differentiates the likelihood window by window, so **peak
memory and runtime scale with the number of fitted windows, not with `cores`**.
Budget roughly 0.1 MB per fitted window and 0.85 s per 1,000 windows at
`cores = 4`. For example, 34,000 fitted windows peak near 3.6 GB; 50,000
fitted windows typically stay within 5-6 GB, and a high-resolution cohort with
the fit capped at 50,000 peaks around 8-12 GB including R and index mappings.

A 10,000-window run over 28 contigs produces 280,000 analysis windows per
sample. Keep `max_fit_windows` at 50,000 or lower, expect tens of minutes per
sample, and provision at least 16 GB. The ShinyProxy spec in `deployment/`
uses `container-memory-limit: 16g` and `container-cpu-limit: 4` for this
reason.

For a batch of samples, run **one process per sample with `cores = 1`** rather
than one wide multi-threaded call, with `OMP_NUM_THREADS=1` in each worker. The
thread count shifts fitted values by around `1e-7`, so fixing it makes a
sample's result independent of the worker it ran on, while throughput comes
from the number of workers. The overview vignette has the full guidance.

## Counting semantics and verification

Reads are counted with the htslib backend at their 5' end: a forward-strand
read belongs to the window containing its leftmost coordinate, a reverse-strand
read to the window containing the alignment end. This is **start-based**, not
overlap-based. `bedtools coverage -counts`, used by the legacy workflow, counts
a read in every window it overlaps, so the two can differ by the few reads
whose alignment spans a window boundary while their 5' end is outside the
window. The difference is confined to those boundary-crossing reads and is
expected, not a counting bug. Both `countReads()` and `analyzeBamInput()` use
the package's htslib counter; there is no bamsignals dependency. The pipeline
counts every BAM, one sample or several, with the same flag and boundary rules.

Both MAPQ and flag filtering follow `samtools`:

```sh
# excluded when any of the 2308 bits is set, MAPQ >= 20
samtools view -c -F 2308 -q 20 sample.bam d1:1-200
```

The command counts reads overlapping `d1:1-200`; to reproduce a window count,
keep reads whose 5' end is inside the window (forward: `POS`; reverse: alignment
end). GC values from `gcrank` match `bedtools nuc` at the six-decimal precision
stored in the legacy BED files. `?countReads` and `?gcrankGc` document the
details.

## Core analysis API

`analyzeBamInput()` returns a `BamAnalysisResult` with the per-sample `RCA` and
`ISA` objects. They can be inspected, subset, reversed, plotted and exported
exactly like the objects built by hand:

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

`RCA()`, `correctReadCounts()`, `segments()`, `writeToBedpe()` and
`exploreIntrogressions()` are the supported table-level API; they are unchanged
by the removal of the directory-layout launcher.

## Migration from the directory layout

The directory-layout launcher and its readers (`readMaMaMIAInput()`,
`listMaMaMIASamples()`, `listMaMaMIAReferences()`, `launcher_app()` and the
`analysis_dir` / `assemblies_dir` / `input` arguments of `runMaMaMIA()`) were
removed. `runMaMaMIA(config = ...)` is the single Shiny entrypoint. To
migrate:

1. Keep the BAMs and the `gcrank` indexes; declare the two genomes as `species:`
   in a registry (see `inst/config/assemblies.example.yml`).
2. Remove the precomputed `analysis_dir` / `assemblies_dir` tables; the app now
   tiles, extracts GC and counts reads itself.
3. Run `analyzeBamInput()` in scripts, or open the app with
   `runMaMaMIA(config = "assemblies.yml")`.

If you already have windowed tables and want to keep them, call the core API
directly: `RCA()`, `correctReadCounts()` and `segments()` accept the same
four-column tables as before, and `buildSampleTables()` can assemble them into
the pipeline contract. `defaultTargetPairs()` remains available for
chromosome-number pairing when an assembly declares no `pairs`.

## Deployment

The reproducible image and ShinyProxy configuration live in `deployment/`.
Build from the repository root with
`docker build -f deployment/Dockerfile -t mamamia:0.99.1 .` and see
[`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) for the data
prerequisites, the local smoke test, read-only mounts, upload/timeouts and
resource limits. The container entrypoint is
`uvr run /app/deployment/start.R`; the registry path comes from
`MAMAMIA_CONFIG`.

## Testing

```sh
uvr sync --frozen
. .uvr/activate
Rscript -e 'testthat::test_local(reporter = "summary")'
R CMD build .
R CMD check --no-manual MaMaMIA_*.tar.gz
```

Real-data tests are gated: set `MAMAMIA_TEST_DATA=1` (and optionally
`MAMAMIA_TEST_DATA_DIR`) and `MAMAMIA_SLOW_TESTS=1` for the golden regression.
The complete validation gate is described in `AGENTIC_PLAN.md` section 6.3.

## Disclaimer

This package is still under active development, the content is therefore subject to change.

## Contact

Suggestions and bug reports: please [open an issue](https://github.com/alermol/mamamia/issues).
