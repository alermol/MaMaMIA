# MaMaMIA

**Max Mapping for Meticulous Introgression Analysis**

MaMaMIA detects interspecific introgressions from Genotyping-by-Sequencing
(GBS) coverage. It measures breadth and read counts in exact windows, computes
a paired donor/recipient log-ratio of per-component depth, corrects GC bias,
masks unreliable windows, segments with CBS, calibrates thresholds by segment
length, and classifies calls as substitutions, additions, deletions, complex,
or other.

The method is ploidy-agnostic. A **component** is an arbitrary label for a set
of contigs (for example `A`, `B`, or `D`); it makes no genome-type or ploidy
assumption, and donor and recipient labels may overlap. The assembly registry
lists role-free **species**; the user assigns a donor and a recipient at run
time.

MaMaMIA is BAM-first. It reads BAM headers, creates CSI indexes, tiles the
reference into windows, reads GC from a `gcrank` index, measures breadth and
counts in one htslib pass, and runs the cohort pipeline. The dashboard is the
main entry point; the same pipeline is available from R.

**Version 0.1.0.** This release is not backward compatible with 0.99.x: the
configuration schema, result object, and function surface changed. See
[`NEWS.md`](NEWS.md).

## Requirements

- R 4.6 or newer and a C++23 compiler (for Rhtslib and the vendored gcrank).
- `uvr` for the pinned environment. The repository commits `uvr.toml`,
  `uvr.lock`, and `.r-version` (R 4.6.1); `uvr` installs the exact package set
  into `.uvr/library/`.
- `samtools` and `gcrank-build` only for preparing BAMs and GC indexes; the
  analysis itself does not call them.

## Install

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
fresh R session. The package vignette walks through the same workflow:
`vignette("Overview", package = "MaMaMIA")`.

## Quick start

Write a version 2 registry, for example `assemblies.yml`:

```yaml
version: 2

defaults:
  windows_per_chromosome: 1000
  allowed_windows: [100, 1000, 5000]
  mapq: 20
  exclude_flags: 1796      # unmapped, secondary, QC-fail, duplicate
  metric: breadth
  score: paired_log_ratio

server:
  max_bams: 1              # dashboard: one BAM per run
  max_total_bytes: 5368709120

species:
  triticum_aestivum:
    label: "Triticum aestivum CS T2T"
    gc_indexes: [/srv/mamamia/gcidx/cs_t2t/index.gcrank]
    contigs:
      - {id: CP155614.1, name: 1A, component: A}
      # ... 21 wheat contigs
  aegilops_speltoides:
    label: "Aegilops speltoides"
    gc_indexes: [/srv/mamamia/gcidx/ae_speltoides/index.gcrank]
    contigs:
      - {id: CH1S, name: 1S}
      # ... 7 S contigs
```

`group` is optional. When omitted, MaMaMIA uses the first number in the contig
`name` (`1A` and `1S` are both group 1). `component` is optional for a
single-component species and defaults to the species `default_component`, or
to `"default"`. A diploid species therefore needs only `id` and `name`; its
groups become 1 and 2, its component becomes `"default"`, and `copy_number`
defaults to 2.

The registry is species-only: there is no `assemblies:` section. The donor and
recipient are chosen at run time, and the assembly is derived from the two
species.

```r
library(MaMaMIA)

config <- readMaMaMIAConfig("assemblies.yml")
pair <- resolveSpeciesAssembly(
    config, donor = "aegilops_speltoides", recipient = "triticum_aestivum"
)
result <- analyzeBamInput(
    "Sq74.bam", pair$config,
    assembly = pair$key,
    windows_per_chromosome = 1000L,
    cores = 8L
)

result$samples[, c("sample", "mode", "n_windows", "n_calls")]
result$calls[result$calls$pass, ]
writeCallsToBedpe(result, "calls.bedpe")
```

`resolveSpeciesAssembly()` accepts explicit `recipient_components` and
`donor_components` for multi-component species, plus
`default_sample_composition`; a runtime `sample_composition` overrides it.
`copy_number: 0` removes a component from the recipient aggregate.

A multi-sample BAM is still one file: samples are read from its `@RG SM`
tags. `server.max_bams` limits the number of BAM files per run, and the
example registry sets it to one.

`windows_per_chromosome` must lie between the smallest and largest
`defaults.allowed_windows` value in steps of 10. The default `1000` gives about
28,000 windows for the wheat + *Ae. speltoides* reference.

The bundled demo runs the same path:

```r
demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
pair <- resolveSpeciesAssembly(
    config, donor = "secale_like", recipient = "triticum_like"
)
result <- analyzeBamInput(
    file.path(demo, "demo.bam"), pair$config,
    assembly = pair$key, windows_per_chromosome = 100L
)
```

## BAM input and registry

- Inputs must be coordinate-sorted BAM. MaMaMIA creates missing CSI indexes
  and reuses existing ones; SAM and CRAM are rejected.
- A multi-sample BAM is one file. Samples come from `@RG SM`; a BAM without
  `@RG` falls back to its file name.
- `server.max_bams` limits the number of BAM files per run. The dashboard
  accepts one BAM; a programmatic cohort run must raise the limit explicitly.
- `server.max_total_bytes` (default 5 GiB) caps dashboard uploads. An
  oversized file is rejected in the browser before transfer; the app then
  links to [bamhusk](https://github.com/alermol/bamhusk) with a command to
  slim the BAM locally.
- `gc_indexes` point to `gcrank` indexes built with `gcrank-build`. They are
  read-only and can be shared across runs.

## Pipeline

```text
preflightBams()                  headers, samples, CSI status
  -> makeAnalysisWindows()       exact-N, group-aware windows
  -> extractWindowGc()           gcrank
  -> countWindowMetrics()        counts + covered bases in one htslib pass
  -> buildWindowMetrics()        breadth, lambda_breadth, lambda_count
  -> scoreIntrogression()        paired log-ratio or donor-only
  -> buildWindowMask()           panel / cohort / blacklist / physical
  -> segmentScores()             CBS + gap merging
  -> calibrateSegments()         length-stratified empirical null
  -> classifyCalls()             class, copy state, dosage
  -> BamAnalysisResult
```

### Scoring

For a chromosome group with donor depth `lambda_d` and recipient depth
`lambda_r`:

```text
lambda = -log(1 - breadth)                    # metric = "breadth"
lambda = count * read_length / length         # metric = "count"
y      = log(lambda_d + eps) - log(lambda_r + eps)
```

Recipient homoeologs are averaged. Donor homoeologs are max-selected by
default: at each window the donor component with the strongest donor/recipient
log-ratio is used, so one strongly supported donor component is not diluted by
the others. `mean` and `sum` are available through
`defaults.donor_aggregation`. A single sample is GC-corrected with a natural
spline and centred; a cohort uses median polish. `score = "donor_only"` is the
CNV-robust fallback.

### Masks and calibration

A mask is indexed by chromosome group and relative window. Sources are tried
in order: a fixed panel mask, a cohort recurrence mask, a BED blacklist, and a
physical terminal fallback. Masked windows are set to `NA` before CBS, so
segment boundaries and gap merging still see the full grid.

Calls are calibrated by segment length. The threshold is the empirical
`1 - fdr` quantile of control segment means, raised to a per-bin floor. A
single sample is its own control; a cohort uses all samples (or `controls`).

### Classification and dosage

Calls combine donor gain, recipient loss, and a relative total-depth residual:
`substitution` (balanced replacement), `addition` (donor gain plus total gain),
`deletion` (recipient loss plus total loss), `complex` (total imbalance), or
`other`. Each call carries `copy_state`, `class_evidence`, `class_confidence`,
and `class_status`; an assigned call also reports a provisional
`recipient_copy_number`, `recipient_observed_breadth`, `affected_fraction`,
`affected_copies`, and `dosage_confidence`.

## Result contract

`analyzeBamInput()` returns a `BamAnalysisResult` with twelve fields:

| Field | Contents |
|---|---|
| `assembly`, `config` | resolved assembly and validated registry |
| `bams`, `samples` | input metadata and per-sample call counts |
| `windows`, `gc` | exact-N windows and reference GC |
| `metrics` | sample × window: count, covered, high-MAPQ fraction, breadth, GC, lambdas |
| `scores` | sample × group × window: score, side residuals and total-depth residual |
| `masks` | window-level mask with `source` and `version` |
| `calibration` | thresholds by mode and segment-length bin |
| `calls` | segments with class, copy state, total evidence, confidence, method, assignment, dosage and `pass` |
| `provenance` | versions, parameters, mask/calibration hashes |

## Dashboard

```r
runMaMaMIA(config = "assemblies.yml")
exploreIntrogressions(result)  # existing result
```

`runMaMaMIA()` is the main entry point. It preflights the BAM, matches its
header to the species registry, stages and indexes it in a private session
directory, and runs the same pipeline as `analyzeBamInput()`. The dashboard
accepts one BAM per run; a multi-sample BAM is one file and adds a **Sample**
selector as soon as its `@RG SM` tags are read.

The **Karyotype** tab shows a vertical, group-aware karyotype with score-heatmap
and call-map modes. Hover a call for class, status, placement, evidence,
confidence, and dosage. The **Calls** tab is a sortable table with configurable
columns. The **Export** tab downloads one ZIP with call tables, BEDPE, BED,
and bedGraph tracks, karyotype SVGs, per-window diagnostics, masks,
calibration, QC plots, `report.html`, provenance, session information, and MD5
checksums.

Four analysis profiles are available: `balanced` (default), `sensitive` (low
coverage), `cnv_robust` (recipient copy-number variation), and `complex`
(multi-component donors). Advanced controls expose donor aggregation,
calibration, resolution, and masking.

The bundled **Example data** button loads eight synthetic 160–260 kb chromosome
pairs covering substitutions, additions, deletions, complex and other segments,
a full substitution, and whole-chromosome gains and losses.

See [`docs-source/karyotype-tooltip.md`](docs-source/karyotype-tooltip.md) for
the tooltip field guide and
[`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) for deployment,
upload limits, and troubleshooting.

## Validation

The paired log-ratio recovers **6/6** reference events on the 24
*Aegilops speltoides* BAMs (Sq74 7S; Sq75 7S/2S/6S/5S; Sq76 7S), with one
control call ≥ 20 Mb. The donor-only fallback also recovers 6/6. On the
development host, a 28,000-window BAM takes about 19 s on one core; the 24-BAM
cohort, including global assignment, takes about 450 s with 30 cores.

Global assignment places five of the six reference events on recipient contigs
(Sq74 7D; Sq75 2B/5B/6D/7D). The Sq76 7A event is detected as an addition
without recipient loss and receives no host. Synthetic copy-neutral events are
never marked `assigned`. This validation is a programmatic cohort run and
raises `server.max_bams` explicitly; the dashboard still accepts one BAM file.

Reproduce it with the gated 24-BAM integration test:

```sh
MAMAMIA_TEST_DATA=1 MAMAMIA_SLOW_TESTS=1 \
  Rscript -e 'testthat::test_local(filter = "integration")'
```

Prior cytogenetic and FISH records for the `Sq52`–`Sq71` lines, with
sample-level caveats and reference mappings, are curated in
[`docs-source/prior-introgression-records.md`](docs-source/prior-introgression-records.md).

## Testing

```sh
Rscript -e 'testthat::test_local()'
R CMD check --no-manual
```

The suite covers htslib breadth parity (`samtools depth` and `bedtools
coverage`), the WindowMetrics/ScoreTable/MaskTable contracts, diploid through
hexaploid scoring, multi-donor max/sum, GC correction, masks, CBS boundaries,
length calibration, BEDPE export, the Shiny launcher, and the bundled example.
Real-data tests are gated because the 24 BAMs and the 1.2 GB gcrank indexes are
not shipped.

## Deployment

The production container mounts the registry read-only and passes its path in
`MAMAMIA_CONFIG`:

```sh
MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml \
  R -e 'MaMaMIA::runMaMaMIA()'
```

Uploaded BAMs are staged in a private per-session directory and removed when
the session ends. The container listens on port 3838. See
[`docs-source/shinyproxy.md`](docs-source/shinyproxy.md) for ShinyProxy,
nginx, upload limits, and troubleshooting.

## License

GPL (>= 3). See [`LICENSE.md`](LICENSE.md).

## Disclaimer

MaMaMIA is research software. Validate calls experimentally before using them
for breeding or publication decisions.

## Contact

Aleksey Ermolaev <alerm@al-erm.ru>

Issues and pull requests: <https://github.com/alermol/MaMaMIA/issues>
