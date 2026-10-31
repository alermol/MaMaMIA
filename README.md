# MaMaMIA

`MaMaMIA` (Max Mapping for Meticulous Introgression Analysis) detects
interspecific introgressions from Genotyping-by-Sequencing (GBS) coverage data.

Version 0.1.0 uses **breadth of coverage and read counts** in exact windows. It
computes a paired donor/recipient log-ratio of per-component depth, corrects GC
bias with a spline (single sample) or median polish (cohort), applies a
window-level recurrence/blacklist mask, segments with CBS, calibrates
thresholds by segment length, and classifies calls as substitutions, additions,
deletions or other. The method is ploidy-agnostic: a donor and a recipient may each
contribute any number of components. A **component** is an arbitrary label for
a set of contigs (for example `A`, `B` or `D`); it carries no genome-type or
ploidy assumption, and donor and recipient labels may overlap.

The BAM-first workflow is the only workflow. `MaMaMIA` reads BAM headers,
creates CSI indexes, tiles the reference into windows, reads GC from a `gcrank`
index, measures breadth and counts in one htslib pass, and runs the cohort
pipeline. The assembly registry lists role-free **species**; the user assigns a
donor and a recipient at run time.

**Upgrading from 0.99.1?** This version changes the configuration schema, the
result object and the function surface. See [`MIGRATION.md`](MIGRATION.md).

## Requirements

- R 4.6 or newer and a C++23 compiler (for Rhtslib and the vendored gcrank).
- `uvr` for the reproducible, pinned environment. The repository commits
  `uvr.toml`, `uvr.lock` and `.r-version` (R 4.6.1); `uvr` installs the exact
  package set into the project library `.uvr/library/`.
- `samtools` and `gcrank-build` only for preparing BAMs and GC indexes; the
  analysis itself does not call them.

## Reproducible environment with uvr

Dependencies are declared in `uvr.toml`, resolved in `uvr.lock`, and installed
into `.uvr/library/`.

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
fresh R session.

## Quick start

Write a version 2 registry, for example `assemblies.yml`:

```yaml
version: 2

defaults:
  windows_per_chromosome: 1000
  allowed_windows: [100, 1000, 5000]
  mapq: 20
  exclude_flags: 1796      # unmapped + secondary + QC-fail + duplicate
  metric: breadth
  score: paired_log_ratio
  donor_aggregation: max   # default; mean | max | sum, Advanced panel

server:
  max_bams: 30
  max_total_bytes: 5368709120

species:
  triticum_aestivum:
    label: "Triticum aestivum CS T2T"
    gc_indexes: [/srv/mamamia/gcidx/triticum-aestivum/index.gcrank]
    contigs:
      - {id: CP155614.1, name: 1A, component: A}
      # ... 21 wheat contigs
  aegilops_speltoides:
    label: "Aegilops speltoides"
    gc_indexes: [/srv/mamamia/gcidx/aegilops-speltoides/index.gcrank]
    contigs:
      - {id: CH1S, name: 1S}
      # ... 7 S contigs
```

`group` is optional. When omitted, MaMaMIA uses the first number in the contig
`name` (`1A` and `1S` are both group 1). Provide an explicit `group` only when
the names do not encode it.

`component` is also optional for a single-component species. It defaults to
the species `default_component`, or to `"default"` when no default is given.
A truly diploid species therefore needs only `id` and `name`:

```yaml
contigs:
  - {id: chr1, name: 1}
  - {id: chr2, name: 2}
```

`group` becomes 1 and 2, `component` becomes `"default"`, and `copy_number`
defaults to 2.

The registry is species-only: there is no `assemblies:` section. The donor and
recipient are chosen at run time, and the assembly is derived from the two
species. Run the pipeline:

```r
library(MaMaMIA)

config <- readMaMaMIAConfig("assemblies.yml")
pair <- resolveSpeciesAssembly(
    config, donor = "aegilops_speltoides", recipient = "triticum_aestivum"
)
result <- analyzeBamInput(
    c("Sq74.bam", "Sq75.bam", "Sq76.bam"),
    pair$config,
    assembly = pair$key,
    windows_per_chromosome = 1000L,
    cores = 8L
)

result$samples[, c("sample", "mode", "n_windows", "n_calls")]
result$calls[result$calls$pass, ]
writeCallsToBedpe(result, "calls.bedpe")
```

`windows_per_chromosome` must lie between the smallest and largest
`defaults.allowed_windows` value in steps of 10. The allowed values are the
slider tick marks (100, 1,000 and 5,000 by default). The default `1000` gives
about 28,000 windows for the wheat + *Ae. speltoides* reference.

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

## Pipeline

```
preflightBams()                  headers, samples, CSI status
  -> makeAnalysisWindows()       exact-N, group-aware windows
  -> extractWindowGc()           gcrank
  -> countWindowMetrics()        counts + covered bases in one htslib pass
  -> buildWindowMetrics()        breadth, lambda_breadth, lambda_count
  -> scoreIntrogression()        paired log-ratio or donor-only
  -> buildWindowMask()           panel / cohort / blacklist / physical
  -> segmentScores()             CBS + gap merging
  -> calibrateSegments()         length-stratified empirical null
  -> classifyCalls()             substitution / addition / deletion / other
  -> BamAnalysisResult
```

### Score

For a chromosome group with donor depth `lambda_d` and recipient depth
`lambda_r`:

```
lambda = -log(1 - breadth)                    # metric = "breadth"
lambda = count * read_length / length         # metric = "count"
y      = log(lambda_d + eps) - log(lambda_r + eps)
```

Recipient homoeologs are averaged. Donor homoeologs are max-selected by
default: at each window the donor component with the strongest donor/recipient
log-ratio is used, so one strongly supported donor component is not diluted by
the other donor homoeologs. `mean` and `sum` remain available in the registry
(`defaults.donor_aggregation`); the dashboard exposes the setting in the
**Advanced** sidebar panel. The same panel selects calibration controls:
chosen samples are used only to estimate the cohort null, and an empty
selection pools every sample. A single sample is GC-corrected with a natural
spline and centred; a cohort uses median polish. `score = "donor_only"` is the
CNV-robust fallback.

### Ploidy

The assembly selects one or more components (or the species default), not
individual contigs, and it is built at run time from two species:

```r
pair <- resolveSpeciesAssembly(
    config,
    recipient = "triticum_aestivum",
    donor = "aegilops_speltoides",
    recipient_components = c("A", "B", "D"),
    donor_components = "S",
    default_sample_composition = c(A = 2, B = 2, D = 2)
)
```

A diploid donor/recipient, a tetraploid recipient, a hexaploid recipient, a
multi-donor design and an autotetraploid donor all use the same code path. When
a species has only one component, `recipient_components` and
`donor_components` can be omitted. A runtime `sample_composition` overrides the
default; `copy_number: 0` removes a component from the recipient aggregate.

### Masks and calibration

A mask is indexed by chromosome group and relative window. Sources are tried in
order: a fixed panel mask, a cohort recurrence mask (multi-sample), a BED
blacklist, and a physical terminal fallback. Masked windows are set to `NA`
before CBS, so segment boundaries and gap merging still see the full grid.

Calls are calibrated by segment length. The threshold is the empirical
`1 - fdr` quantile of control segment means, raised to a per-bin floor. A
single sample is its own control; a cohort uses all samples (or `controls`).

When a call is assigned to a recipient contig, MaMaMIA also estimates a
provisional dosage. Breadth is converted to a per-copy depth scale with
`-log(1 - breadth)`, so losing `k` of `n` copies gives
`lambda_observed / lambda_expected = (n - k) / n`. The call reports
`recipient_copy_number`, `recipient_observed_breadth`, `affected_fraction`,
`affected_copies` and `dosage_confidence`. A donor-only addition uses the
donor-copy estimate, and a deletion without donor support uses the recipient
loss fraction at lower confidence. Unassigned and copy-neutral calls remain
`unavailable`.

Classification uses three axes: donor gain, recipient loss and a relative
total-depth residual. A balanced replacement is `substitution`, donor gain with
a total gain is `addition`, recipient loss with a total loss is `deletion`, and
a total imbalance is `complex`. Donor and recipient cutoffs are calibrated from
control segments when controls are supplied; without controls the fixed
floors are used. `copy_state` is `balanced`, `gain`, `loss`, `compound` or
`neutral`; `class_evidence` records which side(s) supported the call,
`class_confidence` is a class-level confidence, and `class_status` is
`confident`, `ambiguous` or `needs_review`. `component_status` flags
component-heterogeneous calls and `mappability_status` flags low-MAPQ windows.
A recipient loss without donor support is `deletion`; a width mismatch is
`complex`.

## Result contract

`BamAnalysisResult` has twelve fields:

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
# or, for an existing result:
exploreIntrogressions(result)
```

The dashboard analyses every component declared for the selected species and
offers compact pill dropdowns for class, class status, copy state and
confidence, plus a passing-only switch and a width filter; class choices are
substitution, addition, deletion, complex and other. An empty dropdown means
"all".

The main tab is a vertical karyotype with a score heatmap or call-map mode. Chromosomes are
shown in assembly-config order without group or component columns, as diploid
homolog pairs: a small gap inside each pair, a larger gap between chromosomes,
the same `1A`/`1B` label on both copies, and every bar scaled by its configured
length. Every donor segment is compared with every recipient contig; a
per-contig corrected recipient residual, a standard-error-aware candidate
score, a calibrated width-matched null and a same-group homology prior feed a
per-sample global assignment with a "no host" option. Passing calls are matched
before non-passing calls, multiple non-overlapping runs on one contig are
scheduled by weighted interval scheduling, cross-group choices are reported as
`cross_group` translocation candidates rather than confirmed hosts, and a
strong assigned loss reclassifies the call as a substitution. Deleted
chromosomes and added donor chromosomes are drawn separately. Hovering a call
opens a structured card: class and status, placement, evidence residuals and
quality chips, dosage, geometry, and a group-filter hint. A sortable call
table and BEDPE export follow. See
[`docs-source/karyotype-tooltip.md`](docs-source/karyotype-tooltip.md) for the
field-by-field tooltip reference.

## Validation

The paired log-ratio recovers **6/6** reference events on the 24
*Aegilops speltoides* BAMs (Sq74 7S; Sq75 7S/2S/6S/5S; Sq76 7S), with one
control call >= 20 Mb. The donor-only fallback also recovers 6/6. A
28,000-window BAM takes about 19 s on one core; the 24-BAM cohort, including
global assignment, takes about 450 s with 30 cores. Global assignment places
five of the six reference events on recipient contigs (Sq74 7D; Sq75
2B/5B/6D/7D). The Sq76 7A event is detected as an addition without recipient
loss and receives no host. Synthetic copy-neutral events are never marked
`assigned`. The reproducible check is:

```sh
Rscript tools/validate-24.R --out test_data/introgression-study/results/validation_24bam.rds
```

Set `MAMAMIA_TEST_DATA=1` and `MAMAMIA_SLOW_TESTS=1` to run the gated real-data
tests.

## Testing

```sh
Rscript -e 'testthat::test_local()'
R CMD check --no-manual
```

The test suite covers the htslib breadth parity (samtools depth and bedtools
coverage), the WindowMetrics/ScoreTable/MaskTable contracts, diploid through
hexaploid scoring, multi-donor max/sum, GC correction, masks, CBS boundaries,
length calibration, BEDPE export, the Shiny launcher and the bundled example.
Real-data tests are gated because the 24 BAMs and the 1.2 GB gcrank indexes are
not shipped.

## Migration

See [`MIGRATION.md`](MIGRATION.md) for the configuration, result, function and
BEDPE changes since 0.99.1.

## Deployment

The production deployment mounts the registry read-only and passes its path in
`MAMAMIA_CONFIG`:

```sh
MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml R -e 'MaMaMIA::runMaMaMIA()'
```

Uploaded BAMs are staged in a private per-session directory and removed when
the session ends.

## Disclaimer

MaMaMIA is research software. Calls should be validated experimentally before
they are used for breeding or publication decisions.

## Contact

Aleksey Ermolaev <alerm@al-erm.ru>
