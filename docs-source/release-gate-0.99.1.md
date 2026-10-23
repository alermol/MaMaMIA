# Release gate 0.99.1 (T7.4)

Status: **PASS**

Date: 2026-10-21

> Post-gate update: counting was later unified on the in-tree htslib counters
> described in `bam-workflow-api.md`, and the bamsignals dependency was
> removed. This document remains the record of the original gate run. See
> `NEWS.md` for the count-backend change.



Scope: wave 10 / T7.4 of `AGENTIC_PLAN.md` — full validation and release
gate. This run validated the BAM-first workflow against the real Sq54, Sq81
and Sq84 BAMs, the wheat and rye gcrank indexes, the golden legacy tables, the
benchmark harness and the production ShinyProxy image.

No frozen decision was changed. The legacy directory-layout code was still
present during this run because D10 keeps it only until the T7.4 gate; the
T7.5 follow-up at the end of this document records its removal and the
re-validation on the decommissioned tree.

## Environment

| Component | Version |
|---|---|
| MaMaMIA | 0.99.1 |
| R | 4.6.1 |
| uvr | 0.4.7 |
| Packages in `uvr.lock` | 124 |
| Installed in `.uvr/library` | 125 (124 + MaMaMIA) |
| bamsignals | 1.44.1 |
| Rhtslib / htslib | 3.8.0 / 1.1.80 |
| GenomicRanges | 1.64.0 |
| DNAcopy | 1.86.0 |
| glmmTMB / TMB | 1.1.15.2 / 1.9.25 |
| Rcpp | 1.1.2 |
| yaml | 2.3.12 |
| testthat | 3.3.2 |
| shiny | 1.14.0 |

`glmmTMB` 1.1.15.2 loads against TMB 1.9.25 without a version warning. The
mismatch recorded in the plan's baseline is resolved by the frozen lockfile.

## 1. Reproducibility and drift guards

| Check | Command | Result |
|---|---|---|
| Frozen sync | `uvr sync --frozen --unattended` | PASS, 124 packages, "Everything is up to date" |
| Fresh-image sync | `uvr r install 4.6.1` + `uvr sync --frozen --unattended` in the Debian trixie image | PASS; the lockfile resolves from scratch |
| Environment doctor | `uvr doctor` | PASS, all project/R/library checks |
| Manifest drift | scratch copy with an extra `uvr.toml` dependency, then `uvr sync --frozen` | FAILS as required: "Lockfile is out of date with the current manifest" (exit 1) |
| Dependency-set drift | injected an extra `uvr.toml` dependency and ran `tests/testthat/test-build.R` | FAILS as required: `fakepkg` reported as unexpected in `[dev-dependencies]` |
| DESCRIPTION/uvr.toml agreement | `tests/testthat/test-build.R` | PASS |
| Conda scope | `conda list -n mamamia` | 0 R/Bioconductor packages; only `samtools`, `bedtools`, `pandoc` and non-R tools |

Operational note: `uvr sync --frozen` prunes MaMaMIA itself from
`.uvr/library` (it is the project package, not a lockfile dependency). Run
`R CMD INSTALL --preclean .` after every frozen sync; the Dockerfile follows
this order.

## 2. Build and package check

```sh
uvr sync --frozen
source .uvr/activate
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'roxygen2::roxygenise()'
R CMD INSTALL --preclean .
R CMD build .
R CMD check --no-manual --output=check-output MaMaMIA_0.99.1.tar.gz
```

- `compileAttributes()` and `roxygenise()` reproduce the committed generated
  files; no diff was produced.
- `R CMD INSTALL --preclean .` succeeds with the C++23 build and Rhtslib
  linking.
- `R CMD build` succeeds, including the vignette.
- `R CMD check` finishes with **1 WARNING and 1 NOTE, no ERRORs**, and the
  test suite inside the check passes:
  - WARNING: GNU extensions in `src/Makevars` (`OBJECTS`, `$(shell ...)`). This
    is pre-existing and required by the vendored gcrank/Rhtslib build.
  - NOTE: `Rhtslib` is declared in `Imports` but only used through `LinkingTo`.
    This is pre-existing; the package links to Rhtslib's headers and static
    library.

## 3. Test suites

| Suite | Command | Result |
|---|---|---|
| Fast suite | `testthat::test_local(".", reporter="summary")` | All tests pass; only the 2 `MAMAMIA_SLOW_TESTS` golden skips |
| Golden regression | `MAMAMIA_SLOW_TESTS=1 MAMAMIA_TEST_DATA=1 testthat::test_file("tests/testthat/test-golden-bam.R")` | All tests pass |
| Real-data integration | `MAMAMIA_TEST_DATA=1 Rscript tests/integration/run-bam-e2e.R --cores=8` | 33 checks, 0 failures |

## 4. Golden regression evidence

Against the legacy Sq54 `w1000` BED tables:

| Metric | Result | Gate |
|---|---|---|
| Count mismatches | 229 / 28,000 windows (0.818%), recipient 229/21,000, donor 0/7,000 | ≤2% |
| Maximum per-window count delta | 14 reads | ≤20 |
| Total absolute delta | 708 reads (0.020% of legacy counts) | <0.1% |
| Mismatch explanation | 229 / 229 are reads whose 5' end is outside the window; max unexplained 0 | all explained |
| Per-contig counts | equal to `samtools view -q 20 -F 2308 -c` | exact |
| GC | max absolute delta 4.9995e-7 over 28,000 windows | matches at six decimals |
| Window tiling | exact-N tiling, 28,000 windows; at most 963 bp shift vs legacy fixed-width tiling | documented |
| Full pipeline segments | 7 pairs; segment-count delta ≤0; putative delta ≤0; introgression bp relative delta 9.272e-7 | ≤1 segment, ≤2% bp |
| Same-interval correction | 38 vs 38 segments; breakpoints identical; max median delta 0.3906 (0.1222% of scale); 99% relative 2.024e-3, max 0.1427 on near-zero-coverage outliers | documented tolerances |
| Putative calls | 1 vs 1 on the same intervals | exact |

The count difference is the documented start-based (`bamsignals`) vs
overlap-based (`bedtools coverage`) boundary semantics, not a regression.

## 5. Real-data end-to-end scenarios

All analysis scenarios (1–4) used the 28-contig wheat + rye assembly, two
gcrank indexes, `mapq = 20`, `exclude_flags = 2308`, and `cores = 8`. Scenario
0 only resolves the species registry and runs no BAM analysis.

| # | Scenario | Elapsed | Peak RSS | Result |
|---|---|---:|---:|---|
| 0 | Role-free species registry; resolve wheat/rye in both directions | instant | — | 7 pairs each way; swapping reverses every pair |
| 1 | Sq54 single BAM, 1,000 windows | 96.9 s | 5,285 MB | 28,000 windows; 38 segments, 1 putative; BEDPE 1 row, 15 columns |
| 2 | Sq54 + Sq81 + Sq84, median normalization | 267.2 s | 7,525 MB | 7,000 × 3 donor / 21,000 × 3 recipient; finite dose everywhere; 29 segments, 21 putative; combined BEDPE 21 rows from 2 samples |
| 3 | Sq54 as reference sample | 219.0 s | 7,893 MB | reference dose exactly 1; other doses finite; 33 segments |
| 4 | Sq54 high resolution, 10,000 windows, `max_fit_windows = 50000` | 137.0 s | 10,190 MB | 280,000 windows; 271,712 ideal → 50,000 fitted (seed 44,170,674); 26 segments |

Notes:

- Single-BAM runs correctly resolve to normalization `none` and segment on
  `cor.gc`; cohort runs write `dose` and segment on `dose`.
- The high-resolution peak of about 10.2 GB is inside the documented 8–12 GB
  envelope and the 16 GB container limit. The earlier 16.4 GB reading was an
  artefact of the runner retaining three prior results; the runner now releases
  each result before the next scenario, and the external monitor reports a
  process peak of 10,437,416 kB (about 10.0 GiB).
- The median cohort's combined BEDPE contains rows from Sq81 and Sq84 only;
  Sq54 had no putative introgression after median normalization. The export
  correctly includes a `sample` column and does not invent empty rows.

## 6. Benchmark

`tools/benchmark-bam-pipeline.R --bam=test_data/bam_files/Sq54/cs_t2t_se_cereale.bam --resolution=1000 --cores=8`

- 20 / 20 golden checks pass.
- Wall time 324 s; peak RSS 7,741,620 kB (about 7.4 GiB).

| Step | Elapsed |
|---|---:|
| `countReads` on legacy intervals | 54.760 s |
| gcrank GC on legacy intervals | 0.279 s |
| `makeAnalysisWindows` | 0.071 s |
| `countReads` on generated windows | 44.979 s |
| `analyzeBamInput` (full pipeline) | 84.558 s |
| legacy RCA + correction | 29.726 s |
| legacy segments | 3.649 s |
| `buildSampleTables` on legacy intervals | 0.089 s |
| BAM-first core on legacy intervals | 33.261 s |

Fingerprint: package 0.99.1, R 4.6.1, htslib 1.1.80, gcrank commit
`8fdbc32557973c7737496959ecb91acd9bed3ee5`; 229 count mismatches, max delta 14,
total absolute delta 708, max GC delta 4.9995e-7.

## 7. Deployment smoke test

- Image built from `deployment/Dockerfile`:
  `mamamia:0.99.1` (`sha256:47e71fbc3fbfd160ed0c66f537cee61aef7e908a3afe5eaacc6d4554bc927842`).
- `uvr doctor --unattended` inside the image: PASS. (`cargo` is absent by
  design; every R package is installed from `uvr.lock`.)
- `uvr run /app/deployment/smoke-test.R` with the example species registry,
  read-only gcrank indexes and read-only BAM root, under
  `--memory=16g --cpus=4`:
  - species pair `triticum_aestivum__A__secale_cereale__R`;
  - 14,000 windows; 39 segments; 1 putative introgression;
  - result is a `BamAnalysisResult` with sample `Sq54`.
- ShinyProxy itself is not required for the image build or this smoke test; it
  is only needed for the proxied deployment.

## 8. Frozen decision audit

| ID | Decision | Evidence | Status |
|---|---|---|---|
| D1 | Pairwise donor − recipient difference | `R/segmentation.R` computes `diff = <value>.don - <value>.rec`; `test-segmentation.R`; golden segment summary | PASS |
| D2 | One `@RG SM`, basename fallback, duplicate samples rejected | `R/bam-header.R`, `R/preflight.R`; `test-bam-header.R`, `test-preflight.R` | PASS |
| D3 | Always CSI, `min_shift = 14` | `R/bam-index.R` default and `.csi` only; `test-bam-index.R`; staged BAMs indexed in the integration run | PASS |
| D4 | Per-genome index list, route each contig | `R/config.R`, `R/gcrank.R` groups queries by `index_path`; `test-config.R`, `test-gc-extraction.R`; two indexes in the integration run | PASS |
| D5 | Pseudocount ratio, per-side default | `R/normalization.R`: `(C + eps) / (denominator + eps)`, `eps = 0.5`, `by_side = TRUE`; `test-normalization.R`; integration median/reference checks | PASS |
| D6 | Normalize after correction, before CBS; store `dose` | `R/pipeline-bam.R` correction → normalization → segmentation; `test-pipeline-bam.R` stage-order tests | PASS |
| D7 | `segments(value_col = NULL)` uses `dose` if present, else `cor.gc` | `R/segmentation.R` resolves the column and records `ISA$meta$value_col`; integration checks both modes | PASS |
| D8 | No Rsamtools; Rhtslib header/index + bamsignals counts | `DESCRIPTION` imports Rhtslib and bamsignals only; no Rsamtools reference; `test-build.R` | PASS |
| D9 | Deterministic stratified fit subsample, predict all | `R/correction.R`, `R/fit-cache.R`; `test-correction.R`; high-resolution run fitted 50,000 of 271,712 ideal windows with a recorded seed | PASS |
| D10 | Keep legacy code through T7.4 as the regression oracle | `R/app.R`, `R/input.R`, `R/validation.R`, `R/dashboard.R` still present; golden test passes | PASS (see readiness note) |
| D11 | Per-sample ISA, UI sample selector, combined BEDPE | `R/app-bam.R` sample selector; `R/bedpe.R` combined export with `sample`; `test-dashboard-bam.R`, `test-bedpe.R`; integration BEDPE checks | PASS |
| D12 | Uploads staged in session temp; server-side BAM directories removed | `R/preflight.R`, `R/bam-index.R` staging; `test-preflight.R`, `test-dashboard-bam.R` | PASS |
| D13 | uvr (not renv) | `uvr.toml`, `uvr.lock`, `.r-version`; `uvr sync --frozen` and `uvr doctor` pass; drift guards verified | PASS |
| D14 | Conda only for non-R tools; uvr for R packages | Dockerfile has no conda; conda `mamamia` has 0 R packages and only samtools/bedtools/pandoc | PASS |

No approved deviations were required.

## 9. Definition-of-done mapping

| DoD item | Evidence |
|---|---|
| 1–4 BAM ingestion, assembly detection, CSI, 100/1,000/10,000 windows | Fast suite; integration scenarios; real config uses 28,000 and 280,000 windows |
| 5 GC at six-decimal precision | Golden test and benchmark: max delta 4.9995e-7 |
| 6 Counts match up to boundary-crossing reads | Golden test and benchmark: 229/28,000, all explained, max 14 |
| 7 `correctReadCounts()` path | Golden same-interval correction and full pipeline checks |
| 8 Median/reference normalization | Integration scenarios 2 and 3; `test-normalization.R` |
| 9 CBS per sample on the pairwise difference | Integration segment tables; `test-segmentation.R` |
| 10 Dashboard plots and lists introgressions | `test-dashboard-bam.R` passes in the fast suite |
| 11 BEDPE with `sample` column | Integration single and cohort exports; `test-bedpe.R` |
| 12 Legacy removal | Deferred to T7.5 by D10; the readiness note below clears it |
| 13 `R CMD check` | 0 ERRORs; 1 pre-existing WARNING, 1 pre-existing NOTE |
| 14 Real-data regression | Golden test and benchmark |
| 15 ShinyProxy image/config smoke test | Section 7 |
| 16 Documentation | README, vignette, function docs, config reference, runbook; docs fixture examples run |
| 17 10,000-window memory budget | 10.2 GB peak, inside the documented 8–12 GB and 16 GB limit |
| 18 Frozen decisions | Section 8, no deviations |
| 19 Reproducible environment | Section 1; `uvr sync --frozen`, `uvr doctor`, zero R packages in conda |

## 10. D10 readiness note for T7.5

The BAM-first path passes every release gate and covers the legacy
directory-layout workflow's purpose:

- direct BAM ingestion replaces prepared `mamamia_analysis` tables, including
  header-based assembly/sample detection and CSI creation;
- GC correction uses the same `correctReadCounts()` implementation and matches
  the legacy GC at six-decimal precision;
- pairwise CBS reproduces the legacy segment breakpoints, medians and putative
  calls on the same intervals;
- the dashboard exposes the same plots and BEDPE export, now with a per-sample
  ISA selector and a combined cohort BEDPE carrying a `sample` column;
- multi-sample median/reference normalization is an additive capability that
  the legacy workflow did not have.

The only legacy-specific capability not reproduced by BAM-first is reading
prepared directory-layout tables (`readMaMaMIAInput()` and
`launcher_app()`), which is superseded by direct BAM ingestion and was removed
in T7.5. There is no capability loss that blocks decommissioning. **T7.5 may
proceed.**

## Reproduction commands

```sh
uvr sync --frozen --unattended
source .uvr/activate
uvr doctor
Rscript -e 'Rcpp::compileAttributes()'
Rscript -e 'roxygen2::roxygenise()'
R CMD INSTALL --preclean .
NOT_CRAN=true Rscript -e 'testthat::test_local(".", reporter="summary")'
MAMAMIA_SLOW_TESTS=1 MAMAMIA_TEST_DATA=1 \
  Rscript -e 'library(MaMaMIA); testthat::test_file("tests/testthat/test-golden-bam.R")'
MAMAMIA_TEST_DATA=1 Rscript tests/integration/run-bam-e2e.R --cores=8
Rscript tools/benchmark-bam-pipeline.R \
  --bam=test_data/bam_files/Sq54/cs_t2t_se_cereale.bam \
  --resolution=1000 --cores=8 --out=/tmp/bench-bam
R CMD build .
mkdir -p check-output
R CMD check --no-manual --output=check-output MaMaMIA_0.99.1.tar.gz
```

## T7.5 follow-up: legacy decommissioning

Status: **PASS**

Date: 2026-10-21

Precondition: the D10 readiness note in section 10 ("T7.5 may proceed").
No frozen decision was changed by this task.

Removed from the package:

* `launcher_app()`, the `analysis_dir` / `assemblies_dir` arguments and the
  `input = "auto"` / `"bam"` / `"files"` / `"example"` dispatch of
  `runMaMaMIA()`. The entrypoint is now
  `runMaMaMIA(cache, config, max_request_size, host, port,
  launch.browser)` and always opens the BAM-first app.
* `R/input.R`: `readMaMaMIAInput()`, `listMaMaMIASamples()`,
  `listMaMaMIAReferences()`, `list_window_labels()`, `read_window_table()`,
  `read_chromosome_meta()` and the directory discovery helpers.
* `tests/testthat/test-input.R` (replaced by `test-default-pairs.R` for the
  retained pairing helper).

Retained: `RCA()`, `correctReadCounts()`, `segments()`, `ISA`, `writeToBedpe()`,
`writeCohortToBedpe()`, `buildSampleTables()`, `defaultTargetPairs()` and the
BAM-first API. `chr_number()` moved to `R/config.R`; `defaultTargetPairs()`
moved to `R/sample-tables.R` with updated documentation.

The gated golden regression and `tools/benchmark-bam-pipeline.R` keep the
legacy BED oracle by reading the four-column tables locally
(`read_golden_legacy()` in `tests/testthat/helper-bam.R` and
`read_legacy_layout()` in the benchmark). No directory-layout reader ships in
the package.

### Validation on the decommissioned tree

| Gate | Result |
|---|---|
| Fast suite (`test_local`) | All pass; only the two expected gated slow-test skips. |
| Golden (`MAMAMIA_SLOW_TESTS=1 MAMAMIA_TEST_DATA=1`) | All pass: 229 / 28,000 count mismatches (all 5'-end explained, max 14, 0.020 %), GC max 4.9995e-7, 7 pairs, 38 vs 38 segments, identical breakpoints, median max delta 0.3906 (0.1222 %), putative 1 vs 1. |
| Real-data E2E (`run-bam-e2e.R --cores=8`) | 33 checks, 0 failures; single 38 segments / 1 putative, median 29 segments / 21 putative, reference dose exactly 1, high-resolution 280,000 windows with 50,000 fitted (seed 44170674) and 26 segments; reference 218.7 s / 7,903 MB, high-resolution 136.3 s / 10,505 MB. |
| Benchmark (`--resolution=1000 --cores=8`) | 20 / 20 checks pass; 38 vs 38 segments, 1 vs 1 putative, breakpoints identical, max median delta 0.1222 %; `analyzeBamInput` 84.5 s, count 45.1 s. |
| `uvr sync --frozen` / `uvr doctor` | Pass (126 packages in the lockfile; doctor reports no issues). |
| `R CMD build` / `R CMD check --no-manual` | Build OK; check 0 ERRORs, 1 WARNING and 1 NOTE, both pre-existing (`src/Makevars` GNU extensions; unused `Rhtslib` Import). Tests inside check pass. |
| Docker image | `mamamia:0.99.1-t7.5`, image `sha256:82263a13d9c641605d6b21794e60d908d10fcf97f9095df2eaf67f6c78c333f7`, 2.66 GB; `uvr doctor` passes. |
| Container smoke test | Exit 0: species pair `triticum_aestivum__A__secale_cereale__R`, 14,000 windows, 39 segments, 1 putative. |

Definition of done item 12 ("the legacy directory-layout and packaged-example
workflows are removed in T7.5") is satisfied. README, the overview vignette,
`docs-source/bam-workflow-api.md` and `NEWS.md` document the breaking change
and the migration path for users with precomputed tables.

Note on concurrent edits: the `interacCircos` declaration that appeared in
`uvr.toml` during T7.5 was completed in the follow-up below (`DESCRIPTION`,
usage and tests). The transient unrelated `report/` directory is excluded from
the package build via `.Rbuildignore`.

## Circos follow-up (wave 11)

Status: **PASS**

Date: 2026-10-21

* Removed the remaining BioCircos-related code and assets: the disabled
  `R/circos.R` module, `inst/mamamia/mamamia-circos.*`, the inactive explorer
  card and server hooks, and their tests.
* Added an interactive circos overview built with `interacCircos`: chromosome
  ideograms, a CNV track with the per-window donor - recipient difference, an
  arc track marking the visible segments on both chromosomes, and links
  joining each donor and recipient interval. The view follows the explorer's
  chromosome, pair, width and putative filters.
* `interacCircos` and `htmlwidgets` are declared in `DESCRIPTION` `Imports` and
  in `uvr.toml` dependencies; `uvr.lock` was regenerated (126 packages).
* New `R/circos.R` helpers prepare the data independently of the widget;
  `isa_dashboard_circos_output()` and `isa_dashboard_render_circos()` provide
  the Shiny bindings that interacCircos does not ship itself.

### Validation

| Gate | Result |
|---|---|
| Fast suite (`test_local`) | All pass, including the new `test-circos.R` and `test-build.R`; only the two expected gated slow-test skips. |
| `R CMD build` / `R CMD check --no-manual` | Build OK; 0 ERRORs, 1 WARNING and 1 NOTE, both pre-existing (`src/Makevars` GNU extensions; unused `Rhtslib` Import). Tests inside check pass. |
| Widget serialization | Demo ISA at 100 windows: modules `BACKGROUND_cnv_background`, `CNV_difference`, `ARC_segment_spans`, `LINK_introgressions`; 51,720-byte HTML with the `interacCircos` binding. |
| Docker image | `mamamia:0.99.1-interac`, `sha256:f5b4d68747883d0d98034bcc64988f4a7b83641a69e6a24b5f1d6403637b0469`, 2.66 GB; `uvr doctor` passes. The first two build attempts hit truncated P3M binary downloads (`minqa`, `RcppEigen`); the third, uncached attempt succeeded. |
| Container smoke test | Exit 0: species pair `triticum_aestivum__A__secale_cereale__R`, 14,000 windows, 39 segments, 1 putative. |

`interacCircos` is an unused-import-free runtime dependency: the `R CMD check`
NOTE now lists only `Rhtslib`.
