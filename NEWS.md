# MaMaMIA 0.99.1

## Grid dashboard layout

* The BAM-first dashboard now uses the `gridlayout` package with the layout
  `2rem` gap, `200px` sidebar, `150px` header, fluid main area and `120px`
  footer. The application name is shown in the header, and the main area holds
  only the **Introgression circos**, **Coverage / dosage** and **Segments**
  tabs.
* The sidebar starts with a **Single BAM / Multiple BAMs** selector, the upload
  control and an **Example data** button for both modes. Analysis controls
  appear after preflight succeeds and live in single-open collapsible panels
  (**Analysis**, **Assembly**, **Advanced**, **Filters**, **Export**).
* Preflight status, detected assemblies, assembly candidates, problems,
  suggested fixes and every pipeline stage are written to the calculation log
  in the footer. The standalone preflight, problems and detected-window cards
  are no longer part of the main area.
* The bundled example is now loaded with the **Example data** button; pressing
  it again clears the example and restores the registry resolution choices.
  The button is available in both single-BAM and multi-BAM modes.

## Upload-only BAMs and bundled example data

* BAMs are accepted only through the Shiny upload control. `server.bam_root`,
  the server-side BAM picker, the `bam_paths` argument of `runMaMaMIA()` and the
  `/srv/mamamia/bams` container mount are removed. An old registry key is
  ignored; `server.max_bams` and `server.max_total_bytes` still apply.
* Added `inst/extdata/example`: a self-contained registry and pre-computed
  `BamAnalysisResult` tables at 100 and 1,000 windows per chromosome. The live
  sidebar offers the **Example data** button; it replays the normal pipeline
  stage messages (preflight, indexes, windows, GC, counts, correction,
  segmentation) while loading the stored result, so no example BAM is needed in
  the container. Regenerate it with `tools/make-example-data.R`.
* `deployment/smoke-test.R` now uses the bundled demo BAM and registry by
  default; set `MAMAMIA_SMOKE_BAM` and `MAMAMIA_CONFIG` to test a production
  registry instead. Deployment and API documentation reflect the upload-only
  model.

## Deployment fix: fontconfig and cairo text metrics

* The production image now installs `fontconfig` and `fonts-dejavu-core`. uvr's
  R bundles the cairo and fontconfig shared libraries but no system
  configuration or fonts, so the first `renderPlot()` in a container aborted the
  whole R process with `Fontconfig error: Cannot load default config file`
  followed by a segfault in `grid.Call(C_stringMetric, ...)`.
* `deployment/Dockerfile` renders a cairo/ggplot PNG during the build, and
  `deployment/smoke-test.R` renders a text probe and the real ISA plot. A
  missing font setup now fails the build or the smoke test instead of the
  dashboard. `deployment/start.R` warns at startup when cairo is enabled but no
  fontconfig configuration is visible.

## Dashboard layout fix: fillable pages

* The launcher and the standalone explorer now use
  `bslib::page_sidebar(fillable = FALSE)`. The default fillable page makes the
  growing stack of result cards flex items, which shrank the 420 px plot output
  to a few pixels (or zero); Shiny's PNG device then failed with
  `invalid 'height' argument` before the first plot could render. The pages
  scroll normally, and `test-dashboard-bam.R` asserts that the launcher main
  area is not a fill container.

## Documentation and release preparation (wave 9)

* The BAM-first workflow is now the documented primary path. The README has a
  complete quick start (`uvr sync --frozen`, `runMaMaMIA(config = ...)`), the
  species registry reference, resolution and normalization formulas, counting
  semantics, resource requirements and the ShinyProxy pointer.
* Added `inst/extdata/bam-demo`: a self-contained 140 KB synthetic registry with
  two role-free species, a variable-GC gcrank index and a BAM. The complete
  pipeline runs on it and writes a BEDPE, so the README, the vignette and the
  roxygen examples are executable without external data. Regenerate it with
  `tools/make-docs-fixture.R`.
* Added roxygen examples for `readMaMaMIAConfig()`, `bamHeaderInfo()`,
  `ensureBamIndex()`, `matchAssembly()`, `makeAnalysisWindows()`, `gcrankGc()`,
  `countReads()`, `normalizeReadCounts()`, `analyzeBamInput()` and
  `writeCohortToBedpe()`, and regenerated `man/`.
* The vignette starts with a BAM-first walkthrough (species registry, species
  roles and swapping, resolution and normalization, the demo fixture, cohort
  BEDPE) and documents the count/flag defaults, 5'-end boundary semantics and
  how to verify counts with `samtools`. The directory-layout walkthrough was
  marked deprecated at that point; the `RCA()`/`ISA` sections remain the core
  API reference.
* Documented the uvr migration (`uvr.toml`, `uvr.lock`, `.r-version`,
  `.uvr/library/`, `uvr sync --frozen`, no renv) in the README, and the
  deployment runbook in `docs-source/shinyproxy.md` now points at the same
  entrypoint.

## Breaking change: directory-layout workflow removed (wave 11)

* The legacy directory-layout workflow is removed. `launcher_app()`, the
  `analysis_dir` / `assemblies_dir` arguments of `runMaMaMIA()`, the
  `input = "files"` / `"example"` modes, `readMaMaMIAInput()`,
  `listMaMaMIASamples()`, `listMaMaMIAReferences()` and the packaged-example
  input are gone. `runMaMaMIA(config = ..., bam_paths = ...)` is the single
  Shiny entrypoint and `analyzeBamInput()` is the scripted pipeline. `RCA()`,
  `correctReadCounts()`, `segments()`, `buildSampleTables()`,
  `defaultTargetPairs()`, `writeToBedpe()`, `writeCohortToBedpe()` and the
  `ISA` explorers remain supported.
* The minimum supported environment is unchanged: R 4.6.0 or newer (the
  reproducible `uvr.lock` pins R 4.6.1) and a C++23 compiler. `renv` is not
  used; `uvr sync --frozen` reproduces the package set.
* The gated golden regression now reads its legacy BED oracle with a local
  test helper, and `tools/benchmark-bam-pipeline.R` carries its own oracle
  reader, so no directory-layout reader ships in the package.

## Dashboard cleanup (wave 11)

* Removed all BioCircos-related code and assets: `R/circos.R`, the
  `inst/mamamia/mamamia-circos.*` files, the disabled circos card, the server
  hooks that fed it and their tests. The explorer keeps the mirror/difference
  plot, the segment table and the BEDPE export.
* Added an interactive circos overview built with `interacCircos`: chromosome
  ideograms, the per-window donor-recipient difference, called segment spans and
  donor-recipient links. It follows the explorer's chromosome, pair, width and
  putative filters. `interacCircos` is declared in `DESCRIPTION` and
  `uvr.toml`.

# MaMaMIA 0.99.0

## Initial submission

* Initial submission.

## BAM-first workflow (species selection)

* The registry now supports a role-free `species:` section. A species declares
  only its label, gcrank indexes and contigs (with an optional `group` and
  `default_subgenome`); it is not tied to the donor or recipient role. A
  species-only config needs at least two species, and `assemblies` may be
  omitted entirely.
* Added `listMaMaMIASpecies()` and `resolveSpeciesAssembly()`. Resolution
  selects one subgenome per side and pairs contigs by chromosome group (the
  optional `group` field, otherwise the first number in `name`), includes only
  the paired contigs, and returns the derived assembly with the configuration
  that contains it. `unpaired` explains every excluded contig. An explicit
  `pairs` argument remains available as an override.
* Named assemblies may now reference species (`recipient: wheat`, `donor: rye`),
  including optional `recipient_subgenome` / `donor_subgenome` and `pairs`.
  Inline `contigs:` assemblies keep working in the scripted API; the app falls
  back to BAM-header auto-detection when a registry has no species section.
* The BAM-first sidebar is species-only: choose a recipient species and a
  donor species (in either direction), with subgenome selectors only when a
  species has several subgenomes. A "Swap donor and recipient" button exchanges
  the two species and carries their subgenome choices with them. The resolved
  pairs and any unpaired contigs are shown before Run; the previous result is
  invalidated when the pair changes. The predefined-assembly selector has been
  removed.
* Auto-detection is available as a checkbox when the registry declares named
  assemblies, without exposing an assembly list; a registry without species
  falls back to BAM-header detection automatically.
* `inst/config/assemblies.example.yml` now demonstrates the species registry,
  `tools/validate-config.R` reports and checks species indexes, and the BAM API
  note documents the pairing rule and the ShinyProxy registry-growth caveat.

## BAM-first workflow (wave 8)

* The cohort card now offers a per-sample BEDPE export next to the combined
  selected export. It follows the sample selector: the rows selected for the
  current sample, or all of that sample's putative calls when no rows are
  selected. File names are sanitized per sample, and the cohort export keeps
  the `sample` column last.
* Added the production deployment artifacts: a Debian trixie image with a
  pinned uvr release and R 4.6.1, `uvr sync --frozen` from `uvr.lock`, a
  non-root runtime user and `OMP_NUM_THREADS=1`; a ShinyProxy spec with
  read-only registry/index/BAM mounts and 16 GB / 4 CPU limits; a ShinyProxy
  entrypoint that reads `MAMAMIA_CONFIG`; and an nginx example with large
  uploads and long analysis timeouts.
* Added `tools/validate-config.R` for production registries (config, gcrank
  payloads and `server.bam_root`) and documented the data prerequisites, build,
  local smoke test, upload limits, resource sizing and cleanup in
  `docs-source/shinyproxy.md`.

## BAM-first workflow (wave 7)

* The BAM-first progress bar now reports pipeline detail, not only the stage
  name: BAM count and assembly, reused/staged indexes, total and per-chromosome
  window counts, reads per sample, ideal fit windows, normalization mode and
  reference, and segments per sample.
* Every run and preflight error now carries the affected source BAM(s) in
  addition to the stage, message and suggested fix. Sample errors are mapped
  back through the preflight table, so notifications and the Problems card show
  the supplied path instead of the staged session copy. Run failures also stay
  visible in a persistent "Last run error" card; the previous result is kept.
* The Run button is disabled, with the reason shown, when a contig is shorter
  than the selected resolution. The same guard runs inside the handler, and the
  small correction/GC failure and zero-coverage paths remain actionable
  (`Fewer than N ideal windows`, no reads on a side, all-zero reference
  windows). Zero/near-zero coverage warnings are repeated as a notification
  after a successful run.
* Invalid run controls (resolution, MAPQ, exclude flags, cores, fit cap,
  reference sample) are now MaMaMIA conditions, so they use the same
  stage/message/suggestion presentation as pipeline failures.

## BAM-first workflow (wave 6)

* `runMaMaMIA()` now dispatches between the BAM-first app and the legacy
  directory app. `input = "auto"` chooses the BAM app when `config` or
  `bam_paths` is supplied; `"bam"` forces it, and `"files"`/`"example"` keep
  the legacy launcher reachable until T7.5. The new `config`, `bam_paths`,
  `input` and `max_request_size` arguments are documented; `max_request_size`
  sets `options(shiny.maxRequestSize)` (5 GiB by default).
* The BAM-first app (`R/app-bam.R`) shows the configuration summary, uploads and
  a bounded server-side BAM picker, resolution/reference/assembly controls and
  advanced MAPQ, exclude-flag, core and fit-window inputs. Preflight is staged
  into a private session directory (symlink preferred) and displayed as a
  per-BAM table with index action and warnings; the detected-assembly card
  lists donor/recipient contigs, pairs and window counts for every allowed
  resolution; assembly failures show candidate scores. The Run button stays
  disabled until preflight succeeds.
* The run flow reports pipeline stages through `withProgress()`, catches every
  failure with the stage and a suggested fix, keeps the previous result on
  error, and invalidates a result as soon as an input that affects it changes.
  The session blocks during analysis (no background process); the UI says so.
  `server.max_bams` and `server.max_total_bytes` are enforced by preflight, and
  zero/near-zero coverage, all-zero reference windows and failed fits are
  reported without escaping the Shiny server.
* Multi-sample results gain a sample selector, a cohort summary (reads on
  donor and recipient targets, segments, putative introgressions and
  normalization mode), a per-sample selection map and a combined
  "Download selected BEDPE (all samples)" export through
  `writeCohortToBedpe()`. Switching samples restores the rows selected for the
  other sample.
* Uploaded and server-side BAMs are staged read-only into the session temp
  directory; the source BAM and `server.bam_root` are never written to. Small
  files (below 1 MiB) are copied, larger files are symlinked, so
  temporary-directory cleanup cannot touch the source metadata. The session
  directory, staged BAMs and CSI indexes are removed on session end. GC index
  handles stay cached per process; no session clears them.
* `explorer_server()` gained `selection` and `on_selection`, which preselect
  rows from a parent-owned map and report the current selection per sample.
* `correctReadCounts()` gained `max_fit_windows` (default `Inf`) for
  high-resolution data. Above the cap, the GC fit uses a deterministic sample
  stratified by subgenome and GC decile, with at least ten windows per
  non-empty stratum and predictions still covering every window. The RCA
  records `fit_sampling`, the fit cache key covers the cap and the selected
  rows, and `analyzeBamInput()` applies the cap per sample. A 280,000-window
  benchmark fits 50,000 ideal windows in 47 s with a 4.2 GB peak RSS; see
  `inst/benchmarks/`.

## BAM-first workflow (wave 5)

* The introgression explorer now renders the thresholded `ISA` as a mirror or
  donor-recipient difference plot next to the segment table. The plot follows
  the chromosome, pair, minimum-width, call-threshold and putative-only
  filters, and shows a validation message instead of an empty plot when no
  analysis or no chromosome pair is selected. `plot.ISA()` gained `max_points`
  for point thinning; the dashboard caps the drawn points at 20,000 while
  keeping every segment line and the full export data.
* `writeToBedpe()` gained `sample`. A supplied sample name is appended as the
  last column; `sample = NULL` keeps the legacy BEDPE layout byte for byte.
* `writeCohortToBedpe()` combines putative calls, or selected segment rows,
  from several samples of a `BamAnalysisResult` (or a named list of `ISA`
  objects) into one BEDPE with a `sample` column. It keeps the header when no
  row is selected and validates sample names and row ids per sample.

## BAM-first workflow (wave 4)

* `normalizeReadCounts()` turns the GC-corrected coverage of a named list of
  `RCA` objects into dosage. Reference mode uses
  `(C_s + eps) / (C_ref + eps)`; median mode uses
  `(C_s + eps) / (median_s(C_s) + eps)`. `cor.gc` is preserved, `dose` is added
  to the RCA table, and non-finite doses are marked not ideal. Donor and
  recipient windows are normalized independently by default; `by_side = FALSE`
  is an explicit sensitivity option that pools both sides per window position.
* `segments()` gained `value_col`: it consumes `dose` when present and falls
  back to `cor.gc` otherwise. `ISA$meta$value_col` records the choice,
  `reverseWindows()` reverses and re-segments dose-aware results, and
  `plot.ISA()` plots dosage with a `Dosage` axis. Segmentation without a dose
  column is unchanged.
* `analyzeBamInput()` now normalizes multi-BAM cohorts after all corrections
  and before any CBS. `normalization = "median"` (the default for cohorts) and
  `"reference"` are supported; `"reference"` requires `reference_sample` and
  gives that sample a dose of exactly 1. A single BAM still uses `"none"`.
  Results gained a `samples` metadata table, provenance gained `n_samples` and
  `total_windows`, and a warning guards cohorts whose total window count
  exceeds `max_fit_windows`.

## BAM-first workflow (wave 3)

* `buildSampleTables()` joins generated windows, gcrank GC and bamsignals
  counts into the `SampleTables` contract and validates the result through
  `RCA()`. Explicit `assembly$pairs` win over the chromosome-number rule, and
  malformed samples fail with the sample name before the GC-bias fit.
* `analyzeBamInput()` turns one or more BAMs into `BamAnalysisResult`: preflight,
  CSI indexing (staging read-only BAMs in session temp), exact-N windows,
  gcrank GC, read counting, per-sample RCA, GC correction and pairwise CBS. It
  records provenance and reports `preflight`, `index`, `windows`, `gc`,
  `counts`, `correction` and `segmentation` progress stages. A single BAM is
  always normalized with mode `"none"`; multi-sample normalization arrives in
  wave 4.
* `countReads()` now documents bamsignals' actual 5'-end semantics: forward
  reads are assigned by their leftmost coordinate and reverse reads by their
  rightmost coordinate. The Sq54 golden regression explains all 229 legacy
  count differences with this rule, and a reverse-strand boundary test locks it
  down.
* `tools/benchmark-bam-pipeline.R` and the opt-in
  `tests/testthat/test-golden-bam.R` compare legacy-interval counts and GC,
  exact-N windows, and the full segment summary with the legacy
  directory-layout pipeline. The known tiling difference (exact-N windows vs
  fixed-width BED windows, at most 999 bp) is documented rather than hidden.

## Dashboard

* The previous circos view was removed together with its dependency: the
  dashboard shows the segment table with filtering and the BEDPE export. A new
  `interacCircos` overview is documented under the 0.99.1 dashboard cleanup.
* `runMaMaMIA()` is the new main entry point: a locally running Shiny app that
  reads the windowed input layout, runs RCA construction, GC-bias correction
  and segmentation, and hands the result to the interactive explorer. The
  explorer part can still be launched on its own with
  `exploreIntrogressions()`.
* `readMaMaMIAInput()`, `listMaMaMIASamples()` and `listMaMaMIAReferences()`
  read the on-disk layout (per-sample coverage tables, per-reference GC-content
  tables and chromosome metadata), and `defaultTargetPairs()` pairs donor and
  recipient chromosomes by chromosome number within selected subgenomes.
* The explorer gained call threshold and minimum segment width controls that
  re-flag the segments through `subset.ISA()` without re-running CBS, plus
  select-all and clear buttons for the BEDPE export.

## Performance

* `correctReadCounts()` no longer calls `glmmTMB::predict()` twice on the full
  window table. The fitted ZINB response is evaluated once for the distinct
  (GC, subgenome) pairs, with design matrices built from the stored smooth rather
  than via `predict(debug = TRUE)`, dropping that step from ~22 s to ~1.4 s on the
  packaged `triticum` example. Corrected coverage agrees with the previous
  implementation to ~1e-14 relative. The remaining cost is the `glmmTMB` fit
  (~28 s at `cores = 4`).
* A new `cache` argument memoises that fit to an `.rds` file, keyed on the
  fitted windows, the model specification, `cores` and the R, package and
  `glmmTMB` versions, so repeated calls on the same windows reuse it (~16x
  faster in the example).

## Reproducibility

* `segments()` now pins the seed of the permutation-based p-values used by
  `DNAcopy::segment()` through a new `seed` argument (default `1L`). CBS was
  previously unseeded, so the same object could yield 22 or 23 segments from one
  call to the next. The caller's RNG state is restored afterwards.
