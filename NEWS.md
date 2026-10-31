# MaMaMIA 0.1.0

## Breaking changes

* **One pipeline.** `analyzeBamInput()` now runs the breadth-based pipeline. The
  count -> ZINB -> dose -> difference -> CBS path and its `pipeline` switch are
  gone.
* **Configuration version 2 is required.** `version: 1` is rejected with a
  migration error; see `MIGRATION.md`. The `max_fit_windows` and
  `normalization_by_side` defaults were removed.
* **New result contract.** `BamAnalysisResult` now exposes `metrics`, `scores`,
  `masks`, `calibration` and `calls`. The `counts`, `rca`, `isa` and
  `normalization` fields were removed.
* **Call classes are multi-axis.** A recipient loss without donor support is
  now `deletion`; donor gain with recipient loss and a neutral total is
  `substitution`; donor gain with a total gain is `addition`; a total imbalance
  or a width mismatch is `complex`. Donor and recipient cutoffs are calibrated
  from control segments. Calls add `class_evidence`, `class_confidence`,
  `class_status` (`confident`, `ambiguous`, `needs_review`), `copy_state`
  (`balanced`, `gain`, `loss`, `compound`, `neutral`), `total_evidence`,
  `component_status` and `mappability_status`. Deletion-only calls previously
  labelled `substitution` change class.
* **Assignment is calibrated.** Calling positions now use per-contig corrected
  recipient residuals, standard-error-aware depletion scores, a width-matched
  permutation null, contig-fraction reciprocal overlap, and weighted interval
  scheduling for multiple runs on one contig. Masked windows are excluded.
  Cross-group hosts are reported as `cross_group`.
* **Relative total depth.** The score table adds `total_residual`; calls add
  `total_mean`. This separates balanced substitutions from additions,
  deletions, amplifications and compound events without an external
  mappability track.
* **Dosage covers additions and deletions.** Donor-only additions use
  `donor_copies`; recipient losses without donor support use the loss fraction
  at lower confidence. `affected_copies` is no longer limited to donor-backed
  substitutions.
* **MAPQ proxy.** `bam_metrics_rg_cpp()` also returns per-window high-MAPQ read
  counts. `mapq_high_fraction` feeds `mappability_status`; low-MAPQ calls are
  downgraded to `ambiguous`.
* **Components replace genome-type labels.** Every normalized contig carries
  an arbitrary `component` label. The label is used only to aggregate contigs
  and copies; it carries no ploidy or genome-type assumption, and donor and
  recipient labels may overlap. Configuration and result columns use
  `component`, `recipient_component` and `donor_component`. A
  single-component species may omit `component`; it defaults to the species
  `default_component`, or to `"default"` when none is declared. A truly
  diploid species can therefore declare only `id` and `name` per contig.
* **Removed API.** `RCA()`, `correctReadCounts()`, `normalizeReadCounts()`,
  `segments()`, `buildSampleTables()`, `writeCohortToBedpe()` and
  `writeToBedpe()` were removed, together with the `RCA`/`ISA` S3 methods.
  Use `writeCallsToBedpe()` for export.
* **Removed dependencies.** `glmmTMB`, `mgcv`, `dplyr`, `purrr`,
  `ggnewscale`, `htmlwidgets` and `interacCircos` are no longer dependencies.
* **Dashboard.** The RCA/ISA circos and mirror/difference view were replaced by
  a score profile and call table. `exploreIntrogressions()` now accepts a
  `BamAnalysisResult`.

## New method

* **Breadth and counts in one pass.** `bam_metrics_rg_cpp()` returns both the
  number of alignments overlapping each exact-N window and the number of
  covered bases. Breadth is `covered / length`; count depth is
  `count * read_length / length`. No covidx index or binary is used.
* **Ploidy-agnostic scoring.** Assemblies select one or more recipient and
  donor components. The paired donor/recipient log-ratio averages recipient
  homoeologs and supports `mean`, `max` and `sum` donor aggregation. Donor
  aggregation now defaults to `max`: at each window the donor component with
  the strongest donor/recipient log-ratio is used, so a single-component
  signal is not diluted by the other donor homoeologs. The registry supplies
  the default (`defaults.donor_aggregation`); the dashboard exposes the
  setting in the **Advanced** sidebar panel, where calibration controls can
  also be selected. An empty control selection pools all cohort samples. A
  donor-only fallback is available for CNV-rich material.
* **GC correction.** A single sample uses a natural spline and centring; a
  cohort uses median polish. The vendored gcrank library remains the only GC
  backend.
* **Window masks.** Panel, cohort recurrence, BED blacklist and physical
  terminal fallback masks are applied before CBS.
* **Length-stratified calibration.** Segment thresholds are the empirical
  `1 - fdr` quantile of control segments per length bin, raised to a floor.
* **Recipient-copy dosage.** Assigned calls estimate how many recipient
  copies were replaced. Breadth is converted to a per-copy depth scale with
  `-log(1 - breadth)`; a call with donor support reports
  `affected_fraction`, `affected_copies`, `recipient_observed_breadth`,
  `recipient_copy_number` and `dosage_confidence`. Copy-neutral additions,
  unassigned calls and deletion-only signals have no dosage estimate.
* **Two-sided classification.** Calls now use donor gain and recipient loss
  together: both = `substitution`, donor only = `addition`, recipient only =
  `deletion`, neither = `other`. `class_evidence` records the supporting
  side(s), and `class_confidence` is a separate class-level confidence. A
  deletion-only signal is no longer promoted to substitution.
* **BEDPE export.** `writeCallsToBedpe()` writes donor and recipient intervals
  with `score_mean`, `class`, `class_group`, `class_evidence`,
  `class_confidence`, `confidence`, `method` and `mode`.

## Dashboard

* The sidebar accepts exactly one BAM per run, single- or multi-sample. The
  file picker has no `multiple` attribute, a drag-and-drop guard keeps only
  the first dropped file, and the server rejects a multi-file payload with a
  message instead of silently ignoring the extra files. A multi-sample BAM
  shows the **Sample** selector in the result filter bar as soon as its
  `@RG SM` samples pass preflight, before the analysis runs, and keeps it
  available afterwards to drive the karyotype, call table and exports.
* Component selectors are gone. The app resolves every component declared for
  the selected species — for example A+B+D versus S — and analyses them in one
  assembly, with a read-only summary under **Assembly**.
* The resolution picker is a slider from the smallest to the largest
  `allowed_windows` value, in 10-window steps, with labelled ticks at those
  allowed values. An **Exact** number field beside the **Resolution** label
  accepts a typed value in the same range and step; the slider and field stay
  in sync. The default range is 100–5,000 with ticks at 100, 1,000 and 5,000.
  In the bundled example all analysis parameters are fixed: the resolution
  slider and exact field are locked and shown as inactive at the precomputed
  1,000-window value, while the recipient species, donor species and donor
  aggregation are read-only. Clicking the slider explains that a BAM upload
  is needed to change resolution.
* The main tab is now a vertical, interactive karyotype built with
  `htmlwidgets`. Chromosomes are shown in the order they appear in the
  assembly config, with no group or component columns. Recipient and donor
  chromosomes are drawn as diploid homolog pairs: the two copies have a small
  gap inside the pair and a larger gap before the next chromosome. Both copies
  carry the same label (`1A`, `1B`, …), and each bar is scaled by its
  configured length. A deleted recipient chromosome stays as a grey dashed
  bar, and donor chromosomes are added only when a passing call covers at
  least 95% of the group. Added donor chromosomes sit next to their recipient
  group rather than at the end of the layout, and a full substitution is drawn
  as the lost recipient plus the added donor. A full substitution marks only
  the recipient chromosome it is assigned to; in a multi-contig group such as
  wheat 7A/7B/7D it no longer dashes every homoeolog in the group. A toggle
  switches between the introgression-score
  heatmap and the coloured call map; clicking a chromosome filters its group
  and clicking the background clears the filter. The Score profile tab was
  removed. The chromosome-group picker is removed from both the BAM dashboard
  and the standalone explorer; group filtering is driven by the karyotype
  click.
* The **Karyotype** tab adds a **Tooltip help** button beside the
  Heatmap/Calls switch. It opens a viewport-fixed, scrollable reference with
  quick interpretation patterns, placement and evidence guidance, chip
  definitions, dosage caveats, geometry checks, no-call behavior and caveats.
* The dashboard has a new visual system: bundled Inter type, warm paper
  surfaces, a deep-teal accent and a single white result card with underline
  tabs. Result panes are inset to align with the tab labels. The header adds
  a chromosome-crossover mark — the two homolog strokes meet at a central
  chiasma — and the package version, the BAM chooser is a
  click-or-drop zone, **Run analysis** is the only filled action, and an
  empty result shows a centred illustrated state with a short hint instead
  of a bare message. Karyotype labels use larger Inter type, title-case
  legend entries and captions, and a legend that wraps on narrow cards.
  Analysis progress and toast notifications use compact spacing in the
  bottom-right corner instead of bslib's large default margins and padding.
  The progress card shows a short stage label such as **Counting reads** or
  **Calibrating calls**, without counts or sample names.
* The **Export** tab no longer repeats the near-zero coverage warnings; they
  remain in the calculation log.
* The karyotype legend labels the dashed deleted-chromosome state **Deleted
  recipient**; the caption above the chromosome remains **Deleted**.
* The **Karyotype** toolbar downloads the current filtered view or a full
  view with every call for the sample as an SVG.
* The sidebar species dropdowns show the registry labels without the
  internal species keys.
* The Assembly panel no longer repeats the component description.
* The **Calls** and **Export** tabs stay disabled until an analysis result
  is available.
* While an analysis runs, the main plot area shows an animated chromosome
  loader with the current stage instead of a progress toast.
* A BAM whose header does not match a configured assembly now shows a calm
  amber recovery card in the sidebar instead of a red error line.
* The call table drops the inert per-column "All" filter row that server-side
  DT renders disabled. The global search and the pill filters remain.
* The result filter bar is now a compact row of pill dropdowns for class,
  class status, copy state and confidence, plus a **Passing only** switch.
  Empty selections mean "all", so the bar no longer shows four inline
  checkbox groups.
* The bundled example now uses eight chromosome pairs of 160–260 kb. It
  demonstrates every legend category: segment substitution, addition,
  deletion, complex and other calls, a full-chromosome substitution, a whole
  added donor chromosome and a whole deleted recipient chromosome. All five
  copy states and all three class statuses remain visible. The example
  resolves at 1,000 windows per chromosome.
* The karyotype tooltip is now a structured call card. It has a class/status
  header, a placement block (recipient, donor, relation, assignment and group
  class), an evidence block (score, donor, recipient and total residuals plus
  copy/evidence/quality chips), a dosage block, a geometry block and a
  group-filter hint. When no call is under the pointer it falls back to a
  chromosome card. The old flat list of lines is gone.
* The karyotype now performs an all-vs-all assignment of every donor segment
  against every recipient contig. Detection stays group-level, but candidate
  scoring uses a per-contig recipient residual, a standard-error-aware
  depletion score, a width-matched permutation null and a same-group homology
  bonus. Candidates are compared in contig-fraction space with reciprocal
  overlap. A per-sample assignment chooses at most one call per recipient
  contig plus a dummy "no host" option; weighted interval scheduling places
  additional non-overlapping runs on one contig. Masked windows are excluded
  from the baseline and the run. Each call records `recipient_group`,
  `recipient_cross_group`, `recipient_expected_breadth`,
  `recipient_deletion_width_bp`, `recipient_weight`, `recipient_margin`,
  `recipient_depletion_se`, `recipient_width_overlap` and an
  `assigned`/`cross_group`/`ambiguous`/`no_recipient_loss` flag. A strong
  assigned recipient loss reclassifies the call as a substitution, while
  `class_group` keeps the original group-level class. Cross-group hosts are
  reported as `cross_group` translocation candidates, not confirmed hosts.
  The call table and BEDPE export carry the assignment; calls with no host are
  exported with `NA` recipient coordinates.
* Passing calls are matched before non-passing calls, so a non-passing
  segment cannot take the best host and push a passing call onto a worse
  one. A cross-group host is never marked `assigned`; without a reliable
  same-group candidate it is `cross_group`, and with a competing same-group
  candidate it stays `ambiguous`. This keeps a copy-neutral event from being
  confirmed on an unrelated deletion.

## Packaging

* The chromosome-crossover mark now ships as `inst/www/mamamia-logo.svg` and
  is served as the browser favicon. The package, pipeline and widget-asset
  versions are reset to 0.1.0. The version is pinned from here and changes
  only on an explicit request.

## Validation

* The paired log-ratio recovers **6/6** reference events on the 24
  *Aegilops speltoides* BAMs (Sq74 7S; Sq75 7S/2S/6S/5S; Sq76 7S), with one
  control call >= 20 Mb.
* Global assignment places five of the six reference events on recipient
  contigs: Sq74 7D and Sq75 2B/5B/6D/7D. The Sq76 7A event is detected as an
  addition without recipient loss and receives no host. Across synthetic
  copy-neutral sample/group events (donor present, no loss in the call's own
  group), none is marked `assigned`: events receive no host or remain
  ambiguous when the background sample carries a natural deletion elsewhere.
  Synthetic copy-neutral events retain the group-level class.
* The donor-only fallback also recovers 6/6, with two control calls >= 20 Mb.
* Single-sample analysis recovers the same six events, including the distal
  Sq75 5S and Sq76 5S segments.
* The 24-BAM run, including global assignment, takes about 450 s with 30
  cores; one 28,000-window BAM takes about 19 s on one core.
* The 24-sample panel mask masks 16.3% of windows, masks terminal artifacts and
  shortens (does not delete) the Sq76 5S event.
* The five-sample UI fixture exercises same-group assignment, a cross-group
  translocation candidate and a copy-neutral no-host call. Cross-group hosts
  are reported as `cross_group`, independent of the older within-group
  placement rule.

# MaMaMIA 0.99.1

## Unified BAM counting

* `analyzeBamInput()` now counts every BAM with the same htslib backend,
  `.count_reads_by_sample()`: a one-sample BAM counts every read, and a BAM
  with several `@RG SM` samples is split by `RG` tag. The single-sample
  bamsignals branch is gone from the pipeline, so both cases share one flag,
  MAPQ and window-boundary rule. `countReads()` now uses a C++ range counter
  with the same htslib rules, and `bamsignals` is removed from the package
  dependencies entirely.
* Fixed the htslib window assignment at exact window boundaries. The counter
  now uses the same floor boundaries as `makeAnalysisWindows()`, so a read on
  the first base of a window is no longer assigned to the previous window when
  the contig length is not divisible by the window count. Provenance records
  `count_backend = "htslib"` and `count_semantics = "five_prime_start_based"`.

## Vendored gcrank 0.3.1

* Updated the vendored gcrank sources from 0.1.0 to 0.3.1 (upstream commit
  `aa3a32b`). The on-disk `format_version` stays `0.1.0`, so existing indexes
  remain valid, while new indexes record `builder_version = "0.3.1"` and
  provenance records the new commit.
* gcrank 0.3.0 adds `GcRank::stats()` and speeds up index builds (SHA-NI
  hashing and hashing overlapped with parsing) and window queries; 0.3.1 fixes
  SHA-NI detection on Clang 18. MaMaMIA's wrapper keeps using
  `gc_fraction()` and `contigs()`, so the computed GC content is unchanged.

## Grid dashboard layout

* The BAM-first dashboard now uses the `gridlayout` package with a `0.5rem`
  gap, `240px` sidebar, `75px` header, fluid main area and `100px` footer. The
  slimmer header and footer give the result tabs more vertical space; the
  application name is shown in large bold type in the header. The main area
  has a shared filter bar above the **Introgression circos**, **Coverage /
  dosage**, **Segments** and **Export** tabs.
* The sidebar starts with one BAM upload control and an **Example data** button.
  Analysis controls appear after preflight succeeds and live in single-open
  collapsible panels (**Analysis**, **Assembly**). **Run analysis** sits below
  both panels and is disabled until preflight and the selected resolution are
  ready. There is no manual single/multi mode switch: the sample count is
  detected from the BAM header. The **Advanced** panel exposes donor
  aggregation and calibration-control selection; MAPQ, flag masks, fit limits
  and cores come from the validated registry and the BAM is processed as-is.
* The **Filters** and **Export** panels were removed from the sidebar. The main
  area now keeps the chromosome, chromosome-pair and putative-only filters
  above every result tab; call threshold, minimum width and row selection live
  in **Segments**; the cohort summary and BEDPE downloads live in **Export**.
* Multi-sample BAMs now show a **Sample** selector in the main filter bar. It
  switches the circos, coverage plot, segment table and per-sample export
  together; single-sample runs do not show it.
* The first run now precomputes every valid component combination for the
  selected species pair and caches each completed analysis. Returning to a
  component that was analysed earlier restores its result without recalculating
  every sample; combinations that cannot be resolved or fail are logged and
  skipped. The first run takes longer, and changing the resolution or reference
  still needs one new run.
* The circos view is zoomable (scroll to zoom, drag to pan), uses the
  chromosome names from the registry instead of contig IDs, and hides hover
  tooltips when the pointer leaves the plot or a stale widget render leaves a
  tooltip behind.
* The footer calculation log can be collapsed with **Hide log**. The footer
  shrinks to its status line, and the sidebar, circos and coverage plot take the
  freed height; **Show log** restores the log. The log names detected and
  precomputed assemblies by their human-readable labels instead of internal
  registry keys such as `triticum_aestivum__A__secale_cereale__R`.
* A single BAM may contain one or several samples. Samples are read from
  `@RG SM` tags and read groups from `@RG ID`; incomplete or duplicate read
  groups are rejected, and reads are counted per read group for multi-sample
  BAMs. Multi-file upload and the old multiple-BAM mode were removed; the
  single/multi distinction is now inferred from the header.
* The **Assembly** panel contains only the recipient/donor species and
  component selectors plus **Swap donor and recipient**. Assembly detection,
  preflight tables and candidates are no longer shown there.
* Preflight status, detected assemblies, assembly candidates, problems,
  suggested fixes and every pipeline stage are written to the calculation log
  in the footer. The line above the log reports only the number of BAM files.
  The log autoscrolls to new entries. The standalone preflight, problems and
  detected-window cards are no longer part of the main area.
* The bundled example is now loaded with the **Example data** button; pressing
  it again clears the example and restores the registry resolution choices.
  The button is available in both modes and fills the example species and
  resolution fields. Example progress messages are identical to a real BAM
  run; the example itself remains pre-computed and needs no mounted BAM.

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
* The image now installs `pandoc` so the rmarkdown/knitr development
  dependencies and vignette tooling work in the container. The lockfile is
  installed with the normal binary `uvr sync --frozen`; no source-build flag is
  required. On a host that cannot reach `bioconductor.org`, build with
  `--build-arg BIOC_MIRROR=https://bioconductor.statistik.tu-dortmund.de/packages`
  and `--network=host` as documented in the Dockerfile.

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
  `default_component`); it is not tied to the donor or recipient role. A
  species-only config needs at least two species, and `assemblies` may be
  omitted entirely.
* Added `listMaMaMIASpecies()` and `resolveSpeciesAssembly()`. Resolution
  selects one component per side and pairs contigs by chromosome group (the
  optional `group` field, otherwise the first number in `name`), includes only
  the paired contigs, and returns the derived assembly with the configuration
  that contains it. `unpaired` explains every excluded contig. An explicit
  `pairs` argument remains available as an override.
* Named assemblies may now reference species (`recipient: wheat`, `donor: rye`),
  including optional `recipient_component` / `donor_component` and `pairs`.
  Inline `contigs:` assemblies keep working in the scripted API; the app falls
  back to BAM-header auto-detection when a registry has no species section.
* The BAM-first sidebar is species-only: choose a recipient species and a
  donor species (in either direction), with component selectors only when a
  species has several components. A "Swap donor and recipient" button exchanges
  the two species and carries their component choices with them. The resolved
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
  windows). Zero/near-zero coverage warnings are written to the calculation log
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
  stratified by component and GC decile, with at least ten windows per
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
  recipient chromosomes by chromosome number within selected components.
* The explorer gained call threshold and minimum segment width controls that
  re-flag the segments through `subset.ISA()` without re-running CBS, plus
  select-all and clear buttons for the BEDPE export.

## Performance

* `correctReadCounts()` no longer calls `glmmTMB::predict()` twice on the full
  window table. The fitted ZINB response is evaluated once for the distinct
  (GC, component) pairs, with design matrices built from the stored smooth rather
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
