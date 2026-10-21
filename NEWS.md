# MaMaMIA 0.99.0

* Initial submission.

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

* The BioCircos-based circos view is temporarily disabled and the BioCircos
  dependency removed: the dashboard shows the segment table with filtering and
  the BEDPE export. The circos code (in `R/circos.R`, its test, and the
  `inst/mamamia/` assets) is commented out and can be reinstated together with
  the `BioCircos` entry in `Suggests`.
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
