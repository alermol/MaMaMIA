# MaMaMIA BAM-first API contract

**Status:** frozen for implementation
**Frozen on:** 2026-10-20
**Source of truth:** `AGENTIC_PLAN.md` sections 3.3-4
**Applies from:** Wave 0 onward

This note is the contract between the BAM-first implementation tasks. When a
task changes a signature, a data shape or an invariant, it must update this
note in the same change. The plan's decision table (D1-D14) is authoritative;
this note only makes it operational.

---

## 1. Data contracts

All window tables use **0-based, half-open** coordinates
(`chr_id`, `start`, `stop`). `GenomicRanges` objects use **1-based, inclusive**
coordinates. A table and its `GenomicRanges` representation describe the same
windows.

```r
# Parsed YAML configuration; returned by readMaMaMIAConfig()
MaMaMIAConfig <- list(
  version = 1L,
  defaults = list(
    windows_per_chromosome = 1000L,
    allowed_windows = c(100L, 1000L, 10000L),
    mapq = 20L,
    exclude_flags = 2308L,
    max_fit_windows = 50000L,
    gc_mmap = TRUE,
    gc_validate_crc = FALSE,
    normalization_by_side = TRUE
  ),
  assemblies = named list of AssemblyConfig,
  server = list(bam_root = NULL, max_bams = 20L, max_total_bytes = 5 * 1024^3)
)

# AssemblyConfig; contig lengths come from the gcrank index, not from YAML
AssemblyConfig <- list(
  key = "cs_t2t_se_cereale",
  label = "Triticum aestivum CS T2T + Secale cereale",
  description = "Competitive mapping to wheat CS T2T plus rye",
  gc_indexes = c("/srv/.../wheat/index.gcrank", "/srv/.../rye/index.gcrank"),
  contigs = list(
    recipient = data.frame(id, name, subgenome, length, index_path),
    donor     = data.frame(id, name, subgenome, length, index_path)
  ),
  pairs = c("OZ284382.1" = "CP155614.1")
)

# BAM header; returned by bamHeaderInfo()
BamHeaderInfo <- list(
  path = "/abs/path/sample.bam",
  format = "BAM",
  targets = data.frame(name, length),
  samples = c("Sq54"),
  read_groups = data.frame(id, sample),
  programs = data.frame(id, name, version, command),
  header_text = "@HD..."
)

# Generated windows; returned by makeAnalysisWindows()
AnalysisWindows <- list(
  gr = GenomicRanges::GRanges,   # 1-based inclusive
  table = data.frame(chr_id, start, stop, role, chr_name, subgenome), # 0-based
  n_per_chromosome = 1000L
)

# Per-sample tables; returned by buildSampleTables()
SampleTables <- list(
  dcov, dgc, rcov, rgc, meta,
  don_chr_ids, rec_chr_ids,
  sample, reference, window,
  target_pairs
)

# Full BAM analysis result; returned by analyzeBamInput()
BamAnalysisResult <- list(
  assembly = AssemblyConfig,
  config = MaMaMIAConfig,
  bams = data.frame(path, sample, index_path, n_targets),
  windows = AnalysisWindows,
  gc = data.frame(chr_id, start, stop, gc),
  counts = list(donor = matrix, recipient = matrix),
  rca = named list of RCA,
  isa = named list of ISA,
  normalization = list(mode, reference, by_side, pseudocount),
  provenance = list(...)
)
```

### 1.1 Example `SampleTables`

```r
# One sample, one donor contig, one recipient contig, 100 windows each.
# dcov/rcov are half-open coverage tables; dgc/rgc are GC fractions.
SampleTables <- list(
  dcov = data.frame(
    chr_id = rep("d1", 100L),
    start = seq(0, by = 100, length.out = 100L),
    stop = seq(100, by = 100, length.out = 100L),
    cov = rep(8, 100L)
  ),
  dgc = data.frame(
    chr_id = rep("d1", 100L),
    start = seq(0, by = 100, length.out = 100L),
    stop = seq(100, by = 100, length.out = 100L),
    gc = rep(0.40, 100L)
  ),
  rcov = data.frame(
    chr_id = rep("r1", 100L),
    start = seq(0, by = 100, length.out = 100L),
    stop = seq(100, by = 100, length.out = 100L),
    cov = rep(10, 100L)
  ),
  rgc = data.frame(
    chr_id = rep("r1", 100L),
    start = seq(0, by = 100, length.out = 100L),
    stop = seq(100, by = 100, length.out = 100L),
    gc = rep(0.60, 100L)
  ),
  meta = data.frame(
    chr_id = c("d1", "r1"),
    chr_name = c("1D", "1R"),
    subgenome = c("D", "R")
  ),
  don_chr_ids = "d1",
  rec_chr_ids = "r1",
  sample = "sample1",
  reference = "tiny",
  window = 100L,
  target_pairs = c(d1 = "r1")
)
```

`SampleTables` is passed to the existing `RCA()` constructor as
`do.call(RCA, SampleTables[c("dcov", "dgc", "rcov", "rgc", "meta")])`.
The extra fields are provenance for the pipeline and the UI.

### 1.2 Example `BamAnalysisResult`

```r
BamAnalysisResult <- list(
  assembly = config$assemblies$tiny,
  config = config,
  bams = data.frame(
    path = "/abs/tiny.bam",
    sample = "tiny-sample",
    index_path = "/abs/tiny.bam.csi",
    n_targets = 2L
  ),
  windows = list(
    gr = GenomicRanges::GRanges(...),
    table = data.frame(chr_id, start, stop, role, chr_name, subgenome),
    n_per_chromosome = 100L
  ),
  gc = data.frame(chr_id, start, stop, gc),
  counts = list(donor = matrix(0, 100L, 1L), recipient = matrix(0, 100L, 1L)),
  rca = list(`tiny-sample` = rca),
  isa = list(`tiny-sample` = isa),
  normalization = list(
    mode = "median", reference = NULL, by_side = TRUE, pseudocount = 0.5
  ),
  provenance = list(
    package_version = "0.99.0",
    htslib_version = "1.24",
    gcrank_format_version = "0.1.0",
    sample_source = c(`tiny-sample` = "rg"),
    config_path = "/abs/assemblies.yml",
    window_rule = "exact_equal_windows"
  )
)
```

---

## 2. Public function signatures

```r
# Configuration (implemented in R/config.R)
readMaMaMIAConfig(path)
validateMaMaMIAConfig(config)
listMaMaMIAAssemblies(config)

# BAM preflight
bamHeaderInfo(path)
ensureBamIndex(path, min_shift = 14L, threads = 1L)
matchAssembly(header, config, assembly = NULL)
preflightBams(paths, config, assembly = NULL)

# Windows and GC
readContigsTsv(path)
makeAnalysisWindows(assembly, contigs, n_per_chromosome = 1000L)
gcrankContigs(index_dir)
gcrankGc(index_dir, contig, start, end,
         remove_ambiguous = FALSE, mmap = TRUE, validate_crc = FALSE)
gcrankClearCache()

# Counts and pipeline
countReads(bam_path, gr, mapq = 20L, exclude_flags = 2308L, verbose = FALSE)
buildSampleTables(counts, gc, windows, assembly, sample)
analyzeBamInput(bam_paths, config,
                assembly = NULL,
                windows_per_chromosome = 1000L,
                reference_sample = NULL,
                normalization = c("median", "reference", "none"),
                mapq = 20L, exclude_flags = 2308L,
                cores = 1L, max_fit_windows = 50000L,
                cache = NULL, progress = NULL)

# Normalization and export
normalizeReadCounts(rcas, reference = NULL, by_side = TRUE, pseudocount = 0.5)
writeCohortToBedpe(result, file, samples = NULL, selected = NULL)
```

`readMaMaMIAConfig()` and `validateMaMaMIAConfig()` return the normalized
configuration invisibly; `listMaMaMIAAssemblies()` returns one row per assembly
with `key`, `label`, `description`, `n_donor`, `n_recipient`, `n_pairs` and
`n_gc_indexes`.

---

## 3. Invariants

1. **Coordinates.** Tables are 0-based half-open; `GenomicRanges` are 1-based
   inclusive. Convert once at the boundary, never mix conventions inside a
   single object.
2. **Equal windows.** Every chromosome in a sample has exactly
   `n_per_chromosome` windows. Contiguous half-open tiling is required:
   `start[i + 1] == stop[i]` and the final window ends at the contig length.
3. **Disjoint subgenomes.** Donor and recipient contig ID sets and subgenome
   label sets are disjoint. `RCA()` enforces this; generated tables must
   satisfy it rather than bypass it.
4. **One BAM = one sample.** The sample is the single distinct non-empty
   `@RG SM` value; if `@RG`/`SM` is absent, it is `basename(path)` without
   `.bam`. More than one distinct `SM` stops preflight. The resolution source
   (`"rg"` or `"basename"`) is recorded in provenance.
5. **Index.** Every BAM is read through a CSI index with `min_shift = 14`. A
   BAI is never used. The index is written next to the BAM when that directory
   is writable; otherwise it is staged into the session temporary directory.
6. **Normalization order.** `correctReadCounts()` writes `cor.gc`. When
   normalization applies, it reads corrected coverage, writes `dose`, and CBS
   consumes `dose`. `cor.gc` is never overwritten.
7. **Segmentation target.** Pairwise CBS uses `dose.don - dose.rec` when
   `dose` exists, otherwise `cor.gc.don - cor.gc.rec`. `ISA$meta$value_col`
   records the resolved column.
8. **GC fit.** At most `max_fit_windows` ideal windows are fitted, stratified
   by subgenome and GC; prediction covers every window. `max_fit_windows` is
   part of the fit-cache key and the result metadata.
9. **Configuration.** Index paths are resolved relative to the config file.
   Every configured target contig appears exactly once across the assembly's
   GC indexes. Only config-defined index paths are read.
10. **Result containers.** `BamAnalysisResult$rca` and `$isa` are named lists
    keyed by sample; every `ISA` is pairwise/two-ended. The cohort BEDPE has a
    mandatory `sample` column.

---

## 4. Error contract

All MaMaMIA errors are R conditions with class `MaMaMIAError` plus a
subclass. The message is user-facing and names the object that failed (BAM,
assembly, contig, file). Functions must stop before expensive work starts.

| Class | Raised by | Message shape |
|---|---|---|
| `MaMaMIAConfigError` | config loader/validator | `Assembly '<key>': ... '<contig>' ... '<path>'` |
| `MaMaMIABamHeaderError` | `bamHeaderInfo()` | `Could not open/read BAM '<path>': <htslib reason>` |
| `MaMaMIAUnsupportedFormatError` | `bamHeaderInfo()` | `Unsupported input format '<fmt>' for '<path>'; only coordinate-sorted BAM files are accepted` |
| `MaMaMIAAssemblyMatchError` | `matchAssembly()` / `preflightBams()` | `BAM '<path>' does not match any configured assembly` or `BAM '<path>' matches more than one assembly: <keys>` |
| `MaMaMIAIndexError` | `ensureBamIndex()` | `Could not build a CSI index for '<path>': <reason>` |
| `MaMaMIAIndexError` | `ensureBamIndex()` | `BAM '<path>' has a corrupt or unreadable CSI index: <path>` |
| `MaMaMIASampleError` | preflight | `BAM '<path>' has more than one distinct @RG SM value: <values>` |
| `MaMaMIASampleError` | preflight | `Duplicate sample name '<sample>' in input BAMs: <paths>` |
| `MaMaMIAHeaderMismatchError` | preflight | `BAM '<path>' target '<contig>' has length <observed>; assembly '<key>' expects <expected>` |
| `MaMaMIAPathError` | preflight | `BAM path '<path>' is outside server.bam_root '<root>'` |
| `MaMaMIALimitError` | preflight | `Too many BAM files: <n> (limit <max_bams>)` or `Total BAM size <bytes> exceeds <max_total_bytes> bytes` |
| `MaMaMIAAssemblyMatchError` | preflight | `BAM '<path>' covers only one side of assembly '<key>'` |

Config validation messages implemented in Wave 0 include, among others:

```text
`version` is required in the MaMaMIA config
Unsupported MaMaMIA config version <v> (supported: 1)
`assemblies` must contain at least one assembly
Assembly '<key>': gc_indexes must not be empty
Assembly '<key>': gc_index '<path>' does not exist or is not a directory
Assembly '<key>': gc_index '<path>' has no contigs.tsv (<file>)
Assembly '<key>': gc_index '<path>' has no manifest.json (<file>)
Assembly '<key>': gc_index '<path>' uses unsupported gcrank format_version '<v>' (supported: 0.1.0)
Assembly '<key>': contig '<id>' appears in more than one gc_index (<paths>)
Assembly '<key>': contig '<id>' is <bp> bp, shorter than the largest allowed window (<bp> bp)
Assembly '<key>': contig '<id>' is <bp> bp; 32-bit positions require lengths below 2^31
Assembly '<key>': donor and recipient subgenome labels must be disjoint; shared: <labels>
Assembly '<key>': pairs refers to unknown donor contig(s): <ids>
Assembly '<key>': pairs is not one-to-one; a recipient is used more than once
```

---

## 5. Frozen decisions (D1-D14)

| ID | Decision |
|---|---|
| D1 | Pairwise donor-recipient difference per sample. CBS target is `dose.don - dose.rec` when `dose` exists, otherwise `cor.gc.don - cor.gc.rec`. `ISA`, plots and BEDPE stay pairwise/two-ended. |
| D2 | One unique `@RG SM`; otherwise basename without `.bam`. Multiple distinct `SM` values stop preflight. |
| D3 | Always CSI, `min_shift = 14`. No BAI code path. Index in place when writable, otherwise session temp. |
| D4 | `gc_indexes` is a list; every target contig appears exactly once across indexes with consistent lengths; queries are grouped by index. |
| D5 | Reference: `(C_s + eps) / (C_ref + eps)`; median: `(C_s + eps) / (median_s(C_s) + eps)`; per side by default. `eps = pseudocount` (default 0.5). Both numerator and denominator zero gives dose 1. |
| D6 | Normalize after GC correction, before CBS; store `RCA$data$dose`; preserve `cor.gc`. |
| D7 | `segments(value_col = NULL)` uses `dose` if present, else `cor.gc`; `ISA$meta$value_col` records the choice. |
| D8 | No Rsamtools. Use Rhtslib for headers/indexes, bamsignals for counts. |
| D9 | Fit at most `max_fit_windows` ideal windows, stratified by subgenome and GC; predict all; include the limit in the cache key and metadata. |
| D10 | Remove legacy directory-layout code after the T7.4 gate; keep the core analysis API. |
| D11 | Per-sample `ISA` list, sample selector in the UI, combined BEDPE download with a `sample` column. |
| D12 | Support uploads (session temp) and server paths (only under `server.bam_root`); never write under `bam_root`. |
| D13 | Use `uvr`, not `renv`; commit `uvr.toml`, `uvr.lock`, `.r-version`; ignore `.uvr/library/`; CI uses `uvr sync --frozen --unattended`. |
| D14 | R packages come only from uvr; conda `mamamia` is for non-R tools only. Production Docker uses apt/GCC/uvr, not the conda environment. |

Resolved questions are frozen as in the plan: pairwise CBS, window-wise ratio
with pseudocount, per-side normalization, stratified subsampling for high
resolution, one combined BEDPE with a `sample` column, and real GC indexes for
the current assemblies.

---

## 6. bamsignals flag semantics (Wave 0 finding)

`bamsignals::bamCount(filteredFlag = mask)` filters a read only when **all**
bits in `mask` are set. This differs from `samtools view -F mask`, which
excludes a read when **any** bit is set. The real Sq54 BAM happens to contain
no secondary, supplementary or unmapped alignments, so the original spike did
not expose the difference.

`countReads(exclude_flags = ...)` must keep samtools-compatible semantics.
T3.1 must therefore not pass a combined mask such as `2308` straight to
`filteredFlag`. Use `filteredFlag = -1` plus per-bit filtering (for example an
inclusion-exclusion decomposition over `4`, `256` and `2048`), or another
approach that produces the same counts. The tiny fixture deliberately includes
one secondary, one supplementary, one duplicate and one low-MAPQ read so the
wrapper can be tested at single-read resolution.

---

## 7. Wave 0 build decisions

- gcrank is **vendored**, not a git submodule. `tools/vendor-gcrank.R` copies
  `include/` and `src/` plus `LICENSE` and `CHANGELOG.md` into `src/gcrank/`,
  records the upstream commit in `tools/vendor-gcrank.txt`, and is documented
  in `inst/COPYRIGHTS`.
- `src/Makevars` uses `CXX_STD = CXX23`, `Rhtslib::pkgconfig()`, an explicit
  `OBJECTS` list and `-pthread`. The R environment is activated through uvr,
  so `R_LIBS_USER` must point at `.uvr/library` when installing or compiling.
- `DESCRIPTION`, `uvr.toml`, `uvr.lock` and `.r-version` are the dependency
  contract. `tests/testthat/test-build.R` fails when the runtime or development
  dependency sets drift.
- Tiny fixtures live in `tests/testthat/fixtures/` and are generated by
  `tools/make-test-fixtures.R`. They contain two 20 kb contigs (`d1` at 40% GC,
  `r1` at 60% GC), a tiny gcrank index, a coordinate-sorted CSI-indexed BAM,
  and an `assemblies.yml` that loads without external tools.

---

## 8. Wave 1 implementation notes

### 8.1 BAM headers

`bamHeaderInfo(path)` returns the `BamHeaderInfo` contract plus `sample`,
`sample_source` (`"rg"` or `"basename"`), `sort_order` and `n_targets`.
`@RG SM` values are trimmed; missing `SM` values are `NA_character_` in
`read_groups$sample` and are excluded from `samples`. A header with multiple
distinct `SM` values raises `MaMaMIASampleError`.

A truncated or non-BAM file whose format cannot be determined raises
`MaMaMIABamHeaderError`. A readable SAM/CRAM raises
`MaMaMIAUnsupportedFormatError`.

### 8.2 CSI indexes

`ensureBamIndex(path, min_shift = 14L, threads = 1L)` returns the `.csi` path.
The CSI is BGZF-compressed, so validity is checked by decompressing the first
four bytes and comparing them with `CSI\1`. A current valid index is reused;
a stale, empty or corrupt index is rebuilt when the BAM directory is writable.
Permission bits, not root's ability to write, decide writability.

`stageBam(path, session_dir = tempdir())` is internal (`@noRd`). It creates a
private subdirectory per BAM, symlinks the BAM when possible and copies it
otherwise, and returns the staged path **without** `normalizePath()`. Callers
pass that path to `ensureBamIndex()`; the CSI then lands next to the symlink in
session temp. No `.bai` is ever created or consulted.

### 8.3 Assembly matching

`matchAssembly()` returns:

```r
list(
  assembly = "tiny",                  # selected key
  assembly_config = AssemblyConfig,
  header_path = "/abs/sample.bam",
  candidates = data.frame(key, label, n_targets, n_matched, n_missing,
                          n_length_mismatches, score, exact),
  matched = TRUE, exact = TRUE,
  missing = character(0),
  length_mismatches = data.frame(chr_id, role, observed, expected),
  score = 1,
  summary = "BAM '...' matches assembly 'tiny' exactly (2/2 targets)"
)
```

An explicit `assembly` key is validated; a length mismatch raises
`MaMaMIAHeaderMismatchError`, a missing target raises
`MaMaMIAAssemblyMatchError` (with `covers only one side` when exactly one role
is complete). Auto-selection errors carry the candidate table in
`condition$data$candidates`. `matchAssemblyCohort()` is internal and enforces
one assembly key across several headers.

### 8.4 Windows

`readContigsTsv()` enforces exactly three tab-separated columns, unique names,
positive integer lengths, non-negative integer offsets and offsets equal to
the cumulative lengths of preceding contigs.

`.assembly_target_contigs()` (internal) flattens recipient contigs first and
then donor contigs, preserving assembly order, and validates a supplied
contigs table (a path or a data frame with `name`/`length` or `id`/`length`).
`makeAnalysisWindows()` uses `start_i = floor((i - 1) * L / n)` and
`stop_i = floor(i * L / n)`; the GRanges and the 0-based table therefore
describe the same windows. `n_per_chromosome` is recorded in the result.

**Known plan inconsistency (wave 1 finding).** The legacy
`*.w1000.*.bed` tables appear to use a fixed window width of `floor(L / n)`
with the last window stretched to `L`, whereas T2.2 specifies the distributed
remainder formula above. The two tilings differ after the first window (for
`L = 603818960` and `n = 1000`: legacy stops are `603818, 1207636, ...`, the
T2.2 stops are `603818, 1207637, ...`). Wave 1 implements the T2.2 formula.
T2.4/T3.4 must either compare GC/counts at common intervals rather than by row,
or the plan must be amended to legacy tiling before those tasks run.

### 8.5 gcrank

`gcrankContigs()`, `gcrankGc()` and `gcrankClearCache()` are exported and use
an internal LRU cache of at most two open `XPtr<gcrank::GcRank>` handles keyed
by normalized index directory. Queries are vectorized over `contig`, `start`
and `end` (each length 1 or the common length) and raise
`MaMaMIAGcrankError` with contig/interval context. Defaults are `mmap = TRUE`
and `validate_crc = FALSE`; resident memory grows as mapped pages are touched.

New error subclasses introduced in wave 1: `MaMaMIAWindowError` and
`MaMaMIAGcrankError` (both inherit `MaMaMIAError`).
