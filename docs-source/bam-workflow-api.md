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
The extra fields are provenance for the pipeline and the UI. `buildSampleTables()`
returns this list with class `c("SampleTables", "list")`.

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
  samples = data.frame(
    sample = "tiny-sample",
    path = "/abs/tiny.bam",
    assembly = "tiny",
    index_path = "/abs/tiny.bam.csi",
    reference = FALSE,
    normalization = "none",
    n_windows = 200L,
    n_segments = 1L,
    n_putative_introgressions = 0L
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
    mode = "none", reference = NULL, by_side = TRUE, pseudocount = 0.5
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

`analyzeBamInput()` returns this list with class `c("BamAnalysisResult", "list")`.

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
segments(RCA, target_pairs, value_col = NULL, ...)
writeToBedpe(ISA, file, sample = NULL)
writeCohortToBedpe(result, file, samples = NULL, selected = NULL)

# Application entrypoint
runMaMaMIA(analysis_dir = NULL, assemblies_dir = NULL, cache = NULL,
           config = NULL, bam_paths = NULL,
           input = c("auto", "bam", "files", "example"),
           max_request_size = 5 * 1024^3,
           host = "127.0.0.1", port = NULL, launch.browser = interactive())
```

`readMaMaMIAConfig()` and `validateMaMaMIAConfig()` return the normalized
configuration invisibly; `listMaMaMIAAssemblies()` returns one row per assembly
with `key`, `label`, `description`, `n_donor`, `n_recipient`, `n_pairs` and
`n_gc_indexes`.

`runMaMaMIA()` is the only application entrypoint. The legacy directory app is
still reachable through `input = "files"`/`"example"` until T7.5, but the BAM
workflow is selected whenever `config` or `bam_paths` is supplied (or when
`input = "bam"`). `bam_launcher_app()` and `mamamia_app()` are internal.

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
| `MaMaMIAPreflightError` | preflight | malformed `paths` argument |
| `MaMaMIACountError` | `countReads()` | `Could not count reads in BAM '<path>': <reason>` |
| `MaMaMIAAnalysisError` | `analyzeBamInput()` | invalid pipeline argument or unsupported normalization mode |
| `MaMaMIAAnalysisError` | `normalizeReadCounts()` | invalid sample list, incompatible windows, unknown reference or non-positive pseudocount |

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

---

## 9. Wave 2 implementation notes

### 9.1 Preflight

`preflightBams(paths, config, assembly = NULL)` reads file metadata and BAM
headers only; it never builds an index or reads alignments. It returns one row
per BAM with:

```text
path, file_size, sample, sample_source, assembly, n_targets,
n_assembly_targets, index_path, index_status, index_action, warnings, valid
```

`index_status` is one of `valid`, `missing`, `stale` or `corrupt`.
`index_action` says what the pipeline will do: `reuse`, `create`, `rebuild`, or
`stage` when the BAM directory is not writable. A read-only BAM with no valid
CSI is therefore not an error at preflight; it is planned for staging.

Limits come from `config$server` (`max_bams`, `max_total_bytes`). Samples follow
D2 and duplicate resolved names stop preflight. `matchAssemblyCohort()` enforces
one assembly across all inputs; extra header targets are allowed and reported
as a warning.

Path containment is enforced when `server.bam_root` is set: paths are resolved
with `normalizePath()` (so `..` and symlinks cannot escape), and a path is
accepted only under the root or under the session temporary directory. The
upload directory is `getOption("MaMaMIA.session_dir", tempdir())`; the app can
set this option per session. A symlink in session temp that resolves outside
both locations is rejected.

`.list_bam_files(root, max_depth, max_files)` (internal) is the bounded listing
for the server-side picker. It returns regular `.bam` files only, drops
symlinks that resolve outside `root`, and records `attr(x, "truncated")` when
`max_files` is reached.

New error subclasses in this wave: `MaMaMIAPreflightError`,
`MaMaMIAPathError`, `MaMaMIALimitError` and `MaMaMIACountError` (all inherit
`MaMaMIAError`).

### 9.2 GC extraction for all windows

`extractWindowGc(windows, assembly, mmap = TRUE, validate_crc = FALSE)` is
internal (`@noRd`). It accepts an `AnalysisWindows` list or a table with
`chr_id`, `start` and `stop`; builds the contig-to-index map from
`assembly$contigs$index_path`; calls `gcrankGc()` once per index; and returns
`data.frame(chr_id, start, stop, gc)` in input order. Coordinates are validated
against the assembly lengths, and GC values must be non-missing and in `[0, 1]`.

Wave 1's tiling finding still applies: `makeAnalysisWindows()` follows the
T2.2 distributed-remainder formula, while the legacy `*.w1000.*.bed` files use
a fixed width with a stretched last window. The gated real-data GC test
therefore compares `extractWindowGc()` with the legacy BED at the legacy
intervals (all 28,000 rows, six decimals). New windows are covered by the tiny
fixture and by the ordering/validation tests. This avoids changing the plan's
window definition without an explicit decision.

### 9.3 Read counting

`countReads(bam_path, gr, mapq = 20L, exclude_flags = 2308L, verbose = FALSE)`
returns an integer vector aligned to `gr`. It validates the GRanges, its
seqlengths and bounds, and ensures a CSI index with `ensureBamIndex()` before
counting.

`exclude_flags` keeps `samtools view -F` semantics: a read is excluded when
**any** mask bit is set. Because bamsignals excludes only when **all** bits are
set, the mask is decomposed with inclusion-exclusion over its set bits:

```text
pass = sum over non-empty subsets S of (-1)^(|S| + 1) * C(mask(S))
```

For the default `2308 = 4 + 256 + 2048` this is seven bamsignals passes. The
implementation supports at most six bits (63 passes) and raises
`MaMaMIACountError` for larger masks instead of silently distorting counts. The
Sq54 tiny fixture locks the difference down: a direct `filteredFlag = 2308`
call returns `c(5, 8)` in window order, while the correct count is `c(5, 6)`.

**Counting is 5'-end based.** bamsignals assigns a forward-strand read to the
window containing its leftmost coordinate and a reverse-strand read to the
window containing its rightmost coordinate (the alignment end). A reverse read
can therefore belong to the window after the one containing its SAM `POS`.
This differs from `bedtools coverage -counts`, which counts a read in every
window it overlaps. The difference is confined to reads whose 5' end is outside
the window while the alignment span overlaps it; it is not a pure
"alignment start" rule. `tests/testthat/test-counts.R` includes a reverse-strand
boundary test that locks this down.

The gated Sq54 regression counts the legacy BED intervals and matches the
legacy tables within the documented tolerance: 229/28,000 windows differ
(0.82%; all on the recipient side, donor 0%), maximum delta 14, total absolute
delta 708/3,535,298 (0.02%), all deltas non-positive. Every mismatch is
explained exactly by the number of overlapping reads whose 5' end is outside
the window (verified with `samtools view`, 229/229), and per-contig 5'-end
totals equal `samtools view -c` exactly.

Set `options(MaMaMIA.count_chunk_size = n)` to split very large window sets into
chunks. The default `0` disables chunking because every chunk repeats the
bamsignals passes.

---

## 10. Wave 3 implementation notes

### 10.1 Per-sample tables

`buildSampleTables(counts, gc, windows, assembly, sample)` returns a classed
`SampleTables` list with the documented names. `counts` and `gc` may be numeric
vectors aligned to `windows$table`, or a `donor`/`recipient` list; `gc` may
also be the data frame returned by `extractWindowGc()`. The function validates
that the values are aligned to the window table, that every contig has one
role/name/subgenome, that each chromosome has the same number of windows, and
that donor and recipient IDs are disjoint.

Pairs come from `assembly$pairs` when configured; otherwise
`defaultTargetPairs()` derives them from chromosome numbers. Explicit pairs
that name an unknown or missing contig, duplicate a donor, or reuse a
recipient raise `MaMaMIASampleError`. The assembled tables are passed through
`RCA()` and the result must be uncorrected, so malformed data fails before the
GC-bias fit. Validation errors are re-raised with the sample name.

### 10.2 Single-BAM orchestration

`analyzeBamInput()` runs preflight, index, windows, GC, counts, tables, RCA,
correction and segmentation in that order. It stages BAMs whose directory is
not writable before indexing, uses a per-process fit cache in session temp when
`cache = NULL`, and records provenance (config path/hash, BAM sizes/mtimes,
package/htslib/gcrank/uvr versions, `uvr.lock` hash, window rule, filters,
normalization and timings).

The `progress` callback is `function(stage, detail)` and is called after each
completed stage: `preflight`, `index`, `windows`, `gc`, `counts`,
`correction`, `segmentation`. For several samples the last three stages repeat
per sample. A callback error aborts the run (used for cancellation).

For one BAM the normalization mode is forced to `"none"` and CBS consumes
`cor.gc` (D7). For several BAMs, `"median"` (the default), `"reference"` and
`"none"` are supported. `"reference"` requires `reference_sample`; the
reference sample's dose is exactly 1. A sample with zero reads in every
analysis window stops before the fit with `MaMaMIASampleError`; a nearly empty
sample surfaces the correction-stage error with the sample name, never a raw
model failure.

`counts$donor` and `counts$recipient` are integer matrices with one column per
sample (columns named by sample). `rca` and `isa` are named lists keyed by
sample. `samples` is a data frame with one row per sample (`sample`, `path`,
`assembly`, `index_path`, `reference`, `normalization`, `n_windows`,
`n_segments`, `n_putative_introgressions`). The result is classed
`BamAnalysisResult`.

### 10.3 Golden regression

`tools/benchmark-bam-pipeline.R` and the opt-in
`tests/testthat/test-golden-bam.R` (`MAMAMIA_SLOW_TESTS=1`) compare Sq54 with
the legacy directory-layout tables. Results:

| Check | Result |
|---|---|
| counts on legacy intervals | 229/28,000 windows differ, max delta 14, total abs delta 708 (0.02%), all explained by 5' ends |
| GC on legacy intervals | identical at six decimals (28,000/28,000) |
| generated windows | exact-N tiling; per-contig 5'-end totals equal `samtools view -c` |
| full pipeline vs legacy segments | same 7 pairs; segment/putative count deltas 0; introgression bp relative delta 9.3e-7 |
| core on the same intervals | identical breakpoints and putative calls; segment medians within 0.13% of scale |

Both pipelines use the config's explicit pairs (seven rye donors to the seven
A-genome recipients); the legacy chromosome-number rule would leave these pairs
unmatched because every rye chromosome number matches three wheat chromosomes.

The legacy BEDs use fixed-width windows with a stretched last window, whereas
`makeAnalysisWindows()` uses the T2.2 distributed-remainder formula. The two
tilings differ by at most 999 bp (963 bp on this data) and only 126/28,000
intervals are identical, so per-window golden comparisons are made on the
legacy intervals and the generated windows are validated by exact-N tiling plus
per-contig totals. No production or test-only GC fallback is used; the golden
test skips cleanly when `test_data` is absent.

---

## 11. Wave 4 implementation notes

### 11.1 Normalization core

`normalizeReadCounts(rcas, reference = NULL, by_side = TRUE, pseudocount = 0.5)`
takes a named list of corrected `RCA` objects. Every sample must describe the
same `(chr_id, iid)` windows with the same `start`/`stop` coordinates and the
same chromosome sets; rows may be in any order because they are aligned by
`(chr_id, iid)`. Missing `cor.gc` values are allowed; infinite values are
treated as missing.

- Reference mode (`reference` is a sample name):
  `dose_s = (C_s + eps) / (C_ref + eps)`. The reference sample's dose is exactly
  1. A zero reference window with all-zero samples also gives 1; a zero
  reference window with positive samples uses the pseudocount ratio.
- Median mode (`reference = NULL`):
  `dose_s = (C_s + eps) / (median_s(C_s) + eps)`, with the median taken over all
  finite sample values of that window, not restricted to `ideal` windows.
- `by_side = TRUE` (default) computes the denominator separately for donor and
  recipient windows. `by_side = FALSE` is a sensitivity option: donor and
  recipient windows with the same position index (`iid`) are pooled into one
  per-position baseline. Joint mode is only available for median normalization;
  reference mode must keep the sides separate so that the reference dose stays
  exactly 1.

`cor.gc` is never modified. `dose` is added to `RCA$data`, `RCA$normalization`
records `mode`, `reference`, `by_side` and `pseudocount`, and
`ideal <- ideal & is.finite(dose)` so CBS and plots skip non-finite windows.
`pseudocount` must be positive; this guarantees finite ratios whenever the
input is finite.

### 11.2 Dose-aware segmentation and reversal

`segments(RCA, target_pairs, value_col = NULL, ...)` resolves `value_col` to
`"dose"` when the RCA has a dose column and to `"cor.gc"` otherwise (D7). The
`ISA` records the choice in `meta$value_col`; `out$diff` is the donor minus
recipient difference of that column. `segments()` output is identical to the
pre-wave behavior when `dose` is absent. `reverseWindows()` reverses `dose`
together with the other value columns and recomputes `diff` and the affected
pairs with the recorded value column. `plot.ISA()` plots `dose` when the ISA was
segmented on it (axis label `Dosage`), otherwise `cor.gc`. `print.RCA()` and
`summary.RCA()` report the normalization mode and settings.

### 11.3 Cohort orchestration

For more than one BAM, `analyzeBamInput()` counts and corrects every sample
first, then calls `normalizeReadCounts()` once, then segments every sample with
the same target pairs. `normalization = "reference"` requires
`reference_sample`; an unknown name fails with `MaMaMIAAnalysisError`. A single
BAM is still forced to `"none"`. The result gains a `samples` data frame with
`sample`, `path`, `assembly`, `index_path`, `reference`, `normalization`,
`n_windows`, `n_segments` and `n_putative_introgressions`. Provenance adds
`n_samples` and `total_windows`. A `"normalization"` progress stage is reported
between the correction and segmentation stages.

The performance guard warns when `n_samples * n_windows` exceeds
`max_fit_windows`. The existing `max_fit_windows` setting is used because the
configuration schema defines no separate cohort-window limit; the fit cap
itself is applied per sample by `correctReadCounts()` (section 13.6).

### 11.4 Validation

The synthetic variable-GC fixture builds three samples with different read
densities (`step` 4, 6 and 8). Tests cover the hand-computed three-sample
dosage table, reference dose = 1, all-zero windows, missing values, one sample,
shuffled window orders, approximate scale invariance, by-side vs joint mode,
invalid inputs, a 3-BAM median cohort, a 3-BAM reference cohort, duplicate
sample names, a nearly empty sample and the performance warning. The gated
golden regression confirms that segmentation without dose is unchanged.

---

## 12. Wave 5 implementation notes

### 12.1 Dashboard introgression plot

`explorer_results_ui()` gained an "Introgression plot" card with
`radioButtons("plot_type", c("Mirror" = "mirror", "Difference" = "diff"))` and
`plotOutput("isa_plot", height = "420px")`. `explorer_server()` renders the
thresholded ISA through `plot.ISA()`. Before plotting, `filter_dashboard_isa()`
narrows the ISA to the selected chromosomes and chromosome pairs and to segments
that pass the minimum-width filter, so the plot and the table agree on the
visible calls. `putative_only = FALSE` passes `all_segments = TRUE`, which keeps
the plot in step with the "show putative introgressions only" checkbox.

The render validates when no analysis is loaded ("Run the analysis to plot") and
when the filters leave no window ("Select at least one chromosome and pair to
plot"). A failure inside `plot.ISA()` is reported as a warning notification and
a plot-level validation message; it does not break the table or the download.
The plot card is independent of the parked BioCircos view, which stays disabled.

A 280,000-window synthetic ISA rendered in about 7.6 s (`diff`) and 11.1 s
(`mirror`). `plot.ISA()` gained `max_points` (default `Inf`, so all other callers
are unchanged): a finite cap keeps the peak absolute window in each of
`max_points` bins, while the axis and the segment lines still use every window.
The dashboard passes `dashboard_max_plot_points = 20000L`, which renders the
same case in about 2.7 s and 1.6 s. The BEDPE export always uses the full,
unthinned ISA.

### 12.2 BEDPE export with sample identity

`writeToBedpe(ISA, file, sample = NULL)` keeps the legacy 14-column layout when
`sample` is `NULL`; the table, header and numeric formatting are unchanged.
With a sample name, a `sample` column is appended as the last column. The
zero-row path keeps the header (including `sample`) and warns.

`writeCohortToBedpe(result, file, samples = NULL, selected = NULL)` accepts a
`BamAnalysisResult` or a named list of `ISA` objects. `samples` selects and
orders the exported samples. `selected` is a named list of 1-based row ids into
each sample's `ISA$segments` (the explorer's `.row_id`):

- `selected = NULL` exports every putative introgression of every exported
  sample.
- A sample named in `selected` exports exactly the listed rows, even when a row
  is not flagged as a putative introgression.
- A sample absent from `selected` contributes no rows.
- An explicit `NULL` entry exports every putative introgression for that
  sample; `integer(0)` exports none.

Row ids are validated per sample (whole numbers within `1..nrow(segments)`) and
samples are validated against the result. Zero rows still produce a header, with
`sample` last. The explorer passes its `sample` argument to `writeToBedpe()` so a
per-sample download is already tagged; the combined cohort download is wired in
T6.3 with the per-sample selection map.

### 12.3 Validation

Tests cover the exact legacy BEDPE bytes with `sample = NULL`, the appended
sample column, header preservation with zero rows, cohort combination, sample
subsetting, per-sample row selection, samples missing from the selection, and
invalid sample names or row ids. Dashboard tests cover
`filter_dashboard_isa()`, `testServer` rendering of both plot types, the
no-data and no-selection validation paths, and downloads with and without a
sample name. A gated synthetic two-sample cohort run confirms that a real
`BamAnalysisResult` exports a combined BEDPE with `sample` last.

---

## 13. Wave 6 implementation notes

### 13.1 Entrypoint dispatch

`runMaMaMIA(analysis_dir, assemblies_dir, cache, config, bam_paths, input,
max_request_size, host, port, launch.browser)` keeps the legacy arguments and
adds the BAM-first ones. `input` is one of `"auto"`, `"bam"`, `"files"` or
`"example"`. `"auto"` selects the BAM app when `config` or `bam_paths` is
supplied and the legacy app otherwise; `"bam"` always selects the BAM app and
`"files"`/`"example"` keep the legacy launcher available through T7.5.
`max_request_size` (default 5 GiB) is passed to `options(shiny.maxRequestSize)`.
The internal `mamamia_app()` builds the app object without launching it and
records the resolved mode in the `mamamia_input` attribute, which tests use.

### 13.2 BAM-first application

`bam_launcher_app(config, bam_paths = NULL, cache = NULL)` is the internal app
factory. `config` is a YAML path, a validated configuration list or `NULL`; a
`NULL` falls back to `MAMAMIA_CONFIG`. A missing or invalid registry renders a
configuration error card instead of stopping the process.

Sidebar controls: configuration summary; `fileInput("bams")`;
`pickerInput("server_bams")` when `server.bam_root` is configured or initial
`bam_paths` are supplied (bounded to depth 3 and 2,000 files, listed once per
session); resolution radio buttons built from `defaults.allowed_windows` (Low =
100, Medium = 1,000, High = 10,000, default
`defaults.windows_per_chromosome`); the reference selector (hidden for one
BAM); the assembly override; advanced MAPQ, exclude-flags, cores and
`max_fit_windows`; and a Run button that carries the `disabled` attribute and
the Bootstrap `disabled` class until preflight succeeds.

Main cards: a preflight table (`file`, `sample`, `assembly`, `targets`,
`index`, `warnings`), the detected-assembly card (label, description, donor and
recipient contigs, pairs, window counts per allowed resolution) and a problems
card with the error stage, message, suggested fix and, for assembly failures,
the candidate scores. Candidate scores are computed with the same
`.match_candidate()`/`.candidate_table()` helpers as `matchAssembly()`, so they
are available even for header-length mismatches whose condition payload has no
candidate table.

### 13.3 Preflight, staging and run flow

The preflight reactive collects uploaded `datapath`s and server-side picker
paths, stages each BAM (and an existing CSI sibling) into a per-session
directory (symlink for files of at least 1 MiB, copy for smaller files and when
symlinks are unsupported), and calls `preflightBams()` on the staged paths. The
copy path also means temporary-directory cleanup that follows and chmods
symlinks cannot touch the source metadata. Staging normalizes sources to
absolute paths, so relative picker values cannot create dangling symlinks. The
table reports the original path as `source` but keeps the staged path for the
run.

The run handler validates the controls through `bam_run_parameters()` (single
BAM → `none`; cohort with an empty reference → `median`; a named reference →
`reference`; resolution from `allowed_windows`; MAPQ 0-255; flags 0-4095; cores
bounded by `detectCores()`; positive `max_fit_windows`), then calls
`analyzeBamInput()` inside `shiny::withProgress()`. `bam_progress_stages()`
reports `preflight`, `index`, `windows`, `gc`, one `counts`, `correction` and
`segmentation` per sample, plus `normalization` for a cohort.
`bam_progress_label()` adds the detail reported by each stage to the progress
text: the BAM count and detected assembly, reused versus staged indexes, total
and per-chromosome windows, reads per sample, ideal fit windows, the
normalization mode and reference, and segments per sample. Every failure is
caught, classified by `bam_error_stage()`, and shown as
`[stage] message / File: <source> / Suggested fix: ...`; the previous result is
restored and the last error remains visible in the Problems card. A
`MaMaMIASampleError` is mapped back to its source BAM through the preflight
table, so the notification names the file the user supplied rather than the
staged session copy. The result is invalidated when the BAM selection,
resolution, assembly, reference, MAPQ, exclude flags or `max_fit_windows`
change, so a stale BEDPE cannot be downloaded for different inputs.

`bam_resolution_problem()` mirrors the `makeAnalysisWindows()` length check
before the run starts. A configured contig shorter than the selected resolution
keeps the Run button disabled and displays the reason; the same check also runs
inside the handler as defense in depth. The registry validator rejects such a
configuration at startup, so this guard only matters when a server reuses a
validated configuration whose allowed set changed. After a successful run,
zero or near-zero coverage on a side and all-zero windows in a reference sample
are shown both in the cohort card and as a warning notification.

`session$onSessionEnded()` removes the session directory, including every
staged BAM and generated CSI index. `server.bam_root` is never written to
because the pipeline only sees the staged copies. GC index handles are
process-wide and are not cleared on session end. The analysis is synchronous:
the UI states that the session blocks and no cancellation path is offered.

### 13.4 Multi-sample results

`bam_cohort_summary()` reports `sample`, `reads_donor`, `reads_recipient`,
`reads_total`, `n_segments`, `n_putative_introgressions`, `normalization` and
`reference`. The results card shows it with a `selectInput("result_sample")`
for cohorts and a combined `downloadButton("download_cohort")`.
`bam_coverage_warnings()` flags a side with zero reads or fewer than one read
per window on average and reference samples with all-zero windows. The
selection map is a named list of `.row_id` vectors initialized to
`integer(0)` per sample; the explorer reports changes through
`on_selection(rows, sample)` and preselects the stored rows through
`selection`. `download_cohort` calls
`writeCohortToBedpe(result, file, selected = selection_map())`, so every sample
with selected rows contributes to one file with `sample` last.

### 13.5 Validation

`tests/testthat/test-dashboard-bam.R` covers the configuration cards, the
`MAMAMIA_CONFIG` fallback, resolution/assembly helpers, the error taxonomy, run
parameters, the cohort summary and coverage warnings, selection maps, preflight
text and button state, upload staging, assembly candidates, missing paths, both
server limits, a full single-BAM `testServer` run (plot and download), a
two-sample cohort run with selector, selection persistence and combined BEDPE,
a failing run that stays inside the server, and the entrypoint dispatch. The
dashboard tests add preselection and `on_selection` coverage for the module.

### 13.6 High-resolution fit subsampling

`correctReadCounts()` gained `max_fit_windows = Inf`. When the ideal set is
larger, `select_fit_windows()` draws a deterministic sample from the
cross-classification of subgenome and GC decile. Every non-empty stratum keeps
at least ten windows (or all of its windows when it is smaller); remaining
slots are allocated in proportion to the remaining capacity, with the largest
fractional remainders receiving the leftover slots. A cap below the required
per-stratum minimum stops with an error that names the minimum valid value.

Sampling uses a seed derived from the full-set cache key (`fit_sampling_seed()`)
and restores the caller's RNG state. `fit_cache_key()` covers
`max_fit_windows`, the fitted rows and a fingerprint of the selected rows.
Prediction still covers every window. The `RCA$fit_sampling` metadata records
`sampled`, `max_fit_windows`, `n_ideal`, `n_fitted`, `seed` and the per-stratum
allocation; `summary.RCA()` prints the fitted/ideal counts. `analyzeBamInput()`
passes its `max_fit_windows` to every sample, sets
`provenance$max_fit_windows_enforced = TRUE`, and keeps the cohort-wide
performance warning.

Validation: `tests/testthat/test-cache.R` covers the new key fields and seed
stability; `tests/testthat/test-correction.R` covers the sample size, the
per-stratum minimum, reproducibility, RNG isolation, cache invalidation,
finite predictions, invalid caps and a full-versus-subsample comparison on the
packaged `triticum` pair; `tests/testthat/test-pipeline-bam.R` confirms the cap
is applied per sample and recorded in provenance.

`inst/benchmarks/benchmark-fit-subsampling.R` runs 28,000 and 280,000 windows
in fresh processes and records elapsed time and peak RSS in
`inst/benchmarks/fit-subsampling-results.csv`. On an 8-core host, the 280,000
window run fitted 50,000 of 273,028 ideal windows in 47.0 s with a 4.2 GB peak
RSS, below the 8 GB budget; the 28,000-window full fit took 9.8 s and 1.1 GB.

---

## 14. Wave 7 implementation notes

### 14.1 Progress and error reporting

`bam_progress_label(stage, detail)` renders the pipeline detail forwarded by
`analyzeBamInput()`: BAM count and assembly for `preflight`, reused/staged
indexes for `index`, total and per-chromosome counts for `windows`, window
counts for `gc`, reads for `counts`, ideal windows for `correction`, mode and
reference for `normalization`, and segments for `segmentation`. Labels fall
back to the plain stage name when the pipeline calls `progress()` without
detail, so the callback remains optional.

`bam_error_info()` returns `stage`, `class`, `message`, `files` and
`suggestion`. Files are resolved in this order: a `path` in the condition
payload, the sample name in a `MaMaMIASampleError` mapped through the preflight
`sample`/`source` columns, then the caller's fallback list. Notifications and
the Problems card use `File:` for one file and `Files:` for several, capped at
five paths with a `(+N more)` suffix. `bam_error_file_label()` is the single
formatter for notifications, the preflight error card and the persistent
last-run error card.

Run-control errors are MaMaMIA conditions: an invalid resolution raises
`MaMaMIAWindowError`, a bad MAPQ/flag/core/fit-window value raises
`MaMaMIAAnalysisError`, and an unknown reference raises `MaMaMIAAnalysisError`
with a `Reference sample ...` message that `bam_error_stage()` maps to
`normalization`. This keeps every user-facing failure on the same
stage/message/file/suggestion path as pipeline failures.

### 14.2 Limits, guards and cancellation

`server.max_bams` and `server.max_total_bytes` are enforced by `preflightBams()`
before headers are parsed, so an oversized selection never reaches counting.
The error is a `MaMaMIALimitError` and maps to the `limits` stage. Limits are
checked on staged paths as well, so an upload cannot bypass them.

Pathological inputs are handled as follows:

| Input | Behaviour |
|---|---|
| Contig shorter than the selected resolution | `makeAnalysisWindows()` stops with `MaMaMIAWindowError`; the dashboard disables Run before the analysis through `bam_resolution_problem()`. Config validation rejects such a registry at startup. |
| Fewer than 2 valid or 10 ideal windows | `correctReadCounts()` stops; `analyzeBamInput()` wraps it as `MaMaMIASampleError` with the `GC correction` stage; the app maps it to `correction` and suggests a lower resolution or a larger fit cap. |
| Zero or near-zero reads on a side | Normalization uses the pseudocount dosage rule and completes; the cohort card and a warning notification name the sample and side. |
| All-zero reference windows | `(0 + eps) / (0 + eps) = 1` by the frozen D5 rule; the cohort warning reports the all-zero window count for the reference sample. |

No background process or wall-clock timeout is used: the analysis runs in the
session's R process and blocks it. The dashboard states this in the sidebar and
the deployment notes set Shiny, ShinyProxy and reverse-proxy timeouts longer
than the largest expected run. A browser disconnect is handled after the
current pipeline stage returns, when Shiny regains control and
`session$onSessionEnded()` removes the staging directory.

### 14.3 Validation

`tests/testthat/test-dashboard-bam.R` now covers the detailed progress labels,
the file/sample attribution and file-label truncation, the short-contig guard
(unit and `testServer`), a mocked `analyzeBamInput()` that raises every error
category and checks the stage, suggestion, file, persistent card, previous
result and continued session, an invalid reference selection, the real
correction failure stage, and the existing end-to-end, cohort, limit and
cleanup tests. `tests/testthat/test-normalization.R` covers all-zero median and
reference windows; `tests/testthat/test-correction.R` covers the fit-window cap
and its errors.
