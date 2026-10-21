# Fit-subsampling benchmark (T7.1)

`benchmark-fit-subsampling.R` builds a synthetic 28-contig assembly (14 donor
and 14 recipient chromosomes) and runs `correctReadCounts()` at the two
high-resolution sizes from the plan:

| Windows | Windows per contig | Ideal windows | Fitted windows | Mode |
|---:|---:|---:|---:|:---|
| 28,000 | 1,000 | 27,328 | 27,328 | full fit |
| 280,000 | 10,000 | 273,028 | 50,000 | stratified subsample |

Each case runs in a fresh R process, so `VmHWM` from `/proc/self/status` is the
peak resident set size of that run. `peak_heap_mb` is the R heap high-water
mark after `gc(reset = TRUE)`.

## Recorded results

`fit-subsampling-results.csv` (8 cores, `max_fit_windows = 50000`):

| Windows | Elapsed (s) | Peak RSS (MB) | R heap (MB) | Finite `cor.gc` |
|---:|---:|---:|---:|:---:|
| 28,000 | 9.842 | 1,134.1 | 249.0 | yes |
| 280,000 | 46.990 | 4,213.8 | 333.9 | yes |

The 280,000-window run stays within the plan's 8 GB budget. The cap keeps the
fit at 50,000 ideal windows; prediction still covers all 280,000 windows.

## Reproduce

```sh
Rscript inst/benchmarks/benchmark-fit-subsampling.R \
    --cores=8 --max-fit-windows=50000
```

The no-argument form runs both sizes and writes
`inst/benchmarks/fit-subsampling-results.csv`. Use `--out=DIR` to write
elsewhere. The synthetic data are deterministic; timings vary by host.
