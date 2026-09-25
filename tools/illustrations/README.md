# Illustration assets

Optional maintainer script for regenerating static vignette figures in
`vignettes/figures/`. These plotting packages are not package runtime deps.

```sh
/opt/R-4.6.1/bin/Rscript tools/illustrations/render-illustrations.R
/opt/R-4.6.1/bin/Rscript tools/illustrations/render-illustrations.R --install
/opt/R-4.6.1/bin/Rscript tools/illustrations/render-illustrations.R --rebuild-data
```

Uses `circlize`, `gggenes`, `ggraph`, `tidygraph`, `ggdist`, `ggtext`,
`patchwork`, `png`, and `ragg`. Caches fit/data in `/tmp/mamamia-fit.rds` and
`/tmp/mamamia-illus.rds`.

| File | What it shows | Main library |
|---|---|---|
| `01-mechanism.png` | Introgression as donor excess in real `1At/1A` coverage | `gggenes` + `patchwork` |
| `02-workflow.png` | Pipeline from sample to called introgression | `ggraph`/`tidygraph` |
| `03-circos.png` | Chromosome overview, coverage heatmap, introgression ribbons | `circlize` |
| `04-gc-bias.png` | Fitted GC response removed by correction | `ggplot2`/`ggtext` |
| `05-pair-profile.png` | Per-window difference, CBS segments, called interval | `ggplot2`/`ggtext` |
| `06-orientation.png` | Segmentation before/after donor chromosome reversal | `ggplot2`/`ggtext` |
| `07-arrowmap.png` | Chromosome-arrow map of called intervals | `gggenes` |
| `08-distributions.png` | Donor-recipient difference distributions by pair | `ggdist` |
