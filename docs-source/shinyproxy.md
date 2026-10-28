# ShinyProxy deployment

This runbook builds the production image and runs MaMaMIA under
[ShinyProxy](https://shinyproxy.io/). It covers the data prerequisites, the
image build, a local smoke test, the proxy limits, and session cleanup.

Artifacts:

| File | Purpose |
|---|---|
| `deployment/Dockerfile` | Debian trixie image with GCC 14, pinned uvr 0.4.7, R 4.6.1 and the full `uvr.lock` set |
| `deployment/start.R` | Entrypoint: `MaMaMIA::runMaMaMIA(host = "0.0.0.0", port = 3838, config = Sys.getenv("MAMAMIA_CONFIG"))` |
| `deployment/shinyproxy.yml` | ShinyProxy spec, read-only mounts, limits and upload settings |
| `deployment/nginx.conf` | Optional reverse proxy with large bodies and long timeouts |
| `deployment/smoke-test.R` | Runs the full pipeline in the image on the bundled demo BAM (or `MAMAMIA_SMOKE_BAM`) |
| `deployment/gcidx-checksums.txt` | Format, sizes and SHA-256 of the validated gcrank indexes |
| `tools/validate-config.R` | Validates a production registry and its GC indexes before mounting it |

ShinyProxy itself is **not** required to build the image or run the local smoke
test. It is only needed for the final proxied deployment.

## 0. Runbook checklist

1. Build or copy the `gcrank` indexes under `/srv/mamamia/gcidx/` and write the
   species registry to `/srv/mamamia/config/assemblies.yml` (sections 1 and
   1.1). Record the index checksums in `deployment/gcidx-checksums.txt`.
2. Validate the mounted registry: `Rscript tools/validate-config.R
   /srv/mamamia/config/assemblies.yml`. Do not continue on a non-zero exit.
3. Build the image from the repository root (section 2) and run
   `uvr doctor --unattended` inside it.
4. Run the local smoke test (section 3). The script must print a result with
   at least one segment. It uses the bundled demo BAM by default, so no BAM
   mount is needed.
5. Start the image through the ShinyProxy command
   `uvr run /app/deployment/start.R`; `MAMAMIA_CONFIG` must point at the
   mounted registry (section 4).
6. Put the reverse proxy in front, with the upload and timeout limits raised
   together (section 5), and confirm the resource limits (section 6).
7. Check session cleanup and authentication before exposing the app (section
   7), and use section 8 when a run fails.

A registry change is a new deployment: replace the file atomically and start a
fresh app container, because the config is read once per process (section 1.1).
New index directories alone are visible without a restart.

## 1. Data prerequisites

The application reads a YAML registry (see
`inst/config/assemblies.example.yml`) and never writes to it. For every
species the registry must provide:

- a `gcrank` index directory per genome with `contigs.tsv`, `manifest.json`,
  `genome.bv` and `genome.rb` (the `index.gcrank` directory produced by
  `tools/build-gcrank-index.sh`; the format version must be `0.1.0`);
- contig ids that exist in those indexes, with a `component` and a chromosome
  group (from `group` or the first number in `name`);
- a `default_component` for a polyploid species, so the application does not
  have to ask which component to use.

A species has no donor/recipient role. The user chooses the recipient and the
donor in the app, and the assembly is derived from those two species. The
registry declares no `assemblies:` section; a registry with one is rejected at
startup. A registry needs at least two species.

- a largest `allowed_windows` value no longer than the shortest contig that
  can enter an analysis;

Validate the exact file that will be mounted:

```sh
Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml
```

The script calls `readMaMaMIAConfig()` (the same validation as application
startup) and additionally checks the index payload files. It exits non-zero on
the first problem. Do not ship a registry that only works on the recipient
side: every contig that can pair with a homoeolog must be present in an index.

The indexes validated for this deployment are recorded in
`deployment/gcidx-checksums.txt`: format version, contig/bp totals, SHA-256 for
`contigs.tsv`, `genome.bv`, `genome.rb`, `manifest.json` and their CRC sidecars,
plus the source-FASTA digest. `manifest.json`'s `source.sha256` covers the
uncompressed FASTA stream, while the `.sha256` sidecar next to the gzip file
covers the compressed file.

The expected host layout is:

```
/srv/mamamia/config/assemblies.yml
/srv/mamamia/gcidx/<species>/index.gcrank/
```

BAMs are uploaded by users; the container does not mount a BAM directory.
Uploaded BAMs are staged into a per-session temporary directory and the source
upload directory is never modified. The app also offers a bundled synthetic
example whose pre-computed tables are part of the image.

### 1.1 Growing the registry (new species and indexes)

The container mounts the **whole** `/srv/mamamia/gcidx` directory, not a fixed
list of index paths. A new
`/srv/mamamia/gcidx/<species>/index.gcrank/` therefore becomes visible inside
already-running containers as soon as it is complete; no image rebuild, no
`uvr sync`, and no ShinyProxy restart is needed. Indexes are memory-mapped on
demand, so mounting all of them does not load all of them into memory.

The registry itself is different: `bam_launcher_state()` reads and validates
`assemblies.yml` once, when the app process starts (the UI and server share
that state). A newly started ShinyProxy app container sees the new registry;
an already-running container keeps the configuration it started with. Close and
reopen the app so ShinyProxy starts a fresh container after a registry change.
This also applies to pre-initialized/shared containers: sessions in one R
process share the config loaded at process start.

Keep the config as a **directory** mount (`/srv/mamamia/config`), as in
`deployment/shinyproxy.yml`. A single-file bind mount pins the old inode when
the file is replaced atomically. Install a new registry with a rename:

```sh
# 1. build and validate the new index somewhere under /srv/mamamia/gcidx
Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml.new

# 2. atomically replace the registry; running containers see the new file
mv /srv/mamamia/config/assemblies.yml.new /srv/mamamia/config/assemblies.yml
```

An index added under `/srv/mamamia/gcidx` is visible immediately, but the
`species:` entry that references it is part of the registry and therefore
follows the process-lifetime rule above. Adding a species means adding an
index (or reusing one already mounted) and one `species:` entry; no preset
assembly and no contig pairs are required. After the restart, the new species
appears in both the recipient and the donor selector.

### 1.2 Species-pair selection

The registry lists role-free species. The sidebar asks for a **recipient
species** and a **donor species**; the same species can fill either role, so
swapping the two selectors reverses the comparison. There are no per-component
selectors: every component declared for a species is analysed in one assembly,
and the **Assembly** panel shows a read-only summary such as
`recipient A+B+D; donor S`.

Contigs are paired by chromosome group across all selected components. The
group is the optional `group` field or the first number in `name`, so rye `1R`
pairs with wheat `1A`, `2R` with `2A`, and so on. Only paired contigs enter the
analysis; the card shows which contigs were excluded and why.

The sidebar has no predefined-assembly selector. A "Swap donor and recipient"
button reverses the comparison without reselecting the species. These
selectors live in the **Assembly** panel, which appears after preflight
succeeds. The collapsed **Advanced** panel exposes donor aggregation
(`max` by default, or `mean`/`sum`); the registry supplies the initial value,
and changing it clears the current result until the analysis is rerun. The
resolution, upload limits, resource limits and registry-growth caveats are
unchanged.

## 2. Build the image

From the repository root:

```sh
docker build -f deployment/Dockerfile -t mamamia:0.99.2 .
```

If the build host cannot reach `bioconductor.org`, pass a mirror (the package
checksums still come from the committed lockfile):

```sh
docker build --network=host \
  --build-arg BIOC_MIRROR=https://bioconductor.statistik.tu-dortmund.de/packages \
  -f deployment/Dockerfile -t mamamia:0.99.2 .
```

The validation sandbox needed `--network=host` and
this mirror because its Docker daemon could not resolve or reach the default
hosts; the default build path remains the public bioconductor.org URLs.

The build:

1. installs the compiler and system libraries needed by Rhtslib,
   Bioconductor, Matrix, stringi and xml2;
2. downloads uvr 0.4.7 for `x86_64-unknown-linux-gnu` and verifies its SHA-256
   before installing it;
3. creates the unprivileged `mamamia` user and installs R 4.6.1 under its home;
4. runs `uvr sync --frozen --unattended` from the committed lockfile and
   `uvr doctor`;
5. installs `fontconfig` and `fonts-dejavu-core`. uvr's R bundles the cairo and
   fontconfig shared libraries but not `/etc/fonts/fonts.conf` or fonts; without
   them the first `renderPlot()` dies with `Fontconfig error: Cannot load
   default config file` followed by a segfault in `grid.Call(C_stringMetric, ...)`;
6. installs MaMaMIA into `.uvr/library` with the uvr-managed R and renders a
   cairo/ggplot PNG, so a missing font setup fails the build before the image is
   tagged.

The image is not based on a preinstalled R and does not use the local conda
environment. `OMP_NUM_THREADS=1` is set at build and runtime so a sample's
result does not depend on the worker it ran on.

## 3. Local smoke test

Start the image with the real read-only indexes:

```sh
docker run --rm -p 3838:3838 \
  -e MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml \
  -v "$PWD/inst/config/assemblies.example.yml:/srv/mamamia/config/assemblies.yml:ro" \
  -v "$PWD/test_data/gcidx:/srv/mamamia/gcidx:ro" \
  mamamia:0.99.2
```

Open `http://127.0.0.1:3838/`, confirm that the calculation log shows the two
species registry, choose `Triticum aestivum CS T2T` as the recipient and
`Secale cereale` as the donor (or the reverse) in the **Assembly** panel, then
either upload a BAM or press **Example data**. The example needs no BAM mount
and takes about a second; a successful run populates the **Karyotype**,
**Score profile**, **Calls** and **Export** tabs, and the footer log ends with
the analysis summary. Stop the container with Ctrl-C; the session
directory and staged indexes are removed on session end.

To check the same pipeline without the browser, run the bundled smoke script in
a one-off container. It uses the packaged demo BAM and registry, so it needs no
BAM or registry mount (set `MAMAMIA_SMOKE_BAM` and `MAMAMIA_CONFIG` to exercise
a production registry instead):

```sh
docker run --rm mamamia:0.99.2 uvr run /app/deployment/smoke-test.R
```

The script resolves the first two species when an external registry is used.
Override the pair with `MAMAMIA_SMOKE_RECIPIENT` and `MAMAMIA_SMOKE_DONOR` when
the registry order is not the intended comparison. Expected output on the demo
at resolution 1000: 2,000 windows, a non-zero segment count, a
`calibrated calls:` line, and a `plot:` line with a non-zero PNG size. The
script prints every pipeline stage and fails if the result is not a
`BamAnalysisResult` with at least one segment. It renders a cairo text probe
before the pipeline and the real score profile at the end, so an image without
fontconfig or fonts fails in the smoke test rather than in the dashboard.

Check the image without starting the app:

```sh
docker run --rm mamamia:0.99.2 uvr doctor --unattended
```

For a startup-only check with the real mounts, run the full command above and
request the app page with `curl -fsS http://127.0.0.1:3838/`.

## 4. Run under ShinyProxy

1. Copy `deployment/shinyproxy.yml` to the ShinyProxy host as
   `/opt/shinyproxy/application.yml`, or start the jar with
   `--spring.config.location=file:/path/to/shinyproxy.yml`.
2. Create `/srv/mamamia/{config,gcidx}` and mount them read-only for the
   ShinyProxy user. The container mounts in the spec repeat the same paths.
   BAMs are not mounted; users upload them.
3. Make the image available to the Docker daemon that ShinyProxy uses
   (`docker build`, `docker save`/`load`, or a registry).
4. Start ShinyProxy. The app appears as **MaMaMIA** and starts containers with:

```sh
uvr run /app/deployment/start.R
```

The image's `WORKDIR` is `/app`, so `uvr run` finds `uvr.toml`, `uvr.lock` and
`.uvr/library` without relying on `.Rprofile`. `MAMAMIA_CONFIG`,
`OMP_NUM_THREADS`, `R_USER_CACHE_DIR` and `HOME` come from `container-env`.

`deployment/shinyproxy.yml` uses `container-volumes` (the ShinyProxy key for
bind mounts; `volumes` is not a ShinyProxy property). If ShinyProxy itself
runs in a container on the same Docker host, set
`SHINYPROXY_INTERNAL_NETWORKING=true`.

## 5. Reverse proxy, uploads and timeouts

Put `deployment/nginx.conf` in front of ShinyProxy (listener
`127.0.0.1:8080`). It sets:

- `client_max_body_size 6g` (above the 5 GiB registry limit);
- `proxy_request_buffering off` so large BAM uploads stream;
- one-hour read/send timeouts;
- WebSocket upgrade headers for Shiny.

Three limits must stay consistent:

| Layer | Setting | Default |
|---|---|---|
| Shiny | `max_request_size` in `deployment/start.R` | 5 GiB |
| ShinyProxy | `spring.servlet.multipart.*` and `server.undertow.max-http-post-size` | 6 GB file / 7 GB request / unlimited |
| nginx | `client_max_body_size` | 6 GB |

Raise all three together when the registry allows larger cohorts. A
high-resolution multi-sample run can take tens of minutes; the ShinyProxy
`heartbeat-timeout` is 10 minutes and the session timeout is 12 hours. The
analysis is synchronous, so the R session blocks during a run; there is no
in-app cancellation.

## 6. Resources and concurrency

- Medium resolution (28,000 windows): one BAM completes in about 20 s on one
  core; peak memory is dominated by R, htslib and the gcrank index mappings.
- High resolution (280,000 windows): about ten times the window table; the
  per-window measurements are bounded by one bitset per sample and window.
- `container-memory-limit: 16g` and `container-cpu-limit: 4` match the
  documented high-resolution budget.
- `max-instances: 1` allows one analysis per user; `max-total-instances: 2`
  caps the host at two 16 GB containers.
- GC indexes are bind-mounted, not copied into the image, and are mapped
  read-only. Resident memory depends on touched pages and the host page cache
  is shared across user containers, so mounting all indexes does not multiply
  memory by the number of users. Do not load more indexes than the run needs.
- Do not clear the shared GC index cache at session end in a shared R
  process. The LRU bound handles it; `gcrankClearCache()` is for process
  shutdown or an explicit admin action.

## 7. Session cleanup and security

- Uploaded BAMs and generated CSI indexes live in the session temp directory
  and are removed by `session$onSessionEnded()`.
- Uploaded BAMs are staged by symlink into session temp (copied only when
  symlinking is not supported) and are never written in place. The container
  mounts no BAM directory.
- The image runs as the unprivileged `mamamia` user; the registry and indexes
  are mounted read-only; only `/tmp` and the R cache are writable.
- There are no arbitrary server paths: every BAM comes from a Shiny upload and
  is validated by preflight.
- ShinyProxy listens on `127.0.0.1` in the supplied configuration. Enable
  ShinyProxy authentication (LDAP/OIDC/simple) or an authenticating reverse
  proxy before exposing it.

## 8. Troubleshooting

| Symptom | Check |
|---|---|
| Container exits immediately | `docker logs <id>`; `MAMAMIA_CONFIG` must point at a mounted file |
| Configuration error card | run `tools/validate-config.R` on the mounted registry; paths are resolved relative to the YAML file |
| `uvr: ... package(s) not installed` | the build did not finish `uvr sync --frozen`; rebuild without a cached layer |
| Upload rejected before the app | the nginx body limit or `shiny.maxRequestSize`; raise both |
| Container killed during a run | `container-memory-limit`; lower the resolution at run time |
| App stopped while running | `heartbeat-timeout`/session timeout; raise them for long runs |
| `permission denied` reading indexes | the host mounts must be readable by the container's `mamamia` user (UID 1000) |
| `Fontconfig error: Cannot load default config file` followed by a segfault in `grid.Call(C_stringMetric, ...)` | the image has no fontconfig configuration or fonts. Install `fontconfig` and `fonts-dejavu-core` (both are in `deployment/Dockerfile`) and rebuild; the build-time cairo render and `deployment/smoke-test.R` now fail before a user hits it |
