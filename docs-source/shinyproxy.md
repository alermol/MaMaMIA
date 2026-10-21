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
| `deployment/smoke-test.R` | Runs the full pipeline in the image on one mounted BAM |
| `deployment/gcidx-checksums.txt` | Format, sizes and SHA-256 of the validated gcrank indexes |
| `tools/validate-config.R` | Validates a production registry and its GC indexes before mounting it |

ShinyProxy itself is **not** required to build the image or run the local smoke
test. It is only needed for the final proxied deployment.

## 1. Data prerequisites

The application reads a YAML registry (see
`inst/config/assemblies.example.yml`) and never writes to it. For every
assembly the registry must provide:

- a `gcrank` index directory per genome with `contigs.tsv`, `manifest.json`,
  `genome.bv` and `genome.rb` (the `index.gcrank` directory produced by
  `tools/build-gcrank-index.sh`; the format version must be `0.1.0`);
- donor and recipient contig ids that exist in those indexes;
- pairs covering every donor contig that should be analysed;
- a largest `allowed_windows` value no longer than the shortest target contig.

Validate the exact file that will be mounted:

```sh
Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml
```

The script calls `readMaMaMIAConfig()` (the same validation as application
startup) and additionally checks the index payload files and `server.bam_root`.
It exits non-zero on the first problem. Do not ship a registry that only works
on the recipient side: every donor contig used in `pairs` must be present in an
index.

The indexes validated for this deployment are recorded in
`deployment/gcidx-checksums.txt`: format version, contig/bp totals, SHA-256 for
`contigs.tsv`, `genome.bv`, `genome.rb`, `manifest.json` and their CRC sidecars,
plus the source-FASTA digest. `manifest.json`'s `source.sha256` covers the
uncompressed FASTA stream, while the `.sha256` sidecar next to the gzip file
covers the compressed file.

The expected host layout is:

```
/srv/mamamia/config/assemblies.yml
/srv/mamamia/gcidx/<assembly>/index.gcrank/
/srv/mamamia/bams/<sample>/<name>.bam          # optional, read-only
```

`server.bam_root` in the registry must be `/srv/mamamia/bams` when the
server-side picker is used. Uploaded BAMs are staged into a per-session
temporary directory; the source directory is never modified.

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

An assembly entry still describes one donor/recipient pair, including its
contigs, subgenome labels, GC indexes and explicit contig pairs. Adding a
species pair means adding an index (or two) and one `assemblies:` entry; the
current UI selects the predefined assembly, not two independent species. See
the next section for the reason and the options.

### 1.2 Species-pair selection and the current UI

The frozen plan models an analysis as an **assembly** (T1.3, T6.1, D4): a
donor contig set, a recipient contig set, explicit donor-to-recipient pairs,
and the GC indexes that contain those contigs. The app auto-detects the
assembly from the BAM header and offers an override over the configured
assemblies. It does not offer two independent "donor species" and "recipient
species" selectors.

That is deliberate: pairing two arbitrary species needs a rule for which donor
contig corresponds to which recipient contig, and which subgenome labels to
use. In the wheat/rye case the seven rye chromosomes are paired explicitly with
1A, 2A, ..., 7A. Automatically pairing by name or order would be a scientific
decision, not a deployment detail. Adding a two-species selector therefore
requires a plan amendment and a species/pair registry; it is not solved by
mounting more indexes.

## 2. Build the image

From the repository root:

```sh
docker build -f deployment/Dockerfile -t mamamia:0.99.0 .
```

If the build host cannot reach `bioconductor.org`, pass a mirror (the package
checksums still come from the committed lockfile):

```sh
docker build --network=host \
  --build-arg BIOC_MIRROR=https://bioconductor.statistik.tu-dortmund.de/packages \
  -f deployment/Dockerfile -t mamamia:0.99.0 .
```

The validation sandbox needed `--network=host` and
this mirror because its Docker daemon could not resolve or reach the default
hosts; the default build path remains the public bioconductor.org URLs.

The build:

1. installs the compiler and system libraries needed by Rhtslib,
   Bioconductor, TMB/glmmTMB, Matrix, stringi and xml2;
2. downloads uvr 0.4.7 for `x86_64-unknown-linux-gnu` and verifies its SHA-256
   before installing it;
3. creates the unprivileged `mamamia` user and installs R 4.6.1 under its home;
4. runs `uvr sync --frozen --unattended` from the committed lockfile and
   `uvr doctor`;
5. installs MaMaMIA into `.uvr/library` with the uvr-managed R.

The image is not based on a preinstalled R and does not use the local conda
environment. `OMP_NUM_THREADS=1` is set at build and runtime so a sample's
result does not depend on the worker it ran on.

## 3. Local smoke test

Start the image with the real read-only indexes and BAMs:

```sh
docker run --rm -p 3838:3838 \
  -e MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml \
  -v "$PWD/inst/config/assemblies.example.yml:/srv/mamamia/config/assemblies.yml:ro" \
  -v "$PWD/test_data/gcidx:/srv/mamamia/gcidx:ro" \
  -v "$PWD/test_data/bam_files:/srv/mamamia/bams:ro" \
  mamamia:0.99.0
```

Open `http://127.0.0.1:3838/`, confirm that the configuration summary shows the
assembly, select a server-side BAM (for example
`Sq54/cs_t2t_se_cereale.bam`), wait for preflight, choose resolution 1000 and
run. A successful run shows the cohort summary and the explorer plot. Stop the
container with Ctrl-C; the session directory and staged indexes are removed on
session end.

To check the same pipeline without the browser, run the bundled smoke script in
a one-off container (it stages a symlink into `/tmp` exactly as the app does):

```sh
docker run --rm \
  -e MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml \
  -v "$PWD/inst/config/assemblies.example.yml:/srv/mamamia/config/assemblies.yml:ro" \
  -v "$PWD/test_data/gcidx:/srv/mamamia/gcidx:ro" \
  -v "$PWD/test_data/bam_files:/srv/mamamia/bams:ro" \
  mamamia:0.99.0 uvr run /app/deployment/smoke-test.R
```

Expected output on `Sq54` at resolution 1000: 28,000 windows, a non-zero
segment count, and `putative:` calls. The script prints every pipeline stage
and fails if the result is not a `BamAnalysisResult` with at least one segment.

Check the image without starting the app:

```sh
docker run --rm mamamia:0.99.0 uvr doctor --unattended
```

For a startup-only check with the real mounts, run the full command above and
request the app page with `curl -fsS http://127.0.0.1:3838/`.

## 4. Run under ShinyProxy

1. Copy `deployment/shinyproxy.yml` to the ShinyProxy host as
   `/opt/shinyproxy/application.yml`, or start the jar with
   `--spring.config.location=file:/path/to/shinyproxy.yml`.
2. Create `/srv/mamamia/{config,gcidx,bams}` and mount them read-only for the
   ShinyProxy user. The container mounts in the spec repeat the same paths.
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
high-resolution multi-BAM run can take tens of minutes; the ShinyProxy
`heartbeat-timeout` is 10 minutes and the session timeout is 12 hours. The
analysis is synchronous, so the R session blocks during a run; there is no
in-app cancellation.

## 6. Resources and concurrency

- Medium resolution (28,000 windows): about 4-6 GB peak for the GC fit plus R
  and index mappings.
- High resolution (280,000 windows) with `max_fit_windows: 50000`: about
  8-12 GB peak.
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
- Server-side BAMs are staged by symlink (copied only when symlinking is not
  supported) and are never written in place.
- The image runs as the unprivileged `mamamia` user; the registry, indexes and
  BAM root are mounted read-only; only `/tmp` and the R cache are writable.
- There are no arbitrary server paths: the BAM root comes from the registry,
  and every selected upload is validated by preflight.
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
| Container killed during a run | `container-memory-limit`; lower the resolution or `max_fit_windows` |
| App stopped while running | `heartbeat-timeout`/session timeout; raise them for long runs |
| `permission denied` reading indexes | the host mounts must be readable by the container's `mamamia` user (UID 1000) |
