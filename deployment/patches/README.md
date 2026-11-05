# Deployment patches

These patches are applied to the pinned dependency sources while building the
container image. They are not part of the MaMaMIA application sources.

## `httpuv-1.6.17-body-coalescing.patch`

**Why.** Shiny receives uploads through httpuv. In httpuv 1.6.17 the HTTP body
is handed to R once per socket read, which is typically 16-64 KiB. A 1 GiB BAM
therefore causes roughly 16,000 R callbacks and a single core saturates at
about 40-50 MB/s before Shiny even starts writing the file.

**What.** The patch accumulates request-body bytes in a C++ buffer and calls
`onBodyData` with 1 MiB blocks (and flushes the tail when the message
completes). Measured on the validation host, a 1 GiB upload dropped from
about 22 s to about 9 s; a real 814 MB BAM uploaded in 6.6 s with a matching
MD5.

**Pin.** The patch is written for httpuv 1.6.17, the version in `uvr.lock`.
The Dockerfile downloads that exact tarball and verifies its MD5 before
applying the patch. Do not reuse it with another httpuv release without
checking `src/httprequest.cpp` and `src/httprequest.h`.

**License.** httpuv is GPL (>= 2). The image is built from the pinned CRAN
source plus this patch; if you distribute the image, keep that source and patch
available to satisfy the license.

**Removing it.** Delete the `RUN` block that applies this patch from
`deployment/Dockerfile` and delete this directory. Uploads still work, at the
stock httpuv rate.
