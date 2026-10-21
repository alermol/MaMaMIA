# MaMaMIA production image.
#
# Build from the repository root so the image gets the package sources:
#
#   docker build -f deployment/Dockerfile -t mamamia:0.99.0 .
#
# The image pins the uvr release, R and every R package (uvr.lock). It runs as
# an unprivileged user and is meant to be started by ShinyProxy through
# `uvr run /app/deployment/start.R`.

FROM debian:trixie-slim

ARG UVR_VERSION=0.4.7
ARG UVR_SHA256=e064ecbdba342114bd5dcf23158685d7cea8aa936a9654b5da2394e90f8810bf
ARG UVR_TARGET=x86_64-unknown-linux-gnu
ARG R_VERSION=4.6.1
# Optional Bioconductor mirror. The default lockfile points at
# https://bioconductor.org; set this to a mirror root (for example
# https://bioconductor.statistik.tu-dortmund.de/packages) when the build host
# cannot reach bioconductor.org. The package checksums still come from the
# committed lockfile.
ARG BIOC_MIRROR=

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# GCC 14 (Debian trixie) provides C++23. The development libraries are needed
# for source builds of Rhtslib, Bioconductor packages and MaMaMIA itself.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        make \
        gcc \
        g++ \
        gfortran \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libicu-dev \
        libreadline-dev \
        libpcre2-dev \
        libdeflate-dev \
        libtirpc-dev \
        libblas-dev \
        liblapack-dev \
        libncurses-dev \
        libuv1 \
        libjpeg-dev \
        libpng-dev \
        libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# Pin uvr by version and checksum; `latest` is never used in the image.
RUN curl -fsSL -o /tmp/uvr.tar.gz \
        "https://github.com/nbafrank/uvr/releases/download/v${UVR_VERSION}/uvr-${UVR_TARGET}.tar.gz" \
    && echo "${UVR_SHA256}  /tmp/uvr.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/uvr.tar.gz -C /usr/local/bin \
    && chmod 0755 /usr/local/bin/uvr \
    && rm -f /tmp/uvr.tar.gz \
    && uvr --version

# The runtime user owns the uvr-managed R, the package library and the caches.
RUN useradd --create-home --shell /bin/bash --uid 1000 mamamia \
    && mkdir -p /app /tmp/mamamia-cache /srv/mamamia \
    && chown -R mamamia:mamamia /app /tmp/mamamia-cache /srv/mamamia

ENV HOME=/home/mamamia \
    OMP_NUM_THREADS=1 \
    R_USER_CACHE_DIR=/tmp/mamamia-cache \
    MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml

WORKDIR /app

COPY --chown=mamamia:mamamia uvr.toml uvr.lock .r-version ./

USER mamamia

# Install the pinned R under the runtime user's home and resolve the full
# lockfile (runtime and development dependencies; shiny/bslib/DT live in the
# project's dev-dependencies). P3M binaries are used when available and the
# GCC 14 toolchain covers the source fallbacks.
RUN if [ -n "${BIOC_MIRROR}" ]; then \
        sed -i "s#https://bioconductor.org/packages/#${BIOC_MIRROR}/#g" uvr.lock; \
    fi \
    && uvr r install "${R_VERSION}" --unattended \
    && uvr sync --frozen --unattended --jobs 4 \
    && uvr doctor --unattended

COPY --chown=mamamia:mamamia DESCRIPTION NAMESPACE LICENSE.md README.md ./
COPY --chown=mamamia:mamamia R ./R
COPY --chown=mamamia:mamamia src ./src
COPY --chown=mamamia:mamamia inst ./inst
COPY --chown=mamamia:mamamia data ./data
COPY --chown=mamamia:mamamia man ./man
COPY --chown=mamamia:mamamia deployment ./deployment

# Install MaMaMIA into the uvr library. The activation shim puts the
# uvr-managed R on PATH and R_LIBS_USER on .uvr/library.
RUN uvr activate --write-shim \
    && . .uvr/activate \
    && R CMD INSTALL --no-multiarch --with-keep.source . \
    && Rscript -e 'library(MaMaMIA); cat("MaMaMIA", as.character(packageVersion("MaMaMIA")), "installed\n")'

EXPOSE 3838

# ShinyProxy overrides this with the same command; keeping it in the image
# makes local smoke tests realistic.
CMD ["uvr", "run", "/app/deployment/start.R"]
