# MaMaMIA production image.
#
# Build from the repository root so the image gets the application sources:
#
#   docker build -f deployment/Dockerfile -t mamamia:0.1.0 .
#
# On a host that cannot reach bioconductor.org, use host networking and the
# mirror argument instead:
#
#   docker build --network=host \
#     --build-arg BIOC_MIRROR=https://bioconductor.statistik.tu-dortmund.de/packages \
#     -f deployment/Dockerfile -t mamamia:0.1.0 .
#
# The image pins the uvr release, R and every dependency (uvr.lock). It runs as
# an unprivileged user and is started by ShinyProxy or directly through
# `uvr run /app/deployment/start.R`.

FROM debian:trixie-slim

ARG UVR_VERSION=0.4.7
ARG UVR_SHA256=e064ecbdba342114bd5dcf23158685d7cea8aa936a9654b5da2394e90f8810bf
ARG UVR_TARGET=x86_64-unknown-linux-gnu
ARG R_VERSION=4.6.1
# Optional Bioconductor mirror. The default lockfile points at
# https://bioconductor.org; set this to a mirror root (for example
# https://bioconductor.statistik.tu-dortmund.de/packages) when the build host
# cannot reach bioconductor.org. The dependency checksums still come from the
# committed lockfile.
ARG BIOC_MIRROR=

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# GCC 14 (Debian trixie) provides C++23. The development libraries are needed
# for source builds of Rhtslib, Bioconductor packages and the MaMaMIA native
# bridge.
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
        libuv1-dev \
        cmake \
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

# The runtime user owns the uvr-managed R, the dependency library and the caches.
RUN useradd --create-home --shell /bin/bash --uid 1000 mamamia \
    && mkdir -p /app/tools /tmp/mamamia-cache /srv/mamamia \
    && chown -R mamamia:mamamia /app /tmp/mamamia-cache /srv/mamamia

ENV HOME=/home/mamamia \
    OMP_NUM_THREADS=1 \
    R_USER_CACHE_DIR=/tmp/mamamia-cache \
    MAMAMIA_ROOT=/app \
    MAMAMIA_CONFIG=/srv/mamamia/config/assemblies.yml

WORKDIR /app

COPY --chown=mamamia:mamamia uvr.toml uvr.lock .r-version ./
COPY --chown=mamamia:mamamia deployment/patches /tmp/mamamia-patches

USER mamamia

# Install the pinned R under the runtime user's home and resolve the full
# lockfile. P3M binaries are used when available and the GCC 14 toolchain
# covers the source fallbacks.
RUN if [ -n "${BIOC_MIRROR}" ]; then \
        sed -i "s#https://bioconductor.org/packages/#${BIOC_MIRROR}/#g" uvr.lock; \
    fi \
    && uvr r install "${R_VERSION}" --unattended \
    && uvr sync --frozen --unattended --jobs "$(nproc)" \
    && uvr doctor --unattended

# httpuv 1.6.17 invokes an R callback for every 16-64 KiB request-body read.
# A multi-GB BAM upload therefore spends most of its time in tens of thousands
# of R calls before Shiny sees the file. Rebuild the pinned httpuv from source
# with the coalescing patch (1 MiB blocks); the version and uvr.lock entry are
# unchanged. Remove this block and deployment/patches/httpuv-* together.
RUN uvr activate --write-shim \
    && . .uvr/activate \
    && curl -fsSL -o /tmp/httpuv.tar.gz \
        "https://cran.r-project.org/src/contrib/httpuv_1.6.17.tar.gz" \
    && echo "50763765704e547c15aac4be97cec437  /tmp/httpuv.tar.gz" | md5sum -c - \
    && mkdir -p /tmp/httpuv-src \
    && tar -xzf /tmp/httpuv.tar.gz -C /tmp/httpuv-src --strip-components=1 \
    && cd /tmp/httpuv-src \
    && git apply -p1 /tmp/mamamia-patches/httpuv-1.6.17-body-coalescing.patch \
    && MAKEFLAGS="-j$(nproc)" R CMD INSTALL --no-multiarch . \
    && cd / \
    && rm -rf /tmp/httpuv.tar.gz /tmp/httpuv-src

# uvr's R bundles cairo and fontconfig libraries, but not the system
# fontconfig configuration or any fonts. Without /etc/fonts/fonts.conf, the
# first grid string metric (for example an axis label in renderPlot()) aborts R
# with "Fontconfig error: Cannot load default config file" and a segfault.
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        fontconfig \
        fonts-dejavu-core \
        zip \
    && rm -rf /var/lib/apt/lists/*
USER mamamia

COPY --chown=mamamia:mamamia VERSION LICENSE.md COPYRIGHTS README.md app.R ./
COPY --chown=mamamia:mamamia R ./R
COPY --chown=mamamia:mamamia src ./src
COPY --chown=mamamia:mamamia www ./www
COPY --chown=mamamia:mamamia data ./data
COPY --chown=mamamia:mamamia config ./config
COPY --chown=mamamia:mamamia deployment ./deployment
COPY --chown=mamamia:mamamia tools/build-native.R ./tools/build-native.R

# Compile the htslib/gcrank bridge and verify that the source loader, native
# library and cairo/ggplot stack all work before the image is tagged. The
# activation shim puts the uvr-managed R on PATH and R_LIBS_USER on
# .uvr/library.
RUN uvr activate --write-shim \
    && . .uvr/activate \
    && Rscript tools/build-native.R \
    && Rscript -e 'source("R/load.R"); mamamia_load("."); stopifnot(identical(mamamia_version(), readLines("VERSION", n = 1L))); cat("MaMaMIA", mamamia_version(), "ready\n")' \
    && Rscript -e 'library(ggplot2); f <- tempfile(fileext = ".png"); png(f, width = 400, height = 300, type = "cairo"); print(ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()); dev.off(); stopifnot(file.exists(f), file.info(f)$size > 0); cat("cairo text metrics OK\n")'

EXPOSE 3838

# ShinyProxy overrides this with the same command; keeping it in the image
# makes local smoke tests realistic.
CMD ["uvr", "run", "/app/deployment/start.R"]
