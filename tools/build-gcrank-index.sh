#!/bin/sh
# Build a gcrank index for a (optionally gzipped) FASTA file.
#
# Usage:
#   tools/build-gcrank-index.sh <fasta> <out-dir> [index-name]
#
# The builder is looked up in this order:
#   1. $GCRANK_BUILD, if set;
#   2. gcrank-build on PATH;
#   3. test_data/gcidx/gcrank-build, when running from a source checkout.
set -eu

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <fasta> <out-dir> [index-name]" >&2
    exit 1
fi

fasta=$1
out_dir=$2
index_name=${3:-}

if [ ! -r "$fasta" ]; then
    echo "FASTA file is not readable: $fasta" >&2
    exit 1
fi

mkdir -p "$out_dir"

if [ -n "${GCRANK_BUILD:-}" ]; then
    builder=$GCRANK_BUILD
elif command -v gcrank-build >/dev/null 2>&1; then
    builder=gcrank-build
elif [ -x "test_data/gcidx/gcrank-build" ]; then
    builder=test_data/gcidx/gcrank-build
else
    echo "gcrank-build not found; set GCRANK_BUILD or add it to PATH" >&2
    exit 1
fi

if [ -n "$index_name" ]; then
    exec "$builder" --fasta "$fasta" --out "$out_dir" --name "$index_name"
fi
exec "$builder" --fasta "$fasta" --out "$out_dir"
