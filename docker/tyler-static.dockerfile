# Builder-only Dockerfile for producing a static tyler-glb binary.
#
# Uses Alpine (musl-native) to avoid glibc/musl ABI mismatches.
# All C/C++ code is compiled natively against musl — no cross-compilation.
# SQLite is compiled from source automatically by libsqlite3-sys (bundled
# feature).  A small Debian stage provides proj.db (Alpine has no proj-data
# package).
#
# Usage:
#   docker build --output type=local,dest=. -f docker/tyler-static.dockerfile .
#
# This writes ./tyler-glb directly to the current directory.
#
# The resulting ./tyler-glb binary is fully static (musl) and embeds proj.db.
# Copy it to your HPC node and run directly — no containers needed.

# ── Stage 1: grab proj.db from Debian ──
FROM debian:bookworm-slim AS projdb
RUN apt-get update && apt-get install -y --no-install-recommends proj-data \
    && rm -rf /var/lib/apt/lists/*

# ── Stage 2: Alpine build ──
FROM rust:1.93-alpine AS builder

# Install corporate CA certificates so TLS works behind corporate proxies.
# To use: place .crt/.pem files in the certs/ directory at the repo root.
# If certs/ is empty (default), this step is a no-op and the build uses
# the standard Alpine CA bundle.
COPY certs/ /usr/local/share/ca-certificates/
RUN apk add --no-cache ca-certificates && update-ca-certificates

RUN apk add --no-cache \
    build-base \
    cmake \
    clang-dev \
    pkgconf \
    sqlite

# On Alpine, x86_64-unknown-linux-musl is the native target — no cross-compilation.
# Skip bindgen (which needs dlopen for libclang, incompatible with musl static
# builds) and use pre-generated bindings instead.
ENV PROJ_SYS_SKIP_BINDGEN=1

# Copy proj.db from the Debian stage for embedding.
COPY --from=projdb /usr/share/proj/proj.db /usr/share/proj/proj.db
ENV PROJ_DB_PATH=/usr/share/proj/proj.db

WORKDIR /usr/src/tyler
COPY Cargo.toml Cargo.lock ./
COPY build.rs ./
COPY resources ./resources
COPY src ./src
COPY proj ./proj

RUN cargo build --release --target x86_64-unknown-linux-musl && \
    cp target/x86_64-unknown-linux-musl/release/tyler-glb /usr/local/bin/tyler-glb

# Verify the binary is static.
RUN file /usr/local/bin/tyler-glb && ldd /usr/local/bin/tyler-glb 2>&1 || true

# Minimal final stage — just the binary, for easy extraction.
FROM scratch
COPY --from=builder /usr/local/bin/tyler-glb /tyler-glb
