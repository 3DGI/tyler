# Builder-only Dockerfile for producing a static tyler-glb binary (musl).
#
# Cross-compiles from Debian to x86_64-unknown-linux-musl using clang.
# NOTE: This Dockerfile has known issues with glibc/musl C++ header conflicts.
# Prefer tyler-static.dockerfile (Alpine native) for reliable builds.
#
# Usage:
#   docker build --output type=local,dest=. -f docker/tyler-glb.dockerfile .
#   # Or to use corporate CA certs:
#   # Place .crt/.pem files in certs/ directory, then build.

FROM rust:1.88-bookworm AS builder

# ── 1. System packages ──
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake clang libclang-dev lld g++ ca-certificates \
    musl-tools musl-dev \
    proj-data sqlite3 \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Corporate CA certificates (optional) ──
# Place .crt/.pem files in certs/ at the repo root. If empty, this is a no-op.
COPY certs/ /usr/local/share/ca-certificates/
RUN update-ca-certificates
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt

# ── 3. Rust musl target ──
RUN rustup target add x86_64-unknown-linux-musl

# ── 4. Musl cross-compilation environment ──
ENV TARGET=x86_64-unknown-linux-musl \
    CC_x86_64_unknown_linux_musl=clang \
    CXX_x86_64_unknown_linux_musl=clang++ \
    AR_x86_64_unknown_linux_musl=ar \
    CFLAGS_x86_64_unknown_linux_musl="--target=x86_64-unknown-linux-musl -I/usr/include/x86_64-linux-musl" \
    CXXFLAGS_x86_64_unknown_linux_musl="--target=x86_64-unknown-linux-musl -I/usr/include/x86_64-linux-musl" \
    CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=clang \
    RUSTFLAGS="-C linker=clang -C link-arg=--target=x86_64-unknown-linux-musl -C link-arg=-fuse-ld=lld" \
    PROJ_SYS_SKIP_BINDGEN=1 \
    PROJ_DB_PATH=/usr/share/proj/proj.db

WORKDIR /usr/src/tyler

# Copy source
COPY Cargo.toml Cargo.lock build.rs ./
COPY resources ./resources
COPY src ./src
COPY proj ./proj

# ── 5. Build ──
RUN --mount=type=cache,target=/usr/src/tyler/target-docker \
    CARGO_TARGET_DIR=/usr/src/tyler/target-docker \
    cargo build --release --target x86_64-unknown-linux-musl \
    && cp /usr/src/tyler/target-docker/x86_64-unknown-linux-musl/release/tyler-glb \
          /usr/local/bin/tyler-glb

# Verify static linking
RUN file /usr/local/bin/tyler-glb && ldd /usr/local/bin/tyler-glb 2>&1 || true

# ── Minimal final stage — just the binary ──
FROM scratch
COPY --from=builder /usr/local/bin/tyler-glb /tyler-glb
