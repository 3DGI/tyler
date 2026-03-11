#!/usr/bin/env bash
set -euo pipefail

# Build a fully-static Tyler binary for x86_64-linux (musl).
#
# Prerequisites (install via your package manager or the Dockerfile):
#   - rustup target add x86_64-unknown-linux-musl
#   - musl-tools (provides musl-gcc)
#   - cmake, make
#
# Usage:
#   ./build-musl.sh          # release build
#   ./build-musl.sh --debug  # debug build (faster compile)

TARGET="x86_64-unknown-linux-musl"
PROFILE="release"

if [[ "${1:-}" == "--debug" ]]; then
    PROFILE="dev"
fi

# Ensure musl target is installed.
rustup target add "$TARGET" 2>/dev/null || true

# Set linker for musl.
export CC_x86_64_unknown_linux_musl="musl-gcc"
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER="musl-gcc"

# Build.
if [[ "$PROFILE" == "release" ]]; then
    cargo build --release --target "$TARGET"
    BINARY="target/${TARGET}/release/tyler-glb"
else
    cargo build --target "$TARGET"
    BINARY="target/${TARGET}/debug/tyler-glb"
fi

# Verify static linking.
if command -v file &>/dev/null; then
    echo ""
    file "$BINARY"
fi
if command -v ldd &>/dev/null; then
    echo ""
    ldd "$BINARY" 2>&1 || true
fi

echo ""
echo "Binary: $BINARY"
ls -lh "$BINARY"
