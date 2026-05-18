#!/usr/bin/env just --justfile

_default:
    just --list

# Fast compile-time validation across the workspace.
check:
    cargo check --workspace --all-targets --no-default-features --features proj-system

# Build all workspace targets with system-preferred PROJ enabled.
build *args:
    cargo build --workspace --all-targets --no-default-features --features proj-system {{args}}

# Run a perf + Massif profiling session for tyler.
profile *args:
    @bash scripts/profile-tyler.sh {{args}}

# Run clippy with strict lint settings across the workspace.
lint:
    RUSTFLAGS='-Dclippy::all -Dclippy::pedantic' RUSTC_WORKSPACE_WRAPPER="$(command -v clippy-driver)" cargo check --workspace --all-targets --no-default-features --features proj-system

# Format the workspace with rustfmt.
fmt:
    cargo fmt --package tyler --package cityjson-convert

# Check workspace formatting without rewriting files.
fmt-check:
    cargo fmt --package tyler --package cityjson-convert --check

# Run the workspace tests with system-preferred PROJ enabled.
test:
    cargo test --workspace --all-targets --no-default-features --features proj-system

# Collect test coverage with cargo-tarpaulin.
coverage:
    cargo tarpaulin --workspace --all-targets --no-default-features --features proj-system --out Stdout --out Xml

# Clean the workspace by removing all build artifacts and test artifacts.
clean: clean-output
    cargo clean

# Clean the test output directories.
clean-output:
    rm -rf tests/output*
    rm -rf cityjson-convert/tests/output

# Run the full local validation sequence.
ci: fmt lint check build test

# Run the full validation sequence without modifying files.
ci-check: fmt-check lint check build test

# Run full validation using the system-preferred PROJ mode.
ci-check-system: fmt-check lint check build test

# Validate bundled PROJ source builds without running the full test suite.
ci-check-bundled: fmt-check
    cargo check --workspace --all-targets --no-default-features --features proj-bundled
    cargo build --workspace --all-targets --no-default-features --features proj-bundled
