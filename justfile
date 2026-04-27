#!/usr/bin/env just --justfile

_default:
    just --list

# Fast compile-time validation across the workspace.
check:
    cargo check --workspace --all-targets --all-features

# Build all workspace targets with all features enabled.
build *args:
    cargo build --workspace --all-targets --all-features {{args}}

# Run clippy with strict lint settings across the workspace.
lint:
    RUSTFLAGS='-Dclippy::all -Dclippy::pedantic' RUSTC_WORKSPACE_WRAPPER="$(command -v clippy-driver)" cargo check --workspace --all-targets --all-features

# Format the workspace with rustfmt.
fmt:
    cargo fmt --package tyler --package cityjson-convert

# Check workspace formatting without rewriting files.
fmt-check:
    cargo fmt --package tyler --package cityjson-convert --check

# Run the workspace tests with all features enabled.
test:
    cargo test --workspace --all-targets --all-features

# Clean the workspace by removing all build artifacts and test artifacts.
clean: clean-output
    cargo clean

# Clean the test output directories.
clean-output:
    rm -rf tests/output*
    rm -rf cityjson-convert/tests/output

# Run the full local validation sequence.
ci:
    just fmt
    just lint
    just check
    just build
    just test

# Run the full validation sequence without modifying files.
ci-check:
    just fmt-check
    just lint
    just check
    just build
    just test
