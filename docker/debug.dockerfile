FROM 3dgi/geoflow-bundle-builder:2025.09.01 AS base

USER root

ARG GF_PLUGIN_FOLDER="/usr/local/lib/geoflow-plugins"

RUN apt-get update && apt-get install -y unzip curl

# Download Dutch transformation grids
RUN wget https://cdn.proj.org/nl_nsgi_nlgeo2018.tif -O /usr/local/share/proj/nl_nsgi_nlgeo2018.tif && \
    wget https://cdn.proj.org/nl_nsgi_rdtrans2018.tif -O /usr/local/share/proj/nl_nsgi_rdtrans2018.tif

# Needed for the proj-sys bindings
RUN apt-get install -y clang-15

# Install rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

WORKDIR /usr/src/tyler

COPY Cargo.toml Cargo.lock ./
COPY resources ./resources
COPY src ./src
COPY proj ./proj

RUN --mount=type=cache,target=/usr/src/tyler/target-docker CARGO_TARGET_DIR=/usr/src/tyler/target-docker /root/.cargo/bin/cargo build --manifest-path ./Cargo.toml
