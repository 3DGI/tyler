FROM mcr.microsoft.com/devcontainers/rust:2-1-trixie

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        clang \
        cmake \
        curl \
        git \
        just \
        libproj-dev \
        libsqlite3-dev \
        libtiff-dev \
        pkg-config \
        sqlite3 \
    && rm -rf /var/lib/apt/lists/*
