FROM 3dgi/geoflow-bundle-builder:2024.08.09 AS builder

USER root

RUN apt-get update && apt-get install -y unzip curl postgresql

# Download Dutch transformation grids
RUN mkdir -p /usr/local/share/proj
RUN curl -o /usr/local/share/proj/nl_nsgi_nlgeo2018.tif https://cdn.proj.org/nl_nsgi_nlgeo2018.tif && \
    curl -o /usr/local/share/proj/nl_nsgi_rdtrans2018.tif https://cdn.proj.org/nl_nsgi_rdtrans2018.tif

# Needed for the proj-sys bindings
RUN apt-get install -y clang-15

# Install rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

WORKDIR /usr/src/tyler
COPY Cargo.toml Cargo.lock ./
COPY resources ./resources
COPY src ./src
COPY proj ./proj
RUN --mount=type=cache,target=/usr/src/tyler/target /root/.cargo/bin/cargo install --path . --bin tyler-db

COPY docker/strip-docker-image-export ./
RUN rm -rf /export
RUN mkdir /export && \
    bash ./strip-docker-image-export \
    -v \
    -d /export \
    -f /usr/local/share/proj/proj.db \
    -f /root/.cargo/bin/tyler-db

FROM ubuntu:lunar-20230301
ARG VERSION
LABEL org.opencontainers.image.authors="Balázs Dukai <balazs.dukai@3dgi.nl>"
LABEL org.opencontainers.image.vendor="3DGI"
LABEL org.opencontainers.image.title="tyler-db"
LABEL org.opencontainers.image.description=""
LABEL org.opencontainers.image.version=$VERSION
LABEL org.opencontainers.image.license="(APACHE-2.0 AND GPL-3 AND AGPL-3)"

RUN apt-get -y update && apt-get -y install curl postgresql && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/share/proj /usr/local/share/proj
COPY --from=builder /export/lib/ /lib/
COPY --from=builder /export/lib64/ /lib64/
COPY --from=builder /export/usr/ /usr/
COPY --from=builder /export/root/.cargo/bin/tyler-db /usr/local/bin/tyler-db

# Update library links
RUN ldconfig
CMD ["tyler-db"]