FROM balazsdukai/geoflow-bundle-builder:2023.07.17 as builder

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
RUN --mount=type=cache,target=/usr/src/tyler/target /root/.cargo/bin/cargo install --path .

COPY docker/strip-docker-image-export ./
RUN rm -rf /export
RUN mkdir /export && \
    bash ./strip-docker-image-export \
    -v \
    -d /export \
    -f /usr/local/share/proj/proj.db \
    -f /usr/local/bin/geof \
    -f $GF_PLUGIN_FOLDER/gfp_buildingreconstruction.so \
    -f $GF_PLUGIN_FOLDER/gfp_core_io.so \
    -f $GF_PLUGIN_FOLDER/gfp_gdal.so \
    -f $GF_PLUGIN_FOLDER/gfp_val3dity.so \
    -f $GF_PLUGIN_FOLDER/gfp_las.so \
    -f /root/.cargo/bin/tyler

FROM ubuntu:lunar-20230301
ARG VERSION
LABEL org.opencontainers.image.authors="Balázs Dukai <balazs.dukai@3dgi.nl>"
LABEL org.opencontainers.image.vendor="3DGI"
LABEL org.opencontainers.image.title="tyler"
LABEL org.opencontainers.image.description="Create tiles from 3D city objects encoded as CityJSONFeatures."
LABEL org.opencontainers.image.version=$VERSION
LABEL org.opencontainers.image.license="(APACHE-2.0 AND GPL-3 AND AGPL-3)"

ARG GF_PLUGIN_FOLDER="/usr/local/lib/geoflow-plugins"

RUN rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/src/tyler/resources/geof /usr/src/tyler/resources/geof
COPY --from=builder /usr/local/share/proj /usr/local/share/proj
COPY --from=builder $GF_PLUGIN_FOLDER $GF_PLUGIN_FOLDER
COPY --from=builder /export/lib/ /lib/
COPY --from=builder /export/lib64/ /lib64/
COPY --from=builder /export/usr/ /usr/
COPY --from=builder /export/root/.cargo/bin/tyler /usr/local/bin/tyler

# Update library links
RUN ldconfig
CMD ["tyler"]