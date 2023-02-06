FROM balazsdukai/geoflow-bundle-builder-ubuntu:latest as builder

ARG GF_PLUGIN_FOLDER="/usr/local/lib/geoflow-plugins"

RUN apt-get install unzip curl
RUN wget https://github.com/zeux/meshoptimizer/releases/download/v0.18/gltfpack-ubuntu.zip -O /tmp/gltfpack-ubuntu.zip
RUN unzip /tmp/gltfpack-ubuntu.zip -d /tmp && \
    chmod a+x /tmp/gltfpack

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

FROM ubuntu:lunar-20221216
ARG VERSION
LABEL org.opencontainers.image.authors="Balázs Dukai <balazs.dukai@3dgi.nl>"
LABEL org.opencontainers.image.vendor="3DGI"
LABEL org.opencontainers.image.title="tyler"
LABEL org.opencontainers.image.description="Create tiles from 3D city objects encoded as CityJSONFeatures."
LABEL org.opencontainers.image.version=$VERSION

ARG GF_PLUGIN_FOLDER="/usr/local/lib/geoflow-plugins"

RUN rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/src/tyler/resources/geof /usr/src/tyler/resources/geof
COPY --from=builder /tmp/gltfpack /usr/local/bin/gltfpack
COPY --from=builder /usr/local/share/proj /usr/local/share/proj
COPY --from=builder $GF_PLUGIN_FOLDER $GF_PLUGIN_FOLDER
COPY --from=builder /export/lib/ /lib/
COPY --from=builder /export/lib64/ /lib64/
COPY --from=builder /export/usr/ /usr/
COPY --from=builder /export/root/.cargo/bin/tyler /usr/local/bin/tyler

# Update library links
RUN ldconfig
CMD ["tyler"]