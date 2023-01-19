rm -rf /dev/shm/3dtiles
mkdir /dev/shm/3dtiles

CARGO_MANIFEST_DIR=/home/ravi/git/tyler \
RUST_LOG=debug \
# /usr/bin/time -v \
/home/ravi/git/tyler/target/release/tyler \
        --metadata /mnt/Data/LocalData/Kadaster/db3dnl_features/metadata.city.json \
        --features /mnt/Data/LocalData/Kadaster/db3dnl_features/gb2 \
        --output /dev/shm/3dtiles \
        --python-bin "/home/ravi/git/tyler/resources/geof/run_geof.sh" \
        --format 3dtiles \
        --cellsize 500 \
        --quadtree-criteria vertices \
        --quadtree-limit 18000 \
        --object-type building \
        --object-type building-part \
        --limit-minz=-5 \
        --limit-maxz=300 \
        --gltfpack-bin "/home/ravi/git/meshoptimizer/build/gltfpack"
