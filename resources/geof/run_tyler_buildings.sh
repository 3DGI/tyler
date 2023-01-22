rm -rf /dev/shm/3dtiles-buildings
mkdir /dev/shm/3dtiles-buildings

CARGO_MANIFEST_DIR=/home/ravi/git/tyler \
RUST_LOG=debug \
/usr/bin/time -v \
/home/ravi/git/tyler/target/release/tyler \
        --metadata /mnt/Data/LocalData/Kadaster/3dbasisvoorziening_dh/features/metadata.city.json \
        --features /mnt/Data/LocalData/Kadaster/3dbasisvoorziening_dh/features/ \
        --output /dev/shm/3dtiles-buildings \
        --python-bin "/home/ravi/git/tyler/resources/geof/run_geof_buildings.sh" \
        --format 3dtiles \
        --cellsize 500 \
        --quadtree-criteria vertices \
        --quadtree-limit 42000 \
	--object-type building \
	--object-type building-part \
        --limit-minz=-5 \
        --limit-maxz=300 \
        --gltfpack-bin "/home/ravi/git/meshoptimizer/build/gltfpack" \
 &> /dev/shm/3dtiles-buildings/tyler_building_"$(date +"%Y-%m-%d")".log \
 rm -rf /dev/shm/3dtiles-buildings/inputs
