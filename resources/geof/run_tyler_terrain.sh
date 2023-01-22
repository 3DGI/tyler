rm -rf /dev/shm/3dtiles-terrain
mkdir /dev/shm/3dtiles-terrain

CARGO_MANIFEST_DIR=/home/ravi/git/tyler \
RUST_LOG=debug \
/usr/bin/time -v \
/home/ravi/git/tyler/target/release/tyler \
        --metadata /mnt/Data/LocalData/Kadaster/3dbasisvoorziening_dh/features/metadata.city.json \
        --features /mnt/Data/LocalData/Kadaster/3dbasisvoorziening_dh/features/ \
        --output /dev/shm/3dtiles-terrain \
        --python-bin "/home/ravi/git/tyler/resources/geof/run_geof_terrain.sh" \
        --format 3dtiles \
        --cellsize 400 \
        --quadtree-criteria vertices \
        --quadtree-limit 42000 \
	--object-type land-use \
	--object-type plant-cover \
	--object-type water-body \
	--object-type road \
	--object-type generic-city-object \
	--object-type bridge \
        --limit-minz=-5 \
        --limit-maxz=300 \
        --gltfpack-bin "/home/ravi/git/meshoptimizer/build/gltfpack" \
 &> /dev/shm/3dtiles-terrain/tyler_terrain_"$(date +"%Y-%m-%d")".log \
 rm -rf /dev/shm/3dtiles-terrain/inputs
