# tyler

<p align="center">
  <img src="https://github.com/3DGI/tyler/blob/master/tyler.png" />
</p>

*tyler* creates tiles from 3D city objects.

As input, *tyler* takes [CityJSON Features](https://www.cityjson.org/specs/1.1.3/#text-sequences-and-streaming-with-cityjsonfeature), where each feature is stored in a separate file.

As output, *tyler* can create:

- 3D Tiles

## Usage

*tyler* is a command line application.

Use `--help` to see the help menu.

```shell
tyler --help
```

Generate 3D Tiles from a directory (or directory tree) of CityJSONFeatures.
We need information about the coordinate reference system and transformation properties of the feature vertices.
This information is stored in a main `.city.json` file, containing at least the [CRS](https://www.cityjson.org/specs/1.1.3/#referencesystem-crs) and [transform](https://www.cityjson.org/specs/1.1.3/#transform-object) objects.
Pass the main `.city.json` file as the `--metadata` argument.

```shell
tyler \
  --metadata metadata.city.json \
  --features directory_with_cityjsonfeatures \
  --output directory_for_output \
  --format 3dtiles
```

## Algorithm

### Assumptions

1. The feature coordinates are expected to be in a projected, Cartesian coordinate reference system with metric units.
2. The total number of vertices in a feature do not exceed *65535*, so that we can store the value in a `u16`.
3. The quantized coordinates fit in the range of an `i32`, that is `[-2147483648, 2147483647]`.
4. If a feature has multiple children, it is expected that the children are adjacent, or at least very close to each other. This is because the whole feature, including all its children are treated as a unit in the tiling process, in order to make sure that all the children of a feature always end up in the same tile.


### Inputs

1. A main `.city.json` file, containing at least the [CRS](https://www.cityjson.org/specs/1.1.3/#referencesystem-crs) and [transform](https://www.cityjson.org/specs/1.1.3/#transform-object) objects.
2. A directory (or directory tree) of `.city.jsonl` files, each containing one CityJSON Feature, including all its children City Objects.
3. (optional) A GeoPackage `.gpkg` with the `feature_tuple`s, where a `feature_tuple` is *(morton code of feature centroid, nr. vertices, feature file path)*.

### Process

Read the *CRS* and *transform* properties from the `.city.json` file and keep them in memory.
We need these properties for parsing the feature geometries.

If the `.gpkg` is not provided then create the `feature_set` as explained below. 
Else the `.gpkg` is the `feature_set`, and the `extent` (BBOX of all features) is computed from the `.gpkg`.
<-- There is one table/layer in the `.gpkg` that contains the `feature_tuple`s (see below). We don't really need to do any spatial operations directly in the `.gpkg`, thus we can use any SQLite library to read the records from the table. Eg. [sqlite3 in Python](https://docs.python.org/3.8/library/sqlite3.html), maybe [SQLiteCpp](https://srombauts.github.io/SQLiteCpp/), or the [native C SQLite library](https://www.sqlite.org/c3ref/intro.html), [sqlite in Rust](https://crates.io/crates/sqlite).

---

Creating the `feature_set` from the `.city.jsonl` files:

`feature_count` <-- Count the number of `.city.jsonl` files on the filesystem.
Initialize a vector `feature_set` with size `feature_count` to store the `feature_tuple`s.

A `feature_tuple` is *(morton code of feature centroid [`u64`], nr. vertices [`u16`], feature file path [`string`])*.
Therefore, for 10mio features the `feature_set` is expected to use *0.58–1.05*Gb memory (explained below).

Thus, the complete set of `feature_tuple`s contains all the information that is needed to construct a quadtree with all the features.
The morton code is 64-bit integer, as it is composed of a pair of 32-bit x,y coordinates (integers).
The integer coordinates are created by scaling the original floating point coordinates by 100, to retain centimeter precision.
By decoding the morton code and scaling the integer coordinates (0.01), we obtain the original coordinates of the feature centroid.
Thus, we only need to store the morton code of the centroid, instead of the coordinate pair.
Storing the morton code instead of the coordinates allows us to sort the features in z-order with conventional 1-D sorting algorithms.

Depending on the size of the *feature file path* string (UTF-8) in the `feature_tuple`, the `feature_set` might use more or less memory.

```shell
echo "b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl" | wc -m
# 49 characters

echo "/data/3DBAGv2/export/cityjson/v210908_fd2cee53/b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl" | wc -m
# 96 characters
```

Then in Python we have in bytes:

```python
assert(sys.getsizeof("b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl") == 97)

assert(sys.getsizeof("/data/3DBAGv2/export/cityjson/v210908_fd2cee53/b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl") == 144)
```

While in Rust (I guess C++ is same), in bytes:

```rust
assert_eq!(std::mem::size_of_val("b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl"), 48);

assert_eq!(std::mem::size_of_val("/data/3DBAGv2/export/cityjson/v210908_fd2cee53/b3bd7e17c-deb5-11e7-951f-610a7ca84980.city.jsonl"), 95);
```

Therefore, we expect that we can store the whole `feature_set` in memory.

Then for each `.city.jsonl` file, read sequentially: 
<-- superfast and low memory with Rust's serde, okay with Python's *json* or C++ *nlohmann*. Since we practically only need to parse the `vertices` array of the CityJSONFeature, we could write a tiny library in Rust and bind it to Python/C++ (if we use those languages).

1. Compute the feature centroid as the median of x,y in the `vertices` array of the CityJSONFeature.
2. Compute the morton code of the feature centroid.
3. Count the total number of vertices of the feature (simply the length of `vertices`, because shared vertices are irrelevant).
4. Get the file path.
5. Insert the `feature_tuple` into the `feature_set`.
6. Update the `extent` with the feature centroid.

*Optimization A: Parse the .city.jsonl files asynchronously.*

*Optimization B: Use a MapReduce paradigm (single machine) where we have parallel .city.jsonl --> feature_tuple mappers and their result is written into multiple feature_sets. The number of feature_sets depend on the level or parallelism, and a feature_set is identified by the first Nx2 bits of the morton code of a feature. Where N is the level of parallelism (more or less could work like this I think). Need to take care when updating a feature_set from multiple processes. Finally, merge (reduce) the multiple feature_sets into a single feature_set.*

---

*Do we need to sort the features? If so, we can sort in-memory if possible, or do [external sorting](https://en.wikipedia.org/wiki/External_sorting).*

Initialize the quadtree with the `extent`.
<-- Quadtree examples [in Rust 1](https://github.com/snorrwe/morton-table), [in Rust 2](https://docs.rs/quadtree_rs/0.1.2/quadtree_rs/), [in Rust 3](https://github.com/dmac/rust-quadtree), [in Rust 4](https://github.com/fschutt/quadtree-f32), [in Rust 5](https://docs.rs/crate/aabb-quadtree/0.2.0/source/src/lib.rs), [in Python 1](https://scipython.com/blog/quadtrees-2-implementation-in-python/), [in Python 2](https://github.com/CartoDB/QuadGrid), for sure there is a gazillion in C++.

Read sequentially from the `feature_set` and insert the `feature_tuple` ID into the quadtree.

The `feature_tuple` ID is the index of the `feature_tuple` in the `feature_set` (vector or table row).
I think we shouldn't use the morton code as the ID, because there could be multiple features with the same morton code, esp. since we have centimeter precision centroids of the complete feature (possibly a group of CityObjects).

There are multiple possible split criteria for the quadtree:

1. fixed cellsize: We recursively subdivide the extent until we reach the approximate cellsize. So, basically the split criteria is the edge length of the cell. Then some leafs will have many features in them, other leafs will have just a couple, others will have zero, but we will have a tile hierarchy still. To build a tileset with this type of subdivision we only need to know the morton code of the centroids. This spatial window query can be done efficiently with the morton codes.
2. max. nr. of features (aka 3D BAG)
3. max. nr. of vertices

The quadtree is built from bottom up. 
Eg. define a minimum tree depth/cellsize, then assign all the centroids to a node on that depth (simple grid indexing super fast). 
Then merge nodes that contain too few data (4 nodes -> 1 node one level shallower in the tree). 

Thus, the quadtree is constructed sequentially, in memory.
Its leafs contain a set of `feature_tuple` IDs.

*Optimization: Prune the quadtree to remove empty nodes. Although, I'm not sure how this works.*

Create the `tileset.json` from the quadtree.
Where at least for each leaf a path to a `.glb` file is generated and added to the `tileset.json`.
But probably we need `.glb` assets for higher level nodes too, that contain less detailed models.

For each node/leaf, done in parallel:

1. merge the features into a citymodel
2. reproject to ECEF crs (required by Cesium)
3. convert to gltf with cjio
4. write the gltf to the `.glb` path that is set in `tileset.json`

#### Struct of arrays to store features

Instead of having an array of structs for the features, we should consider a struct of arrays.
And therefore benefit from cache-locality when looping over the features in order.
For this to work, the features probably should be in some spatial order.

Array of structs:

```rust
struct Feature {
    centroid_quantized: [i32;2],
    nr_vertices: u16,
    path_jsonl: PathBuf
}

type FeatureSet = Vec<Feature>;
```

Struct of arrays:

```rust
struct FeatureSet {
    centroid_quantized: Vec<[i32;2]>,
    nr_vertices: Vec<u16>,
    path_jsonl: Vec<PathBuf>
}

type FeatureID = usize;

impl FeatureSet {
    fn get(&self, id: &FeatureID) -> (&[i32;2], &u16, &PathBuf) {
        (&self.centroid_quantized[*id], &self.nr_vertices[*id], &self.path_jsonl[*id])
    }
}
```

#### Estimating memory use for a country like Netherlands

The Netherlands extent: `Polygon ((13565.3984375 306846.1875, 278026.125 306846.1875, 278026.125 619315.6875, 13565.3984375 619315.6875, 13565.3984375 306846.1875))`

y-side = 312469.5
x-side = 264460.7

Square grid: 312469.5 x 312469.5 meters.
Assuming a cellsize of 600m (~600m is the smallest tile in the 3D BAG).
Nr. cells = `ceil(312469.5/600) = 521^2 = 271441`.

We store the square grid in 3D vector, where the 3rd dimension is the cell contents.
So we have one vector for the x-side, one vector for each y "column", and one vector for each cell that stores the pointers to the features.
The features are stored in a separate container and grid cell stores pointers to them.
So we have 1 + 521 + 271441 vectors, plus number of features x `usize`.
An empty Vec is 128bit, and let's assume 10mio features which are all the buildings in the Netherlands.

Then a square grid with 600m cells, covering the Netherlands and indexing 10mio features is `(1 + 521 + 271441) × 128 + 10000000 × 64` bits = 84.35Mb.

#### Merging CityJSONFeatures in Python and writing gltf

A simple cjio-based script merges a list of `.city.jsonl` files and exports a `.glb` from it.

```python
import json
from sys import argv
from pathlib import Path

from cjio.cityjson import CityJSON


def merge(cityjson_path, features_dir):
    lcount = 1
    #-- read first line
    with Path(cityjson_path).resolve().open("r") as fo:
        j1 = json.load(fo)
    cm = CityJSON(j=j1)
    if "CityObjects" not in cm.j:
        cm.j["CityObjects"] = {}
    if "vertices" not in cm.j:
        cm.j["vertices"] = []
    for child in Path(features_dir).iterdir():
        if child.suffix == ".jsonl":
            with child.open("r") as fo:
                j1 = json.load(fo)
            if not( "type" in j1 and j1["type"] == 'CityJSONFeature'):
               raise IOError("Line {} is not of type 'CityJSONFeature'.".format(lcount))
            cm.add_cityjsonfeature(j1)
    return cm

if __name__ == "__main__":
    cm = merge(argv[1], argv[2])
    glb = cm.export2glb()
    glb.seek(0)
    with Path(argv[3]).open("wb") as bo:
        bo.write(glb.getvalue())
```

Testing on *3dbag_v21031_7425c21b_5910.json*, which is 2.7Mb on disk, triangulated faces, 436 CityObjects, ~30k vertices.

`/usr/bin/time -v` reports 56 seconds, 127Mb memory use.

```
User time (seconds): 56.34
System time (seconds): 0.64
Percent of CPU this job got: 101%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:56.08
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 127740
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 38536
Voluntary context switches: 20
Involuntary context switches: 977
Swaps: 0
File system inputs: 0
File system outputs: 9376
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```
