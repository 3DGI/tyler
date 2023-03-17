# Design document

The goal of Tyler is to create tiles for large areas from 3D city objects that are stored as CityJSON. 
In order to tile up large areas efficiently, it only loads the minimum required information from the 3D city objects for creating the tiles and keeps the remaining data on disk. 
For fast, city object access, Tyler relies on CityJSONFeature-s. 
CityJSONFeature is part of the CityJSON specification and it allows to store each feature (or city object) in a separate JSON document, or separate file. 
In case of Tyler,  each CityJSONFeature is stored in a separate file.

Tyler uses a square grid for indexing and counting the features. 
The extent of the grid is determined by visiting each feature. 
In case the features contain outliers, the vertical extent of the grid can be limited (`--grid-min/maxz`). 
The 2D cell size of the grid is an argument set by the user (`--grid-cellsize`).

After the grid, a quadtree is generated bottom-up, by merging each 4 grid cells in Morton-order (a quadrant) if the number of vertices in the quadrant are less than or equal to the configured capacity. 
The capacity can be set with the `--qtree-capacity` argument.
Thus, the grid cellsize (`--grid-cellsize`) determines the size of the smallest possible node in the quadtree, the capacity (`--qtree-capacity`) determines the maximum vertex capacity of the leaves of the quadtree.
However, if the grid cellsize is set too large in relation to the quadtree capacity, then the grid cells will not be merged, since they will already contain more vertices than the configured capacity. 
On the other hand, if the grid cellsize is too small, then the grid will contain too many cells and the tiling process will take significantly longer and consume more resources. 
The default cellsize is 250m and capacity is 42000 vertices.

The 3D Tiles tileset is created from the quadtree, by transforming the quadtree to the schema that is required by 3D Tiles. 
By default, an explicit tiling is created. 
It is possible to create implicit tiling with the `--3dtiles-implicit` argument.

## Assumptions set by Tyler

1. The feature coordinates are expected to be in a projected, Cartesian coordinate reference system with metric units.
2. If a feature has multiple children, it is expected that the children are adjacent, or at least very close to each other. This is because the whole feature, including all its children are treated as a unit in the tiling process, in order to make sure that all the children of a feature always end up in the same tile.

## Input data

1. A main `.city.json` file, containing at least the `CRS` and `transform` objects.
2. A directory (or directory tree) of `.city.jsonl` files, each containing one CityJSON Feature, including all its children City Objects.


## Misc. notes

Depending on the size of the *feature file path* string (UTF-8) in the `Feature`, the `FeatureSet` might use more or less memory.

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

*Optimization A: Parse the .city.jsonl files asynchronously.*

*Optimization B: Use a MapReduce paradigm (single machine) where we have parallel .city.jsonl --> feature_tuple mappers and their result is written into multiple feature_sets. The number of feature_sets depend on the level or parallelism, and a feature_set is identified by the first Nx2 bits of the morton code of a feature. Where N is the level of parallelism (more or less could work like this I think). Need to take care when updating a feature_set from multiple processes. Finally, merge (reduce) the multiple feature_sets into a single feature_set.*

#### Struct of arrays to store features

Instead of having an array of structs for the features, we should consider a struct of arrays.
And therefore benefit from cache-locality when looping over the features in order.
For this to work, the features probably should be in some spatial order.

Array of structs:

```rust
struct Feature {
    centroid_qc: [i64; 2],
    nr_vertices: u16,
    path_jsonl: PathBuf,
    bbox_qc: BboxQc,
}

type FeatureSet = Vec<Feature>;
```

Struct of arrays:

```rust
struct FeatureSet {
    centroid_qc: Vec<[i64; 2]>,
    nr_vertices: Vec<u16>,
    path_jsonl: Vec<PathBuf>,
    bbox_qc: Vec<BboxQc>,
}

type FeatureID = usize;

impl FeatureSet {
    fn get(&self, id: &FeatureID) -> (&[i64;2], &u16, &PathBuf, &BboxQc) {
        (&self.centroid_qc[*id], &self.nr_vertices[*id], &self.path_jsonl[*id], &self.bbox_qc[*id])
    }
}
```

#### Estimating memory use for a country like Netherlands

The Netherlands extent: `Polygon ((13565.3984375 306846.1875, 278026.125 306846.1875, 278026.125 619315.6875, 13565.3984375 619315.6875, 13565.3984375 306846.1875))`

y-side = 312469.5
x-side = 264460.7

Square grid: 312469.5 x 312469.5 meters.
Assuming a cellsize of 250m.
Nr. cells = `ceil(312469.5/250) = 12499^2 = 156225001`.

We store the square grid in 3D vector, where the 3rd dimension is the cell contents.
So we have one vector for the x-side, one vector for each y "row", and one vector for each cell that stores the pointers to the features.
The features are stored in a separate container and grid cell stores pointers to them.
So we have 1 + 521 + 156225001 vectors, plus number of features x `usize`.
An empty Vec is 128bit, and let's assume 10mio features which are all the buildings in the Netherlands.

Then a square grid with 250m cells, covering the Netherlands and indexing 10mio features is `(1 + 521 + 156225001) × 128 + 10000000 × 64` bits = 2.6Gb.

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
