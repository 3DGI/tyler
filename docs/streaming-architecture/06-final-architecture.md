# 06 — Final Architecture

This document describes the unified final architecture that combines
all accepted decisions:

- Pre-declared extent + per-cell spill files (from [02](./02-alternatives-considered.md))
- Progressive per-feature conversion + progressive assembly (from [03](./03-progressive-tile-generation.md))
- Morton-ordered input + low-water mark (from [04](./04-morton-ordered-input.md))
- Native Rust geof library (from [05](./05-geof-library.md))

## System topology

```
┌──────────┐    CityJSONSeq     ┌─────────────┐    morton-ordered    ┌─────────────────────┐
│  roofer  │ ─────stdout─────▶ │ morton-sort │ ───CityJSONSeq─────▶ │  tyler (--stdin)    │
└──────────┘                    └─────────────┘                      │                     │
                                                                     │  ┌─────────────┐   │
                                                                     │  │ geof crate  │◀──┤ linked library
                                                                     │  └─────────────┘   │
                                                                     │                     │
                                                                     │  ┌───────────────┐ │
                                                                     │  │ ingestion +   │ │
                                                                     │  │ progressive   │ │
                                                                     │  │ assembly      │ │
                                                                     │  └───────────────┘ │
                                                                     │                     │
                                                                     │  output/            │
                                                                     │    t/*/*/*.glb      │
                                                                     │    tileset.json     │
                                                                     └─────────────────────┘

Separately:  geof input.city.jsonl output.glb      ← standalone CLI
             (same library, one-shot tile mode)
```

## Components

### 1. The producer (roofer or any CityJSONSeq source)

Emits CityJSONSeq to stdout. Line 0 is the CityJSON metadata header;
subsequent lines are CityJSONFeature documents (one per line) or
optional `+`-prefixed control lines.

**Required:** valid CityJSONSeq with CRS and transform in the header.

**Preferred:** the header should include `metadata.geographicalExtent`
so `morton-sort` and tyler can derive the extent without a CLI argument.

**Optional:** approximate spatial coherence in output order (reduces
`morton-sort` buffer requirements) and/or `+RegionComplete` control
lines (enables explicit region finalization without morton ordering).

### 2. The `morton-sort` tool

Standalone binary. Reads CityJSONSeq from stdin, writes CityJSONSeq to
stdout in approximate morton order. Detailed in
[04-morton-ordered-input](./04-morton-ordered-input.md).

**Skippable** if the producer emits in morton order natively:

```
roofer --morton-output | tyler --stdin
```

### 3. The geof Rust crate

Linked into tyler (no subprocess). Also built as a standalone CLI for
ad-hoc conversions. Detailed in [05-geof-library](./05-geof-library.md).

### 4. Tyler (`--stdin` mode)

Consumes the morton-ordered CityJSONSeq stream, builds tiles
progressively, emits output files continuously. Detailed below.

## Tyler's internals in streaming mode

### Initialization

1. Parse CLI arguments: `--stdin`, `--format`, `--extent`, output
   directory, grid cellsize, quadtree capacity, etc.
2. Read line 0 from stdin (the CityJSONSeq header). Extract CRS and
   transform. Validate.
3. Determine the spatial extent:
   - From `--extent` CLI argument if provided.
   - Otherwise from the header's `metadata.geographicalExtent`.
   - Otherwise: fall back to two-pass mode (spill to temp file, scan for
     extent, then process the temp file via the same streaming loop).
4. Allocate the `SquareGrid` from the extent (unchanged from current
   code).
5. Create the output directory layout (`output/t/`, `output/spill/`).
6. Initialize a Rayon thread pool for assembly workers.
7. Initialize a live quadtree tracker (see below).

### Data structures

**In-memory:**

```rust
struct StreamingWorld {
    grid: SquareGrid,                    // unchanged
    features: Vec<Feature>,              // 84 bytes × N; reference = SpillRef
    live_quadtree: LiveQuadtree,         // sealing state tracker
    low_water_mark: u128,                // max morton code seen
    cell_spill_writers: HashMap<CellId, BufWriter<File>>,
    format: Format,
    geof_encoder: Box<dyn FeatureEncoder<Intermediate = ...>>,
}

pub struct Feature {
    pub centroid: [f64; 2],
    pub bbox: [f64; 6],
    pub reference: FeatureReference,
}

pub enum FeatureReference {
    LegacyPath(PathBuf),
    CjIndexId(String),
    SpillRef { cell_id: CellId, byte_offset: u64, byte_len: u32 },
}
```

**On-disk:**

```
output/
├── spill/                          # per-cell intermediates (temporary)
│   ├── 0_0.bin                     # binary (format-specific intermediates)
│   ├── 0_1.bin
│   └── ...
├── t/                              # finalized tile outputs
│   ├── 0/0/0.glb
│   ├── 0/0/1.glb
│   └── ...
└── tileset.json                    # written at EOF
```

### The live quadtree

A lightweight mutable structure tracking the sealing state. Can be
implemented as a sparse tree (only nodes touched by features are
allocated) or dense (arrays indexed by level + morton code).

```rust
struct LiveQuadtree {
    levels: Vec<LiveLevel>,
}

struct LiveLevel {
    /// For each node at this level: vertex count, sealed flag.
    nodes: HashMap<NodeId, LiveNode>,
}

struct LiveNode {
    nr_vertices: usize,
    sealed: bool,                    // sum > capacity → won't merge
    siblings_complete: u8,           // 0..4; how many of the 4 siblings
                                     // have been fully drained by LWM
}

impl LiveQuadtree {
    fn on_feature_added(&mut self, cell_id: CellId, vtx: usize) -> Vec<LeafTileId> {
        // Walk up the tree, updating counts, sealing groups that cross
        // the capacity threshold, and returning any newly-finalized
        // leaf tile IDs.
    }

    fn on_lwm_advanced(&mut self, lwm: u128) -> Vec<LeafTileId> {
        // Mark cells whose morton code < lwm as drained; check whether
        // their quadrants are now complete; recursively advance the
        // completion state up the tree. Return newly-finalized leaves.
    }
}
```

### The ingestion loop

```rust
fn streaming_loop(
    stdin: impl BufRead,
    world: &mut StreamingWorld,
    assembly_tx: &Sender<AssemblyTask>,
) -> Result<()> {
    for line in stdin.lines() {
        let line = line?;
        let bytes = line.as_bytes();

        // Control line?
        if bytes.starts_with(b"{\"type\":\"+") {
            handle_control_line(bytes, world, assembly_tx)?;
            continue;
        }

        // Parse feature
        let feature = cjlib::parse_feature(bytes)?;

        // Extract stats, assign to grid cell(s)
        let stats = selected_geometry_stats(&feature, world.cityobject_types.as_ref());
        let Some(bbox) = stats.bbox else { continue; };
        let centroid = stats.centroid.unwrap();
        let cell_vtx_cnt = count_vertices_in_grid(&feature, &stats.selected_vertices, &world.grid, &bbox);
        let cell_assignment = feature_to_cells(&cell_vtx_cnt, &stats.selected_object_types);

        // Encode feature per the selected format (triangulate for GLB,
        // passthrough for CityJSON, etc.)
        let intermediate = world.geof_encoder.encode_feature(&feature, &world.options)?;

        // Append to cell spill file(s)
        for (cell_id, vtx_in_cell) in &cell_assignment {
            let writer = world.cell_spill_writers.entry(*cell_id)
                .or_insert_with(|| open_cell_spill(&world.output_dir, *cell_id));
            let (byte_offset, byte_len) = write_intermediate(writer, &intermediate)?;
            // (For unique-assignment types like Building, only write once)

            // Update grid
            let cell = world.grid.cell_mut(cell_id);
            cell.feature_ids.push(world.features.len());
            cell.nr_vertices += vtx_in_cell;

            // Update live quadtree
            let newly_finalized = world.live_quadtree.on_feature_added(*cell_id, *vtx_in_cell);
            for leaf_id in newly_finalized {
                // Leaf is sealed, but its feature set may still grow.
                // Actual dispatch to assembly happens when LWM passes
                // the last cell in the leaf.
            }
        }

        // Store Feature (84 bytes)
        world.features.push(Feature {
            centroid,
            bbox,
            reference: FeatureReference::SpillRef { /* ... */ },
        });

        // Advance low-water mark
        let morton = morton_of(&centroid, &world.grid);
        if morton > world.low_water_mark {
            let newly_drained = world.live_quadtree.on_lwm_advanced(morton);
            for leaf_id in newly_drained {
                assembly_tx.send(AssemblyTask::Tile(leaf_id))?;
            }
            world.low_water_mark = morton;
        }
    }

    // EOF: finalize all outstanding work
    let remaining = world.live_quadtree.finalize_all();
    for leaf_id in remaining {
        assembly_tx.send(AssemblyTask::Tile(leaf_id))?;
    }
    drop(assembly_tx); // signal workers to finish
    Ok(())
}
```

### The assembly workers

A Rayon thread pool consumes `AssemblyTask`s from the channel:

```rust
fn assembly_worker(
    rx: Receiver<AssemblyTask>,
    world: Arc<StreamingWorld>,
) -> Result<()> {
    while let Ok(task) = rx.recv() {
        match task {
            AssemblyTask::Tile(leaf_id) => {
                let tile_cells = world.live_quadtree.cells_of_leaf(&leaf_id);
                let intermediates = read_intermediates(&world.output_dir, tile_cells)?;
                let output_path = tile_output_path(&world.output_dir, &leaf_id);
                let mut output = BufWriter::new(File::create(&output_path)?);
                world.geof_encoder.assemble_tile(
                    intermediates,
                    &tile_metadata_for(leaf_id, &world),
                    &mut output,
                    &world.options,
                )?;
            }
        }
    }
    Ok(())
}
```

### EOF handling

1. Flush any `BufWriter`s on cell spill files.
2. For each quadrant that hasn't been sealed and isn't yet complete
   (sparse regions): compute the final merge decision with complete
   vertex counts.
3. Assemble all remaining leaves.
4. Build the final `Tileset` structure from the completed quadtree.
5. Write `tileset.json`.
6. (Optionally) delete the `output/spill/` directory.

## Memory profile

| Data | Size | Notes |
|---|---|---|
| Reorder buffer (in `morton-sort`) | ~20MB | 10K × 2KB avg feature |
| Grid (`Vec<Vec<Cell>>`) | O(cells) | Unchanged from current tyler |
| Feature set | N × 84 bytes | Lightweight index; ~840MB for 10M features |
| Live quadtree | O(non-empty nodes) × small | Sparse; small constant per node |
| Cell spill write buffers | O(active cells) × 64KB | Typical `BufWriter` default |
| Geometry data | **0 bytes in RAM** | On disk in per-cell spill files |

Total memory for a Netherlands-scale run: on the order of 1–2 GB,
bounded by the grid + feature set, independent of the total geometry
volume.

## Disk profile

- Per-cell spill files total approximately the size of the input
  (format-specific intermediates may be slightly more or less compact
  than raw CityJSON).
- Final output files are the 3D Tiles / OBJ / CityJSON / GPKG output.
- Spill files can be deleted after EOF if tiles are fully assembled.

## Progressive behavior

- **Dense urban tiles** (cells exceeding capacity individually): sealed
  on cell threshold crossing; assembled as soon as the LWM passes.
  Typically completed well before EOF.
- **Moderate-density tiles**: sealed when quadrant sum crosses capacity
  OR when LWM completes the quadrant. Assembled progressively.
- **Sparse rural tiles**: quadrant sums may stay below capacity until
  EOF. Merge decision and assembly deferred to EOF.
- **`tileset.json`**: always written at EOF (requires complete tree
  structure). Writing is fast (small JSON file).

For the typical Netherlands-scale dataset, the majority of assembly
work (the expensive triangulation + packing) happens concurrently with
ingestion. EOF work is limited to a small fraction of sparse tiles plus
the tileset JSON.

## Fallback modes

The architecture supports graceful degradation:

1. **No morton-sort in pipe.** Tyler still works. LWM advances
   non-monotonically (with warnings on out-of-order features); sealed
   leaves are still assembled progressively; sparse regions are
   handled at EOF. The progressive benefit is reduced but correctness
   is preserved.
2. **No extent provided.** Tyler falls back to two-pass mode: spills
   all features to a temp file, scans for extent, then runs the
   streaming loop against the temp file.
3. **geof library integration failure.** The crate also ships as a CLI;
   tyler can invoke it per-tile as a subprocess (matching the current
   architecture) as an emergency fallback. This is not the designed-for
   path.

## Backward compatibility

The existing tyler invocation modes (legacy feature-file directories,
cjindex datasets) continue to work. Streaming is an additional input
mode (`--stdin`), not a replacement. The existing grid, quadtree,
tileset, and tile-export code paths are reused with minimal changes
(primarily: the `FeatureReference` enum gains a `SpillRef` variant;
`write_inputs` learns to read from spill files instead of or in
addition to legacy/cjindex sources).
