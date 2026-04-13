# 05 — Native Rust geof Library and CLI

This document specifies the geof crate that replaces the geof
subprocess. The crate exposes both a high-level tile-conversion API
(used by the standalone CLI and for simple cases) and a low-level
feature-grained API that enables tyler's progressive pipeline.

## Motivation

Today tyler calls `geof` as a subprocess (`src/main.rs:696-728`). Each
tile spawns a process, passes parameters as command-line strings,
communicates input via an NDJSON file, and reads the output GLB from
disk. This has several costs:

- **Subprocess spawn overhead** per tile. Tyler spawns one process per
  leaf tile, of which there are typically thousands.
- **Serialization roundtrips.** Tyler writes NDJSON; geof parses it.
  Tyler ends up with a GLB on disk; never inspects it.
- **No feature-granular calls.** Calling geof per-feature would
  multiply the subprocess overhead. Progressive per-feature conversion
  is therefore impractical with the current design.
- **Deployment complexity.** Tyler depends on `geof` being installed and
  discoverable; errors manifest as subprocess failures rather than
  typed Rust errors.

Integrating geof as a Rust library eliminates all of these costs.

## Scope of the geof crate

The crate handles everything currently performed by the external `geof`
executable:

1. Parse CityJSONFeature input.
2. Apply attribute filtering / selection.
3. Apply LoD filtering.
4. Triangulate surfaces (for output formats that require it).
5. Compute normals.
6. Produce output in the target format:
   - 3D Tiles (GLB / glTF binary)
   - CityJSON (merged document)
   - OBJ
   - GeoPackage (GPKG)

The crate is usable both by tyler (in-process, with feature-grained API)
and as a standalone CLI (tile-grained, for one-shot conversions and
testing).

## API design

### Layer 1: high-level tile conversion

The simple API; the CLI wraps this directly. Also usable by tyler when
progressive behavior is not required.

```rust
pub enum Format {
    Gltf,
    Glb,
    CityJson,
    Obj,
    Gpkg,
}

pub struct Options {
    pub attribute_spec: Option<Vec<AttributeSelector>>,
    pub cotypes: Option<Vec<CityObjectType>>,
    pub metadata_class: Option<String>,
    pub geometric_error: Option<f64>,
    pub lod_filter: LodFilter,
    pub color_spec: ColorSpec,
    // ...
}

pub struct TileMetadata {
    pub bbox: [f64; 6],
    pub path_metadata: Vec<u8>,
    // ...
}

pub fn convert_tile(
    features: impl Iterator<Item = Result<CityJsonFeature>>,
    tile_metadata: &TileMetadata,
    output: impl Write,
    format: Format,
    options: &Options,
) -> Result<ConversionStats, Error>;
```

**Semantics:**
- Consumes features lazily from the iterator.
- Writes output to `output: impl Write` (file, buffer, pipe, etc.).
- Returns statistics (vertex count, feature count, conversion time) on
  success.
- Errors propagate through `Result`.

### Layer 2: feature-grained pipeline for progressive processing

Used by tyler's streaming loop to do per-feature work during ingestion
and per-tile work during assembly.

```rust
pub trait FeatureEncoder: Send + Sync {
    type Intermediate: Encode + Decode + Send;

    /// Encode a single feature into a format-specific intermediate
    /// representation. This is the expensive per-feature work
    /// (triangulation, normal computation, etc.) when required.
    fn encode_feature(
        &self,
        feature: &CityJsonFeature,
        options: &Options,
    ) -> Result<Self::Intermediate>;
}

pub trait TileAssembler: Send + Sync {
    type Intermediate;

    /// Assemble a tile from pre-encoded per-feature intermediates.
    /// This is mostly buffer layout work: concatenation, index
    /// rebasing, header packing.
    fn assemble_tile(
        &self,
        intermediates: impl Iterator<Item = Self::Intermediate>,
        tile_metadata: &TileMetadata,
        output: impl Write,
        options: &Options,
    ) -> Result<ConversionStats>;
}

// Concrete implementations:
pub struct GlbConverter;    // Intermediate = TriangulatedMesh
pub struct CityJsonConverter;// Intermediate = RawFeatureBytes (passthrough)
pub struct ObjConverter;    // Intermediate = ObjFragment
pub struct GpkgConverter;   // Intermediate = WkbGeometry + Attributes

impl FeatureEncoder for GlbConverter { type Intermediate = TriangulatedMesh; /* ... */ }
impl TileAssembler   for GlbConverter { type Intermediate = TriangulatedMesh; /* ... */ }
// and so on for each format
```

### Relationship between layers

Layer 1 is implemented on top of layer 2:

```rust
pub fn convert_tile<C>(
    features: impl Iterator<Item = Result<CityJsonFeature>>,
    tile_metadata: &TileMetadata,
    output: impl Write,
    format: Format,
    options: &Options,
) -> Result<ConversionStats, Error>
where
    C: FeatureEncoder + TileAssembler<Intermediate = <C as FeatureEncoder>::Intermediate>,
{
    let converter = converter_for(format);
    let intermediates: Vec<_> = features
        .map(|f| converter.encode_feature(&f?, options))
        .collect::<Result<_>>()?;
    converter.assemble_tile(intermediates.into_iter(), tile_metadata, output, options)
}
```

## The per-feature intermediates

Each format has its own intermediate type. The intermediate must be
serializable (the `Encode + Decode` bounds) so tyler can spill it to a
file and read it back during assembly.

### `TriangulatedMesh` (GLB / 3D Tiles)

```rust
pub struct TriangulatedMesh {
    pub feature_id: String,
    pub semantic_surfaces: Vec<SemanticSurface>,
    pub attributes: AttributeSet,
}

pub struct SemanticSurface {
    pub semantic: SemanticType,      // Roof, Wall, Ground, etc.
    pub positions: Vec<[f32; 3]>,    // triangulated vertices
    pub normals: Vec<[f32; 3]>,
    pub indices: Vec<u32>,           // triangle list
}
```

Wire format: binary, length-prefixed. Written with a fast encoder
(e.g., `bincode` or a bespoke zero-copy format).

### `RawFeatureBytes` (CityJSON)

```rust
pub struct RawFeatureBytes(pub Vec<u8>);
```

Passthrough. The per-feature "encoding" is just the raw CityJSON bytes.
Assembly concatenates features and reindexes vertices into a single
CityJSON document.

### `WkbGeometry + Attributes` (GPKG)

```rust
pub struct GpkgIntermediate {
    pub wkb: Vec<u8>,       // Well-Known Binary geometry
    pub attributes: Vec<AttributeValue>,
}
```

Triangulation is typically not required; tessellation may be applied if
the target GPKG consumer expects it. Assembly batches INSERTs into the
GeoPackage SQLite file.

### `ObjFragment` (OBJ)

```rust
pub struct ObjFragment {
    pub feature_id: String,
    pub vertices: Vec<[f32; 3]>,
    pub faces: Vec<Vec<u32>>,
}
```

The converter may triangulate or may keep n-gons depending on the OBJ
variant chosen. Assembly concatenates fragments with vertex index
rebasing.

## Streaming encoder ergonomics (optional third layer)

For tyler's streaming loop, the per-feature API is called inline in a
hot loop. To avoid allocating a `TriangulatedMesh` per feature in a
tight loop, the encoder can expose a "push-style" API that writes into
a caller-supplied output buffer:

```rust
pub trait StreamingEncoder {
    fn encode_feature_into(
        &self,
        feature: &CityJsonFeature,
        options: &Options,
        output: &mut impl Write,
    ) -> Result<()>;
}
```

Tyler then passes its cell-spill-file writer directly as `output`, and
no intermediate `TriangulatedMesh` is allocated. The binary format of
the spill file is defined by the encoder.

This is an optimization; the base API is sufficient for correctness.

## The geof CLI

The CLI is a thin wrapper around `convert_tile`:

```
geof --format glb input.city.jsonl output.glb [--tile-bbox x1,y1,z1,x2,y2,z2] [options...]
geof --format obj input.city.jsonl output.obj
geof --format cityjson input.city.jsonl output.city.json
geof --format gpkg input.city.jsonl output.gpkg
```

**Behavior:**
- Reads CityJSONFeatures from the input file (or stdin with `-`).
- Writes a single output file in the chosen format.
- Mirrors tyler's current geof-subprocess invocation options (colors,
  LoD filters, attribute spec, etc.) for backward compatibility with
  existing workflows.

**Purpose:**
- Testing and debugging the library without running tyler.
- Ad-hoc conversions (convert a single tile's worth of features).
- Backward compatibility — if tyler's library integration has a bug in
  some edge case, the same crate can be invoked as a subprocess as a
  fallback.

## Error model

Errors are typed:

```rust
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("failed to parse CityJSONFeature: {0}")]
    Parse(#[from] cjlib::Error),

    #[error("triangulation failed for feature {feature_id}: {reason}")]
    Triangulation { feature_id: String, reason: String },

    #[error("attribute {name} not found in feature {feature_id}")]
    MissingAttribute { name: String, feature_id: String },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    // ...
}
```

Tyler's streaming loop can decide per error whether to skip the feature,
abort the run, or log and continue.

## Dependencies

The geof crate depends on:

- `cjlib` — CityJSON parsing (already in tyler's workspace).
- A triangulation library (e.g., `earcut-rs` for 2D, or a 3D polygon
  triangulator).
- `gltf` or a hand-rolled GLB writer for 3D Tiles output.
- `rusqlite` + `gpkg-rs` (or similar) for GPKG output.
- `thiserror` for typed errors.

Tyler's `Cargo.toml` adds `geof = { path = "../geof" }` and drops the
`subprocess` crate dependency.

## Migration from geof subprocess

Once the library exists, tyler's migration is:

1. Replace the subprocess invocation in `src/main.rs:696-728` with
   `geof::convert_tile(...)`.
2. Remove the NDJSON write step in `write_inputs`
   (`src/main.rs:171-238`) — `convert_tile` takes an iterator of
   features directly, so the intermediate file is no longer needed for
   the in-process call.
3. Remove the `subprocess` crate dependency.
4. Deprecate the `--exe-geof` CLI flag (or keep it as a fallback for
   users who want to run the external CLI for reproducibility).

After migration, tyler's tile-export phase collapses from a
subprocess-orchestration to a library call inside the Rayon parallel
loop.

## Testability

Each layer is independently testable:

- **Feature encoders** — unit tests with representative CityJSONFeatures
  and golden-file comparisons on the encoded intermediate.
- **Tile assemblers** — unit tests with synthetic intermediates and
  golden-file comparisons on the output.
- **`convert_tile`** — integration tests with full CityJSONSeq files
  and golden-file comparisons on the final output.
- **The CLI** — end-to-end tests invoking the binary on fixture files.
