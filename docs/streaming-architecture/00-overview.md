# Streaming 3D Tile Generation — Design Overview

## Purpose

This document set captures the design for a streaming tile-generation
architecture in tyler. It enables the following pipeline:

```
roofer | morton-sort | tyler
```

Where `roofer` (the building reconstruction tool) emits reconstructed 3D
models as CityJSONSeq to stdout, and tyler generates 3D Tiles (or OBJ,
CityJSON, GPKG) progressively — without collecting the whole input into
memory and without waiting for EOF to begin conversion.

## Goals

1. **Memory efficiency.** Tyler's memory footprint must be bounded by the
   spatial grid plus a small per-feature index entry — **not** by the total
   size of the input geometry data.
2. **Elegance and clarity.** A clean architecture where the streaming
   concern is isolated from the tiling concern, and each tool in the
   pipeline has a single, well-defined responsibility.
3. **Performance.** Single-pass ingestion, zero re-serialization, parallel
   assembly. The expensive per-feature work (triangulation, format
   conversion) happens progressively during streaming, spread over the full
   runtime rather than concentrated at EOF.
4. **Multi-format output.** The same streaming architecture supports
   triangulated (3D Tiles / GLB) and non-triangulated (CityJSON, OBJ, GPKG)
   output formats.
5. **Optional progressive output.** When input arrives in (approximate)
   morton order, tiles are assembled and written to disk as the stream
   progresses — not at EOF.

## Non-goals

- Replacing the `SquareGrid` + `QuadTree` spatial structure. The existing
  design is preserved; streaming is layered on top.
- Handling unbounded streams (tyler still produces a finite tileset; EOF
  is required to write `tileset.json` and finalize sparse regions).
- Perfect morton ordering of the input. An approximate ordering is
  sufficient; the reorder buffer smooths local disorder.

## Document index

| # | Document | Topic |
|---|---|---|
| 00 | this file | Overview and goals |
| 01 | [current-architecture](./01-current-architecture.md) | Analysis of tyler's current pipeline and the streaming bottlenecks |
| 02 | [alternatives-considered](./02-alternatives-considered.md) | Streaming input alternatives evaluated |
| 03 | [progressive-tile-generation](./03-progressive-tile-generation.md) | Quadtree monotonicity, sealing, and layered progressive work |
| 04 | [morton-ordered-input](./04-morton-ordered-input.md) | The `morton-sort` tool and the low-water mark mechanism |
| 05 | [geof-library](./05-geof-library.md) | Native Rust geof library and CLI design |
| 06 | [final-architecture](./06-final-architecture.md) | The unified architecture and all its components |
| 07 | [implementation-phases](./07-implementation-phases.md) | Layered rollout plan |

## The pipeline at a glance

```
┌──────────┐    CityJSONSeq     ┌─────────────┐    morton-ordered    ┌─────────────────────┐
│  roofer  │ ─────stdout─────▶ │ morton-sort │ ───CityJSONSeq─────▶ │  tyler (--stdin)    │
└──────────┘                    └─────────────┘                      │                     │
                                                                     │  ┌─────────────┐   │
                                                                     │  │ geof crate  │◀──┤ linked library
                                                                     │  └─────────────┘   │
                                                                     │                     │
                                                                     │  streaming loop +   │
                                                                     │  progressive        │
                                                                     │  assembly workers   │
                                                                     └─────────────────────┘

Separately:  geof input.city.jsonl output.glb     ← standalone CLI
             (same library, one-shot tile mode)
```

## Summary of decisions

1. **Stream CityJSONSeq from stdin** with a pre-declared or header-provided
   extent. Features are spilled to per-cell files on disk as they arrive;
   only a lightweight `Feature` struct (centroid, bbox, reference) is kept
   in memory.
2. **Maintain a live quadtree** that tracks cell vertex counts. The
   monotonicity property of the quadtree merge decision (`sum > capacity`
   is permanent) lets us seal groups as soon as they exceed capacity.
3. **Standalone `morton-sort` tool** reorders the CityJSONSeq stream into
   approximate morton order using a bounded min-heap buffer. Tyler's
   low-water mark mechanism then uses the spatial ordering as an implicit
   completion signal for cells.
4. **Native Rust geof library** replaces the geof subprocess. Exposes a
   two-layer API: high-level `convert_tile()` for the CLI and simple
   cases, and low-level `FeatureEncoder` / `TileAssembler` traits for
   tyler's progressive pipeline.
5. **Format-specific per-feature encoders.** Each output format has its
   own per-feature intermediate (triangulated mesh for GLB, raw bytes for
   CityJSON, WKB for GPKG, OBJ fragment for OBJ). Streaming produces
   intermediates; assembly concatenates and packs them.
6. **Progressive assembly.** Leaf tiles are assembled as soon as their
   quadrant is finalized (via capacity seal or via morton frontier
   passing the quadrant). Output files grow continuously during
   processing; only `tileset.json` and sparse-region tiles are written at
   EOF.
7. **Optional `+RegionComplete` control signals** in CityJSONSeq provide
   an escape hatch for explicit completion from upstream tools that can
   afford it. Forward-compatible; tyler ignores unknown control lines.
