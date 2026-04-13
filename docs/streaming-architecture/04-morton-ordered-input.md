# 04 — Morton-Ordered Input and the `morton-sort` Tool

This document specifies the `morton-sort` tool and the low-water mark
mechanism that uses spatial ordering as an implicit completion signal.

## The problem

In [03-progressive-tile-generation](./03-progressive-tile-generation.md)
we established that a sealed leaf tile is assemblable only when its
feature set is complete. Without an ordering guarantee, we cannot know
when more features might arrive in a cell.

## The insight

If features arrive in (approximate) morton order of their centroid,
tyler can maintain a **low-water mark** — the maximum morton code seen
so far. The invariant:

> **Low-water mark invariant:** all cells with morton code strictly less
> than the low-water mark have received all their features.

This is sufficient to finalize the tile structure for completed
quadrants.

## How it plays with the quadtree

Tyler's quadtree is built from morton-ordered cells. Groups of 4
consecutive morton codes form quadrants at each level:

```
Morton order:  0  1  2  3 │ 4  5  6  7 │ 8  9  10 11 │ 12 13 14 15
               ╰────┬────╯ ╰────┬────╯   ╰────┬─────╯  ╰────┬─────╯
Level N-1:       Q0           Q1              Q2             Q3
                 ╰──────┬──────╯              ╰──────┬───────╯
Level N-2:           QQ0                          QQ1
```

As the low-water mark advances:

1. Low-water mark passes cell 3 → cells {0,1,2,3} are complete → Q0's
   merge decision is final.
2. Decision: if `sum(Q0) <= capacity` → merged leaf tile → **assemble
   immediately**. If `sum(Q0) > capacity` → each cell is its own leaf
   tile → **assemble each immediately**.
3. Low-water mark passes cell 7 → Q1 is also final → assemble its
   leaves.
4. Q0 and Q1 both final → QQ0's merge decision is also final → assemble
   if merged.

The quadtree builds itself incrementally, bubbling up. Leaf tiles are
assembled as soon as their quadrant is decided. No regeneration, no
speculation, no dirty tracking.

## Combined with capacity sealing

The two progressive mechanisms cooperate:

| Condition | Boundary finalized | Feature set finalized |
|---|---|---|
| Cell exceeds capacity | Immediately | When low-water mark passes cell |
| Quadrant complete (via LWM) | When LWM passes last sibling | Simultaneously |

So dense cells (exceeding capacity individually) are assemblable the
moment the LWM passes them. Sparse cells are assemblable when their full
quadrant has been traversed by the LWM.

## The `morton-sort` tool

A standalone, pipeable, single-purpose tool that sits in the pipe
between the producer (roofer) and the consumer (tyler):

```
roofer | morton-sort | tyler --stdin
```

### Responsibility

Reorder a CityJSONSeq stream into approximate morton order using a
bounded-size reorder buffer (min-heap).

### Invariant

When `morton-sort` emits a feature with morton code M, no future output
will have a code strictly less than M. The buffer size determines the
maximum local disorder it can absorb; beyond that, the tool must either
(a) reject out-of-order features, (b) emit them anyway with a warning
(downstream treats them as stragglers), or (c) grow the buffer (with a
memory cap).

### Interface

```
morton-sort [--buffer-size 10000] [--extent x1,y1,x2,y2] [--order yx|xy] [--on-late warn|error|emit]
```

- `--buffer-size` (default: 10000): maximum number of features buffered.
- `--extent`: spatial extent used to quantize coordinates into integer
  grid indices before bit-interleaving. If omitted, read from the
  CityJSONSeq header's `metadata.geographicalExtent` field.
- `--order`: `yx` (Z-curve, matches tyler's current morton encoding at
  `src/spatial_structs.rs:49-51`) or `xy` (N-curve).
- `--on-late`: behavior when a feature arrives with a morton code less
  than the last emitted code.

### Behavior

1. **Header passthrough.** The first line of the CityJSONSeq stream is
   the metadata document. Pass it through unchanged immediately; begin
   buffering from line 1 onward.
2. **Per-feature processing.**
   - Parse just enough of the feature to extract its centroid (e.g.,
     via `cjlib` or a minimal JSON scan).
   - Compute morton code from the centroid using the extent to quantize.
   - Insert `(morton_code, raw_bytes)` into a `BinaryHeap<Reverse<(u128, Vec<u8>)>>`.
   - If the heap has reached `buffer-size`, pop the minimum and emit it.
3. **Control line passthrough.** Lines starting with `{"type": "+...` are
   control signals (e.g., `+RegionComplete`). These must not be
   reordered relative to features. On encountering a control line:
   drain the heap (emit all buffered features in morton order), emit the
   control line, resume buffering.
4. **EOF.** Drain the heap in morton order and close.

### Pseudo-code

```rust
struct ReorderBuffer {
    heap: BinaryHeap<Reverse<(u128, Vec<u8>)>>,
    capacity: usize,
    extent: Extent,
    last_emitted: u128,
    on_late: LatePolicy,
}

impl ReorderBuffer {
    fn push(&mut self, bytes: Vec<u8>, stdout: &mut impl Write) -> Result<()> {
        let centroid = extract_centroid(&bytes)?;
        let morton = morton_code(centroid, &self.extent);
        if morton < self.last_emitted {
            match self.on_late {
                LatePolicy::Error => return Err("feature arrived out of order beyond buffer"),
                LatePolicy::Warn  => eprintln!("late feature: morton {} < {}", morton, self.last_emitted),
                LatePolicy::Emit  => {},
            }
        }
        self.heap.push(Reverse((morton, bytes)));
        if self.heap.len() > self.capacity {
            self.emit_min(stdout)?;
        }
        Ok(())
    }

    fn emit_min(&mut self, stdout: &mut impl Write) -> Result<()> {
        if let Some(Reverse((morton, bytes))) = self.heap.pop() {
            self.last_emitted = morton;
            stdout.write_all(&bytes)?;
            if !bytes.ends_with(b"\n") { stdout.write_all(b"\n")?; }
        }
        Ok(())
    }

    fn drain(&mut self, stdout: &mut impl Write) -> Result<()> {
        while !self.heap.is_empty() {
            self.emit_min(stdout)?;
        }
        Ok(())
    }
}
```

### Memory cost

`buffer_size × avg_feature_bytes`. For 10K features averaging ~2KB each,
~20MB. Negligible.

### Design rationale

**Why a separate tool and not built into tyler?**

- Single responsibility — each tool does one thing well.
- Independently testable. `morton-sort` can be tested in isolation with
  synthetic input.
- Composable. Users can drop in a different sorter (e.g., Hilbert
  curve) without modifying tyler.
- Pipeable. Standard Unix tooling idioms.
- Debuggable. You can observe the sorted stream by piping `morton-sort`
  to a file.
- Future-proof. If other tools want ordered CityJSONSeq input, they can
  use the same sorter.

**Why morton order specifically (not Hilbert, not row-major)?**

Tyler's existing quadtree is built in morton order
(`src/spatial_structs.rs:52-55`). Matching the sorter's output order to
the consumer's consumption order is what gives us the low-water mark
guarantee on the quadtree.

**Why a reorder buffer (not a strict sort)?**

A strict sort of an unbounded stream is impossible — it would require
buffering everything. A reorder buffer trades a small amount of memory
for a bounded-disorder guarantee that is sufficient in practice.

Roofer's output has strong spatial coherence (it processes neighborhoods
sequentially), so the natural local disorder is small. A buffer of
10K–50K features is enough to produce effectively-sorted output for
real-world inputs.

### When to skip morton-sort

If roofer (or another producer) already emits features in morton order
natively — e.g., because it processes its input in morton-indexed order —
the `morton-sort` tool can be omitted:

```
roofer --morton-output | tyler --stdin
```

Tyler's low-water mark mechanism works identically whether the sorter
exists or not. The sorter is an adapter for producers that don't
natively order their output.

## Tyler's side: the low-water mark

Tyler's streaming loop maintains a single `u128` counter: the maximum
morton code of any feature seen so far (approximately — see caveat
below). For each finalized cell (morton code < LWM), tyler checks if the
cell's quadrant is complete; if so, it triggers assembly.

### Caveat on monotonicity

The LWM should be advanced only as a monotonically increasing counter.
If `morton-sort`'s buffer absorbs all local disorder, the stream tyler
receives is perfectly ordered and the LWM is simply the last-seen code.
If a late feature slips through (with `--on-late emit`), tyler must
treat it as belonging to an already-finalized cell; depending on
severity, this could cause a tile regeneration or a warning.

With a sufficiently large buffer, late features are rare or impossible,
and the LWM logic is trivial.

## Optional: `+RegionComplete` control lines

For producers that can afford it, an explicit completion signal can be
added to the CityJSONSeq stream:

```
{"type": "+RegionComplete", "bbox": [x1, y1, x2, y2]}
```

Tyler treats this as an explicit LWM advancement: all cells fully
contained in `bbox` are marked complete, regardless of morton ordering.
The `+` prefix makes it a tyler-specific extension; readers that don't
understand it can skip.

This is strictly additive — it does not replace the LWM mechanism.
Tyler uses whichever mechanism finalizes a cell first. Producers that
cannot produce ordered output may prefer completion signals over a sort
step.
