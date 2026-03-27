# Use Segment-to-Grid Traversal for Feature Cell Ownership

## Status

Accepted

## Related Commits

- none yet

## Context

Tyler assigns some feature classes, especially `Building` and
`BuildingPart`, to exactly one grid cell. The per-cell "vertex count" used
for that assignment is not a true geometry metric. It is an ownership score
used to answer a simpler question:

"Which cell should own this feature?"

The legacy implementation built that score from two inputs:

- feature vertices located in cells
- a fallback `+1` for every cell intersecting the feature bbox

The bbox fallback existed for a real reason. A large polygon can overlap
multiple cells without having a stored vertex inside each overlapped cell. If
ownership uses only vertex locations, those crossed cells disappear from the
candidate set.

However, bbox expansion is a coarse proxy for actual geometry coverage:

- it includes cells touched only by the bbox, not the geometry
- it can over-credit elongated or rotated features
- it makes the score depend on a rough envelope instead of topology

After the move to `cjlib` and `cityjson-rs`, Tyler now works with real-world
coordinates directly. That makes it practical to traverse geometry segments in
world space and map them onto grid cells without preserving the legacy
quantized-coordinate logic.

## Decision

Tyler will replace the bbox-based ownership fallback with segment-to-grid
traversal.

For feature classes that require unique cell assignment:

1. Tyler will continue to use a per-cell ownership score rather than a strict
   geometric measure such as area.
2. Candidate cells will be derived from actual geometry traversal:
   - cells containing selected geometry vertices
   - cells crossed by selected geometry boundary segments
3. The old bbox-based `+1` fallback will be removed from the ownership score.
4. Tie-breaking should remain deterministic. If two cells receive the same
   ownership score, Tyler should prefer a stable geometric fallback such as the
   centroid cell or a bbox-center based rule.

The first implementation target is segment traversal in the XY plane over the
selected geometry boundaries. This keeps the heuristic aligned with tiling
ownership while avoiding full polygon overlay work.

## Consequences

Good:

- overlapped cells without interior vertices are still discovered
- candidate cells come from actual geometry coverage instead of bbox inflation
- the ownership score becomes less sensitive to feature storage artifacts
- the approach matches the `cjlib` / `cityjson-rs` real-world coordinate model

Trade-offs:

- Tyler needs a robust grid traversal algorithm for 2D segments
- the score is still a heuristic, not an exact measure of area overlap
- boundary traversal costs more than bbox intersection alone
- curved or highly detailed boundaries may touch many cells and need careful
  implementation to stay fast

## Rejected Alternatives

- Keep the bbox fallback.
  This preserves legacy behavior but keeps the ownership score tied to a coarse
  proxy that can add cells the geometry never touches.

- Use only vertex locations.
  This misses cells crossed by long segments or large polygons with sparse
  vertices.

- Move directly to polygon-cell overlap area.
  This is more semantically precise, but it is a much larger geometry problem
  and not the right first step for restoring robust ownership assignment.

- Reintroduce legacy quantized-coordinate counting.
  Tyler no longer needs to preserve quantized-coordinate handling now that
  `cjlib` and `cityjson-rs` provide real-world coordinates directly.

## Notes

This ADR is about cell ownership scoring, not about reported geometry
statistics. The implementation should eventually rename the current
`nr_vertices`-style counters to reflect that they are ownership scores rather
than literal vertex totals.
