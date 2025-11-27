# Changelog

<!-- markdownlint-disable MD024 -->

All notable changes to this project will be documented in this file. The
format roughly follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.3.0] - 2025-11-27

### Added

- **Spike Orientation Refinement**: Implemented a Coulomb-like relaxation solver
  that redistributes spike orientations to maximise angular separation before
  geometry generation.
- Added `refinedOrientation` and `refinedOrientationThreshold` parameters to
  control the relaxation process.
- Introduced a momentum-based update strategy to ensure stable convergence without
  oscillation, even for high spike counts.
- Added an `URCHIN_REFINE_VIS` environment toggle to visualise the relaxation
  process in real-time.
- **Cone Seam Solver**: Redefined the calculation logic for solving cone seams in
  `solveConeSeams.m` to robustly handle the tangency condition between conical
  bodies and spherical tips.

### Changed

- Updated `urchin.m` to invoke the refinement pass when requested.

## [4.2.0] - 2025-10-10

<!-- markdownlint-disable-next-line MD024 -->
### Added

- Derived per-spike maximum base radii from orientation dot products and exposed
  them via `Metrics.SpikeBaseMaxima`.
- Extended automated tests to cover the new per-spike maxima output.

<!-- markdownlint-disable-next-line MD024 -->
### Changed

- Used the per-spike maxima bounds during cone synthesis so short spikes retain
  stable bases without exceeding their angular neighbourhoods.
- Refined documentation and examples to describe the per-spike maxima metric.

<!-- markdownlint-disable-next-line MD024 -->
### Fixed

- Removed redundant base radius tightening when per-spike maxima already bound
  the footprint, yielding consistent results across dense orientation sets.

## [4.1.0] - 2025-10-04

### Added

- Introduced per-spike spherical tip clamping so short spikes keep their tip
  seams tangent without collapsing into the core.
- Added a spacing-driven cone sampling routine that adapts ring placement to
  each spike's axial length while allowing adjacent rings to carry different
  polygon counts.

### Changed

- Tightened the short-spike base fallback so cones maintain a finite footprint
  even when the tangential base limit degenerates.
- Refined documentation and examples to highlight the new geometry safeguards.

### Fixes

- Prevented dense ring stacking on short spikes by decoupling axial ring count
  from per-ring segment interpolation.

## [4.0.1] - 2025-10-01

### Fixed

- Corrected the default spike fluctuation factor (`sf`) to 0 so meshes are
  deterministic unless jitter is explicitly requested.
- Wrapped the surface mesh output in a struct so auxiliary data (volume mask,
  octree metadata) can be attached without breaking MATLAB's `surfaceMesh`
  class.
- Updated MATLAB regression tests to target the B-Rep pipeline, covering
  volume exports and the fluctuation toggle.

## [4.0.0] - 2025-10-01

### Highlights

- Replaced the voxel-oriented pipeline with a watertight B-Rep surface mesh
  generator that stitches cone, cap, and core patches on-the-fly.
- Added deterministic per-spike trimming and seam bridging, eliminating the
  duplicate-face and leak issues present in the 3.x series.
- Introduced built-in mesh diagnostics (watertight, manifold, orientable,
  self-intersection checks) powered by MATLAB Surface Mesh functionality.

### Breaking changes

- The primary output of `urchin` is now a `surfaceMesh` object. The function
  returns `[mesh, diagnostics]` instead of the previous `[mesh, mask,
  threshold, eqRadius]` signature. Callers that relied on the old outputs must
  update their code.
- Voxel masks are no longer produced automatically. Dense and adaptive volume
  exports are still available, but only when the corresponding name-value
  parameters are enabled.

### Improvements

- Seam vertices are now incorporated into subsequent trims, ensuring that
  adjacent spikes share rings without manual bookkeeping.
- Added robust fallback handling for minimal spike footprints, preventing trim
  failures on high-density orientation sets.
- Simplified seam bridging logic to avoid duplicate ring connections and reduce
  numerical instability near the core surface.

## [3.5.0] - 2025-08-07

- Refined voxel mask generation workflow and documentation tweaks ahead of the
  B-Rep rewrite.

## [3.0.0] - 2025-05-27

- Initial public release after the project rename to **urchin**, featuring the
  voxel-based meshing pipeline and MATLAB GUI assets.
