# Urchin Project — AI Coding Agent Instructions

## Architecture Overview

This is a **MATLAB-based B-Rep surface mesh generator** for nano-urchin particles (spiky spherical cores). The v4.x pipeline stitches **watertight `surfaceMesh` objects** by composing geometric patches (core sphere, conical spike bodies, spherical tip caps, optional toroidal fillets) with deterministic seam alignment.

**Key modules:**
- `matlab/src/urchin.m` — Main orchestrator: handles parameter parsing, spike orientation/length distribution, per-spike geometry derivation, and patch assembly
- `matlab/src/Urchin_Creator.m` — High-level batch driver with structured parameter groups and automated export workflows
- `matlab/src/meshDiagnostics.m` — Quality validation using MATLAB's Lidar Toolbox (`isWatertight`, `isEdgeManifold`, etc.)
- Helper functions (`solveConeSeams.m`, `refineSpikeOrientations.m`, `triangulate_icosphere.m`) — Geometry solvers and mesh primitives

## Critical Developer Workflows

### Running Tests
```matlab
cd matlab/tests
results = runtests;
table(results)
```
Tests cover deterministic defaults, volume exports, per-spike parameter derivation, and mesh diagnostics. Always run tests after modifying core geometry logic.

### Generating an Urchin
```matlab
addpath("matlab/src");
urchinStruct = urchin('rCore', 30, 'spikeLength', 15, 'spikeCount', 100, ...
    'spikeTip', 3, 'spikeConicality', 0.6);
diagnostics = meshDiagnostics(urchinStruct);
fprintf('Watertight: %d\n', diagnostics.IsWatertight);
```
The **return struct** contains `.SurfaceMesh` (the core output), `.Parameters` (echoes inputs), and `.Metrics` (derived quantities like `SpikeBaseMaxima`, `MinimumSpacing`, `TotalVolume`).

### Visualizing Refinement Process
Set `URCHIN_REFINE_VIS=1` in your environment before launching MATLAB to see the live Coulomb-based spike orientation relaxation solver (`refineSpikeOrientations.m`).

## Project-Specific Conventions

### Geometry Parameter Hierarchy
1. **Inputs** → `rCore`, `spikeLength`, `spikeCount`, `spikeTip` (diameter), `spikeConicality` (∈[−1,1])
2. **Intermediate** → `scale = 2*(rCore + spikeLength)`, `minSpacing = scale/resolution`
3. **Per-spike derived** → `spikeLengths` (with fluctuations), `spikeOrientations`, `rBaseMaxs`, `rTipMax`, `alphaCone`
4. **Seam tangency** → `solveConeSeams(rBases, rTips, zBases, zTipCenters)` returns `zTipSeams`, `rTipSeams`, `alphaCones` by enforcing cone-sphere tangency

**Document updates:** All geometry parameters and their dependencies are catalogued in `documents/urchin_geometry_parameters.md` (categorized into Core/General and Spike Configuration).

### Determinism and Stochasticity
- **Default `flucFactor=0`** → deterministic spike lengths (v4.0.1 fix)
- `flucMethod='uniform'` uses Sobol sequences; `'random'` uses seeded RNG
- `distMethod='uniform'` applies Fibonacci sphere packing; `'random'` uses seeded orientations
- `refinedOrientation=true` enables electrostatic-like relaxation to maximize angular separation

### Breaking Changes (v3.x → v4.x)
- **Old signature:** `[mesh, mask, threshold, eqRadius] = urchin(...)`
- **New signature:** `urchinStruct = urchin(...)` returns a struct with `.SurfaceMesh`, `.Vertices`, `.Faces`, `.Parameters`, `.Metrics`
- Volume outputs (`.VolumeMask`, `.VolumeOctree`) are **opt-in** via `genVolume=true`

### Mesh Stitching Patterns
- `append_vertices(V, newVerts)` deduplicates globally
- `bridge_rings_idx(ring1Idx, ring2Idx)` creates quad strips between circular seam rings
- Core trimming: compute per-spike footprint, remove core vertices inside circular regions, regenerate boundary loops
- **Watertightness:** seam rings must share identical precomputed vertices to avoid gaps

### Helper Function Roles
- `plane_vectors(orientation)` → orthonormal basis for circular spike cross-sections
- `ring_segment_count(radius, minSpacing)` → polygon count for curvature-aware sampling
- `seam_ring_points(center, orientation, radius, nSegments)` → deterministic ring vertices
- `computeMeshVolume(V, F)` → signed volume via divergence theorem

## Integration Points

### External Dependencies
- **Statistics and Machine Learning Toolbox** → `sobolset` for quasi-random sequences
- **Lidar Toolbox** → `surfaceMesh`, `surfaceMeshShow`, `isWatertight`, mesh diagnostics
- **Image Processing Toolbox** (optional) → `volshow` for dense volume visualization
- **Parallel Computing Toolbox** (optional) → GPU experiments (not actively used)

### Export Formats
- `writeSurfaceMesh(mesh, 'file.stl', Encoding='binary')` → STL for COMSOL/FEM
- `Urchin_Creator.m` automates JSON config exports, MAT file archives, and volume masks

## Common Pitfalls

### Short Spike Handling
When `spikeLength` < `rTip`, the code **clamps tip radii per spike** to maintain tangency. Do not assume `rTip = spikeTip/2` always holds—check `urchinStruct.Metrics.SpikeTipRadius`.

### Resolution Tuning
`resolution` is dimensionless; effective spacing = `2*(rCore+spikeLength)/resolution`. High `resolution` → dense mesh but slower generation. Test with `resolution=50` first.

### Zero Spike Counts
`spikeCount=0` is valid and produces only the core sphere. Tests explicitly cover this edge case.

### Collapsed Tips
When `spikeTip` falls below `minSpacing`, tips may collapse to points. Check `urchinStruct.Metrics.SpikeTipCollapsed` flag.

## Citation and Licensing
Developed for *"Influence of Spike Geometry on Electromagnetic Field Enhancement..."* (ACS Applied Nano Materials, DOI: 10.1021/acsanm.5c03297). MIT licensed. See `CITATION.cff` for BibTeX metadata.
