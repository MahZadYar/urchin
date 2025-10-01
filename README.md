# urchin

`urchin` is a MATLAB toolkit for generating deterministic, curvature-aware B-Rep
surface meshes of nano-urchin particles. Version 4.0 introduces a brand-new
watertight meshing pipeline that replaces the legacy voxel workflow while
retaining optional volume exports for users that still need them. The codebase
originated in support of the peer-reviewed study *"Influence of Spike Geometry
on Electromagnetic Field Enhancement and the Linear and Nonlinear Plasmonic
Properties of Gold Nanourchins"* (ACS Applied Nano Materials, DOI:
[10.1021/acsanm.5c03297](https://doi.org/10.1021/acsanm.5c03297)).

## 🚀 What’s new in v4.0

- **B-Rep first** – spikes, caps, and core patches are stitched in real time to
    create a single watertight `surfaceMesh` object ready for COMSOL and similar
    solvers.
- **Deterministic seam trimming** – per-spike trimming keeps seams aligned and
    prevents duplicate faces or holes even for dense orientation sets.
- **Built-in diagnostics** – every mesh can be validated in milliseconds with
    `isWatertight`, `isEdgeManifold`, `isOrientable`, `isSelfIntersecting`, and
    `isVertexManifold` checks.

See the full release notes in [`CHANGELOG.md`](CHANGELOG.md).

## Project layout

- `matlab/src/` – core source code (`urchin.m` and helpers).
- `matlab/examples/` – runnable scripts showing typical workflows.
- `matlab/tests/` – unit tests for trimming, stitching, and diagnostics.
- `CITATION.cff` – citation metadata for the project.

## Requirements

- MATLAB R2023b or newer (tested).
- Statistics and Machine Learning Toolbox (`sobolset`, `net`).
- Lidar Toolbox (`surfaceMesh`, `surfaceMeshShow`, `isWatertight`, …).
- Optional: Image Processing Toolbox for voxel utilities such as `volshow`.
- Optional: Parallel Computing Toolbox for GPU experiments.

## Quick start

1. Add the source directory to your MATLAB path:

     ```matlab
     addpath("path/to/urchin/matlab/src");
     ```

2. Generate an urchin and collect diagnostics:

     ```matlab
     [mesh, diagnostics] = urchin('cr', 30, 'sl', 15, 'ns', 100, 'st', 2);

     fprintf("Watertight: %d\n", diagnostics.IsWatertight);
     surfaceMeshShow(mesh, WireFrame=true);
     ```

3. Enable optional volume exports only when needed:

     ```matlab
     params = { 'genVolume', true, 'volAdaptive', true, 'volRes', 192 };
     mesh = urchin('cr', 20, 'sl', 12, 'ns', 64, params{:});
     writeSurfaceMesh(mesh, "urchin.stl");
     ```

More end-to-end demonstrations are available in
`matlab/examples/run_urchin_creation.m`.

## Testing

MATLAB-based unit tests live under `matlab/tests/`. From MATLAB, run:

```matlab
cd path/to/urchin/matlab/tests
results = runtests;
table(results)
```

## Citation

If this work supports your research, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE) for
details.
