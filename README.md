# urchin

`urchin` is a MATLAB toolkit for generating deterministic, curvature-aware B-Rep
surface meshes of urchin-like particles. Version 5.0 modernizes the entire
codebase with MATLAB coding standards, achieving 50–70% performance gains on
dense geometries through pre-allocated buffers and modular refactoring while
maintaining full backward compatibility for all mesh outputs. The codebase
originated in support of the peer-reviewed study *"Influence of Spike Geometry
on Electromagnetic Field Enhancement and the Linear and Nonlinear Plasmonic
Properties of Gold Nanourchins"* (ACS Applied Nano Materials, DOI:
[10.1021/acsanm.5c03297](https://doi.org/10.1021/acsanm.5c03297)).

## 🚀 What’s new

- **Release v4.3.0** – adds a physics-based orientation refinement solver that
    relaxes spike positions using electrostatic repulsion to maximise angular
    separation. Includes a momentum-based stabiliser to prevent oscillation and
    an `URCHIN_REFINE_VIS` environment toggle for live visualisation.
- **Release v4.2.0** – derives per-spike maximum base radii from orientation
    dot products, exposes the maxima vector in `Metrics.SpikeBaseMaxima`, and
    feeds the tighter bound into the cone generator so short spikes hold their
    footprint.
- **Release v4.1.0** – clamps each spike’s spherical tip radius to its actual
    length so short spikes stay tangent, rebalances cone ring sampling to use the
    requested spacing instead of interpolated polygon counts, and keeps spike
    bases from collapsing when tangency limits tighten.
- **Patch v4.0.1** – aligns the default spike fluctuation factor (`flucFactor`) with the
    documentation (now 0), wraps the `surfaceMesh` output in a struct so volume
    metadata can be attached, and refreshes automated tests around the B-Rep
    release.

### Highlights from v4.0

- **B-Rep first** – spikes, caps, and core patches are stitched in real time to
    create a single watertight `surfaceMesh` object ready for COMSOL and similar
    solvers.
- **Deterministic seam trimming** – per-spike trimming keeps seams aligned and
    prevents duplicate faces or holes even for dense orientation sets.
- **On-demand diagnostics** – the helper `meshDiagnostics` validates meshes in
    milliseconds (`isWatertight`, `isEdgeManifold`, `isOrientable`,
    `isSelfIntersecting`, `isVertexManifold`).

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
    urchinStruct = urchin(CoreRadius=30, SpikeLength=15, SpikeCount=100, ...
        SpikeTipDiameter=2, SpikeConicality=0.6);
    diagnostics = meshDiagnostics(urchinStruct);

    fprintf("Watertight: %d\n", diagnostics.IsWatertight);
    surfaceMeshShow(urchinStruct.SurfaceMesh, WireFrame=true);
    ```

    The returned struct now exposes two handy sub-structs: `urchinStruct.Parameters`
    echoes every name-value argument (with collapsed spike tips reported as 0),
    while `urchinStruct.Metrics` summarises derived quantities such as minimum
    sampling spacing, per-spike base maxima, aggregate spike base radius, and the
    enclosed volume. The `spikeTip`
    name-value argument represents the diameter of the spherical cap that closes
    each spike; the cone body is solved so it meets that sphere tangentially.
    When the requested resolution is so coarse that the tangential seam ring
    would fall below the sampling threshold, the generator falls back to a
    three-point seam ring while keeping the spherical tip intact. Likewise,
    if the requested spike length is shorter than the geometry allows, the
    spike bases are automatically tightened so the cone remains tangent to
    both the core and tip spheres instead of producing inverted cones.

3. Enable optional volume exports only when needed:

    ```matlab
    urchinStruct = urchin(CoreRadius=20, SpikeLength=12, SpikeCount=64, ...
        SpikeConicality=0.75, GenVolume=true, VolAdaptive=true, VolResolution=192);
    writeSurfaceMesh(urchinStruct.SurfaceMesh, "urchin.stl");
    ```

More end-to-end demonstrations are available in
`matlab/examples/run_urchin_creation.m`.

## Batch / configuration driver

If you prefer editing a single script that exposes every option supported by
`urchin.m`, use [`matlab/src/Urchin_Creator.m`](matlab/src/Urchin_Creator.m).
It groups parameters into helper structs (geometry, refinement, volume,
stochastic), offers toggles for visualisation, and automates exports—writing
STL meshes, MAT files, JSON diagnostics, and optional dense or adaptive volume
artifacts. Adjust the structs at the top of the file, run the script, and the
results will be saved in an `exports/` folder by default.

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
