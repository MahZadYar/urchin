%% urchin v4.0.1 example
%
% This script demonstrates the B-Rep workflow introduced in version 4.0 and the
% deterministic defaults restored in patch 4.0.1.

% Add the source directory to the MATLAB path when running the script directly.
% addpath(fullfile(fileparts(mfilename("fullpath")), "..", "src"));

%% Parameter block
cr = 30;          % Core radius (nm)
sl = 15;          % Spike length measured from the core surface (nm)
ns = 128;         % Number of spikes packed on the core
st = 3;           % Seam diameter at the distal cut plane (nm)
sc = 0.6;         % Conicality (0=cylindrical, 1=widest non-overlapping base)
sf = 0.35;        % Spike length fluctuation factor (set to zero to disable)

refine = 1.2;     % Unified refinement multiplier for all patch samplings
distMethod = "uniform"; % "uniform" or "random" spike orientations
flucMethod = "uniform"; % "uniform", "random", or "gaussian" fluctuations

fprintf("Generating urchin with %d spikes...\n", ns);

%% Core mesh generation with diagnostics
[mesh, diagnostics] = urchin( ...
    "cr", cr, "sl", sl, "ns", ns, "st", st, ...
    "sc", sc, "sf", sf, ...
    "meshRefine", refine, ...
    "distMethod", distMethod, ...
    "flucMethod", flucMethod ...
);

fprintf("Watertight: %d\n", diagnostics.IsWatertight);
fprintf("Edge manifold: %d\n", diagnostics.IsEdgeManifold);
fprintf("Self-intersections: %d\n", diagnostics.IsSelfIntersecting);

%% Visualize the B-Rep surface
viewer = viewer3d;
surfaceMeshShow(mesh.SurfaceMesh, Parent=viewer, WireFrame=true);
surfaceMeshShow(mesh.SurfaceMesh, Parent=viewer, Alpha=0.6);

%% Optional: adaptive sparse volume (turned off by default)
% volumeMesh = urchin( ...
%     "cr", cr, "sl", sl, "ns", ns, "st", st, ...
%     "sc", sc, "sf", sf, ...
%     "volAdaptive", true, "genVolume", true, ...
%     "volRes", 160 ...
% );
% leaves = volumeMesh.VolumeOctree.Leaves; %#ok<NASGU>

%% Optional: export STL
% writeSurfaceMesh(mesh.SurfaceMesh, "urchin_v4_example.stl");

fprintf("Done.\n");
