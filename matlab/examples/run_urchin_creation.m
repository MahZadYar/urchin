%% urchin v4.3.0 example
%
% This script demonstrates the B-Rep workflow introduced in version 4.0, the
% deterministic defaults restored in patch 4.0.1, the short-spike safeguards added
% in release 4.1.0 (per-spike tip clamping, spacing-driven cone sampling), the
% per-spike base maxima introduced in release 4.2.0, and the new spike orientation
% refinement solver introduced in release 4.3.0.

% Add the source directory to the MATLAB path when running the script directly.
% addpath(fullfile(fileparts(mfilename("fullpath")), "..", "src"));

%% Parameter block
coreRadius = 30;          % Core radius (nm)
spikeLength = 15;         % Spike length measured from the core surface (nm)
spikeCount = 128;         % Number of spikes packed on the core
spikeTip = 3;             % Diameter of the spherical tip cap (nm)
spikeConicality = 0.6;    % Conicality (0=cylindrical, 1=widest non-overlapping base)
flucFactor = 0.35;        % Spike length fluctuation factor (set to zero to disable)

refine = 1.2;     % Unified refinement multiplier for all patch samplings
distMethod = "uniform"; % "uniform" or "random" spike orientations
flucMethod = "uniform"; % "uniform", "random", or "gaussian" fluctuations

fprintf("Generating urchin with %d spikes...\n", spikeCount);

%% Core mesh generation with diagnostics
urchinStruct = urchin( ...
    "coreRadius", coreRadius, "spikeLength", spikeLength, "spikeCount", spikeCount, "spikeTip", spikeTip, ...
    "spikeConicality", spikeConicality, "flucFactor", flucFactor, ...
    "resolution", 100 * refine, ...
    "distMethod", distMethod, ...
    "flucMethod", flucMethod ...
);

diagnostics = meshDiagnostics(urchinStruct);

fprintf("Watertight: %d\n", diagnostics.IsWatertight);
fprintf("Edge manifold: %d\n", diagnostics.IsEdgeManifold);
fprintf("Self-intersections: %d\n", diagnostics.IsSelfIntersecting);

%% Visualize the B-Rep surface
viewer = viewer3d;
surfaceMeshShow(urchinStruct.SurfaceMesh, Parent=viewer, WireFrame=true);
surfaceMeshShow(urchinStruct.SurfaceMesh, Parent=viewer, Alpha=0.6);

%% Optional: adaptive sparse volume (turned off by default)
% volumeMesh = urchin( ...
%     "coreRadius", coreRadius, "spikeLength", spikeLength, "spikeCount", spikeCount, "spikeTip", spikeTip, ...
%     "spikeConicality", spikeConicality, "flucFactor", flucFactor, ...
%     "volAdaptive", true, "genVolume", true, ...
%     "volRes", 160 ...
% );
% leaves = volumeMesh.VolumeOctree.Leaves; %#ok<NASGU>

%% Optional: export STL
% writeSurfaceMesh(urchinStruct.SurfaceMesh, "urchin_v4_example.stl");

fprintf("Done.\n");
