classdef test_urchin < matlab.unittest.TestCase
    % Regression tests for the B-Rep based urchin generator

    properties (Constant)
        % Use a lightweight configuration to keep runtimes short in CI
        SmallSpikeCount = 24;
    end

    methods (TestClassSetup)
        function addProjectPaths(~)
            thisDir = fileparts(mfilename('fullpath'));
            srcDir = fullfile(thisDir, '..', 'src');
            addpath(srcDir);
        end
    end

    methods (Test)
        function testDefaultRunProducesSurfaceMesh(testCase)
            % Default call should yield a surfaceMesh; diagnostics are computed on demand
            urchinStruct = urchin(SpikeCount=testCase.SmallSpikeCount);
            diagnostics = meshDiagnostics(urchinStruct);

            testCase.verifyTrue(isstruct(urchinStruct) && isfield(urchinStruct, "SurfaceMesh"), ...
                "urchin should return a struct exposing the underlying surface mesh.");
            testCase.verifyTrue(isa(urchinStruct.SurfaceMesh, "surfaceMesh"), ...
                "SurfaceMesh field must be a surfaceMesh instance.");
            testCase.verifyGreaterThan(size(urchinStruct.Vertices, 1), 0, ...
                "Default run should produce vertices.");
            testCase.verifyGreaterThan(size(urchinStruct.Faces, 1), 0, ...
                "Default run should produce faces.");

            testCase.verifyTrue(isfield(urchinStruct, "Parameters"), ...
                "Returned struct should list the input parameters.");
            params = urchinStruct.Parameters;
            testCase.verifyTrue(isfield(params, "SpikeConicality"), ...
                "Parameters must include SpikeConicality.");
            testCase.verifyGreaterThanOrEqual(params.SpikeTipDiameter, 0, ...
                "Spike tip diameter should be non-negative and zero when collapsed.");

            testCase.verifyTrue(isfield(urchinStruct, "Metrics"), ...
                "Returned struct should expose derived metrics.");
            metrics = urchinStruct.Metrics;
            testCase.verifyGreaterThan(metrics.MinimumSpacing, 0, ...
                "Minimum spacing metric should be positive.");
            testCase.verifyTrue(isfield(metrics, "SpikeBaseMaxima"), ...
                "Metrics should include per-spike maximum base radii.");
            testCase.verifySize(metrics.SpikeBaseMaxima, [testCase.SmallSpikeCount, 1], ...
                "Per-spike maxima vector should match the spike count.");
            testCase.verifyGreaterThanOrEqual(min(metrics.SpikeBaseMaxima), 0, ...
                "Per-spike base maxima should be non-negative.");
            testCase.verifyGreaterThanOrEqual(metrics.SpikeBaseRadius, 0, ...
                "Spike base radius metric should be non-negative.");
            testCase.verifyTrue(isfield(metrics, "SpikeSeamRadius"), ...
                "Metrics should expose the cone-to-sphere seam radius.");
            testCase.verifyGreaterThanOrEqual(metrics.SpikeSeamRadius, 0, ...
                "Seam radius metric should be non-negative.");
            tipRadiusExpected = 0.5 * urchinStruct.Parameters.SpikeTipDiameter;
            testCase.verifyGreaterThanOrEqual(metrics.SpikeTipRadius, 0, ...
                "Spherical tip radius metric should be non-negative.");
            testCase.verifyEqual(metrics.SpikeTipRadius, tipRadiusExpected, "AbsTol", 1e-9, ...
                "Reported tip radius should match half of the spherical tip diameter.");
            testCase.verifyGreaterThanOrEqual(metrics.TotalVolume, 0, ...
                "Total volume metric should be non-negative.");

            testCase.verifyTrue(isstruct(diagnostics), ...
                "Diagnostics output must be a struct with quality flags.");
            testCase.verifyTrue(isfield(diagnostics, "IsWatertight"), ...
                "Diagnostics should report watertight status.");
        end

        function testDefaultIsDeterministic(testCase)
            % Two identical calls should produce identical geometry
            meshA = urchin(SpikeCount=testCase.SmallSpikeCount);
            meshB = urchin(SpikeCount=testCase.SmallSpikeCount);

            delta = max(abs(meshA.Vertices(:) - meshB.Vertices(:)));
            testCase.verifyLessThanOrEqual(delta, 1e-9, ...
                "Deterministic configuration should produce identical vertices.");
        end

        function testSpikeFluctuationDisabledByDefault(testCase)
            % Default flucFactor must be 0 so there is no jitter unless requested
            meshDefault = urchin(SpikeCount=testCase.SmallSpikeCount);
            meshFluct = urchin(SpikeCount=testCase.SmallSpikeCount, FlucFactor=0.5, FlucMethod="random");

            bboxDefault = [min(meshDefault.Vertices, [], 1); max(meshDefault.Vertices, [], 1)];
            bboxFluct   = [min(meshFluct.Vertices, [], 1); max(meshFluct.Vertices, [], 1)];
            bboxDelta = max(abs(bboxDefault(:) - bboxFluct(:)));

            testCase.verifyGreaterThan(bboxDelta, 1e-6, ...
                "Enabling spike fluctuation should perturb the geometry compared to defaults.");
        end

        function testVolumeMaskExport(testCase)
            % Dense volume export should populate VolumeMask metadata
            mesh = urchin(SpikeCount=12, GenVolume=true, VolResolution=48);

            diagnostics = meshDiagnostics(mesh);
            testCase.verifyTrue(isstruct(diagnostics), "Diagnostics helper should return a struct.");
            testCase.verifyTrue(isfield(mesh, 'VolumeMask'), ...
                "VolumeMask field must be present when genVolume=true.");
            testCase.verifyTrue(~isempty(mesh.VolumeMask), ...
                "VolumeMask should contain voxel data.");
            testCase.verifyEqual(ndims(mesh.VolumeMask), 3, ...
                "VolumeMask must be a 3-D array.");
            testCase.verifyTrue(islogical(mesh.VolumeMask), ...
                "VolumeMask should be logical to represent occupancy.");
            testCase.verifyEqual(mesh.Parameters.VolResolution, 48, ...
                "Parameters struct should reflect the requested VolResolution.");
        end

        function testAdaptiveVolumeLeaves(testCase)
            % Sparse adaptive volume export should expose octree leaves metadata
            mesh = urchin(SpikeCount=12, GenVolume=true, VolAdaptive=true, VolResolution=64);

            testCase.verifyTrue(isfield(mesh, 'VolumeOctree'), ...
                "VolumeOctree field must be populated for adaptive voxelizations.");
            testCase.verifyTrue(isfield(mesh.VolumeOctree, 'Leaves'), ...
                "VolumeOctree should expose Leaves metadata.");
            testCase.verifyGreaterThan(numel(mesh.VolumeOctree.Leaves), 0, ...
                "Adaptive voxelization should produce at least one leaf block.");
            testCase.verifyTrue(mesh.Parameters.VolAdaptive, ...
                "Parameters struct should record VolAdaptive flag.");
        end

        function testMinimalSeamFallbackKeepsGeometryStable(testCase)
            % When tangential seams fall below sampling limits, we still expect
            % a valid mesh thanks to the three-point seam fallback.
            mesh = urchin(SpikeCount=12, SpikeLength=1e-3, Resolution=32);

            testCase.verifyGreaterThan(size(mesh.Vertices, 1), 0, ...
                "Minimal seam fallback should still produce geometry.");

            minSpacing = 2 * (mesh.Parameters.CoreRadius + mesh.Parameters.SpikeLength) / mesh.Parameters.Resolution;
            seamDiameter = 2 * mesh.Metrics.SpikeSeamRadius;

            testCase.verifyLessThan(seamDiameter, minSpacing, ...
                "Fallback scenario should correspond to seams smaller than the sampling threshold.");

            fallbackBase = 0.5 * minSpacing;
            testCase.verifyGreaterThanOrEqual(mesh.Metrics.SpikeBaseRadius, fallbackBase - 1e-9, ...
                "Fallback should keep the base radius at least at the spacing-driven minimum.");
            testCase.verifyGreaterThan(mesh.Metrics.SpikeBaseRadius, 0, ...
                "Base radius should stay finite even when tangency degenerates.");
        end

        function testShortSpikeMaintainsFiniteBaseAndSeam(testCase)
            % Versions 4.1.0 and later clamp tip spheres and keep short spikes from collapsing.
            mesh = urchin(SpikeCount=testCase.SmallSpikeCount, ...
                SpikeLength=0.05, ...
                SpikeTipDiameter=0.5, ...
                Resolution=80);

            minSpacing = 2 * (mesh.Parameters.CoreRadius + mesh.Parameters.SpikeLength) / mesh.Parameters.Resolution;
            minBaseRadius = 0.5 * minSpacing;

            baseRadius = mesh.Metrics.SpikeBaseRadius;
            seamRadius = mesh.Metrics.SpikeSeamRadius;

            testCase.verifyGreaterThan(baseRadius, 0, ...
                "Short spikes should keep a finite base radius.");
            testCase.verifyGreaterThanOrEqual(baseRadius, minBaseRadius - 1e-9, ...
                "Base radius should respect the spacing-driven fallback.");

            testCase.verifyGreaterThanOrEqual(seamRadius, 0, ...
                "Seam radius metric should never be negative.");
            testCase.verifyLessThanOrEqual(seamRadius, baseRadius + 1e-9, ...
                "Seam radius should not exceed the base radius.");
        end
    end
end
