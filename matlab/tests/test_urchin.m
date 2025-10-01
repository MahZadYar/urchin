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
            % Default call should yield a surfaceMesh with diagnostics
            [mesh, diagnostics] = urchin('ns', testCase.SmallSpikeCount);

            testCase.verifyTrue(isstruct(mesh) && isfield(mesh, "SurfaceMesh"), ...
                "urchin should return a struct exposing the underlying surface mesh.");
            testCase.verifyTrue(isa(mesh.SurfaceMesh, "surfaceMesh"), ...
                "SurfaceMesh field must be a surfaceMesh instance.");
            testCase.verifyGreaterThan(size(mesh.Vertices, 1), 0, ...
                "Default run should produce vertices.");
            testCase.verifyGreaterThan(size(mesh.Faces, 1), 0, ...
                "Default run should produce faces.");

            testCase.verifyTrue(isstruct(diagnostics), ...
                "Diagnostics output must be a struct with quality flags.");
            testCase.verifyTrue(isfield(diagnostics, "IsWatertight"), ...
                "Diagnostics should report watertight status.");
        end

        function testDefaultIsDeterministic(testCase)
            % Two identical calls should produce identical geometry
            [meshA, ~] = urchin('ns', testCase.SmallSpikeCount);
            [meshB, ~] = urchin('ns', testCase.SmallSpikeCount);

            delta = max(abs(meshA.Vertices(:) - meshB.Vertices(:)));
            testCase.verifyLessThanOrEqual(delta, 1e-9, ...
                "Deterministic configuration should produce identical vertices.");
        end

        function testSpikeFluctuationDisabledByDefault(testCase)
            % Default sf must be 0 so there is no jitter unless requested
            [meshDefault, ~] = urchin('ns', testCase.SmallSpikeCount);
            [meshFluct, ~] = urchin('ns', testCase.SmallSpikeCount, 'sf', 0.5, 'flucMethod', 'random');

            delta = max(abs(meshDefault.Vertices(:) - meshFluct.Vertices(:)));
            testCase.verifyGreaterThan(delta, 1e-6, ...
                "Enabling spike fluctuation should perturb the geometry compared to defaults.");
        end

        function testVolumeMaskExport(testCase)
            % Dense volume export should populate VolumeMask metadata
            [mesh, diagnostics] = urchin('ns', 12, 'genVolume', true, 'volRes', 48);

            testCase.verifyTrue(isstruct(diagnostics), "Diagnostics should still be returned.");
            testCase.verifyTrue(isfield(mesh, 'VolumeMask'), ...
                "VolumeMask field must be present when genVolume=true.");
            testCase.verifyTrue(~isempty(mesh.VolumeMask), ...
                "VolumeMask should contain voxel data.");
            testCase.verifyEqual(ndims(mesh.VolumeMask), 3, ...
                "VolumeMask must be a 3-D array.");
            testCase.verifyTrue(islogical(mesh.VolumeMask), ...
                "VolumeMask should be logical to represent occupancy.");
        end

        function testAdaptiveVolumeLeaves(testCase)
            % Sparse adaptive volume export should expose octree leaves metadata
            mesh = urchin('ns', 12, 'genVolume', true, 'volAdaptive', true, 'volRes', 64);

            testCase.verifyTrue(isfield(mesh, 'VolumeOctree'), ...
                "VolumeOctree field must be populated for adaptive voxelizations.");
            testCase.verifyTrue(isfield(mesh.VolumeOctree, 'Leaves'), ...
                "VolumeOctree should expose Leaves metadata.");
            testCase.verifyGreaterThan(numel(mesh.VolumeOctree.Leaves), 0, ...
                "Adaptive voxelization should produce at least one leaf block.");
        end
    end
end
