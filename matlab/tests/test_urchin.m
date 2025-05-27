classdef test_urchin < matlab.unittest.TestCase
    % Test cases for the urchin function

    properties
        TestDataPath = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'test_outputs'); % Path to save test outputs
    end

    methods (TestClassSetup)
        function createTestDataPath(testCase)
            if ~exist(testCase.TestDataPath, 'dir')
                mkdir(testCase.TestDataPath);
            end
        end
    end

    methods (Test)
        function testDefaultValues(testCase)
            % Test urchin with default parameters
            [mesh, mask, threshold, eqRadius] = urchin();
            testCase.verifyNotEmpty(mesh.Vertices, 'Mesh vertices should not be empty for default run.');
            testCase.verifyNotEmpty(mesh.Faces, 'Mesh faces should not be empty for default run.');
            testCase.verifyTrue(isnumeric(mask) && ndims(mask) == 3, 'Mask should be a 3D numeric array.');
            testCase.verifyTrue(isnumeric(threshold) && isscalar(threshold), 'Threshold should be a numeric scalar.');
            testCase.verifyTrue(isnumeric(eqRadius) && isscalar(eqRadius) && eqRadius > 0, 'Equivalent radius should be a positive numeric scalar.');
            
            % Optional: Save output for inspection
            % save(fullfile(testCase.TestDataPath, 'default_urchin.mat'), 'mesh', 'mask', 'threshold', 'eqRadius');
        end

        function testCustomCoreRadius(testCase)
            % Test urchin with a custom core radius
            cr_val = 5;
            [mesh, ~, ~, eqRadius] = urchin('cr', cr_val);
            testCase.verifyNotEmpty(mesh.Vertices);
            % A simple check: equivalent radius should be related to core radius
            testCase.verifyGreaterThan(eqRadius, cr_val * 0.5, 'Equivalent radius seems too small for the given core radius.');
        end
        
        function testNoSpikes(testCase)
            % Test urchin with zero spikes (should generate just the core)
            ns_val = 0;
            [mesh, mask, ~, eqRadius] = urchin('ns', ns_val, 'sl', 5); % sl is needed to define overall size if ns=0
            testCase.verifyNotEmpty(mesh.Vertices, 'Mesh vertices should not be empty even with no spikes.');
            % Check if the equivalent radius is close to the core radius
            % This assumes 'cr' default is 1. If urchin has a different default, adjust.
            testCase.verifyEqual(eqRadius, feval(@(x) x, urchin('cr')), 'AbsTol', 0.1, ...
                 'Equivalent radius should be close to core radius when no spikes are present.');
        end

        function testHollowCore(testCase)
            % Test with a hollow core
            cr_val = 10;
            hr_val = 5;
            [mesh, mask, ~, eqRadius] = urchin('cr', cr_val, 'hr', hr_val, 'ns', 10); % Added ns for stability
            testCase.verifyNotEmpty(mesh.Vertices);
            testCase.verifyTrue(eqRadius < cr_val && eqRadius > hr_val, 'Equivalent radius for hollow core is out of expected range.');
            % Further checks could involve inspecting the mask center
        end

        function testSpikeLengthFluctuation(testCase)
            % Test with spike length fluctuation
            [mesh_no_fluc, ~, ~, eqR_no_fluc] = urchin('sf', 0, 'ns', 20);
            [mesh_fluc, ~, ~, eqR_fluc] = urchin('sf', 0.8, 'ns', 20);
            testCase.verifyNotEmpty(mesh_fluc.Vertices);
            % It's hard to give a precise test, but fluctuation might change eqRadius slightly
            % or visual inspection of saved meshes would be needed.
            testCase.log(1, sprintf('EqR without fluctuation: %f, with fluctuation: %f', eqR_no_fluc, eqR_fluc));
        end

        function testDifferentDistributionMethods(testCase)
            % Test 'random' and 'uniform' distribution
            [mesh_uniform, ~] = urchin('distMethod', 'uniform', 'ns', 30);
            [mesh_random, ~]  = urchin('distMethod', 'random',  'ns', 30);
            testCase.verifyNotEmpty(mesh_uniform.Vertices);
            testCase.verifyNotEmpty(mesh_random.Vertices);
            % Check that vertex counts are not identical (probabilistically)
            testCase.verifyNotEqual(size(mesh_uniform.Vertices,1), size(mesh_random.Vertices,1), ...
                'Random and Uniform distribution produced meshes with the exact same number of vertices, which is unlikely.');
        end
        
        function testAntialiasingEffect(testCase)
            % Test antialiasing effect (indirectly, by checking mask size or computation time)
            tic;
            [~, mask_no_aa] = urchin('antialiasing', false, 'res', 1, 'cr', 5, 'sl', 2, 'ns', 10);
            time_no_aa = toc;
            tic;
            [~, mask_aa] = urchin('antialiasing', true, 'res', 1, 'cr', 5, 'sl', 2, 'ns', 10);
            time_aa = toc;
            
            testCase.verifyTrue(all(size(mask_aa) > size(mask_no_aa)), ...
                'Antialiased mask should be larger internally before downscaling if antialiasing involves upsampling.');
            testCase.log(1, sprintf('Time without AA: %.2fs, Time with AA: %.2fs', time_no_aa, time_aa));
        end

        function testResolutionImpact(testCase)
            % Test different resolutions
            [mesh_low_res, ~, ~, eqR_low] = urchin('res', 2, 'cr', 10, 'sl', 5, 'ns', 20);
            [mesh_high_res, ~, ~, eqR_high] = urchin('res', 0.5, 'cr', 10, 'sl', 5, 'ns', 20);
            testCase.verifyTrue(size(mesh_high_res.Vertices,1) > size(mesh_low_res.Vertices,1),
                'Higher resolution should result in more mesh vertices.');
            testCase.verifyAlmostEqual(eqR_low, eqR_high, 'RelTol', 0.1, ...
                'Equivalent radius should be similar despite resolution changes.');
        end

        function testSmoothingEffect(testCase)
            % Test smoothing
            [mesh_no_smooth, ~] = urchin('smth', 0, 'cr', 8, 'sl', 4, 'ns', 15);
            [mesh_smooth, ~] = urchin('smth', 1, 'cr', 8, 'sl', 4, 'ns', 15);
            % Smoothing might change the number of vertices or their positions.
            % A robust test would be complex, possibly involving curvature analysis.
            % For now, just ensure it runs and produces a mesh.
            testCase.verifyNotEmpty(mesh_smooth.Vertices);
            testCase.verifyNotEqual(size(mesh_no_smooth.Vertices,1), size(mesh_smooth.Vertices,1), ...
                'Smoothing should ideally alter the mesh, unlikely to have same vertex count.');
        end

        function testEdgeCaseLowSpikeCount(testCase)
            % Test with a very low number of spikes
            [mesh, ~] = urchin('ns', 1, 'cr', 5, 'sl', 5);
            testCase.verifyNotEmpty(mesh.Vertices, 'Should produce a mesh even with 1 spike.');
        end

        function testEdgeCaseHighSpikeCount(testCase)
            % Test with a high number of spikes (might be slow)
            % Mark as slow if there's a mechanism or keep ns moderate for CI
            [mesh, ~] = urchin('ns', 150, 'cr', 10, 'sl', 3, 'res', 1.5); % Adjusted res for speed
            testCase.verifyNotEmpty(mesh.Vertices, 'Should produce a mesh with many spikes.');
        end

        function testParameterUrchinSize_us(testCase)
            % Test 'us' parameter for defining overall size
            us_val = 30;
            [mesh, ~, ~, eqRadius] = urchin('us', us_val, 'ns', 10); % cr and sl should be derived
            testCase.verifyNotEmpty(mesh.Vertices);
            % Expected eqRadius should be somewhat less than us_val/2
            testCase.verifyLessThan(eqRadius, us_val/2, 'Equivalent radius should be less than us/2.');
            testCase.verifyGreaterThan(eqRadius, us_val/8, 'Equivalent radius seems too small for given us.'); % Heuristic lower bound
        end

    end
end
