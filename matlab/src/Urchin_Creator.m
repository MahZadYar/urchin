%% Urchin Creator — comprehensive driver for urchin.m
% This orchestration script exposes every public parameter accepted by
% `urchin.m`, captures both surface and volume outputs, and offers a single
% place to configure visualization and export behaviour (STL, MAT, JSON, …).
%
% Edit the configuration structs below to explore different geometries or
% batch-export datasets. Leave any field set to [] to let `urchin.m` fall back
% to its internal default for that parameter.

%% Geometry and spike layout -------------------------------------------------
geometry = struct( ...
    'rCore', 1, ...           % Core radius (nm)
    'spikeLength', 0.2, ...          % Spike length measured from the core surface (nm)
    'spikeCount', 100, ...         % Number of spikes
    'spikeTip', 0.2, ...           % [] → use default spherical tip diameter (spikeLength/10)
    'spikeConicality', 1 ...          % Conicality: 0 = cylinder, 1 = widest base
    );

%% Surface quality and refinement controls ----------------------------------
surfaceControls = struct( ...
    'filletRatio', 0.25, ...   % Toroidal fillet radius relative to seam diameter
    'resolution', 50, ...    % Dimensionless resolution controlling mesh spacing
    'useFillet', true, ...     % Include toroidal core↔spike fillet (future toggles)
    'includeCore', true ...    % Keep the trimmed core surface in the final mesh
    );

%% Spike stochasticity -------------------------------------------------------
stochastic = struct( ...
    'flucFactor', 1, ...          % Spike length fluctuation factor (0 = no fluctuation)
    'flucMethod', 'uniform', ... % 'uniform' | 'random' | 'gaussian'
    'distMethod', 'uniform', ...  % 'uniform' | 'random'
    'refinedOrientation', true, ... % Relax spike orientations using electrostatic settling
    'refinedOrientationThreshold', 0.1 ... % Minimum angular separation target (degrees)
    );

%% Volume representation controls -------------------------------------------
volumeControls = struct( ...
    'genVolume', false, ...    % Enable volume generation (dense or adaptive)
    'volRes', 128, ...         % Target voxels along longest axis (dense volume)
    'volPadding', 0.05, ...    % Fractional AABB padding for voxelization
    'volAlpha', [], ...        % AlphaShape fallback parameter (auto when [])
    'volAdaptive', false, ...  % true = adaptive octree, false = dense grid
    'volDxMax', [], ...        % Max voxel size for adaptive volume (auto when [])
    'volDxMin', [], ...        % Min voxel size for adaptive volume (auto when [])
    'volBlockSize', 8, ...     % Leaf block resolution for adaptive volume
    'volCriterion', 'boundary' ... % Refinement criterion: 'boundary' | 'distance' | 'curvature' | 'hybrid'
    );

%% Consolidate parameter structs --------------------------------------------
urchinParams = mergeStructs(geometry, surfaceControls, stochastic, volumeControls);

%% Visualization preferences -------------------------------------------------
visualize = struct( ...
    'printDiagnostics', false, ...       % Display diagnostic flags and mesh stats
    'showSurfaceMesh', true, ...        % Visualize the surface mesh
    'surfaceMeshWireframe', true, ...  % Render surface mesh in wireframe mode
    'surfaceMeshAlpha', 0.75, ...       % Opacity for solid rendering (0-1)
    'showVolumeMask', false, ...        % Visualize dense volume mask (requires genVolume=true & volAdaptive=false)
    'annotateSpikeInfo', true, ...          % Show spike indices and lengths on the mesh
    'volumeRenderingStyle', "CinematicRendering" ... % volshow rendering style
    );

%% Export preferences --------------------------------------------------------
exports = struct( ...
    'outputDir', fullfile(pwd, "exports"), ... % Output directory (created automatically)
    'label', [], ...                  % Optional base filename label (auto when [])
    'saveSurfaceMesh', true, ...      % Write surface mesh to disk
    'surfaceMeshFormat', "stl", ...  % 'stl' | 'ply' (limited by writeSurfaceMesh)
    'surfaceMeshEncoding', "binary", ... % STL encoding ('ascii' or 'binary')
    'saveMeshMat', false, ...          % Save mesh struct + diagnostics as MAT file
    'saveVolumeMask', false, ...      % Save dense volume mask (requires genVolume & ~volAdaptive)
    'saveVolumeMaskFormat', "mat", ... % Currently 'mat' supported
    'saveVolumeOctree', false, ...    % Save adaptive octree (requires genVolume & volAdaptive)
    'saveDiagnosticsJson', false, ...  % Export diagnostics to JSON
    'saveDiagnosticsMat', false, ...  % Export diagnostics to MAT
    'saveConfigurationJson', true, ...% Export run configuration to JSON
    'saveConfigurationMat', false ... % Export run configuration to MAT
    );

%% Derived configuration checks ---------------------------------------------
if exports.saveVolumeMask && exports.saveVolumeOctree
    error('Urchin_Creator:conflictExports', ...
        'Dense volume masks and adaptive octrees require separate runs. Disable one of the export toggles.');
end

if (visualize.showVolumeMask || exports.saveVolumeMask) && ~urchinParams.genVolume
    warning('Urchin_Creator:enablingDenseVolume', ...
        'Enabling dense volume generation to honour volume mask visualization/export request.');
    urchinParams.genVolume = true;
end

if exports.saveVolumeOctree && (~urchinParams.genVolume || ~urchinParams.volAdaptive)
    warning('Urchin_Creator:enablingAdaptiveVolume', ...
        'Enabling adaptive volume generation to honour octree export request.');
    urchinParams.genVolume = true;
    urchinParams.volAdaptive = true;
end

if exports.saveVolumeMask && urchinParams.volAdaptive
    error('Urchin_Creator:denseMaskRequiresNonAdaptive', ...
        'Dense volume masks are only produced when volAdaptive=false. Adjust configuration.');
end

if visualize.showVolumeMask && urchinParams.volAdaptive
    warning('Urchin_Creator:noDenseVolume', ...
        'volAdaptive=true -> no dense mask to visualise. Disable showVolumeMask or switch volAdaptive=false.');
end

%% Invoke urchin -------------------------------------------------------------
nvPairs = struct2NameValue(urchinParams);

fprintf('Launching urchin with %d spikes, rCore=%.3g nm, spikeLength=%.3g nm, spikeTip=%.3g nm, spikeConicality=%.3g...\n', ...
    urchinParams.spikeCount, urchinParams.rCore, urchinParams.spikeLength, urchinParams.spikeTip, urchinParams.spikeConicality);
urchinStruct = urchin(nvPairs{:});

needDiagnostics = visualize.printDiagnostics || ...
    exports.saveMeshMat || exports.saveDiagnosticsJson || exports.saveDiagnosticsMat;
if needDiagnostics
    diagnostics = meshDiagnostics(urchinStruct);
else
    diagnostics = [];
end

%% Diagnostics ---------------------------------------------------------------
if visualize.printDiagnostics
    fprintf('\nSurface mesh summary:\n');
    fprintf('  Vertices: %d\n', size(urchinStruct.Vertices, 1));
    fprintf('  Faces   : %d\n', size(urchinStruct.Faces, 1));
    if isfield(urchinStruct, 'VolumeMask') && ~isempty(urchinStruct.VolumeMask)
        fprintf('  Volume mask size: %s\n', mat2str(size(urchinStruct.VolumeMask)));
    end
    if isfield(urchinStruct, 'VolumeOctree') && ~isempty(urchinStruct.VolumeOctree)
        fprintf('  Adaptive volume leaves: %d\n', numel(urchinStruct.VolumeOctree.Leaves));
    end
    fprintf('\nQuality diagnostics:\n');
    disp(struct2table(diagnostics));
end

%% Visualisation ------------------------------------------------------------
if visualize.showSurfaceMesh && isfield(urchinStruct, 'SurfaceMesh')
    if exist('viewer3d', 'file') == 2 && exist('surfaceMeshShow', 'file') == 2
        viewer = viewer3d;
        if visualize.surfaceMeshAlpha
            surfaceMeshShow(urchinStruct.SurfaceMesh, Parent=viewer, Alpha=visualize.surfaceMeshAlpha);
        end
        if visualize.surfaceMeshWireframe
            surfaceMeshShow(urchinStruct.SurfaceMesh, Parent=viewer, WireFrame=true);
        end
        if visualize.annotateSpikeInfo
            annotateSpikeApexPoints(viewer, urchinStruct);
        end
    else
        figure('Name', 'Urchin Surface Mesh');
        if visualize.surfaceMeshWireframe
            edgeColor = [0 0 0];
            faceAlpha = visualize.surfaceMeshAlpha;
        else
            edgeColor = 'none';
            faceAlpha = visualize.surfaceMeshAlpha;
        end
    trisurf(urchinStruct.Faces, urchinStruct.Vertices(:,1), urchinStruct.Vertices(:,2), urchinStruct.Vertices(:,3), ...
            'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'FaceColor', [0.8 0.8 0.9]);
        axis equal off; lighting gouraud; camlight;
        title('Urchin Surface Mesh (viewer3d unavailable)');
    end
end

if visualize.showVolumeMask
    if isfield(urchinStruct, 'VolumeMask') && ~isempty(urchinStruct.VolumeMask)
        if exist('volshow', 'file') == 2
            volshow(urchinStruct.VolumeMask, Colormap=parula, RenderingStyle=visualize.volumeRenderingStyle);
        else
            warning('Urchin_Creator:volshowUnavailable', ...
                'volshow is unavailable. Save the mask and inspect it manually (e.g., sliceViewer).');
        end
    end
end

%% Exports ------------------------------------------------------------------
exportRequests = [exports.saveSurfaceMesh, exports.saveMeshMat, exports.saveVolumeMask, ...
    exports.saveVolumeOctree, exports.saveDiagnosticsJson, exports.saveDiagnosticsMat, ...
    exports.saveConfigurationJson, exports.saveConfigurationMat];

if any(exportRequests)
    if ~isfolder(exports.outputDir)
        mkdir(exports.outputDir);
    end
    if isempty(exports.label)
        exports.label = autoLabel(urchinParams);
    end
end

savedFiles = strings(0, 1);

if exports.saveSurfaceMesh
    meshFile = fullfile(exports.outputDir, sprintf('%s.%s', exports.label, exports.surfaceMeshFormat));
    writeSurfaceMesh(urchinStruct.SurfaceMesh, meshFile, "Encoding", exports.surfaceMeshEncoding);
    savedFiles(end+1,1) = string(meshFile); %#ok<SAGROW>
end

if exports.saveMeshMat
    matPath = fullfile(exports.outputDir, sprintf('%s_mesh.mat', exports.label));
    meshStruct = urchinStruct; %#ok<NASGU>
    diagnosticsStruct = diagnostics; %#ok<NASGU>
    urchinParameters = urchinParams; %#ok<NASGU>
    save(matPath, 'meshStruct', 'diagnosticsStruct', 'urchinParameters');
    savedFiles(end+1,1) = string(matPath); %#ok<SAGROW>
end

if exports.saveVolumeMask
    if isfield(urchinStruct, 'VolumeMask') && ~isempty(urchinStruct.VolumeMask)
        switch lower(exports.saveVolumeMaskFormat)
            case 'mat'
                maskPath = fullfile(exports.outputDir, sprintf('%s_volumeMask.mat', exports.label));
                maskPayload = struct( ...
                    'VolumeMask', urchinStruct.VolumeMask, ...
                        'VoxelGrid', urchinStruct.VoxelGrid, ...
                        'VoxelSize', urchinStruct.VoxelSize, ...
                        'Bounds', urchinStruct.Bounds);
                save(maskPath, 'maskPayload');
                savedFiles(end+1,1) = string(maskPath); %#ok<SAGROW>
            otherwise
                error('Urchin_Creator:unsupportedMaskFormat', ...
                    'Unsupported volume mask format: %s', exports.saveVolumeMaskFormat);
        end
    else
        warning('Urchin_Creator:noDenseMask', ...
            'saveVolumeMask enabled but urchinStruct.VolumeMask is empty. Ensure genVolume=true and volAdaptive=false.');
    end
end

if exports.saveVolumeOctree
    if isfield(urchinStruct, 'VolumeOctree') && ~isempty(urchinStruct.VolumeOctree)
        octreePath = fullfile(exports.outputDir, sprintf('%s_volumeOctree.mat', exports.label));
    volumeOctree = urchinStruct.VolumeOctree; %#ok<NASGU>
        save(octreePath, 'volumeOctree');
        savedFiles(end+1,1) = string(octreePath); %#ok<SAGROW>
    else
        warning('Urchin_Creator:noOctree', ...
            'saveVolumeOctree enabled but no VolumeOctree data was produced.');
    end
end

if exports.saveDiagnosticsJson
    jsonPath = fullfile(exports.outputDir, sprintf('%s_diagnostics.json', exports.label));
    fid = fopen(jsonPath, 'w');
    fprintf(fid, '%s', jsonencode(diagnostics));
    fclose(fid);
    savedFiles(end+1,1) = string(jsonPath); %#ok<SAGROW>
end

if exports.saveDiagnosticsMat
    diagMatPath = fullfile(exports.outputDir, sprintf('%s_diagnostics.mat', exports.label));
    diagnosticsStruct = diagnostics; %#ok<NASGU>
    save(diagMatPath, 'diagnosticsStruct');
    savedFiles(end+1,1) = string(diagMatPath); %#ok<SAGROW>
end

if exports.saveConfigurationJson
    cfgPath = fullfile(exports.outputDir, sprintf('%s_configuration.json', exports.label));
    cfgPayload = struct('urchinParams', urchinParams, 'visualize', visualize, 'exports', exports);
    fid = fopen(cfgPath, 'w');
    fprintf(fid, '%s', jsonencode(cfgPayload));
    fclose(fid);
    savedFiles(end+1,1) = string(cfgPath); %#ok<SAGROW>
end

if exports.saveConfigurationMat
    cfgMatPath = fullfile(exports.outputDir, sprintf('%s_configuration.mat', exports.label));
    cfgPayload = struct('urchinParams', urchinParams, 'visualize', visualize, 'exports', exports); %#ok<NASGU>
    save(cfgMatPath, 'cfgPayload');
    savedFiles(end+1,1) = string(cfgMatPath); %#ok<SAGROW>
end

if ~isempty(savedFiles)
    fprintf('\nArtifacts written to %s:\n', exports.outputDir);
    for idx = 1:numel(savedFiles)
        fprintf('  %s\n', savedFiles(idx));
    end
end

%% Helper functions ---------------------------------------------------------
% Utility helpers are now provided as standalone functions:
%   mergeStructs.m, struct2NameValue.m, autoLabel.m, valueToken.m

function annotateSpikeApexPoints(viewer, urchinStruct)
%ANNOTATESPIKEApexPOINTS Place 3D point annotations at spike apices when requested.
    if isempty(viewer) || ~isgraphics(viewer) || nargin < 2
        return;
    end
    if exist('images.ui.graphics3d.roi.Point', 'class') ~= 8
        warning('Urchin_Creator:pointAnnotationUnavailable', ...
            'images.ui.graphics3d.roi.Point is unavailable. Spike annotations skipped.');
        return;
    end
    if ~isfield(urchinStruct, 'Parameters')
        return;
    end
    params = urchinStruct.Parameters;
    requiredFields = {'rCore', 'spikeLengths', 'spikeOrientations'};
    for idx = 1:numel(requiredFields)
        if ~isfield(params, requiredFields{idx}) || isempty(params.(requiredFields{idx}))
            return;
        end
    end
    orientations = params.spikeOrientations;
    lengths = params.spikeLengths(:);
    if size(orientations, 2) ~= 3 || isempty(lengths)
        return;
    end
    nSpikes = min(size(orientations,1), numel(lengths));
    orientations = orientations(1:nSpikes, :);
    lengths = lengths(1:nSpikes);
    apexRadius = params.rCore + lengths;
    apexPositions = orientations .* apexRadius;
    pointArray = images.ui.graphics3d.roi.Point.empty;
    pointArray(1, nSpikes) = images.ui.graphics3d.roi.Point;
    for k = 1:nSpikes
        label = sprintf('Spike #%d (%.3g nm)', k, lengths(k));
        pointArray(k) = images.ui.graphics3d.roi.Point( ...
            Parent=viewer, Label=label, Position=apexPositions(k, :));
    end
    viewer.Annotations = pointArray;
end

