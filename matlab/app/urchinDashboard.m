function urchinDashboard()
% URCHINDASHBOARD Interactive dashboard to generate, preview, and export nano-urchin meshes.
%   Launches a uifigure with an HTML control panel (uihtml) and a 3D
%   preview pane. The panel exposes all urchin parameters, triggers mesh
%   generation via the urchin() function, and supports STL/export actions.

    rootDir = fileparts(mfilename('fullpath'));
    htmlPath = fullfile(rootDir, 'urchin-dashboard.html');

    fig = uifigure('Name', 'Urchin Dashboard', ...
        'Position', [100 100 1400 900]);
    mainLayout = uigridlayout(fig, [1, 2]);
    mainLayout.ColumnWidth = {440, '1x'};
    mainLayout.RowHeight = {'1x'};
    mainLayout.Padding = [12 12 12 12];
    mainLayout.RowSpacing = 12;
    mainLayout.ColumnSpacing = 12;

    ctrl = uihtml(mainLayout, 'HTMLSource', htmlPath);
    ctrl.Layout.Row = 1;
    ctrl.Layout.Column = 1;
    ctrl.UserData = struct();

    viewerSlot = uipanel(mainLayout);
    viewerSlot.BackgroundColor = [0.08 0.09 0.11];
    viewerSlot.BorderType = 'none';
    viewerSlot.Layout.Row = 1;
    viewerSlot.Layout.Column = 2;
    if isprop(viewerSlot, 'AutoResizeChildren')
        viewerSlot.AutoResizeChildren = 'off';
    end

    % Don't create viewer yet - will be created on first render
    % Store empty viewer reference and panel in figure UserData
    fig.UserData = struct('viewer', [], 'panel', viewerSlot);
    viewerSlot.SizeChangedFcn = @(src, evt) layoutViewerFromFig(fig);
    fig.SizeChangedFcn = @(src, evt) layoutViewerFromFig(fig);

    defaults = defaultParameters();
    ranges = parameterRanges();
    options = parameterOptions();

    uiConfig = struct('defaults', defaults, 'ranges', ranges, 'options', options);
    ctrl.UserData = struct('config', uiConfig, 'lastResult', [], 'volHandle', [], 'meshHandle', []);

    viewHandles = struct('panel', viewerSlot, 'viewer', [], 'fig', fig);
    ctrl.HTMLEventReceivedFcn = @(src, event) handleHtmlEvent(src, event, viewHandles);

    % Don't force initial layout - viewer will be created on first render
    drawnow();

    ctrl.Data = uiConfig; % Makes defaults available to JS via DataChanged
end

function layoutViewerFromFig(fig)
% Helper that reads viewer and panel from fig.UserData and positions viewer
    if isempty(fig) || ~isvalid(fig)
        return;
    end
    if ~isstruct(fig.UserData) || ~isfield(fig.UserData, 'viewer') || ~isfield(fig.UserData, 'panel')
        return;
    end
    view3d = fig.UserData.viewer;
    slotPanel = fig.UserData.panel;
    layoutViewerSlot(slotPanel, view3d);
end

function scheduleViewerLayout(fig)
% Schedule repeated repositioning to counteract async viewer resize
    t = timer('ExecutionMode', 'fixedRate', ...
              'Period', 0.1, ...
              'TasksToExecute', 5, ...
              'TimerFcn', @(~,~) safeLayoutViewer(fig), ...
              'StopFcn', @(obj,~) delete(obj));
    start(t);
end

function safeLayoutViewer(fig)
    try
        if isvalid(fig)
            layoutViewerFromFig(fig);
        end
    catch
        % Ignore errors during async layout
    end
end

function layoutViewerSlot(slotPanel, view3d)
    if isempty(view3d) || ~isvalid(view3d)
        return;
    end
    if isempty(slotPanel) || ~isvalid(slotPanel)
        return;
    end

    prevUnitsSlot = slotPanel.Units;
    slotPanel.Units = 'pixels';
    slotPos = slotPanel.Position; % [x y w h] in pixels relative to figure
    slotPanel.Units = prevUnitsSlot;

    if numel(slotPos) < 4
        return;
    end

    if view3d.Parent == slotPanel
        if isprop(view3d, 'Units')
            view3d.Units = 'normalized';
        end
        if isprop(view3d, 'Position')
            view3d.Position = [0 0 1 1];
        end
        return;
    end

    if isprop(view3d, 'Units')
        view3d.Units = 'pixels';
    end
    if isprop(view3d, 'Position')
        view3d.Position = slotPos;
    end
end

function handleHtmlEvent(src, event, viewHandles)
    state = src.UserData;
    config = state.config;
    debugFlag = debugEnabled(event.HTMLEventData);

    try
        switch string(event.HTMLEventName)
            case "Ready"
                sendEventToHTMLSource(src, "Boot", config);
                % Auto-generate a default urchin on startup for immediate preview
                try
                    params = config.defaults;
                    tic;
                    result = callUrchin(params); %#ok<NASGU>
                    elapsed = toc;
                    state.lastResult = result;
                    state.lastParams = params;
                    src.UserData = state;

                    renderVisualization(src, viewHandles, result);
                    metrics = summarizeResult(result, elapsed);
                    sendEventToHTMLSource(src, "Result", metrics);
                catch autoErr
                    reportError(src, autoErr, debugFlag);
                end

            case "Generate"
                params = mergeParams(config.defaults, event.HTMLEventData);
                tic;
                result = callUrchin(params); %#ok<NASGU>
                elapsed = toc;
                state.lastResult = result;
                state.lastParams = params;
                src.UserData = state;

                renderVisualization(src, viewHandles, result);
                metrics = summarizeResult(result, elapsed);
                sendEventToHTMLSource(src, "Result", metrics);

            case "ExportSTL"
                state = exportStl(state, event, debugFlag);
                src.UserData = state;

            case "ExportConfig"
                exportConfig(event, config.defaults, debugFlag);

            otherwise
                sendEventToHTMLSource(src, "Error", "Unsupported event received.");
        end
    catch ME
        reportError(src, ME, debugFlag);
    end
end

function renderVisualization(src, viewHandles, result)
    state = src.UserData;

    % Delete old handles
    if isfield(state, 'volHandle') && ~isempty(state.volHandle) && isvalid(state.volHandle)
        delete(state.volHandle);
    end
    state.volHandle = [];

    if isfield(state, 'meshHandle') && ~isempty(state.meshHandle) && isvalid(state.meshHandle)
        delete(state.meshHandle);
    end
    state.meshHandle = [];

    % Delete old viewer from fig.UserData (the actual active viewer)
    if isstruct(viewHandles.fig.UserData) && isfield(viewHandles.fig.UserData, 'viewer')
        oldViewer = viewHandles.fig.UserData.viewer;
        if ~isempty(oldViewer) && isvalid(oldViewer)
            delete(oldViewer);
        end
    end

    volData = pickVolumeData(result);
    if ~isempty(volData)
        state.volHandle = volshow(volData, ...
            'Parent', viewHandles.viewer, ...
            'Colormap', [0 0.65 1], ...
            'RenderingStyle', 'MaximumIntensityProjection');
        src.UserData = state;
        return;
    end

    mesh = pickSurfaceMesh(result);

    if isempty(mesh)
        sendEventToHTMLSource(src, "Error", "No surface mesh to display.");
        src.UserData = state;
        return;
    end

    if isa(mesh, 'surfaceMesh')
        mesh = struct('Vertices', mesh.Vertices, 'Faces', mesh.Faces);
    end

    if isstruct(mesh) && isfield(mesh, 'Vertices') && isfield(mesh, 'Faces') && ~isempty(mesh.Faces)
        sm = surfaceMesh(mesh.Vertices, mesh.Faces);
        % Recreate viewer3d to avoid Name property errors with surfaceMeshShow
        if ~isempty(viewHandles.viewer) && isvalid(viewHandles.viewer)
            delete(viewHandles.viewer);
        end
        newViewer = viewer3d(viewHandles.fig, 'BackgroundColor', [0.08 0.09 0.11]);
        if isprop(newViewer, 'Units')
            newViewer.Units = 'pixels';
        end
        layoutViewerSlot(viewHandles.panel, newViewer);
        viewHandles.viewer = newViewer;
        % Update figure UserData so SizeChangedFcn uses the new viewer
        viewHandles.fig.UserData.viewer = newViewer;
        viewHandles.fig.UserData.panel = viewHandles.panel;
                % Match visualization style from urchin.m
                try
                    % Compute scale from parameters if available, else from vertices
                    scale = [];
                    if isstruct(result) && isfield(result, 'Parameters')
                        p = result.Parameters;
                        if isfield(p, 'CoreRadius') && isfield(p, 'SpikeLength')
                            scale = 2 * (double(p.CoreRadius) + double(p.SpikeLength));
                        end
                    end
                    if isempty(scale) && isstruct(mesh) && isfield(mesh, 'Vertices') && ~isempty(mesh.Vertices)
                        r = vecnorm(double(mesh.Vertices), 2, 2);
                        scale = 2 * max(r);
                        if ~isfinite(scale) || scale <= 0
                            scale = [];
                        end
                    end
                    if ~isempty(scale) && isprop(newViewer, 'CameraPosition')
                        newViewer.CameraPosition = [0 0 -scale];
                    end
                    if isprop(newViewer, 'CameraUpVector')
                        newViewer.CameraUpVector = [0 1 0];
                    end
                    if isprop(newViewer, 'Lighting')
                        newViewer.Lighting = "off";
                    end
                catch
                    % Continue without styling if properties not available
                end

        % Call surfaceMeshShow
        try
            surfaceMeshShow(sm, 'Parent', newViewer);
        catch
            % Fallback: call without Parent
            surfaceMeshShow(sm);
        end
        % Reposition after surfaceMeshShow in case it changed size
                % Call surfaceMeshShow with translucent fill and wireframe overlay
                try
                    surfaceMeshShow(sm, 'Parent', newViewer, 'Alpha', 0.5);
                catch
                    surfaceMeshShow(sm, 'Alpha', 0.5);
                end
                try
                    surfaceMeshShow(sm, 'Parent', newViewer, 'WireFrame', true);
                catch
                    surfaceMeshShow(sm, 'WireFrame', true);
                end
    end

    src.UserData = state;
                scheduleViewerLayout(viewHandles.fig);
end

function volData = pickVolumeData(result)
    volData = [];
    if ~isstruct(result)
        return;
    end

    if isfield(result, 'VolumeMask') && ~isempty(result.VolumeMask)
        volData = result.VolumeMask;
    end

end

function mesh = pickSurfaceMesh(result)
    if isstruct(result) && isfield(result, 'SurfaceMesh')
        mesh = result.SurfaceMesh;
        return;
    end

    mesh = struct();
    if isstruct(result) && all(isfield(result, {'Vertices', 'Faces'}))
        mesh.Vertices = result.Vertices;
        mesh.Faces = result.Faces;
        if isempty(mesh.Vertices) || isempty(mesh.Faces)
            mesh = [];
        end
    end
end

function metrics = summarizeResult(result, elapsed)
    metrics = struct();
    metrics.Status = "ok";
    metrics.ElapsedSeconds = elapsed;

    mesh = pickSurfaceMesh(result);
    metrics.Vertices = safeSize(mesh, 'Vertices');
    metrics.Faces = safeSize(mesh, 'Faces');

    if isfield(result, 'Parameters')
        params = result.Parameters;
        metrics.SpikeCount = getFieldOr(params, 'SpikeCount', NaN);
        metrics.CoreRadius = getFieldOr(params, 'CoreRadius', NaN);
        metrics.SpikeLength = getFieldOr(params, 'SpikeLength', NaN);
        metrics.SpikeTipDiameter = getFieldOr(params, 'SpikeTipDiameter', NaN);
        metrics.SpikeConicality = getFieldOr(params, 'SpikeConicality', NaN);
    end

    if isfield(result, 'Metrics')
        metrics.MinSpacing = getFieldOr(result.Metrics, 'MinSpacing', NaN);
        metrics.MaxBaseRadius = getFieldOr(result.Metrics, 'SpikeBaseMaxima', NaN);
        metrics.TotalVolume = getFieldOr(result.Metrics, 'TotalVolume', NaN);
        metrics.TipCollapsed = getFieldOr(result.Metrics, 'SpikeTipCollapsed', false);
    end
end

function n = safeSize(mesh, fieldName)
    if isempty(mesh) || ~isfield(mesh, fieldName) || ~isnumeric(mesh.(fieldName))
        n = 0;
        return;
    end
    n = size(mesh.(fieldName), 1);
end

function val = getFieldOr(s, name, default)
    if isstruct(s) && isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end

function tf = debugEnabled(payload)
    tf = false;
    if nargin < 1 || isempty(payload)
        return;
    end
    if isstruct(payload) && isfield(payload, 'Debug')
        tf = logical(payload.Debug);
    end
end

function reportError(src, ME, debugFlag)
    msg = ME.message;
    if debugFlag && ~isempty(ME.stack)
        trace = arrayfun(@(s) sprintf('%s:%d', s.name, s.line), ME.stack, 'UniformOutput', false);
        msg = sprintf('%s | stack: %s', ME.message, strjoin(trace, ' > '));
    end
    sendEventToHTMLSource(src, "Error", msg);
end

function result = callUrchin(params)
    % Convert struct to name/value cell list for urchin()
    names = fieldnames(params);
    values = struct2cell(params);
    nv = reshape([names.'; values.'], 1, []);
    result = urchin(nv{:});
end

function params = mergeParams(defaults, payload)
    params = defaults;
    if isempty(payload) || ~isstruct(payload)
        return;
    end

    userFields = fieldnames(payload);
    for k = 1:numel(userFields)
        key = userFields{k};
        if isfield(defaults, key)
            params.(key) = coerceValue(payload.(key), defaults.(key));
        end
    end
end

function value = coerceValue(inputValue, defaultValue)
    if islogical(defaultValue)
        value = logical(inputValue);
    elseif isnumeric(defaultValue)
        if isempty(inputValue)
            value = defaultValue;
        else
            value = double(inputValue);
        end
    elseif isstring(defaultValue) || ischar(defaultValue)
        value = string(inputValue);
    else
        value = inputValue;
    end
end

function defaults = defaultParameters()
    defaults = struct(...
        'CoreRadius', 30, ...
        'SpikeLength', 15, ...
        'SpikeCount', 120, ...
        'SpikeTipDiameter', 3, ...
        'SpikeConicality', 0.5, ...
        'FilletRatio', 0.25, ...
        'Resolution', 100, ...
        'UseFillet', true, ...
        'IncludeCore', true, ...
        'FlucFactor', 0.5, ...
        'FlucMethod', "uniform", ...
        'DistMethod', "uniform", ...
        'RefinedOrientation', true, ...
        'RefinedOrientationThresholdDeg', 0.1, ...
        'GenVolume', false, ...
        'VolResolution', 128, ...
        'VolPadding', 0.05, ...
        'VolAlpha', [], ...
        'VolAdaptive', false, ...
        'VolDxMax', [], ...
        'VolDxMin', [], ...
        'VolBlockSize', 8, ...
        'VolCriterion', "boundary" ...
    );
end

function ranges = parameterRanges()
    ranges = struct();
    ranges.CoreRadius = [0.1, 300];
    ranges.SpikeLength = [0.1, 300];
    ranges.SpikeCount = [0, 500];
    ranges.SpikeTipDiameter = [0, 50];
    ranges.SpikeConicality = [-1, 1];
    ranges.FilletRatio = [0, 0.5];
    ranges.Resolution = [10, 400];
    ranges.FlucFactor = [0, 1];
    ranges.RefinedOrientationThresholdDeg = [0, 10];
    ranges.VolResolution = [16, 512];
    ranges.VolPadding = [0, 0.5];
    ranges.VolDxMax = [0.001, 50];
    ranges.VolDxMin = [0.0001, 10];
    ranges.VolBlockSize = [4, 32];
end

function opts = parameterOptions()
    opts = struct();
    opts.FlucMethod = {"uniform", "random", "gaussian"};
    opts.DistMethod = {"uniform", "random"};
    opts.VolCriterion = {"boundary", "distance", "curvature", "hybrid"};
end

function state = exportStl(state, event, debugFlag)
    if ~isfield(state, 'lastResult') || isempty(state.lastResult)
        sendEventToHTMLSource(event.Source, "Error", "Generate an urchin before exporting.");
        return;
    end

    [fileName, filePath] = uiputfile('*.stl', 'Export urchin surface as STL');
    if isequal(fileName, 0)
        return;
    end

    mesh = pickSurfaceMesh(state.lastResult);
    if isempty(mesh)
        sendEventToHTMLSource(event.Source, "Error", "No surface mesh available for export.");
        return;
    end

    target = fullfile(filePath, fileName);
    try
        writeSurfaceMesh(mesh, target, 'Encoding', 'binary');
        sendEventToHTMLSource(event.Source, "Exported", struct('path', target));
    catch ME
        reportError(event.Source, ME, debugFlag);
    end
end

function exportConfig(event, defaults, debugFlag)
    params = mergeParams(defaults, event.HTMLEventData);
    [fileName, filePath] = uiputfile('*.json', 'Save urchin parameter set');
    if isequal(fileName, 0)
        return;
    end

    try
        payload = jsonencode(params, 'PrettyPrint', true);
        target = fullfile(filePath, fileName);
        fid = fopen(target, 'w');
        if fid < 0
            error('exportConfig:IO', 'Unable to write configuration file.');
        end

        cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fwrite(fid, payload, 'char');
        sendEventToHTMLSource(event.Source, "Exported", struct('path', target));
    catch ME
        reportError(event.Source, ME, debugFlag);
    end
end
