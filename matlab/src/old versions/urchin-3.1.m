function [mesh, mask, threshold, eqRadius] = urchin(varargin)
    % URCHIN  Three-Dimensional Urchin Model Creator
    % =========================================================================
    % This function generates a 3D model of a spherical urchin particle with
    % conical spikes. The spikes' length, tip radius, and orientation can be
    % customized or randomized.
    %
    %   Usage (Name-Value Pairs):
    %       [mesh, mask, threshold, eqRadius] = urchin('cr',10,'sl',5,'ns',50,'st',0.5, ...
    %                                       'sc',0.7,'sf',0.2,'res',1, ...
    %                                       'smth',0.5,'flucMethod','random','distMethod','uniform');
    %
    %   If you call the function without output arguments, it automatically
    %   visualizes the urchin mask and mesh using the viewer3d tool.
    %       urchin(); 
    %       urchin('cr',10,'sl',5);
    %
    % INPUT (Name-Value Pairs):
    %   cr          - Core radius. (default 10)
    %   sl          - Nominal spike length. (default 5)
    %   us          - Urchin size. us = 2*(cr + sl) if provided, missing cr or sl will be derived. if both missing, cr = us/6 and sl = us/3. if both provided, us is ignored. (default [])
    %   ns          - Number of spikes. (default 20)
    %   st          - Spikes tip thickness. Must satisfy st <= cr. (default 0)
    %   sc          - Conicality factor, 0= cylindrical, 1= maximally conical when spike bases filling the entire core. (default 1)
    %   sf          - Spike fluctuation range [0,1]. 0 means no fluctuation. (default 0)
    %   res         - Voxel dimension (smallest discretization unit). (default (cr+sl)/49)
    %   smth        - Smoothing factor (0 for no smoothing | 1 for maximum voxel smoothing preserving spikes | >1 if you want to melt the spikes away!). (default 1)
    %   flucMethod  - 'uniform' (Sobol set) or 'random' for spike length flucs. (default 'uniform')
    %   distMethod  - 'uniform' (Fibonacci) or 'random' for spike orientation. (default 'uniform')
    %
    % OUTPUT:
    %   mesh       - A struct containing the surface mesh faces (F) and vertices (V).
    %   mask       - A 3D float array representing the density mask.
    %   eqRadius   - Equivalent sphere radius of the final structure.
    %   threshold  - Threshold value for mask used to define particle surface.
    %
    % Dependencies: Image Processing Toolbox, Viewer3D
    %
    % Version: 3.1
    % Created by: Maziar Moussavi
    % Enhanced via: Github Copilot
    % Date: 2025-02-07
    % =========================================================================
    
    %% 1) Parse Name-Value Inputs
    p = inputParser;
    p.addParameter('cr',          1,      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('sl',          2,      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('us',         [],      @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('ns',         50,      @(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
    p.addParameter('st',         [],      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('sc',          1,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    p.addParameter('sf',          0,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    p.addParameter('res',        [],      @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('smth',        1,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
    p.addParameter('flucMethod','uniform',@(x)any(validatestring(x,{'uniform','random'})));
    p.addParameter('distMethod','uniform',@(x)any(validatestring(x,{'uniform','random'})));
    
    p.parse(varargin{:});
    
    cr = p.Results.cr;
    sl = p.Results.sl;
    ns = p.Results.ns;
    st = p.Results.st;
    sc = p.Results.sc;
    sf = p.Results.sf;
    res = p.Results.res;
    smth = p.Results.smth;
    flucMethod = p.Results.flucMethod;
    distMethod = p.Results.distMethod;
    us = p.Results.us;
    
    % If 'us' is provided, derive missing cr or sl values.
    if ~isempty(us)
        defaultsUsed = p.UsingDefaults;
        if any(strcmp('cr', defaultsUsed)) && ~any(strcmp('sl', defaultsUsed))
            cr = us/2 - sl;
        elseif ~any(strcmp('cr', defaultsUsed)) && any(strcmp('sl', defaultsUsed))
            sl = us/2 - cr;
        elseif any(strcmp('cr', defaultsUsed)) && any(strcmp('sl', defaultsUsed))
            cr = (1/3) * us/2;
            sl = (2/3) * us/2;
        elseif ~any(strcmp('cr', defaultsUsed)) && ~any(strcmp('sl', defaultsUsed))
            warning('Both core radius and spike length are provided. Ignoring us.');
        end
    end
    
    % If res is empty, set default based on cr + sl.
    if isempty(res)
        res = (cr + sl) / 49;
    end
    % If st is empty, set default based on res.
    if isempty(st)
        st = res;
    end
    %% 2) Validate Some Conditions
    
    % If cr < res, clamp cr to res
    if cr < res
        warning('Core diameter must be >= voxel dimension. Resetting cr to voxel dimension.');
        cr = res/2;
    end
    
    % If st > cr, clamp st to cr
    if st > cr
        warning('Spikes tip radius (st) must be <= core diameter (cr). Resetting st to core diameter.');
        st = cr/2;
    end
    
    % If smth > 1, clamp it.
    if smth > 1
        warning('Smoothing factor > 1. Spike features may be lost.');
    end
    
    %% 3) Distribution of Spikes
    fprintf('Computing spikes geometry...\n');
    tic;
    
    switch distMethod
        case 'uniform'
            % Fibonacci spiral distribution for uniform coverage
            Phi = (1 + sqrt(5)) / 2;  % Golden ratio
            i = (1:ns)';
            theta = acos(1 - 2 .* i ./ (ns + 1));
            phi = mod(2 * pi * i / Phi, 2 * pi);
    
        case 'random'
            % Random distribution on the sphere
            phi = 2 * pi * rand(ns, 1);
            cosTheta = 2 * rand(ns, 1) - 1;  % Ensures no pole clustering
            theta = acos(cosTheta);
    end
    
    % Orientation vectors
    sp = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    
    %% 4) Fluctuations in Spike Length
    switch flucMethod
        case 'uniform'
            flucs = -net(sobolset(1), ns);  % Sobol points in ~[-1, 0]
        case 'random'
            flucs = -rand(ns, 1);          % Random points in [-1,0]
    end
    % Scale them to full spike length and apply fluctuation factor
    flucs = sf * sl * flucs;
    
    %% 5) Spike Geometry Calculations
    % Spherical cap approx to guess max base diameter
    h = 2 * cr / ns;
    SB = 2 * sqrt(2 * cr * h - h^2);  % Max cone base diameter
    
    if st > SB
        warning('Spike tip radius > possible base. Setting st = SB and sc = 0');
        st = SB;
        sc = 0;
    end
    
    Alpha = atan(0.5*(SB - st) / (sl + h));
    alpha = max(sc * Alpha, 1e-3);
    
    RPE = cr + sl + flucs + 0.5*st / tan(alpha); % position of the spike tips
    
    %% 6) Build the 3D Mask
    N = 2 * ceil((cr + sl + smth * 0.5 * st) / res) + 1;
    centerCoord = res * (N + 1) / 2;

    xa = (1:N) * res; xa = xa - mean(xa);
    ya = (1:N) * res; ya = ya - mean(ya);
    za = (1:N) * res; za = za - mean(za);
    [Y, X, Z] = meshgrid(ya, xa, za);
    
    if license('test', 'Distrib_Computing_Toolbox') && gpuDeviceCount > 0
        X = gpuArray(X);
        Y = gpuArray(Y);
        Z = gpuArray(Z);
    else
        warning('Parallel Toolbox not available or no GPU found. Continuing on CPU arrays.');
    end
    
    fprintf('Creating the core...\n');
    mask = (X.^2 + Y.^2 + Z.^2) <= cr^2;
    coreVolume = sum(mask(:)) * res^3;
    
    %% 7) Add Spikes
    fprintf('Adding %d spikes...\n', ns);  % Corrected the format specifier
    for ii = 1:ns
        rpe = RPE(ii);
    
        % Angle between point in space and spike axis
        angleCheck = acos(complex((sp(ii,1)*(X - rpe*sp(ii,1)) + ...
                       sp(ii,2)*(Y - rpe*sp(ii,2)) + ...
                       sp(ii,3)*(Z - rpe*sp(ii,3))) ./ ...
                       sqrt((X - rpe*sp(ii,1)).^2 + (Y - rpe*sp(ii,2)).^2 + (Z - rpe*sp(ii,3)).^2)));
    
        spike = (angleCheck >= (pi - alpha)) & ...                      % within spike half-angle
                ((X - rpe*sp(ii,1)).^2 + (Y - rpe*sp(ii,2)).^2 + ...
                 (Z - rpe*sp(ii,3)).^2) <= (rpe - cr + 2*h)^2 & ...    % limit spike base
                ((X.^2 + Y.^2 + Z.^2) <= (cr + sl + flucs(ii))^2);      % cut excess tip
        
        mask = mask | spike;
    end
    
    volume = sum(mask(:)) * res^3;

    %% 8) Smooth if Requested

    mask = single(mask); % Convert to single for smoothing

    if smth > 0
        fprintf('Applying smoothing...\n');
        kernelSize = 2 * floor(smth * (st / res)) + 1; % Ensure an odd kernel size
        mask = smooth3(mask, 'gaussian', kernelSize); 
    end

    %% 9) Extract Surface Mesh
    threshold = 0.44946787; % default

        fprintf('Shifting mesh border to match volume... \n');
        targetVolume = volume;
        tolerance = 1e-5 * targetVolume;
        thresholdMin = 0; 
        thresholdMax = threshold * 2; 
        maxIter = 50;
        
        for i = 1:maxIter
            disp(['    Iteration: ', num2str(i), ' | Threshold: ', num2str(threshold, '%.5f')]);
            threshold = (thresholdMin + thresholdMax) / 2;
            [F, V] = isosurface(mask, threshold);
            V = V * res - centerCoord;

            % Compute the volume of the mesh
            volume_mesh = 0;
            for j = 1:size(F,1)
                v1 = V(F(j,1), :);
                v2 = V(F(j,2), :);
                v3 = V(F(j,3), :);
                volume_mesh = volume_mesh + dot(v1, cross(v2, v3));
            end
            volume_mesh = abs(volume_mesh) / 6; % Divide by 6 to get the volume of the tetrahedron
            
            if abs(volume_mesh - targetVolume) < tolerance
            break;
            elseif volume_mesh > targetVolume
            thresholdMin = threshold;
            else
            thresholdMax = threshold;
            end
        end

        if nargout ~= 0
            % Create the final mesh with the final threshold
            mesh = surfaceMesh(V, F);
        end
    
    %% 10) Compute Stats
    fprintf('Urchin created successfully.\n');
    elapsedTime = toc;
    fprintf('Elapsed time: %.3f seconds\n', elapsedTime);
    
    eqRadius = (3 * volume / (4 * pi))^(1/3);
    
    fprintf('Urchin Particle Properties:\n');
    fprintf('   Mask dimensions       : %s voxels\n', mat2str(size(mask)));
    fprintf('   Equivalent Sphere R   : %.3f\n', eqRadius);
    fprintf('   Mask Volume           : %.3f\n', volume);
    fprintf('   Mesh Volume           : %.3f\n', volume_mesh); % Added to output mesh volume
    
    alphaDeg = alpha * 180 / pi;
    fprintf('   Spike Half-Angle      : %.3f degrees\n', alphaDeg);
    
    meanSpikeLength = sl + mean(flucs);
    fprintf('   Mean Spike Length     : %.3f\n', meanSpikeLength);
    
    spikeVolume = volume - coreVolume;
    eqSpikeRadius = sqrt(spikeVolume / (ns * pi * meanSpikeLength));
    fprintf('   Mean Spike Thickness  : %.3f\n', 2*eqSpikeRadius);
    
    %% 11) Visualize If No Output
    if nargout == 0
        % create a sliced mesh
        [F, V] = isosurface(mask .* single(X <= 0), threshold); % Use the threshold variable
        V = V * res - centerCoord;  % Ensure correct transformation
        mesh_sliced = surfaceMesh(V, F);

        % create a sliced mask
        tform = affinetform3d([res 0 0 -centerCoord; 0 res 0 -centerCoord; 0 0 res -centerCoord; 0 0 0 1]);
        mask_sliced = gather(mask .* single(X >= 0));

        % Visualize
        viewer = viewer3d;
        surfaceMeshShow(mesh_sliced, Parent=viewer);
        volshow(mask_sliced, 'Parent', viewer, 'RenderingStyle', 'CinematicRendering', Colormap=parula, Transformation=tform);
    end

    clearvars -except mesh mask eqRadius threshold; % To clean up the workspace if running without output
end