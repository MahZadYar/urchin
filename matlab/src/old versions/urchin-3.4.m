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
    %   hr          - Hollow radius. (default 0)    
    %   sl          - Nominal spike length from core. (default 5)
    %   us          - Urchin size. us = 2*(cr + sl) if provided, missing cr or sl will be derived. if both missing, cr = us/6 and sl = us/3. if both provided, us is ignored. (default [])
    %   ns          - Number of spikes. (default 20)
    %   st          - Spikes tip thickness. Must satisfy st <= cr. (default 0)
    %   sc          - Conicality factor, 0= cylindrical, 1= maximally conical when spike bases filling the entire core. (default 1)
    %   sf          - Spike fluctuation range [0,1]. 0 means no fluctuation. (default 0)
    %   res         - Voxel dimension (smallest discretization unit). (default (cr+sl)/49)
    %   smth        - Smoothing factor (0 for no smoothing | 1 for maximum voxel smoothing preserving spikes | >1 if you want to melt the spikes away!). (default 1)
    %   flucMethod  - 'uniform' (Sobol set) or 'random' for spike length flucs. (default 'uniform')
    %   distMethod  - 'uniform' (Fibonacci) or 'random' for spike orientation. (default 'uniform')
    %   antialising - If true, upscales the resolution to avoid aliasing. (default true)
    %
    % OUTPUT:
    %   mesh       - A struct containing the surface mesh faces (F) and vertices (V).
    %   mask       - A 3D float array representing the density mask.
    %   eqRadius   - Equivalent sphere radius of the final structure.
    %   threshold  - Threshold value for mask used to define particle surface.
    %
    % Dependencies: Image Processing Toolbox, Viewer3D
    %
    % Version: 3.4
    % Created by: Maziar Moussavi
    % Enhanced via: Github Copilot
    % Date: 2025-05-07
    % =========================================================================

    tic;

    %% 1) Parse Name-Value Inputs
    p = inputParser;
    p.addParameter('cr',          1,      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('hr',         [],      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('sl',          2,      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('us',         [],      @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('ns',         100,      @(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
    p.addParameter('st',         [],      @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('sc',          0.75,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',-1,'<=',1}));
    p.addParameter('sf',          0.75,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    p.addParameter('res',        [],      @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('smth',        1,      @(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
    p.addParameter('flucMethod','uniform',@(x)any(validatestring(x,{'uniform','random'})));
    p.addParameter('distMethod','uniform',@(x)any(validatestring(x,{'uniform','random'})));
    p.addParameter('antialiasing', true,  @(x)validateattributes(x,{'logical'},{'scalar'}));
    
    p.parse(varargin{:});
    
    cr = p.Results.cr;
    hr = p.Results.hr;
    sl = p.Results.sl;
    ns = p.Results.ns;
    st = p.Results.st;
    sc = p.Results.sc;
    sf = p.Results.sf;
    res = p.Results.res;
    smth = p.Results.smth;
    flucMethod = p.Results.flucMethod;
    distMethod = p.Results.distMethod;
    antialiasing = p.Results.antialiasing; % Updated variable name
    us = p.Results.us;

    %% 2) Validate and initialize parameters    

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
        st = 2*res;
    end

    if isempty(antialiasing) % Updated variable name to match the new parameter
        antialiasing = true; % Default to no antialiasing
    end
    
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
    
    fprintf('Creating urchin with the following parameters:\n');
    fprintf('   Core radius (cr)      : %.3f\n', cr);
    fprintf('   Hollow radius (hr)    : %.3f\n', hr);
    fprintf('   Spike length (sl)     : %.3f\n', sl);
    fprintf('   Urchin size (us)      : %.3f\n', us);
    fprintf('   Num spikes (ns)       : %d\n', ns);
    fprintf('   Spike tip (st)        : %.3f\n', st);
    fprintf('   Conicality (sc)       : %.3f\n', sc);
    fprintf('   Fluctuation (sf)      : %.3f\n', sf);
    fprintf('   Resolution (res)      : %.3f\n', res);
    fprintf('   Smoothing (smth)      : %.3f\n', smth);
    fprintf('   Fluc method           : %s\n', flucMethod);
    fprintf('   Dist method           : %s\n', distMethod);
    fprintf('   Antialiasing          : %s\n', mat2str(logical(antialiasing)));

    %% 5) Setup Coordinates

    if antialiasing
        antialiasing = 2; % Upscaling factor for antialiasing
        res = res / antialiasing; % upscaling the resolution to avoid aliasing
    end

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

    %% 6) Geometry Calculations
    fprintf('Computing spikes geometry...\n');

    % Distribution of Spikes
    switch distMethod
        case 'uniform'
            % Fibonacci spiral distribution for uniform coverage
            Phi = (1 + sqrt(5)) / 2;  % Golden ratio (can be replaced with any other irrational value)
            i = (1:ns)';
            theta = acos(1 - 2 .* i ./ (ns + 1));
            phi = mod(2 * pi * i / Phi, 2 * pi);
    
        case 'random'
            % Random distribution on the sphere
            phi = 2 * pi * rand(ns, 1);
            cosTheta = 2 * rand(ns, 1) - 1;  % Ensures no pole clustering
            theta = acos(cosTheta);
    end
    
    % Spikes' orientation vectors
    sp = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    
    % Fluctuations in Spike Length
    switch flucMethod
        case 'uniform'
            flucs = -net(sobolset(1), ns);  % Sobol points in ~[-1, 0]
        case 'random'
            % Reproducible seed based only on the number of spikes
            seed = ns;
            rng(seed, 'twister');    % seed RNG deterministically from ns
            flucs = -rand(ns, 1);    % random points in [-1,0]
    end
    % Scale them to full spike length and apply fluctuation factor
    flucs = sf * sl * flucs;
    
    % Spike Cone Calculation
    % ns spherical caps approximately fill the core surface. 
    if ns > 1
        % Sum of CapArea = 2 * pi * cr * hCap * ns = 4 * pi * cr^2 (sum of caps surface area = core surface area)
        hCap = 2 * cr / ns;   % Height of the spherical cap 
    else
        hCap = cr;           % If only one spike, its cap can only fill half of the core
    end
    rCap = sqrt(2 * cr * hCap - hCap^2);  % Max cone base radius is the radius of the spherical cap
    rTip = st/2;
    x = sl - rTip + hCap;

    if sc >= 0 % Conical spikes        
        if st > 2*rCap
            warning('Spike tip is bigger than possible base. neglecting sc parameter.');
            sc = 1;
        end

        % Maximum possible spike half-angle (the cone base encapsulates the maximum sppheical cap and tangent to the tip sphere)
        % sin(Alpha) = rTip/d   &    tan(Alpha) = rCap/(x+rTip)   =>
        Alpha = asin((-x*rTip + rCap*sqrt(x^2 + rCap^2 - rTip^2))/(x^2 + rCap^2)); 
    else % reversed conical spikes
        % Maximum possible spike half-angle (the cone's tip is on the core and tangent to the tip sphere)
        Alpha = asin(rTip / x); 
    end

    alpha = sc * Alpha; % Spike half-angle

    if abs(alpha) < 1e-3
        alpha = 1e-3; % Avoid numerical issues with very small angles
    end

    RPE = cr + sl + flucs + 0.5 * st * (-1 + 1 / sin(alpha)); % distances of the extended cones' tips from the center    
    SCR = cr + sl + flucs + rTip * (sin(alpha) - 1);  % distances from the center to the tangent point of the cone and the tip sphere
    TipPoints = sp .* RPE; % Cone tip points of the spikes
    TipCenters = sp .* (cr + sl + flucs - rTip); % sphere center moved rTip back along the spike axis from spike tip    

    %% 7) Creating the Particle Mask
    fprintf('Creating the core...\n');
    mask = (X.^2 + Y.^2 + Z.^2) <= cr^2;    
    
    fprintf('Adding %d spikes...\n', ns);
    for ii = 1:ns
        % 1) create the cone
        % Angle between spike axis and point in space relative to cone tip 
        angleCheck = acos(complex( ...
            ( sp(ii,1)*(X - TipPoints(ii, 1)) + ...
              sp(ii,2)*(Y - TipPoints(ii, 2)) + ...  % Updated index from tipPoint(1) to tipPoint(2)
              sp(ii,3)*(Z - TipPoints(ii, 3)) ) ./ ...  % Updated index from tipPoint(1) to tipPoint(3)
            sqrt((X - TipPoints(ii, 1)).^2 + (Y - TipPoints(ii, 2)).^2 + (Z - TipPoints(ii, 3)).^2) ...  % Updated rpe*sp(ii,1) to tipPoint(1)
        )); 
        coneMask = (angleCheck >= (pi - abs(alpha))) | (angleCheck <= abs(alpha));        % Include both cones in either direction along the spike axis

        % 2) create tip sphere
        tipMask  = (X - TipCenters(ii,1)).^2 + (Y - TipCenters(ii,2)).^2 + (Z - TipCenters(ii,3)).^2 ...
                   <= (rTip)^2;

        % 3) limit spike base to core
        baseCut = (X - TipCenters(ii,1)).^2 + (Y - TipCenters(ii,2)).^2 + (Z - TipCenters(ii,3)).^2 ...
                  <= rCap^2 + (hCap + sl + flucs(ii) - rTip)^2; % calculating spike base cut radius from tip center to core intersection
        
        % 4) cut off cone tip beyound the tangent point of the tip sphere
        tipCut  = (X.^2 + Y.^2 + Z.^2) ...
                  <= SCR(ii)^2;

        % 5) combine masks
        spike = (coneMask & baseCut & tipCut) | tipMask ;
        mask  = mask | spike;
    end
    
    %% 8) Hollow Core if Requested
    if ~isempty(p.Results.hr) && p.Results.hr > 0
        fprintf('Creating hollow core...\n');
        mask = mask & (X.^2 + Y.^2 + Z.^2 >= (hr)^2);
    end

    %% 9) Smooth if Requested

    mask = single(mask); % Convert to single for smoothing

    if smth > 0
        fprintf('Applying smoothing...\n');
        kernelSize = 2 * floor(smth) + 1; % Ensure an odd kernel size
        mask = smooth3(mask, 'gaussian', kernelSize); % options: 'gaussian', 'box', 'average'
    else
        mask = gather(mask); % Ensure mask is on CPU for further processing
    end

    %% Calculate Volume
    fprintf('Calculating volumes...\n');
    volume = sum(mask(:)) * res^3;
    coreVolume = sum(mask((X.^2 + Y.^2 + Z.^2) <= (cr)^2)) * res^3;

    %% Scaledown to original resolution
    if antialiasing
        fprintf('Downscaling to original resolution...\n');
        mask = imresize3(mask, 1/antialiasing, 'linear'); % options: 'linear', 'nearest', 'cubic'
        X = imresize3(gather(X), 1/antialiasing, 'linear');
        res = res * antialiasing; % downscale the resolution
    end

    %% 10) Extract Surface Mesh
    fprintf('Shifting mesh border to match volume... \n');
    targetVolume = volume;
    tolerance = 1e-5 * targetVolume;
    threshold = 0.44946787; % default
    thresholdMin = 0; 
    thresholdMax = threshold * 2; 
    maxIter = 50;
    SmoothIterations = 5; % Number of iterations for smoothing

    for i = 1:maxIter
        disp(['    Iteration: ', num2str(i), ' | Threshold: ', num2str(threshold, '%.5f')]);
        threshold = (thresholdMin + thresholdMax) / 2;
        [F, V] = isosurface(mask, threshold); % options: 'linear', 'nearest', 'cubic'
        V = V * res - centerCoord;

        % smooth the mesh        
        if SmoothIterations > 0
            mesh = surfaceMesh(V, F);
            mesh = smoothSurfaceMesh(mesh, SmoothIterations, Method="Average"); % Use numIterations for smoothing
            V = mesh.Vertices;
            F = mesh.Faces;
        end

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
        mesh = surfaceMesh(V, F);
    end

    % %% Extract point cloud from the mesh
    % pc = pointCloud(V);

    %% 11) Compute Stats
    eqRadius = (3 * volume / (4 * pi))^(1/3);
    alphaDeg = alpha * 180 / pi;
    meanSpikeLength = sl + mean(flucs);
    spikeVolume = volume - coreVolume;
    eqSpikeRadius = sqrt(spikeVolume / (ns * pi * meanSpikeLength));
    elapsedTime = toc;

    fprintf('Urchin created successfully!\n');
    fprintf('Elapsed time: %.3f seconds\n', elapsedTime);    
    fprintf('Urchin Particle Properties:\n');
    fprintf('   Mask dimensions       : %s voxels\n', mat2str(size(mask)));
    fprintf('   Equivalent Sphere R   : %.3f\n', eqRadius);
    fprintf('   Mask Volume           : %.3f\n', volume);
    fprintf('   Mesh Volume           : %.3f\n', volume_mesh);    
    fprintf('   Spike Half-Angle      : %.3f degrees\n', alphaDeg);    
    fprintf('   Mean Spike Length     : %.3f\n', meanSpikeLength);
    fprintf('   Mean Spike Thickness  : %.3f\n', 2*eqSpikeRadius);

    
    %% 12) Visualize If No Output
    if nargout == 0
        % create a sliced mesh
        [F, V] = isosurface(mask .* single(X <= 0), threshold); % Use the threshold variable
        V = V * res - centerCoord;  % Ensure correct transformation
        mesh_sliced = surfaceMesh(V, F);
        if SmoothIterations > 0
            mesh_sliced = smoothSurfaceMesh(mesh_sliced, SmoothIterations, Method="Average"); % Use numIterations for smoothing
        end

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