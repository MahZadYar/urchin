function urchin = urchin(varargin)
% URCHIN  Deterministic, curvature-aware B-Rep nano-urchin mesher
% =========================================================================
% Builds a COMSOL-ready surface by stitching parametric patches:
%   - Core sphere with circular holes at spike footprints
%   - Conical frusta for spike bodies (smoothly controlled by spikeConicality∈[0,1])
%   - Spherical caps for spike tips (equal-angle ring spacing)
%   - Toroidal fillets blending core↔cone seams (G1-tangent, geometric)
%
% Point sampling is deterministic and curvature-aware:
%   - Spherical patches via Fibonacci lattice or geodesic icosphere
%   - Cones/fillets via structured grids in parameter space
%   - Seam edges share identical, precomputed rings for watertight joins
%
% New functionality:
%   - Spike fluctuations: per-spike length jitter (flucFactor) with uniform/random (Sobol or seeded RNG).
%   - Orientation distributions: Fibonacci-like uniform or random (seeded).
%   - Resolution-driven sampling: a single dimensionless knob controls
%     all curvature-aware feature counts while keeping seams watertight.
%   - Volume representations:
%       • Dense uniform voxel grid via inpolyhedron/alphaShape fallback
%       • Adaptive sparse octree (VDB-style) with user control of voxel
%         size range (volDxMin/Max), leaf block size, and refinement
%         criterion (boundary/distance/curvature/hybrid). Stored in
%         mesh.VolumeOctree for sparse workflows; optional densification.
% Usage:
%   urchinStruct = urchin('rCore',30,'spikeLength',15,'spikeCount',50,'spikeTip',5);
%   urchinStruct = urchin('rCore',25,'spikeLength',10,'spikeCount',75);
%   diagnostics = meshDiagnostics(urchinStruct.SurfaceMesh);
%   fprintf('Watertight: %d\n', diagnostics.IsWatertight);
%
% Inputs (name-value):
%   rCore   Core radius (nm). Default 1
%   spikeLength  Spike length from core surface (nm). Default 1
%   spikeCount   Number of spikes. Default 100
%   spikeTip     Diameter of the spherical tip cap (nm).
%                Defaults to spikeLength/10 when omitted.
%   spikeConicality Conicality [-1..1]. -1 collapses base, 0=cylinder, 1=widest base.
%   resolution   Dimensionless mesh resolution (default 100). Smallest feature
%                length is ~2*(rCore+spikeLength)/resolution.
%   useFillet    Enable toroidal fillet patch. Default true (fillet radius
%                defaults to spikeTip/2 and is clamped relative to base radius)
%   flucFactor   Spike length fluctuation factor [0..1]. Default 0.5
%   flucMethod   'uniform' (Sobol) or 'random' (seeded by spikeCount). Default 'uniform'
%   distMethod   'uniform' (Fibonacci-like) or 'random' (seeded by spikeCount). Default 'uniform'
%   refinedOrientation Toggle iterative Coulomb-like relaxation (default true)
%   refinedOrientationThreshold Minimum allowable angular separation in degrees (default 0.1)
%
% Volume (dense):
%   genVolume    If true, generate dense voxel mask. Default false
%   volRes       Target voxels along largest dimension. Default 128
%   volPadding   Fractional AABB padding for voxelization. Default 0.05
%   volAlpha     Alpha for alphaShape fallback (auto if empty)
%
% Volume (adaptive sparse / VDB-style):
%   volAdaptive  If true, build sparse octree volume instead of dense mask
%   volDxMax     Max voxel size (auto: maxDim/128)
%   volDxMin     Min voxel size (auto: volDxMax/4)
%   volBlockSize Leaf block resolution (voxel dimension), e.g., 8 or 16
%   volCriterion 'boundary' | 'distance' | 'curvature' | 'hybrid' (boundary by default)
%
% Output:
%   mesh struct with fields: Vertices [n×3], Faces [m×3]
%     + Parameters (struct echoing all name-value inputs; collapsed spike tips report 0)
%     + Metrics (struct with min spacing, spike base radius, enclosed volume, etc.)
%     + VolumeMask (logical 3D) and VoxelGrid axes if genVolume=true & !volAdaptive
%     + VolumeOctree (sparse leaves) if genVolume=true & volAdaptive=true
%
% Dependencies: 
%   sobolset requires Statistics and Machine Learning ToolboxR
%   gpuDeviceCount requires Parallel Computing Toolbox
%   wrapToPi requires Mapping Toolbox
%   surfaceMesh requires Lidar Toolbox
%
% =========================================================================

    %% 1) Input Parser
    p = inputParser;
    addParameter(p, 'rCore', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'spikeLength', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'spikeCount', 100, @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0}));
    addParameter(p, 'spikeTip', [],  @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    addParameter(p, 'spikeConicality', 0.5, @(x)validateattributes(x,{'numeric'},{'scalar','>=',-1,'<=',1}));
    addParameter(p, 'filletRatio', 0.25, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',0.5}));
    addParameter(p, 'resolution', 100, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'useFillet', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'includeCore', true, @(x)islogical(x) && isscalar(x));
    % Spike fluctuations and orientation distribution
    addParameter(p, 'flucFactor', 0.5, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    addParameter(p, 'flucMethod', 'uniform', @(x)any(validatestring(x,{'uniform','random','gaussian'})));
    addParameter(p, 'distMethod', 'uniform', @(x)any(validatestring(x,{'uniform','random'})));
    addParameter(p, 'refinedOrientation', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'refinedOrientationThreshold', 0.1, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',180}));
    % Optional volume representation controls (voxel grid from mesh)
    addParameter(p, 'genVolume', false, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'volRes', 128, @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',16}));
    addParameter(p, 'volPadding', 0.05, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',0.5}));
    addParameter(p, 'volAlpha', [], @(x)isempty(x) || (isscalar(x) && isnumeric(x) && x>0));
    % Adaptive (VDB-style) sparse octree volume options
    addParameter(p, 'volAdaptive', false, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'volDxMax', [], @(x)isempty(x) || (isscalar(x) && isnumeric(x) && x>0));
    addParameter(p, 'volDxMin', [], @(x)isempty(x) || (isscalar(x) && isnumeric(x) && x>0));
    addParameter(p, 'volBlockSize', 8, @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',4}));
    addParameter(p, 'volCriterion', 'boundary', @(x)any(validatestring(x,{'boundary','distance','curvature','hybrid'})));
    parse(p, varargin{:});

    % Assign parsed inputs to variables
    rCore = p.Results.rCore;
    spikeLength = p.Results.spikeLength;
    spikeCount = p.Results.spikeCount;
    tipDiameter = p.Results.spikeTip;
    spikeTipRequested = tipDiameter;
    spikeConicality = p.Results.spikeConicality;
    filletRatio = p.Results.filletRatio;
    resolution = p.Results.resolution;
    useFillet   = p.Results.useFillet;
    includeCore = p.Results.includeCore;
    flucFactor  = p.Results.flucFactor;
    flucMethod  = p.Results.flucMethod;
    distMethod  = p.Results.distMethod;
    refinedOrientation = p.Results.refinedOrientation;
    refinedOrientationThresholdDeg = p.Results.refinedOrientationThreshold;
    genVolume   = p.Results.genVolume;
    volRes      = p.Results.volRes;
    volPadding  = p.Results.volPadding;
    volAlphaInp = p.Results.volAlpha;
    volAdaptive = p.Results.volAdaptive;
    volDxMax    = p.Results.volDxMax;
    volDxMin    = p.Results.volDxMin;
    volBlockSz  = p.Results.volBlockSize;
    volCriterion= p.Results.volCriterion;

    fprintf('Starting B-Rep Urchin Generation...\n');
    tic;

    %% 2) Initial Geometric Calculations

    scale = 2 * (rCore + spikeLength); % urchin scale for normalization
    minSpacing = max(1e-9,  scale/resolution);
    coreSpacingFactor = 2.0;
    azimuthSpacingFactor = 1.2;
    axialSpacingFactor = 1.0;
    capSpacingFactor = 1.0;
    collapseSpacingTol = minSpacing;
    % conicalityFlatThreshold = 0.05;
    conicalityEffective = spikeConicality;
    % forceCylinder = false;
    PointyTip = false;

    % Tip sphere diameter (default: spikeLength/10 unless specified)
    if isempty(tipDiameter)
        tipDiameter = spikeLength / 10;
    end

    if tipDiameter <= collapseSpacingTol & tipDiameter > 0
        tipDiameter = collapseSpacingTol;
        if conicalityEffective < 0
            % negative conicality handled later via tapered base rule
            conicalityEffective = 0;
        end
    elseif tipDiameter == 0
        % Only collapse tip if explicitly requested
        PointyTip = true;
    end

    rTip = tipDiameter / 2;

    forceCylinder = abs(conicalityEffective) <= 1e-12;
    pendingPointyTip = PointyTip;

    tipCollapseTol = max(0.5 * minSpacing, 1e-9);
    if forceCylinder
        tipCollapseTol = min(tipCollapseTol, 0.25 * minSpacing);
    end

    %% Spike Orientations
    switch distMethod
        case 'uniform'
            % Golden spiral distribution across sphere via (theta,phi)
            Phi = (1 + sqrt(5)) / 2;  % golden ratio
            i = (1:spikeCount)';
            theta = acos(1 - 2 .* i ./ (spikeCount + 1));
            phi   = mod(2 * pi * i / Phi, 2 * pi);
        case 'random'
            rng(spikeCount,'twister');
            phi = 2 * pi * rand(spikeCount, 1);
            cosTheta = 2 * rand(spikeCount, 1) - 1;
            theta = acos(cosTheta);
    end

    spikeOrientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)]; % #KeyParameter

    if refinedOrientation && spikeCount > 1
        spikeOrientations = refineSpikeOrientations(spikeOrientations, refinedOrientationThresholdDeg);
    end
    %% Spike length fluctuations
    switch flucMethod
        case 'uniform'
            flucs = net(sobolset(1), spikeCount);   % values in ~[0,1]
        case 'random'
            rng(spikeCount,'twister');
            flucs = rand(spikeCount,1);            % values in ~[0,1]
        case 'gaussian'
            rng(spikeCount,'twister');
            flucs = randn(spikeCount,1);            % values in (-∞,∞)
    end

    flucs = flucs - min(flucs);     % shift to [0,∞)
    flucs = flucs / mean(flucs);   % normalize to mean on 1
    flucs = flucs -1;               % shift to [-1,∞) with mean 0
    flucs = flucFactor * spikeLength * flucs;        % scale to full spike length range [-flucFactor*spikeLength, +∞) with mean 0
    
    spikeLengths = spikeLength + flucs;
    spikeLengths = max(spikeLengths, 0);
    zTipApexs = rCore + spikeLengths; % #KeyParameter
    

    %% Spike maximum base radii (non-overlapping cones & tangency limits)
    if spikeCount == 0
        thetaMins = zeros(0,1);
    elseif spikeCount == 1
        thetaMins = pi;
    else
        % Compute maximum cosine (excluding self) to find minimal angular sep per spike
        Ddot = spikeOrientations * spikeOrientations.'; % spikeCount x spikeCount dot products
        Ddot(1:spikeCount+1:end) = -Inf; % ignore diagonal
        maxDotPerSpike = max(Ddot, [], 2);
        maxDotPerSpike(~isfinite(maxDotPerSpike)) = -1;
        maxDotPerSpike = min(1, max(-1, maxDotPerSpike));
        thetaMins = acos(maxDotPerSpike); % smallest angle between each spike and its nearest neighbour
    end

    if isempty(thetaMins)
        thetaMinGlobal = pi;
    else
        thetaMinGlobal = min(thetaMins);
    end

    alphaBaseMaxSpikes = 0.5 * thetaMins;           % tangent half-angle per spike

    alphaBaseMaxTangentRatios = (rCore - rTip) / (rCore + spikeLengths - rTip);
    alphaBaseMaxTangentRatios = min(1, max(-1, alphaBaseMaxTangentRatios));
    alphaBaseMaxTangents = acos(alphaBaseMaxTangentRatios);

    alphaBaseMaxs = min(alphaBaseMaxSpikes, alphaBaseMaxTangents);      % clamp to tangent angle per spike  % #KeyParameter

    rBaseMaxs = rCore * sin(alphaBaseMaxs);      % maximal non-overlapping base radius per spike #KeyParameter

    if ~isempty(rBaseMaxs)
        rBaseMaxs = min(rBaseMaxs, rCore - 0.5 * minSpacing);
        rBaseMaxs = max(rBaseMaxs, minSpacing);
    end

    % rBaseMax = rCore * sin(0.5 * thetaMinGlobal);
    % rBaseMax = min(rBaseMax, rCore - 0.5 * minSpacing);
    % rBaseMax = max(rBaseMax, minSpacing);

    %% Spike tip radii and centers
    rTips = ones(spikeCount,1) * rTip;
    rTipMaxs = zTipApexs./(1+1./sin(alphaBaseMaxs)); % max tip radius for tangency
    rTips = min(rTipMaxs, rTips);  % #KeyParameter
    zTipCenters = zTipApexs - rTips;  % #KeyParameter
    
    % tipCollapseTol = max(0.5 * minSpacing, 1e-9);
    % if forceCylinder
    %     tipCollapseTol = min(tipCollapseTol, 0.25 * minSpacing);
    % end
    % spikeTipEffectiveDiameter = rTip*2;
    % spikeTipSphereRadius = rTip;

    %% Spike base radii and centers
    baseCollapseTol = 0.5 * minSpacing; % collapse tolerance for base radius
    if forceCylinder
        baseCollapseTol = min(baseCollapseTol, 0.25 * minSpacing);
    end

    if conicalityEffective >= 0
        % Positive conicality: linear interpolation between tip and max base
        rBases = rTips + conicalityEffective * (rBaseMaxs - rTips);   % #KeyParameter
    else
        % Negative conicality: tapered base rule (base radius < tip radius)
        rBases = (1 + conicalityEffective) * rTips; % linear tapering to zero at conicality=-1
    end

    rBases = min(rBases, rBaseMaxs);
    rBases = max(rBases, 0);
    % collapseBase = rBases <= baseCollapseTol;
    % rBases(collapseBase) = 0;

    alphaBases = asin(rBases/rCore);  % #KeyParameter
    
    zBases = sqrt(max(0, rCore^2 - rBases.^2));  % #KeyParameter

    %% Tip Seam Radii and centers
    [zTipSeams, rTipSeams, alphaCones] = solveConeSeams(rBases, rTips, zBases, zTipCenters);   % #KeyParameter

    %% more preparations before mesh generation if needed?
    
    % e.g. number of rings/segments along seams and tip caps and number of radial segments, etc.

    %% 3) Build core sphere
    % Core sphere subdivision target driven purely by resolution
    targetCoreSpacing = coreSpacingFactor * minSpacing;
    targetCoreSpacing = max(targetCoreSpacing, 1e-9*scale);
    coreSurfaceArea = 4 * pi * rCore^2;
    targetVertices = max(12, ceil(coreSurfaceArea / (targetCoreSpacing^2)));
    if targetVertices <= 12
        coreSubdiv = 0;
    else
        coreSubdiv = ceil(0.5 * log2((targetVertices - 2) / 10));
        coreSubdiv = max(coreSubdiv, 0);
    end
    V = zeros(0,3); F = zeros(0,3);
    if includeCore
        coreSubdiv = max(floor(coreSubdiv), 0);
        [coreVFull, coreFFull] = triangulate_icosphere(rCore, coreSubdiv);
        [V, coreIdx] = append_vertices(V, coreVFull); % core vertices global indices
        coreFCurrent = coreFFull + coreIdx(1) - 1;   % active core faces (global indices)
        coreDirsFull = coreVFull / rCore;                % unit directions for initial removal test
        coreMask = false(size(V,1),1);
        coreMask(coreIdx) = true;
    else
        coreIdx = zeros(0,1);
        coreFCurrent = zeros(0,3);
        coreDirsFull = zeros(0,3);
        coreMask = false(0,1);
    end
    %% 4) Generate spikes one-by-one with per-spike core trimming & stitching
    fprintf('Generating & stitching %d spikes...\n', spikeCount);
    minBaseRadiusMetric = Inf;
    minSeamRadiusMetric = Inf;
    tipCollapsedAll = true;
    producedSpikeCount = 0;
    seamIndices = cell(spikeCount,1);

    if isempty(rTips)
        spikeTipSphereRadius = 0;
    else
        finiteTips = rTips(isfinite(rTips));
        if isempty(finiteTips)
            spikeTipSphereRadius = 0;
        else
            spikeTipSphereRadius = max(finiteTips);
        end
    end
    spikeTipEffectiveDiameter = 2 * spikeTipSphereRadius;
    for i = 1:spikeCount
        orientation = spikeOrientations(i, :);
        % Local orthonormal frame:
        %  - orientation: spike axis (unit vector)
        %  - u, v: span orthogonal plane (used for circular rings)
        [u, v] = plane_vectors(orientation);

        coreLoop = [];
        seamIndices{i} = [];

        spikeLength = spikeLengths(i);
        if ~isfinite(spikeLength)
            spikeLength = 0;
        end
        spikeLength = max(spikeLength, 0);

        rTip = rTips(i);
        if ~isfinite(rTip)
            rTip = 0;
        end
        rTip = max(0, min(rTip, spikeLength));

        zTipCenter = zTipCenters(i);
        if ~isfinite(zTipCenter)
            zTipCenter = rCore + spikeLength - rTip;
        end

        if rTip > spikeTipSphereRadius
            spikeTipSphereRadius = rTip;
            spikeTipEffectiveDiameter = 2 * spikeTipSphereRadius;
        end

        rBaseMax = rBaseMaxs(i);
        if ~isfinite(rBaseMax)
            rBaseMax = 0;
        end
        rBaseMax = max(rBaseMax, 0);

        rBase = rBases(i);
        if ~isfinite(rBase)
            rBase = 0;
        end
        rBase = max(0, min(rBase, rBaseMax));

        minBaseRadius = baseCollapseTol;
        if rBase > 0 && rBase < minBaseRadius
            rBase = min(minBaseRadius, max(rBase, rBaseMax));
        end
        if rBase >= rCore
            rBase = 0.999 * rCore;
        end

        collapseBaseI = rBase <= baseCollapseTol;
        if collapseBaseI
            rBase = 0;
        end

        cosCutoffI = 1;
        if ~collapseBaseI && rCore > 0
            ratio = min(1, max(rBase / rCore, 0));
            cosCutoffI = sqrt(max(0, 1 - ratio.^2));
        end

        zBase = zBases(i);
        if ~isfinite(zBase)
            zBase = sqrt(max(0, rCore^2 - rBase^2));
        end
        zBase = max(zBase, 0);

        if zBase < minSpacing
            tipCollapsedAll = false;
            continue;
        end

        if includeCore && ~collapseBaseI
            coreIdx = find(coreMask);
            if isempty(coreIdx)
                coreDirsFull = zeros(0,3);
            else
                dirs = V(coreIdx,:);
                norms = vecnorm(dirs,2,2);
                norms(norms < 1e-12) = 1;
                coreDirsFull = dirs ./ norms;
            end

            coreDot = coreDirsFull * orientation';
            newRemoveMask = coreDot > cosCutoffI + 1e-12;
            if ~any(newRemoveMask) && ~isempty(coreIdx)
                [~, fallbackPos] = max(coreDot);
                if ~isnan(fallbackPos) && fallbackPos >= 1 && fallbackPos <= numel(newRemoveMask)
                    newRemoveMask(fallbackPos) = true;
                end
            end
            if any(newRemoveMask)
                vertsRemoveGlobal = coreIdx(newRemoveMask);
                coreMask(vertsRemoveGlobal) = false;
                faceRemoveMask = any(ismember(coreFCurrent, vertsRemoveGlobal), 2);
                removedFaces = coreFCurrent(faceRemoveMask, :);
                coreFCurrent(faceRemoveMask, :) = [];

                vertsFromRemovedFaces = unique(removedFaces(:));
                coreLoop = setdiff(vertsFromRemovedFaces, vertsRemoveGlobal);

                candidates = {};
                if numel(coreLoop) >= 3
                    candidates{end+1} = coreLoop(:)'; %#ok<AGROW>
                end
                if ~isempty(vertsRemoveGlobal)
                    for j = 1:i-1
                        ring = seamIndices{j};
                        if numel(ring) < 3, continue; end
                        if all(ismember(ring, vertsRemoveGlobal))
                            candidates{end+1} = ring(:)'; %#ok<AGROW>
                        end
                    end
                end
                if ~isempty(candidates)
                    bestDot = -Inf;
                    bestLoop = [];
                    for c = 1:numel(candidates)
                        loop = candidates{c};
                        ctr = mean(V(loop, :), 1);
                        nrm = norm(ctr);
                        if nrm < 1e-12, continue; end
                        dval = dot(ctr / nrm, orientation);
                        if dval > bestDot
                            bestDot = dval;
                            bestLoop = loop;
                        end
                    end
                    coreLoop = bestLoop;
                end
            end
        end

        seamBaseCenter = zBase * orientation;
        baseSegmentsEffective = 3;
        if collapseBaseI
            [V, seamIndices{i}] = append_vertices(V, seamBaseCenter);
        else
            baseSegmentsI = ring_segment_count(rBase, minSpacing, azimuthSpacingFactor);
            if forceCylinder
                baseSegmentsI = max(baseSegmentsI, 3);
            end
            baseSegmentsEffective = max(3, baseSegmentsI);
            anglesBase = circle_angles(baseSegmentsEffective);
            seamBaseRing = seam_ring_points(seamBaseCenter, u, v, rBase, anglesBase);
            [V, seamIndices{i}] = append_vertices(V, seamBaseRing);
            if includeCore
                coreMask(seamIndices{i}) = true;
            end

            if includeCore && numel(coreLoop) >= 3
                seamRing = seamIndices{i}(:)';
                newBridgeFaces = bridge_rings_idx(coreLoop, seamRing, V, orientation);
                coreFCurrent = [coreFCurrent; newBridgeFaces]; %#ok<AGROW>
                coreMask(coreLoop) = true;
            end

            if includeCore
                coreIdx = find(coreMask);
                if isempty(coreIdx)
                    coreDirsFull = zeros(0,3);
                else
                    dirs = V(coreIdx,:);
                    norms = vecnorm(dirs,2,2);
                    norms(norms < 1e-12) = 1;
                    coreDirsFull = dirs ./ norms;
                end
            end
        end

        if ~includeCore && ~collapseBaseI
            baseCenter = seamBaseCenter;
            [V, baseCenterIdx] = append_vertices(V, baseCenter);
            baseFaces = connect_ring_to_apex(seamIndices{i}, baseCenterIdx);
            if ~isempty(baseFaces)
                tri = baseFaces(1,:);
                v1 = V(tri(1),:);
                v2 = V(tri(2),:);
                v3 = V(tri(3),:);
                baseNormal = cross(v2 - v1, v3 - v1);
                if dot(baseNormal, -orientation) < 0
                    baseFaces = baseFaces(:, [1 3 2]);
                end
            end
            F = [F; baseFaces]; %#ok<AGROW>
        end

        forcePointyTipI = pendingPointyTip && ~collapseBaseI && (rBase > 2 * minSpacing);
        if forceCylinder
            forcePointyTipI = false;
        end

        minimalSeam = false;
        alphaCone = alphaCones(i);
        if ~isfinite(alphaCone)
            alphaCone = alphaBases(i);
        end
        alphaCone = max(1e-6, alphaCone);

        zTipSeam = zTipSeams(i);
        if ~isfinite(zTipSeam)
            zTipSeam = zTipCenter + rTip * cos(alphaCone);
        end

        rTipSeam = rTipSeams(i);
        if ~isfinite(rTipSeam)
            rTipSeam = min(max(min(rBase, rTip), 0), rTip);
        end
        rTipSeam = max(rTipSeam, 0);

        if rTip <= 0
            rTipSeam = 0;
            alphaCone = pi/2;
            zTipSeam = max(zBase + 1e-6, zTipCenter);
        else
            % Re-solve the seam with the sanitized radii to preserve tangency.
            if ~collapseBaseI && rBase > 0
                try
                    [zTipSeamLoc, rTipSeamLoc, alphaConeLoc] = solveConeSeams(rBase, rTip, zBase, zTipCenter);
                catch
                    zTipSeamLoc = NaN;
                    rTipSeamLoc = NaN;
                    alphaConeLoc = NaN;
                end
                if isfinite(zTipSeamLoc) && isfinite(rTipSeamLoc) && isfinite(alphaConeLoc)
                    zTipSeam = zTipSeamLoc;
                    rTipSeam = rTipSeamLoc;
                    alphaCone = alphaConeLoc;
                end
            end

            rTipSeam = min(rTipSeam, min(rBase, rTip));
            ratio = min(1, max(rTipSeam / max(rTip, 1e-12), 0));
            alphaCone = max(1e-6, asin(ratio));
            zTipSeam = max(zBase + 1e-6, zTipCenter + rTip * cos(alphaCone));

            if rTipSeam > 0 && (2 * rTipSeam < minSpacing)
                minimalSeam = true;
            end
        end

        rTipSeam = min(rTipSeam, rBase);
        tipCapHeight = rTip * max(0, cos(alphaCone));
        tipCapRadius = rTip * max(0, sin(alphaCone));

        collapseTipI = (forcePointyTipI || (((rTipSeam <= tipCollapseTol) && ~minimalSeam) || ((tipCapHeight < minSpacing) && (tipCapRadius < minSpacing))));
        if forceCylinder
            collapseTipI = false;
        end
        if collapseTipI
            rTipSeam = 0;
        end

        minBaseRadiusMetric = min(minBaseRadiusMetric, rBase);
        minSeamRadiusMetric = min(minSeamRadiusMetric, rTipSeam);
        tipCollapsedAll = tipCollapsedAll && collapseTipI;

        tipCenter = zTipCenter * orientation;

        producedSpikeCount = producedSpikeCount + 1;

        coneHeight = max(1e-6, zTipSeam - zBase);
        baseCountSeq = max(3, baseSegmentsEffective);
        if collapseTipI
            seamSegmentsEffective = max(3, baseCountSeq);
        elseif minimalSeam
            seamSegmentsEffective = 3;
        else
            seamSegmentsEffective = ring_segment_count(rTipSeam, minSpacing, azimuthSpacingFactor);
            seamSegmentsEffective = max(3, seamSegmentsEffective);
        end
        axialSteps = max(1, ceil(coneHeight / (minSpacing * axialSpacingFactor)));
        numConeSteps = axialSteps;

        prevIdx = seamIndices{i};
        seamIsPoint = numel(prevIdx) == 1;

        for k = 1:numConeSteps
            t = k / numConeSteps;              % 0 at base, 1 at seam
            ringRadius = (1 - t) * rBase + t * rTipSeam;
            ringCenter = ((1 - t) * zBase + t * zTipSeam) * orientation;
            if ((forcePointyTipI && k == numConeSteps) || (ringRadius <= tipCollapseTol && k == numConeSteps))
                [V, idxPoint] = append_vertices(V, ringCenter);
                if numel(prevIdx) > 1
                    F = [F; connect_ring_to_apex(prevIdx, idxPoint)]; %#ok<AGROW>
                elseif numel(prevIdx) == 1
                    % nothing to connect; two points on axis
                end
                prevIdx = idxPoint;
                seamIsPoint = true;
                break;
            end

            if k == numConeSteps
                segmentsCurr = seamSegmentsEffective;
            else
                segmentsCurr = ring_segment_count(ringRadius, minSpacing, azimuthSpacingFactor);
            end
            segmentsCurr = max(3, segmentsCurr);

            anglesCurr = circle_angles(segmentsCurr);
            ringPts = seam_ring_points(ringCenter, u, v, ringRadius, anglesCurr);
            [V, idxRing] = append_vertices(V, ringPts);
            if includeCore
                coreMask(idxRing) = false;
            end

            if numel(prevIdx) == 1
                F = [F; connect_ring_to_apex(idxRing, prevIdx)]; %#ok<AGROW>
            else
                facesBridge = bridge_rings_idx(prevIdx(:)', idxRing(:)', V, orientation);
                F = [F; facesBridge]; %#ok<AGROW>
            end

            prevIdx = idxRing;
        end

        seamRingIdx = prevIdx;

        % --- Patch B: Spherical tip cap (ring-based) ---
        apex = tipCenter + rTip * orientation;
        if seamIsPoint
            % Seam already collapsed to a point; cone-to-tip faces created in loop.
        elseif collapseTipI
            [V, idxApex] = append_vertices(V, apex);
            if includeCore
                coreMask(idxApex) = false;
            end
            if numel(seamRingIdx) > 1
                F = [F; connect_ring_to_apex(seamRingIdx, idxApex)]; %#ok<AGROW>
            end
        else
            capArcLength = rTip * max(alphaCone, 1e-6);
            if minimalSeam
                nCapRings = 0;
            else
                nCapRings = max(1, ceil(capArcLength / (minSpacing * capSpacingFactor)));
            end
            prevCapIdx = seamRingIdx;
            for j = 1:nCapRings
                tCap = (nCapRings - j + 1) / (nCapRings + 1);
                psi = alphaCone * tCap;
                ringRadius = rTip * sin(psi);
                ringCenter = tipCenter + (rTip * cos(psi)) * orientation;
                segmentsCap = ring_segment_count(ringRadius, minSpacing, capSpacingFactor);
                segmentsCap = max(3, segmentsCap);
                anglesCap = circle_angles(segmentsCap);
                ringPts = seam_ring_points(ringCenter, u, v, ringRadius, anglesCap);
                [V, idxRing] = append_vertices(V, ringPts);
                if includeCore
                    coreMask(idxRing) = false;
                end
                facesBridge = bridge_rings_idx(prevCapIdx(:)', idxRing(:)', V, orientation);
                F = [F; facesBridge]; %#ok<AGROW>
                prevCapIdx = idxRing;
            end

            [V, idxApex] = append_vertices(V, apex);
            if includeCore
                coreMask(idxApex) = false;
            end
            F = [F; connect_ring_to_apex(prevCapIdx, idxApex)]; %#ok<AGROW>
        end

        % --- Patch C: Toroidal fillet (suppressed for now) ---
        % Intentionally disabled to focus on core sphere + spike body + tip
    end

    if producedSpikeCount == 0
        tipCollapsedAll = false;
    end

    if isinf(minBaseRadiusMetric)
        minBaseRadiusMetric = 0;
    end
    if isinf(minSeamRadiusMetric)
        minSeamRadiusMetric = 0;
    end
    spikeSeamRadius = minSeamRadiusMetric;
    spikeBaseRadiusMetric = minBaseRadiusMetric;
    spikeTipCollapsedMetric = tipCollapsedAll;
    if tipCollapsedAll
        spikeTipEffectiveDiameter = 0;
        spikeTipSphereRadius = 0;
    end

    % After all spikes processed, append remaining active core faces
    F = [coreFCurrent; F];
    
    % Final weld then construct surface mesh
    fprintf('Final welding and surface mesh construction...\n');
    [Vw, ~, ic] = uniquetol(V, 1e-9, 'ByRows', true);
    Fw = ic(F);
    meshSurface = surfaceMesh(Vw, Fw);
    totalVolume = computeMeshVolume(Vw, Fw);

    % Wrap the surface mesh so we can attach auxiliary data (volume, bounds, etc.)
    fprintf('Wrapping up urchin structure...\n');
    urchin = struct();
    urchin.SurfaceMesh = meshSurface;
    urchin.Vertices = meshSurface.Vertices;
    urchin.Faces = meshSurface.Faces;

    urchin.Parameters = struct( ...
        'rCore', rCore, ...
        'coreRadius', rCore, ...
        'spikeLength', spikeLength, ...
        'spikeCount', spikeCount, ...
        'spikeTipRequested', spikeTipRequested, ...
        'spikeTip', spikeTipEffectiveDiameter, ...
        'spikeConicality', spikeConicality, ...
        'filletRatio', filletRatio, ...
        'resolution', resolution, ...
        'useFillet', useFillet, ...
        'includeCore', includeCore, ...
        'flucFactor', flucFactor, ...
        'flucMethod', flucMethod, ...
        'distMethod', distMethod, ...
        'refinedOrientation', refinedOrientation, ...
        'refinedOrientationThreshold', refinedOrientationThresholdDeg, ...
        'genVolume', genVolume, ...
        'volRes', volRes, ...
        'volPadding', volPadding, ...
        'volAlpha', volAlphaInp, ...
        'volAdaptive', volAdaptive, ...
        'volDxMax', volDxMax, ...
        'volDxMin', volDxMin, ...
        'volBlockSize', volBlockSz, ...
        'volCriterion', volCriterion, ...
        'spikeOrientations', spikeOrientations, ...
        'spikeLengths', spikeLengths, ...
        'spikeBaseRadii', rBases, ...
        'spikeTipRadii', rTips ...
    );

    if isempty(rBaseMaxs)
        spikeBaseMaximaVector = rBaseMaxs;
    elseif size(rBaseMaxs, 2) > 1
        spikeBaseMaximaVector = rBaseMaxs(:, 1);
    else
        spikeBaseMaximaVector = rBaseMaxs;
    end
    spikeBaseMaximaVector = spikeBaseMaximaVector(:);

    urchin.Metrics = struct( ...
        'MinimumSpacing', minSpacing, ...
        'SpikeBaseRadius', spikeBaseRadiusMetric, ...
        'SpikeBaseMaxima', spikeBaseMaximaVector, ...
        'SpikeSeamRadius', spikeSeamRadius, ...
        'SpikeTipRadius', spikeTipSphereRadius, ...
        'SpikeTipDiameter', spikeTipEffectiveDiameter, ...
        'SpikeTipCollapsed', spikeTipCollapsedMetric, ...
        'TotalVolume', totalVolume ...
    );


    %% 5) Build a volumetric mask from the surface mesh by voxelizing
    % the shape using either inpolyhedron (if available) or alphaShape fallback.
    if genVolume && ~volAdaptive
        fprintf('Generating dense volumetric mask via voxelization...\n');
        Vmin = min(urchin.Vertices,[],1);
        Vmax = max(urchin.Vertices,[],1);
        span = Vmax - Vmin;
        Vmin = Vmin - volPadding * span;
        Vmax = Vmax + volPadding * span;
        % Build cubic voxels sized to largest dimension / volRes
        maxDim = max(span);
        dx = maxDim / volRes;
        xs = Vmin(1):dx:Vmax(1);
        ys = Vmin(2):dx:Vmax(2);
        zs = Vmin(3):dx:Vmax(3);
        [XX,YY,ZZ] = ndgrid(xs,ys,zs);
        P = [XX(:), YY(:), ZZ(:)];
        inside = false(size(P,1),1);
        insideHasValue = false;
        % Try inpolyhedron (File Exchange) if available
        try
            inside = inpolyhedron(urchin.Faces, urchin.Vertices, P);
            insideHasValue = true;
        catch
            % Fallback: alphaShape from mesh vertices, choose alpha automatically or use input
            warning('inpolyhedron failed; falling back to alphaShape.');
            if isempty(volAlphaInp)
                % heuristic alpha: small fraction of bounding radius
                Rb = 0.5 * maxDim;
                alphaVol = 0.08 * Rb;
            else
                alphaVol = volAlphaInp;
            end
            try
                shpVol = alphaShape(urchin.Vertices, alphaVol);
                inside = inShape(shpVol, P(:,1), P(:,2), P(:,3));
                insideHasValue = true;
            catch
                % As a last resort, mark nothing inside
                inside = false(size(P,1),1);
                insideHasValue = true;
            end
        end
        if ~insideHasValue
            inside = false(size(P,1),1);
        end
        VolumeMask = reshape(logical(inside), size(XX));
        urchin.VolumeMask = VolumeMask;
        urchin.VoxelGrid.X = xs;
        urchin.VoxelGrid.Y = ys;
        urchin.VoxelGrid.Z = zs;
        urchin.VoxelSize = dx;
        urchin.Bounds = [Vmin; Vmax];
    elseif genVolume && volAdaptive
        % Adaptive sparse octree voxelization (VDB-style)
        fprintf('Generating adaptive sparse octree volumetric representation...\n');
        Vmin = min(urchin.Vertices,[],1);
        Vmax = max(urchin.Vertices,[],1);
        span = Vmax - Vmin;
        Vmin = Vmin - volPadding * span;
        Vmax = Vmax + volPadding * span;
        maxDim = max(span);
        if isempty(volDxMax), volDxMax = maxDim / 128; end
        if isempty(volDxMin), volDxMin = volDxMax / 4; end
    insideFn = make_inside_tester(urchin, volAlphaInp);
        leaves = build_adaptive_octree(Vmin, Vmax, volDxMax, volDxMin, volBlockSz, insideFn, volCriterion);
        urchin.VolumeOctree.Leaves = leaves;
        urchin.VolumeOctree.BlockSize = volBlockSz;
        urchin.VolumeOctree.DxMin = volDxMin;
        urchin.VolumeOctree.DxMax = volDxMax;
        urchin.VolumeOctree.Bounds = [Vmin; Vmax];
    end

    elapsedTime = toc;
    fprintf('Urchin created successfully in %.2f seconds.\n', elapsedTime);
    fprintf('Generated mesh with %d vertices and %d faces.\n', size(urchin.Vertices, 1), size(urchin.Faces, 1));

    urchin.Metrics.GenerationTimeSeconds = elapsedTime;
    urchin.Metrics.VertexCount = size(urchin.Vertices, 1);
    urchin.Metrics.FaceCount = size(urchin.Faces, 1);

    %% 6) Visualize if no output is requested
    if nargout == 0
        viewer = viewer3d;
        viewer.CameraPosition = [0 0 -scale];
        viewer.CameraUpVector = [0 1 0];
        viewer.Lighting = 'off';
        surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, Alpha=0.5);
        surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, Wireframe=true);
        % surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, VerticesOnly=true);

    end
end

