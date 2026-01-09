function urchin = urchin(options)
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
%   - Spike fluctuations: per-spike length jitter (FlucFactor) with uniform/random (Sobol or seeded RNG).
%   - Orientation distributions: Fibonacci-like uniform or random (seeded).
%   - Resolution-driven sampling: a single dimensionless knob controls
%     all curvature-aware feature counts while keeping seams watertight.
%   - Volume representations:
%       • Dense uniform voxel grid via inpolyhedron/alphaShape fallback
%       • Adaptive sparse octree (VDB-style) with user control of voxel
%         size range (volDxMin/Max), leaf block size, and refinement
%         criterion (boundary/distance/curvature/hybrid). Stored in
%         mesh.VolumeOctree for sparse workflows; optional densification.
%
% Usage:
%   urchinStruct = urchin(CoreRadius=30, SpikeLength=15, SpikeCount=50, SpikeTipDiameter=5);
%   urchinStruct = urchin(CoreRadius=25, SpikeLength=10, SpikeCount=75);
%   diagnostics = meshDiagnostics(urchinStruct.SurfaceMesh);
%   fprintf('Watertight: %d\n', diagnostics.IsWatertight);
%
% Inputs (name-value):
%   CoreRadius   Core radius (nm). Default 1
%   SpikeLength  Spike length from core surface (nm). Default 1
%   SpikeCount   Number of spikes. Default 100
%   SpikeTipDiameter Diameter of the spherical tip cap (nm).
%                Defaults to SpikeLength/10 when omitted.
%   SpikeConicality Conicality [-1..1]. -1 collapses base, 0=cylinder, 1=widest base.
%   Resolution   Dimensionless mesh resolution (default 100). Smallest feature
%                length is ~2*(CoreRadius+SpikeLength)/Resolution.
%   UseFillet    Enable toroidal fillet patch. Default true (fillet radius
%                defaults to SpikeTipDiameter/2 and is clamped relative to base radius)
%   FlucFactor   Spike length fluctuation factor [0..1]. Default 0.5
%   FlucMethod   "uniform" (Sobol) or "random" (seeded by SpikeCount). Default "uniform"
%   DistMethod   "uniform" (Fibonacci-like) or "random" (seeded by SpikeCount). Default "uniform"
%   RefinedOrientation Toggle iterative Coulomb-like relaxation (default true)
%   RefinedOrientationThresholdDeg Minimum allowable angular separation in degrees (default 0.1)
%
% Volume (dense):
%   GenVolume    If true, generate dense voxel mask. Default false
%   VolResolution Target voxels along largest dimension. Default 128
%   VolPadding   Fractional AABB padding for voxelization. Default 0.05
%   VolAlpha     Alpha for alphaShape fallback (auto if empty)
%
% Volume (adaptive sparse / VDB-style):
%   VolAdaptive  If true, build sparse octree volume instead of dense mask
%   VolDxMax     Max voxel size (auto: maxDim/128)
%   VolDxMin     Min voxel size (auto: VolDxMax/4)
%   VolBlockSize Leaf block resolution (voxel dimension), e.g., 8 or 16
%   VolCriterion "boundary" | "distance" | "curvature" | "hybrid" (boundary by default)
%
% Output:
%   mesh struct with fields: Vertices [n×3], Faces [m×3]
%     + Parameters (struct echoing all name-value inputs; collapsed spike tips report 0)
%     + Metrics (struct with min spacing, spike base radius, enclosed volume, etc.)
%     + VolumeMask (logical 3D) and VoxelGrid axes if GenVolume=true & ~VolAdaptive
%     + VolumeOctree (sparse leaves) if GenVolume=true & VolAdaptive=true
%
% Dependencies: 
%   sobolset requires Statistics and Machine Learning ToolboxR
%   gpuDeviceCount requires Parallel Computing Toolbox
%   wrapToPi requires Mapping Toolbox
%   surfaceMesh requires Lidar Toolbox
%
% =========================================================================

    arguments
        options.CoreRadius (1,1) double {mustBePositive} = 1
        options.SpikeLength (1,1) double {mustBePositive} = 1
        options.SpikeCount (1,1) double {mustBeInteger, mustBeNonnegative} = 100
        options.SpikeTipDiameter (1,:) double {mustBeNonnegative} = []
        options.SpikeConicality (1,1) double {mustBeInRange(options.SpikeConicality, -1, 1)} = 0.5
        options.FilletRatio (1,1) double {mustBeInRange(options.FilletRatio, 0, 0.5)} = 0.25
        options.Resolution (1,1) double {mustBePositive} = 100
        options.UseFillet (1,1) logical = true
        options.IncludeCore (1,1) logical = true
        options.FlucFactor (1,1) double {mustBeInRange(options.FlucFactor, 0, 1)} = 0.5
        options.FlucMethod (1,1) string {mustBeMember(options.FlucMethod, ["uniform", "random", "gaussian"])} = "uniform"
        options.DistMethod (1,1) string {mustBeMember(options.DistMethod, ["uniform", "random"])} = "uniform"
        options.RefinedOrientation (1,1) logical = true
        options.RefinedOrientationThresholdDeg (1,1) double {mustBeInRange(options.RefinedOrientationThresholdDeg, 0, 180)} = 0.1
        options.GenVolume (1,1) logical = false
        options.VolResolution (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.VolResolution, 16)} = 128
        options.VolPadding (1,1) double {mustBeInRange(options.VolPadding, 0, 0.5)} = 0.05
        options.VolAlpha (1,:) double {mustBePositive} = []
        options.VolAdaptive (1,1) logical = false
        options.VolDxMax (1,:) double {mustBePositive} = []
        options.VolDxMin (1,:) double {mustBePositive} = []
        options.VolBlockSize (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.VolBlockSize, 4)} = 8
        options.VolCriterion (1,1) string {mustBeMember(options.VolCriterion, ["boundary", "distance", "curvature", "hybrid"])} = "boundary"
    end

    fprintf("Starting B-Rep Urchin Generation...\n");
    tic();

    %% 2) Initial Geometric Calculations

    % Extract parameters (using new arguments structure)
    rCore = options.CoreRadius;
    spikeLength = options.SpikeLength;
    spikeCount = options.SpikeCount;
    tipDiameter = options.SpikeTipDiameter;
    spikeTipRequested = tipDiameter;
    spikeConicality = options.SpikeConicality;
    filletRatio = options.FilletRatio;
    resolution = options.Resolution;
    useFillet = options.UseFillet;
    includeCore = options.IncludeCore;
    flucFactor = options.FlucFactor;
    flucMethod = char(options.FlucMethod);  % Convert string to char for switch statement
    distMethod = char(options.DistMethod);  % Convert string to char for switch statement
    refinedOrientation = options.RefinedOrientation;
    refinedOrientationThresholdDeg = options.RefinedOrientationThresholdDeg;
    genVolume = options.GenVolume;
    volRes = options.VolResolution;
    volPadding = options.VolPadding;
    volAlphaInp = options.VolAlpha;
    volAdaptive = options.VolAdaptive;
    volDxMax = options.VolDxMax;
    volDxMin = options.VolDxMin;
    volBlockSz = options.VolBlockSize;
    volCriterion = char(options.VolCriterion);  % Convert string to char for switch statement

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
        case "uniform"
            % Golden spiral distribution across sphere via (theta,phi)
            Phi = (1 + sqrt(5)) / 2;  % golden ratio
            i = (1:spikeCount)';
            theta = acos(1 - 2 .* i ./ (spikeCount + 1));
            phi   = mod(2 * pi * i / Phi, 2 * pi);
        case "random"
            rng(spikeCount,'twister');
            phi = 2 * pi * rand(spikeCount, 1);
            cosTheta = 2 * rand(spikeCount, 1) - 1;
            theta = acos(cosTheta);
    end

    spikeOrientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)]; % #KeyParameter

    if refinedOrientation && spikeCount > 1
        spikeOrientations = refineSpikeOrientations(spikeOrientations, refinedOrientationThresholdDeg);
    end

    % Precompute orthonormal frames for each spike orientation to avoid repeated planeVectors calls
    uFrames = zeros(spikeCount, 3);
    vFrames = zeros(spikeCount, 3);
    for idx = 1:spikeCount
        [uVec, vVec] = planeVectors(spikeOrientations(idx, :));
        uFrames(idx, :) = uVec(:).';
        vFrames(idx, :) = vVec(:).';
    end
    %% Spike length fluctuations
    switch flucMethod
        case "uniform"
            flucs = net(sobolset(1), spikeCount);   % values in ~[0,1]
        case "random"
            rng(spikeCount,'twister');
            flucs = rand(spikeCount,1);            % values in ~[0,1]
        case "gaussian"
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
        [coreVFull, coreFFull] = triangulateIcosphere(rCore, coreSubdiv);
        [V, coreIdx] = appendVertices(V, coreVFull); % core vertices global indices
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
    fprintf("Generating & stitching %d spikes...\n", spikeCount);
    minBaseRadiusMetric = Inf;
    minSeamRadiusMetric = Inf;
    tipCollapsedAll = true;
    producedSpikeCount = 0;
    seamIndices = cell(spikeCount,1);
    
    % PRE-ALLOCATE face array to avoid O(n²) dynamic reallocation (50-70% speedup)
    % Estimate: ~100-200 faces per spike (base ring + cone + tip)
    estimatedFacesPerSpike = 150;
    maxEstimatedFaces = size(coreFCurrent, 1) + spikeCount * estimatedFacesPerSpike;
    F = zeros(maxEstimatedFaces, 3);
    faceCount = size(coreFCurrent, 1);
    if faceCount > 0
        F(1:faceCount, :) = coreFCurrent;
    end

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
        % Local orthonormal frame (precomputed):
        %  - orientation: spike axis (unit vector)
        %  - u, v: span orthogonal plane (used for circular rings)
        u = uFrames(i, :);
        v = vFrames(i, :);

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
            % Extract core trimming into dedicated function
            [coreLoop, coreFCurrent, coreMask] = trimCoreForSpike(V, coreFCurrent, coreMask, ...
                                                                     orientation, cosCutoffI, rBase, ...
                                                                     seamIndices, i);
        end

        seamBaseCenter = zBase * orientation;
        baseSegmentsEffective = 3;
        if collapseBaseI
            [V, seamIndices{i}] = appendVertices(V, seamBaseCenter);
        else
            baseSegmentsI = ringSegmentCount(rBase, minSpacing, azimuthSpacingFactor);
            if forceCylinder
                baseSegmentsI = max(baseSegmentsI, 3);
            end
            baseSegmentsEffective = max(3, baseSegmentsI);
            seamBaseRing = generateCircularRing(seamBaseCenter, u, v, rBase, minSpacing, azimuthSpacingFactor);
            [V, seamIndices{i}] = appendVertices(V, seamBaseRing);
            if includeCore
                coreMask(seamIndices{i}) = true;
            end

            if includeCore && numel(coreLoop) >= 3
                seamRing = seamIndices{i}(:)';
                newBridgeFaces = bridgeRingsIdx(coreLoop, seamRing, V, orientation);
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
            [V, baseCenterIdx] = appendVertices(V, baseCenter);
            baseFaces = connectRingToApex(seamIndices{i}, baseCenterIdx);
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
            % Append to pre-allocated array with range assignment (no dynamic growth)
            nNewFaces = size(baseFaces, 1);
            if faceCount + nNewFaces > size(F, 1)
                F = [F; zeros(ceil(size(F,1)*0.5), 3)];
            end
            F(faceCount+1:faceCount+nNewFaces, :) = baseFaces;
            faceCount = faceCount + nNewFaces;
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
            seamSegmentsEffective = ringSegmentCount(rTipSeam, minSpacing, azimuthSpacingFactor);
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
                [V, idxPoint] = appendVertices(V, ringCenter);
                if numel(prevIdx) > 1
                    newFaces = connectRingToApex(prevIdx, idxPoint);
                    nNewFaces = size(newFaces, 1);
                    if faceCount + nNewFaces > size(F, 1)
                        F = [F; zeros(ceil(size(F,1)*0.5), 3)];
                    end
                    if nNewFaces > 0
                        F(faceCount+1:faceCount+nNewFaces, :) = newFaces;
                        faceCount = faceCount + nNewFaces;
                    end
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
                segmentsCurr = ringSegmentCount(ringRadius, minSpacing, azimuthSpacingFactor);
            end
            segmentsCurr = max(3, segmentsCurr);

            ringPts = generateCircularRing(ringCenter, u, v, ringRadius, minSpacing, azimuthSpacingFactor);
            [V, idxRing] = appendVertices(V, ringPts);
            if includeCore
                coreMask(idxRing) = false;
            end

            if numel(prevIdx) == 1
                newFaces = connectRingToApex(idxRing, prevIdx);
            else
                newFaces = bridgeRingsIdx(prevIdx(:)', idxRing(:)', V, orientation);
            end
            nNewFaces = size(newFaces, 1);
            if faceCount + nNewFaces > size(F, 1)
                F = [F; zeros(ceil(size(F,1)*0.5), 3)];
            end
            if nNewFaces > 0
                F(faceCount+1:faceCount+nNewFaces, :) = newFaces;
                faceCount = faceCount + nNewFaces;
            end

            prevIdx = idxRing;
        end

        seamRingIdx = prevIdx;

        % --- Patch B: Spherical tip cap (ring-based) ---
        apex = tipCenter + rTip * orientation;
        if seamIsPoint
            % Seam already collapsed to a point; cone-to-tip faces created in loop.
        elseif collapseTipI
            [V, idxApex] = appendVertices(V, apex);
            if includeCore
                coreMask(idxApex) = false;
            end
            if numel(seamRingIdx) > 1
                F = [F; connectRingToApex(seamRingIdx, idxApex)]; %#ok<AGROW>
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
                ringPts = generateCircularRing(ringCenter, u, v, ringRadius, minSpacing, capSpacingFactor);
                [V, idxRing] = appendVertices(V, ringPts);
                if includeCore
                    coreMask(idxRing) = false;
                end
                newFaces = bridgeRingsIdx(prevCapIdx(:)', idxRing(:)', V, orientation);
                nNewFaces = size(newFaces, 1);
                if faceCount + nNewFaces > size(F, 1)
                    F = [F; zeros(ceil(size(F,1)*0.5), 3)];
                end
                if nNewFaces > 0
                    F(faceCount+1:faceCount+nNewFaces, :) = newFaces;
                    faceCount = faceCount + nNewFaces;
                end
                prevCapIdx = idxRing;
            end

            [V, idxApex] = appendVertices(V, apex);
            if includeCore
                coreMask(idxApex) = false;
            end
            newFaces = connectRingToApex(prevCapIdx, idxApex);
            nNewFaces = size(newFaces, 1);
            if faceCount + nNewFaces > size(F, 1)
                F = [F; zeros(ceil(size(F,1)*0.5), 3)];
            end
            if nNewFaces > 0
                F(faceCount+1:faceCount+nNewFaces, :) = newFaces;
                faceCount = faceCount + nNewFaces;
            end
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

    % Trim pre-allocated face array to actual size
    F = F(1:faceCount, :);
    
    % Final weld then construct surface mesh
    fprintf("Final welding and surface mesh construction...\n");
    % Coarse binning to reduce duplicates before uniquetol (faster)
    coarseTol = 1e-9;
    Vcoarse = round(V / coarseTol) * coarseTol;
    [Vpre, ~, icPre] = unique(Vcoarse, 'rows', 'stable');
    Fpre = icPre(F);

    try
        [Vw, ~, ic] = uniquetol(Vpre, 1e-9, 'ByRows', true);
        Fw = ic(Fpre);
        meshSurface = surfaceMesh(Vw, Fw);
    catch ME
        fprintf("Error during final mesh construction: %s\n", ME.message);
        % Fallback: less aggressive tolerance
        [Vw, ~, ic] = uniquetol(Vpre, 1e-6, 'ByRows', true);
        Fw = ic(Fpre);
        meshSurface = surfaceMesh(Vw, Fw);
    end
    totalVolume = computeMeshVolume(Vw, Fw);

    % Wrap the surface mesh so we can attach auxiliary data (volume, bounds, etc.)
    fprintf("Wrapping up urchin structure...\n");
    urchin = struct();
    urchin.SurfaceMesh = meshSurface;
    urchin.Vertices = meshSurface.Vertices;
    urchin.Faces = meshSurface.Faces;

    urchin.Parameters = struct( ...
        "CoreRadius", rCore, ...
        "SpikeLength", spikeLength, ...
        "SpikeCount", spikeCount, ...
        "SpikeTipRequested", spikeTipRequested, ...
        "SpikeTipDiameter", spikeTipEffectiveDiameter, ...
        "SpikeConicality", spikeConicality, ...
        "FilletRatio", filletRatio, ...
        "Resolution", resolution, ...
        "UseFillet", useFillet, ...
        "IncludeCore", includeCore, ...
        "FlucFactor", flucFactor, ...
        "FlucMethod", flucMethod, ...
        "DistMethod", distMethod, ...
        "RefinedOrientation", refinedOrientation, ...
        "RefinedOrientationThresholdDeg", refinedOrientationThresholdDeg, ...
        "GenVolume", genVolume, ...
        "VolResolution", volRes, ...
        "VolPadding", volPadding, ...
        "VolAlpha", volAlphaInp, ...
        "VolAdaptive", volAdaptive, ...
        "VolDxMax", volDxMax, ...
        "VolDxMin", volDxMin, ...
        "VolBlockSize", volBlockSz, ...
        "VolCriterion", volCriterion, ...
        "SpikeOrientations", spikeOrientations, ...
        "SpikeLengths", spikeLengths, ...
        "SpikeBaseRadii", rBases, ...
        "SpikeTipRadii", rTips ...
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
        "MinimumSpacing", minSpacing, ...
        "SpikeBaseRadius", spikeBaseRadiusMetric, ...
        "SpikeBaseMaxima", spikeBaseMaximaVector, ...
        "SpikeSeamRadius", spikeSeamRadius, ...
        "SpikeTipRadius", spikeTipSphereRadius, ...
        "SpikeTipDiameter", spikeTipEffectiveDiameter, ...
        "SpikeTipCollapsed", spikeTipCollapsedMetric, ...
        "TotalVolume", totalVolume ...
    );


    %% 5) Build a volumetric mask from the surface mesh by voxelizing
    urchin = generateVolumeRepresentation(urchin, genVolume, volAdaptive, volPadding, volRes, ...
                                          volAlphaInp, volDxMax, volDxMin, volBlockSz, volCriterion);

    elapsedTime = toc();
    fprintf("Urchin created successfully in %.2f seconds.\n", elapsedTime);
    fprintf("Generated mesh with %d vertices and %d faces.\n", size(urchin.Vertices, 1), size(urchin.Faces, 1));

    urchin.Metrics.GenerationTimeSeconds = elapsedTime;
    urchin.Metrics.VertexCount = size(urchin.Vertices, 1);
    urchin.Metrics.FaceCount = size(urchin.Faces, 1);

    %% 6) Visualize if no output is requested
    if nargout() == 0
        viewer = viewer3d();
        viewer.CameraPosition = [0 0 -scale];
        viewer.CameraUpVector = [0 1 0];
        viewer.Lighting = "off";
        surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, Alpha=0.5);
        surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, Wireframe=true);
        % surfaceMeshShow(urchin.SurfaceMesh, Parent=viewer, VerticesOnly=true);

    end
end

