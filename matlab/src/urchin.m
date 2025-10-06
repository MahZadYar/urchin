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
%   urchinStruct = urchin('coreRadius',30,'spikeLength',15,'spikeCount',50,'spikeTip',5);
%   urchinStruct = urchin('coreRadius',25,'spikeLength',10,'spikeCount',75);
%   diagnostics = meshDiagnostics(urchinStruct.SurfaceMesh);
%   fprintf('Watertight: %d\n', diagnostics.IsWatertight);
%
% Inputs (name-value):
%   coreRadius   Core radius (nm). Default 1
%   spikeLength  Spike length from core surface (nm). Default 1
%   spikeCount   Number of spikes. Default 100
%   spikeTip     Diameter of the spherical tip cap (nm).
%                Defaults to spikeLength/10 when omitted.
%   spikeConicality Conicality [-1..1]. -1 collapses base, 0=cylinder, 1=widest base.
%   resolution   Dimensionless mesh resolution (default 100). Smallest feature
%                length is ~2*(coreRadius+spikeLength)/resolution.
%   useFillet    Enable toroidal fillet patch. Default true (fillet radius
%                defaults to spikeTip/2 and is clamped relative to base radius)
%   flucFactor   Spike length fluctuation factor [0..1]. Default 0.5
%   flucMethod   'uniform' (Sobol) or 'random' (seeded by spikeCount). Default 'uniform'
%   distMethod   'uniform' (Fibonacci-like) or 'random' (seeded by spikeCount). Default 'uniform'
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
    addParameter(p, 'coreRadius', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
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
    coreRadius = p.Results.coreRadius;
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

    %% 2) Initial Geometric Calculations & Spike Distribution

    scale = 2 * (coreRadius + spikeLength); % urchin scale for normalization
    min_spacing = max(1e-9,  scale/resolution);
    core_spacing_factor = 2.0;
    azimuth_spacing_factor = 1.2;
    axial_spacing_factor = 1.0;
    cap_spacing_factor = 1.0;
    collapse_st_tol = min_spacing;
    conicalityFlatThreshold = 0.05;
    conicalityEffective = spikeConicality;
    force_cylinder = false;
    pending_pointy_tip = false;
    % Tip sphere diameter (default: spikeLength/10 unless specified)
    if isempty(tipDiameter)
        tipDiameter = spikeLength / 10;
    end
    tipDiameter = max(tipDiameter, 0);
    tipSphereRadius = tipDiameter / 2;

    if tipDiameter <= collapse_st_tol
        tipDiameter = collapse_st_tol;
        tipSphereRadius = tipDiameter / 2;
        if conicalityEffective < 0
            % negative conicality handled later via tapered base rule
        elseif conicalityEffective <= conicalityFlatThreshold
            conicalityEffective = 0;
            force_cylinder = true;
        else
            pending_pointy_tip = true;
        end
    end


    % Note: Seam axial location and cone tip are computed per-spike (fluctuations)

    % Spike height unused explicitly; seam defined at z=coreRadius+spikeLength

    % Distribute spike orientations (uniform Fibonacci or random)
    switch distMethod
        case 'uniform'
            % Golden spiral distribution across sphere via (theta,phi)
            Phi = (1 + sqrt(5)) / 2;  % golden ratio
            i = (1:spikeCount)';
            theta = acos(1 - 2 .* i ./ (spikeCount + 1));
            phi   = mod(2 * pi * i / Phi, 2 * pi);
            spike_orientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
        case 'random'
            rng(spikeCount,'twister');
            phi = 2 * pi * rand(spikeCount, 1);
            cosTheta = 2 * rand(spikeCount, 1) - 1;
            theta = acos(cosTheta);
            spike_orientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    end

    % Spike length fluctuations
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

    % =====================================================================
    % New maximum base radius rationale (non-overlapping spherical footprints)
    % ---------------------------------------------------------------------
    % Previous approach: equal-area spherical caps -> could yield overlapping
    % circular footprints on the core sphere because sum of max circle areas
    % (π r_base^2 each) is strictly less than total sphere area, so area-based
    % allocation is not a tight packing constraint.
    % New approach: find the smallest angular separation θ_min between any two
    % spike orientation vectors (after distribution). Circular footprints on
    % the core are spherical discs with angular (geodesic) radius α. Two discs
    % are just tangent when 2α = θ_sep. To avoid overlap for all pairs we set
    % α_max = θ_min / 2. The Euclidean (chord) radius of that footprint circle
    % on a sphere of radius coreRadius is:
    %   r_base_max = coreRadius * sin(α_max)
    % This guarantees non-overlapping base circles (they may be tangent for at
    % least one closest pair) while being as large as analytically possible
    % under the given orientation set.
    if spikeCount == 0
        theta_min_per_spike = zeros(0,1);
    elseif spikeCount == 1
        theta_min_per_spike = pi;
    else
        % Compute maximum cosine (excluding self) to find minimal angular sep per spike
        Ddot = spike_orientations * spike_orientations.'; % spikeCount x spikeCount dot products
        Ddot(1:spikeCount+1:end) = -Inf; % ignore diagonal
        maxDotPerSpike = max(Ddot, [], 2);
        maxDotPerSpike(~isfinite(maxDotPerSpike)) = -1;
        maxDotPerSpike = min(1, max(-1, maxDotPerSpike));
        theta_min_per_spike = acos(maxDotPerSpike); % smallest angle between each spike and its nearest neighbour
    end

    if isempty(theta_min_per_spike)
        theta_min_global = pi;
    else
        theta_min_global = min(theta_min_per_spike);
    end

    alpha_base_max_per_spike = 0.5 * theta_min_per_spike;           % tangent half-angle per spike

    theta_ratio = coreRadius / (coreRadius + spikeLength + tipSphereRadius);
    theta_ratio = min(1, max(-1, theta_ratio));
    theta_max = acos(theta_ratio);
    alpha_base_max_per_spike = min(alpha_base_max_per_spike, theta_max);      % clamp to tangent angle per spike

    r_base_max_per_spike = coreRadius * sin(alpha_base_max_per_spike);      % maximal non-overlapping base radius per spike
    if ~isempty(r_base_max_per_spike)
        r_base_max_per_spike = min(r_base_max_per_spike, 0.999 * coreRadius);
        r_base_max_per_spike = max(r_base_max_per_spike, 1e-9);
    end

    r_base_max = coreRadius * sin(0.5 * theta_min_global);
    r_base_max = min(r_base_max, 0.999 * coreRadius);
    r_base_max = max(r_base_max, 1e-9);


    base_collapse_tol = max(0.5 * min_spacing, 1e-9);
    if force_cylinder
        base_collapse_tol = min(base_collapse_tol, 0.25 * min_spacing);
    end

    tip_collapse_tol = max(0.5 * min_spacing, 1e-9);
    if force_cylinder
        tip_collapse_tol = min(tip_collapse_tol, 0.25 * min_spacing);
    end

    spikeTipEffectiveDiameter = tipDiameter;
    spikeTipSphereRadius = tipSphereRadius;

    % Points container not used; direct append to V

    % Core sphere subdivision target driven purely by resolution
    target_core_spacing = core_spacing_factor * min_spacing;
    target_core_spacing = max(target_core_spacing, 1e-9);
    core_surface_area = 4 * pi * coreRadius^2;
    target_vertices = max(12, ceil(core_surface_area / (target_core_spacing^2)));
    if target_vertices <= 12
        core_subdiv = 0;
    else
        core_subdiv = ceil(0.5 * log2((target_vertices - 2) / 10));
        core_subdiv = max(core_subdiv, 0);
    end

    %% 3) Build full core sphere FIRST (no trimming yet)
    V = zeros(0,3); F = zeros(0,3);
    if includeCore
        core_subdiv = max(floor(core_subdiv), 0);
    [core_V_full, core_F_full] = triangulate_icosphere(coreRadius, core_subdiv);
        [V, core_idx] = append_vertices(V, core_V_full); % core vertices global indices
        core_F_current = core_F_full + core_idx(1) - 1;   % active core faces (global indices)
    core_dirs_full = core_V_full / coreRadius;                % unit directions for initial removal test
        core_mask = false(size(V,1),1);
        core_mask(core_idx) = true;
    else
        core_idx = zeros(0,1);
        core_F_current = zeros(0,3);
        core_dirs_full = zeros(0,3);
        core_mask = false(0,1);
    end

    minBaseRadiusMetric = Inf;
    minSeamRadiusMetric = Inf;
    tipCollapsedAll = true;
    produced_spike_count = 0;

    %% 4) Generate spikes one-by-one WITH per-spike core trimming & stitching
    fprintf('Generating & stitching %d spikes...\n', spikeCount);
    seam_indices = cell(spikeCount,1);
    for i = 1:spikeCount
        orientation = spike_orientations(i, :);
        % Local orthonormal frame:
        %  - orientation: spike axis (unit vector)
        %  - u, v: span orthogonal plane (used for circular rings)
        [u, v] = plane_vectors(orientation);

        core_loop = [];
        seam_indices{i} = [];

        spike_length_i = max(spikeLength + flucs(i), 0);
        tipSphereRadius_i = max(0, tipSphereRadius);
        if tipSphereRadius_i > 0
            tipSphereRadius_i = min(tipSphereRadius_i, spike_length_i);
        end
        z_apex_i = coreRadius + spike_length_i;
        z_tip_center_i = z_apex_i - tipSphereRadius_i;

        if isempty(r_base_max_per_spike)
            r_base_max_i = r_base_max;
        else
            r_base_max_i = r_base_max_per_spike(i);
        end

        [r_base_nominal_i, r_tip_nominal_i, alpha_nominal_default_i, z_seam_nominal_default_i] = compute_nominal_spike_profile( ...
            coreRadius, spike_length_i, tipSphereRadius_i, conicalityEffective, r_base_max_i, scale);

        r_base_i = r_base_nominal_i;
        min_base_radius = 0.5 * min_spacing;
        if tipSphereRadius_i > 0
            tangent_limit = tangent_base_limit(coreRadius, tipSphereRadius_i, spike_length_i);
            if tangent_limit > 1e-9
                r_base_i = min(r_base_i, tangent_limit);
            else
                fallback_base = min(r_base_max_i, max(min_base_radius, 1e-9));
                if fallback_base > 0
                    r_base_i = max(r_base_i, fallback_base);
                end
            end
        end
        r_base_i = min(r_base_i, r_base_max_i);
        r_base_i = max(r_base_i, 0);
        if r_base_i > 0 && r_base_i < min_base_radius
            r_base_i = min(min_base_radius, max(r_base_i, r_base_max_i));
        end
        if r_base_i >= coreRadius
            r_base_i = 0.999 * coreRadius;
        end

        collapse_base_i = r_base_i <= 1e-9;
        if collapse_base_i
            r_base_i = 0;
        end

        cos_cutoff_i = 1;
        if ~collapse_base_i && coreRadius > 0
            ratio = min(1, max(r_base_i / coreRadius, 0));
            cos_cutoff_i = sqrt(max(0, 1 - ratio.^2));
        end

        z_base_i = sqrt(max(1e-12, coreRadius^2 - r_base_i^2));

        if z_base_i < min_spacing
            tipCollapsedAll = false;
            continue;
        end

        if includeCore && ~collapse_base_i
            core_idx = find(core_mask);
            if isempty(core_idx)
                core_dirs_full = zeros(0,3);
            else
                dirs = V(core_idx,:);
                norms = vecnorm(dirs,2,2);
                norms(norms < 1e-12) = 1;
                core_dirs_full = dirs ./ norms;
            end

            core_dot = core_dirs_full * orientation';
            newRemoveMask = core_dot > cos_cutoff_i + 1e-12;
            if ~any(newRemoveMask) && ~isempty(core_idx)
                [~, fallbackPos] = max(core_dot);
                if ~isnan(fallbackPos) && fallbackPos >= 1 && fallbackPos <= numel(newRemoveMask)
                    newRemoveMask(fallbackPos) = true;
                end
            end
            if any(newRemoveMask)
                verts_remove_global = core_idx(newRemoveMask);
                core_mask(verts_remove_global) = false;
                faceRemoveMask = any(ismember(core_F_current, verts_remove_global), 2);
                removedFaces = core_F_current(faceRemoveMask, :);
                core_F_current(faceRemoveMask, :) = [];

                verts_from_removedFaces = unique(removedFaces(:));
                core_loop = setdiff(verts_from_removedFaces, verts_remove_global);

                candidates = {};
                if numel(core_loop) >= 3
                    candidates{end+1} = core_loop(:)'; %#ok<AGROW>
                end
                if ~isempty(verts_remove_global)
                    for j = 1:i-1
                        ring = seam_indices{j};
                        if numel(ring) < 3, continue; end
                        if all(ismember(ring, verts_remove_global))
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
                    core_loop = bestLoop;
                end
            end
        end

        seam_base_center = z_base_i * orientation;
        base_segments_effective = 3;
        if collapse_base_i
            [V, seam_indices{i}] = append_vertices(V, seam_base_center);
        else
            base_segments_i = ring_segment_count(r_base_i, min_spacing, azimuth_spacing_factor);
            if force_cylinder
                base_segments_i = max(base_segments_i, 3);
            end
            base_segments_effective = max(3, base_segments_i);
            angles_base = circle_angles(base_segments_effective);
            seam_base_ring = seam_ring_points(seam_base_center, u, v, r_base_i, angles_base);
            [V, seam_indices{i}] = append_vertices(V, seam_base_ring);
            if includeCore
                core_mask(seam_indices{i}) = true;
            end

            if includeCore && numel(core_loop) >= 3
                seam_ring = seam_indices{i}(:)';
                new_bridge_faces = bridge_rings_idx(core_loop, seam_ring, V, orientation);
                core_F_current = [core_F_current; new_bridge_faces]; %#ok<AGROW>
                core_mask(core_loop) = true;
            end

            if includeCore
                core_idx = find(core_mask);
                if isempty(core_idx)
                    core_dirs_full = zeros(0,3);
                else
                    dirs = V(core_idx,:);
                    norms = vecnorm(dirs,2,2);
                    norms(norms < 1e-12) = 1;
                    core_dirs_full = dirs ./ norms;
                end
            end
        end

        if ~includeCore && ~collapse_base_i
            base_center = seam_base_center;
            [V, base_center_idx] = append_vertices(V, base_center);
            base_faces = connect_ring_to_apex(seam_indices{i}, base_center_idx);
            if ~isempty(base_faces)
                tri = base_faces(1,:);
                v1 = V(tri(1),:);
                v2 = V(tri(2),:);
                v3 = V(tri(3),:);
                base_normal = cross(v2 - v1, v3 - v1);
                if dot(base_normal, -orientation) < 0
                    base_faces = base_faces(:, [1 3 2]);
                end
            end
            F = [F; base_faces]; %#ok<AGROW>
        end

        force_pointy_tip_i = pending_pointy_tip && ~collapse_base_i && (r_base_i > 2 * min_spacing);
        if force_cylinder
            force_pointy_tip_i = false;
        end

    minimal_seam = false;
    alpha_i = alpha_nominal_default_i;
    z_seam_i = z_seam_nominal_default_i;
    seam_radius_i = r_tip_nominal_i;

        if tipSphereRadius_i > 0
            [seam_candidate, alpha_candidate, z_seam_candidate, tangentValid_i] = solve_tip_tangent(max(r_base_i, 0), coreRadius, z_tip_center_i, tipSphereRadius_i);
            if tangentValid_i
                seam_radius_i = seam_candidate;
                alpha_i = alpha_candidate;
                z_seam_i = z_seam_candidate;
            else
                seam_radius_i = min(max(min(r_base_i, tipSphereRadius_i), 0), tipSphereRadius_i);
                ratio = min(1, max(seam_radius_i / max(tipSphereRadius_i, 1e-12), 0));
                alpha_i = max(1e-6, asin(ratio));
                z_seam_i = max(z_base_i + 1e-6, z_tip_center_i + tipSphereRadius_i * cos(alpha_i));
            end

            if z_seam_i <= z_base_i + 1e-9
                seam_radius_i = min(max(min(r_base_i, tipSphereRadius_i), 0), tipSphereRadius_i);
                ratio = min(1, max(seam_radius_i / max(tipSphereRadius_i, 1e-12), 0));
                alpha_i = max(1e-6, asin(ratio));
                z_seam_i = max(z_base_i + 1e-6, z_tip_center_i + tipSphereRadius_i * cos(alpha_i));
            end

            if seam_radius_i > 0 && (2 * seam_radius_i < min_spacing)
                minimal_seam = true;
                desired_radius = min(max(min_spacing / 2, seam_radius_i), tipSphereRadius_i);
                seam_radius_i = desired_radius;
                ratio = min(1, max(seam_radius_i / max(tipSphereRadius_i, 1e-12), 0));
                alpha_i = max(1e-6, asin(ratio));
                z_seam_i = max(z_base_i + 1e-6, z_tip_center_i + tipSphereRadius_i * cos(alpha_i));
            end
        else
            seam_radius_i = 0;
            alpha_i = pi/2;
            z_seam_i = max(z_base_i + 1e-6, z_tip_center_i);
        end

        seam_radius_i = min(seam_radius_i, r_base_i);
        tip_cap_height = tipSphereRadius_i * max(0, cos(alpha_i));
        tip_cap_radius = tipSphereRadius_i * max(0, sin(alpha_i));

        collapse_tip_i = (force_pointy_tip_i || (((seam_radius_i <= tip_collapse_tol) && ~minimal_seam) || ((tip_cap_height < min_spacing) && (tip_cap_radius < min_spacing))));
        if force_cylinder
            collapse_tip_i = false;
        end
        if collapse_tip_i
            seam_radius_i = 0;
        end

        minBaseRadiusMetric = min(minBaseRadiusMetric, r_base_i);
        minSeamRadiusMetric = min(minSeamRadiusMetric, seam_radius_i);
        tipCollapsedAll = tipCollapsedAll && collapse_tip_i;

    tip_center = z_tip_center_i * orientation;

    produced_spike_count = produced_spike_count + 1;

        L_cone = max(1e-6, z_seam_i - z_base_i);
        base_count_seq = max(3, base_segments_effective);
        if collapse_tip_i
            seam_segments_effective = max(3, base_count_seq);
        elseif minimal_seam
            seam_segments_effective = 3;
        else
            seam_segments_effective = ring_segment_count(seam_radius_i, min_spacing, azimuth_spacing_factor);
            seam_segments_effective = max(3, seam_segments_effective);
        end
        axial_steps = max(1, ceil(L_cone / (min_spacing * axial_spacing_factor)));
        num_cone_steps = axial_steps;

        prev_idx = seam_indices{i};
        seam_is_point = numel(prev_idx) == 1;

        for k = 1:num_cone_steps
            t = k / num_cone_steps;              % 0 at base, 1 at seam
            ring_radius = (1 - t) * r_base_i + t * seam_radius_i;
            ring_center = ((1 - t) * z_base_i + t * z_seam_i) * orientation;
            if ((force_pointy_tip_i && k == num_cone_steps) || (ring_radius <= tip_collapse_tol && k == num_cone_steps))
                [V, idx_point] = append_vertices(V, ring_center);
                if numel(prev_idx) > 1
                    F = [F; connect_ring_to_apex(prev_idx, idx_point)]; %#ok<AGROW>
                elseif numel(prev_idx) == 1
                    % nothing to connect; two points on axis
                end
                prev_idx = idx_point;
                seam_is_point = true;
                break;
            end

            if k == num_cone_steps
                segments_curr = seam_segments_effective;
            else
                segments_curr = ring_segment_count(ring_radius, min_spacing, azimuth_spacing_factor);
            end
            segments_curr = max(3, segments_curr);

            angles_curr = circle_angles(segments_curr);
            ring_pts = seam_ring_points(ring_center, u, v, ring_radius, angles_curr);
            [V, idx_ring] = append_vertices(V, ring_pts);
            if includeCore
                core_mask(idx_ring) = false;
            end

            if numel(prev_idx) == 1
                F = [F; connect_ring_to_apex(idx_ring, prev_idx)]; %#ok<AGROW>
            else
                faces_bridge = bridge_rings_idx(prev_idx(:)', idx_ring(:)', V, orientation);
                F = [F; faces_bridge]; %#ok<AGROW>
            end

            prev_idx = idx_ring;
        end

        seam_ring_idx = prev_idx;

        % --- Patch B: Spherical tip cap (ring-based) ---
        apex = tip_center + tipSphereRadius_i * orientation;
        if seam_is_point
            % Seam already collapsed to a point; cone-to-tip faces created in loop.
        elseif collapse_tip_i
            [V, idx_apex] = append_vertices(V, apex);
            if includeCore
                core_mask(idx_apex) = false;
            end
            if numel(seam_ring_idx) > 1
                F = [F; connect_ring_to_apex(seam_ring_idx, idx_apex)]; %#ok<AGROW>
            end
        else
            cap_arc_length = tipSphereRadius_i * max(alpha_i, 1e-6);
            if minimal_seam
                n_cap_rings = 0;
            else
                n_cap_rings = max(1, ceil(cap_arc_length / (min_spacing * cap_spacing_factor)));
            end
            prev_cap_idx = seam_ring_idx;
            for j = 1:n_cap_rings
                t_cap = (n_cap_rings - j + 1) / (n_cap_rings + 1);
                psi = alpha_i * t_cap;
                ring_radius = tipSphereRadius_i * sin(psi);
                ring_center = tip_center + (tipSphereRadius_i * cos(psi)) * orientation;
                segments_cap = ring_segment_count(ring_radius, min_spacing, cap_spacing_factor);
                segments_cap = max(3, segments_cap);
                angles_cap = circle_angles(segments_cap);
                ring_pts = seam_ring_points(ring_center, u, v, ring_radius, angles_cap);
                [V, idx_ring] = append_vertices(V, ring_pts);
                if includeCore
                    core_mask(idx_ring) = false;
                end
                faces_bridge = bridge_rings_idx(prev_cap_idx(:)', idx_ring(:)', V, orientation);
                F = [F; faces_bridge]; %#ok<AGROW>
                prev_cap_idx = idx_ring;
            end

            [V, idx_apex] = append_vertices(V, apex);
            if includeCore
                core_mask(idx_apex) = false;
            end
            F = [F; connect_ring_to_apex(prev_cap_idx, idx_apex)]; %#ok<AGROW>
        end

        % --- Patch C: Toroidal fillet (suppressed for now) ---
        % Intentionally disabled to focus on core sphere + spike body + tip
    end

    if produced_spike_count == 0
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
    F = [core_F_current; F];
    
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
        'coreRadius', coreRadius, ...
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
        'genVolume', genVolume, ...
        'volRes', volRes, ...
        'volPadding', volPadding, ...
        'volAlpha', volAlphaInp, ...
        'volAdaptive', volAdaptive, ...
        'volDxMax', volDxMax, ...
        'volDxMin', volDxMin, ...
        'volBlockSize', volBlockSz, ...
        'volCriterion', volCriterion ...
    );

    urchin.Metrics = struct( ...
        'MinimumSpacing', min_spacing, ...
        'SpikeBaseRadius', spikeBaseRadiusMetric, ...
        'SpikeBaseMaxima', r_base_max_per_spike, ...
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
        inside_has_value = false;
        % Try inpolyhedron (File Exchange) if available
        try
            inside = inpolyhedron(urchin.Faces, urchin.Vertices, P);
            inside_has_value = true;
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
                inside_has_value = true;
            catch
                % As a last resort, mark nothing inside
                inside = false(size(P,1),1);
                inside_has_value = true;
            end
        end
        if ~inside_has_value
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
    fprintf('B-Rep Urchin created successfully in %.2f seconds.\n', elapsedTime);
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

function [r_base_nominal, r_tip_nominal, alpha_nominal, z_seam_nominal] = compute_nominal_spike_profile(coreRadius, spikeLength, tipSphereRadius, conicalityEffective, r_base_max_local, scale)
%COMPUTE_NOMINAL_SPIKE_PROFILE Baseline spike profile honoring tip tangency and conicality.

    r_base_max_local = max(0, min(r_base_max_local, 0.999 * coreRadius));
    spikeLength = max(spikeLength, 0);
    z_apex = coreRadius + spikeLength;
    tipSphereRadius = max(0, min(tipSphereRadius, spikeLength));
    z_tip_center = z_apex - tipSphereRadius;

    if tipSphereRadius > 0 && r_base_max_local > 0
        r_base = min(max(min(tipSphereRadius, r_base_max_local), 0), r_base_max_local);
        for iter = 1:8
            [r_tip_candidate, ~, ~, tangentValid] = solve_tip_tangent(r_base, coreRadius, z_tip_center, tipSphereRadius);
            if ~tangentValid
                r_tip_candidate = min(r_base, tipSphereRadius);
            end
            if conicalityEffective >= 0
                r_base_target = r_tip_candidate + conicalityEffective * (r_base_max_local - r_tip_candidate);
            else
                r_base_target = max(0, (1 + conicalityEffective) * r_tip_candidate);
            end
            r_base_target = min(r_base_target, r_base_max_local);
            if abs(r_base_target - r_base) <= 1e-9 * max(scale, 1)
                r_base = r_base_target;
                break;
            end
            r_base = r_base_target;
        end

        [r_tip_nominal, alpha_nominal, z_seam_nominal, tangentValid] = solve_tip_tangent(r_base, coreRadius, z_tip_center, tipSphereRadius);
        if ~tangentValid
            r_tip_nominal = min(r_base, tipSphereRadius);
            ratio = min(1, max(r_tip_nominal / max(tipSphereRadius, 1e-12), 0));
            alpha_nominal = max(1e-6, asin(ratio));
            z_seam_nominal = z_tip_center + tipSphereRadius * cos(alpha_nominal);
        end
        r_base_nominal = r_base;
    else
        r_base_nominal = min(r_base_max_local, max(0, min(tipSphereRadius, r_base_max_local)));
        r_tip_nominal = 0;
        alpha_nominal = pi/2;
        z_seam_nominal = z_apex;
    end

    r_base_nominal = max(r_base_nominal, 0);
    r_tip_nominal = max(min(r_tip_nominal, tipSphereRadius), 0);
    alpha_nominal = min(max(alpha_nominal, 1e-6), pi/2);
    z_seam_nominal = max(z_seam_nominal, coreRadius);
end

function r_base_limit = tangent_base_limit(coreRadius, tipSphereRadius, spikeLength)
%TANGENT_BASE_LIMIT Maximum base radius for a spike length to remain tangent.

    if coreRadius <= 0
        r_base_limit = 0;
        return;
    end

    if tipSphereRadius <= 0
        r_base_limit = coreRadius;
        return;
    end

    spikeLength = max(spikeLength, 0);
    total_height = coreRadius + spikeLength;
    tip_total = coreRadius + tipSphereRadius;

    if total_height <= tip_total
        r_base_limit = 0;
        return;
    end

    sin_alpha = tip_total / total_height;
    sin_alpha = min(max(sin_alpha, 0), 1);
    cos_alpha = sqrt(max(0, 1 - sin_alpha^2));
    r_base_limit = coreRadius * cos_alpha;
end

function volume = computeMeshVolume(vertices, faces)
%COMPUTEMESHVOLUME Calculate enclosed volume of a triangular surface mesh.
    if isempty(vertices) || isempty(faces)
        volume = 0;
        return;
    end
    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);
    tetraVol = dot(v1, cross(v2, v3, 2), 2) / 6;
    volume = abs(sum(tetraVol));
end

