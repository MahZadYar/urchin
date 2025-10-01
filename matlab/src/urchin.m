function [mesh, diagnostics] = urchin(varargin)
% URCHIN  Deterministic, curvature-aware B-Rep nano-urchin mesher
% =========================================================================
% Builds a COMSOL-ready surface by stitching parametric patches:
%   - Core sphere with circular holes at spike footprints
%   - Conical frusta for spike bodies (smoothly controlled by sc∈[0,1])
%   - Spherical caps for spike tips (equal-angle ring spacing)
%   - Toroidal fillets blending core↔cone seams (G1-tangent, geometric)
%
% Point sampling is deterministic and curvature-aware:
%   - Spherical patches via Fibonacci lattice or geodesic icosphere
%   - Cones/fillets via structured grids in parameter space
%   - Seam edges share identical, precomputed rings for watertight joins
%
% New functionality:
%   - Spike fluctuations: per-spike length jitter (sf) with uniform/random (Sobol or seeded RNG).
%   - Orientation distributions: Fibonacci-like uniform or random (seeded).
%   - Unified meshing refinement (meshRefine) with three tunables in code
%     (CF_CORE, CF_BASE, CF_RINGS) that scale together; curvature-aware
%     counts keep quality where curvature is high (tips/fillets/seams).
%   - Volume representations:
%       • Dense uniform voxel grid via inpolyhedron/alphaShape fallback
%       • Adaptive sparse octree (VDB-style) with user control of voxel
%         size range (volDxMin/Max), leaf block size, and refinement
%         criterion (boundary/distance/curvature/hybrid). Stored in
%         mesh.VolumeOctree for sparse workflows; optional densification.
% Usage:
%   mesh = urchin('cr',30,'sl',15,'ns',50,'st',5);
%   [mesh, diagnostics] = urchin('cr',25,'sl',10,'ns',75);
%   fprintf('Watertight: %d\n', diagnostics.IsWatertight);
%
% Inputs (name-value):
%   cr           Core radius (nm). Default 1
%   sl           Spike length from core surface (nm). Default 2
%   ns           Number of spikes. Default 100
%   st           Spike seam diameter (nm) at the distal cut plane located at
%                distance (cr+sl) from center. A spherical cap is added beyond
%                the cut, tangent to the cone/cylinder at that seam. Default
%                st = sl/10.
%   sc           Conicality [0..1]. 0=cylinder; 1=full taper away from core.
%   density      Base point density (pts / nm^2). Default 10
%   useFillet    Enable toroidal fillet patch. Default true (fillet radius
%                defaults to st/2 and is clamped relative to base radius)
%
%   meshRefine   Unified refinement scalar scaling all feature counts. Default 1.0
%   sf           Spike length fluctuation factor [0..1]. Default 0.0
%   flucMethod   'uniform' (Sobol) or 'random' (seeded by ns). Default 'uniform'
%   distMethod   'uniform' (Fibonacci-like) or 'random' (seeded by ns). Default 'uniform'
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
    addParameter(p, 'cr', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'sl', 2, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'ns', 50, @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0}));
    addParameter(p, 'st', [],  @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    addParameter(p, 'sc', 1, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    addParameter(p, 'filletRatio', 0.25, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',0.5}));
    addParameter(p, 'density', 10, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'useFillet', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'includeCore', true, @(x)islogical(x) && isscalar(x));
    % Unified meshing refinement scalar (scales all counts proportionally)
    addParameter(p, 'meshRefine', 1.0, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    % Spike fluctuations and orientation distribution
    addParameter(p, 'sf', 1, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
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
    % Deterministic resolution controls (geometry-independent)
    addParameter(p, 'nAzimuth', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',3}));
    addParameter(p, 'nConeRings', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',3}));
    addParameter(p, 'nCapRings', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',3}));
    addParameter(p, 'nFilletRings', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',3}));
    addParameter(p, 'nCorePoints', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',100}));
    parse(p, varargin{:});

    % Assign parsed inputs to variables
    cr = p.Results.cr;
    sl = p.Results.sl;
    ns = p.Results.ns;
    st = p.Results.st;
    sc = p.Results.sc;
    base_density = p.Results.density;
    includeCore = p.Results.includeCore;
    refine      = p.Results.meshRefine;
    sf          = p.Results.sf;
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

    scale = 2 * (cr^2 + sl + (isempty(st) * (sl/20)/2)); % urchin scale for normalization
    % Seam diameter at distal cut plane (default: st = sl/10 unless specified)
    if isempty(st)
        st = sl / 10;
    end
    r_tip = st / 2;


    % Note: Seam axial location and cone tip are computed per-spike (fluctuations)

    % Spike height unused explicitly; seam defined at z=cr+sl

    % Distribute spike orientations (uniform Fibonacci or random)
    switch distMethod
        case 'uniform'
            % Golden spiral distribution across sphere via (theta,phi)
            Phi = (1 + sqrt(5)) / 2;  % golden ratio
            i = (1:ns)';
            theta = acos(1 - 2 .* i ./ (ns + 1));
            phi   = mod(2 * pi * i / Phi, 2 * pi);
            spike_orientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
        case 'random'
            rng(ns,'twister');
            phi = 2 * pi * rand(ns, 1);
            cosTheta = 2 * rand(ns, 1) - 1;
            theta = acos(cosTheta);
            spike_orientations = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    end

    % Spike length fluctuations
    switch flucMethod
        case 'uniform'
            flucs = net(sobolset(1), ns);   % values in ~[0,1]
        case 'random'
            rng(ns,'twister');
            flucs = rand(ns,1);            % values in ~[0,1]
        case 'gaussian'
            rng(ns,'twister');
            flucs = randn(ns,1);            % values in (-∞,∞)
    end

    flucs = flucs - min(flucs);     % shift to [0,∞)
    flucs = flucs / mean(flucs);   % normalize to mean on 1
    flucs = flucs -1;               % shift to [-1,∞) with mean 0
    flucs = sf * sl * flucs;        % scale to full spike length range [-sf*sl, +∞) with mean 0

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
    % on a sphere of radius cr is:
    %   r_base_max = cr * sin(α_max)
    % This guarantees non-overlapping base circles (they may be tangent for at
    % least one closest pair) while being as large as analytically possible
    % under the given orientation set.
    if ns == 0
        theta_min = pi; % arbitrary; no spikes
    elseif ns == 1
        theta_min = pi; % single spike: allow up to hemisphere
    else
        % Compute maximum cosine (excluding self) to find minimal angular sep
        Ddot = spike_orientations * spike_orientations.'; % ns x ns dot products
        Ddot(1:ns+1:end) = -Inf; % ignore diagonal
        maxDot = max(Ddot(:));   % largest cosine between distinct spikes
        maxDot = min(1, max(-1, maxDot));
        theta_min = acos(maxDot); % smallest angle between two spikes
    end
    alpha_base_max = 0.5 * theta_min;           % tangent half-angle

    theta_max = acos(cr/(cr+sl+st));
    alpha_base_max = min(theta_max, alpha_base_max);      % clamp to tangent angle
    r_base_max = cr * sin(alpha_base_max);      % maximal non-overlapping base radius
    % numeric safety clamps
    r_base_max = min(r_base_max, 0.999 * cr);
    r_base_max = max(r_base_max, 1e-9);

    % Interpolate actual base radius using conicality sc:
    %   sc=0 -> r_base = r_tip (cylindrical)
    %   sc=1 -> r_base = r_base_max (widest non-overlapping base)
    %   sc~>-1 -> r_base ~ 0 (base smaller than tip for negative sc)
    r_base = max(1e-6, sc * r_base_max + (1 - sc) * r_tip);
    z_base  = sqrt(max(1e-6, cr^2 - r_base^2));

    % Cone design base (independent of any fillet) for taper control
    r_cone_design = max(r_tip, r_tip + sc * (r_base_max - r_tip));

    % Points container not used; direct append to V

    % Curvature-aware density multipliers
    % Tunable base factors (adjust manually if desired)
    CF_CORE  = 1.0;   % baseline for core sphere points
    CF_BASE  = 1.0;   % baseline for base ring azimuthal segmentation
    CF_RINGS = 1.0;   % baseline for cap/fillet ring counts
    k_sphere = 1 / cr;
    k_tip    = 2 / max(st, 1e-9);   % 1 / r_tip proxy using seam diameter
    dens_mult_sphere = 0.8 + 0.2 * (k_sphere / (k_sphere + k_tip));
    dens_mult_tip    = 3.0;         % higher
    dens_mult_cone   = 1.25;        % higher to keep spikes

    % Edge segment counts (also used for seam alignment)
    % Global azimuthal segmentation based on the largest base ring (r_base_max),
    % scaled by unified refinement factor. Individual spikes use a fraction.
    n_theta_base = max(24, ceil(CF_BASE * refine * 2 * pi * r_base_max * sqrt(base_density * dens_mult_cone)));
    % Deterministic override for azimuthal segmentation
    if ~any(strcmp('nAzimuth', p.UsingDefaults)) && ~isempty(p.Results.nAzimuth)
        n_theta_base = p.Results.nAzimuth;
    end
    % Note: Actual per-spike azimuth angles are set via 'angles_uniform'

    %% 3) Build full core sphere FIRST (no trimming yet)
    V = zeros(0,3); F = zeros(0,3);
    core_subdiv = max(3, ceil(0.25 * sqrt(base_density * dens_mult_sphere))); % NEEDED TO BE ADAPTED TO RES
    [core_V_full, core_F_full] = triangulate_icosphere(cr, core_subdiv);
    [V, core_idx] = append_vertices(V, core_V_full); % core vertices global indices
    core_F_current = core_F_full + core_idx(1) - 1;   % active core faces (global indices)
    core_dirs_full = core_V_full / cr;                % unit directions for initial removal test
    core_mask = false(size(V,1),1);
    core_mask(core_idx) = true;

    % Precompute cutoff parameters (same for all spikes currently)
    base_angle   = asin(min(1, r_base / cr));
    cos_cutoff   = cos(base_angle);

    %% 4) Generate spikes one-by-one WITH per-spike core trimming & stitching
    fprintf('Generating spikes with on-the-fly core trimming & stitching for %d spikes...\n', ns);
    seam_indices = cell(ns,1);
    for i = 1:ns
        orientation = spike_orientations(i, :);
        % Local orthonormal frame:
        %  - orientation: spike axis (unit vector)
        %  - u, v: span orthogonal plane (used for circular rings)
        [u, v] = plane_vectors(orientation);

        % Determine azimuthal segmentation for this spike (must match all rings)
        n_seg = max(5, round(0.6 * n_theta_base));
        if ~any(strcmp('nAzimuth', p.UsingDefaults)) && ~isempty(p.Results.nAzimuth)
            n_seg = p.Results.nAzimuth;
        end
        angles_uniform = linspace(0, 2*pi, n_seg+1); angles_uniform(end) = [];

        % Refresh active core vertex set and their directions for this spike
        core_idx = find(core_mask);
        if isempty(core_idx)
            core_dirs_full = zeros(0,3);
        else
            dirs = V(core_idx,:);
            norms = vecnorm(dirs,2,2);
            norms(norms < 1e-12) = 1;
            core_dirs_full = dirs ./ norms;
        end

        core_loop = [];
        seam_indices{i} = [];

        % === Per-spike core trimming ===
        % Determine newly removed core vertices for this spike.
        % Angular test over current active core vertices
        core_dot = core_dirs_full * orientation'; % cos(angle) per active core vertex
        newRemoveMask = core_dot > cos_cutoff + 1e-12; % inside cutoff
        if ~any(newRemoveMask) && ~isempty(core_idx)
            % Fallback: remove the single vertex whose direction is closest to the base center
            [~, fallbackPos] = max(core_dot);
            if ~isnan(fallbackPos) && fallbackPos >= 1 && fallbackPos <= numel(newRemoveMask)
                newRemoveMask(fallbackPos) = true;
            end
        end
        if any(newRemoveMask)
            % Remove faces referencing any removed vertex
            verts_remove_global = core_idx(newRemoveMask);
            core_mask(verts_remove_global) = false;
            faceRemoveMask = any(ismember(core_F_current, verts_remove_global), 2);

            % Save removed faces for boundary-vertex inference BEFORE deleting them
            removedFaces = core_F_current(faceRemoveMask, :);

            % Now delete those faces from the current core face list
            core_F_current(faceRemoveMask, :) = [];

            % Determine candidate boundary vertices:
            % vertices that participated in removed faces but were NOT removed
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

        % --- Create seam ring AFTER trimming so it participates in next spike's core ---
        seam_base_center = z_base * orientation;
        seam_base_ring = seam_ring_points(seam_base_center, u, v, r_base, angles_uniform);
        [V, seam_indices{i}] = append_vertices(V, seam_base_ring);
        core_mask(seam_indices{i}) = true;

        if numel(core_loop) >= 3
            seam_ring = seam_indices{i}(:)';
            new_bridge_faces = bridge_rings_idx(core_loop, seam_ring, V, orientation);
            core_F_current = [core_F_current; new_bridge_faces]; %#ok<AGROW>
            core_mask(core_loop) = true;
        end

        % Refresh stored directions for debugging / subsequent logic
        core_idx = find(core_mask);
        if isempty(core_idx)
            core_dirs_full = zeros(0,3);
        else
            dirs = V(core_idx,:);
            norms = vecnorm(dirs,2,2);
            norms(norms < 1e-12) = 1;
            core_dirs_full = dirs ./ norms;
        end

        % Per-spike seam and tip geometry with fluctuations
        % z_seam_i: axial position of cone/tip seam; alpha_i: tangent half-angle; r_tip_i: tip sphere radius
        z_seam_i = cr + sl + flucs(i);
        alpha_i  = solve_alpha_conicality(r_base_max, r_tip, cr, sl + flucs(i), sc);
        r_tip_i  = r_tip / max(1e-9, sin(alpha_i));
        z_tip_center = z_seam_i - r_tip_i * cos(alpha_i);
        tip_center   = z_tip_center * orientation;

        L_cone = max(1e-6, z_seam_i - z_base);
        num_cone_rings = max(5, ceil(CF_RINGS * refine * 0.5 * (1 + k_sphere * L_cone)));
        if ~any(strcmp('nConeRings', p.UsingDefaults)) && ~isempty(p.Results.nConeRings)
            num_cone_rings = p.Results.nConeRings;
        end
        prev_idx = seam_indices{i};
        for k = 1:num_cone_rings
            t = k / num_cone_rings;              % 0 at base, 1 at seam
            % Smooth conicality: sc=0 → cylinder (r=r_tip), sc=1 → full taper
            dr = max(0, r_cone_design - r_tip);
            ring_radius = r_tip + (1 - t) * dr;
            ring_center = ((1 - t) * z_base + t * z_seam_i) * orientation; % z_base_local removed (fillet suppressed)
            ring_pts = seam_ring_points(ring_center, u, v, ring_radius, angles_uniform);
            [V, idx_ring] = append_vertices(V, ring_pts);
            core_mask(idx_ring) = false;
            if ~isempty(prev_idx)
                F = [F; connect_rings(prev_idx, idx_ring)]; %#ok<AGROW>
            end
            prev_idx = idx_ring;
        end

        % --- Patch B: Spherical tip cap (ring-based) ---
        % 4 rings at equal polar angle spacing from seam (ψ≈α) to apex (ψ≈0).
        apex = tip_center + r_tip_i * orientation;
        n_cap_rings = max(3, round(CF_RINGS * refine * 4)); % default ~4 rings
        if ~any(strcmp('nCapRings', p.UsingDefaults)) && ~isempty(p.Results.nCapRings)
            n_cap_rings = p.Results.nCapRings;
        end
        prev_cap_idx = prev_idx; % start from cone seam indices
        for j = 1:n_cap_rings
            % Equal angular spacing from the spherical cap center (linear in psi)
            % j=1 near seam (psi≈alpha), j=n near apex (psi≈0)
            t = (n_cap_rings - j + 1) / (n_cap_rings + 1);
            psi = alpha_i * t;
            ring_radius = r_tip_i * sin(psi);
            ring_center = tip_center + (r_tip_i * cos(psi)) * orientation;
            ring_pts = seam_ring_points(ring_center, u, v, ring_radius, angles_uniform);
            [V, idx_ring] = append_vertices(V, ring_pts);
            core_mask(idx_ring) = false;
            F = [F; connect_rings(prev_cap_idx, idx_ring)]; %#ok<AGROW>
            prev_cap_idx = idx_ring;
        end
        % Apex fan
        [V, idx_apex] = append_vertices(V, apex);
        core_mask(idx_apex) = false;
        F = [F; connect_ring_to_apex(prev_cap_idx, idx_apex)]; %#ok<AGROW>

        % --- Patch C: Toroidal fillet (suppressed for now) ---
        % Intentionally disabled to focus on core sphere + spike body + tip
    end

    % After all spikes processed, append remaining active core faces
    % Append remaining (trimmed) core faces
    F = [core_F_current; F];

    % Final weld then construct surface mesh
    [Vw, ~, ic] = uniquetol(V, 1e-9, 'ByRows', true);
    Fw = ic(F);
    mesh = surfaceMesh(Vw, Fw);

    % Built-in mesh quality diagnostics
    diagnostics = struct( ...
        'IsWatertight',       isWatertight(mesh), ...
        'IsEdgeManifold',     isEdgeManifold(mesh, false), ...
        'IsOrientable',       isOrientable(mesh), ...
        'IsSelfIntersecting', isSelfIntersecting(mesh), ...
        'IsVertexManifold',   isVertexManifold(mesh) ...
    );

    if ~diagnostics.IsWatertight
        warning('Urchin mesh is not watertight. Consider adjusting parameters.');
    end
    if ~diagnostics.IsEdgeManifold
        warning('Urchin mesh has non-manifold edges.');
    end
    if ~diagnostics.IsVertexManifold
        warning('Urchin mesh has non-manifold vertices.');
    end
    if diagnostics.IsSelfIntersecting
        warning('Urchin mesh contains self-intersections.');
    end

    if nargout < 2
        diagnostics = [];
    end
    %% 5) Build a volumetric mask from the surface mesh by voxelizing
    % the shape using either inpolyhedron (if available) or alphaShape fallback.
    if genVolume && ~volAdaptive
        Vmin = min(mesh.Vertices,[],1);
        Vmax = max(mesh.Vertices,[],1);
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
        % Try inpolyhedron (File Exchange) if available
        try
            inside = inpolyhedron(mesh.Faces, mesh.Vertices, P);
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
                shpVol = alphaShape(mesh.Vertices, alphaVol);
                inside = inShape(shpVol, P(:,1), P(:,2), P(:,3));
            catch
                % As a last resort, mark nothing inside
                inside = false(size(P,1),1);
            end
        end
        VolumeMask = reshape(logical(inside), size(XX));
        mesh.VolumeMask = VolumeMask;
        mesh.VoxelGrid.X = xs;
        mesh.VoxelGrid.Y = ys;
        mesh.VoxelGrid.Z = zs;
        mesh.VoxelSize = dx;
        mesh.Bounds = [Vmin; Vmax];
    elseif genVolume && volAdaptive
        % Adaptive sparse octree voxelization (VDB-style)
        Vmin = min(mesh.Vertices,[],1);
        Vmax = max(mesh.Vertices,[],1);
        span = Vmax - Vmin;
        Vmin = Vmin - volPadding * span;
        Vmax = Vmax + volPadding * span;
        maxDim = max(span);
        if isempty(volDxMax), volDxMax = maxDim / 128; end
        if isempty(volDxMin), volDxMin = volDxMax / 4; end
        insideFn = make_inside_tester(mesh, volAlphaInp);
        leaves = build_adaptive_octree(Vmin, Vmax, volDxMax, volDxMin, volBlockSz, insideFn, volCriterion);
        mesh.VolumeOctree.Leaves = leaves;
        mesh.VolumeOctree.BlockSize = volBlockSz;
        mesh.VolumeOctree.DxMin = volDxMin;
        mesh.VolumeOctree.DxMax = volDxMax;
        mesh.VolumeOctree.Bounds = [Vmin; Vmax];
    end

    elapsedTime = toc;
    fprintf('B-Rep Urchin created successfully in %.2f seconds.\n', elapsedTime);
    fprintf('Generated mesh with %d vertices and %d faces.\n', size(mesh.Vertices, 1), size(mesh.Faces, 1));

    %% 6) Visualize if no output is requested
    if nargout == 0
        viewer = viewer3d;
        surfaceMeshShow(mesh, Parent=viewer, WireFrame=true);
        surfaceMeshShow(mesh, Parent=viewer, Alpha=0.5);
    end
end

%% --- Helper Functions ---

function insideFn = make_inside_tester(mesh, volAlpha)
    % Build a robust point-in-mesh tester closure. Prefers inpolyhedron if available.
    F = mesh.Faces; V = mesh.Vertices;
    hasInpoly = false; %#ok<NASGU>
    try
        % Probe inpolyhedron availability without calling it on large sets
        test = inpolyhedron(F, V, mean(V,1)); %#ok<NASGU>
        hasInpoly = true; %#ok<NASGU>
        insideFn = @(P) inpolyhedron(F, V, P);
        return;
    catch
        % Fallback to alphaShape
        if isempty(volAlpha)
            bbox = max(V,[],1) - min(V,[],1);
            volAlpha = 0.08 * max(bbox); % heuristic
        end
        shp = alphaShape(V(:,1), V(:,2), V(:,3), volAlpha);
        insideFn = @(P) inShape(shp, P(:,1), P(:,2), P(:,3));
    end
end

function leaves = build_adaptive_octree(bmin, bmax, dxMax, dxMin, blockSize, insideFn, criterion)
    % Recursively build a sparse octree of voxel blocks within [bmin,bmax].
    % Each leaf contains: origin [1x3], dx scalar, occupancy [blockSize^3 logical].
    leaves = recurse_node(bmin, bmax, dxMax);

    function outLeaves = recurse_node(minC, maxC, dx)
        center = 0.5*(minC + maxC);
        span = maxC - minC;
        % Sample corners to determine trivial inside/outside
        [CX,CY,CZ] = ndgrid([minC(1), maxC(1)], [minC(2), maxC(2)], [minC(3), maxC(3)]);
        Pc = [CX(:), CY(:), CZ(:)];
        inCorners = insideFn(Pc);
        if all(inCorners)
            % Full block leaf
            occ = true(blockSize, blockSize, blockSize);
            outLeaves = struct('origin',minC,'dx',dx,'occupancy',occ);
            return;
        elseif ~any(inCorners) && ~needs_refine(minC, maxC, insideFn, criterion)
            % Empty block leaf
            occ = false(blockSize, blockSize, blockSize);
            outLeaves = struct('origin',minC,'dx',dx,'occupancy',occ);
            return;
        end
        if dx/2 < dxMin
            % Reached finest resolution: build occupancy at this leaf
            [occ, leafDx] = rasterize_block(minC, maxC, blockSize, insideFn);
            outLeaves = struct('origin',minC,'dx',leafDx,'occupancy',occ);
            return;
        end
        % Subdivide into 8 children
        outLeaves = repmat(struct('origin',[],'dx',[],'occupancy',[]), 0, 1);
        for ix = 0:1
            for iy = 0:1
                for iz = 0:1
                    childMin = minC + (span/2).* [ix, iy, iz];
                    childMax = childMin + span/2;
                    outLeaves = [outLeaves; recurse_node(childMin, childMax, dx/2)]; %#ok<AGROW>
                end
            end
        end
    end
end

function tf = needs_refine(minC, maxC, insideFn, criterion)
    % Decide if a node needs refinement based on boundary/distance/curvature heuristics.
    % For speed, we implement a boundary band test only (others are stubs/placeholders).
    switch criterion
        case 'boundary'
            % Boundary test: if corners disagree, refine
            [CX,CY,CZ] = ndgrid([minC(1), maxC(1)], [minC(2), maxC(2)], [minC(3), maxC(3)]);
            Pc = [CX(:), CY(:), CZ(:)];
            inCorners = insideFn(Pc);
            tf = any(inCorners) && ~all(inCorners);
        case 'distance'
            % Placeholder: use boundary as proxy
            tf = needs_refine(minC, maxC, insideFn, 'boundary');
        case 'curvature'
            % Placeholder: use boundary as proxy
            tf = needs_refine(minC, maxC, insideFn, 'boundary');
        case 'hybrid'
            tf = needs_refine(minC, maxC, insideFn, 'boundary');
        otherwise
            tf = needs_refine(minC, maxC, insideFn, 'boundary');
    end
end

function [occ, dx] = rasterize_block(minC, maxC, blockSize, insideFn)
    % Rasterize a leaf block uniformly at resolution 'blockSize'.
    span = maxC - minC;
    dx = max(span) / blockSize;
    xs = minC(1) + (0.5:1:blockSize) * dx;
    ys = minC(2) + (0.5:1:blockSize) * dx;
    zs = minC(3) + (0.5:1:blockSize) * dx;
    [XX,YY,ZZ] = ndgrid(xs,ys,zs);
    P = [XX(:), YY(:), ZZ(:)];
    inside = insideFn(P);
    occ = reshape(logical(inside), size(XX));
end

function alpha = solve_alpha_conicality(r_base_max, r_tip, cr, sl, sc)
    % Solve for cone half-angle alpha given conicality sc ∈ [0,1]
    % sc=0 → cylinder with hemispherical tip (alpha=pi/2)
    % sc=1 → base radius equals r_base_max with G1 tangency at tip seam
    if sc <= 0
        alpha = pi/2; % perfect cylinder, hemispherical cap
        return;
    end
    if sc >= 1
        sc = 1; % clamp
    end

    z_seam = cr + sl; % seam axial position
    function e = err(a)
        rs = r_tip;                          % seam radius fixed by st
        rb = sc * r_base_max + (1 - sc) * rs; % base radius from conicality rule
        zb = sqrt(max(1e-12, cr^2 - rb^2));   % base z on core sphere
        zs = z_seam;                           % seam z on axis
        denom = max(1e-9, (zs - zb));
        slope_cone = (rs - rb) / denom;       % radial slope of frustum
        e = slope_cone + cot(a);              % G1 tangency condition
    end

    aL = 2*pi/180; aU = 89*pi/180; % search bounds (avoid 0 and 90 exactly)
    fL = err(aL); fU = err(aU);
    if sign(fL) == sign(fU)
        % Fallback: bias towards cylinder for small sc, sharper for large sc
        alpha = (1 - 0.6*sc) * (pi/2); % in [~0.2*pi, 0.5*pi]
        return;
    end
    for it = 1:60
        aM = 0.5*(aL + aU);
        fM = err(aM);
        if sign(fM) == sign(fL)
            aL = aM; fL = fM;
        else
            aU = aM; fU = fM;
        end
        if abs(fM) < 1e-9 || (aU - aL) < 1e-7
            break;
        end
    end
    alpha = 0.5*(aL + aU);
end

function [V, idx] = append_vertices(V, new_pts)
    % Append new points to vertex list and return their indices
    if isempty(new_pts)
        idx = zeros(0,1);
        return;
    end
    startIdx = size(V,1) + 1;
    V = [V; new_pts];
    idx = (startIdx:size(V,1)).';
end

function F = connect_rings(idxA, idxB)
    % Create quad strips between two rings and split into triangles
    n = numel(idxA);
    assert(numel(idxB) == n, 'Rings must have same segmentation');
    F = zeros(2*n, 3);
    for k = 1:n
        k2 = mod(k, n) + 1;
        a1 = idxA(k); a2 = idxA(k2);
        b1 = idxB(k); b2 = idxB(k2);
        F(2*k-1, :) = [a1, b1, b2];
        F(2*k,   :) = [a1, b2, a2];
    end
end

function F = connect_ring_to_apex(idxRing, idxApex)
    % Fan triangulation from ring to apex
    n = numel(idxRing);
    F = zeros(n, 3);
    for k = 1:n
        k2 = mod(k, n) + 1;
        F(k, :) = [idxRing(k), idxRing(k2), idxApex];
    end
end

function [V, F] = triangulate_icosphere(R, subdivisions)
    % Build a geodesic icosphere with given subdivisions and radius R
    t = (1 + sqrt(5)) / 2;
    verts = [
        -1,  t,  0;
         1,  t,  0;
        -1, -t,  0;
         1, -t,  0;
         0, -1,  t;
         0,  1,  t;
         0, -1, -t;
         0,  1, -t;
         t,  0, -1;
         t,  0,  1;
        -t,  0, -1;
        -t,  0,  1
    ];
    faces = [
        1,12,6; 1,6,2; 1,2,8; 1,8,11; 1,11,12;
        2,6,10; 6,12,5; 12,11,3; 11,8,7; 8,2,9;
        4,10,5; 4,5,3; 4,3,7; 4,7,9; 4,9,10;
        5,10,6; 3,5,12; 7,3,11; 9,7,8; 10,9,2
    ];
    % Normalize to unit sphere
    verts = verts ./ vecnorm(verts, 2, 2);
    % Subdivide
    for s = 1:subdivisions
        midCache = containers.Map('KeyType','char', 'ValueType','int32');
        newFaces = zeros(size(faces,1)*4, 3);
        newVerts = verts;
        for i = 1:size(faces,1)
            a = faces(i,1); b = faces(i,2); c = faces(i,3);
            [ab, newVerts, midCache] = midpointIndex(a,b,newVerts,midCache);
            [bc, newVerts, midCache] = midpointIndex(b,c,newVerts,midCache);
            [ca, newVerts, midCache] = midpointIndex(c,a,newVerts,midCache);
            idx4 = (i-1)*4;
            newFaces(idx4+1,:) = [a, ab, ca];
            newFaces(idx4+2,:) = [b, bc, ab];
            newFaces(idx4+3,:) = [c, ca, bc];
            newFaces(idx4+4,:) = [ab, bc, ca];
        end
        faces = newFaces;
        verts = newVerts;
    end
    % Scale to radius R
    V = R * (verts ./ vecnorm(verts,2,2));
    F = faces;
end

function [idx, newVerts, midCache] = midpointIndex(a,b,newVerts,midCache)
    % Return index of midpoint vertex between a and b, caching results
    if a > b
        key = sprintf('%d_%d', b,a);
    else
        key = sprintf('%d_%d', a,b);
    end
    if isKey(midCache, key)
        idx = midCache(key);
        return;
    end
    p = (newVerts(a,:) + newVerts(b,:)) / 2;
    p = p / norm(p);
    newVerts = [newVerts; p];
    idx = size(newVerts,1);
    midCache(key) = idx;
end

function [u, v] = plane_vectors(normal_vec)
    % Orthonormal basis (u,v) perpendicular to normal_vec
    n = normal_vec(:) / norm(normal_vec);
    if abs(n(1)) > 0.9
        other = [0;1;0];
    else
        other = [1;0;0];
    end
    u = cross(n, other); u = u / max(1e-12, norm(u));
    v = cross(n, u);     v = v / max(1e-12, norm(v));
end

function ring_pts = seam_ring_points(center, u, v, radius, angles)
    % Deterministic points on a circle with given center, basis {u,v}, radius
    ca = cos(angles); sa = sin(angles);
    ring_pts = [center(1) + radius * (u(1)*ca + v(1)*sa)', ...
                center(2) + radius * (u(2)*ca + v(2)*sa)', ...
                center(3) + radius * (u(3)*ca + v(3)*sa)'];
end

function Fpatch = bridge_rings_idx(inner_idx, outer_idx, V, orientation)
% BRIDGE_RINGS_IDX  Zipper triangulation between two index loops.
% This function is robust to unordered inputs, jagged/non-planar rings,
% misaligned starting points, and opposite winding orders.
% It guarantees all face normals point towards the given orientation.
%
% inner_idx, outer_idx: row vectors of vertex indices (can be unordered)
% V: The N-by-3 vertex coordinate matrix.
% orientation: A 1x3 vector normal to the plane of the inner ring.
% Returns faces referencing original V; does not add vertices.

    % --- Input Validation & Formatting ---
    inner_idx = inner_idx(:)'; outer_idx = outer_idx(:)';
    Ni = numel(inner_idx); No = numel(outer_idx);
    if Ni < 2 || No < 2
        Fpatch = zeros(0,3); return;
    end

    % --- Define Sorting Coordinate System for Inner Ring ---
    center_point = mean(V(inner_idx, :), 1);
    w = orientation(:)' / norm(orientation);
    vecs_from_center = bsxfun(@minus, V(inner_idx,:), center_point);
    [~, max_dist_idx] = max(sum(vecs_from_center.^2, 2));
    ref_vec = V(inner_idx(max_dist_idx), :) - center_point;
    u = ref_vec - dot(ref_vec, w) * w; u = u / norm(u);
    v = cross(w, u);
    
    % --- Angular Sorting of Inner Ring (Assumed Planar) ---
    inner_vecs = bsxfun(@minus, V(inner_idx, :), center_point);
    angles_inner = atan2(inner_vecs * v', inner_vecs * u');
    [~, sort_order_inner] = sort(angles_inner);
    inner_idx = inner_idx(sort_order_inner);

    % --- 3D Nearest-Neighbor Sort for Outer Ring ---
    if No > 2
        sorted_outer_idx = zeros(1, No);
        dists_from_center_sq = sum(bsxfun(@minus, V(outer_idx,:), center_point).^2, 2);
        [~, start_pos] = max(dists_from_center_sq);
        sorted_outer_idx(1) = outer_idx(start_pos);
        remaining_indices = 1:No;
        remaining_indices(start_pos) = [];
        for k = 2:No
            last_added_v_idx = sorted_outer_idx(k-1);
            remaining_v_indices = outer_idx(remaining_indices);
            coords_last_added = V(last_added_v_idx,:);
            coords_remaining = V(remaining_v_indices,:);
            distances_sq = sum(bsxfun(@minus, coords_remaining, coords_last_added).^2, 2);
            [~, closest_rem_pos] = min(distances_sq);
            sorted_outer_idx(k) = remaining_v_indices(closest_rem_pos);
            remaining_indices(closest_rem_pos) = [];
        end
        outer_idx = sorted_outer_idx;
    end

    % --- Align Start Point ---
    v_inner_start = V(inner_idx(1), :);
    distances_sq = sum(bsxfun(@minus, V(outer_idx,:), v_inner_start).^2, 2);
    [~, best_start_j_pos] = min(distances_sq);
    outer_idx = circshift(outer_idx, [0, 1 - best_start_j_pos]);

    % --- Correct Outer Ring Winding Direction (GLOBAL PROJECTION FIX) ---
    V_outer_on_plane = V(outer_idx,:) * [u', v'];
    signed_area = sum( V_outer_on_plane(1:end-1, 1) .* V_outer_on_plane(2:end, 2) ) - sum( V_outer_on_plane(2:end, 1) .* V_outer_on_plane(1:end-1, 2) );
    signed_area = signed_area + (V_outer_on_plane(end, 1) * V_outer_on_plane(1, 2) - V_outer_on_plane(1, 1) * V_outer_on_plane(end, 2));

    if signed_area < 0
        outer_idx = [outer_idx(1), fliplr(outer_idx(2:end))];
    end

    % --- Triangulation (Greedy Shortest-Diagonal Algorithm) ---
    Ftmp = zeros(Ni + No, 3); c = 0;
    i = 1; j = 1; ai = 0; aj = 0;
    
    while ai < Ni || aj < No
        curr_i_idx = inner_idx(i); curr_j_idx = outer_idx(j);
        inext = mod(i, Ni) + 1; jnext = mod(j, No) + 1;
        next_i_idx = inner_idx(inext); next_j_idx = outer_idx(jnext);
        new_face = zeros(1,3);

        if ai == Ni
            new_face = [curr_i_idx, next_j_idx, curr_j_idx];
            j = jnext; aj = aj + 1; 
        elseif aj == No
            new_face = [curr_i_idx, next_i_idx, curr_j_idx];
            i = inext; ai = ai + 1;
        else
            % --- Greedy Choice Logic (THE FIX) ---
            % This logic now solely relies on choosing the shorter of the
            % two possible diagonals to form the next triangle. This
            % avoids creating geometrically unsound, overlapping faces.
            d1 = sum((V(curr_j_idx, :) - V(next_i_idx, :)).^2);
            d2 = sum((V(curr_i_idx, :) - V(next_j_idx, :)).^2);
            
            if d1 < d2
                new_face = [curr_i_idx, next_i_idx, curr_j_idx];
                i = inext; ai=ai+1;
            else
                new_face = [curr_i_idx, next_j_idx, curr_j_idx];
                j = jnext; aj=aj+1;
            end
        end
        
        c = c + 1;
        
        % --- Enforce Outward-Facing Normal (Safety Check) ---
        p1 = V(new_face(1),:);
        p2 = V(new_face(2),:);
        p3 = V(new_face(3),:);
        face_normal = cross(p2 - p1, p3 - p1);

        if dot(face_normal, orientation) < 0
            Ftmp(c,:) = [new_face(1), new_face(3), new_face(2)];
        else
            Ftmp(c,:) = new_face;
        end
    end
    Fpatch = Ftmp(1:c,:);
end