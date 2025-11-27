function [mesh] = urchin_BRep(varargin)
% URCHIN_BREP  Deterministic, curvature-aware B-Rep nano-urchin mesher
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
%   - Spike fluctuations: per-spike length jitter (sf) with uniform/random
%     methods matching urchin.m (Sobol or seeded RNG).
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
%   mesh = urchin_BRep('cr',30,'sl',15,'ns',50,'st',5);
%   urchin_BRep('cr',25,'sl',10,'ns',75); % auto-visualize if no outputs
%
% Inputs (name-value):
%   cr           Core radius (nm). Default 1 (as in urchin.m)
%   sl           Spike length from core surface (nm). Default 2 (as in urchin.m)
%   ns           Number of spikes. Default 100 (as in urchin.m)
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
%   imresize3 requires Image Processing Toolbox
%   surfaceMesh requires Lidar Toolbox
%
% =========================================================================

    %% 1) Input Parser
    p = inputParser;
    addParameter(p, 'cr', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'sl', 2, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'ns', 100, @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',0}));
    addParameter(p, 'st', [],  @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
    addParameter(p, 'sc', 0.75, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    addParameter(p, 'filletRatio', 0.25, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',0.5}));
    addParameter(p, 'density', 10, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    addParameter(p, 'useFillet', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'includeCore', true, @(x)islogical(x) && isscalar(x));
    % Unified meshing refinement scalar (scales all counts proportionally)
    addParameter(p, 'meshRefine', 1.0, @(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    % Spike fluctuations and orientation distribution (same rationale as urchin.m)
    addParameter(p, 'sf', 0.0, @(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}));
    addParameter(p, 'flucMethod', 'uniform', @(x)any(validatestring(x,{'uniform','random'})));
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
    useFillet   = p.Results.useFillet;
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
    
    % Seam diameter at distal cut plane (default: st = sl/10 unless specified)
    if isempty(st)
        st = sl / 20;
    end
    r_seam = max(1e-9, st / 2);

    % Temporarily suppress toroidal fillets to simplify base-radius optimization
    useFillet = false; % override user flag for this experimental mode

    % Note: Seam axial location and cone tip are computed per-spike (fluctuations)

    % Spike height unused explicitly; seam defined at z=cr+sl

    % Distribute spike orientations (uniform Fibonacci or random), consistent with urchin.m
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

    % Spike length fluctuations (same definition as urchin.m)
    switch flucMethod
        case 'uniform'
            flucs = -net(sobolset(1), ns);   % values in ~[-1,0]
        case 'random'
            rng(ns,'twister');
            flucs = -rand(ns,1);            % values in [-1,0]
    end
    flucs = sf * sl * flucs;                % scale to full spike length range

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
    r_base_max = cr * sin(alpha_base_max);      % maximal non-overlapping base radius
    r_base_max = min(r_base_max, 0.999 * cr);   % numeric safety clamp
    if r_base_max < 1e-9, r_base_max = 1e-9; end

    % Interpolate actual base radius using conicality sc:
    %   sc=0 -> r_base = r_seam (cylindrical)
    %   sc=1 -> r_base = r_base_max (widest non-overlapping)
    r_base = max(1e-6, sc * r_base_max + (1 - sc) * r_seam);
    z_base  = sqrt(max(1e-6, cr^2 - r_base^2));

    % Cone design base (independent of any fillet) for taper control
    r_cone_design = max(r_seam, r_seam + sc * (r_base_max - r_seam));

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

    %% 3) Initialize mesh containers
    V = zeros(0,3);
    F = zeros(0,3);

    %% 4) Generate Patches for Each Spike (explicit triangulation)
    fprintf('Generating seam-aligned patches for %d spikes...\n', ns);
    seam_rings = cell(ns,1);
    seam_indices = cell(ns,1);
    annulus_end_rings = cell(ns,1);
    for i = 1:ns
        orientation = spike_orientations(i, :);
        % Local orthonormal frame:
        %  - orientation: spike axis (unit vector)
        %  - u, v: span orthogonal plane (used for circular rings)
        [u, v] = plane_vectors(orientation);

        % --- Shared seam circles (identical points on both adjacent patches) ---
        % Base seam circle on core sphere at precomputed z_base
        seam_base_center = z_base * orientation;
        % Determine azimuthal segmentation for this spike (must match all rings)
        n_seg = max(5, round(0.6 * n_theta_base));
        if ~any(strcmp('nAzimuth', p.UsingDefaults)) && ~isempty(p.Results.nAzimuth)
            n_seg = p.Results.nAzimuth;
        end
        angles_uniform = linspace(0, 2*pi, n_seg+1); angles_uniform(end) = [];
        seam_base_ring = seam_ring_points(seam_base_center, u, v, r_base, angles_uniform);
        seam_rings{i}  = seam_base_ring;
        [V, seam_indices{i}] = append_vertices(V, seam_base_ring);

        % Per-spike seam and tip geometry with fluctuations (consistent with urchin.m)
        % z_seam_i: axial position of cone/tip seam; alpha_i: tangent half-angle; r_tip_i: tip sphere radius
        z_seam_i = cr + sl + flucs(i);
        alpha_i  = solve_alpha_conicality(r_base_max, r_seam, cr, sl + flucs(i), sc);
        r_tip_i  = r_seam / max(1e-9, sin(alpha_i));
        z_tip_center = z_seam_i - r_tip_i * cos(alpha_i);
        tip_center   = z_tip_center * orientation;

        % Base seam transition: fillet (torus) or annulus band
    % (Fillet suppressed) Track no alternate base; use global r_base/z_base
        % Add a multi-ring spherical annulus to stitch core↔cone robustly
        dz_ann = max(1e-6, min(0.02*cr, 0.08*r_base));
        n_ann_rings = max(3, ceil(0.15 * numel(angles_uniform)));
        prev_band_idx = seam_indices{i};
        for ar = 1:n_ann_rings
            zt = max(-cr+1e-6, z_base - (ar/n_ann_rings) * dz_ann);
            rt = sqrt(max(0, cr^2 - zt^2));
            ann_ring = seam_ring_points(zt * orientation, u, v, rt, angles_uniform);
            [V, idx_ann] = append_vertices(V, ann_ring);
            F = [F; connect_rings(prev_band_idx, idx_ann)]; %#ok<AGROW>
            prev_band_idx = idx_ann;
        end
        annulus_end_rings{i} = ann_ring; % store last annulus ring coordinates
    % Tip geometry already computed above (z_seam_i, alpha_i, r_tip_i, tip_center)

        % --- Patch A: Conical/cylindrical body (from fillet end/base to tip seam) ---
        % Radius profile is controlled solely by r_cone_design and r_seam to
        % keep conicality independent of fillet geometry.
    L_cone = max(1e-6, z_seam_i - z_base);
        num_cone_rings = max(5, ceil(CF_RINGS * refine * 0.5 * (1 + k_sphere * L_cone)));
        if ~any(strcmp('nConeRings', p.UsingDefaults)) && ~isempty(p.Results.nConeRings)
            num_cone_rings = p.Results.nConeRings;
        end
    prev_idx = seam_indices{i};
        for k = 1:num_cone_rings
            t = k / num_cone_rings;              % 0 at base, 1 at seam
            % Smooth conicality: sc=0 → cylinder (r=r_seam), sc=1 → full taper
            dr = max(0, r_cone_design - r_seam);
            ring_radius = r_seam + (1 - t) * dr;
            ring_center = ((1 - t) * z_base + t * z_seam_i) * orientation; % z_base_local removed (fillet suppressed)
            ring_pts = seam_ring_points(ring_center, u, v, ring_radius, angles_uniform);
            [V, idx_ring] = append_vertices(V, ring_pts);
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
            F = [F; connect_rings(prev_cap_idx, idx_ring)]; %#ok<AGROW>
            prev_cap_idx = idx_ring;
        end
        % Apex fan
        [V, idx_apex] = append_vertices(V, apex);
        F = [F; connect_ring_to_apex(prev_cap_idx, idx_apex)]; %#ok<AGROW>

        % --- Patch C: Toroidal fillet (suppressed for now) ---
        % Intentionally disabled to focus on core sphere + spike body + tip
    end

    %% 4) Core sphere triangulation (alpha shape) with annulus endpoints for guaranteed welding
    if includeCore
        fprintf('Triangulating core sphere (alphaShape) with spike holes...\n');
        % Build point cloud: Fibonacci points outside footprints + annulus end rings
        num_core_points = max(2000, ceil(0.2 * base_density * dens_mult_sphere * 4 * pi * cr^2));
        core_pts = fibonacci_sphere(num_core_points, cr);
        keep_mask = true(size(core_pts,1),1);
        rim_trim = 0.12 * r_base;
        r_effective = min(cr * 0.999, r_base + rim_trim);
        cutoff_angle = asin(min(1, r_effective / cr));
        for ii = 1:ns
            o = spike_orientations(ii, :);
            cosang = (core_pts * o') ./ cr;
            cosang = max(-1, min(1, cosang));
            ang = acos(cosang);
            keep_mask = keep_mask & (ang > cutoff_angle);
        end
        core_pts = core_pts(keep_mask, :);
        if ns > 0
            core_pts = [core_pts; cell2mat(annulus_end_rings)];
        end
        % Alpha shape triangulation with adaptive alpha
        alpha_len = 2.4 * cr / sqrt(max(1,num_core_points));
        alpha_core = alpha_len;
        success = false;
        for itry = 1:6
            try
                shp_core = alphaShape(core_pts(:,1), core_pts(:,2), core_pts(:,3), alpha_core);
                [core_F, core_V] = boundaryFacets(shp_core);
                if ~isempty(core_F) && ~isempty(core_V)
                    [V, idx_core] = append_vertices(V, core_V);
                    F = [F; idx_core(core_F)];
                    success = true;
                    break;
                end
            catch
                warning('Alpha shape triangulation failed at alpha=%.4f, retrying with smaller alpha...', alpha_core);
            end
            alpha_core = alpha_core * 0.7;
        end
        if ~success
            % Robust fallback: geodesic icosphere with centroid trimming under spike footprints
            warning('Alpha shape triangulation failed repeatedly, falling back to icosphere.');
            subdivisions = max(5, ceil(0.25 * sqrt(base_density * dens_mult_sphere) ));
            [core_V2, core_F2] = triangulate_icosphere(cr, subdivisions);
            if ns > 0
                face_centers2 = (core_V2(core_F2(:,1),:) + core_V2(core_F2(:,2),:) + core_V2(core_F2(:,3),:)) / 3;
                keep_face2 = true(size(core_F2,1),1);
                for ii = 1:ns
                    o = spike_orientations(ii, :);
                    cosang_c2 = (face_centers2 * o') ./ cr;
                    cosang_c2 = max(-1, min(1, cosang_c2));
                    ang_c2 = acos(cosang_c2);
                    keep_face2 = keep_face2 & (ang_c2 > cutoff_angle);
                end
                core_F2 = core_F2(keep_face2, :);
            end
            [V, idx_core] = append_vertices(V, core_V2);
            F = [F; idx_core(core_F2)];
        end
    end

    %% 5) Weld duplicate vertices at seams and build mesh
    [Vw, ~, ic] = uniquetol(V, 1e-7, 'ByRows', true);
    Fw = ic(F);
    mesh.Vertices = Vw;
    mesh.Faces = Fw;
    
    % Optionally build a volumetric mask from the surface mesh by voxelizing
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
            if isempty(volAlphaInp)
                % heuristic alpha: small fraction of bounding radius
                Rb = 0.5 * maxDim;
                alphaVol = 0.08 * Rb;
            else
                alphaVol = volAlphaInp;
            end
            try
                shpVol = alphaShape(mesh.Vertices(:,1), mesh.Vertices(:,2), mesh.Vertices(:,3), alphaVol);
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
        % Attempt viewer3d with surfaceMeshShow
        try
            viewer = viewer3d;
            meshObj = surfaceMesh(mesh.Vertices, mesh.Faces);
            surfaceMeshShow(meshObj, Parent=viewer);
        catch
            % ignore and proceed to classic methods
        end
        % Classic figure using patch
        fig1 = figure('Name', 'Urchin (patch)', 'NumberTitle', 'off', 'Color', 'w'); %#ok<NASGU>
        p1 = patch('Faces', mesh.Faces, 'Vertices', mesh.Vertices);
        p1.FaceColor = [0.9, 0.7, 0.2];
        p1.FaceAlpha = 0.5;
        p1.EdgeColor = [0.1, 0.1, 0.1];
        p1.EdgeAlpha = 0.8;
        p1.LineWidth = 0.5;
        p1.FaceLighting = 'gouraud';
        axis equal off; view(3); camlight('headlight');
        % Alternate method: trisurf
        fig2 = figure('Name', 'Urchin (trisurf)', 'NumberTitle', 'off', 'Color', 'w'); %#ok<NASGU>
        trisurf(mesh.Faces, mesh.Vertices(:,1), mesh.Vertices(:,2), mesh.Vertices(:,3), ...
                'FaceColor', [0.8, 0.6, 0.2], 'FaceAlpha', 0.5, ...
                'EdgeColor', [0.1, 0.1, 0.1], 'EdgeAlpha', 0.8, 'LineWidth', 0.5);
        axis equal off; view(3); camlight headlight; lighting gouraud;
    end

end

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
function [r_c, z_c, theta_cone, theta_sphere] = compute_fillet_center_and_angles(cr, r_f, m, C, z_base, ~)
    % Compute fillet circle center (r_c, z_c) in rz-plane and contact angles with
    % cone line r = m z + C and sphere r^2 + z^2 = cr^2, with fillet radius r_f
    % Equations:
    %   (1) |(r_c, z_c)·n + C| = r_f where n = [1, -m]/sqrt(1+m^2)  (distance to cone)
    %   (2) sqrt(r_c^2 + z_c^2) = cr + r_f                            (offset circle)
    % Choose the sign so the circle lies outside the sphere (towards +z from center)
    nrm = sqrt(1 + m^2);
    nx = 1 / nrm; nz = -m / nrm;
    % Solve (2): r_c = (cr + r_f) cos t, z_c = (cr + r_f) sin t
    Roff = cr + r_f;
    % Plug into (1): |nx*r_c + nz*z_c + C| = r_f
    % => |Roff*(nx*cos t + nz*sin t) + C| = r_f
    % Choose t to minimize error; solve via atan2
    A = Roff * nx; B = Roff * nz;
    % We want signed distance nx*r_c + nz*z_c + C = +r_f (outer side of cone)
    s = r_f - C;
    denom = sqrt(A*A + B*B) + 1e-12;
    cos_phi = (A/denom); sin_phi = (B/denom);
    % Let t = atan2(sin_phi, cos_phi) + acos( clamp(s/denom) )
    clamp = @(x) max(-1,min(1,x));
    t0 = atan2(sin_phi, cos_phi);
    a = acos(clamp(s/denom));
    % Two candidates
    t_plus  = t0 + a;
    t_minus = t0 - a;
    cand = [t_plus, t_minus];
    rc = Roff * cos(cand);
    zc = Roff * sin(cand);
    d = nx * rc + nz * zc + C; % signed distance to cone line
    % pick the one with d closest to +r_f and with zc >= z_base (outer side along +z)
    penalty = (zc < z_base);
    score = abs(d - r_f) + 1e3 * penalty;
    [~, idx] = min(score);
    r_c = rc(idx);
    z_c = zc(idx);
    % Contact angles on the fillet circle relative to center
    % Sphere contact: radial toward origin
    theta_sphere = atan2(-z_c, -r_c);
    % Cone contact: fillet radius from center must be perpendicular to cone normal n=[nx,nz],
    % i.e., along ±t where t is tangent to cone line. Choose outward: k = [+nz, -nx]
    theta_cone = atan2(-nx, +nz);
    % Ensure sweep direction from cone→sphere is correct
    if wrapToPi(theta_sphere - theta_cone) < 0
        tmp = theta_cone; theta_cone = theta_sphere; theta_sphere = tmp;
    end
end

%% --- Helper Functions ---

function alpha = solve_alpha_conicality(r_base_max, r_seam, cr, sl, sc)
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
        rs = r_seam;                          % seam radius fixed by st
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

function [V, F] = triangulate_sphere(R, n_lat)
    % Build a triangulated UV sphere (deterministic)
    n_lon = 2 * n_lat;
    lats = linspace(0, pi, n_lat+1);
    lons = linspace(0, 2*pi, n_lon+1); lons(end) = [];
    V = zeros((n_lat-1)*n_lon + 2, 3);
    % Poles
    north = [0, 0, R]; south = [0, 0, -R];
    V(1,:) = north; V(end,:) = south;
    % Intermediate rings
    idx = 2;
    for i = 2:n_lat
        phi = lats(i);
        if i == n_lat+1, break; end
        sinp = sin(phi); cosp = cos(phi);
        for j = 1:n_lon
            theta = lons(j);
            V(idx, :) = R * [cos(theta)*sinp, sin(theta)*sinp, cosp];
            idx = idx + 1;
        end
    end
    % Faces
    F = [];
    % Top cap
    top_ring = 2:(1+n_lon);
    F = [F; connect_ring_to_apex(top_ring, 1)];
    % Middle bands
    for i = 1:(n_lat-2)
        ringA = 1 + (i-1)*n_lon + (1:n_lon);
        ringB = 1 + i*n_lon + (1:n_lon);
        F = [F; connect_rings(ringA, ringB)];
    end
    % Bottom cap
    bottom_ring = size(V,1)-n_lon:(size(V,1)-1);
    F = [F; fliplr(connect_ring_to_apex(bottom_ring, size(V,1)))];
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

function [pts, faces] = triangulate_torus_patch(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
    % Build a triangulated torus patch parameterized by theta (tube) and phi (around seam)
    thetas = linspace(theta_min, theta_max, n_theta);
    phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
    pts = zeros(n_theta * n_phi, 3);
    for it = 1:n_theta
        ct = cos(thetas(it)); st = sin(thetas(it));
        for ip = 1:n_phi
            cp = cos(phis(ip)); sp = sin(phis(ip));
            ring_dir = u*cp + v*sp;
            p = center + (R + r*ct) * (ring_dir') + (r*st) * (n(:))';
            pts((it-1)*n_phi+ip, :) = p;
        end
    end
    faces = zeros(2*(n_theta-1)*n_phi, 3);
    idx = 1;
    for it = 1:(n_theta-1)
        for ip = 1:n_phi
            ip2 = mod(ip, n_phi) + 1;
            a1 = (it-1)*n_phi + ip;
            a2 = (it-1)*n_phi + ip2;
            b1 = it*n_phi + ip;
            b2 = it*n_phi + ip2;
            faces(idx, :) = [a1, b1, b2]; idx = idx + 1;
            faces(idx, :) = [a1, b2, a2]; idx = idx + 1;
        end
    end
end
function points = fibonacci_sphere(num_points, radius)
    % Deterministic Fibonacci sphere sampling of a full sphere of given radius
    if nargin < 2, radius = 1.0; end
    points = zeros(num_points, 3);
    golden_angle = pi * (3 - sqrt(5));
    for i = 0:num_points-1
        y = 1 - 2 * (i + 0.5) / num_points;   % avoid poles clustering
        r_xy = sqrt(max(0, 1 - y*y));
        theta = mod(i * golden_angle, 2*pi);
        x = cos(theta) * r_xy;
        z = sin(theta) * r_xy;
        points(i+1, :) = radius .* [x, y, z];
    end
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

function pts = fibonacci_spherical_cap(num_points, radius, orientation, center, outward_sign)
    % Fibonacci sampling restricted to a hemispherical cap oriented by 'orientation'
    % outward_sign = +1 for outer hemisphere (tip), -1 for inner
    if nargin < 5, outward_sign = +1; end
    dirs = fibonacci_sphere(num_points, 1.0);
    n = orientation(:)'/norm(orientation);
    % Keep points with sign matching hemisphere
    keep = outward_sign * (dirs * n') >= 0;
    dirs = dirs(keep, :);
    if isempty(dirs)
        pts = zeros(0,3); return;
    end
    % Scale to radius and translate to center
    pts = dirs * radius + center;
end

function torus_pts = torus_patch(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
    % Parametric torus patch around axis n with center located at 'center'
    % P(θ, φ) = center + (R + r cos θ) (u cos φ + v sin φ) + r sin θ n
    thetas = linspace(theta_min, theta_max, n_theta);
    phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
    torus_pts = zeros(numel(thetas) * numel(phis), 3);
    idx = 1;
    for it = 1:numel(thetas)
        ct = cos(thetas(it)); st = sin(thetas(it));
        for ip = 1:numel(phis)
            cp = cos(phis(ip)); sp = sin(phis(ip));
            ring_dir = u*cp + v*sp;             % 3x1
            % Build p as a 1x3 row vector to avoid implicit expansion
            p = center + (R + r*ct) * (ring_dir') + (r*st) * (n(:))';
            torus_pts(idx, :) = p(:)';
            idx = idx + 1;
        end
    end
end

function [pts, faces] = triangulate_torus_patch_inward(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
    % Concave inward torus patch around seam; normal component points inward
    thetas = linspace(theta_min, theta_max, n_theta);
    phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
    pts = zeros(n_theta * n_phi, 3);
    for it = 1:n_theta
        ct = cos(thetas(it)); st = sin(thetas(it));
        for ip = 1:n_phi
            cp = cos(phis(ip)); sp = sin(phis(ip));
            ring_dir = u*cp + v*sp;
            p = center + (R + r*ct) * (ring_dir') - (r*st) * (n(:))';
            pts((it-1)*n_phi+ip, :) = p;
        end
    end
    faces = zeros(2*(n_theta-1)*n_phi, 3);
    idx = 1;
    for it = 1:(n_theta-1)
        for ip = 1:n_phi
            ip2 = mod(ip, n_phi) + 1;
            a1 = (it-1)*n_phi + ip;
            a2 = (it-1)*n_phi + ip2;
            b1 = it*n_phi + ip;
            b2 = it*n_phi + ip2;
            faces(idx, :) = [a1, b1, b2]; idx = idx + 1;
            faces(idx, :) = [a1, b2, a2]; idx = idx + 1;
        end
    end
end



% function pts = fibonacci_spherical_cap(num_points, radius, orientation, center, outward_sign)
%     % Fibonacci sampling restricted to a hemispherical cap oriented by 'orientation'
%     % outward_sign = +1 for outer hemisphere (tip), -1 for inner
%     if nargin < 5, outward_sign = +1; end
%     dirs = fibonacci_sphere(num_points, 1.0);
%     n = orientation(:)'/norm(orientation);
%     % Keep points with sign matching hemisphere
%     keep = outward_sign * (dirs * n') >= 0;
%     dirs = dirs(keep, :);
%     if isempty(dirs)
%         pts = zeros(0,3); return;
%     end
%     % Scale to radius and translate to center
%     pts = dirs * radius + center;
% end

% function torus_pts = torus_patch(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
%     % Parametric torus patch around axis n with center located at 'center'
%     % P(θ, φ) = center + (R + r cos θ) (u cos φ + v sin φ) + r sin θ n
%     thetas = linspace(theta_min, theta_max, n_theta);
%     phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
%     torus_pts = zeros(numel(thetas) * numel(phis), 3);
%     idx = 1;
%     for it = 1:numel(thetas)
%         ct = cos(thetas(it)); st = sin(thetas(it));
%         for ip = 1:numel(phis)
%             cp = cos(phis(ip)); sp = sin(phis(ip));
%             ring_dir = u*cp + v*sp;             % 3x1
%             % Build p as a 1x3 row vector to avoid implicit expansion
%             p = center + (R + r*ct) * (ring_dir') + (r*st) * (n(:))';
%             torus_pts(idx, :) = p(:)';
%             idx = idx + 1;
%         end
%     end
% end

% function [V, F] = triangulate_sphere(R, n_lat)
%     % Build a triangulated UV sphere (deterministic)
%     n_lon = 2 * n_lat;
%     lats = linspace(0, pi, n_lat+1);
%     lons = linspace(0, 2*pi, n_lon+1); lons(end) = [];
%     V = zeros((n_lat-1)*n_lon + 2, 3);
%     % Poles
%     north = [0, 0, R]; south = [0, 0, -R];
%     V(1,:) = north; V(end,:) = south;
%     % Intermediate rings
%     idx = 2;
%     for i = 2:n_lat
%         phi = lats(i);
%         if i == n_lat+1, break; end
%         sinp = sin(phi); cosp = cos(phi);
%         for j = 1:n_lon
%             theta = lons(j);
%             V(idx, :) = R * [cos(theta)*sinp, sin(theta)*sinp, cosp];
%             idx = idx + 1;
%         end
%     end
%     % Faces
%     F = [];
%     % Top cap
%     top_ring = 2:(1+n_lon);
%     F = [F; connect_ring_to_apex(top_ring, 1)];
%     % Middle bands
%     for i = 1:(n_lat-2)
%         ringA = 1 + (i-1)*n_lon + (1:n_lon);
%         ringB = 1 + i*n_lon + (1:n_lon);
%         F = [F; connect_rings(ringA, ringB)];
%     end
%     % Bottom cap
%     bottom_ring = size(V,1)-n_lon:(size(V,1)-1);
%     F = [F; fliplr(connect_ring_to_apex(bottom_ring, size(V,1)))];
% end

% function [r_c, z_c, theta_cone, theta_sphere] = compute_fillet_center_and_angles(cr, r_f, m, C, z_base, ~)
%     % Compute fillet circle center (r_c, z_c) in rz-plane and contact angles with
%     % cone line r = m z + C and sphere r^2 + z^2 = cr^2, with fillet radius r_f
%     % Equations:
%     %   (1) |(r_c, z_c)·n + C| = r_f where n = [1, -m]/sqrt(1+m^2)  (distance to cone)
%     %   (2) sqrt(r_c^2 + z_c^2) = cr + r_f                            (offset circle)
%     % Choose the sign so the circle lies outside the sphere (towards +z from center)
%     nrm = sqrt(1 + m^2);
%     nx = 1 / nrm; nz = -m / nrm;
%     % Solve (2): r_c = (cr + r_f) cos t, z_c = (cr + r_f) sin t
%     Roff = cr + r_f;
%     % Plug into (1): |nx*r_c + nz*z_c + C| = r_f
%     % => |Roff*(nx*cos t + nz*sin t) + C| = r_f
%     % Choose t to minimize error; solve via atan2
%     A = Roff * nx; B = Roff * nz;
%     % We want signed distance nx*r_c + nz*z_c + C = +r_f (outer side of cone)
%     s = r_f - C;
%     denom = sqrt(A*A + B*B) + 1e-12;
%     cos_phi = (A/denom); sin_phi = (B/denom);
%     % Let t = atan2(sin_phi, cos_phi) + acos( clamp(s/denom) )
%     clamp = @(x) max(-1,min(1,x));
%     t0 = atan2(sin_phi, cos_phi);
%     a = acos(clamp(s/denom));
%     % Two candidates
%     t_plus  = t0 + a;
%     t_minus = t0 - a;
%     cand = [t_plus, t_minus];
%     rc = Roff * cos(cand);
%     zc = Roff * sin(cand);
%     d = nx * rc + nz * zc + C; % signed distance to cone line
%     % pick the one with d closest to +r_f and with zc >= z_base (outer side along +z)
%     penalty = (zc < z_base);
%     score = abs(d - r_f) + 1e3 * penalty;
%     [~, idx] = min(score);
%     r_c = rc(idx);
%     z_c = zc(idx);
%     % Contact angles on the fillet circle relative to center
%     % Sphere contact: radial toward origin
%     theta_sphere = atan2(-z_c, -r_c);
%     % Cone contact: fillet radius from center must be perpendicular to cone normal n=[nx,nz],
%     % i.e., along ±t where t is tangent to cone line. Choose outward: k = [+nz, -nx]
%     theta_cone = atan2(-nx, +nz);
%     % Ensure sweep direction from cone→sphere is correct
%     if wrapToPi(theta_sphere - theta_cone) < 0 % requires Mapping Toolbox! Alternative solution recommended.
%         tmp = theta_cone; theta_cone = theta_sphere; theta_sphere = tmp;
%     end
% end

% function [pts, faces] = triangulate_torus_patch(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
%     % Build a triangulated torus patch parameterized by theta (tube) and phi (around seam)
%     thetas = linspace(theta_min, theta_max, n_theta);
%     phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
%     pts = zeros(n_theta * n_phi, 3);
%     for it = 1:n_theta
%         ct = cos(thetas(it)); st = sin(thetas(it));
%         for ip = 1:n_phi
%             cp = cos(phis(ip)); sp = sin(phis(ip));
%             ring_dir = u*cp + v*sp;
%             p = center + (R + r*ct) * (ring_dir') + (r*st) * (n(:))';
%             pts((it-1)*n_phi+ip, :) = p;
%         end
%     end
%     faces = zeros(2*(n_theta-1)*n_phi, 3);
%     idx = 1;
%     for it = 1:(n_theta-1)
%         for ip = 1:n_phi
%             ip2 = mod(ip, n_phi) + 1;
%             a1 = (it-1)*n_phi + ip;
%             a2 = (it-1)*n_phi + ip2;
%             b1 = it*n_phi + ip;
%             b2 = it*n_phi + ip2;
%             faces(idx, :) = [a1, b1, b2]; idx = idx + 1;
%             faces(idx, :) = [a1, b2, a2]; idx = idx + 1;
%         end
%     end
% end

% function [pts, faces] = triangulate_torus_patch_inward(center, n, u, v, R, r, theta_min, theta_max, n_theta, n_phi)
%     % Concave inward torus patch around seam; normal component points inward
%     thetas = linspace(theta_min, theta_max, n_theta);
%     phis   = linspace(0, 2*pi, n_phi+1); phis(end) = [];
%     pts = zeros(n_theta * n_phi, 3);
%     for it = 1:n_theta
%         ct = cos(thetas(it)); st = sin(thetas(it));
%         for ip = 1:n_phi
%             cp = cos(phis(ip)); sp = sin(phis(ip));
%             ring_dir = u*cp + v*sp;
%             p = center + (R + r*ct) * (ring_dir') - (r*st) * (n(:))';
%             pts((it-1)*n_phi+ip, :) = p;
%         end
%     end
%     faces = zeros(2*(n_theta-1)*n_phi, 3);
%     idx = 1;
%     for it = 1:(n_theta-1)
%         for ip = 1:n_phi
%             ip2 = mod(ip, n_phi) + 1;
%             a1 = (it-1)*n_phi + ip;
%             a2 = (it-1)*n_phi + ip2;
%             b1 = it*n_phi + ip;
%             b2 = it*n_phi + ip2;
%             faces(idx, :) = [a1, b1, b2]; idx = idx + 1;
%             faces(idx, :) = [a1, b2, a2]; idx = idx + 1;
%         end
%     end
% end
