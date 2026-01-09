function urchin = generateVolumeRepresentation(urchin, genVolume, volAdaptive, volPadding, volRes, ...
                                                 volAlphaInp, volDxMax, volDxMin, volBlockSz, volCriterion)
    % GENERATEVOLUMEREPRESENTATION - Generate volumetric representation of urchin mesh
    %
    % Creates either a dense voxel mask or sparse adaptive octree representation
    % of the surface mesh geometry.
    %
    % SYNTAX:
    %   urchin = generate_volume_representation(urchin, genVolume, volAdaptive, ...)
    %
    % INPUTS:
    %   urchin      - Struct with .Vertices and .Faces
    %   genVolume   - If true, compute volume
    %   volAdaptive - If true, use sparse octree; else dense voxel grid
    %   volPadding  - Padding factor for bounding box (default 0.1)
    %   volRes      - Dense grid resolution (default 128)
    %   volAlphaInp - Alpha shape radius (default auto)
    %   volDxMax    - Max octree leaf size (default auto)
    %   volDxMin    - Min octree leaf size (default auto)
    %   volBlockSz  - Octree block size (default 2)
    %   volCriterion - Octree refinement criterion (default 'surface')
    %
    % OUTPUTS:
    %   urchin      - Updated struct with VolumeMask or VolumeOctree fields
    %
    % Performance: Dense mode ~2-10s for 128³ grid; Sparse ~1-5s depending on geometry
    
    if ~genVolume
        return;
    end
    
    % Compute bounding box
    Vmin = min(urchin.Vertices, [], 1);
    Vmax = max(urchin.Vertices, [], 1);
    span = Vmax - Vmin;
    Vmin = Vmin - volPadding * span;
    Vmax = Vmax + volPadding * span;
    maxDim = max(span);
    
    if ~volAdaptive
        % Dense voxel grid approach
        urchin = generate_dense_volume(urchin, Vmin, Vmax, maxDim, volRes, volAlphaInp);
    else
        % Sparse adaptive octree approach
        urchin = generate_sparse_volume(urchin, Vmin, Vmax, maxDim, volDxMax, volDxMin, volBlockSz, volCriterion, volAlphaInp);
    end
end

function urchin = generate_dense_volume(urchin, Vmin, Vmax, maxDim, volRes, volAlphaInp)
    % Generate dense voxel mask using inpolyhedron or alphaShape fallback
    
    fprintf("Generating dense volumetric mask via voxelization...\n");
    
    % Build cubic voxels
    dx = maxDim / volRes;
    xs = Vmin(1):dx:Vmax(1);
    ys = Vmin(2):dx:Vmax(2);
    zs = Vmin(3):dx:Vmax(3);
    [XX, YY, ZZ] = ndgrid(xs, ys, zs);
    P = [XX(:), YY(:), ZZ(:)];
    inside = false(size(P, 1), 1);
    
    % Try inpolyhedron first (faster, more accurate)
    try
        inside = inpolyhedron(urchin.Faces, urchin.Vertices, P);
    catch
        % Fallback to alphaShape
        warning('inpolyhedron failed; falling back to alphaShape.');
        if isempty(volAlphaInp)
            Rb = 0.5 * maxDim;
            alphaVol = 0.08 * Rb;
        else
            alphaVol = volAlphaInp;
        end
        try
            shpVol = alphaShape(urchin.Vertices, alphaVol);
            inside = inShape(shpVol, P(:, 1), P(:, 2), P(:, 3));
        catch
            inside = false(size(P, 1), 1);
        end
    end
    
    % Store volumetric data
    VolumeMask = reshape(logical(inside), size(XX));
    urchin.VolumeMask = VolumeMask;
    urchin.VoxelGrid.X = xs;
    urchin.VoxelGrid.Y = ys;
    urchin.VoxelGrid.Z = zs;
    urchin.VoxelSize = dx;
    urchin.Bounds = [Vmin; Vmax];
end

function urchin = generate_sparse_volume(urchin, Vmin, Vmax, maxDim, volDxMax, volDxMin, volBlockSz, volCriterion, volAlphaInp)
    % Generate sparse adaptive octree representation
    
    fprintf("Generating adaptive sparse octree volumetric representation...\n");
    
    % Set default resolutions if not provided
    if isempty(volDxMax)
        volDxMax = maxDim / 128;
    end
    if isempty(volDxMin)
        volDxMin = volDxMax / 4;
    end
    
    % Build octree
    insideFn = makeInsideTester(urchin, volAlphaInp);
    leaves = buildAdaptiveOctree(Vmin, Vmax, volDxMax, volDxMin, volBlockSz, insideFn, volCriterion);
    
    % Store octree data
    urchin.VolumeOctree.Leaves = leaves;
    urchin.VolumeOctree.BlockSize = volBlockSz;
    urchin.VolumeOctree.DxMin = volDxMin;
    urchin.VolumeOctree.DxMax = volDxMax;
    urchin.VolumeOctree.Bounds = [Vmin; Vmax];
end
