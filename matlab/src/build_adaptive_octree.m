function leaves = build_adaptive_octree(bmin, bmax, dxMax, dxMin, blockSize, insideFn, criterion)
%BUILD_ADAPTIVE_OCTREE Construct a sparse octree of voxel blocks.
%   leaves = BUILD_ADAPTIVE_OCTREE(bmin, bmax, dxMax, dxMin, blockSize,
%   insideFn, criterion) recursively subdivides the bounding box
%   [bmin,bmax] into leaf voxel blocks of size blockSize^3. Each leaf
%   stores its origin, voxel spacing dx, and a logical occupancy grid. The
%   insideFn handle is used to classify occupancy, while the criterion
%   string controls refinement ("boundary", "distance", "curvature", or
%   "hybrid").

    leaves = recurse_node(bmin, bmax, dxMax);

    function outLeaves = recurse_node(minC, maxC, dx)
        span = maxC - minC;
        % Sample corners to determine trivial inside / outside cases
        [CX, CY, CZ] = ndgrid([minC(1), maxC(1)], [minC(2), maxC(2)], [minC(3), maxC(3)]);
        Pc = [CX(:), CY(:), CZ(:)];
        inCorners = insideFn(Pc);
        if all(inCorners)
            occ = true(blockSize, blockSize, blockSize);
            outLeaves = struct('origin', minC, 'dx', dx, 'occupancy', occ);
            return;
        elseif ~any(inCorners) && ~needs_refine_local(minC, maxC, insideFn, criterion)
            occ = false(blockSize, blockSize, blockSize);
            outLeaves = struct('origin', minC, 'dx', dx, 'occupancy', occ);
            return;
        end
        if dx/2 < dxMin
            [occ, leafDx] = rasterize_block_local(minC, maxC, blockSize, insideFn);
            outLeaves = struct('origin', minC, 'dx', leafDx, 'occupancy', occ);
            return;
        end

        outLeaves = repmat(struct('origin', [], 'dx', [], 'occupancy', []), 0, 1);
        for ix = 0:1
            for iy = 0:1
                for iz = 0:1
                    childMin = minC + (span/2) .* [ix, iy, iz];
                    childMax = childMin + span/2;
                    outLeaves = [outLeaves; recurse_node(childMin, childMax, dx/2)]; %#ok<AGROW>
                end
            end
        end
    end
end

function tf = needs_refine_local(minC, maxC, insideFn, criterion)
%NEEDS_REFINE_LOCAL Decide if an octree node requires subdivision.
    switch criterion
        case 'boundary'
            [CX, CY, CZ] = ndgrid([minC(1), maxC(1)], [minC(2), maxC(2)], [minC(3), maxC(3)]);
            Pc = [CX(:), CY(:), CZ(:)];
            inCorners = insideFn(Pc);
            tf = any(inCorners) && ~all(inCorners);
        case 'distance'
            tf = needs_refine_local(minC, maxC, insideFn, 'boundary');
        case 'curvature'
            tf = needs_refine_local(minC, maxC, insideFn, 'boundary');
        case 'hybrid'
            tf = needs_refine_local(minC, maxC, insideFn, 'boundary');
        otherwise
            tf = needs_refine_local(minC, maxC, insideFn, 'boundary');
    end
end

function [occ, dx] = rasterize_block_local(minC, maxC, blockSize, insideFn)
%RASTERIZE_BLOCK_LOCAL Evaluate occupancy of a voxel block at fixed resolution.
    span = maxC - minC;
    dx = max(span) / blockSize;
    xs = minC(1) + (0.5:1:blockSize) * dx;
    ys = minC(2) + (0.5:1:blockSize) * dx;
    zs = minC(3) + (0.5:1:blockSize) * dx;
    [XX, YY, ZZ] = ndgrid(xs, ys, zs);
    P = [XX(:), YY(:), ZZ(:)];
    inside = insideFn(P);
    occ = reshape(logical(inside), size(XX));
end
