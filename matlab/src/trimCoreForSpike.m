function [coreLoop, coreFCurrent, coreMask] = trimCoreForSpike(V, coreFCurrent, coreMask, ...
                                                                     orientation, cosCutoffI, rBase, ...
                                                                     seamIndices, i)
    % TRIMCOREFORSPIKE - Extract core vertices that occlude spike placement
    %
    % Removes core sphere vertices and faces that would collide with spike,
    % and identifies the trim boundary loop for stitching.
    %
    % SYNTAX:
    %   [coreLoop, coreFCurrent, coreMask] = trim_core_for_spike(V, coreFCurrent, coreMask, ...)
    %
    % INPUTS:
    %   V           - Vertex array (Nx3)
    %   coreFCurrent - Active core face indices (Mx3)
    %   coreMask    - Boolean mask of active core vertices
    %   orientation - Unit spike direction vector (1x3)
    %   cosCutoffI  - Cosine of max angle before occlusion
    %   rBase       - Spike base radius
    %   seamIndices - Cell array of seam boundary rings
    %   i           - Current spike index
    %
    % OUTPUTS:
    %   coreLoop    - Boundary vertex indices after trim (or empty)
    %   coreFCurrent - Updated active core faces
    %   coreMask    - Updated active vertex mask
    %
    % Performance: O(num_core_vertices * num_previous_spikes) for candidate search
    
    coreLoop = [];
    
    % Find core vertices that occlude this spike
    coreIdx = find(coreMask);
    if isempty(coreIdx)
        coreDirsFull = zeros(0,3);
    else
        dirs = V(coreIdx,:);
        norms = vecnorm(dirs, 2, 2);
        norms(norms < 1e-12) = 1;
        coreDirsFull = dirs ./ norms;
    end

    % Compute dot products to identify occlusion region
    coreDot = coreDirsFull * orientation';
    newRemoveMask = coreDot > cosCutoffI + 1e-12;
    
    % Fallback: always remove at least one vertex (nearest to spike axis)
    if ~any(newRemoveMask) && ~isempty(coreIdx)
        [~, fallbackPos] = max(coreDot);
        if ~isnan(fallbackPos) && fallbackPos >= 1 && fallbackPos <= numel(newRemoveMask)
            newRemoveMask(fallbackPos) = true;
        end
    end
    
    if any(newRemoveMask)
        % Mark vertices for removal
        vertsRemoveGlobal = coreIdx(newRemoveMask);
        coreMask(vertsRemoveGlobal) = false;
        
        % Remove faces containing removed vertices
        faceRemoveMask = any(ismember(coreFCurrent, vertsRemoveGlobal), 2);
        removedFaces = coreFCurrent(faceRemoveMask, :);
        coreFCurrent(faceRemoveMask, :) = [];

        % Identify trim boundary from removed faces
        vertsFromRemovedFaces = unique(removedFaces(:));
        coreLoop = setdiff(vertsFromRemovedFaces, vertsRemoveGlobal);

        % Find best candidate loop for stitching (prefer collinear with spike axis)
        candidates = {};
        if numel(coreLoop) >= 3
            candidates{end+1} = coreLoop(:)'; %#ok<AGROW>
        end
        
        % Check previous spike seams for existing loop candidates
        if ~isempty(vertsRemoveGlobal)
            for j = 1:i-1
                ring = seamIndices{j};
                if numel(ring) < 3, continue; end
                if all(ismember(ring, vertsRemoveGlobal))
                    candidates{end+1} = ring(:)'; %#ok<AGROW>
                end
            end
        end
        
        % Select best candidate (highest dot product with spike direction)
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
