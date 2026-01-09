function orientations = refineSpikeOrientations(orientations, thresholdDeg)
%REFINESPIKEORIENTATIONS Relax spike directions on the unit sphere.
%   orientations = refineSpikeOrientations(orientations, thresholdDeg)
%   perturbs the input directions (each row is a unit vector) using a
%   Coulomb-like repulsion model until the nearest-neighbour angular
%   separation meets the requested threshold (degrees) or a fixed maximum
%   iteration count is reached.
%
%   Inputs:
%       orientations  N-by-3 unit vectors.
%       thresholdDeg  Minimum allowable angle in degrees (>=0).
%
%   Output:
%       orientations  Refined unit vectors.
%
%   The routine is intentionally conservative: it exits early when the
%   requested threshold is satisfied or when the step updates no longer
%   increase the minimum separation.
%
% See also: acos, dot

    if nargin < 2 || ~isfinite(thresholdDeg)
        thresholdDeg = 0;
    end

    if isempty(orientations) || size(orientations, 1) < 2 || thresholdDeg <= 0
        orientations = normalizeRows(orientations);
        return;
    end

    pts = normalizeRows(orientations);
    n = size(pts, 1);
    maxIterations = 1000;
    % The solver now uses a damped momentum update (similar to gradient-descent
    % optimisers) instead of raw force application. A small learning rate keeps
    % the system responsive, while the momentum accumulator reuses a fraction
    % of the previous step to avoid oscillations when neighbouring spikes react
    % to each other simultaneously.
    step = 0.005;
    momentum = 0.8;
    velocity = zeros(n, 3);
    clampCos = @(x) min(1, max(-1, x));
    thresholdRad = max(0, thresholdDeg) * pi / 180;
    [initialMinAngleRad, ~] = pairwiseAngleExtents(pts);
    initialMinAngleDeg = initialMinAngleRad * 180 / pi;
    fprintf('Orientation refinement: min-\angle %.3f°, change threshold %.3g° (max %d iterations)\n', initialMinAngleDeg, thresholdDeg, maxIterations);
    % Set URCHIN_REFINE_VIS=1 in the environment before launching MATLAB to
    % visualise the relaxation (otherwise this stays false to avoid GUI work).
    visRefinement = ~isempty(getenv('URCHIN_REFINE_VIS'));
    visHandles = struct('enabled', visRefinement, 'fig', [], 'ax', [], ...
        'scatter', [], 'stepLines', gobjects(0), 'cmap', parula(128), 'maxIter', maxIterations);
    if visRefinement
        visHandles = initRefinementVisualization(visHandles, pts);
    end
    iterCount = 0;
    maxDelta = Inf;
    for iter = 1:maxIterations
        prevPts = pts;
        % Vectorised Coulomb-like repulsion (avoids per-spike loop)
        % diffs: n x n x 3 tensor of pairwise deltas
        diffs = reshape(pts, [n, 1, 3]) - reshape(pts, [1, n, 3]);
        distSq = sum(diffs.^2, 3);
        distSq(1:n+1:end) = Inf;           % ignore self-distances
        invDist3 = distSq.^(-1.5);
        invDist3(~isfinite(invDist3)) = 0;
        forces = squeeze(sum(diffs .* invDist3, 2));

        % Remove radial component to keep pure tangential motion
        radialComp = sum(forces .* pts, 2);
        forces = forces - radialComp .* pts;
        
        % Update velocity with momentum
        velocity = momentum * velocity + step * forces;
        
        % Remove radial component from velocity as well to stay on surface
        velRadial = sum(velocity .* pts, 2);
        velocity = velocity - velRadial .* pts;

        [minAngleRad, ~] = pairwiseAngleExtents(pts);
        limitRad = minAngleRad / 5;
        
        % Check if velocity is too large
        velMags = sqrt(sum(velocity.^2, 2));
        maxVel = max(velMags);
        
        if limitRad > 0 && isfinite(limitRad) && maxVel > limitRad
            scaleFactor = limitRad / maxVel;
            velocity = velocity * scaleFactor;
        end
        
        candidate = applyStep(prevPts, velocity, 1.0);
        [maxDelta, candidate] = measureDelta(prevPts, candidate, clampCos);

        pts = candidate;
        iterCount = iter;

        visHandles = updateRefinementVisualization(visHandles, prevPts, pts, iter);

        if maxDelta < thresholdRad
            break;
        end
    end

    [finalMinAngleRad, ~] = pairwiseAngleExtents(pts);
    finalMinAngleDeg = finalMinAngleRad * 180 / pi;
    finalDeltaDeg = maxDelta * 180 / pi;
    thresholdMet = finalDeltaDeg < thresholdDeg - 1e-6;
    fprintf('  iterations %d, final min %.3f°, max Δ %.4f° (threshold %.3g°) -> %s\n', ...
        iterCount, finalMinAngleDeg, finalDeltaDeg, thresholdDeg, mat2str(thresholdMet));

    orientations = pts;
end

function ptsNew = applyStep(prevPts, forces, scaledStep)
%APPLYSTEP Apply scaled step and renormalize orientations.
    ptsNew = prevPts + scaledStep * forces;
    ptsNew = normalizeRows(ptsNew);
end

function [maxDelta, ptsNew] = measureDelta(prevPts, ptsCurrent, clampCos)
%MEASUREDELTA Compute angular change per spike between iterations.
    deltaCos = sum(prevPts .* ptsCurrent, 2);
    deltaCos = clampCos(deltaCos);
    deltaAngles = acos(deltaCos);
    maxDelta = max(deltaAngles);
    ptsNew = ptsCurrent;
end

function vis = initRefinementVisualization(vis, pts)
%INITREFINEMENTVISUALIZATION Prepare 3D plot of spike orientations.
    vis.enabled = true;
    vis.fig = figure('Name','Spike Orientation Refinement','NumberTitle','off');
    vis.ax = axes('Parent', vis.fig); %#ok<LAXES>
    hold(vis.ax, 'on');
    [sx, sy, sz] = sphere(50);
    surf(vis.ax, sx, sy, sz, 'FaceAlpha', 0.5, 'EdgeColor', [0.1 0.1 0.1], 'FaceColor', [0 0 0]);
    vis.scatter = scatter3(vis.ax, pts(:,1), pts(:,2), pts(:,3), 36, 'filled', 'MarkerFaceColor', [0.1 0.45 0.9]);
    vis.stepLines = gobjects(0);
    axis(vis.ax, 'equal');
    grid(vis.ax, 'on');
    xlabel(vis.ax, 'X'); ylabel(vis.ax, 'Y'); zlabel(vis.ax, 'Z');
    title(vis.ax, 'Spike Orientation Refinement');
    view(vis.ax, 135, 30);
end

function vis = updateRefinementVisualization(vis, prevPts, pts, iter)
%UPDATEREFINEMENTVISUALIZATION Refresh scatter plot and draw displacement traces.
    if isempty(vis) || ~isfield(vis, 'enabled') || ~vis.enabled
        return;
    end
    if isempty(vis.scatter) || ~isgraphics(vis.scatter)
        return;
    end
    set(vis.scatter, 'XData', pts(:,1), 'YData', pts(:,2), 'ZData', pts(:,3));
    if nargin < 4
        iter = numel(vis.stepLines) + 1;
    end
    cmap = vis.cmap;
    idx = max(1, min(size(cmap,1), round((iter-1)/(max(1, vis.maxIter-1)) * (size(cmap,1)-1)) + 1));
    lineColor = cmap(idx, :);
    newLines = plot3(vis.ax, [prevPts(:,1) pts(:,1)]', [prevPts(:,2) pts(:,2)]', [prevPts(:,3) pts(:,3)]', ...
        '-', 'Color', lineColor, 'LineWidth', 0.5);
    vis.stepLines = [vis.stepLines; newLines(:)];
    drawnow limitrate;
end

function v = normalizeRows(v)
%NORMALIZEROWS Project vectors back to unit sphere.
    if isempty(v)
        return;
    end
    norms = sqrt(sum(v.^2, 2));
    norms(norms == 0) = 1;
    v = v ./ norms;
end

function [minAngleRad, maxAngleRad] = pairwiseAngleExtents(pts)
%PAIRWISEANGLEEXTENTS Smallest and largest inter-spike angles.
    n = size(pts,1);
    if n < 2
        minAngleRad = 0;
        maxAngleRad = pi;
        return;
    end
    dots = pts * pts.';
    dots = max(-1, min(1, dots));

    dotsForNearest = dots;
    dotsForNearest(1:n+1:end) = -Inf;
    cosClosest = max(dotsForNearest, [], 2);
    cosClosest(~isfinite(cosClosest)) = -1;

    dotsForFarthest = dots;
    dotsForFarthest(1:n+1:end) = Inf;
    cosFarthest = min(dotsForFarthest, [], 2);
    cosFarthest(~isfinite(cosFarthest)) = 1;

    minAngleRad = acos(max(-1, min(1, max(cosClosest))));
    maxAngleRad = acos(max(-1, min(1, min(cosFarthest))));
end
