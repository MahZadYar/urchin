function rBaseLimit = tangentBaseLimit(coreRadius, tipSphereRadius, spikeLength)
%TANGENTBASELIMIT Maximum cone base radius that preserves spherical tip tangency.
%   rBaseLimit = TANGENTBASELIMIT(coreRadius, tipSphereRadius, spikeLength)
%   returns the largest allowable base radius on the core sphere such that a
%   straight conical flank can remain tangent to the spherical tip. The
%   result is bounded to the geometric limits of the core sphere radius and
%   degrades gracefully when a closed-form tangency is impossible (e.g.,
%   degenerate tip or extremely short spikes).

    coreRadius = max(coreRadius, 1e-9);
    spikeLength = max(spikeLength, 0);
    tipSphereRadius = max(tipSphereRadius, 0);

    if tipSphereRadius <= 0
        % No spherical tip -- tangency constraint does not apply.
        rBaseLimit = min(0.999 * coreRadius, coreRadius);
        return;
    end

    % Axial location of the spherical tip centre.
    zTipCenter = coreRadius + spikeLength - tipSphereRadius;
    if zTipCenter <= 0
        % Tip centre is at or inside the core origin; tangency collapses.
        rBaseLimit = 0;
        return;
    end

    % Solve for the core intersection height using the closed-form tangency
    % condition derived from the quadratic cone-sphere system.
    numerator = zTipCenter^2 + coreRadius^2 - tipSphereRadius^2;
    denominator = 2 * zTipCenter;

    if denominator <= 0
        rBaseLimit = 0;
        return;
    end

    zBaseCandidate = numerator / denominator;
    zBaseCandidate = min(max(zBaseCandidate, 0), coreRadius);

    rSquared = max(coreRadius^2 - zBaseCandidate^2, 0);
    rBaseLimit = sqrt(rSquared);
    rBaseLimit = min(rBaseLimit, 0.999 * coreRadius);

    if ~isfinite(rBaseLimit)
        rBaseLimit = 0;
    end
end
