function segments = ring_segment_count(radius, spacing, spacingFactor)
%RING_SEGMENT_COUNT Number of segments required to sample a circular ring.
%   segments = RING_SEGMENT_COUNT(radius, spacing, spacingFactor) returns the
%   minimum number of evenly spaced points needed to respect the requested
%   linear spacing. Radii below half the spacing collapse to a single point.

    if radius <= spacing * 0.5
        segments = 0;
        return;
    end
    circumference = 2 * pi * radius;
    segments = max(3, ceil(circumference / max(spacing * spacingFactor, 1e-9)));
end
