function ringPts = generateCircularRing(ringCenter, u, v, ringRadius, minSpacing, spacingFactor)
    % GENERATECIRCULARRING - Consolidated ring point generation
    %
    % Generates a uniform circular ring with adaptive segment count based on radius.
    % Handles both regular rings and minimal rings (3 points).
    %
    % SYNTAX:
    %   ringPts = generateCircularRing(ringCenter, u, v, ringRadius, minSpacing, spacingFactor)
    %
    % INPUTS:
    %   ringCenter  - Ring center position (1x3)
    %   u, v        - Orthonormal frame vectors (1x3 each)
    %   ringRadius  - Ring radius in local frame
    %   minSpacing  - Minimum spacing for segment count
    %   spacingFactor - Spacing factor (typically 1.0 for body, <1.0 for caps)
    %
    % OUTPUTS:
    %   ringPts     - Ring vertex positions (Nx3)
    %
    % Performance: ~1-2ms per ring
    
    % Determine segment count
    segments = ringSegmentCount(ringRadius, minSpacing, spacingFactor);
    segments = max(3, segments);  % Minimum 3 vertices per ring
    
    % Generate angles around the ring
    angles = circleAngles(segments);
    
    % Generate ring points in 3D space
    ringPts = seamRingPoints(ringCenter, u, v, ringRadius, angles);
end
