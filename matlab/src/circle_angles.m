function angles = circle_angles(segments)
%CIRCLE_ANGLES Evenly spaced angular samples around a circle.
%   angles = CIRCLE_ANGLES(segments) returns SEGMENTS angles in radians
%   covering [0, 2*pi) without duplicating the endpoint. When segments < 3
%   the function returns an empty array so that rings collapse to a point.

    if segments < 3
        angles = [];
        return;
    end
    angles = linspace(0, 2*pi, segments + 1);
    angles(end) = [];
end
