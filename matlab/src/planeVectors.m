function [u, v] = planeVectors(normal_vec)
%PLANEVECTORS Generate orthonormal basis vectors perpendicular to normal.
%   [u, v] = PLANEVECTORS(normal_vec) returns two orthonormal vectors that
%   span the plane orthogonal to normal_vec.

    n = normal_vec(:) / norm(normal_vec);
    if abs(n(1)) > 0.9
        other = [0; 1; 0];
    else
        other = [1; 0; 0];
    end
    u = cross(n, other);
    u = u / max(1e-12, norm(u));
    v = cross(n, u);
    v = v / max(1e-12, norm(v));
end
