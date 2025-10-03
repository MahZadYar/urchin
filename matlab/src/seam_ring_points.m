function ring_pts = seam_ring_points(center, u, v, radius, angles)
%SEAM_RING_POINTS Deterministic points on a seam ring given basis vectors.
%   ring_pts = SEAM_RING_POINTS(center, u, v, radius, angles) returns the
%   Cartesian coordinates of points located on a circle defined by the
%   orthonormal basis {u, v}.

    ca = cos(angles);
    sa = sin(angles);
    ring_pts = [center(1) + radius * (u(1) * ca + v(1) * sa)', ...
                center(2) + radius * (u(2) * ca + v(2) * sa)', ...
                center(3) + radius * (u(3) * ca + v(3) * sa)'];
end
