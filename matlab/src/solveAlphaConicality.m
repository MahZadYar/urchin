function alpha = solveAlphaConicality(r_base, r_tip, z_base, z_seam)
%SOLVEALPHACONICALITY Compute the cone half-angle ensuring G1 continuity.
%   alpha = SOLVEALPHACONICALITY(r_base, r_tip, z_base, z_seam) solves for
%   the half-angle of the conical frustum that smoothly connects the base
%   ring at (z_base, r_base) to the spherical tip seam at (z_seam, r_tip).

    if abs(z_seam - z_base) < 1e-9
        alpha = pi/2;
        return;
    end

    rs = r_tip;
    rb = r_base;
    slope_cone = (rs - rb) / (z_seam - z_base); % dr/dz along cone

    alpha = atan2(1, -slope_cone);
    if alpha <= 0
        alpha = alpha + pi;
    end

    % Clamp to avoid singular configurations exactly at 0 or pi
    alpha = min(max(alpha, 2*pi/180), pi - 2*pi/180);
end
