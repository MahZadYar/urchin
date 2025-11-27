function [zTipSeams, rTipSeams, alphaCones] = solveConeSeams(rBases, rTips, zBases, zTipCenters)
%SOLVECONESEAMS Calculate the tangent seam between a cone body and a spherical tip.
%   [zTipSeams, rTipSeams, alphaCones] = solveConeSeams(rBases, rTips, zBases, zTipCenters)
%   computes the geometric properties of the circular seam where a conical spike
%   body meets its spherical tip cap tangentially.
%
%   Inputs:
%       rBases       - Radius of the cone base.
%       rTips        - Radius of the spherical tip.
%       zBases       - Axial position (Z) of the cone base.
%       zTipCenters  - Axial position (Z) of the center of the tip sphere.
%
%   Outputs:
%       zTipSeams    - Axial position (Z) of the tangent seam ring.
%       rTipSeams    - Radius of the tangent seam ring.
%       alphaCones   - Half-angle of the cone (radians).
%
%   The solver enforces tangency by ensuring the cone surface is tangent to the
%   tip sphere at the seam. It solves the system of equations derived from
%   similar triangles formed by the cone slope and the sphere normal.

    % We have two triangles that must be similar:
    % [1] tan(alphaCones) = (rBases - rTipSeams) / (zTipSeams - zBases)
    % [2] tan(alphaCones) = (zTipSeams - zTipCenters) / rTipSeams
    % [3] rTips^2 = (zTipSeams - zTipCenters)^2 + rTipSeams^2

    % Rearranging gives:
    % [1]+[2]=> [4] (rBases - rTipSeams) / (zTipSeams - zBases) = (zTipSeams - zTipCenters) / rTipSeams
  
    % [3] + [4] => quadratic in zTipSeams
    zTipSeams = (rTips.^2.*zBases + rBases.^2*zTipCenters - rTips.^2.*zTipCenters ...
               - 2*zBases.*zTipCenters.^2 + zBases.^2*zTipCenters + zTipCenters.^3 ...
               + rBases.*rTips.*(rBases.^2 - rTips.^2 + zBases.^2 - 2*zBases.*zTipCenters + zTipCenters.^2).^0.5) ...
               /(rBases.^2 + zBases.^2 - 2*zBases.*zTipCenters + zTipCenters.^2);
    
    % [3] =>
    rTipSeams = sqrt(rTips.^2 - (zTipSeams - zTipCenters).^2);
    alphaCones = atan((rBases - rTipSeams) ./ (zTipSeams - zBases));
    zTipSeams = real(zTipSeams);
    rTipSeams = real(rTipSeams);
    alphaCones = real(alphaCones);
end
