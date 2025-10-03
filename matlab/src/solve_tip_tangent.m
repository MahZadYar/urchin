function [r_seam, alpha, z_seam, isValid] = solve_tip_tangent(r_base, coreRadius, z_tip_center, tipSphereRadius)
%SOLVE_TIP_TANGENT Compute tangency between a cone and spherical tip.
%   [r_seam, alpha, z_seam, isValid] = solve_tip_tangent(r_base, coreRadius,
%   z_tip_center, tipSphereRadius) returns the seam radius, the half-angle
%   at the tangent point, the axial seam location, and a validity flag. The
%   flag indicates whether the closed-form tangency was feasible without
%   clamping (true) or whether fallback heuristics were applied (false).

	r_base = max(r_base, 0);
	if coreRadius <= 0
		coreRadius = 1e-9;
	end
	if r_base >= coreRadius
		r_base = min(r_base, 0.999 * coreRadius);
	end

	R = tipSphereRadius;
	if R <= 0
		z_base_local = sqrt(max(0, coreRadius^2 - r_base^2));
		r_seam = 0;
		alpha = pi/2;
		z_seam = z_base_local;
		isValid = false;
		return;
	end

	z_base_local = sqrt(max(0, coreRadius^2 - r_base^2));
	A = z_tip_center - z_base_local;
	if A <= 0
		r_seam = max(0, min(r_base, R));
		ratio = min(1, max(0, r_seam / max(R, 1e-12)));
		alpha = asin(ratio);
		if alpha < 1e-6
			alpha = 1e-6;
		end
		z_seam = max(z_base_local, z_tip_center);
		isValid = false;
		return;
	end

	D2 = A^2 + r_base^2;
	term = D2 - R^2;
	if term < 0
		term = 0;
		isValid = false;
	else
		isValid = true;
	end

	sqrtTerm = sqrt(term);
	denom = max(D2, 1e-12);
	r_seam = (r_base * R^2 + R * A * sqrtTerm) / denom;
	r_seam = max(0, min(r_seam, R));

	ratio = min(1, max(0, r_seam / max(R, 1e-12)));
	alpha = asin(ratio);
	if alpha < 1e-6
		alpha = 1e-6;
	end

	C = (r_base * r_seam - R^2) / max(A, 1e-12);
	z_seam = z_tip_center + C;
	if z_seam <= z_base_local
		z_seam = z_base_local + max(1e-6, min(R, A));
		isValid = false;
	end
end
