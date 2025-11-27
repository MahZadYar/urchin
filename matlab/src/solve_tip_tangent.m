function [rTipSeam, alphaCone, zTipSeam, isValid] = solve_tip_tangent(rBase, coreRadius, zTipCenter, tipSphereRadius)
%SOLVE_TIP_TANGENT Compute tangency between a cone and spherical tip.
%   [rTipSeam, alphaCone, zTipSeam, isValid] = solve_tip_tangent(rBase, coreRadius,
%   zTipCenter, tipSphereRadius) returns the seam radius, the half-angle
%   at the tangent point, the axial seam location, and a validity flag. The
%   flag indicates whether the closed-form tangency was feasible without
%   clamping (true) or whether fallback heuristics were applied (false).

	rBase = max(rBase, 0);
	if coreRadius <= 0
		coreRadius = 1e-9;
	end
	if rBase >= coreRadius
		rBase = min(rBase, 0.999 * coreRadius);
	end

	rTipSphere = tipSphereRadius;
	if rTipSphere <= 0
		zBaseLocal = sqrt(max(0, coreRadius^2 - rBase^2));
		rTipSeam = 0;
		alphaCone = pi/2;
		zTipSeam = zBaseLocal;
		isValid = false;
		return;
	end

	zBaseLocal = sqrt(max(0, coreRadius^2 - rBase^2));
	A = zTipCenter - zBaseLocal;
	if A <= 0
		rTipSeam = max(0, min(rBase, rTipSphere));
		ratio = min(1, max(0, rTipSeam / max(rTipSphere, 1e-12)));
		alphaCone = asin(ratio);
		if alphaCone < 1e-6
			alphaCone = 1e-6;
		end
		zTipSeam = max(zBaseLocal, zTipCenter);
		isValid = false;
		return;
	end

	D2 = A^2 + rBase^2;
	term = D2 - rTipSphere^2;
	if term < 0
		term = 0;
		isValid = false;
	else
		isValid = true;
	end

	sqrtTerm = sqrt(term);
	denom = max(D2, 1e-12);
	rTipSeam = (rBase * rTipSphere^2 + rTipSphere * A * sqrtTerm) / denom;
	rTipSeam = max(0, min(rTipSeam, rTipSphere));

	ratio = min(1, max(0, rTipSeam / max(rTipSphere, 1e-12)));
	alphaCone = asin(ratio);
	if alphaCone < 1e-6
		alphaCone = 1e-6;
	end

	C = (rBase * rTipSeam - rTipSphere^2) / max(A, 1e-12);
	zTipSeam = zTipCenter + C;
	if zTipSeam <= zBaseLocal
		zTipSeam = zBaseLocal + max(1e-6, min(rTipSphere, A));
		isValid = false;
	end
end
