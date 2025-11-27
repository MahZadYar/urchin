function [rBaseNominal, rTipNominal, alphaConeNominal, zTipSeamNominal] = compute_nominal_spike_profile(coreRadius, spikeLength, tipSphereRadius, conicalityEffective, rBaseLimitLocal, scale)
%COMPUTE_NOMINAL_SPIKE_PROFILE Baseline spike profile honoring tip tangency and conicality.
%
%   [rBaseNominal, rTipNominal, alphaConeNominal, zTipSeamNominal] =
%   COMPUTE_NOMINAL_SPIKE_PROFILE(coreRadius, spikeLength, tipSphereRadius,
%   conicalityEffective, rBaseLimitLocal, scale) returns the baseline cone
%   base radius, seam radius, seam angle, and seam height for a spike. The
%   profile honours the requested conicality while staying within the
%   tangency limits imposed by the spherical tip and core geometry.

    rBaseLimitLocal = max(0, min(rBaseLimitLocal, 0.999 * coreRadius));
    spikeLength = max(spikeLength, 0);
    zApex = coreRadius + spikeLength;
    tipSphereRadius = max(0, tipSphereRadius);
    zTipCenter = zApex - tipSphereRadius;

    if tipSphereRadius > 0 && rBaseLimitLocal > 0
        rBase = min(max(min(tipSphereRadius, rBaseLimitLocal), 0), rBaseLimitLocal);
        for iter = 1:8
            [rTipCandidate, ~, ~, tangentValid] = solve_tip_tangent(rBase, coreRadius, zTipCenter, tipSphereRadius);
            if ~tangentValid
                rTipCandidate = min(rBase, tipSphereRadius);
            end
            if conicalityEffective >= 0
                rBaseTarget = rTipCandidate + conicalityEffective * (rBaseLimitLocal - rTipCandidate);
            else
                rBaseTarget = max(0, (1 + conicalityEffective) * rTipCandidate);
            end
            rBaseTarget = min(rBaseTarget, rBaseLimitLocal);
            if abs(rBaseTarget - rBase) <= 1e-9 * max(scale, 1)
                rBase = rBaseTarget;
                break;
            end
            rBase = rBaseTarget;
        end

        [rTipNominal, alphaConeNominal, zTipSeamNominal, tangentValid] = solve_tip_tangent(rBase, coreRadius, zTipCenter, tipSphereRadius);
        if ~tangentValid
            rTipNominal = min(rBase, tipSphereRadius);
            ratio = min(1, max(rTipNominal / max(tipSphereRadius, 1e-12), 0));
            alphaConeNominal = max(1e-6, asin(ratio));
            zTipSeamNominal = zTipCenter + tipSphereRadius * cos(alphaConeNominal);
        end
        rBaseNominal = rBase;
    else
        rBaseNominal = min(rBaseLimitLocal, max(0, min(tipSphereRadius, rBaseLimitLocal)));
        rTipNominal = 0;
        alphaConeNominal = pi/2;
        zTipSeamNominal = zApex;
    end

    rBaseNominal = max(rBaseNominal, 0);
    rTipNominal = max(min(rTipNominal, tipSphereRadius), 0);
    alphaConeNominal = min(max(alphaConeNominal, 1e-6), pi/2);
    zTipSeamNominal = max(zTipSeamNominal, coreRadius);
end
