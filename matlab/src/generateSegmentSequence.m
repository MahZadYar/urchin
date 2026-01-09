function counts = generateSegmentSequence(startCount, endCount, steps)
%GENERATESEGMENTSEQUENCE Interpolate segment counts along the cone axis.
%   counts = GENERATESEGMENTSEQUENCE(startCount, endCount, steps) returns a
%   vector of length steps+1 with integer segment counts transitioning
%   between startCount and endCount using a greedy shortest-diagonal rule.

    if steps <= 0
        counts = [startCount; endCount];
        return;
    end

    counts = zeros(steps + 1, 1);
    counts(1) = startCount;
    counts(end) = endCount;
    diffVal = endCount - startCount;
    numChanges = abs(diffVal);
    if numChanges == 0
        for idx = 2:(steps + 1)
            counts(idx) = startCount;
        end
        return;
    end

    changeIdx = round(linspace(1, steps, numChanges));
    changeMask = false(steps, 1);
    changeMask(changeIdx) = true;
    stepSign = sign(diffVal);
    for s = 1:steps
        counts(s+1) = counts(s);
        if changeMask(s)
            counts(s+1) = counts(s+1) + stepSign;
        end
    end
end
