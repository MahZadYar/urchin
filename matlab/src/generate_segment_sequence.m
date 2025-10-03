function counts = generate_segment_sequence(startCount, endCount, steps)
%GENERATE_SEGMENT_SEQUENCE Interpolate segment counts along the cone axis.
%   counts = GENERATE_SEGMENT_SEQUENCE(startCount, endCount, steps) returns a
%   vector of length steps+1 with integer segment counts transitioning
%   between startCount and endCount. The interpolation is spacing-driven so it
%   can accommodate large differences even with a small number of axial steps.

    if steps <= 0
        counts = [startCount; endCount];
        return;
    end

    counts = round(linspace(startCount, endCount, steps + 1));
    counts(1) = startCount;
    counts(end) = endCount;
    counts = max(counts, 1);
end
