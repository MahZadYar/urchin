function [V, idx] = append_vertices(V, new_pts)
%APPEND_VERTICES Concatenate new vertices and report their indices.
%   [V, idx] = APPEND_VERTICES(V, new_pts) appends the rows in new_pts to
%   the existing vertex array V, returning the updated matrix and the
%   indices of the appended vertices.

    if isempty(new_pts)
        idx = zeros(0, 1);
        return;
    end
    startIdx = size(V, 1) + 1;
    V = [V; new_pts]; %#ok<AGROW>
    idx = (startIdx:size(V, 1)).';
end
