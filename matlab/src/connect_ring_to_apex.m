function F = connect_ring_to_apex(idxRing, idxApex)
%CONNECT_RING_TO_APEX Triangulate a fan between a ring and apex vertex.
%   Faces are returned with consistent winding following the ordering of
%   idxRing.

    n = numel(idxRing);
    F = zeros(n, 3);
    for k = 1:n
        k2 = mod(k, n) + 1;
        F(k, :) = [idxRing(k), idxRing(k2), idxApex];
    end
end
