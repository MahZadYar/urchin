function volume = computeMeshVolume(vertices, faces)
%COMPUTEMESHVOLUME Calculate enclosed volume of a triangular surface mesh.
%
%   volume = COMPUTEMESHVOLUME(vertices, faces) returns the absolute volume
%   enclosed by a triangular mesh defined by vertices and face indices. The
%   mesh is assumed to be consistently oriented.

    if isempty(vertices) || isempty(faces)
        volume = 0;
        return;
    end

    v1 = vertices(faces(:,1), :);
    v2 = vertices(faces(:,2), :);
    v3 = vertices(faces(:,3), :);
    tetraVol = dot(v1, cross(v2, v3, 2), 2) / 6;
    volume = abs(sum(tetraVol));
end
