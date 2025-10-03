function [V, F] = triangulate_icosphere(R, subdivisions)
%TRIANGULATE_ICOSPHERE Generate a geodesic icosphere mesh.
%   [V, F] = TRIANGULATE_ICOSPHERE(R, subdivisions) creates a geodesic
%   sphere with radius R and subdivision depth. The output vertices lie on
%   the sphere of radius R, and faces contains triangle indices.

    t = (1 + sqrt(5)) / 2;
    verts = [
        -1,  t,  0;
         1,  t,  0;
        -1, -t,  0;
         1, -t,  0;
         0, -1,  t;
         0,  1,  t;
         0, -1, -t;
         0,  1, -t;
         t,  0, -1;
         t,  0,  1;
        -t,  0, -1;
        -t,  0,  1
    ];
    faces = [
        1,12,6; 1,6,2; 1,2,8; 1,8,11; 1,11,12;
        2,6,10; 6,12,5; 12,11,3; 11,8,7; 8,2,9;
        4,10,5; 4,5,3; 4,3,7; 4,7,9; 4,9,10;
        5,10,6; 3,5,12; 7,3,11; 9,7,8; 10,9,2
    ];
    verts = verts ./ vecnorm(verts, 2, 2);

    for s = 1:subdivisions
        midCache = containers.Map('KeyType','char', 'ValueType','int32');
        newFaces = zeros(size(faces, 1) * 4, 3);
        newVerts = verts;
        for i = 1:size(faces, 1)
            a = faces(i, 1); b = faces(i, 2); c = faces(i, 3);
            [ab, newVerts, midCache] = midpointIndex(a, b, newVerts, midCache);
            [bc, newVerts, midCache] = midpointIndex(b, c, newVerts, midCache);
            [ca, newVerts, midCache] = midpointIndex(c, a, newVerts, midCache);
            idx4 = (i - 1) * 4;
            newFaces(idx4 + 1, :) = [a, ab, ca];
            newFaces(idx4 + 2, :) = [b, bc, ab];
            newFaces(idx4 + 3, :) = [c, ca, bc];
            newFaces(idx4 + 4, :) = [ab, bc, ca];
        end
        faces = newFaces;
        verts = newVerts;
    end

    V = R * (verts ./ vecnorm(verts, 2, 2));
    F = faces;
end

function [idx, newVerts, midCache] = midpointIndex(a, b, newVerts, midCache)
%MIDPOINTINDEX Return the vertex index of the edge midpoint, caching results.
    if a > b
        key = sprintf('%d_%d', b, a);
    else
        key = sprintf('%d_%d', a, b);
    end
    if isKey(midCache, key)
        idx = midCache(key);
        return;
    end
    p = (newVerts(a, :) + newVerts(b, :)) / 2;
    p = p / norm(p);
    newVerts = [newVerts; p]; %#ok<AGROW>
    idx = size(newVerts, 1);
    midCache(key) = idx;
end
