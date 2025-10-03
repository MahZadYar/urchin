function insideFn = make_inside_tester(mesh, volAlpha)
%MAKE_INSIDE_TESTER Build a robust point-in-mesh tester closure.
%   insideFn = MAKE_INSIDE_TESTER(mesh, volAlpha) returns a function handle
%   that evaluates whether query points lie inside the mesh. The helper
%   prefers INPOLYHEDRON when available and falls back to an ALPHASHAPE
%   implementation when needed. volAlpha controls the alpha radius used by
%   the fallback. If empty, a heuristic based on the mesh bounding box is
%   chosen.

    F = mesh.Faces;
    V = mesh.Vertices;

    try
        % Probe inpolyhedron availability without processing large point sets
        inpolyhedron(F, V, mean(V, 1));
        insideFn = @(P) inpolyhedron(F, V, P);
        return;
    catch
        % Fallback: alphaShape
        if isempty(volAlpha)
            bbox = max(V, [], 1) - min(V, [], 1);
            volAlpha = 0.08 * max(bbox); % heuristic
        end
        shp = alphaShape(V(:,1), V(:,2), V(:,3), volAlpha);
        insideFn = @(P) inShape(shp, P(:,1), P(:,2), P(:,3));
    end
end
