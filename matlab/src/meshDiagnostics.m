function diagnostics = meshDiagnostics(surfaceOrStruct, warnOnIssues)
% MESHDIAGNOSTICS Compute standard quality checks for a surface mesh.
%   diagnostics = meshDiagnostics(surfaceMesh) evaluates watertightness,
%   manifoldness, orientability, and self-intersection flags for the provided
%   surface mesh. The input may be either a surfaceMesh object or any struct
%   containing a SurfaceMesh field (as returned by URCHIN).
%
%   diagnostics = meshDiagnostics(surfaceMesh, warnOnIssues) allows toggling
%   warning emission when issues are detected. Set warnOnIssues=false to
%   suppress warnings and retrieve the struct silently.
%
%   The diagnostics struct includes the following logical fields:
%       IsWatertight
%       IsEdgeManifold
%       IsOrientable
%       IsSelfIntersecting
%       IsVertexManifold
%
%   See also: URCHIN, SURFACEMESH

    if nargin < 1
        error('meshDiagnostics:MissingInput', 'A surface mesh or URCHIN struct is required.');
    end
    if nargin < 2 || isempty(warnOnIssues)
        warnOnIssues = true;
    else
        warnOnIssues = logical(warnOnIssues);
    end

    surfaceMeshObj = parseSurface(surfaceOrStruct);

    diagnostics = struct( ...
        'IsWatertight',       isWatertight(surfaceMeshObj), ...
        'IsEdgeManifold',     isEdgeManifold(surfaceMeshObj, false), ...
        'IsOrientable',       isOrientable(surfaceMeshObj), ...
        'IsSelfIntersecting', isSelfIntersecting(surfaceMeshObj), ...
        'IsVertexManifold',   isVertexManifold(surfaceMeshObj) ...
    );

    if warnOnIssues
        emitWarnings(diagnostics);
    end
end

function meshSurface = parseSurface(surfaceOrStruct)
    if isa(surfaceOrStruct, 'surfaceMesh')
        meshSurface = surfaceOrStruct;
        return;
    end
    if isstruct(surfaceOrStruct) && isfield(surfaceOrStruct, 'SurfaceMesh')
        meshSurface = surfaceOrStruct.SurfaceMesh;
        if isa(meshSurface, 'surfaceMesh')
            return;
        end
    end
    error('meshDiagnostics:InvalidInput', ...
        'Input must be a surfaceMesh or struct containing a SurfaceMesh field.');
end

function emitWarnings(diagnostics)
    if ~diagnostics.IsWatertight
        warning('Urchin mesh is not watertight. Consider adjusting parameters.');
    end
    if ~diagnostics.IsEdgeManifold
        warning('Urchin mesh has non-manifold edges.');
    end
    if ~diagnostics.IsVertexManifold
        warning('Urchin mesh has non-manifold vertices.');
    end
    if diagnostics.IsSelfIntersecting
        warning('Urchin mesh contains self-intersections.');
    end
end
