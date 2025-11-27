function out = mergeStructs(varargin)
%MERGESTRUCTS Merge multiple scalar structs into one.
%
%   out = MERGESTRUCTS(s1, s2, ...) combines the fields of each scalar
%   struct argument into a single struct. Later structs overwrite earlier
%   values when field names collide. Empty inputs are ignored.

    out = struct();
    for k = 1:nargin
        s = varargin{k};
        if isempty(s)
            continue;
        end
        fields = fieldnames(s);
        for i = 1:numel(fields)
            out.(fields{i}) = s.(fields{i});
        end
    end
end
