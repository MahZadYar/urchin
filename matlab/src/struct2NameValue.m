function args = struct2NameValue(s)
%STRUCT2NAMEVALUE Convert struct to name-value pair cell array.
%
%   args = STRUCT2NAMEVALUE(s) returns a 1x(2N) cell array suitable for
%   calling functions that accept name-value pairs. Fields with empty values
%   are skipped.

    f = fieldnames(s);
    tmp = cell(1, numel(f) * 2);
    idx = 0;
    for i = 1:numel(f)
        value = s.(f{i});
        if isempty(value)
            continue;
        end
        idx = idx + 1;
        tmp{2*idx-1} = f{i};
        tmp{2*idx} = value;
    end
    args = tmp(1:2*idx);
end
