function tok = valueToken(value)
%VALUETOKEN Convert numeric/logical/string values to safe filename tokens.
%
%   tok = VALUETOKEN(value) returns a lower-case string containing only
%   [a-z0-9] characters, suitable for use in filenames.

    if isnumeric(value)
        tok = regexprep(sprintf('%.3g', double(value)), '[^0-9A-Za-z]', '');
        if isempty(tok)
            tok = "0";
        else
            tok = string(tok);
        end
    elseif islogical(value)
        tok = string(double(value));
    else
        tok = regexprep(lower(string(value)), '[^a-z0-9]', '');
    end
end
