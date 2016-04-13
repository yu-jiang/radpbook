function out = get(a,property)
% function Out = get(A,Property)
%
% DESCRIPTION
%   Get polynomial properties.
%
% INPUTS
%   A: polynomial
%   Property: polynomial property
%
% OUTPUTS
%   Out: Value of the property
%
% SYNTAX
%   Out = get(A,Property);

% 6/7/2002: Written by Packard
% 10/30/2002 PJS  Make PVget act on pmat

GSp = PVget(a,'PropNames');
if nargin==2
    % Use CaseInsensitivePartialMatch to find the correct property name
    [idx,fstr] = PVcipmatch(GSp.GPropNames,property);
    if isempty(idx)
        error(['Property ' property ' is not gettable']);
    else
        % Extract value with PVGET
        out = PVget(a,fstr);
    end
else
    % Create an empty STRUCT, with all the entries from GProp as
    % fieldsnames
    nprop = length(GSp.GPropNames);
    plist = [GSp.GPropNames cell(nprop,1)]';
    tout = struct(plist{:});
    % Fill each field with the appropriate value
    for i=1:nprop
        tout = setfield(tout,GSp.GPropNames{i}, ...
            PVget(a,GSp.GPropNames{i}));
    end
    if nargout==1
        out = tout;
    else
        % If there are no output arguments, simply display the result
        disp(tout);
    end
end


