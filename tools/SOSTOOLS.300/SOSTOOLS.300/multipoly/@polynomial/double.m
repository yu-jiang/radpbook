function b=double(a)
%
% DESCRIPTION
%   Convert constant polynomial object to double.
%
% INPUTS
%   A: constant polynomial object
%
% OUTPUTS
%   B: double
%
% SYNTAX
%   B= double(A)
%

% 3/15/03: SP     Initial coding
% 3/25/03: PJS    Remove "get" within methods

if isempty(a.varname)
    sza = size(a);
    b = reshape(a.coefficient,sza(1),sza(2));
    b = full(b);
else
    error('Input contains independent variables.');
end

