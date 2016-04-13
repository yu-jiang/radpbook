function b = isempty(a)
% function B = isempty(A)
%
% DESCRIPTION
%   Returns true for empty polynomials.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: true if A is an empty polynomial and false otherwise.
%
% SYNTAX
%   B = isempty(A);

% 6/14/2002: PJS  Initial Coding

if min(a.matdim)==0
    b = true;
else
    b = false;
end

