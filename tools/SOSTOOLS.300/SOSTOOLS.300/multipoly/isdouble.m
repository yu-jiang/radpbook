function b = isdouble(a)
% function B = isdouble(A)
%
% DESCRIPTION
%   Returns true if A is a double.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: treu if A is a double or a vector or matrix of doubles.  Returns
%      false otherwise.
%
% SYNTAX
%   B = isdouble(A);

% 12/9/2010: PJS  Initial Coding

if isa(a,'double')
    b = true;
elseif isa(a,'polynomial')
    b = (a.nvar==0);
else
    b = false;
end
