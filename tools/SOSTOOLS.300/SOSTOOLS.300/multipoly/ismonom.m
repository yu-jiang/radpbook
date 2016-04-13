function b = ismonom(a)
% function B = ismonom(A)
%
% DESCRIPTION
%   Returns true for a monomial or a vector or matrix of monomials.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: 1 if A is a monomial or a vector or matrix of monomials.  Returns
%      0 otherwise.
%
% SYNTAX
%   B = ismonom(A);

% 1/29/2008: PJS  Initial Coding

if isa(a,'double') && a==1
    b = true;
    return;
elseif ~isa(a,'polynomial')
    b = false;
    return;
end

% acoef will be Nterms-by-(Nrows*Ncols)
acoef = a.coefficient;

if all(nonzeros(acoef)==1) && all(sum(acoef,1)==1)
    b = true;
else
    b = false;
end
