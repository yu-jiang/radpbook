function b = transpose(a)
% function B=transpose(A)
%
% DESCRIPTION
%   Returns the non-conjugate transpose of a polynomial.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: transpose of A.
%
% SYNTAX
%   B= transpose(A);
%   B= A.';
%
% See also ctranspose

% 6/7/2002: PJS  Initial Coding

b=a;

% Get Dimensions
[nrb,ncb]=size(b);
bcoef = b.coefficient;

idx = [];
tempidx = 0:nrb:(ncb-1)*nrb;
for i1 = 1:nrb
    idx = [idx i1+tempidx];
end
bcoef = bcoef(:,idx);
b.coefficient = bcoef;
b.matdim = [ncb nrb];

