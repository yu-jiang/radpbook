function b = tril(a,k)
% function B = tril(A,K)
%
% DESCRIPTION
%   Extract lower triangular part of a polynomial matrix
%
% INPUTS
%   A: polynomial
%   K: integer
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=tril(A,K)
%     B contains the elements on and below the K-th diagonal of A. K=0
%     is the main diagonal, K>0 is above the main diagonal, and K<0 is
%     below the main diagonal.
%   B=tril(A)
%     This is the same as B=tril(A,0).
%
% See also triu, diag

% 11/24/2009: PJS  Initial Coding

if nargin==1
    k=0;
end

sza = size(a);
% a is a matrix --> b is a vector
if k>=sza(2)
    b = a;
    return
elseif k <= -sza(1)
    b = polynomial(zeros(sza));
    return;
end

% Use double/tril to find indices to zero out
idx = find( ~tril(ones(sza),k) );

% Assign values
b = a;
L.type = '()';
L.subs = {idx};
b = subsasgn(b,L,0);


