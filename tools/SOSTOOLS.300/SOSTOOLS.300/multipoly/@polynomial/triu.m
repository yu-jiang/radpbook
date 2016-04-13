function b = triu(a,k)
% function B = triu(A,K)
%
% DESCRIPTION
%   Extract upper triangular part of a polynomial matrix
%
% INPUTS
%   A: polynomial
%   K: integer
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=triu(A,K)
%     B contains the elements on and above the K-th diagonal of A. K=0
%     is the main diagonal, K>0 is above the main diagonal, and K<0 is
%     below the main diagonal.
%   B=triu(A)
%     This is the same as B=triu(A,0).
%
% See also tril, diag

% 11/24/2009: PJS  Initial Coding

if nargin==1
    k=0;
end

sza = size(a);
% a is a matrix --> b is a vector
if k>=sza(2)
    b = polynomial(zeros(sza));
    return
elseif k <= -sza(1)
    b = a;
    return;
end

% Use double/triu to find indices to zero out
idx = find( ~triu(ones(sza),k) );

% Assign values
b = a;
L.type = '()';
L.subs = {idx};
b = subsasgn(b,L,0);


