function out = end(a,k,n)
% function Out=end(A,K,N)
%
% DESCRIPTION
%   Called for indexing expressions involving the polynomial A
%   when END is part of the K-th index out of N indices.
%
% INPUTS
%   A: polynomial
%   K,N: indices (natural numbers)
%
% OUTPUTS
%   Out: The size of A along dimension K.
%
% SYNTAX
%   Out=end(A,K,N)
%     Returns the size of A along dimension K.

% 7/2/2002: PJS  Initial Coding

if n ==1
    sza = size(a);
    out = sza(1)*sza(2);
elseif n==2
    out = size(a,k);
end
