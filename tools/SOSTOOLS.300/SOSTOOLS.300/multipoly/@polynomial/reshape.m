function b = reshape(a,m,n)
% function B = reshape(A,M,N)
%
% DESCRIPTION
%   Returns an M-by-N polynomial matrix whose elements are taken
%   columnwise from A.  The input, A, must have M*N elements.
%
% INPUTS
%   A: polynomial
%   M: Row dimension
%   N: Column dimension
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=reshape(A,M,N)
%   B=reshape(A,[M N])

% 3/30/2003: PJS   Initial Coding


[nra,nca]=size(a);
if nargin==2
    n = m(2);
    m = m(1);
end
if (nra*nca~=m*n)
    error('To RESHAPE the number of elements must not change.');
else
    b = a;
    b.matdim = [m n];
end


