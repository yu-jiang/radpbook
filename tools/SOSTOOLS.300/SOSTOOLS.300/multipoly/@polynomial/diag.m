function b = diag(a,k)
% function B = diag(A,K)
%
% DESCRIPTION
%   Diagonal polynomial matrices and diagonals of polynomial matrices.
%
% INPUTS
%   A: polynomial
%   K: integer
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=diag(A,K)
%     If A is a polynomial vector of length N then B is an M-by-M
%     square matrix with M=N+ABS(K) and the elements of A on the K^th
%     diagonal.  K=0 is the main diagonal, K>0 is above the main diagonal
%     and K<0 is below the main diagonal
%   B=diag(A,K)
%     If A is a polynomial matrix then B is a column vector whose elements
%     are the entries of the K^th diagonal of A.
%   B=diag(A)
%     This is the same as B=diag(A,0).
%
% See also triu, tril, blkdiag

% 12/3/2009: PJS  Initial Coding

if nargin==1
    k=0;
end

sza = size(a);
if any(sza==1)
    % a is a vector --> b is a matrix
    
    % Use double/diag to find right dimensions and indices
    la = length(a);
    bmat = diag(1:la,k);
    idx = find(bmat);
    
    % Assign values
    b = polynomial( zeros(size(bmat)) );
    L.type = '()';
    L.subs = {idx};
    b = subsasgn(b,L,a);
else
    % a is a matrix --> b is a vector
    if k>=sza(2) || k <= -sza(1)
        b = polynomial(zeros(0,1));
        return
    end
    
    % Use double/diag to find right indices
    amat = reshape( 1:(sza(1)*sza(2)), sza );
    idx = diag(amat,k);
    
    % Assign values
    L.type = '()';
    L.subs = {idx};
    b = subsref(a,L);
end

