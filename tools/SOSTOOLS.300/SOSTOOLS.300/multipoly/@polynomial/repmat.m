function B = repmat(A,M,N)
% function B = repmat(A,M,N)
%
% DESCRIPTION
%   Replicate and tile an array of polynomials
%
% INPUTS
%   A: polynomial matrix
%   M: number of tiles in the row dimension
%   N: number of tiles in the column dimension
%
% OUTPUTS
%   B: polynomial matrix
%
% SYNTAX
%   B = repmat(A,M,N)
%       Creates a large matrix B consisting of an M-by-N tiling of copies
%       of A. The size of B is [size(A,1)*M, size(A,2)*N].
%   B = repmat(A,N)
%       Creates an N-by-N tiling of A.
%   B = REPMAT(A,[M N])
%       Equivalent to repmat(A,M,N).
%
% EXAMPLE
%   pvar x1 x2;
%   p = [3*x2, x1^2+x2];
%   M = repmat(p,[3 2])

% 10/25/2010 PJS   Initial Coding

% Argument checking
error(nargchk(2,3,nargin))
if nargin==2
    if isscalar(M)
        N = M;
    elseif isvector(M) && length(M)==2
        N = M(2);
        M = M(1);
    else
        error('For two argument syntax the second input must be a scalar or 2-by-1 vector.');
    end
end
if ~( isscalar(M) && isa(M,'double') && ceil(M)==floor(M) )
    error('Row tiling dimension must be a scalar integer.');
end
if ~( isscalar(N) && isa(N,'double') && ceil(N)==floor(N) )
    error('Column tiling dimension must be a scalar integer.');
end

% Tile in column direction
% This is the easier direction due to format of coef matrix
B = A;
if N>1
    szB = size(B);
    B.coefficient = repmat( B.coefficient , [1, N] );
    B.matdim = [szB(1), szB(2)*N];
end

% Tile in row direction
% Transpose, perform tile in col dimension, and then transpose back
if M>1
    Bt = B';
    szBt = size(Bt);
    Bt.coefficient = repmat( Bt.coefficient , [1, M] );
    Bt.matdim = [szBt(1), szBt(2)*M];
    B = Bt';
end


