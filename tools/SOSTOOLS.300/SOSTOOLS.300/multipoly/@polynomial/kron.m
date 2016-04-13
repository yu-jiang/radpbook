function b = kron(x,y)
% function B = kron(X,Y)
%
% DESCRIPTION
%   Kronecker tensor product
%
% INPUTS
%   X: rx-by-cx polynomial matrix
%   Y: ry-by-cy polynomial matrix
%
% OUTPUTS
%   B: rx*ry-by-cx*cy polynomial matrix
%
% SYNTAX
%   B=kron(X,Y)
%     B is the Kronecker product of X and Y.
%

% 9/9/2013: PJS  Initial Coding

% Check # of input arguments
if nargin~=2
    error('Syntax is B=kron(X,Y)');
    %narginchk(2,2);
end

% Get matrix dimensions
[rx,cx] = size(x);
[ry,cy] = size(y);

% Inialize variable for parentheses subsref
L.type = '()';
L.subs = {};

% Perform KRON operation 
% Note: This is probably not the most efficient implementation for polys
b = polynomial( zeros(rx*ry,cx*cy) );
for i = 1:rx
    for j = 1:cx
        L.subs = {i,j};
        xij = subsref(x,L);
        
        ridx = (1:ry)+(i-1)*ry;
        cidx = (1:cy)+(j-1)*cy;
        L.subs = {ridx,cidx};
        b = subsasgn(b,L,xij*y);
    end
end
