function b = toeplitz(c,r)
% function B = toeplitz(C,R)
%
% DESCRIPTION
%   Construct Toeplitz Matrices
%
% INPUTS
%   C: polynomial
%   R: integer
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=toeplitz(C,R)
%     B is a Toeplitz matrix having C as its first column and R
%     as its first row.
%   B=toeplitz(C)
%     This is the same as B=toeplitz(C,C), i.e. B is a Toeplitz
%     matrix with C as the first row and first column.
%
% See also hankel

% 12/22/2011: PJS  Initial Coding

if nargin==1
    r = c;
end

% Warning if r(1) and c(1) are not the same
L.type = '()';
L.subs = {1};
r1 = subsref(r,L);
c1 = subsref(c,L);
if ~isequal(r1,c1)
    toeplitz(1,2);
end

% Get indices from double/toeplitz
Nr = length(r);
Nc = length(c);
T = toeplitz(1:Nc,[1 Nc+(1:Nr-1)]);

% Create polynomial Toeplitz matrix
% This first converts r to a column and then grabs 2:end
L.subs = {':'};
ccol = subsref(c,L);
rcol = subsref(r,L);
L.subs = {2:Nr};
rcol = subsref(rcol,L);
v = [ccol; rcol];

L.subs = { T(:) };
vT = subsref(v,L);

L.subs = {':'};
b = polynomial(zeros(Nc,Nr));
b = subsasgn(b,L,vT);


