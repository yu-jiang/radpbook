function b = hankel(c,r)
% function B = hankel(C,R)
%
% DESCRIPTION
%   Construct Hankel Matrices
%
% INPUTS
%   C: polynomial
%   R: integer
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B=hankel(C,R)
%     B is a Hankel matrix having C as its first column and R
%     as its last row.
%   B=hankel(C)
%     This is the same as B=hankel(C,C), i.e. B is a Hankel
%     matrix with C as the first column and last row.
%
% See also toeplitz

% 12/22/2011: PJS  Initial Coding

if nargin==1
    r = c;
end

% Warning if r(1) and c(1) are not the same
L.type = '()';
L.subs = {1};
r1 = subsref(r,L);

L.subs = {length(c)};
cend = subsref(c,L);
if ~isequal(r1,cend)
    hankel(1,2);
end

% Get indices from double/toeplitz
Nc = length(c);
Nr = length(r);
T = hankel(1:Nc,[Nc Nc+(1:Nr-1)]);

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


