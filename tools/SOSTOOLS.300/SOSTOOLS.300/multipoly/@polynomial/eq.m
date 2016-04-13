function pco = eq(p,q)
% function pco = eq(p,q)
%
% DESCRIPTION
%   Creates a polynomial constraint object representing the 
%   constraint p==q.  sosopt interprets this constraint as
%           p-q = zero polynomial
%
% INPUTS
%   p: polynomial
%   q: polynomial
%
% OUTPUT
%   pco: polynomial constraint object
%
% SYNTAX
%   p==q 
%     Creates a polynomial constraint object.  If p is a scalar and
%     q is a vector then p will be expanded to have the same dimension as
%     q.  Similarly, if q is a scalar and p is a vector then q will be 
%     expanded.
%

% 10/22/2010:   PJS  Initial Coding

% Create polynomial constraint object
pco = polyconstr(p,q,'==');

