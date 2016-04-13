function b = ctranspose(a)
% function B=ctranspose(A)
%
% DESCRIPTION
%   Returns the complex conjugate transpose of a polynomial.
%   This currently assumes that the polynomial coefficients and
%   variables are real. Thus ctranspose actually returns the non-conjugate
%   transpose.  ctranspose is the same as transpose.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: transpose of A.
%
% SYNTAX
%   B = ctranspose(A);
%   B = A';
%
% See also transpose

% 6/14/2002: PJS  Initial Coding

b=transpose(a);


