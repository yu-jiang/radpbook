function b = diff(a,x)
% function B=diff(A,X)
%
% DESCRIPTION
%   Element-by-element differentiation of a polynomial with respect
%   to a single variable.
%
% INPUTS
%   A: polynomial
%   X: Differentiate with respect to the (single) variable X.
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B = diff(A,X);
%     Differentiate the polynomial, A, with to X.  A should be a
%     polynomial and X should be a polynomial variable or a string.
%     Differentiation is done element-by-element if A is a matrix.
%
% EXAMPLE
%   pvar x y z;
%   f = 2*x^3+5*y*z-2*x*z^2;
%   df = diff(f,x)
%
% See also: jacobian, int

% 11/5/2002 PJS  Initial Coding
% 10/21/2010 PJS diff calls jacobian (reduces redundant code)

% Error Checking
if nargin~=2
    error('Error in calling diff');
end

if ispvar(x) && length(x)==1
    x = char(x);
    x = x{1};
elseif ~ischar(x)
    error('X must be a single polynomial variable or a string');
end

% Call jacobian to perform differentiation
% jacobian works on a column vector and returns a column vector
sza = size(a);
b = jacobian(a(:),x);

% Reshape derivative back to size of a
b = reshape(b,sza);

