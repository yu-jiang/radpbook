function b = uminus(a)
% function B=uminus(A)
%
% DESCRIPTION 
%   Negate the elements of a polynomial matrix.
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: unary minus of input, -A
%  
% SYNTAX 
%   B = -A;
%     Negates the elements of A. 
%   B = uminus(A)
%     Function-call form for negation.

% 6/8/2002: PJS  Initial Coding  

b = a;
b.coefficient = -b.coefficient;

