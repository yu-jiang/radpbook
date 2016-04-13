function out = length(a)
% function Out=length(A)
%
% DESCRIPTION
%   Length of a polynomial vector.  For a matrix, it is equivalent
%   to max(size(A)) for non-empty arrays and 0 for empty ones.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   Out: Length of A
%
% SYNTAX
%   Out=length(A);


% 6/10/2002: PJS  Initial Coding

out = max(size(a));




