function b = det(a)
% function b=det(a)
%
% DESCRIPTION
%   Determinant of a matrix polynomial
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: polynomial, the determinant of A
%
% SYNTAX
%   B = det(A);


% 6/14/2002: PJS  Initial Coding

sza=size(a);
if sza(1)~=sza(2)
    error('Matrix must be square');
end

if isempty(a)
    b = polynomial(1);
elseif sza(1)==1
    b = a;
else
    L.type = '()';
    b = polynomial(0);
    for i1 = 1:sza(1);
        % Recursive cofactor expansion for det
        % XXX Faster algos for computing det exist
        L.subs = {[2:sza(1)] [1:i1-1, i1+1:sza(1)]};
        M = subsref(a,L);
        cofactor = (-1)^(1+i1)*det(M);
        
        L.subs = {1, i1};
        a1_i1 = subsref(a,L);
        b = b + a1_i1*cofactor;
    end
end




