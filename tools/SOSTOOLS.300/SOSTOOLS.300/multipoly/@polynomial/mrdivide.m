function c = mrdivide(a,b)
% function C=mrdivide(A,B)
%
% DESCRIPTION
%   Right matrix divide
%
% INPUTS
%   A: polynomial
%   B: matrix of constants
%
% OUTPUTS
%   C: polynomial, the result of division.
%
% SYNTAX
%   C= A/B
%   C = mrdivide(A,B)
%     Function-call form

% 10/22/2002: PJS  Initial Coding

% Promote a to polynomial
a = polynomial(a);
sza = size(a);

% Promote b to polynomial
b = polynomial(b);
szb = size(b);


if isempty(a) || isempty(b)
    
    if isempty(a) && all(szb==[1 1])
        % empty/scalar = empty(sza)
        c=polynomial(zeros(sza));
        return;
    elseif sza(2)==szb(2)
        c=polynomial([]);
        return;
    else
        error('Matrix dimensions must agree.');
    end
    
elseif all(szb==[1 1])
    
    % Matrix/Scalar:
    degidx = find(b.degmat);
    if ~isempty(degidx)
        error('B must be a constant.');
    end
    b = combine(b);
    bcoef = b.coefficient;
    c = mtimes(a,1/bcoef);
    
elseif sza(2)==szb(2)
    
    % Matrix/Matrix
    degidx = find(b.degmat);
    if ~isempty(degidx)
        error('B must be a matrix of constants.');
    end
    b = combine(b);
    bcoef = b.coefficient;
    c = mtimes(a,inv(bcoef));
    
else
    error('Matrix dimensions must agree');
end









