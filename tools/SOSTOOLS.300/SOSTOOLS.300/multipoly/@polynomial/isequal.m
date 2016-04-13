function c = isequal(a,b)
% function C=isequal(A,B)
%
% DESCRIPTION
%   Element-by-element comparisons of A and B.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C: Logical array of same size as A and B with the i^th element of C
%      set to true if the i^th elements of A and B are the same. Otherwise
%      the i^th element of C is set to false.
%
% SYNTAX
%   C= isequal(A,B)
%     Element-by-element comparison of A and B.  A and B must have the
%     same dimensions unless one is is a scalar. Scalars can be compared
%     to anything.


% Note: This syntax for isequal is different from the standard syntax.
% Typically EQ does element-by-element comparisons and ISEQUAL does an
% overall comparison of two arrays.  For polynomials EQ is overloaded to
% create poly constraint objects and ISEQUAL does element-by-element
% comparisons.

% 6/14/2002: PJS  Initial Coding
% 11/17/2010: PJS Modified == code

a = polynomial(a);
sza = size(a);

b = polynomial(b);
szb = size(b);

if isempty(a) || isempty(b)
    
    if isempty(a) && all(szb==[1 1])
        % empty==scalar returns empty(sza)
        c= false(sza);
        return;
    elseif isempty(b) && all(sza==[1 1])
        % scalar==empty returns empty(szb)
        c= false(szb);
        return;
    elseif all(sza==szb)
        % empty==empty returns empty
        c= false(sza);
        return;
    else
        error('Matrix dimensions must agree.');
    end
elseif all(sza==szb)
    % Matrix == Matrix case
    d = a-b;
    [nrd,ncd] =size(d);
    dcoef = d.coefficient;
    dcoef = sum(abs(dcoef),1);
    c = full(dcoef)==0;
    c = reshape(c,nrd,ncd);
    
elseif all(sza==[1 1]) || all(szb==[1 1])
    % Scalar == Matrix:
    
    % Make first term be the scalar.
    if ~all(sza==[1 1])
        temp = a;
        a = b;
        b = temp;
        szb = sza;
    end
    
    % Turn into matrix==matrix
    a.coefficient = repmat(a.coefficient,[1 szb(1)*szb(2)]);
    a.matdim = szb;
    c = isequal(a,b);
else
    error('Matrix dimensions must agree');
end



