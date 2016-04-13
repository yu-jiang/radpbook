function c = plus(a,b)
% function C=plus(A,B)
%
% DESCRIPTION
%   Add two polynomials and combine common terms.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C: polynomial, the result of addition A+B
%
% SYNTAX
%   C= A+B
%     Adds A and B.  A and B must have the same dimensions unless one is
%     is a scalar.  A scalar can be added to anything.
%   C = plus(A,B)
%     Function-call form for addition.

% 6/7/2002: PJS  Initial Coding
% 6/8/2002: PJS  Allow matrices of polynomials

% Promote a to polynomial
a = polynomial(a);
sza = size(a);

% Promote b to polynomial
b = polynomial(b);
szb = size(b);

if isempty(a) || isempty(b)
    
    if isempty(a) && all(szb==[1 1])
        % empty+scalar = empty(sza)
        c=polynomial(zeros(sza));
        return;
    elseif isempty(b) && all(sza==[1 1])
        % scalar+empty = empty(szb)
        c=polynomial(zeros(szb));
        return;
    elseif all(sza==szb)
        c=polynomial(zeros(sza));
        return;
    else
        error('Matrix dimensions must agree.');
    end
    
elseif all(sza==szb)
    % Matrix + Matrix
    
    % Get Dimensions
    nta = size(a.degmat,1);
    nva = length(a.varname);
    ntb = size(b.degmat,1);
    nvb = length(b.varname);
    
    % Stack up coefficients
    acoef = a.coefficient;
    bcoef = b.coefficient;
    coefficient = [acoef; bcoef];
    
    if nva==0 && nvb==0
        degmat = zeros(nta+ntb,0);
        varname = {};
    else
        adeg = a.degmat;
        bdeg = b.degmat;
        degmat = blkdiag(adeg,bdeg);
        varname = [a.varname(:); b.varname(:)];
    end
    
    % Form polynomial and combine terms
    chkval = 0; % skip validity check
    c = polynomial(coefficient,degmat,varname,sza,chkval);
    c = combine(c);
    
elseif all(sza==[1 1]) || all(szb==[1 1])
    % Scalar + Matrix / Matrix + Scalar
    
    % Make first term of sum be the scalar.
    if ~all(sza==[1 1])
        temp = a;
        a = b;
        b = temp;
        szb = sza;
    end
    
    % Turn into matrix+matrix
    a.coefficient = repmat(a.coefficient,1,szb(1)*szb(2));
    a.matdim = szb;
    c = plus(a,b);
else
    error('Matrix dimensions must agree');    
end
