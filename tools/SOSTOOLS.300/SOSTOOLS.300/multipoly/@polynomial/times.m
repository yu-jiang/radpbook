function c = times(a,b)
% function C=times(A,B)
%
% DESCRIPTION
%   Element-by-element multiply of two polynomial matrices.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C: polynomial, the result of multiplication.
%
% SYNTAX
%   C= A.*B
%     Multiplies A and B element-by-element. A and B must have the same
%     dimensions unless one is a scalar. A scalar can be multiplied by 
%     anything.
%   C = times(A,B)
%     Function-call form for times.

% 10/22/2002: PJS  Initial Coding

% Promote a to polynomial
a = polynomial(a);
sza = size(a);

% Promote b to polynomial
b = polynomial(b);
szb = size(b);

if isempty(a) || isempty(b)
    
    if isempty(a) && all(szb==[1 1])
        % empty.*scalar = empty(sza)
        c=polynomial(zeros(sza));
        return;
    elseif isempty(b) && all(sza==[1 1])
        % scalar.*empty = empty(szb)
        c=polynomial(zeros(szb));
        return;
    elseif all(sza==szb)
        c=polynomial(zeros(sza));
        return;
    else
        error('Matrix dimensions must agree.');
    end
    
elseif all(sza==szb)
    % Matrix .* Matrix
    
    % Get Dimensions
    nta = size(a.degmat,1);
    nva = length(a.varname);
    ntb = size(b.degmat,1);
    nvb = length(b.varname);
    
    if nva==0 && nvb==0
        % Handle constant.*constant
        acoef = sum(a.coefficient,1);
        bcoef = sum(b.coefficient,1);
        ccoef = reshape(acoef.*bcoef,sza);
        c = polynomial(ccoef);
        return;
    else
        % Get Coefficients
        acoef = a.coefficient;
        acoef = repmat(acoef,ntb,1);
        
        temp = b.coefficient;
        bcoef =spalloc(ntb*nta,szb(1)*szb(2),nta*nnz(temp));
        for i1 = 1:ntb
            idx = (1:nta)+(i1-1)*nta;
            bcoef(idx,:) = repmat(temp(i1,:),nta,1);
        end
        coefficient = acoef.*bcoef;
        
        % Form Degmat and Varname
        adeg = a.degmat;
        bdeg = b.degmat;
        bdegcol = [];
        for i1 = 1:ntb
            bdegcol = [bdegcol; repmat(bdeg(i1,:),nta,1)];
        end
        degmat = [repmat(adeg,ntb,1) bdegcol];
        varname = [a.varname(:); b.varname(:)];
        
        % Form polynomial and combine terms
        chkval = 0; % skip validity check
        c = polynomial(coefficient,degmat,varname,sza,chkval);
        c = combine(c);
        return
    end
    
elseif all(sza==[1 1]) || all(szb==[1 1])
    % Scalar .* Matrix or Matrix .* Scalar
    
    % Make first term be the scalar.
    if ~all(sza==[1 1])
        temp = a;
        a = b;
        b = temp;
        szb = sza;
    end
    
    nva = length(a.varname);
    if nva == 0
        % Handle constant scalar .* poly matrix
        acoef = sum(a.coefficient);
        c = b;
        c.coefficient = c.coefficient*acoef;
    else
        % Turn into matrix.*matrix
        a.coefficient = repmat(a.coefficient,[1 szb(1)*szb(2)]);
        a.matdim = szb;
        c = times(a,b);
    end
else
    error('Matrix dimensions must agree');
end












