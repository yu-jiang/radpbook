function c = mtimes(a,b)
% function C=mtimes(A,B)
%
% DESCRIPTION
%   Multiply two polynomial matrices and combine terms.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C: polynomial, the result of the multiplication A*B
%
% SYNTAX
%   C= A*B
%     Multiplies A and B. A scalar can be multiplied by anything. Otherwise
%     the number of columns of A must be equal to the number of rows of B.
%   C = mtimes(A,B)
%     Function-call form for mtimes.

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
        % empty*scalar = empty(sza)
        c=polynomial(zeros(sza));
        return;
    elseif isempty(b) && all(sza==[1 1])
        % scalar*empty = empty(szb)
        c=polynomial(zeros(szb));
        return;
    elseif sza(2)==szb(1)
        c=polynomial(zeros(sza(1),szb(2)));
        return;
    else
        error('Matrix dimensions must agree.');
    end
    
elseif all(sza==[1 1]) || all(szb==[1 1])
    % Scalar * Matrix or Matrix * Scalar
    
    % Make first term be the scalar.
    if ~all(sza==[1 1])
        temp = a;
        a = b;
        
        b = temp;
        szb = sza;
    end
    
    nva = length(a.varname);
    if nva == 0
        % Handle constant scalar * poly matrix
        acoef = sum(a.coefficient);
        c = b;
        c.coefficient = c.coefficient*acoef;
    else
        % Turn into matrix.*matrix
        a.coefficient = repmat(a.coefficient,[1 szb(1)*szb(2)]);
        a.matdim = szb;
        c = times(a,b);
    end
    
elseif sza(2)==szb(1)
    
    % Matrix * Matrix
    
    % Get Dimensions
    nta = size(a.degmat,1);
    nva = length(a.varname);
    ntb = size(b.degmat,1);
    nvb = length(b.varname);
    
    if nva==0 && nvb==0
        % Handle constant case
        acoef = reshape(sum(a.coefficient,1),sza);
        bcoef = reshape(sum(b.coefficient,1),szb);
        c = polynomial(acoef*bcoef);
        return;
        
    else
        % Get Coefficients:
        %  acoefcol: coef. matrices stacked vertically
        %  bcoefcol: coef. matrices stacked horizontally
        
        %         acoef = a.coefficient;
        %         acoefcol = spalloc(nta*sza(1),sza(2),nnz(acoef));
        %         for i1 = 1:nta
        %             ridx = (i1-1)*sza(1) + (1:sza(1));
        %             acoefcol(ridx,:) = reshape(acoef(i1,:),sza(1),sza(2));
        %         end
        %
        %         bcoef = b.coefficient;
        %         bcoefcol = spalloc(szb(1),ntb*szb(2),nnz(bcoef));
        %         for i1 = 1:ntb
        %             cidx = (i1-1)*szb(2) + (1:szb(2));
        %             bcoefcol(:,cidx) = reshape(bcoef(i1,:),szb(1),szb(2));
        %         end
        %         tempcoef = acoefcol*bcoefcol;
        %
        %         coefficient = spalloc(sza(1),nta*ntb*szb(2),nnz(tempcoef));
        %         for i1 = 1:nta
        %             ridx = (i1-1)*sza(1) + (1:sza(1));
        %             cidx = (i1-1)*ntb*szb(2) + (1:ntb*szb(2));
        %             coefficient(:,cidx) = tempcoef(ridx,:);
        %         end
        %         coefficient = reshape(coefficient,sza(1)*szb(2),nta*ntb);
        %         coefficient = coefficient';
        
        % Vectorized code to compute coef matrix
        idx1 = reshape(1:sza(1)*sza(2),[sza(1),sza(2)]);
        idx1 = repmat(idx1,[nta 1]);
        idx1 = idx1(:);
        
        idx2 = repmat(1:nta,[sza(1) 1]);
        idx2 = repmat(idx2(:),[sza(2) 1]);
        
        idx = sub2ind([sza(1)*sza(2) nta],idx1,idx2);
        acoef = a.coefficient;
        acoefcol = full(acoef');
        acoefcol = reshape(acoefcol(idx),[nta*sza(1) sza(2)]);
        
        bcoef = b.coefficient;
        bcoefcol = reshape(bcoef',[szb(1) szb(2)*ntb]);
        
        tempcoef = acoefcol*bcoefcol;
        
        idx1 = reshape(1:sza(1)*nta,[sza(1),nta]);
        idx1 = repmat(idx1,[ntb*szb(2) 1]);
        idx1 = idx1(:);
        
        idx2 = repmat(1:(ntb*szb(2)),[sza(1) 1]);
        idx2 = repmat(idx2(:),[nta 1]);
        
        idx = sub2ind([nta*sza(1) ntb*szb(2)],idx1,idx2);
        coefficient = reshape(tempcoef(idx),sza(1)*szb(2),nta*ntb);
        coefficient = coefficient';
        
        % Form Degmat and Varname
        if nva==0
            bdeg = b.degmat(:,1:nvb);
            bdeg = repmat(bdeg,nta,1);
            degmat = bdeg;
            varname = b.varname(:);
        else
            tempdeg = a.degmat(:,1:nva);
            %             adeg = spalloc(ntb*nta,nva,ntb*nnz(tempdeg));
            %             for i1 = 1:nta
            %                 ridx = (i1-1)*ntb + (1:ntb);
            %                 adeg(ridx,:) = repmat(tempdeg(i1,:),[ntb,1]);
            %             end
            
            % PJS 4/23/09
            % Copying rows of sparse matrices can greatly increase the
            % memory usage over copying columns.  There is a technical
            % solution description at mathworks.com:
            % "Why does MATLAB version 6.5 (R13) allocate large amounts of
            % memory when I copy a row of data from a sparse matrix?"
            % The vectorized code below does a copy of columns and
            % then transposes to ensure less memory usage.  I also
            % discovered that double transposing a sparse matrix (as in
            % the commented code below) is a fast trick to reduce
            % the memory usage of a sparse matrix.  What is this actually
            % doing to the sparse matrix object?
            allidx=repmat(1:nta,[ntb,1]);
            allidx = allidx(:);
            adeg = tempdeg';
            adeg = (adeg(:,allidx))';
            
            bdeg = b.degmat(:,1:nvb);
            bdeg = repmat(bdeg,nta,1);
            degmat = [adeg bdeg];
            varname = [a.varname(:); b.varname(:)];
        end
        
        % Form polynomial and combine terms
        chkval = 0; % skip validity check
        c = polynomial(coefficient,degmat,varname,[sza(1) szb(2)],chkval);
        c = combine(c);
    end
    
else
    error('Inner matrix dimensions must agree');
end

