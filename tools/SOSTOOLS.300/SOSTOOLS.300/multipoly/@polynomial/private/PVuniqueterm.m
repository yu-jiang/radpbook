function [newcoef,uniquedeg] = PVuniqueterm(coef,degmat,matdim)
% function [NewCoef,UniqueDeg] = PVuniqueterm(Coef,Degmat,MatDim);
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   Find the unique set of terms (i.e. unique rows of Degmat
%   and form coef in terms of this unique set.
%
% INPUTS
%   Coef: coefficient matrix for a polynomial
%   Degmat: degree matrix for a polynomial
%   MatDim  dimensions of the polynomial matrix
%
% OUTPUTS
%   NewCoef: coefficient matrix in terms of new degree matrix
%   UniqueDegmat: degree matrix with unique rows
%
% SYNTAX
%   [NewCoef,UniqueDeg] = PVuniqueterm(Coef,Degmat)

% 11/25/2002: PJS  Initial Coding
%  1/30/2002: PJS  Combine rows using sparse multiply
% 11/11/2008: PJS  Try to use integer sortrows for speed


% Find repeated monomials
% Sort done by monominal degree first and then lexicographic
%   sortidx: index such that degsort = degmat(sortidx)
%   repeatidx: nonzero indices indicate repeats
%   rvec:  map rows of degsort to numbers 1,..,nut
maxdeg = max(max(degmat)); %max(degmat(:));
if maxdeg==0
    nbits = 8;
else
    nbits = log2(maxdeg);
end
if  0
    try
        % Integer sorts are faster but converting degmat to a full matrix
        % might cause memory problems.
        M = [sum(degmat,2) degmat];
        if nbits < 8
            [degsort,sortidx] = sortrows( uint8(full(M)) );
        elseif nbits < 16
            [degsort,sortidx] = sortrows( uint16(full(M)) );
        elseif nbits < 32
            [degsort,sortidx] = sortrows( uint32(full(M)) );
        elseif nbits < 64
            [degsort,sortidx] = sortrows( uint64(full(M)) );
        else
            [degsort,sortidx] = sortrows(M);
        end
    catch
        [degsort,sortidx] = sortrows(M);
    end
else
    % PJS 4/24/2009-- I tried multiplying Z(x)*M*Z(x) where Z is all monoms
    % in 4 vars of deg=2:3.  Sparse sort seemed to be about 50% faster on
    % this multiply.  Is sparse sort now faster in R2009a or is it faster
    % only on this problem?
    M = [sum(degmat,2) degmat];
    [degsort,sortidx] = sortrows(M);
end
degsort = degsort(:,2:end);
% PJS 11/10/2009 -- The following line causes an out of memory error when
% adding large polys.  degsort is very sparse so computing repeatidx with
% == results in a matrix with many non-zero entries.
%repeatidx = all(degsort(2:end,:) == degsort(1:end-1,:),2);

repeatidx = ~any(degsort(2:end,:) ~= degsort(1:end-1,:),2);
rvec = cumsum([1; ~repeatidx]);

% Create unique set of monomials
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniquedeg = degmat(uniqueidx,:);

% Sum repeated terms
%  (Use sparse multiply, summat*coef, to sum the rows)
nt = size(degmat,1);
nut = size(uniquedeg,1);
summat = sparse(rvec,sortidx,1,nut,nt,nt);
newcoef = summat*coef;

% Put terms in "alphabetical" order
% (not necessary, but might look better)
uniquedeg = flipud(uniquedeg);
newcoef = flipud(newcoef);




