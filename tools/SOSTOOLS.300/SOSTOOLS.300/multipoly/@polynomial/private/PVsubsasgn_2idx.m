function b = PVsubsasgn_2idx(a,L,RHS)
% function B =  PVsubsasgn_2idx(A,L,RHS)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a 1x1 structure with the fields:
%    type -- string containing '()'
%    subs -- 1x2 cell array containing indices
%   RHS: Value to be assigned.
%
% OUTPUTS
%   B: object after subsassignment
%
% SYNTAX
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding
% 11/11/2008: PJS Mods to reduce computation
% 1/13/2011: PJS  Fixed bug in sub2ind conversion

% Check Indices
ridx = L(1).subs{1};
lr = length(ridx);
cidx = L(1).subs{2};
lc = length(cidx);
if min(ridx)<1 || min(cidx)<1
    error('Index into matrix is negative or zero.');
end

% Get info about polynomials
b = a;
[nrb,ncb]=size(b);
ntb = size(b.degmat,1);

temp = RHS;
[nrt,nct]=size(temp);
ntt = size(temp.degmat,1);

% Fan out scalar RHS to correct dimension
if all([nrt nct]==[1 1])
    temp.coefficient = repmat(temp.coefficient,1,lr*lc);
    nrt = lr;
    nct = lc;
    temp.matdim = [nrt nct];
end

% Check dimensions
if lr*lc~=nrt*nct
    error(['In an assignment  A(I) = B, the number of' ...
        ' elements in B and I must be the same.']);
end

% Determine proper dimensions for poly after the assignment
nr = max([nrb; ridx(:)]);
nc = max([ncb; cidx(:)]);

% Convert from row/col indices to single index
if isempty(ridx) || isempty(cidx)
    % No substitution
    b = a;
    return;
end
% XXX PJS If ridx and cidx are N-by-1 then the commented code produces
% subsidx of size N-by-1. However, it should be N^2-by-1 (all combos)
% elseif length(cidx)==1
%     cidx = repmat(cidx,length(ridx),1);
% elseif length(ridx)==1
%     ridx = repmat(ridx,length(cidx),1);
% end
% subsidx = sub2ind([nr nc],ridx(:),cidx(:));
tmp = reshape(1:nr*nc,[nr,nc]);
subsidx = tmp(ridx,cidx);
subsidx = subsidx(:);

% Perform the subsassignment
if isempty(b)
    tempcoef = sparse(ntt,nr*nc);
    tempcoef(:,subsidx) = temp.coefficient;
    temp.coefficient = tempcoef;
    temp.matdim = [nr nc];
    
    b = temp;
else
    % Place coefs of b in proper location, zero out terms of b to be
    % assigned, stack in new terms, and then combine
    coef = sparse(ntb,nr*nc);
    [ii,jj]=ind2sub([nrb ncb],1:(nrb*ncb));
    bidx = sub2ind([nr nc],ii,jj);
    coef(:,bidx) = b.coefficient;
    
    coef(:,subsidx) = 0;
    coef(end+1:end+ntt,subsidx) = temp.coefficient;
    
    degmat = blkdiag(b.degmat,temp.degmat);
    varname = [b.varname; temp.varname];
    matdim = [nr nc];
    chkval = 0; % skip validity check
    b = polynomial(coef,degmat,varname,matdim,chkval);
    b = combine(b);
end







