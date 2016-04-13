function b = PVsubsasgn_1idx(a,L,RHS)
% function B =  PVsubsasgn_1idx(A,L,RHS)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a 1x1 structure with the fields:
%    type -- string containing '()'
%    subs -- 1x1 cell array containing indices
%   RHS: Value to be assigned.
%
% OUTPUTS
%   B: object after subsassignment
%
% SYNTAX
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding
% 11/11/2008: PJS Mods to reduce computation

% Check Indices
subsidx = L(1).subs{1};
if isempty(subsidx)
    b = a;
    return;
end

ls = length(subsidx);
if isa(subsidx,'logical')
    subsidx = find(subsidx);
elseif min(subsidx)<1
    error('Index into matrix is negative or zero.');
end

% Get info about polynomials
b = a;
[nrb,ncb]=size(b);

temp = RHS;
[nrt,nct]=size(temp);
ntt = size(temp.degmat,1);

if all([nrt nct]==[1 1])
    % Fan out scalar RHS to correct dimension
    temp.coefficient = repmat(temp.coefficient,1,ls);
    nrt = 1;
    nct = ls;
    temp.matdim = [nrt nct];
elseif all([nrt nct]==[0 0])
    % Assignment with empty removes entries from b
    b.coefficient(:,subsidx) = [];
    if all(size(b)==[1 1])
        b.matdim = [1 0];
    elseif size(b,2)==1
        b.matdim = [ size(b.coefficient,2) 1];
    else
        b.matdim = [ 1 size(b.coefficient,2)];
    end
    b = combine(b);
    return;
end

% Check dimensions
if ls~=nrt*nct
    error(['In an assignment  A(I) = B, the number of' ...
        ' elements in B and I must be the same.']);
end
if (min([nrb ncb])>1) &&  (max(subsidx) > nrb*ncb)
    error('In an assignment  A(I) = B, a matrix A cannot be resized.');
end

% Determine proper dimensions for poly after the assignment
if isempty(b)
    nc = max(subsidx(:));
    nr = 1;
elseif min([nrb ncb])==1 && ncb >= nrb
    nc = max([ncb; subsidx(:)]);
    nr = 1;
elseif min([nrb ncb])==1
    nc = 1;
    nr = max([nrb; subsidx(:)]);
else
    nc = ncb;
    nr = nrb;
end

% Perform the subsassignment
if isempty(b)
    tempcoef = sparse(ntt,nr*nc);
    tempcoef(:,subsidx) = temp.coefficient;
    temp.coefficient = tempcoef;
    temp.matdim = [nr nc];
    
    b=temp;
else
    % Zero out terms of b, stack in terms of temp, and then combine
    coef = b.coefficient;
    coef(:,subsidx) = 0;
    coef(end+1:end+ntt,subsidx) = temp.coefficient;
    
    degmat = blkdiag(b.degmat,temp.degmat);
    varname = [b.varname; temp.varname];
    matdim = [nr nc];
    chkval = 0; % skip validity check
    b = polynomial(coef,degmat,varname,matdim,chkval);
    b = combine(b);
end


