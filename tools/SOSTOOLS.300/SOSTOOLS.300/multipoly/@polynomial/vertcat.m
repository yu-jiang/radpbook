function c = vertcat(varargin)
% function C=vertcat(A,B);
%
% DESCRIPTION
%   Vertical concatenation of polynomial objects.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C:  vertical concatenation of input matrices.
%
% SYNTAX
%   [A; B]
%      Vertical concatenation of polynomial matrices A and B.
%      A and B must have the same number of columns.
%   [A1; A2; A3; ...]
%      Vertical concatenation of several polynomial matrices.
%   C = vertcat(A1,A2,A3,..);
%      Function-call form of vertical concatenation.
%
% See also horzcat

% 6/8/2002: PJS  Initial Coding

if nargin==1
    c = varargin{1};
else
    % Promote a to polynomial
    a = polynomial(varargin{1});
    [nra,nca] = size(a);
    
    % Promote b to polynomial
    b = polynomial(varargin{2});
    [nrb,ncb] = size(b);
    
    if isempty(b);
        c = a;
    elseif isempty(a);
        c = b;
    elseif nca == ncb
        % Get Dimensions
        nta = size(a.degmat,1);
        nva = length(a.varname);
        ntb = size(b.degmat,1);
        nvb = length(b.varname);
        
        if nva==0 && nvb==0
            % Combine constant terms
            ar = combine(a);
            coef1 = reshape(ar.coefficient,[nra nca]);
            br = combine(b);
            coef2 = reshape(br.coefficient,[nrb ncb]);
            
            % Stack up Coefficients and Form Polynomial
            coefficient = [coef1; coef2];
            c = polynomial(coefficient);
        else
            % Form Degmat, Varname, and Matdim
            adeg = a.degmat;
            bdeg = b.degmat;
            degmat = blkdiag(adeg,bdeg);
            varname = [a.varname(:); b.varname(:)];
            matdim = [nra+nrb nca];
            
            % Stack up Coefficients
            idx1 = [];
            idx2 = [];
            for i1 = 0:(nca-1);
                idx1 = [idx1 (1:nra)+i1*(nra+nrb)];
                idx2 = [idx2 nra+(1:nrb)+i1*(nra+nrb)];
            end;
            coef1 = spalloc(nta,(nra+nrb)*nca,nnz(a.coefficient));
            coef1(:,idx1) = a.coefficient;
            coef2 = spalloc(ntb,(nra+nrb)*ncb,nnz(b.coefficient));
            coef2(:,idx2) = b.coefficient;
            coefficient = [coef1; coef2];
            
            % Form Polynomial and combine terms
            chkval = 0; % skip validity check
            c = polynomial(coefficient,degmat,varname,matdim,chkval);
            c = combine(c);
        end
        
    else
        error('All columns must have the same row dimension')
    end
    
    if nargin>2
        c = vertcat(c,varargin{3:end});
    end
end


