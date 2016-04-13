function varargout=prod(A,dim)
% function C=prod(A,Dim)
%
% DESCRIPTION
%   Product of the elements of a polynomial vector or matrix.
%
% INPUTS
%   A: polynomial vector or matrix
%   Dim: dimension along which product is to be performed
%      (Default: Dim = 1)
%
% OUTPUTS
%   C: polynomial or polynomial vector, the result of the product
%
% SYNTAX
%   C = prod(A)
%     If A is a vector then C is the product of the elements in A.  If A is
%     a matrix then C is a row vector with the products of columns of A.
%   C = sum(A,1)
%     C is a row vector with the product over each column of A.
%   C = sum(A,2)
%     C is a column vector with the product over each row of A.

% 12/11/10 PJS  Initial Coding

% Promote A
A = polynomial(A);
szA = size(A);

% Handle empty A
if isempty(A)
    if nargin==1 || isempty(dim)
        C = polynomial(1);
    elseif dim==1
        C = polynomial(zeros(1,0));
    else
        C = polynomial(zeros(0,1));
    end    
    varargout{1} = C;
    return;
end

% Handle non-empty A
if nargin==2 && isempty(dim)
    dim = 1;
end
zerop = polynomial(0);
if (any(szA == 1) && nargin==1) || (szA(2)==1 && dim==1 ) ...
        || (szA(1)==1 && dim==2 )
    % Prod of a vector
    C = polynomial(1);
    L.type = '()';
    for i1=1:length(A)
        L.subs = {i1};
        Ai = subsref(A,L);
        Nt = size(Ai.coefficient,1);
        if isequal(Ai,zerop)
            varargout{1} = polynomial(0);
            return
        else
            C = C*Ai;
        end
    end
    C = combine(C);
elseif nargin == 1 || dim == 1
    % Prod for a matrix along cols of A (dim=1)
    C = prod(A',2)';
else
    % Prod for a matrix along rows of A (dim=2)    
    % XXX Code is not vectorized so computation will be slow.
    C = polynomial( zeros(szA(1),1) );
    for i1 = 1:size(A,1)
        s.type = '()';
        s.subs = {i1,':'};
        Ai = subsref(A,s);
        
        s.subs = {i1};
        C =  subsasgn(C,s,prod(Ai));
    end
end

varargout{1} = C;

