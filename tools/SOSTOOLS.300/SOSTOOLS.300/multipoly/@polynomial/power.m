function b = power(a,n)
% function B=power(A,N)
%
% DESCRIPTION
%   Element-by-element powers
%
% INPUTS
%   A: polynomial
%   N: matrix of natural numbers
%
% OUTPUTS
%   B: polynomial, result of element-by-element power A.^N.
%
% SYNTAX
%   B = A.^N
%     Element-by-element powers. A and B must have the same dimensions
%     unless one is a scalar. A scalar can operate into anything.
%   B=power(A,N)
%     Function-call form of power.

% 6/7/2002 PJS  Initial Coding
% 12/12/2010 PJS Call matrix.^matrix code for scalar expansion cases

% Check number of inputs
error(nargchk(2,2,nargin));
sza=size(a);
szn=size(n);

if isempty(a) || isempty(n)
    if isempty(a) && all(szn==[1 1])
        % empty.^scalar = empty(sza)
        b=polynomial(zeros(sza));
        return;
    elseif isempty(n) && all(sza==[1 1])
        % scalar.^empty = empty(szn)
        b=polynomial(zeros(szn));
        return;
    elseif all(sza==szn)
        b=polynomial(zeros(sza));
        return;
    else
        error('Matrix dimensions must agree.');
    end
    
elseif sza==szn
    % Matrix .^ Matrix
    b = zeros(sza);
    
    % exponent = 0
    idx = find( n==0 );
    b(idx) = 1;
    b = polynomial(b);
    
    % exponent = 1
    idx = find( n==1 );
    L.type = '()';
    L.subs = {idx};
    b = subsasgn(b,L,subsref(a,L));
    
    % exponent = 2;
    idx = find( n>1 );
    idx = idx(:)';
    for i1=idx
        L.subs = {i1};
        ai = subsref(a,L);
        Nt = size(ai.coefficient,1);
        if Nt==1
            bi = ai;
            bi.coefficient = bi.coefficient^n(i1);
            bi.degmat = bi.degmat*n(i1);
        else
            bi = mpower(ai,n(i1));
        end
        b = subsasgn(b,L,bi);
    end
    
    %         % XXX This is pretty slow.
    %         b = a;
    %         L.type = '()';
    %         for i1 = 1:sza(1);
    %             for i2 = 1:sza(2);
    %                 L.subs = {i1,i2};
    %                 aij = subsref(a,L);
    %                 bij = mpower(aij,n(i1,i2));
    %                 b = subsasgn(b,L,bij);
    %             end
    %         end
    %     end
elseif all(szn==[1 1])
    % Matrix.^Scalar:
    n = repmat(n,sza);
    b = power(a,n);
    
    %     b = polynomial(1);
    %     for i1 = 1:n
    %         b = b.*a;
    %     end
elseif all(sza==[1 1])
    % Scalar .^ Matrix:
    a = repmat(a,szn);
    b = power(a,n);
    
    %     % XXX This is pretty slow.
    %     b = polynomial(zeros(szn));
    %     L.type = '()';
    %     for i1 = 1:szn(1);
    %         for i2 = 1:szn(2);
    %             L.subs = {i1,i2};
    %             bij = mpower(a,n(i1,i2));
    %             b = subsasgn(b,L,bij);
    %         end
    %     end
else
    error('Matrix dimensions must agree');
end





