function b = mpower(a,n)
% function B=mpower(A,N)
%
% DESCRIPTION
%   Matrix power, i.e. A to the power N
%
% INPUTS
%   A: polynomial
%   N: natural number
%
% OUTPUTS
%   B: polynomial, result of matrix power A^N
%
% SYNTAX
%   B = A^N
%     A to the power N
%   B=mpower(A,N)
%     Function-call form for mpower

% 6/7/2002: PJS  Initial Coding

% Check number of inputs
error(nargchk(2,2,nargin));

if size(a,1)~=size(a,2)
    error('Matrix must be square.');
end

if isempty(a) || isempty(n)
    if isempty(a) && max(size(n))==1
        % empty^scalar = empty
        b=polynomial;
        return;
    elseif isempty(n) && max(size(a))==1
        % scalar^empty = empty
        b=polynomial;
        return;
    else
        % matrix^empty  or empty^matrix  or empty^empty
        error('At least one operand must be scalar.');
    end
elseif ndims(n)==2 && max(size(n)) == 1 && ...
        floor(n)==ceil(n) && n >= 0
    
    % Call mpower recursively
    % Note: A for loop on mtimes requires n calls to mtimes. A recursion
    % only requires roughly log2(n) calls to mtimes.
    if n==0
        b = polynomial(eye(size(a,1)));
        return;
    elseif n==1
        b=a;
        return;
    elseif floor(n/2)==ceil(n/2)
        b = mpower(a,n/2);
        b = mtimes(b,b);
    else
        b = mpower(a,(n-1)/2);
        b = mtimes(b,b);
        b = mtimes(b,a);
    end
else
    error('Invalid use of power for polynomials');
end
